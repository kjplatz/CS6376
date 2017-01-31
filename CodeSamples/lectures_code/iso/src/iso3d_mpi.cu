/*
 * Copyright (c) 2012, NVIDIA CORPORATION. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  * Neither the name of NVIDIA CORPORATION nor the names of its
 *    contributors may be used to endorse or promote products derived 
 *    from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 * OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * Paulius Micikevicius (pauliusm@nvidia.com)
 * Max Grossman (jmaxg3@gmail.com)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <unistd.h>
#include <mpi.h>
#include "common.h"
#include "common3d.h"

#define BDIMX   32
#define BDIMY   16
#define SHAREDX(radius) (BDIMX + 2 * (radius))
#define SHAREDY(radius) (BDIMY + 2 * (radius))
#define CACHE_INDEX(y, x, radius)   ((y) * SHAREDX(radius) + (x))

#define CHECK_MPI(call) {                                                      \
    int mpi_err = call;                                                        \
    if (mpi_err != MPI_SUCCESS) {                                              \
        fprintf(stderr, "MPI Error at %s:%d - %d\n", __FILE__, __LINE__,       \
                mpi_err);                                                      \
        exit(1);                                                               \
    }                                                                          \
}

#define LOW_HALO_START_Z(gpu, nz_per_gpu) \
    ((gpu) * (nz_per_gpu))

__constant__ TYPE const_c_coeff[NUM_COEFF];

__global__ void fwd_kernel(TYPE *next, TYPE *curr, TYPE *vsq,
        int nx, int ny, int nz, int dimx, int dimy, int radius) {
    extern __shared__ TYPE cache[];

    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    const int this_y = radius + threadIdx.y;
    const int this_x = radius + threadIdx.x;

    // radius must be <= TRANSACTION_LEN
    TYPE zneighbors[2 * TRANSACTION_LEN + 1];

    for (int z = -radius, index = 1; z < radius; z++, index++) {
        zneighbors[index] = curr[POINT_OFFSET(x, y, z, dimy, dimx, radius)];
    }

    for (int z = 0; z < nz; z++) {
        int this_offset = POINT_OFFSET(x, y, z, dimy, dimx, radius);

        for (int i = 0; i <= 2 * radius; i++) {
            zneighbors[i] = zneighbors[i + 1];
        }
        zneighbors[2 * radius] =
            curr[POINT_OFFSET(x, y, z + radius, dimy, dimx, radius)];

        __syncthreads();

        cache[CACHE_INDEX(this_y, this_x, radius)] = curr[POINT_OFFSET(x, y, z, dimy, dimx, radius)];
        if (threadIdx.y < radius) {
            cache[CACHE_INDEX(threadIdx.y, this_x, radius)] =
                curr[POINT_OFFSET(x, y - radius, z, dimy, dimx, radius)];
        }
        if (threadIdx.y >= radius && threadIdx.y < 2 * radius) {
            cache[CACHE_INDEX(threadIdx.y + blockDim.y, this_x, radius)] =
                curr[POINT_OFFSET(x, y - radius + blockDim.y, z, dimy, dimx, radius)];
        }
        if (threadIdx.x < radius) {
            cache[CACHE_INDEX(this_y, threadIdx.x, radius)] =
                curr[POINT_OFFSET(x - radius, y, z, dimy, dimx, radius)];
        }
        if (threadIdx.x >= radius && threadIdx.x < 2 * radius) {
            cache[CACHE_INDEX(this_y, threadIdx.x + blockDim.x, radius)] =
                curr[POINT_OFFSET(x - radius + blockDim.x, y, z, dimy, dimx, radius)];
        }

        __syncthreads();

        TYPE temp = 2.0f * cache[CACHE_INDEX(this_y, this_x, radius)] - next[this_offset];
        TYPE div = const_c_coeff[0] * cache[CACHE_INDEX(this_y, this_x, radius)];
        for (int d = radius; d >= 1; d--) {
            div += const_c_coeff[d] * (zneighbors[radius + d] +
                    zneighbors[radius - d] + cache[CACHE_INDEX(this_y + d, this_x, radius)] +
                    cache[CACHE_INDEX(this_y - d, this_x, radius)] + cache[CACHE_INDEX(this_y, this_x + d, radius)] +
                    cache[CACHE_INDEX(this_y, this_x - d, radius)]);
        }
        next[this_offset] = temp + div * vsq[this_offset];
    }
}

static void enable_peer_access(int src, int target) {
    int can_access;
    
    CHECK(cudaSetDevice(src));
    CHECK(cudaDeviceCanAccessPeer(&can_access, src, target));
    if (can_access) {
        CHECK(cudaDeviceEnablePeerAccess(target, 0));
    } else {
        printf("WARNING: peer access from %d -> %d is unsupported\n", src,
                target);
    }
}

static void exchange_data(int rank, int nprocs, MPI_Datatype mpi_data_type,
        TYPE *high_send_buff, TYPE *high_recv_buff, TYPE *low_send_buff,
        TYPE *low_recv_buff, TYPE **d_next, int ngpus_per_node, int nz_per_gpu,
        int dimy, int dimx, cudaStream_t *streams, int radius) {
    MPI_Request low_send_request, low_recv_request, high_send_request,
                high_recv_request;

    if (rank < nprocs - 1) {
        // copy halo up to the next rank and receive from it
        CHECK(cudaMemcpyAsync(high_send_buff,
                    d_next[ngpus_per_node - 1] + (nz_per_gpu * dimy *
                        dimx),
                    radius * dimy * dimx * sizeof(TYPE),
                    cudaMemcpyDeviceToHost,
                    streams[ngpus_per_node - 1]));
        CHECK(cudaStreamSynchronize(streams[ngpus_per_node - 1]));
        CHECK_MPI(MPI_Isend(high_send_buff, radius * dimy * dimx,
                    mpi_data_type, rank + 1, rank, MPI_COMM_WORLD,
                    &high_send_request));
        CHECK_MPI(MPI_Irecv(high_recv_buff, radius * dimy * dimx,
                    mpi_data_type, rank + 1, rank + 1, MPI_COMM_WORLD,
                    &high_recv_request));
    }
    if (rank > 0) {
        // copy halo down to the previous rank and receive from it
        CHECK(cudaMemcpyAsync(low_send_buff,
                    d_next[0] + (radius * dimy * dimx),
                    radius * dimy * dimx * sizeof(TYPE),
                    cudaMemcpyDeviceToHost, streams[0]));
        CHECK(cudaStreamSynchronize(streams[0]));
        CHECK_MPI(MPI_Isend(low_send_buff, radius * dimy * dimx,
                    mpi_data_type, rank - 1, rank, MPI_COMM_WORLD,
                    &low_send_request));
        CHECK_MPI(MPI_Irecv(low_recv_buff, radius * dimy * dimx,
                    mpi_data_type, rank - 1, rank - 1, MPI_COMM_WORLD,
                    &low_recv_request));
    }

    for (int g = 0; g < ngpus_per_node; g++) {
        if (g < ngpus_per_node - 1) {
            // copy up
            CHECK(cudaMemcpyAsync(d_next[g + 1],
                        d_next[g] + (nz_per_gpu * dimy * dimx),
                        radius * dimy * dimx * sizeof(TYPE),
                        cudaMemcpyDeviceToDevice, streams[g]));
        }
        if (g > 0) {
            // copy down
            CHECK(cudaMemcpyAsync(d_next[g - 1] +
                            ((radius + nz_per_gpu) * dimy * dimx),
                        d_next[g] + (radius * dimy * dimx),
                        radius * dimy * dimx * sizeof(TYPE),
                        cudaMemcpyDeviceToDevice, streams[g]));
        }
    }

    if (rank < nprocs - 1) {
        CHECK_MPI(MPI_Wait(&high_send_request, MPI_STATUS_IGNORE));
        CHECK_MPI(MPI_Wait(&high_recv_request, MPI_STATUS_IGNORE));
        CHECK(cudaMemcpyAsync(d_next[ngpus_per_node - 1] +
                        ((radius + nz_per_gpu) * dimy * dimx),
                    high_recv_buff, radius * dimy * dimx * sizeof(TYPE),
                    cudaMemcpyHostToDevice,
                    streams[ngpus_per_node - 1]));
    }
    if (rank > 0) {
        CHECK_MPI(MPI_Wait(&low_send_request, MPI_STATUS_IGNORE));
        CHECK_MPI(MPI_Wait(&low_recv_request, MPI_STATUS_IGNORE));
        CHECK(cudaMemcpyAsync(d_next[0], low_recv_buff,
                    radius * dimy * dimx * sizeof(TYPE),
                    cudaMemcpyHostToDevice, streams[0]));
    }
}

int main( int argc, char *argv[] ) {
    char name[MPI_MAX_PROCESSOR_NAME];
    int rank, nprocs, namelen;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(name, &namelen);
    MPI_Datatype mpi_data_type;
    if (sizeof(TYPE) == sizeof(float)) {
        mpi_data_type = MPI_FLOAT;
    } else if (sizeof(TYPE) == sizeof(double)) {
        mpi_data_type = MPI_DOUBLE;
    } else {
        fprintf(stderr, "Unrecognized TYPE\n");
        return 1;
    }

    fprintf(stderr, "MPI proc %d/%d starting on %s...\n", rank + 1, nprocs, name);

    config conf;
    setup_config(&conf, argc, argv);
    if (rank == 0) {
        init_progress(conf.progress_width, conf.nsteps, conf.progress_disabled);
    }

#ifndef PADDING
    fprintf(stderr, "Must be compiled with -DPADDING\n");
    return 1;
#endif

    if (conf.nx % BDIMX != 0) {
        fprintf(stderr, "Invalid nx configuration, must be an even multiple of "
                "%d\n", BDIMX);
        return 1;
    }
    if (conf.ny % BDIMY != 0) {
        fprintf(stderr, "Invalid ny configuration, must be an even multiple of "
                "%d\n", BDIMY);
        return 1;
    }
    if (conf.radius > TRANSACTION_LEN) {
        fprintf(stderr, "Radius must be less than TRANSACTION_LEN to include "
                "it in dimx padding\n");
        return 1;
    }

    // Assume that all nodes have the same number of GPUs
    int total_ngpus = nprocs * conf.ngpus;
    int node_base_gpu = rank * conf.ngpus;
    
    if (conf.nz % total_ngpus != 0) {
        fprintf(stderr, "nz must be evenly divisible by ngpus (%d)\n",
                conf.ngpus);
        return 1;
    }
    int nz_per_gpu = conf.nz / total_ngpus;

    TYPE dx = 20.f;
    TYPE dt = 0.002f;

    // compute the pitch for perfect coalescing
    size_t dimx = TRANSACTION_LEN + conf.nx + conf.radius;
    dimx += (TRANSACTION_LEN - (dimx % TRANSACTION_LEN));
    size_t dimy = conf.ny + 2*conf.radius;
    size_t dimz_total = conf.nz + 2*conf.radius;
    size_t dimz_per_gpu = nz_per_gpu + 2 * conf.radius;
    size_t nbytes_total = dimx * dimy * dimz_total * sizeof(TYPE);
    size_t nbytes_per_gpu = dimx * dimy * dimz_per_gpu * sizeof(TYPE);

    if (conf.verbose) {
        printf("x = %zu, y = %zu, z = %zu\n", dimx, dimy, dimz_total);
        printf("nsteps = %d\n", conf.nsteps);
        printf("radius = %d\n", conf.radius);
    }

    TYPE c_coeff[NUM_COEFF];
    TYPE *curr, *next, *vsq;
    CHECK(cudaMallocHost((void **)&curr, nbytes_total));
    CHECK(cudaMallocHost((void **)&next, nbytes_total));
    CHECK(cudaMallocHost((void **)&vsq, nbytes_total));

    TYPE *low_send_buff, *low_recv_buff, *high_send_buff, *high_recv_buff;
    CHECK(cudaMallocHost((void **)&low_send_buff,
                conf.radius * dimy * dimx * sizeof(TYPE)));
    CHECK(cudaMallocHost((void **)&low_recv_buff,
                conf.radius * dimy * dimx * sizeof(TYPE)));
    CHECK(cudaMallocHost((void **)&high_send_buff,
                conf.radius * dimy * dimx * sizeof(TYPE)));
    CHECK(cudaMallocHost((void **)&high_recv_buff,
                conf.radius * dimy * dimx * sizeof(TYPE)));

    config_sources(&conf.srcs, &conf.nsrcs, conf.nx, conf.ny, conf.nsteps);
    TYPE **srcs = sample_sources(conf.srcs, conf.nsrcs, conf.nsteps, dt);

    init_data(curr, next, vsq, c_coeff, dimx, dimy, dimz_total, dx, dt);

    TYPE **d_curr, **d_next, **d_vsq;
    d_curr = (TYPE **)malloc(sizeof(TYPE *) * conf.ngpus);
    d_next = (TYPE **)malloc(sizeof(TYPE *) * conf.ngpus);
    d_vsq = (TYPE **)malloc(sizeof(TYPE *) * conf.ngpus);
    cudaStream_t *streams = (cudaStream_t *)malloc(sizeof(cudaStream_t) *
            conf.ngpus);
    for (int g = 0; g < conf.ngpus; g++) {
        CHECK(cudaSetDevice(g));
        CHECK(cudaMalloc((void **)(d_curr + g), nbytes_per_gpu));
        CHECK(cudaMalloc((void **)(d_next + g), nbytes_per_gpu));
        CHECK(cudaMalloc((void **)(d_vsq + g), nbytes_per_gpu));
        CHECK(cudaStreamCreate(streams + g));

        int device_below = g - 1;
        int device_above = g + 1;
        if (device_below >= 0) {
            enable_peer_access(g, device_below);
        }
        if (device_above < conf.ngpus) {
            enable_peer_access(g, device_above);
        }
    }

    dim3 block(BDIMX, BDIMY);
    dim3 grid(conf.nx / block.x, conf.ny / block.y);

    double mem_start = seconds();

    for (int g = 0; g < conf.ngpus; g++) {
        int low_halo_start = LOW_HALO_START_Z(node_base_gpu + g, nz_per_gpu);
        CHECK(cudaSetDevice(g));
        CHECK(cudaMemcpyToSymbol(const_c_coeff, c_coeff,
                    NUM_COEFF * sizeof(TYPE)));

        CHECK(cudaMemcpyAsync(d_curr[g], curr + (low_halo_start * dimy * dimx),
                    nbytes_per_gpu, cudaMemcpyHostToDevice, streams[g]));
        CHECK(cudaMemcpyAsync(d_next[g], next + (low_halo_start * dimy * dimx),
                    nbytes_per_gpu, cudaMemcpyHostToDevice, streams[g]));
        CHECK(cudaMemcpyAsync(d_vsq[g], vsq + (low_halo_start * dimy * dimx),
                    nbytes_per_gpu, cudaMemcpyHostToDevice, streams[g]));
    }

    double start = seconds();
    for (int step = 0; step < conf.nsteps; step++) {
        if (rank == 0) {
            for (int src = 0; src < conf.nsrcs; src++) {
                if (conf.srcs[src].t > step) continue;
                int src_offset = POINT_OFFSET(conf.srcs[src].x,
                        conf.srcs[src].y, 0, dimy, dimx, conf.radius);
                CHECK(cudaMemcpy(d_curr[0] + src_offset, srcs[src] + step,
                            sizeof(TYPE), cudaMemcpyHostToDevice));
            }
        }

        for (int g = 0; g < conf.ngpus; g++) {
            CHECK(cudaSetDevice(g));
            fwd_kernel<<<grid, block, SHAREDY(conf.radius) * SHAREDX(conf.radius) *
                sizeof(TYPE), streams[g]>>>(d_next[g], d_curr[g], d_vsq[g],
                        conf.nx, conf.ny, nz_per_gpu, dimx, dimy, conf.radius);
        }

        if (step < conf.nsteps - 1) {
            exchange_data(rank, nprocs, mpi_data_type, high_send_buff,
                    high_recv_buff, low_send_buff, low_recv_buff, d_next,
                    conf.ngpus, nz_per_gpu, dimy, dimx, streams, conf.radius);

            for (int g = 0; g < conf.ngpus; g++) {
                CHECK(cudaSetDevice(g));
                CHECK(cudaStreamSynchronize(streams[g]));
            }
        }

        TYPE **tmp = d_next;
        d_next = d_curr;
        d_curr = tmp;

        if (rank == 0) {
            update_progress(step + 1);
        }
    }

    for (int g = 0; g < conf.ngpus; g++) {
        CHECK(cudaSetDevice(g));
        CHECK(cudaDeviceSynchronize());
    }
    double compute_s = seconds() - start;

    for (int g = 0; g < conf.ngpus; g++) {
        CHECK(cudaMemcpy(curr + ((conf.radius + ((node_base_gpu + g) * nz_per_gpu)) *
                        dimy * dimx),
                    d_curr[g] + (conf.radius * dimy * dimx),
                    nz_per_gpu * dimy * dimx * sizeof(TYPE),
                    cudaMemcpyDeviceToHost));
    }

    for (int p = 1; p < nprocs; p++) {
        int z_offset = conf.radius + (p * conf.ngpus * nz_per_gpu);
        int offset = z_offset * dimy * dimx;
        if (rank == 0) {
            CHECK_MPI(MPI_Recv(curr + offset, conf.ngpus * nz_per_gpu * dimy *
                        dimx, mpi_data_type, p, p, MPI_COMM_WORLD,
                        MPI_STATUS_IGNORE));
        } else if (rank == p) {
            CHECK_MPI(MPI_Send(curr + offset, conf.ngpus * nz_per_gpu * dimy *
                        dimx, mpi_data_type, 0, p, MPI_COMM_WORLD));
        }
    }

    double total_s = seconds() - mem_start;

    if (rank == 0) {
        finish_progress();

        float point_rate = (float)conf.nx * conf.ny * conf.nz /
            (compute_s / conf.nsteps);
        fprintf(stderr, "iso_r4_2x:   %8.10f s total, %8.10f s/step, %8.2f "
                "cells/s/step\n", total_s, compute_s / conf.nsteps, point_rate);

        if (conf.save_text != -1) {
            save_layer_text(curr, conf.save_text, dimx, dimy, conf.ny, conf.nx,
                    "snap.text", conf.radius);
        }
    }

    CHECK(cudaFreeHost(curr));
    CHECK(cudaFreeHost(next));
    CHECK(cudaFreeHost(vsq));
    for (int i = 0; i < conf.nsrcs; i++) {
        free(srcs[i]);
    }
    free(srcs);

    for (int g = 0; g < conf.ngpus; g++) {
        CHECK(cudaFree(d_curr[g]));
        CHECK(cudaFree(d_next[g]));
        CHECK(cudaFree(d_vsq[g]));
        CHECK(cudaStreamDestroy(streams[g]));
    }
    free(d_curr);
    free(d_next);
    free(d_vsq);
    free(streams);

    CHECK(cudaFreeHost(low_send_buff));
    CHECK(cudaFreeHost(low_recv_buff));
    CHECK(cudaFreeHost(high_send_buff));
    CHECK(cudaFreeHost(high_recv_buff));

    MPI_Finalize();

    return 0;
}
