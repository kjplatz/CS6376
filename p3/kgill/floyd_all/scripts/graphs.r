library(gridExtra)

omp_all = read.csv("../data/floyd_omp_8192_all.txt")
omp_16 = read.csv("../data/floyd_omp_28_16384.txt")
oacc = read.csv("../data/floyd_acc_28_16384.txt")
oacc_mpi = read.csv("../data/floyd_acc_mpi_p100_4_8192.txt")
mpi = read.csv("../data/floyd_mpi_8_8192.txt")
omp_mpi = read.csv("../data/floyd_omp_mpi_all.txt")

# Convert FLOPS to GFLOPS
omp_all$flops = omp_all$flops /1000000000
omp_16$flops = omp_16$flops/1000000000
oacc$flops = oacc$flops/1000000000
oacc_mpi$flops = oacc_mpi$flops/1000000000
mpi$flops = mpi$flops/1000000000
omp_mpi$flops = omp_mpi$flops/1000000000

# Delete cores from oacc
oacc = within(oacc, rm(cores))
oacc_mpi = within(oacc_mpi, rm(cores))

# Aggregate the cores and mean flops
omp_all_ag = aggregate(flops~cores, omp_all, mean)
omp_16_ag = aggregate(flops~cores, omp_16, mean)
omp_mpi_ag = aggregate(flops~nodes, omp_mpi, mean)

# Add the speedup
omp_all_ag$speedup = omp_all_ag$flops/omp_all_ag$flops[1]

# Add the serial portion
omp_all_ag$serial = (omp_all_ag$speedup - omp_all_ag$cores)/(1 - omp_all_ag$cores)

# Add the efficiency 
omp_all_ag$karp = (((1/omp_all_ag$speedup)-(1/omp_all_ag$cores))/(1 - (1/omp_all_ag$cores)))

# Create plots for the omp_all flops
png("omp_all.png", width=600, height=500)
plot(flops~cores, omp_all_ag[1:18,], main="Cores vs FLOPS (10 iterations)", xlab="Cores", 
	ylab="Giga FLOPS", col="red", pch=16)
dev.off()

# Create plots for the omp_all flops
png("omp_mpi_box.png", width=600, height=500)
plot(flops~nodes, omp_mpi, main="Nodes vs FLOPS (10 iterations)", xlab="Nodes", 
	ylab="Giga FLOPS", col="green", pch=8)
dev.off()


# Create plot for the karp-flatt metic
png("karp_flatt_graph.png")
plot(karp~cores, omp_all_ag[1:18,], col="blue", pch=16, main="Karp-Flatt metric", xlab="Cores")
dev.off()

png("oacc.png", width=400, height=400) 
grid.table(oacc, rows=NULL)
dev.off()

png("oacc_mpi.png", width=400, height=400) 
grid.table(oacc_mpi, rows=NULL)
dev.off()

png("mpi.png", width=400, height=400) 
grid.table(mpi, rows=NULL)
dev.off()

png("omp_mpi.png", width=400, height=400) 
grid.table(omp_mpi_ag, rows=NULL)
dev.off()





