import numpy as np

n = 10

with open('mp_mat', 'wb') as f:
    A = np.empty((n,n))
    A.fill(99999.)

    # Make diagonal 0
    for i in range(n):
        A[i][i] = 0.

    # Make off diag 1
    for i in range(n-1):
        A[i+1][i] = 1.
        A[i][i+1] = 1.

    # Write file as binary
    A.astype('float32').tofile(f)

    #print A

'''
with open('mat', 'rb') as f:
    A = np.fromfile(f, dtype=np.float64)

    print A
'''
