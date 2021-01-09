import numpy as np
import math

with open("internal_dists.txt") as fileObj:
    i = 0
    j = 0

    dist_matrix_from_file = np.zeros((164, 164))
    D_matrix = np.zeros((164, 164))

    for row in fileObj.readlines():
        for dist in row.split():
            dist_matrix_from_file[i, j] = dist
            j = j + 1
        j = 0
        i = i + 1

    for n in range(len(dist_matrix_from_file)):
        for m in range(len(dist_matrix_from_file)):
            D_matrix[n, m] = (-(dist_matrix_from_file[n, m] ** 2) + dist_matrix_from_file[n, 0] ** 2 +
                              dist_matrix_from_file[0, m] ** 2) / 2

    u, sigma, v = np.linalg.svd(D_matrix)
    sigma_1d = []
    for s in sigma:
        sigma_1d.append(math.sqrt(s))

    diag_sigma = np.diag(sigma_1d)

    X = np.matmul(u, diag_sigma)
    XT = (np.matmul(v, diag_sigma)).transpose()

    coord_mat = np.matmul(X, XT)

    rank = np.linalg.matrix_rank(D_matrix)
    print(rank)
