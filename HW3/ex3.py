import numpy as np
from scipy.stats import rankdata
import math
from Bio import PDB

with open("internal_dists.txt") as fileObj:
    i = 0
    j = 0

    dist_matrix_from_file = np.zeros((164, 164))
    dist_matrix_square = np.zeros((164, 164))
    D_matrix = np.zeros((164, 164))
    coord_mat = np.zeros((164, 3))
    sigma_1d = []

    for row in fileObj.readlines():
        for dist in row.split():
            dist_matrix_from_file[i, j] = float(dist)
            j = j + 1
        j = 0
        i = i + 1

    for row in range(len(dist_matrix_square)):
        for col in range(len(dist_matrix_square)):
            dist_matrix_square[row, col] = pow(
                dist_matrix_from_file[row, 0] - dist_matrix_from_file[0, col],
                2)

    for n in range(len(dist_matrix_from_file)):
        for m in range(len(dist_matrix_from_file)):
            D_matrix[n, m] = (
                                     (dist_matrix_square[n, m]) ** 2 - dist_matrix_square[n, 0] ** 2 -
                                     dist_matrix_square[0, m] ** 2) / -2

    u, sigma, v = np.linalg.svd(D_matrix)

    # sqrt of sigma
    for s in sigma:
        sigma_1d.append(math.sqrt(s))

    X = np.matmul(u, np.diag(sigma_1d))
    XT = (np.matmul(v, np.diag(sigma_1d))).transpose()

    for row in range(len(coord_mat)):
        for col in range(3):
            coord_mat[row, col] = "{:.5f}".format((X[row, col]))

    parser = PDB.PDBParser()
    io = PDB.PDBIO()
    struct = parser.get_structure("2lyz", "2lyz.pdb")  # for reference - creating pdb file

    for atom in struct.get_atoms():
        atom_C1 = atom.coord.copy()
        break

    i = 0
    for model in struct:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if (i < 164):
                        atom.set_coord(coord_mat[i])
                    else:
                        atom.set_coord([0, 0, 0])
                    i = i + 1
    io.set_structure(struct)
    io.save("dist_mat_2_coord.pdb")  # need to delete all atoms with (0,0,0) coords (row 164 and above)
