from Bio import PDB
import numpy as np
import math


def quaternion_rotation_matrix(Q):
    q0 = Q[0]
    q1 = Q[1]
    q2 = Q[2]
    q3 = Q[3]

    r00 = (q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3)
    r01 = 2 * (q1 * q2 - q0 * q3)
    r02 = 2 * (q1 * q3 + q0 * q2)

    r10 = 2 * (q2 * q1 - q0 * q3)
    r11 = (q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3)
    r12 = 2 * (q2 * q3 - q0 * q1)

    r20 = 2 * (q3 * q1 + q0 * q2)
    r21 = 2 * (q3 * q2 - q0 * q1)
    r22 = (q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3)

    # 3x3 rotation matrix
    rot_matrix = np.array([[r00, r01, r02],
                           [r10, r11, r12],
                           [r20, r21, r22]])

    return rot_matrix


parser = PDB.PDBParser()
io = PDB.PDBIO()
struct = parser.get_structure('1lyd', '1lyd.pdb')

Q = [math.cos(15), 0, math.sin(15) * (2 / math.sqrt(5)), math.sin(15) * (1 / math.sqrt(5))]

rotation_matrix = quaternion_rotation_matrix(Q)

for atom in struct.get_atoms():
    atom_C1 = atom.coord.copy()
    break

for model in struct:
    for chain in model:
        for residue in chain:
            for atom in residue:
                atom.transform(rotation_matrix, -atom_C1)

io.set_structure(struct)
io.save('1lyd_coord.pdb')