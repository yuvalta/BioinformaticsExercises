import os
from HW1.AtomClass import AtomClass
import numpy as np
import matplotlib.pyplot as plt
import math

protein_coords_folder = "HW1/protein_coords"


def normalize(x):
    length = length_of_vector(x)
    norm_array = [x[0] / length, x[1] / length, x[2] / length]
    return norm_array


def length_of_vector(x):
    return np.sqrt(np.square(x[0]) + np.square(x[1]) + np.square(x[2]))


def degrees(rad_angle):
    angle = rad_angle * (180 / math.pi)
    return angle


def vector_subtract(v1_x, v1_y, v1_z, v2_x, v2_y, v2_z):
    return [v1_x - v2_x, v1_y - v2_y, v1_z - v2_z]


for filename in os.listdir(protein_coords_folder):

    with open(os.path.join(protein_coords_folder, filename)) as protein_file:

        index_list = []

        # create list of all N, CA, C
        for row in protein_file.readlines():
            first_row_split_list = row.split()

            if first_row_split_list[0] == "ATOM":
                protein_name = first_row_split_list[2]
                if protein_name == "CA" or protein_name == "C" or protein_name == "N":
                    first_row = AtomClass(first_row_split_list[5], first_row_split_list[2],
                                          float(first_row_split_list[6]),
                                          float(first_row_split_list[7]), float(first_row_split_list[8]))
                    index_list.append(first_row)
            elif first_row_split_list[0] == "TER":
                break

        angles_dict = {}

        for i in range(0, len(index_list), 3):
            if i + 6 < len(index_list):
                p1_phi = index_list[i + 2]  # C_0
                p2_phi = index_list[i + 3]  # N_1
                p3_phi = index_list[i + 4]  # CA_1
                p4_phi = index_list[i + 5]  # C_1

                r1_phi = vector_subtract(p3_phi.x_coord, p3_phi.y_coord, p3_phi.z_coord, p2_phi.x_coord,
                                         p2_phi.y_coord, p2_phi.z_coord)
                r2_phi = vector_subtract(p1_phi.x_coord, p1_phi.y_coord, p1_phi.z_coord, p2_phi.x_coord,
                                         p2_phi.y_coord, p2_phi.z_coord)
                r3_phi = vector_subtract(p4_phi.x_coord, p4_phi.y_coord, p4_phi.z_coord, p3_phi.x_coord,
                                         p3_phi.y_coord, p3_phi.z_coord)

                n1_phi = np.cross(r1_phi, r2_phi)
                n2_phi = np.cross(r1_phi, r3_phi)

                phi = degrees(math.acos(np.dot(n1_phi, n2_phi) / (length_of_vector(n1_phi) * length_of_vector(n2_phi))))

                arr_pos_d = np.cross(r1_phi, n1_phi)

                if np.dot(n2_phi, arr_pos_d) < 0:  # if dot product is negative - vectors are in opposite directions זווית כהה
                    phi = -phi

                p1_psi = index_list[i + 3]  # N_1
                p2_psi = index_list[i + 4]  # CA_1
                p3_psi = index_list[i + 5]  # C_1
                p4_psi = index_list[i + 6]  # N_2

                r1_psi = vector_subtract(p3_psi.x_coord, p3_psi.y_coord, p3_psi.z_coord, p2_psi.x_coord,
                                         p2_psi.y_coord, p2_psi.z_coord)
                r2_psi = vector_subtract(p1_psi.x_coord, p1_psi.y_coord, p1_psi.z_coord, p2_psi.x_coord,
                                         p2_psi.y_coord, p2_psi.z_coord)
                r3_psi = vector_subtract(p4_psi.x_coord, p4_psi.y_coord, p4_psi.z_coord, p3_psi.x_coord,
                                         p3_psi.y_coord, p3_psi.z_coord)

                n1_psi = np.cross(r1_psi, r2_psi)
                n2_psi = np.cross(r1_psi, r3_psi)

                psi = degrees(math.acos(np.dot(n1_psi, n2_psi) / (length_of_vector(n1_psi) * length_of_vector(n2_psi))))

                arr_pos_d = np.cross(r1_psi, n1_psi)

                if np.dot(n2_psi, arr_pos_d) < 0:  # if dot product is negative - vectors are in opposite directions זווית כהה
                    psi = -psi

                angles_dict[phi] = psi

    plt.scatter(list(angles_dict.keys()), list(angles_dict.values()), s=1)
    plt.title(filename)
    plt.axis([-180, 180, -180, 180])
    plt.xlabel('Phi')
    plt.ylabel('Psi')
    plt.show()
