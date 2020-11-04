import os
from AtomClass import AtomClass
import numpy as np
import matplotlib.pyplot as plt
import math

protein_coords_folder = "./protein_coords"


def normalize(x):
    length_of_vector = np.sqrt(np.square(x[0]) + np.square(x[1]) + np.square(x[2]))
    norm_array = [x[0] / length_of_vector, x[1] / length_of_vector, x[2] / length_of_vector]
    return norm_array


def degrees(rad_angle):
    if rad_angle is None:
        return None
    angle = rad_angle * 180 / math.pi
    while angle > 180:
        angle = angle - 360
    while angle < -180:
        angle = angle + 360
    return angle


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
                    # first_row.convert_to_string()
                    index_list.append(first_row)
            elif first_row_split_list[0] == "TER":
                break

        angles_dict = {}

        for i in range(0, len(index_list), 3):
            if i + 5 < len(index_list):
                p1_phi = index_list[i + 2]  # C_0
                p2_phi = index_list[i]  # N_0
                p3_phi = index_list[i + 1]  # CA_0
                p4_phi = index_list[i + 5]  # C_1

                r1_phi = normalize(np.array([p1_phi.x_coord - p2_phi.x_coord, p1_phi.y_coord - p2_phi.y_coord,
                                             p1_phi.z_coord - p2_phi.z_coord]))

                r2_phi = normalize(np.array([p3_phi.x_coord - p2_phi.x_coord, p3_phi.y_coord - p2_phi.y_coord,
                                             p3_phi.z_coord - p2_phi.z_coord]))

                r3_phi = normalize(np.array([p2_phi.x_coord - p3_phi.x_coord, p2_phi.y_coord - p3_phi.y_coord,
                                             p2_phi.z_coord - p3_phi.z_coord]))

                r4_phi = normalize(np.array([p4_phi.x_coord - p3_phi.x_coord, p4_phi.y_coord - p3_phi.y_coord,
                                             p4_phi.z_coord - p3_phi.z_coord]))

                n1_phi = np.cross(r1_phi, r2_phi)
                n2_phi = np.cross(r3_phi, r4_phi)

                phi = degrees(np.arccos(np.dot(n1_phi, n2_phi)))

                p1_psi = index_list[i]  # N_0
                p2_psi = index_list[i + 1]  # CA_0
                p3_psi = index_list[i + 5]  # C_1
                p4_psi = index_list[i + 3]  # N_2

                r1_psi = normalize(np.array([p1_psi.x_coord - p2_psi.x_coord, p1_psi.y_coord - p2_psi.y_coord,
                                             p1_psi.z_coord - p2_psi.z_coord]))

                r2_psi = normalize(np.array([p3_psi.x_coord - p2_psi.x_coord, p3_psi.y_coord - p2_psi.y_coord,
                                             p3_psi.z_coord - p2_psi.z_coord]))

                r3_psi = normalize(np.array([p2_psi.x_coord - p3_psi.x_coord, p2_psi.y_coord - p3_psi.y_coord,
                                             p2_psi.z_coord - p3_psi.z_coord]))

                r4_psi = normalize(np.array([p4_psi.x_coord - p3_psi.x_coord, p4_psi.y_coord - p3_psi.y_coord,
                                             p4_psi.z_coord - p3_psi.z_coord]))

                n1_psi = np.cross(r1_psi, r2_psi)
                n2_psi = np.cross(r3_psi, r4_psi)

                psi = degrees(np.arccos(np.dot(n1_psi, n2_psi)))

                angles_dict[phi] = psi

    plt.scatter(list(angles_dict.keys()), list(angles_dict.values()), s=5)
    plt.title(filename)
    # plt.axis([-180, 180, -180, 180])
    plt.xlabel('Phi')
    plt.ylabel('Psi')
    plt.show()
