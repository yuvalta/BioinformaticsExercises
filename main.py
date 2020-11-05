import os
from AtomClass import AtomClass
import numpy as np
import matplotlib.pyplot as plt
import math

protein_coords_folder = "./protein_coords"


def normalize(x):
    length = length_of_vector(x)
    norm_array = [x[0] / length, x[1] / length, x[2] / length]
    return norm_array


def length_of_vector(x):
    return np.sqrt(np.square(x[0]) + np.square(x[1]) + np.square(x[2]))


def degrees(rad_angle):
    angle = rad_angle * (180 / math.pi)
    return angle


####

def VectorSubtract(v1_x, v1_y, v1_z, v2_x, v2_y, v2_z):
    # Subtract V2 from V1
    return [v1_x - v2_x, v1_y - v2_y, v1_z - v2_z]


def CrossProduct(v1_x, v1_y, v1_z, v2_x, v2_y, v2_z):
    return [v1_y * v2_z - v1_z * v2_y, v1_z * v2_x - v1_x * v2_z, v1_x * v2_y - v1_y * v2_x]


def Distance3d(x1, y1, z1, x2, y2, z2):
    dx = (x2 - x1);
    dy = (y2 - y1);
    dz = (z2 - z1);
    tmp = dx * dx + dy * dy + dz * dz;
    return math.sqrt(tmp);


def DotProduct3d(v1_x, v1_y, v1_z, v2_x, v2_y, v2_z):
    # Return Dot product of two vectors
    return v1_x * v2_x + v1_y * v2_y + v1_z * v2_z;


def CosAfromDotProduct(aDP, aLine1Len, aLine2Len):
    # Return CosA of Angle between lines aLine1 and aLine2 originating from 0
    # of Lengths aLine1Len and aLine2Len respectively
    return aDP / (aLine1Len * aLine2Len);


def RadToDeg(Radians):
    return Radians * (180 / math.pi);


def Angle(a0_x, a0_y, a0_z, a2_x, a2_y, a2_z, a3_x, a3_y, a3_z):
    # Three Atoms: a1 is Center atom (picked 1st), a2 and a3 are subsequently picked atoms
    # Lines are A and B which join at 0
    # dx, dy, dz;   # Difference in position (to get vector values)
    # lenA, lenB;   # Lengths
    # DProd, CosA;
    lenA = Distance3d(a0_x, a0_y, a0_z, a2_x, a2_y, a2_z);
    lenB = Distance3d(a0_x, a0_y, a0_z, a3_x, a3_y, a3_z);
    dx = a2_x - a0_x;
    dy = a2_y - a0_y;
    dz = a2_z - a0_z;
    vA_x = dx;
    vA_y = dy;
    vA_z = dz;
    dx = a3_x - a0_x;
    dy = a3_y - a0_y;
    dz = a3_z - a0_z;
    vB_x = dx;
    vB_y = dy;
    vB_z = dz;
    DProd = DotProduct3d(vA_x, vA_y, vA_z, vB_x, vB_y, vB_z);
    CosA = CosAfromDotProduct(DProd, lenA, lenB);
    return RadToDeg(math.acos(CosA));


####


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
            if i + 6 < len(index_list):
                p1_phi = index_list[i + 2]  # C_0
                p2_phi = index_list[i + 3]  # N_1
                p3_phi = index_list[i + 4]  # CA_1
                p4_phi = index_list[i + 5]  # C_1


                # arr_dist21 = VectorSubtract(p3_phi.x_coord, p3_phi.y_coord, p3_phi.z_coord, p2_phi.x_coord,
                #                             p2_phi.y_coord, p2_phi.z_coord);
                # arr_dist01 = VectorSubtract(p1_phi.x_coord, p1_phi.y_coord, p1_phi.z_coord, p2_phi.x_coord,
                #                             p2_phi.y_coord, p2_phi.z_coord);
                # arr_dist32 = VectorSubtract(p4_phi.x_coord, p4_phi.y_coord, p4_phi.z_coord, p3_phi.x_coord,
                #                             p3_phi.y_coord, p3_phi.z_coord);
                # dist21_x = arr_dist21[0];
                # dist21_y = arr_dist21[1];
                # dist21_z = arr_dist21[2];
                # dist01_x = arr_dist01[0];
                # dist01_y = arr_dist01[1];
                # dist01_z = arr_dist01[2];
                # dist32_x = arr_dist32[0];
                # dist32_y = arr_dist32[1];
                # dist32_z = arr_dist32[2];
                # arr_dd1 = CrossProduct(dist21_x, dist21_y, dist21_z, dist01_x, dist01_y, dist01_z);
                # arr_dd3 = CrossProduct(dist21_x, dist21_y, dist21_z, dist32_x, dist32_y, dist32_z);
                # dd1_x = arr_dd1[0];
                # dd1_y = arr_dd1[1];
                # dd1_z = arr_dd1[2];
                # dd3_x = arr_dd3[0];
                # dd3_y = arr_dd3[1];
                # dd3_z = arr_dd3[2];
                # r = Angle(0, 0, 0, dd1_x, dd1_y, dd1_z, dd3_x, dd3_y, dd3_z);
                # arr_pos_d = CrossProduct(dist21_x, dist21_y, dist21_z, dd1_x, dd1_y, dd1_z);
                # pos_d_x = arr_pos_d[0];
                # pos_d_y = arr_pos_d[1];
                # pos_d_z = arr_pos_d[2];
                # if (DotProduct3d(dd3_x, dd3_y, dd3_z, pos_d_x, pos_d_y, pos_d_z) < 0):
                #     r = -r;
                #
                # phi = r

                r1_phi = (np.array([p1_phi.x_coord - p2_phi.x_coord, p1_phi.y_coord - p2_phi.y_coord,
                                    p1_phi.z_coord - p2_phi.z_coord]))

                r2_phi = (np.array([p3_phi.x_coord - p2_phi.x_coord, p3_phi.y_coord - p2_phi.y_coord,
                                    p3_phi.z_coord - p2_phi.z_coord]))

                # r3_phi = (np.array([p2_phi.x_coord - p3_phi.x_coord, p2_phi.y_coord - p3_phi.y_coord,
                #                     p2_phi.z_coord - p3_phi.z_coord]))

                r4_phi = (np.array([p4_phi.x_coord - p3_phi.x_coord, p4_phi.y_coord - p3_phi.y_coord,
                                    p4_phi.z_coord - p3_phi.z_coord]))

                n1_phi = np.cross(r1_phi, r2_phi)
                n2_phi = np.cross(r2_phi, r4_phi)

                phi = degrees(math.acos(np.dot(n1_phi, n2_phi) / (length_of_vector(n1_phi) * length_of_vector(n2_phi))))
                #
                p1_psi = index_list[i + 3]  # N_1
                p2_psi = index_list[i + 4]  # CA_1
                p3_psi = index_list[i + 5]  # C_1
                p4_psi = index_list[i + 6]  # N_2


                # arr_dist21 = VectorSubtract(p3_psi.x_coord, p3_psi.y_coord, p3_psi.z_coord, p2_psi.x_coord,
                #                             p2_psi.y_coord, p2_psi.z_coord)
                # arr_dist01 = VectorSubtract(p1_psi.x_coord, p1_psi.y_coord, p1_psi.z_coord, p2_psi.x_coord,
                #                             p2_psi.y_coord, p2_psi.z_coord)
                # arr_dist32 = VectorSubtract(p4_psi.x_coord, p4_psi.y_coord, p4_psi.z_coord, p3_psi.x_coord,
                #                             p3_psi.y_coord, p3_psi.z_coord)
                # dist21_x = arr_dist21[0];
                # dist21_y = arr_dist21[1];
                # dist21_z = arr_dist21[2];
                # dist01_x = arr_dist01[0];
                # dist01_y = arr_dist01[1];
                # dist01_z = arr_dist01[2];
                # dist32_x = arr_dist32[0];
                # dist32_y = arr_dist32[1];
                # dist32_z = arr_dist32[2];
                # arr_dd1 = CrossProduct(dist21_x, dist21_y, dist21_z, dist01_x, dist01_y, dist01_z);
                # arr_dd3 = CrossProduct(dist21_x, dist21_y, dist21_z, dist32_x, dist32_y, dist32_z);
                # dd1_x = arr_dd1[0];
                # dd1_y = arr_dd1[1];
                # dd1_z = arr_dd1[2];
                # dd3_x = arr_dd3[0];
                # dd3_y = arr_dd3[1];
                # dd3_z = arr_dd3[2];
                # r = Angle(0, 0, 0, dd1_x, dd1_y, dd1_z, dd3_x, dd3_y, dd3_z);
                # arr_pos_d = CrossProduct(dist21_x, dist21_y, dist21_z, dd1_x, dd1_y, dd1_z);
                # pos_d_x = arr_pos_d[0];
                # pos_d_y = arr_pos_d[1];
                # pos_d_z = arr_pos_d[2];
                # if (DotProduct3d(dd3_x, dd3_y, dd3_z, pos_d_x, pos_d_y, pos_d_z) < 0):
                #     r = -r;
                #
                # psi = r


                r1_psi = (np.array([p1_psi.x_coord - p2_psi.x_coord, p1_psi.y_coord - p2_psi.y_coord,
                                    p1_psi.z_coord - p2_psi.z_coord]))

                r2_psi = (np.array([p3_psi.x_coord - p2_psi.x_coord, p3_psi.y_coord - p2_psi.y_coord,
                                    p3_psi.z_coord - p2_psi.z_coord]))

                # r3_psi = (np.array([p2_psi.x_coord - p3_psi.x_coord, p2_psi.y_coord - p3_psi.y_coord,
                #                     p2_psi.z_coord - p3_psi.z_coord]))

                r4_psi = (np.array([p4_psi.x_coord - p3_psi.x_coord, p4_psi.y_coord - p3_psi.y_coord,
                                    p4_psi.z_coord - p3_psi.z_coord]))

                n1_psi = np.cross(r1_psi, r2_psi)
                n2_psi = np.cross(r2_psi, r4_psi)

                psi = degrees(np.arccos(np.dot(n1_psi, n2_psi) / (length_of_vector(n1_psi) * length_of_vector(n2_psi))))

                angles_dict[phi] = psi

    plt.scatter(list(angles_dict.keys()), list(angles_dict.values()), s=5)
    plt.title(filename)
    # plt.axis([-180, 180, -180, 180])
    plt.xlabel('Phi')
    plt.ylabel('Psi')
    plt.show()
