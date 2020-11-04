import os

protein_coords_folder = "./protein_coords"

for filename in os.listdir(protein_coords_folder):
    print(os.path.join(protein_coords_folder, filename))