class AtomClass:
    def __init__(self,number, atom, x_coord, y_coord, z_coord):
        self.number = number
        self.atom = atom
        self.x_coord = x_coord
        self.y_coord = y_coord
        self.z_coord = z_coord

    def convert_to_string(self):
        print( "Index: {:4} Atom: {:2s}  Coords: ({:6s}, {:6s}, {:6s})\n".format(
            str(self.number), str(self.atom), str(self.x_coord), str(self.y_coord), str(self.z_coord)))
