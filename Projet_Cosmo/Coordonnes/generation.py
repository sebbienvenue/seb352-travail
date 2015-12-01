__author__ = 'seb'

from op_vector import *

class Gen_graphene(Vectors_operations):
    def __init__(self):
        self.lattice_vector =[(3,0), (0, sqrt(3))]
        self.motif = []
        self.lattice = []

    def gen_graphene(self):
        tri = (0.5, sqrt(3)/2)
        self.motif = [(0,0), tri, (tri[0]+1, tri[1]), (2,0)]
        self.Bohr_radius = 5.2917721e-11
        chaine = []
        for i in list(range(-3, 4)):
            for j in list(range(-3, 4)):
                chaine.append(self.somme_vector(self.mult_const(self.lattice_vector[0], i), self.mult_const(self.lattice_vector[1], j)))
        for i in list(range(0, len(chaine))):
            for j in [0, 1, 2, 3]:
                self.lattice.append(self.somme_vector(chaine[i], self.motif[j]))

    def store(self):
        file_xyz= open("coord0.xyz", 'w')
        file_dat= open("coord0.dat", 'w')
        file_xyz.write("{}\n".format(len(self.lattice)))
        file_dat.write("{}\n".format(len(self.lattice)))
        for i in list(range(0, len(self.lattice))):
            #self.lattice[i][2]=0
            file_xyz.write("C {0} {1} {2}\n".format(1.5*self.lattice[i][0], 1.5*self.lattice[i][1], str(0)))
            file_dat.write("C {0} {1} {2}\n".format(1.5*self.lattice[i][0], 1.5*self.lattice[i][1], str(0)))
        file_xyz.close()
        file_dat.close()

    def affiche(self):
        print(list(range(0, 4)))
        print(list(range(-3, 4)))
        print(self.lattice)
        print(len(self.lattice))
        print()

    def change_file_relaxed(self, str_path):
        #976 ////////////////   17
        file_relaxed = open()




if __name__ == '__main__':

    V=Gen_graphene()
    V.gen_graphene()
    V.affiche()
    V.store()
