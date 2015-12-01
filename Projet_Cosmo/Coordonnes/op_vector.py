__author__ = 'seb'

from math import*

class Vectors_operations(object):

    def mult_const(self, vector, lamda):
        chaine=[]
        for i in range(0, len(vector)):
            chaine.append(vector[i]*lamda)
        return chaine


    def scalaire_product(self, a,b):
        somme =0
        if len(a)!= len(b):
            print("erreur de dimension entre les deux vecteur\n")
            return
        else:
            for i in range(0, len(a), 1):
                somme += a[i]*b[i]
            return somme

    def somme_vector(self, a, b):
        chaine=[]
        if len(a)!= len(b):
            print("erreur de dimension entre les deux vecteur\n")
            return
        else:
            for i in range(0, len(a), 1):
                chaine.append(a[i]+b[i])
            addition = tuple(chaine)
            return addition

    def soustraction_vector(self, a, b):
        """Soustraction du vecteur a par le vecteur b"""""
        chaine=[]
        if len(a)!= len(b):
            print("erreur de dimension entre les deux vecteur\n")
            return
        else:
            for i in range(0, len(a), 1):
                chaine.append(a[i]-b[i])
                print(i)
                print(a[i]-b[i])
            soustraction = tuple(chaine)
            return soustraction

if __name__=='__main__':
    Vec=Vectors_operations()

    Vec.a, Vec.b= (1, 1, 1), (-1, 2, 3)
    print(Vec.somme_vector(Vec.a, Vec.b))
    print(Vec.soustraction_vector(Vec.a, Vec.b))
    print(Vec.scalaire_product(Vec.a, Vec.b))
    print(Vec.a[2]-Vec.b[2])
