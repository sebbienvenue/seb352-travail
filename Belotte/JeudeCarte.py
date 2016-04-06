__author__ = 'seb'

from random import *
from tkinter import *

class JeudeCarte(object):
    """Faire un jeu de carte de 52 pour tirer et melanger
    et lire la premiere carte du paquet"""
    coul=['Pique', 'Trèfle', 'Coeur', 'Carreau']
    nn= ['7', '8', '9', '10','valet','dame','roi', 'as']
    valeur = [0,0,0, 10, 2, 3, 4, 11, 14, 20]

    def __init__(self):
        #creation du jeu de 52
        self.carte=[]
        for coul in range(4):
            for val in valeur[0:8]:
                self.carte.append((val, coul))

    def nom_carte(self, c):
        #c objet carte
        return "{0} de {1}".format(self.nn[c[0]], self.coul[c[1]])

    def battre(self):
        #bat les carte en fesant 2 fois des permutation pour autant dee carte dna sle paquet
        t = len(self.carte)
        for i in range(t):
            h1, h2 = randrange(t), randrange(t)
            self.carte[h1], self.carte[h2] = self.carte[h2], self.carte[h1]

    def tirer(self):
        #tire la premiere carte du paquet
        t = len(self.carte)
        if t > 0:
            carte = self.carte[0]
            del(self.carte[0])
            return carte
        else:
            return None
