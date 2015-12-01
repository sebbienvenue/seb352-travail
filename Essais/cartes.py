__author__ = 'seb'

from random import *
from tkinter import *

class JeudeCarte(object):
    """Faire un jeu de carte de 52 pour tirer et melanger
    et lire la premiere carte du paquet"""
    coul=['Pique', 'Trèfle', 'Coeur', 'Carreau']
    valeur= [0,0,2,3,4,5,6,7,8,9,10,'valet','dame','roi', 'as']

    def __init__(self):
        #creation du jeu de 52
        self.carte=[]
        for coul in range(4):
            for val in range(13):
                self.carte.append((val+2, coul))

    def nom_carte(self, c):
        #c objet carte
        return "{0} de {1}".format(self.valeur[c[0]], self.coul[c[1]])

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

class Bataille(JeudeCarte):
    """Fais un je de bataille"""
    def __init__(self):
        self.A=JeudeCarte()
        self.B=JeudeCarte()
        self.A.battre()
        self.B.battre()
        self.a, self.b = 0,0 #pts respectifs

    def jeu(self):
        for i in range(52):
            self.ca, self.cb = self.A.tirer(), self.B.tirer()
            self.val_a, self.val_b=self.ca[0], self.cb[0]
            if self.val_a > self.val_b:
                self.a+=1
            elif self.val_a < self.val_b:
                self.b+=1
            print("A: {0} / B: {1}\n{2} / {3}".format(self.A.nom_carte(self.ca), self.B.nom_carte(self.cb), self.a, self.b))

class Tapis(JeudeCarte):
    """Creer un tapis de jeu avec les cartes en GUI"""
    def __init__(self, height=600, width=800, bg= 'green',):
        self.fen=Tk()
        JeudeCarte.__init__(self)
        self.can = Canvas(self.fen, height=height, width=width, bg= bg)
        self.can.grid(row=0, column=0)
        self.can.create_rectangle(30, 30, 800-30, 600-30, width=2, fill=None)
        Button(self.can, text="Quitter", command=quit).grid(row=1, column =0)

        self.fen.mainloop()

    def carte_dos(self, x, y):
        """Creer une carte de dos"""
        #x et y en argument le centre de la carte
        # DIM de la carte 80*100
        # Attention pour eviter un depacementajuster avec les valeurs du Canvas prédéfinies
        a, b, c, d= x-80/2, y-100/2-1, x+80/2, y+100/2-1
        self.can.create_rectangle(a, b, c, d, fill='red', width=1)
        a, c= x-16/2, x+16/2
        self.can.create_rectangle(a, b, c, d, fill='ivory', width=0)

    def bataille(self):
        #lance un jeu de bataille

    #def affiche_pts(self, a, b):
        #affiche les pts dans une fenetre sur le Canvas









################################################ boucle test
if __name__=='__main__':
    J=Bataille()
    #J.jeu()
    Q=Tapis()

