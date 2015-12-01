__author__ = 'seb'
#ajouter des coéplement a Newton et faire une animation
#non fini
from tkinter import *
from math import *
from random import randrange

def F_g(m1, m2, d):
    return 6.67E-11 * m1 * m2 / (d*d)

def rand_coul():
    c=["red", "blue", "light blue", "pink", "purple", "yellow", "green", "grey"]
    return c[randrange(0, len(c)-1, step =1)]

def sphere(x, y, r):
    "Creer un objet planete de rayon r en x, y"
    return can.create_oval(x, y, x+r, y+r, width=0, fill=rand_coul())

def bouge(objet, X, Y, x, y):
    "Deplace un objet des coordones indiquees"
    #X, Y viennent de parametres indiquqnts la position des objets
    can.coords(objet, X, Y, X+x, Y+y)

def haut1():
    global x1, y1
    bouge(p1, x1, y1, 0, -10)
def haut2():
    global x2, y2
    bouge(p2, x2, y2, 0, -10)
def bas1():
    global x1, y1
    bouge(p1, x1, y1, 0, 10)
def bas2():
    global x2, y2
    bouge(p2, x2, y2, 0, 10)
def gauche1():
    global x1, y1
    bouge(p1, x1, y1, -10, 0)
def gauche2():
    global x2, y2
    bouge(p2, x2, y2, -10, 0)
def droite1():
    global x1, y1
    bouge(p1, x1, y1, 10, 0)
def droite2():
    global x2, y2
    bouge(p2, x2, y2, 10, 0)


#initialisation
m1, m2 , d = 1E10, 2E5, 1000000
x1, y1, x2, y2, r1, r2 = 100, 120, 400, 350, 80, 120

fen = Tk()
fen.title("Newton et ses planetes")

txt1=Label(fen , text="Planete 1")
txt2=Label(fen , text="Planete 2")
tex3=Label(fen, text= "M1 = " + str(m1) + " kg")
tex4=Label(fen, text= "M2 = " + str(m2) + " kg")
text5=Label(fen, text="Distance" + str(distance))
text6=Label(fen, text="Force")
text3.grid(row= 1, column=0)
text4.grid(row= 1, column=1)
can = Canvas(fen , height=600, width = 750, bg='ivory')
can.grid(row=2, column=0, columnspan=2)
text5.grid(row=3, column=0)
text6.grid(row=3, column=1)

#creation des deux objets qui sont des sphere
p1 = sphere(x1, y1, r1)
p2 = sphere(x2, y2, r2)

#planete 1
f1=Frame(fen)
f1.grid(row=, column=0, sticky= W, padx = 10)
bou1= Button(f1, text="^", command=haut1).pack(side = LEFT)
bou2= Button(f1, text="v", command=bas1).pack(side = LEFT)
bou3= Button(f1, text="<-", command=gauche1).pack(side = LEFT)
bou4= Button(f1, text="->", command=droite1).pack(side = LEFT)

#planete 2
f2=Frame(fen)
f2.grid(row=, column=1, sticky= E, padx = 10)
bou5= Button(f2, text="^", command=droite2).pack(side = RIGHT)
bou6= Button(f2, text="<-", command=gauche2).pack(side = RIGHT)
bou7= Button(f2, text="v", command=bas2).pack(side = RIGHT)
bou8= Button(f2, text="^", command=haut2).pack(side = RIGHT)


fen.mainloop()