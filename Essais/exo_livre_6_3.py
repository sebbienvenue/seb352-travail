#Calculer les dimension d un triangle
#6.2

def Perimetre_tri(a, b, c, d):
#a , b, c triangle
#d longueur du demi perimetre
#retourne le perimetre
    return a+b+c

def Aire_tri(a, b, c, d):
#a , b, c triangle
#d longueur du demi perimetre
    from math import sqrt
    return sqrt(d*(d-a)*(d-b)*(d-c))

def freq_pedulus(L):
#retourne la frequence d un pendule simple 
    from math import sqrt
    from math import pi
    return 2*pi*sqrt(L/9.81)

