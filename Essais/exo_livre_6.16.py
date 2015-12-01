__author__ = 'seb'
#force garvitationnelle
#6.16

def F_grav(m1, m2, d):
    return 6.67*10E-11 *m1*m2/(float(d)*float(d))

def main():
    D=[]    #liste de distance
    F=[]
    ch="stqrt"
    while ch :
        ch = input("Entrer la distance: ")
        if ch != 0 and ch != "":
            D.append(float(ch))
        else :
            ch = 0


    m1=m2=1e4
    i=0
    while i<len(D):
        F.append(F_grav(m1, m2, D[i]))
        i=i+1
        print("d = ", D[i], "m", " laforce vaut ", F[i], 'N')

main()