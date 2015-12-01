#! /home/seb/Documents/Python
#premiers pas
def Pi():
    return 3.14
    
def cube(a):
    return a*a*a
    
def VolumeSphere(r):    
    return 4*Pi()*cube(r)/3



def main():
    liste = [17, 38, 10, 25, 72]  
    
    #trier la liste ou liste.sort()
    sorted(liste)
    #ajout d un element a la fin
    liste.append(12)    
    #pour reverser
    liste.reverse()
    #donner l indie d un element
    a= liste.index(17)
    #le num 17 est a la place
    print "L enplacement du 17 est:",  liste.index(17),  '\n'
    #enlever l element 38 et afficher la liste 
    print liste
    liste.remove(38)  #ou liste.pop([liste.index(38)])
    print liste
    #afficher la sous liste 2e et 3e element !!!!!!!!!!!!!!!ATTENTION comence a zero
    print liste[1]
    print liste[2]
    #afficher sous liste deouis 2e element
    a=len(liste)
    print liste[1:]
    #afficher sous liste du 3e element a la fin de la liste
    print liste[:2]
    #afficher last things avec negative index !!!!!!!!!!!!!!! ce lit depuis la droite
    print liste[-1]
    
main()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
