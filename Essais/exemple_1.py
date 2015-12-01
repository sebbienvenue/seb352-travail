#! /home/seb/Documents/Python
#premiers pas
def Pi():
    return 3.14
    
def cube(a):
    return a*a*a
    
def VolumeSphere(r):    
    return 4*Pi()*cube(r)/3

def main():
    print ("Fuck")  
    r=3
    V=VolumeSphere(r)
    print "Voila le volume de la sphere: ",   V
    print(type(r))
    
main()



