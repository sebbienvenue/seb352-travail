#travailler sur une chaine de caracteres 
#5.8

#fonction qui valide le type de chaine 
def validation(msg):
    ok=0
    A=q
    while ok == 0:
        if type(msg) != type(A):
            print("Erreur dans le format de la chaine, il faut des caracteres !!!!")
            msg=input("Entrer une chaine de caractere: ")
        else:
            ok=1

def main():
    msg=input("Entrer une chaine de caractere: ")
    valiadation(msg)
    etoile='*'
    i=0
    long_msg=""
    
    while i<=len(msg):
        long_msg[2*i]=msg[i]
        long_msg[2*i+1]=etoile
   
    print(long_message)
    
    
main()
