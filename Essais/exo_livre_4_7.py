#affiche n etoile a la suite
def etoile(a):
    m='*'
    msg = '*'
    i=1
    while i<a:
        msg=msg +  m
        i=i+1
    return msg



#affiche tous les multiples de table de / mais pas multiple de trois
def main():
    i=1
    while i<=20:
        a=i*7
        if(a%3 ==0):
            a="*"
        print(a)
        i=i+1
        
    print("\n\n")
#piramide d etoile
    a=input("Nombre de steps :")
    i=1
    while i<=a:
        print(etoile(i))
        i=i+1
    
   
    
main()
