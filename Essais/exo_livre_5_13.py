#generation de liste
#5.11

def listes():
    t1=[31,28,31,30,31,30,31,30,31,30,31]
    t2=['Janvier', 'Fevrier', 'Mars', 'Avril', 'Mai', 'Juin', 'Juillet', 'Aout', 'Septembre', 'Octobre', 'Novembre', 'Decembre']
    a=t1
    b=t2
    return t1,  t2
    
    
def melangeur_listes(t1, t2):
    n1=len(t1)
    n2=len(t2)
    i=0
    t3=''
    
    if n1 != n2 :
        print("Les deux listes n\'ont pas la meme longueuret pas compatibles !!!!!!!!!")
    elif n1-n2 ==1:
        while i<n2:
            t3 = t3 + t1[i] + t2[i]
            i=i+1
            
        t3 =t3 + t1[i]
    elif n1-n2 == -1:
        while i<n1:
            t3 = t3 + t1[i] + t2[i]
            i=i+1
            
        t3 =t3 + t2[i]
    else:
        while i<n2:
            t3 = t3 + t1[i] + t2[i]
            i=i+1
    return t3
    
    
def main():
    t1,  t2 = liste()
    c = melageur_listes(a, b)
    print(c)
    
    
    
main()
