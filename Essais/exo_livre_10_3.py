__author__ = 'seb'
#travail sur des chaines de caracteres

def trouve(ch, ca):
    "Retourne la position du premier caractere dans la chaine, -1 si ne trouve rien"
    n = len(ch)
    i=0
    while i < n:
        if ch[i] == ca:
            return i
        else:
            i += 1
    return -1
#fonction ameliorée qui renvoit a partir d une ligne donnée pour la recherche
def trouve2(ch, ca, indice):
    "Retourne la position du premier caractere dans la chaine, -1 si ne trouve rien"
    n = len(ch)
    i = indice
    while i < n:
        if ch[i] == ca:
            return i
        else:
            i += 1
    return -1

def compteCar(ch, ca):
    "Compte le nombre de caractere ca trouve dans ch"
    n = len(ch)
    i = 0
    ok=0
    while i < n:
        if ch[i] == ca:
            ok += 1
        i += 1
    return ok




#main
print("ncsjncjscsnojcascsacnasocoasc")
print(compteCar("ananas au jus", "a"))
print(compteCar("Gédéon est déjà là", "é"))
