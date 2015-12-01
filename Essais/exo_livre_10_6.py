__author__ = 'seb'
#creer une chaine de caracteres

def n_mots(phrase):
    "Revoit le nombre de mots dans la phrase"
    #exo 10_7
    ok = 0
    for e in phrase:
        if e == " ":
            ok += 1
        elif e == ".":
            return ok
    return ok

def tale():
    #10_6
    for e in 'JKLMNOP':
        print(e+"ack")

def estUnChiffre(a):
    #10_9
    "Indique si a est un chiffre, retourne 1 si vrai, 0 sinon"
    if a in '1234567890':
        return 1
    else:
        return 0

def estMaj(a):
    "Indique si la lettre a est une majuscule, retourne 1 si vrai, 0 sinon"
    MAJ="QWERTZUIOPLKJHGFDSAYXCVBNM"
    if a in MAJ:
        return 1
    else:
        return  0

def  chain_list(ch):
    "converti une phrase en liste de mots"
    #on ajoute une chaine temporaire que
    #l on reinitialise une fois ajoutée dans la liste
    #utilisation de for in
    #fas les lettres unes par unes
    L, ct = [], ""

    for c in ch:
        if c == " " or c == ".":
            L.append(ct)
            ct = ""
        else:
            ct += c
    return L



#main
msg = "ijdw odhsd aoid hsdahs dsaodio dihs aspdo das q."
L = chain_list(msg)

for e in L:
    print(e)
