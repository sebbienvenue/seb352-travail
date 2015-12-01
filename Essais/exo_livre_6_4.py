#cree une liste que l on change dans une boucle
#6.4
def message():
    return "Que voulez vous faire avec votre liste ?\nTapez 1 pour ajouter un element a la liste.\nTapez 2 pour afficher la liste.\nTapez 3 pour arreter\n"

def g_liste(table):
    ok = 1
    while ok == 1:
        a = input(message())
        print(a)
        if a == 1:
            b=input("Saisir ce qui doit etre ajoute: ")
            table.append(b)
        elif a == 2:
            print(t)
        elif a == 3:
            print("Donner rentree non recevable recommencer !!")
        else:
            ok=0




def main():
    t=[]
    ch = "start"
    while ch != "":
        ch=input("Veuiller entrer une valeur : ")
        if ch != "":
            t.append(float(ch))

    print(t)


main()
