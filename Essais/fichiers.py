__author__ = 'seb'

def existe(file_name):
    #script d une fonction pour eviter tout plantage si n existe pas
    "Retourne 1 si fichier existe et 0 si pas"
    try:
        f = open(file_name, 'r')
        f.close()
        return 1
    except:
        return 0

#pickle pour les fichiers binaires
