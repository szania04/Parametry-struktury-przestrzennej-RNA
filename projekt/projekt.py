from __future__ import (absolute_import, division, print_function, unicode_literals)
import csv
from Bio.PDB import calc_dihedral, Vector, PDBList, PDBParser
from os import remove, getcwd

#funkcja pobiera plik pdb o nazwie podanej jako argiment
from os.path import exists


def download_file(name):
    pdbl = PDBList()
    pdbl.retrieve_pdb_file(name.upper(), pdir = getcwd(), file_format='pdb')

#funkcja tworzy gotowy plik csv z wartosciami katow
def csv_file():
    fieldnames = ['NR', 'ALFA', 'BETA', 'GAMMA', 'DELTA', 'EPSILON', 'ZETA']
    with open('csv_file_angles.csv', 'w') as file:
        angles = open('angles.txt', 'r').read().split()
        for x in range(0, int(len(angles)/3)):
            file.write('%s;' % angles[x])
            if angles[x+1] == 'ALFA':
                file.write('%s;' % angles[x+1])
            else: file.write('%s;' % '0')
            if angles[x+2] == 'BETA':
                file.write('%s;' % angles[x+2])
            else: file.write('%s;' % '0')
            if angles[x+3] == 'GAMMA':
                file.write('%s;' % angles[x+3])
            else: file.write('%s;' % '0')
            if angles[x+4] == 'DELTA':
                file.write('%s;' % angles[x+4])
            else: file.write('%s;' % '0')
            if angles[x+5] == 'EPSILON':
                file.write('%s;' % angles[x+5])
            else: file.write('%s;' % '0')
            if angles[x+6] == 'ZETA':
                file.write('%s;' % angles[x+6])
            else: file.write('%s;' % '0')
            file.write('\n')

#funkcja obrabia plik .pdb i przygotowywuje w odpowiedniej formie do znajdywania i liczenia katow torsyjnych
def prepere_atom_file(file):
    param = open(file, 'r').read().split('\n')
    atom = open('atom.txt', 'w')
    for x in range(0 , len(param)):
        #znalezienie jedynie sekcji atomow i przepisanie jedynie czesc zawierajacymi nukleotydy
        if param[x][:4] == 'ATOM' and \
                (param[x][18:20] == ' A' or param[x][18:20] == ' C' or param[x][18:20] == ' T' or
                         param[x][18:20] == ' G' or param[x][18:20] == ' U'):
            #wpisanie do pliku wybranych koumn
            atom.write(
                param[x][6:11] + ' ' + param[x][12:16] + ' ' + param[x][17:20] + ' ' + param[x][21:22] + ' '
                + param[x][22:26] + ' ' +
                param[x][30:38] + ' ' + param[x][38:46] + ' ' + param[x][46:54] + '\n')

#funkcja usuwa utorzone w czsaie dzialania programu pliki
def remove_temp_files(infile):
    if exists('pdb' + infile + '.ent'): #czy plik istnieje (czy mozna go usunac)
        remove('pdb' + infile + '.ent') #usuniecie pliku
    if exists('atom.txt'):
        remove('atom.txt')
    if exists('angles.txt'):
        remove('angles.txt')

#funkcaja przesiewajaca wczesniej utworzony plik i znajdujaca katy oraz zapisujaca je do nowego pliku
def all_angles():
    atom = open('atom.txt', 'r')
    angles = open('angles.txt', 'w')
    parameters = atom.read().split() #lista utworzona z zawartosci przesianego pliku pdb
    end = int(len(parameters)/8) #zakres w jakim bedziemy sie poruszac (ilosc lini w orginalnym pliku pdb po przesianiu)
    for x in range(0, end):
        # zmienne logiczne ktore odpowiednio oznaczaja prawde gdy w danym momencie 'poszukujemy' danego kata torsyjnego
        # find pozwala przejsc do kolejnej petli gdy znajdziemy dany atom stanowiacy czesc kata torsyjnego
        # nr nukleotydu w ktorym aktualnie sie znajdujemy
        alfa, beta, gamma, delta, epsilon, zeta, find, nr = False, False, False, False, False, False, False, 0
        if parameters[x * 8 + 1] == 'O3\'':
            alfa = True
            find = True
            nr = 1
        elif parameters[x * 8 + 1] == 'P':
            beta = True
            find = True
        elif parameters[x * 8 + 1] == 'O5\'':
            gamma = True
            find = True
        elif parameters[x * 8 + 1] == 'C5\'':
            delta = True
            find = True
        elif parameters[x * 8 + 1] == 'C4\'':
            epsilon = True
            find = True
        elif parameters[x * 8 + 1] == 'C3\'':
            zeta = True
            find = True
        if find:
            nr = nr + int(parameters[x * 8 + 1 + 3])
            for y in range(x, end):
                #zmiana find na false, bo znow zostanie zminiona na true jak uda sie znalezc kolejny atom
                find = False
                if alfa == True and parameters[y * 8 + 1] == 'P':
                    if int(parameters[y * 8 + 1 + 3]) != nr:
                        break
                    else:
                        find = True
                elif beta == True and parameters[y * 8 + 1] == 'O5\'':
                    #jako ze wiemy ze kat torsyjny beta jest tworzony w obrebie jednego nr nukl to jak wartosc ta jest inna to znaczy ze juz nie znajdziemy go
                    if int(parameters[y * 8 + 1 + 3]) != nr:
                        break
                    else:
                        find = True
                elif gamma == True and parameters[y * 8 + 1] == 'C5\'':
                    if int(parameters[y * 8 + 1 + 3]) != nr:
                        break
                    else:
                        find = True
                elif delta == True and parameters[y * 8 + 1] == 'O4\'':
                    if int(parameters[y * 8 + 1 + 3]) != nr:
                        break
                    else:
                        find = True
                elif epsilon == True and parameters[y * 8 + 1] == 'C3':
                    if int(parameters[y * 8 + 1 + 3]) != nr:
                        break
                    else:
                        find = True
                elif zeta == True and parameters[y * 8 + 1] == 'O3\'':
                    if int(parameters[y * 8 + 1 + 3]) != nr:
                        break
                    else:
                        find = True
                if find:
                    for z in range(y, end):
                        find = False
                        if alfa == True and parameters[z * 8 + 1] == 'O5\'':
                            if int(parameters[z * 8 + 1 + 3]) != nr:
                                break
                            else:
                                find = True
                        elif beta == True and parameters[z * 8 + 1] == 'C5\'':
                            if int(parameters[z * 8 + 1 + 3]) != nr:
                                break
                            else:
                                find = True
                        elif gamma == True and parameters[z * 8 + 1] == 'C4':
                            if int(parameters[z * 8 + 1 + 3]) != nr:
                                break
                            else:
                                find = True
                        elif delta == True and parameters[z * 8 + 1] == 'C3\'':
                            if int(parameters[z * 8 + 1 + 3]) != nr:
                                break
                            else:
                                find = True
                        elif epsilon == True and parameters[z * 8 + 1] == 'O3\'':
                            if int(parameters[z * 8 + 1 + 3]) != nr:
                                break
                            else:
                                find = True
                        elif zeta == True and parameters[z * 8 + 1] == 'P':
                            if int(parameters[z * 8 + 1 + 3]) != nr + 1:
                                break
                            else:
                                find = True
                        if find:
                            for v in range(z, end):
                                if alfa == True and parameters[v * 8 + 1] == 'C5\'':
                                    if int(parameters[v * 8 + 1 + 3]) != nr:
                                        break
                                    else:
                                        # gdy ostatni atom jest znaleziony to zczytywane sa wektory kazdego jednego atomu
                                        # a nastepnie obliczany jest kat torsyjny
                                        #nr, nazwa kata oraz wartosc jego zostaja zapisane do pliku
                                        vector1 = Vector(parameters[x * 8 + 5], parameters[x * 8 + 6], parameters[x * 8 + 7])
                                        vector2 = Vector(parameters[y * 8 + 5], parameters[y * 8 + 6], parameters[y * 8 + 7])
                                        vector3 = Vector(parameters[z * 8 + 5], parameters[z * 8 + 6], parameters[z * 8 + 7])
                                        vector4 = Vector(parameters[v * 8 + 5], parameters[v * 8 + 6], parameters[v * 8 + 7])
                                        angle = calc_dihedral(vector1, vector2, vector3, vector4)
                                        angles.write(str(nr) + ' ALPHA: ' + str(angle) + '\n')
                                elif beta == True and parameters[v * 8 + 1] == 'C4\'':
                                    if int(parameters[v * 8 + 1 + 3]) != nr:
                                        break
                                    else:
                                        vector1 = Vector(parameters[x * 8 + 5], parameters[x * 8 + 6], parameters[x * 8 + 7])
                                        vector2 = Vector(parameters[y * 8 + 5], parameters[y * 8 + 6], parameters[y * 8 + 7])
                                        vector3 = Vector(parameters[z * 8 + 5], parameters[z * 8 + 6], parameters[z * 8 + 7])
                                        vector4 = Vector(parameters[v * 8 + 5], parameters[v * 8 + 6], parameters[v * 8 + 7])
                                        angle = calc_dihedral(vector1, vector2, vector3, vector4)
                                        angles.write(str(nr) + ' BETA: ' + str(angle) + '\n')
                                elif gamma == True and parameters[v * 8 + 1] == 'C3\'':
                                    if int(parameters[v * 8 + 1 + 3]) != nr:
                                        break
                                    else:
                                        vector1 = Vector(parameters[x * 8 + 5], parameters[x * 8 + 6], parameters[x * 8 + 7])
                                        vector2 = Vector(parameters[y * 8 + 5], parameters[y * 8 + 6], parameters[y * 8 + 7])
                                        vector3 = Vector(parameters[z * 8 + 5], parameters[z * 8 + 6], parameters[z * 8 + 7])
                                        vector4 = Vector(parameters[v * 8 + 5], parameters[v * 8 + 6], parameters[v * 8 + 7])
                                        angle = calc_dihedral(vector1, vector2, vector3, vector4)
                                        angles.write(str(nr) + ' GAMMA: ' + str(angle) + '\n')
                                elif delta == True and parameters[v * 8 + 1] == 'O3\'':
                                    if int(parameters[v * 8 + 1 + 3]) != nr:
                                        break
                                    else:
                                        vector1 = Vector(parameters[x * 8 + 5], parameters[x * 8 + 6], parameters[x * 8 + 7])
                                        vector2 = Vector(parameters[y * 8 + 5], parameters[y * 8 + 6], parameters[y * 8 + 7])
                                        vector3 = Vector(parameters[z * 8 + 5], parameters[z * 8 + 6], parameters[z * 8 + 7])
                                        vector4 = Vector(parameters[v * 8 + 5], parameters[v * 8 + 6], parameters[v * 8 + 7])
                                        angle = calc_dihedral(vector1, vector2, vector3, vector4)
                                        angles.write(str(nr) + ' DELTA: ' + str(angle) + '\n')
                                elif epsilon == True and parameters[v * 8 + 1] == 'P':
                                    if int(parameters[v * 8 + 1 + 3]) != nr + 1:
                                        break
                                    else:
                                        vector1 = Vector(parameters[x * 8 + 5], parameters[x * 8 + 6], parameters[x * 8 + 7])
                                        vector2 = Vector(parameters[y * 8 + 5], parameters[y * 8 + 6], parameters[y * 8 + 7])
                                        vector3 = Vector(parameters[z * 8 + 5], parameters[z * 8 + 6], parameters[z * 8 + 7])
                                        vector4 = Vector(parameters[v * 8 + 5], parameters[v * 8 + 6], parameters[v * 8 + 7])
                                        angle = calc_dihedral(vector1, vector2, vector3, vector4)
                                        angles.write(' EPSILON: ' + str(angle) + '\n')
                                elif zeta == True and parameters[v * 8 + 1] == 'O5\'':
                                    if int(parameters[v * 8 + 1 + 3]) != nr + 1:
                                        break
                                    else:
                                        vector1 = Vector(parameters[x * 8 + 5], parameters[x * 8 + 6], parameters[x * 8 + 7])
                                        vector2 = Vector(parameters[y * 8 + 5], parameters[y * 8 + 6], parameters[y * 8 + 7])
                                        vector3 = Vector(parameters[z * 8 + 5], parameters[z * 8 + 6], parameters[z * 8 + 7])
                                        vector4 = Vector(parameters[v * 8 + 5], parameters[v * 8 + 6], parameters[v * 8 + 7])
                                        angle = calc_dihedral(vector1, vector2, vector3, vector4)
                                        angles.write(str(nr) + ' ZETA: ' + str(angle) + '\n')

#glowna funkcja main
def main():
    infile = input("PDB file: ")
    #infile = raw_input("PDB file: ")
    download_file(infile)
    prepere_atom_file('pdb' + infile + '.ent')
    all_angles()
    csv_file()
    #remove_temp_files(infile)

main()



