from __future__ import (absolute_import, division, print_function, unicode_literals)
import csv
from Bio.PDB import calc_dihedral, Vector
'''
fieldnames = ['model', 'nr_modelu', 'oznaczenie_lancucha', 'nr_nukleotydu']
projekt = open('5no3.pdb', 'r')
only_ATOM = open('only_ATOM.txt', 'w')
n = projekt.read().split('\n')
for x in range(0 , len(n)):
    if n[x][:4] == 'ATOM' and (n[x][18:20] == ' A' or n[x][18:20] == ' C' or n[x][18:20] == ' T' or n[x][18:20] == ' G' or n[x][18:20] == ' U'):
        only_ATOM.write(n[x])
        only_ATOM.write('\n')
only_ATOM.close()
'''

'''
only_ATOM = open('only_ATOM.txt', 'r')
only_NUC =  open('only_NUC.txt', 'w')
n = only_ATOM.read().split('\n')
for x in range(0 , len(n)):
    if n[x][:4] == 'ATOM':
        only_NUC.write(n[x])
        only_NUC.write('\n')
only_ATOM.close()
only_NUC.close()
'''
'''
projekt2 =  open('5no3.pdb', 'r')
only_NUC2 =  open('only_NUC2.txt', 'w')
n3 = projekt2.read().split('\n')
for x in range(0 , len(n3)):
    if n3[x][:4] == 'ATOM' and (n3[x][18:20] == ' A' or n3[x][18:20] == ' C' or n3[x][18:20] == ' T' or n3[x][18:20] == ' G' or n3[x][18:20] == ' U'):
        only_NUC2.write(n3[x] + '\n')
projekt2.close()
only_NUC2.close()
'''
'''
only_NUC2 =  open('only_NUC2.txt', 'r')
α,β,γ,δ,ε,ζ = 0, 0, 0, 0, 0, 0
tor = open('tor.txt', 'w')
n4 = only_NUC2.read().split('\n')
for x in range(0, len(n4)):
    tor.write(n4[x][6:11] + ' ' + n4[x][12:16] + ' ' + n4[x][17:20] + ' ' + n4[x][21:22] + ' ' + n4[x][22:26] + ' ' + n4[x][30:38] + ' ' + n4[x][38:46] + ' ' + n4[x][46:54] + '\n')

tor.close()
only_NUC2.close()
'''

'''
tor = open('tor.txt', 'r')
with open('tor1.txt', 'w') as tor1:
    k2 = tor.read().split()
    ccc = int(len(k2)/8)
    #ALFA
    tor1.write('\n')
    tor1.write('ALFA' + '\n' + '\n')
    for x in range(0, ccc):
        if k2[x*8+1] == 'O3\'':
            nr = int(k2[x * 8 + 1 + 3]) + 1
            for y in range(x, ccc):
                if k2[y*8 + 1] == 'P':
                    if int(k2[y * 8 + 1 + 3]) != nr:
                        break
                    else:
                        for z in range(y, ccc):
                            if k2[z*8 +1] == 'O5\'':
                                if int(k2[z * 8 + 1 + 3]) != nr:
                                    break
                                else:
                                    for v in range(z, ccc):
                                        if k2[v * 8 + 1] == 'C5\'' and  int(k2[v * 8 + 1 + 3]) == nr:
                                            # plik.write(k2[x*8] + ' ' + k2[x*8 +1] + ' ' + k2[x*8 +1 + 3] + " " +
                                            # k2[y * 8] + ' ' + k2[y*8 +1] + " " + k2[y *8 +1+ 3] + ' ' +
                                            # k2[z * 8] + ' ' + " " + k2[z*8 +1] + ' ' + k2[z *8 +1+ 3] + ' ' +
                                            # k2[v * 8] + ' ' +' ' + k2[v*8 +1] + ' ' + k2[v *8 +1+ 3] + '\n')
                                            vector1 = Vector(k2[x * 8 + 5], k2[x * 8 + 6], k2[x * 8 + 7])
                                            vector2 = Vector(k2[y * 8 + 5], k2[y * 8 + 6], k2[y * 8 + 7])
                                            vector3 = Vector(k2[z * 8 + 5], k2[z * 8 + 6], k2[z * 8 + 7])
                                            vector4 = Vector(k2[v * 8 + 5], k2[v * 8 + 6], k2[v * 8 + 7])
                                            angle = calc_dihedral(vector1, vector2, vector3, vector4)
                                            tor1.write(str(angle) + '\n')
    #beta
        tor1.write('\n')
        tor1.write('BETA' + '\n' + '\n')
    for x in range(0, ccc):
        if k2[x*8+1] == 'P':
            nr = int(k2[x * 8 + 1 + 3])
            for y in range(x, ccc):
                if k2[y*8 + 1] == 'O5\'':
                    if int(k2[y * 8 + 1 + 3]) != nr:
                        break
                    else:
                        for z in range(y, ccc):
                            if k2[z*8 +1] == 'C5\'':
                                if int(k2[z * 8 + 1 + 3]) != nr:
                                    break
                                else:
                                    for v in range(z, ccc):
                                        if k2[v * 8 + 1] == 'C4\'' and  int(k2[v * 8 + 1 + 3]) == nr:
                                            vector1 = Vector(k2[x * 8 + 5], k2[x * 8 + 6], k2[x * 8 + 7])
                                            vector2 = Vector(k2[y * 8 + 5], k2[y * 8 + 6], k2[y * 8 + 7])
                                            vector3 = Vector(k2[z * 8 + 5], k2[z * 8 + 6], k2[z * 8 + 7])
                                            vector4 = Vector(k2[v * 8 + 5], k2[v * 8 + 6], k2[v * 8 + 7])
                                            angle = calc_dihedral(vector1, vector2, vector3, vector4)
                                            tor1.write(str(angle) + '\n')
    #gamma
    tor1.write('\n')
    tor1.write('GAMMA' + '\n' + '\n')
    for x in range(0, ccc):
        if k2[x*8+1] == 'O5\'':
            nr = int(k2[x * 8 + 1 + 3])
            for y in range(x, ccc):
                if k2[y*8 + 1] == 'C5\'':
                    if int(k2[y * 8 + 1 + 3]) != nr:
                        break
                    else:
                        for z in range(y, ccc):
                            if k2[z*8 +1] == 'C4':
                                if int(k2[z * 8 + 1 + 3]) != nr:
                                    break
                                else:
                                    for v in range(z, ccc):
                                        if k2[v * 8 + 1] == 'C3\'' and  int(k2[v * 8 + 1 + 3]) == nr:
                                            # plik.write(k2[x*8] + ' ' + k2[x*8 +1] + ' ' + k2[x*8 +1 + 3] + " " +
                                            # k2[y * 8] + ' ' + k2[y*8 +1] + " " + k2[y *8 +1+ 3] + ' ' +
                                            # k2[z * 8] + ' ' + " " + k2[z*8 +1] + ' ' + k2[z *8 +1+ 3] + ' ' +
                                            # k2[v * 8] + ' ' +' ' + k2[v*8 +1] + ' ' + k2[v *8 +1+ 3] + '\n')
                                            vector1 = Vector(k2[x * 8 + 5], k2[x * 8 + 6], k2[x * 8 + 7])
                                            vector2 = Vector(k2[y * 8 + 5], k2[y * 8 + 6], k2[y * 8 + 7])
                                            vector3 = Vector(k2[z * 8 + 5], k2[z * 8 + 6], k2[z * 8 + 7])
                                            vector4 = Vector(k2[v * 8 + 5], k2[v * 8 + 6], k2[v * 8 + 7])
                                            angle = calc_dihedral(vector1, vector2, vector3, vector4)
                                            tor1.write(str(angle) + '\n')
    #delta
    tor1.write('\n')
    tor1.write('DELTA' + '\n' + '\n')
    for x in range(0, ccc):
        if k2[x*8+1] == 'C5\'':
            nr = int(k2[x * 8 + 1 + 3])
            for y in range(x, ccc):
                if k2[y*8 + 1] == 'O4\'':
                    if int(k2[y * 8 + 1 + 3]) != nr:
                        break
                    else:
                        for z in range(y, ccc):
                            if k2[z*8 +1] == 'C3\'':
                                if int(k2[z * 8 + 1 + 3]) != nr:
                                    break
                                else:
                                    for v in range(z, ccc):
                                        if k2[v * 8 + 1] == 'O3\'' and  int(k2[v * 8 + 1 + 3]) == nr:
                                            vector1 = Vector(k2[x * 8 + 5], k2[x * 8 + 6], k2[x * 8 + 7])
                                            vector2 = Vector(k2[y * 8 + 5], k2[y * 8 + 6], k2[y * 8 + 7])
                                            vector3 = Vector(k2[z * 8 + 5], k2[z * 8 + 6], k2[z * 8 + 7])
                                            vector4 = Vector(k2[v * 8 + 5], k2[v * 8 + 6], k2[v * 8 + 7])
                                            angle = calc_dihedral(vector1, vector2, vector3, vector4)
                                            tor1.write(str(angle) + '\n')
    #EPSILON
    tor1.write('\n')
    tor1.write('EPSILON' + '\n' + '\n')
    for x in range(0, ccc):
        if k2[x*8+1] == 'C4\'':
            nr = int(k2[x * 8 + 1 + 3])
            for y in range(x, ccc):
                if k2[y*8 + 1] == 'C3':
                    if int(k2[y * 8 + 1 + 3]) != nr:
                        break
                    else:
                        for z in range(y, ccc):
                            if k2[z*8 +1] == 'O3\'':
                                if int(k2[z * 8 + 1 + 3]) != nr:
                                    break
                                else:
                                    for v in range(z, ccc):
                                        if k2[v * 8 + 1] == 'P' and int(k2[v * 8 + 1 + 3]) == nr+1:
                                            vector1 = Vector(k2[x * 8 + 5], k2[x * 8 + 6], k2[x * 8 + 7])
                                            vector2 = Vector(k2[y * 8 + 5], k2[y * 8 + 6], k2[y * 8 + 7])
                                            vector3 = Vector(k2[z * 8 + 5], k2[z * 8 + 6], k2[z * 8 + 7])
                                            vector4 = Vector(k2[v * 8 + 5], k2[v * 8 + 6], k2[v * 8 + 7])
                                            angle = calc_dihedral(vector1, vector2, vector3, vector4)
                                            tor1.write(str(angle) + '\n')
    #ZETA
    tor1.write('\n')
    tor1.write('ZETA' + '\n' + '\n')
    for x in range(0, ccc):
        if k2[x*8+1] == 'C3\'':
            nr = int(k2[x * 8 + 1 + 3])
            for y in range(x, ccc):
                if k2[y*8 + 1] == 'O3\'':
                    if int(k2[y * 8 + 1 + 3]) != nr:
                        break
                    else:
                        for z in range(y, ccc):
                            if k2[z*8 +1] == 'P':
                                if int(k2[z * 8 + 1 + 3]) != nr+1:
                                    break
                                else:
                                    for v in range(z, ccc):
                                        if k2[v * 8 + 1] == 'O5\'' and  int(k2[v * 8 + 1 + 3]) == nr+1:
                                            # plik.write(k2[x*8] + ' ' + k2[x*8 +1] + ' ' + k2[x*8 +1 + 3] + " " +
                                            # k2[y * 8] + ' ' + k2[y*8 +1] + " " + k2[y *8 +1+ 3] + ' ' +
                                            # k2[z * 8] + ' ' + " " + k2[z*8 +1] + ' ' + k2[z *8 +1+ 3] + ' ' +
                                            # k2[v * 8] + ' ' +' ' + k2[v*8 +1] + ' ' + k2[v *8 +1+ 3] + '\n')
                                            vector1 = Vector(k2[x * 8 + 5], k2[x * 8 + 6], k2[x * 8 + 7])
                                            vector2 = Vector(k2[y * 8 + 5], k2[y * 8 + 6], k2[y * 8 + 7])
                                            vector3 = Vector(k2[z * 8 + 5], k2[z * 8 + 6], k2[z * 8 + 7])
                                            vector4 = Vector(k2[v * 8 + 5], k2[v * 8 + 6], k2[v * 8 + 7])
                                            angle = calc_dihedral(vector1, vector2, vector3, vector4)
                                            tor1.write(str(angle) + '\n')

tor.closw()
tor1.closw()
'''

'''
tor = open('tor.txt', 'r')
with open('tor2.txt', 'w') as tor2:
    k2 = tor.read().split()
    ccc = int(len(k2)/8)

    def kat(r1,n1,r2,n2,r3,n3,r4,n4, ccc, k2):
        for x in range(0, ccc):
            if k2[x*8+1] == r1:
                nr = int(k2[x * 8 + 1 + 3]) + n1
                for y in range(x, ccc):
                    if k2[y*8 + 1] == r2:
                        if int(k2[y * 8 + 1 + 3]) != nr+n2:
                            break
                        else:
                            for z in range(y, ccc):
                                if k2[z*8 +1] == r3:
                                    if int(k2[z * 8 + 1 + 3]) != nr+n3:
                                        break
                                    else:
                                        for v in range(z, ccc):
                                            if k2[v * 8 + 1] == r4 and  int(k2[v * 8 + 1 + 3]) == nr + n4:
                                                vector1 = Vector(k2[x * 8 + 5], k2[x * 8 + 6], k2[x * 8 + 7])
                                                vector2 = Vector(k2[y * 8 + 5], k2[y * 8 + 6], k2[y * 8 + 7])
                                                vector3 = Vector(k2[z * 8 + 5], k2[z * 8 + 6], k2[z * 8 + 7])
                                                vector4 = Vector(k2[v * 8 + 5], k2[v * 8 + 6], k2[v * 8 + 7])
                                                angle = calc_dihedral(vector1, vector2, vector3, vector4)
                                                tor2.write(str(angle) + '\n')

    #ALFA
    tor2.write('\n')
    tor2.write('ALFA' + '\n' + '\n')
    kat('O3\'', 1, 'P', 0 ,'O5\'', 0, 'C5\'', 0, ccc, k2)
    #BETA
    tor2.write('\n')
    tor2.write('BETA' + '\n' + '\n')
    kat('P', 0, 'O5\'', 0 ,'C5\'', 0, 'C4\'', 0, ccc, k2)
    #GAMMA
    tor2.write('\n')
    tor2.write('GAMMA' + '\n' + '\n')
    kat('O5\'', 0, 'C5\'', 0 ,'C4', 0, 'C3\'', 0, ccc, k2)
    #DELTA
    tor2.write('\n')
    tor2.write('DELTA' + '\n' + '\n')
    kat('C5\'', 0, 'O4\'', 0 ,'C3\'', 0, 'O3\'', 0, ccc, k2)
    #EPSILON
    tor2.write('\n')
    tor2.write('EPSILON' + '\n' + '\n')
    kat('C4\'', 0, 'C3', 0 ,'O3\'', 0, 'P', 1, ccc, k2)
    #ZETA
    tor2.write('\n')
    tor2.write('ZETA' + '\n' + '\n')
    kat('C3\'', 0, 'O3\'', 0 ,'P', 1, 'O5\'', 1, ccc ,k2)
'''

tor = open('tor.txt', 'r')
with open('tor3.txt', 'w') as tor3:
    k2 = tor.read().split()
    ccc = int(len(k2)/8)

    def wszystkiekaty(ccc, k2):
        for x in range(0, ccc):
            alfa, beta, gamma, delta, epsilon, zeta, find, nr = False, False, False, False, False, False, False, 0
            if k2[x*8+1] == 'O3\'':
                alfa = True
                find = True
                nr = 1
            elif k2[x*8+1] == 'P':
                beta = True
                find = True
            elif k2[x * 8 + 1] == 'O5\'':
                gamma = True
                find = True
            elif k2[x * 8 + 1] == 'C5\'':
                delta = True
                find = True
            elif k2[x * 8 + 1] == 'C4\'':
                epsilon = True
                find = True
            elif k2[x * 8 + 1] == 'C3\'':
                zeta = True
                find = True
            if find:
                nr = nr + int(k2[x * 8 + 1 + 3])
                for y in range(x, ccc):
                    find = False
                    if alfa == True and k2[y*8 + 1] == 'P':
                        if int(k2[y * 8 + 1 + 3]) != nr:
                            break
                        else: find = True
                    elif beta == True and k2[y*8 + 1] == 'O5\'':
                        if int(k2[y * 8 + 1 + 3]) != nr:
                            break
                        else:
                            find = True
                    elif gamma == True and k2[y * 8 + 1] == 'C5\'':
                        if int(k2[y * 8 + 1 + 3]) != nr:
                            break
                        else:
                            find = True
                    elif delta == True and k2[y * 8 + 1] == 'O4\'':
                        if int(k2[y * 8 + 1 + 3]) != nr:
                            break
                        else:
                            find = True
                    elif epsilon == True and k2[y * 8 + 1] == 'C3':
                        if int(k2[y * 8 + 1 + 3]) != nr:
                            break
                        else:
                            find = True
                    elif zeta == True and k2[y * 8 + 1] == 'O3\'':
                        if int(k2[y * 8 + 1 + 3]) != nr:
                            break
                        else:
                            find = True
                    if find:
                        for z in range(y, ccc):
                            find = False
                            if alfa == True and k2[z * 8 + 1] == 'O5\'':
                                if int(k2[z * 8 + 1 + 3]) != nr:
                                    break
                                else:
                                    find = True
                            elif beta == True and k2[z * 8 + 1] == 'C5\'':
                                if int(k2[z * 8 + 1 + 3]) != nr:
                                    break
                                else:
                                    find = True
                            elif gamma == True and k2[z * 8 + 1] == 'C4':
                                if int(k2[z * 8 + 1 + 3]) != nr:
                                    break
                                else:
                                    find = True
                            elif delta == True and k2[z * 8 + 1] == 'C3\'':
                                if int(k2[z * 8 + 1 + 3]) != nr:
                                    break
                                else:
                                    find = True
                            elif epsilon == True and k2[z * 8 + 1] == 'O3\'':
                                if int(k2[z * 8 + 1 + 3]) != nr:
                                    break
                                else:
                                    find = True
                            elif zeta == True and k2[z * 8 + 1] == 'P':
                                if int(k2[z * 8 + 1 + 3]) != nr + 1:
                                    break
                                else:
                                    find = True
                            if find:
                                for v in range(z, ccc):
                                    if alfa == True and k2[v * 8 + 1] == 'C5\'':
                                        if int(k2[v * 8 + 1 + 3]) != nr:
                                            break
                                        else:
                                            vector1 = Vector(k2[x * 8 + 5], k2[x * 8 + 6], k2[x * 8 + 7])
                                            vector2 = Vector(k2[y * 8 + 5], k2[y * 8 + 6], k2[y * 8 + 7])
                                            vector3 = Vector(k2[z * 8 + 5], k2[z * 8 + 6], k2[z * 8 + 7])
                                            vector4 = Vector(k2[v * 8 + 5], k2[v * 8 + 6], k2[v * 8 + 7])
                                            angle = calc_dihedral(vector1, vector2, vector3, vector4)
                                            tor3.write('ALPHA: ' + str(angle) + '\n')
                                    elif beta == True and k2[v * 8 + 1] == 'C4\'':
                                        if int(k2[v * 8 + 1 + 3]) != nr:
                                            break
                                        else:
                                            vector1 = Vector(k2[x * 8 + 5], k2[x * 8 + 6], k2[x * 8 + 7])
                                            vector2 = Vector(k2[y * 8 + 5], k2[y * 8 + 6], k2[y * 8 + 7])
                                            vector3 = Vector(k2[z * 8 + 5], k2[z * 8 + 6], k2[z * 8 + 7])
                                            vector4 = Vector(k2[v * 8 + 5], k2[v * 8 + 6], k2[v * 8 + 7])
                                            angle = calc_dihedral(vector1, vector2, vector3, vector4)
                                            tor3.write('BETA: ' + str(angle) + '\n')
                                    elif gamma == True and k2[v * 8 + 1] == 'C3\'':
                                        if int(k2[v * 8 + 1 + 3]) != nr:
                                            break
                                        else:
                                            vector1 = Vector(k2[x * 8 + 5], k2[x * 8 + 6], k2[x * 8 + 7])
                                            vector2 = Vector(k2[y * 8 + 5], k2[y * 8 + 6], k2[y * 8 + 7])
                                            vector3 = Vector(k2[z * 8 + 5], k2[z * 8 + 6], k2[z * 8 + 7])
                                            vector4 = Vector(k2[v * 8 + 5], k2[v * 8 + 6], k2[v * 8 + 7])
                                            angle = calc_dihedral(vector1, vector2, vector3, vector4)
                                            tor3.write('GAMMA: ' + str(angle) + '\n')
                                    elif delta == True and k2[v * 8 + 1] == 'O3\'':
                                        if int(k2[v * 8 + 1 + 3]) != nr:
                                            break
                                        else:
                                            vector1 = Vector(k2[x * 8 + 5], k2[x * 8 + 6], k2[x * 8 + 7])
                                            vector2 = Vector(k2[y * 8 + 5], k2[y * 8 + 6], k2[y * 8 + 7])
                                            vector3 = Vector(k2[z * 8 + 5], k2[z * 8 + 6], k2[z * 8 + 7])
                                            vector4 = Vector(k2[v * 8 + 5], k2[v * 8 + 6], k2[v * 8 + 7])
                                            angle = calc_dihedral(vector1, vector2, vector3, vector4)
                                            tor3.write('DELTA: ' + str(angle) + '\n')
                                    elif epsilon == True and k2[v * 8 + 1] == 'P':
                                        if int(k2[v * 8 + 1 + 3]) != nr + 1:
                                            break
                                        else:
                                            vector1 = Vector(k2[x * 8 + 5], k2[x * 8 + 6], k2[x * 8 + 7])
                                            vector2 = Vector(k2[y * 8 + 5], k2[y * 8 + 6], k2[y * 8 + 7])
                                            vector3 = Vector(k2[z * 8 + 5], k2[z * 8 + 6], k2[z * 8 + 7])
                                            vector4 = Vector(k2[v * 8 + 5], k2[v * 8 + 6], k2[v * 8 + 7])
                                            angle = calc_dihedral(vector1, vector2, vector3, vector4)
                                            tor3.write('EPSILON: ' + str(angle) + '\n')
                                    elif zeta == True and k2[v * 8 + 1] == 'O5\'':
                                        if int(k2[v * 8 + 1 + 3]) != nr + 1:
                                            break
                                        else:
                                            vector1 = Vector(k2[x * 8 + 5], k2[x * 8 + 6], k2[x * 8 + 7])
                                            vector2 = Vector(k2[y * 8 + 5], k2[y * 8 + 6], k2[y * 8 + 7])
                                            vector3 = Vector(k2[z * 8 + 5], k2[z * 8 + 6], k2[z * 8 + 7])
                                            vector4 = Vector(k2[v * 8 + 5], k2[v * 8 + 6], k2[v * 8 + 7])
                                            angle = calc_dihedral(vector1, vector2, vector3, vector4)
                                            tor3.write('ZETA: ' + str(angle) + '\n')

    wszystkiekaty(ccc, k2)
#'''