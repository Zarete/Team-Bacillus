# -*- coding: utf-8 -*-

import myBio as bio
import myProject as proj

import os

#seq = proj.readGenBank('sequences/sequencelol.gb')

"""
for a in seq['genes']:
    for b, c in a.items():
        print(b)
        print(c)
        print("\n\n")
"""

lst_fichier = os.listdir('sequences/')
"""
#txt = proj.readFlatFile('sequence.gb')
#print(txt) # ← 'LOCUS       NM_000518 ' ...

#entry = proj.readFlatFile('sequence.gb')
#features = proj.getFeatures(entry)
#print(features) # '      source          1..626 '...

entry = proj.readFlatFile('test.gb')
features = proj.getFeatures(entry)
genes = proj.getGenes(features)
 # Number of genes described in the features section
"""
datas = []
for elem in lst_fichier:
    datas.append(proj.readGenBank('sequences/'+elem))

"""print(datas[0]['genes']["name"])"""

for elem in datas:
    print(elem['genes'])


    print("\n[+++++]", elem, "[+++++]\n")
    
    print('ID :',seq['ID'])
    print('#GENES :',seq['length'])
    print('TYPE :',seq['type'])
    print('ORGANISM :',seq['organism'])
    print('DESCRIPTION :',seq['description'])
    print('GBTYPE :',seq['gbtype'])
    print('DATA :',seq['data'])
    print('\n####################\n')

    for a in seq['genes']:
        print('NAME :', a['name'])
        print('START :', a['start'])
        print('STOP :', a['stop'])
        print('LENGTH :', a['length'])
        print('FRAME :', a['frame'])
        print('PRODUCT :', a['product'])
        print('PROTEIN :', a['protein'])
        print()

    print("\n##############################################\n")

print("\n\n\n ALL WENT GOOD !! \n\n")
"""
seq = proj.readGenBank('sequences/sequence.gb')
print(seq['ID'])
print(seq['length'])
print(seq['type'])
print(seq['organism'])
print(len(seq['genes']) ) # Number of genes
print(seq['genes'][0]['start'])
print(seq['genes'][0]['stop'])
print(seq['genes'][0]['frame']) # 1, 2, 3, -1, -2, or -3
print(seq['genes'][0]['length']) # expressed in bp
print(seq['genes'][0]['protein']) # M.....*
print(seq['genes'][0]['product']) # protein name
print(seq['genes'][0]['name']) # gene name
print(seq['genes'][1]['start'])
"""
