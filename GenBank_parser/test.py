# -*- coding: utf-8 -*-

import myBio as bio
import myProject as proj

#txt = proj.readFlatFile('sequence.gb')
#print(txt) # ← 'LOCUS       NM_000518 ' ...

#entry = proj.readFlatFile('sequence.gb')
#features = proj.getFeatures(entry)
#print(features) # '      source          1..626 '...

#entry = proj.readFlatFile('sequence.gb')
#features = proj.getFeatures(entry)
#genes = proj.getGenes(features)
#print(len(genes)) # Number of genes described in the features section

seq = proj.readGenBank('sequence.gb')
print(seq['ID'])
print(seq['length'])
print(seq['type'])
print(seq['organism'])
print(seq['description'])
print()
print(len(seq['genes']) ) # Number of genes
print(seq['genes'][0]['start'])
print(seq['genes'][0]['stop'])
print(seq['genes'][0]['frame']) # 1, 2, 3, -1, -2, or -3
print(seq['genes'][0]['length']) # expressed in bp
print(seq['genes'][0]['protein']) # M.....*
print(seq['genes'][0]['product']) # protein name
print(seq['genes'][0]['name']) # gene name
print(seq['genes'][1]['start'])
