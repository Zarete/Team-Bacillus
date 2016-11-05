# -*- coding: utf-8 -*-

import myBio as bio
import myProject as proj

#txt = proj.readFlatFile('sequence.gb')
#print(txt) # ← 'LOCUS       NM_000518 ' ...

#entry = proj.readFlatFile('sequence.gb')
#features = proj.getFeatures(entry)
#print(features) # '      source          1..626 '...

entry = proj.readFlatFile('sequence.gb')
features = proj.getFeatures(entry)
genes = proj.getGenes(features)
#print(len(genes)) # Number of genes described in the features section
