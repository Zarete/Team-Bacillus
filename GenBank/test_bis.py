# -*- coding: utf-8 -*-

import myBio as bio
import myProject as proj

import os

print('Test_bis.py')

seq = proj.readGenBank('sequences/sequence(7).gb')

print(len(seq['genes']))

somme = 0
i = 0
for elem in seq['genes']:
    somme += elem['length']
    i += 1
    print(elem['protein'])

print(somme/i)