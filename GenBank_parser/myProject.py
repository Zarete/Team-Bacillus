# -*- coding: utf-8 -*-

# Thomas Maucourt

import myBio as bio

def clean_size(txt):
    """ Gene size cleaner

    Cleans the size information provided to the function,
    this function is written by Maucourt Thomas

    Args :
        txt : text extracted from the line containing the gene size

    Return :
        return a string containing numbers only
    """

    if txt.isdigit():
        return txt

    else:
        nb = ''
        for char in txt.split('\n')[-1]:
            if char.isdigit():
                nb += char
        return nb

def readGenBank(filename):
    """GenBank flat file parser

    Extract all the informations about a GenBank flat file,
    this function is written by Maucourt Thomas

    Args :
        filename : GenBank flat file
    Return :
        return a dictionnary containing all the informations about the
        sequence in the file
    """

    result = {"description" : '', "type" : '', "data" : '', "gbtype" : '',
              "ID" : '', 'length' : '', 'type' : '', 'organism' : '', "codeTableID" : '',
              "genes" : []}

    with open(filename, 'r') as file_content:
        data = file_content.read()

    features = data[data.index('FEATURES'): data.index('ORIGIN')]
    header = data[data.index('LOCUS'): data.index('FEATURES')]

    # header parsing and extraction
    splitted_header = header.split('\n')

    for elem in splitted_header:
        if 'DEFINITION' in elem:
            tmp = elem.strip().split('DEFINITION')
            result['description'] = tmp[1].strip()
        elif 'LOCUS' in elem:
            tmp = elem.strip().split()

    # features parsing and extraction
    splitted_features = features.split('     gene     ')

    lst_pos = []
    lst_genes = []

    """for i in range(len(splitted_features)):
        if '     gene     ' in splitted_features[i]:
            lst_pos.append(i)

    for i in range(0,len(lst_pos)-1,1):
        tmp = splitted_features[lst_pos[i]:lst_pos[i+1]]
        lst_genes.append(tmp)
    lst_genes.append(splitted_features[-1:])"""

    for elem in splitted_features[1:]:

        dic_result = {'start' : 0,'stop' : 0, 'frame' : 0, 'length' : 0, 'name' : 'unknown', 'protein' : 'xxx', 'product' : 'unknown'}

        if "product=" in elem:
            product = elem.split("/product=")[1].replace(" "*20, "").replace('\n', '').split('"')[1]
            dic_result['product'] = product

        section = elem.split('/')

        for part in section:

            if '..' in part:
                start = clean_size(part.strip().split('..')[0])
                stop = clean_size(part.strip().split('..')[1])
                dic_result['start'] = int(start)
                dic_result['stop'] = int(stop)
                dic_result['length'] = int(stop) - int(start)

            elif 'codon_start=' in part:
                dic_result['frame'] = part.strip().replace('\n', '').split('=')[1].strip()

            elif 'gene=' in part:
                name = part.split("=")
                dic_result['name'] = name[1].replace('\n', '').strip()[1:-1]
            elif 'translation=' in part:
                protein = part.strip().replace("\n", "").split('\"')[1]
                dic_result['protein'] = protein.replace(' ', '')

        result['genes'].append(dic_result)

    for elem in result['genes']:
        print("start    ",elem['start'])
        print("stop     ",elem['stop'])
        print("length   ",elem['length'])
        print("frame    ",elem['frame'])
        print("name     ",elem['name'])
        print("product  ",elem['product'])
        print("protein  ",elem['protein'])
        print('\n\n')

    return result
