# -*- coding: utf-8 -*-

# Thomas Maucourt

import myBio as bio

def clean_size(txt):
    if txt.isdigit():
        return txt

    else:
        nb = ''
        for char in txt:
            if char.isdigit():
                nb += char
        return nb

def readFlatFile(filename):
    """ Flat file reader (format GenBank)

    Loads the datas from the file into memory,
    this function is written by Maucourt Thomas

    Args :
        filename : file to extract information from

    Return :
        return the datas extracted as a string
    """

    with open(filename, 'r') as file_content:
        result = file_content.read()
    return str(result)

def getFeatures(txt, specific = "auto"):
    """Features extraction

    Extract the features of the text given in
    parameter,
    this function is written by Maucourt Thomas

    Args :
        txt : string containing the sequence and it's informations

    Return :
        return a string containg the features of the sequence
    """

    pos1 = txt.index("FEATURES")
    pos2 = txt.index('ORIGIN')

    result = txt[pos1:pos2]

    return result

def getGenes(txt):
    """Extract informations about genes & CDS

    Extract all the informations about genes & CDS
    of the parameter txt,
    this function is written by Maucourt Thomas

    Args :
        txt : string containing the sequence and it's informations

    Return :
        return a list of dictionnaries containing all the informations
        about gene & CDS sections
    """

    result = []

    section = txt.split('\n')
    lst_pos = []
    lst_genes = []

    for i in range(len(section)):
        if '     gene     ' in section[i]:
            lst_pos.append(i)

    for i in range(0, len(lst_pos)-1, 2):
        tmp = section[lst_pos[i]:lst_pos[i+1]]
        lst_genes.append(tmp)

    for elem in lst_genes:
        dic_result = {'start' : 0,'stop' : 0, 'frame' : 0, 'length' : 0, 'name' : 'unknown', 'protein' : 'xxx', 'product' : 'unknown'}
        for part in elem:
            if '     gene     ' in part:
                part = part.split('gene')
            elif '     tRNA     ' in part:
                part = part.split('tRNA')
            elif '     tRNA     ' in part:
                part = part.split('tRNA')

            part[1] = part[1].strip().split('..')
            part[1][0] = clean_size(part[1][0])
            part[1][1] = clean_size(part[1][1])
            dic_result['start'] = int(part[1][0])
            dic_result['stop'] = int(part[1][1])
            dic_result['length'] = int(part[1][1]) - int(part[1][0])

            if 'product' in part:
                part = part.split('=')
                dic_result['product'] = part[1][1:-1]
            elif 'codon_start' in part:
                part = part.split('=')
                dic_result['frame'] = part[1]

        result.append(dic_result)
    print(result)
    print(len(result))
