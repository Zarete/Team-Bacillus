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

    result = {"description" : '', "type" : '', "data" : 'xxx', "gbtype" : '',
              "ID" : '', 'length' : '', 'type' : '', 'organism' : '', "codeTableID" : '',
              "genes" : []}

    with open(filename, 'r') as file_content:
        data = file_content.read()

    features = data[data.index('FEATURES'): data.index('ORIGIN')]
    header = data[data.index('LOCUS'): data.index('     gene     ')]

    ########## Header parsing and extraction ##########

    splitted_header = header.split('\n')

    if '##Genome-Annotation-Data-START##' in header and '##Genome-Annotation-Data-END##' in header:
        result['data'] = header[header.index('##Genome-Annotation-Data-START##')+33: header.index('##Genome-Annotation-Data-END##')]

    for elem in splitted_header:
        if 'DEFINITION' in elem:
            tmp = elem.strip().split('DEFINITION')
            result['description'] = tmp[1].strip()
        elif '/mol_type' in elem:
            moltype = elem.split('"')[1].strip()
            result['gbtype'] = moltype

            if 'RNA' in moltype.upper():
                otype = 'rna'
            elif 'DNA' in moltype.upper():
                otype = 'dna'
            elif 'PROTEIN' in moltype.upper():
                otype = 'protein'
            else:
                otype = ''
            result['type'] = otype
        elif 'ACCESSION' in elem:
            ID = elem.split("ACCESSION")[1].strip()
            result['ID'] = ID
        elif '/organism' in elem:
            organism = elem.split('"')[1].strip()
            result['organism'] = organism
        elif '     source     ' in elem:
            size = elem.split('..')[1].strip()
            result['length'] = size

    ########## Features parsing and extraction ##########

    splitted_features = features.split('     gene     ')

    for elem in splitted_features[1:]:

        dic_result = {'start' : 0,'stop' : 0, 'frame' : 0, 'length' : 0, 'name' : 'unknown', 'protein' : 'xxx', 'product' : 'unknown'}

        if "product=" in elem:
            product = elem.split("/product=")[1].replace(" "*20, "").replace('\n', '').split('"')[1]
            dic_result['product'] = product

        sub_section = elem.split('\n')[0]
        if '..' in sub_section and '/' not in sub_section:
            if ',' in sub_section:
                sub, final = [], []
                under = sub_section.split(',')
                for elem in under:
                    sub.extend(elem.split('..'))

                for elem in sub:
                    final.append(clean_size(elem))

                dic_result['start'] = int(final[1])
                dic_result['stop'] = int(final[3])
                dic_result['length'] = (int(final[1]) - int(final[0])) + (int(final[3]) - int(final[2]))

            else:
                start = clean_size(sub_section.strip().split('..')[0])
                stop = clean_size(sub_section.strip().split('..')[1])
                dic_result['start'] = int(start)
                dic_result['stop'] = int(stop)
                dic_result['length'] = int(stop) - int(start)

        section = elem.split('/')

        for part in section:

            if 'codon_start=' in part:
                dic_result['frame'] = part.strip().replace('\n', '').split('=')[1].strip()

            elif 'gene=' in part:
                name = part.split("=")
                dic_result['name'] = name[1].replace('\n', '').strip()[1:-1]
            elif 'translation=' in part:
                protein = part.strip().replace("\n", "").split('\"')[1]
                dic_result['protein'] = protein.replace(' ', '')

        result['genes'].append(dic_result)

    return result

    """if '..' in part.split('\n')[0]:
        start, stop, sub = 0, 0, []
        if 'join(' in part or 'complement(' in part:
            sub = part.replace('\n', '').split('join(')[1].split(',')
            under = []
            for elem in sub:
                print(elem)
                a, b, tmp = 0, 0, []
                tmp = elem.split('..')
                a = clean_size(tmp[0])
                b = clean_size(tmp[1])
                under.append(a)
                under.append(b)
            under = sorted(under)
            start = under[3]
            stop = under[1]
            dic_result['start'] = start
            dic_result['stop'] = stop
            dic_result['length'] = (under[3] - under[2]) + (under[1] - under[0])

        else:
            start = part[0]
            print('f',start)
            start = clean_size(start)
            print('s',start)
            #print(clean_size(start[0]))
            stop = part[1]
            print('f',stop)
            stop = clean_size(stop)
            print('s',stop)
            #print(clean_size(stop[1]))
            #start = clean_size(part.replace("\n", "").strip().split('\n')[0].split('..')[1])
            #stop = clean_size(part.replace("\n", "").strip().split('\n')[0].split('..')[1])
            dic_result['start'] = start
            dic_result['stop'] = stop
            dic_result['length'] = stop - start"""
