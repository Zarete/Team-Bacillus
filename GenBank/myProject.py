# -*- coding: utf-8 -*-

# Thomas Maucourt

import myBio as bio

def clean_numbers(txt):
    """ Gene size cleaner
    
    Cleans the size information provided to the function,
    this function is written by Maucourt Thomas
    Args :
        txt : text extracted from the line containing the gene size
    Return :
        return a string containing numbers only
    """

    if txt.isdigit():
        return int(txt)

    else:
        nb = ''
        for char in txt:
            if char.isdigit():
                nb += char

        return int(nb)

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

    ########## File loading ##########

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
            result['description'] = elem.strip().split('DEFINITION')[1].strip()

        elif '/mol_type' in elem:
            result['gbtype'] = elem.split('"')[1].strip()

            if 'RNA' in result['gbtype']:
                otype = 'rna'

            elif 'DNA' in result['gbtype']:
                otype = 'dna'

            elif 'PROTEIN' in result['gbtype']:
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
            dic_result['product'] = elem.split("/product=")[1].replace(" "*20, "").replace('\n', '').split('"')[1]

        # Extract start, stop and length of the gene
        sub_section = elem.split('\n')[0]
        if '..' in sub_section:
            start, stop, length = 0, 0, 0

            if ',' in sub_section:
                tmp, final = [], []
                for elem in sub_section.split(','):
                    tmp.extend(elem.split('..'))

                for elem in tmp:
                    final.append(clean_numbers(elem))

                final = sorted(final)

                start = final[1]
                stop = final[3]
                length = (final[1] - final[0]) + (final[3] - final[2])

            else:
                start = clean_numbers(sub_section.strip().split('..')[0])
                stop = clean_numbers(sub_section.strip().split('..')[1])
                length = stop - start

            dic_result['start'] = start
            dic_result['stop'] = stop
            dic_result['length'] = length

        # Extract frame, gene name and translated protein
        section = elem.split('/')
        for part in section:
            if 'codon_start=' in part:
                dic_result['frame'] = part.strip().replace('\n', '').split('=')[1].strip()

            elif 'gene=' in part:
                name = part.split("=")
                dic_result['name'] = name[1].replace('\n', '').strip()[1:-1]

            elif 'translation=' in part:
                dic_result['protein'] = part.strip().replace("\n", "").split('\"')[1].replace(' ', '')

        result['genes'].append(dic_result)

    return result
