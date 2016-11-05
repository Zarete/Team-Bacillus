# -*- coding: utf-8 -*-

# Thomas Maucourt

import myBio as bio

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
    section = txt.split('     gene     ')

    for i in range(1, 2, 1):
        dic_result = {'start' : 0,'stop' : 0, 'frame' : 0, 'length' : 0, 'name' : 'unknown', 'protein' : 'xxx', 'product' : 'unknown'}
        for elem in section[i]:
            print(elem)
            if elem[0].isdigit():
                part = elem.split('..')
                dic_result['start'] = part[0]
                dic_result['stop'] = part[1]
                dic_result['length'] = int(part[1])-int(part[0])
            elif 'product' in section[i]:
                part = section[i].split('=')
                dic_result['product'] = part[1]
        result.append(dic_result)
    print("\n", result)


    """for i in range(1, 2, 1):

        dic_result = {'start' : 0,'stop' : 0, 'frame' : 0, 'length' : 0, 'name' : 'unknown', 'protein' : 'xxx', 'product' : 'unknown'}

        section_tmp = section[i].replace('\n', '').replace(' ','').split("/")
        print(section_tmp)
        start_stop = section_tmp[0].split('..')
        dic_result['start'] = start_stop[0]
        dic_result['stop'] = start_stop[1]
        #dic_result['frame'] =
        dic_result['length'] = int(start_stop[1])-int(start_stop[0])
        #dic_result['name'] =
        #dic_result['protein'] =
        #dic_result['product'] = section_tmp[]
        result.append(dic_result)"""
