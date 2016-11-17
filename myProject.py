# -*- coding: utf-8 -*-
#Problem 4.6.1.1
import myBio as bio
import copy
import random

##### ORF Parser #####

#Recupere la table standard ou celle du mycoplasme
def getGeneticCode(NCBI_ID):
    
    """Get the Genetic Code wanted 

    Two genetics code are implemented in the function, the standard code and the code for mycoplasma.
    This function is written by Julien Benetti.

    Args:
        NCBI_ID: Genetic ID from the NCBI as an integer (1=Standard, 4=Mycoplasma)

    Returns:
        The function getGeneticCode returns the codon Table in a dictionnary with the codon
        as keys and the amino acid as value.
    """
    
    if NCBI_ID==1:
        base1='TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
        base2='TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
        base3='TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'
        aas  ='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
        table = {}
        for i in range(0,len(aas)):
            codon = base1[i]+base2[i]+base3[i]
            aa = aas[i]
            table[codon] = aa
        return table
    
    if NCBI_ID==4:
        base1='TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
        base2='TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
        base3='TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'
        aas  ='FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
        table = {}
        for i in range(0,len(aas)):
            codon = base1[i]+base2[i]+base3[i]
            aa = aas[i]
            table[codon] = aa
        return table

def translate(seq,codonTable=1):
    codonT=getGeneticCode(codonTable)
    seq_word=[]
    seq_translated=''
    for i in range(0,len(seq),3):
        seq_word.append(seq[i:i+3])
    for i in seq_word:
        if len(i)==3:
            seq_translated += codonT[i]
    return seq_translated


def findORF(seq, threshold,codeTable):
        
    """Find the ORF of the sequence

    The algorithm searches all the start and stop codons and keeps track of their locations and frames
    (locations stored in different lists according to their frames). Each Start codons are associated with
    Stop codons of the same frame in two lists (one for the start codons, one for the stop codons).
    Each start codon between a start and a stop of the same frame are not added to the list. Each stop
    codon between a stop codon and a start codon of the same frame are not added to the list either.
    The lists contains now all the ORFs (start and stop). A threshold is applied to delete each ORFs with
    a short length. All the information about the ORF is stocked in a dictionnary.
    This function is written by Julien Benetti.

    Args:
        seq: Sequence of nucleobases as a dictionnary (see fromFasta in myBio).
        threshold: the minimum ORF length expressed in base pairs as a integer.
        codeTable: Genetic code Table (see getGeneticCode)

    Returns:
        The function findORF returns a list of dictionnaries. Each dictionnary represents an ORF with 
        information stocked in different keys : start, stop, length, frame, protein
    """
    codons=getGeneticCode(codeTable)
    codonstop=[]
    codonstart=[]
    for i in codons:
        if codons.get(i) == "*":
            codonstop.append(i)
        if codons.get(i) =="M":
            codonstart.append(i)

    #Data doublée afin de trouver les séquences débordantes de l'ADN circulaire
    seq['data']+=seq['data']
    revseq=copy.deepcopy(seq)
    revseq=bio.complement_reverse(revseq)

    #Stock TOUS les stops dans 3 listes, une pour chaque frame. i correspond à la séquence reverse complement
    PosF=[[],[],[]]
    PosFi=[[],[],[]]
    #Stock le premier start entre 2 stops dans 3 listes aussi.
    PosD=[[],[],[]]
    PosDi=[[],[],[]]
    #Stock le stop correspondant.
    PosF2=[[],[],[]]
    PosF2i=[[],[],[]]
    #Boucle pour les stops.
    for k in range(3): #Chaque frame
        for i in range(k,len(seq['data']),3): #Tous les 3 pas
            if seq['data'][i:i+3] in codonstop:
                PosF[k].append(i)
        for i in range(k,len(revseq['data']),3): #boucle revseq
            if revseq['data'][i:i+3] in codonstop:
                PosFi[k].append(i)

    #Boucle pour les 2 autres listes.
    for j in range(3):
        for k in range(0,len(PosF[j])-1): #Avance d'un stop
            a=0 #Le compteur 'a' permet d'arreter le stockage des starts entre 2 stops apres le premier start trouvé
            for i in range(PosF[j][k],PosF[j][k+1],3): #Entre deux stops
                if seq['data'][i:i+3] in codonstart and a==0:
                    PosD[j].append(i)
                    PosF2[j].append(PosF[j][k+1])
                    a=1
        for k in range(0,len(PosFi[j])-1): #revseq
            a=0
            for i in range(PosFi[j][k],PosFi[j][k+1],3):
                if revseq['data'][i:i+3] in codonstart and a==0:
                    PosDi[j].append(i)
                    PosF2i[j].append(PosFi[j][k+1])
                    a=1
                    
    #Nouvelles variables pour stocker uniquement les ORFs dépassant le seuil.
    NewD,NewF=[[],[],[]],[[],[],[]]
    NewDi,NewFi=[[],[],[]],[[],[],[]]
    for k in range(3):
        for i in range(len(PosD[k])):
            if PosF2[k][i]-PosD[k][i] > threshold:
                NewD[k].append(PosD[k][i])
                NewF[k].append(PosF2[k][i])
        for i in range(len(PosDi[k])):
            if PosF2i[k][i]-PosDi[k][i] > threshold:
                NewDi[k].append(PosDi[k][i])
                NewFi[k].append(PosF2i[k][i])

    #Données stockées dans une liste de dictionnaire
    ORFs=[]
    ORFsi=[]
    for k in range(3):
        for i in range(len(NewD[k])):
            #remplissage du dico, start, stop, frame, longueur, proteine traduite
            ORFs.append({'start':NewD[k][i]+1,'stop':NewF[k][i]+3,'frame':k+1,'length':NewF[k][i]-NewD[k][i]+3,'protein':translate(seq['data'][NewD[k][i]:NewF[k][i]],codeTable)})
        for i in range(len(NewDi[k])):
            ORFsi.append({'start':len(seq['data'])-NewDi[k][i],'stop':len(seq['data'])-NewFi[k][i]+1,'frame':-k-1,'length':NewFi[k][i]-NewDi[k][i],'protein':translate(revseq['data'][NewDi[k][i]:NewFi[k][i]],codeTable)})

    #Suppression doublons. Boucles pour comparer les items entre eux
    for i in range(len(ORFs)):
        for j in range(len(ORFs)):
            if ORFs[i]['protein'] == ORFs[j]['protein'] and i!=j:
                ORFs[j]={'protein':''}
    for i in range(len(ORFsi)):
        for j in range(len(ORFsi)):
            if ORFsi[i]['protein'] == ORFsi[j]['protein'] and i!=j:
                ORFsi[j]={'protein':''}
                
    for i in range(len(ORFs)-1,0,-1):
        if ORFs[i]=={'protein':''}:
            del ORFs[i]
    for i in range(len(ORFsi)-1,0,-1):
        if ORFsi[i]=={'protein':''}:
            del ORFsi[i]
    
    #Correction de l'exception : si le condon start se trouve AVANT le tout premier stop et sur le début de la séquence
    for i in range(len(ORFs)):
        #Il est stocké mais n'est trouvé que dans le seq['data'] doublé.
        if ORFs[i]['start']>len(seq['data'])/2:
            #Il faut donc modifier la position start et stop
            ORFs[i]['start']=int(ORFs[i]['start']-len(seq['data'])/2)
            ORFs[i]['stop']=int(ORFs[i]['stop']-len(seq['data'])/2)

    for i in range(len(ORFsi)):
        if ORFsi[i]['stop']>len(revseq['data'])/2:
            ORFsi[i]['start']=int(ORFsi[i]['start']-len(revseq['data'])/2)
            ORFsi[i]['stop']=int(ORFsi[i]['stop']-len(revseq['data'])/2)

    ORFT=ORFs+ORFsi #Fusion liste simple et liste reverse comp

    return ORFT


def getLengths(orf_list):
        
    """Get the lengths of all the ORF in the list. 

    This function is written by Julien Benetti.

    Args:
        orf_list: list of ORF from findORF (a list of dictionnaries)

    Returns:
        The function getLengths returns a list with all the lengths of the ORFs.
    """
    
    result_list=[]
    for i in range(len(orf_list)):
        result_list.append(orf_list[i]['length'])
    return result_list

def getLongestORF(orflist):
        
    """Get the longest ORF of the list.

    This function is written by Julien Benetti.

    Args:
        orflist: list of ORF from findORF (a list of dictionnaries)

    Returns:
        The function getLongestORF returns the longestORF of the list as an integer.
    """
    
    a=getLengths(orflist)
    max=a[0]
    for i in a:
        if i>max:
            max=i
    return max

def getTopLongestORF(orflist,value):
        
    """Get the top %value of longest ORFs of the list

    This function is written by Julien Benetti.

    Args:
        orflsit: list of ORF from findORF (a list of dictionnaries)
        value: percentage value (between 0.0 and 1.0)

    Returns:
        The function getTopLongestORF returns longest ORFs of the list as a list.
    """
    
    a=getLengths(orflist)
    a.sort()
    seuil=a[round(-(len(a)*value))]
    toplist=[]
    for i in a:
        if i>=seuil:
            toplist.append(i)
    return toplist

"""sequence=bio.readFASTA('./seq.txt')

a=findORF(sequence,0,4)
b=getLengths(a)
c=getLongestORF(a)
d=getTopLongestORF(a,0.1)

for i in a:
    print(i)
print(len(a))
print(b)
print(c)
print(d)"""

##### Stats #####

def getRandomSequence(seqLength):
    """Returns a random sequence of the length seqLength

    This function is written by Binet Martin.

    Args:
        seqLength : the length of the desired sequence as an integer

    Returns:
        The function returns a dictionary containing three keys:
            - "description" : contains a description of the sequence similar to "Random sequence | 10bp"
            - "type" : contains the nature of the data, "dna"
            - "data" : contains the sequence as a list of randomly generated A, C, T and G
    """
    nucleotides = ("A", "C", "G", "T")
    seq = ""
    for i in range(seqLength):
        seq += random.choice(nucleotides)
    
    dictionary = {"description": "Random sequence | " + str(seqLength) + "bp", "type": "dna", "data" : seq}
    
    return dictionary


def shuffle(seq):
    """Returns the sequence seq after having shuffled it and permuted the symbols
    
    This function uses the Modern version of the Fisher-Yates shuffle, also known as Durstenfeld's version.
    The algorithm works as follow:
    -- To shuffle an array a of n elements (indices 0..n-1):
        for i from 0 to n−2 do
            j ← random integer such that i ≤ j < n
            exchange a[i] and a[j]
            
    This function is written by Binet Martin.

    Args:
        seq : a nucleic or protein sequence as a dictionary containing "description", "data", and "type"

    Returns:
        The function returns the initial dictionary, with the permuted sequence instead of the original.
        The "description" is also updated by adding "shuffle" at the end of it
    """

    # Transforms the string data into a list:
    listseq = list(seq["data"])
    
    # Fisher–Yates shuffle - Modern Algorithm
    
    for i in range(0, len(listseq)-2):
        n = len(listseq) - i - 1
        position = random.randint(0,n)
        listseq[n], listseq[position] = listseq[position], listseq[n]
        
    #The dictionary is modified to replace the data and add the shuffle keyword
        
    seq["description"] += " | shuffle"
    seq["data"] = str(listseq)
        
    return seq


def writeCSV(filename, separator, data):
    """Writes the data into a new file in the CSV format, using the separator argument

    This function is written by Binet Martin.

    Args:
        filename : the name of the file to be created
        separator : the separator that will be used in the file
        data : a dictionary type object

    Returns:
        This function does not return anything
    """
    
    filetowrite = open(filename, "w")
    values = []
    i = 0  #Count the number of objects already written
    for item in data:
        filetowrite.write(item)
        i += 1
        if i < len(data.keys()):
            filetowrite.write(separator)
        values.append(data[item])
    filetowrite.write("\n")
    i = 0
    for value in values:
        filetowrite.write(str(value))
        i += 1
        if i < len(values):
            filetowrite.write(separator)
    
    filetowrite.close()


def readCSV(filename, separator):
    """Read the data from a file in the CSV format, using the separator argument

    This function is written by Binet Martin.

    Args:
        filename : the name of the file to be read
        separator : the separator that will be used to parse the date from the file

    Returns:
        This function returns a dictionary where keys are strings from the first line of the file and values are from the second line of the file
    """
    
    filetoread = open(filename, "r")
    lines = []
    for line in filetoread:
        line = line.replace("\n", "").split(separator)
        lines.append(line)
    keys, values = lines[0], lines[1]
    dictionary = {}
    for i in range(0,len(keys)):
        try:
            dictionary[keys[i]] = int(values[i])
        except:
            dictionary[keys[i]] = values[i]
    return dictionary
        

def compare(orflist1, orflist2):
    """Compares two lists of ORFs/genes to find all identical genes between the two sets
    
    The function considers that two genes are identital if they share at least 99% of similarities between each other
    If they are of different lengths, the shortest sequence is saved.
    If no genes are found with 99% threshold, the comparison is ran again with a lower threshold until at least one similar ORF is found

    This function is written by Binet Martin.

    Args:
        orflist1 : the first list of ORFs
        orflist2 : the second list of ORFs

    Returns:
        This function returns a list of the identical ORFs between the two sets.
    """
    identicals = []
    threshold = 0   # in percentage
    while len(identicals) == 0:     # if no identical ORF is found, we lower the threshold
        threshold += 1
        for orf1 in orflist1:
            for orf2 in orflist2:
                if orf1 == orf2 or orf1 in orf2 or orf2 in orf1:
                    if len(orf1) > len(orf2):
                        identicals.append(orf1)
                    elif len(orf2) > len(orf1):
                        identicals.append(orf2)
                    else:
                        identicals.append(orf1)
                else:
                    same = 0
                    diff = 0
                    for i in range(0, min(len(orf1), len(orf2))):
                        if orf1[i] == orf2[i]:
                            same += 1
                        else:
                            diff += 1
                        if diff / min(len(orf1), len(orf2)) * 100 > threshold:
                            break
                    percent = same / (same+diff) * 100
                    if percent >= (100 - threshold):
                        if len(orf1) > len(orf2):
                            identicals.append(orf1)
                        elif len(orf2) > len(orf1):
                            identicals.append(orf2)
                        else:
                            identicals.append(orf1)
    print("Sequences identical at " + str(100 - threshold) + "%")
                    
    return identicals

##### GenBank Parser #####

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
