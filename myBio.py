# -*- coding: utf-8 -*-

# Thomas Maucourt
# Bioinformatics Tools

# Modify code to use dict in functions


import copy

def version():
    '''Returns the version number of the myBio module'''
    return '3.3.1'

def author():
    return 'Thomas MAUCOURT'

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

################### F U N C T I O N S ###################

amino_acids_one = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L',
                   'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'X']

amino_acids_three = ['Arg', 'His', 'Lys', 'Asp', 'Glu', 'Ser', 'Thr', 'Asn', 'Gln', 'Cys', 'Sec', 'Gly', 'Pro', 'Ala', 'Val', 'Ile', 'Leu', 'Met', 'Phe', 'Tyr', 'Trp']

bases = ['T', 'G', 'C', 'I', 'A', 'X','Q', 'R', 'Y', 'N', 'U']

amino_acids_props = {'tiny' : ['a', 'c', 'g', 's', 't'],
                    'small' : ['a', 'b', 'd', 'g', 'n', 'p', 's', 't', 'v'],
                    'aliphatic' : ['a','i', 'l', 'v'],
                    'aromatic' : ['f', 'h', 'w', 'y'],
                    'non-polar' : ['a', 'c', 'f', 'g', 'i', 'l', 'm', 'p', 'v', 'w', 'y'],
                    'polar' : ['d', 'e', 'h', 'k', 'n', 'q', 'r', 's', 't', 'z'],
                    'charged' : ['b', 'd', 'e', 'h', 'k', 'r', 'z'],
                    'basic' : ['h', 'k', 'r'],
                    'acidic' : ['b', 'd', 'e', 'z']}

# Parse FASTA sequence and extract data
def fromFASTA(txt, mol_type = 'auto'):
    seq = {'description' : 'mybiosequence', 'data':'None', 'type':mol_type, 'AC': '', 'DB' : '', 'ID' : ''}

    pos1 = txt.index('>')
    pos2 = txt.index('\n', pos1)

    seq['description'] = txt[pos1+1:pos2]
    seq['data'] = txt[pos2:].replace('\n', '').upper()

    # Extraction données de la description
    header = seq['description'].split('|')
    if len(header) == 3:
        seq['DB'] = header[0]
        seq['AC'] = header[1]
        seq['ID'] = header[2]

    if mol_type == 'auto':
        if isDNA(seq):
            mol_type = 'dna'
        if isRNA(seq):
            mol_type = 'rna'
        if isProtein(seq):
            mol_type = 'protein'
        else:
            mol_type = 'none'
    elif mol_type != 'auto':
        moltype = mol_type

    else:
        mol_type = 'none'

    seq['type'] = mol_type

    return seq

# Return length of the given sequence
def length(seq):
    return len(seq['data'])

# Check if the given sequence is DNA
def isDNA(seq):
    return 'M' not in seq['data'] and 'U' not in seq['data'] and 'P' not in seq['data'] and 'V' not in seq['data']

# Check if the given sequence is RNA
def isRNA(seq):
    return 'U' in seq['data'] and 'M' not in seq['data']

# Check if the given sequence is a nucleic acid sequence
def isNucleic(seq):
    return isDNA(seq) or isRNA(seq)

# Check if the given sequence is a protein
def isProtein(seq):
    return not isNucleic(seq)


# Convert the given text into a FASTA sequence
def fromString(txt):
    seq = ''
    dictionnary = {}
    for carac in txt:
        if carac.upper() in amino_acids_one or carac in bases or carac == '*':
            seq += carac

    dictionnary['description'] = ">Mybio Sequence"
    dictionnary['data'] = seq.upper()

    return dictionnary

# Display the Genetic Code table
def getStandardCode():
    bases = ['T', 'C', 'A', 'G']
    codons = [a+b+c for a in bases for b in bases for c in bases]
    amino_acids_one = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    codon_table = dict(zip(codons, amino_acids_one))
    return codon_table

# Transcribe DNA into RNA
def transcribe(dna_seq):
    if dna_seq['type'] == 'dna':
        dna_seq['data'] = dna_seq['data'].replace("T", "U")
        dna_seq['type'] = 'rna'
        return dna_seq
    else:
        return None

def countAll(seq):
    dictionnary = {}
    for elem in seq['data']:
        if elem in dictionnary.keys():
            dictionnary[elem] += 1
        else:
            dictionnary[elem] = 1
    return dictionnary

# Return the complement of the given sequence
def complement(seq):
    compl = {}
    result = ''
    co_bases = {'A' : 'T', 'G' : 'C', 'C' : 'G', 'T' : 'A', 'R' : 'Y', 'Y' : 'R'}
    for letter in seq['data']:
        if letter in co_bases.keys():
            result += co_bases[letter]
        else:
            result += letter
    compl['data'] = result
    return compl

# Return the reverse of the given sequence
def reverse(seq):
    result = {}
    reverse_data = seq['data'][::-1]
    result['data'] = reverse_data
    return result

# Return the reverse complement of the given sequence
def complement_reverse(seq):
    rev_compl = {}
    rev_compl['data'] = seq['data']
    rev_compl['description'] = seq['description']
    rev_compl['ID'] = seq['ID']
    rev_compl['DB'] = seq['DB']
    rev_compl['AC'] = seq['AC']
    rev_compl['type'] = seq['type']
    rev_compl['data'] = complement(seq)['data'][::-1]

    return rev_compl

# Return sequences of wlength size without overlapping
def getWords(seq, wlength):
    words = []
    for i in range(0,len(seq['data']), wlength):
        words.append(seq['data'][i:i+wlength])
    return words

# Return sequences of wlength size with overlapping
def getOverlapWords(seq, wlength):
    words_bis = []
    for i in range(0, len(seq['data'])-(wlength-1), 1):
        words_bis.append(seq['data'][i:i+wlength])
    return words_bis

# Read FASTA files and return the datas fromFASTA.
def readFASTA(filename):
    with open(filename, 'r') as fichier:
        fasta = fichier.read()
    return fromFASTA(fasta)

# Create a sequence in FASTA and write it in a file
def writeFASTA(seq, filename):
    sequence = fromFASTA(seq)
    with open(filename, 'w') as fichier:
        fichier.write('>'+sequence['description']+'\n')
        fichier.write(sequence['data'])

# Take a string sequence and return its description and its sequence in 1 letter aa
def threeToOne(str_seq):
    returned_seq = ''
    seq = {}
    str_seq = str_seq.replace("\n", "")
    for i in range(0,len(str_seq),3):
        returned_seq += amino_acids_one[amino_acids_three.index(str_seq[i:i+3].lower().capitalize())]
    seq =fromString(returned_seq)
    seq['data'] = returned_seq.upper()
    return seq

def oneToThree(str_seq):
    returned_seq = ''
    seq = {}
    str_seq = str_seq.replace('\n', '')
    for elem in str_seq:
        if elem in amino_acids_one:
            returned_seq += amino_acids_three[amino_acids_one.index(elem)]
        else:
            returned_seq += '*'
    return returned_seq

# Return datas about the given sequence : size, composition
def compseq(seq, wlength):

    seq['words'] = {}
    seq['word%'] = {}
    word_seq = getOverlapWords(seq, wlength)

    seq['wordcount'] = len(word_seq)

    for item in word_seq:
        seq['words'][item.upper()] = word_seq.count(item)
        seq['word%'][item.upper()] = seq['words'][item.upper()] / seq['wordcount']

    return seq

# Return stats about the given sequence
def pepstats(prot_seq):
    prot_seq['props'] = {x : 0 for x in amino_acids_props.keys()}
    prot_seq['prop%'] = {}
    prot_seq['residues'] = len(prot_seq['data'])

    for item in amino_acids_one:
        prot_seq[item.upper()] = prot_seq['data'].count(item.upper())

    for item in prot_seq['data']:
        for group, aa in amino_acids_props.items():
            if item.lower() in aa:
                prot_seq['props'][group] += 1
    for cle, valeur in prot_seq['props'].items():
        prot_seq['prop%'][cle] = valeur / len(prot_seq['data']) * 100
    return prot_seq

# basic translation
def _translate(seq, codon_table = getStandardCode()):
    protein = ''
    dic = {'data' : seq}
    sequence = getWords(dic, 3)
    for elem in sequence:
        for cd, aa in codon_table.items():
            if elem.upper() == cd:
                protein += aa
    return protein

# Translation and frames
def translate(nuc_seq, frame, codon_table = getStandardCode()):
    if isRNA(nuc_seq):
        nuc_seq['data'].replace('U', 'T')

    reverse_seq = complement_reverse(nuc_seq)

    result = {'frame1': {'data' : ''}, 'frame2': {'data' : ''}, 'frame3': {'data' : ''},
              'frame-1': {'data' : ''}, 'frame-2': {'data' : ''}, 'frame-3': {'data' : ''}}

    if frame == 1:
        result['frame1']['data'] = _translate(nuc_seq['data'])
    elif frame == 2:
        result['frame2']['data'] = _translate(nuc_seq['data'][1:])
    elif frame == 3:
        result['frame3']['data'] = _translate(nuc_seq['data'][2:])
    elif frame == -1:
        result['frame-1']['data'] = _translate(reverse_seq['data'])
    elif frame == -2:
        result['frame-2']['data'] = _translate(reverse_seq['data'][1:])
    elif frame == -3:
        result['frame-3']['data'] = _translate(reverse_seq['data'][2:])
    elif frame == 6:
        result['frame1']['data'] = _translate(nuc_seq['data'])
        result['frame2']['data'] = _translate(nuc_seq['data'][1:])
        result['frame3']['data'] = _translate(nuc_seq['data'][2:])
    elif frame == -6:
        result['frame-1']['data'] = _translate(reverse_seq['data'])
        result['frame-2']['data'] = _translate(reverse_seq['data'][1:])
        result['frame-3']['data'] = _translate(reverse_seq['data'][2:])
    elif frame == 12:
        result['frame1']['data'] = _translate(nuc_seq['data'])
        result['frame2']['data'] = _translate(nuc_seq['data'][1:])
        result['frame3']['data'] = _translate(nuc_seq['data'][2:])
        result['frame-1']['data'] = _translate(reverse_seq['data'])
        result['frame-2']['data'] = _translate(reverse_seq['data'][1:])
        result['frame-3']['data'] = _translate(reverse_seq['data'][2:])
    else:
        return None

    return result
"""
# Return the different sequences given in multiFASTA format in a dictionnary
def fromMultiFASTA(txt, moltype):
    sequences = txt.split(">")
    pass

# Load sequences from a multiFASTA file
def readMultiFASTA(filename, moltype):
    sequences = ''
    file_content = []
    isolated_seq = []

    with open(filename, 'r') as fichier:
        file_content = fichier.readlines()
    for elem in file_content:
        sequences += elem.replace("\n", "")
    isolated_seq = sequences.split('>')
    results = []

    for elem in sequences:
        elem = '>' + elem + '\n'
        results.append(elem)
    return results
"""
