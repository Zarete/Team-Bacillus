# -*- coding: utf-8 -*-
#Problem 4.6.1.1
import myBio as bio
import copy

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

sequence=bio.readFASTA('./seq.txt')

a=findORF(sequence,0,4)
b=getLengths(a)
c=getLongestORF(a)
d=getTopLongestORF(a,0.1)

for i in a:
    print(i)
print(len(a))
print(b)
print(c)
print(d)
