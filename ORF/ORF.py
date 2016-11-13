# -*- coding: utf-8 -*-
#Problem 4.6.1.1
import myBio as bio

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
        bases = ['T', 'C', 'A', 'G']
        codons = [a+b+c for a in bases for b in bases for c in bases]
        amino_acids_one = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
        codon_table = dict(zip(codons, amino_acids_one))
        return codon_table
    
    if NCBI_ID==4:
        bases = ['T', 'C', 'A', 'G']
        codons = [a+b+c for a in bases for b in bases for c in bases]
        amino_acids_one = 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
        codon_table = dict(zip(codons, amino_acids_one))
        return codon_table

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
    POS_F_1=[0,0,0]
    POS_F_2=[0,0,0]
    POS_D,POS_F,a=[[],[],[]],[[],[],[]],0
    #Boucle pour les 3 frames
    for C in range(3):
        #Trouve le premier STOP
        K=C
        while seq['data'][K:K+3]!='TAA' and seq['data'][K:K+3]!='TAG' and seq['data'][K:K+3]!='TGA':
            POS_F_1[C]=K
            K += 3
        #-(len(seq['data']-C)%3 c'est pour pas etre out of range, K c'est pour commencer apres le premier STOP
        for i in range(K,len(seq['data'])-(len(seq['data'])-C)%3,3):
            #le a c'est un compteur pour avoir autant de start et stop
            if seq['data'][i:i+3]=='ATG' and a==0:
                POS_D[C].append(i)
                a=1
            if (seq['data'][i:i+3]=='TAA' or seq['data'][i:i+3]=='TAG' or seq['data'][i:i+3]=='TGA') and a==1:
                POS_F[C].append(i)
                a=0
    for C in range(3):
        #S'il y a plus de start que de stop c'est que l'ORF finit sur le début de la séquence (ADN circulaire), on rajoute le stop de départ.
        if len(POS_D[C]) > len(POS_F[C]):
            POS_F[C].append(POS_F_1[C])
    rev_seq=bio.complement_reverse(seq)
    #meme chose pour les 3 autres frames
    for C in range(3):
        K=C
        while seq['data'][K:K+3]!='TAA' and seq['data'][K:K+3]!='TAG' and seq['data'][K:K+3]!='TGA':
            POS_F_2[C]=K
            K += 3
        for i in range(K,len(rev_seq['data'])-(len(seq['data'])-C)%3,3):
            if rev_seq['data'][i:i+3]=='ATG' and a==0:
                #positions notés en négatif
                POS_D[C].append(-(len(rev_seq['data'])-i))
                a=1
            if (rev_seq['data'][i:i+3]=='TAA' or rev_seq['data'][i:i+3]=='TAG' or rev_seq['data'][i:i+3]=='TGA') and a==1:
                POS_F[C].append(-(len(rev_seq['data'])-i))
                a=0
    for C in range(3):
        if len(POS_D[C]) > len(POS_F[C]):
            POS_F[C].append(POS_F_2[C])
                
    #application du seuil, seul les orf assez long iront dans les nouvelles listes.
    LD,LF= POS_D[0]+POS_D[1]+POS_D[2], POS_F[0]+POS_F[1]+POS_F[2]
    newD,newF=[],[]
    for i in range(len(LD)):
        if LF[i]-LD[i] > threshold:
            newD.append(LD[i])
            newF.append(LF[i])
    listORFs=[]
    for i in range(len(newD)):
        #remplissage du dico, start, stop, frame, longueur, proteine traduite
        listORFs.append({'start':abs(newD[i])+1,'stop':abs(newF[i]),'frame':'','length':newF[i]-newD[i],'protein':''})
        #boucle pour completer ['frame']
        for C in range(3):
            if newF[i]%3==C and newF[i]>=0:
                listORFs[i]['frame']=C+1
            if newF[i]%3==C and newF[i]<0:
                listORFs[i]['frame']=-C-1
        #Completion de la proteine
        orf=seq['data'][newD[i]:newF[i]]
        listORFs[i]['protein']=translate(orf,codeTable)
    return listORFs


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
    seuil=a[round(-(len(a)*value))]
    toplist=[]
    for i in a:
        if i>seuil:
            toplist.append(i)
    return toplist

#sequence=bio.readFASTA('./seq.txt')

#a=findORF(sequence,180,4)
#b=getLengths(a)
#c=getLongestORF(a)
#d=getTopLongestORF(a,0.1)
