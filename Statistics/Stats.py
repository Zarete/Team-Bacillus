import random 

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
