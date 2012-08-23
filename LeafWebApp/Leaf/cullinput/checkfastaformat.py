'''
Created on 2 Feb 2011

@author: Simon Bull
'''

def main(fileToCheck, minLength=-1, maxLength=-1):
    """Determines whether fileToCheck is appropriately formatted as a fasta file.
    
    An appropriately formatted fasta file is returned in correctFormatting.
    
    Fasta files are accepted if they have the format:
    >PID1
    letters
    >PID2
    letters
    
    Where PID1 and PID2 can be anything, and letters are a (possibly multiline) sequence of alphabetic letters.
    The letters can be upper or lower case, and each letter is interpreted as one amino acid.
    The correctly formatted fasta file is returned with the sqeuence only going over one line, and all letters
    in upper case.
    
    If a protein (i.e. a fasta information line) appears in the file more than one time, then the final appearance is
    taken to be the correct one. Prior appearances are discarded.
    
    @param fileToCheck: The input file to check for appropriate FASTA formatting.
    @type fileToCheck: string
    @param correctFormatting: The output location for the correctly formatted FASTA file.
    @type correctFormatting: string
    
    """
    
    lineCount = 1
    proteinCount = 0
    protDescription = True
    firstLine = True
    
    proteinsInFile = {}
    
    # The first line of a well formatted FASTA file will start with a '>'.
    # Following this will be a single line of uppercase letters denoting the amino acid sequence.
    # This alternation of lines beginning with '>' and lines of all uppercase letters will continue until
    # the file ends with a line of all uppercase letters.
    checking = fileToCheck.rstrip()
    checking = checking.lstrip()
    checking = checking.split('\n')
    for line in checking:
        line = line.rstrip()
        if firstLine:
            if line[0] == '>':
                currentProt = line
                protDescription = False
                firstLine = False
                currentSeq = ''
            else:
                # If the line was not the beginning of a protein record, terminate the program.
                errorMessage = "Expected line " + str(lineCount) + " to start with a >, but instead got: " + line
                return 1, errorMessage
        elif protDescription:
            # This is true only if a line beginning with a > is expected.
            if line[0] == '>':
                if minLength == -1:
                    if maxLength == -1:
                        proteinsInFile[currentProt] = currentSeq
                    elif len(currentSeq) <= maxLength:
                        proteinsInFile[currentProt] = currentSeq
                elif len(currentSeq) >= minLength:
                    if maxLength == -1:
                        proteinsInFile[currentProt] = currentSeq
                    elif len(currentSeq) <= maxLength:
                        proteinsInFile[currentProt] = currentSeq
                currentProt = line
                protDescription = False
            else:
                # If the line does not begin with a >, and it is expected to, it is possible that
                # the amino acid sequence is split over multiple lines.
                if line.isalpha():
                    # Check if every character on the line is a letter.
                    # If every character is then write out the line in upper case letters.
                    currentSeq += line.upper()
                else:
                    # If the line was not formatted as an amino acid sequence, terminate the program.
                    errorMessage = "Expected line " + str(lineCount) + " to start with a >, but instead got: " + line
                    return 1, errorMessage
        else:
            # If an amino acid sequence is expected
            if line.isalpha():
                # If the line is all alphabetic characters, write the line out and indicate that we are expecting a
                # protein description line next (i.e. one beginning with a >).
                currentSeq = line.upper()
                protDescription = True
            else:
                # If the line did not contain only letters, terminate the program.
                errorMessage = "Expected line " + str(lineCount) + " to contain only letters, but instead got: " + line
                return 2, errorMessage
        
        lineCount += 1
    
    # Catch the final protein from the file.
    if minLength == -1:
        if maxLength == -1:
            proteinsInFile[currentProt] = currentSeq
        elif len(currentSeq) <= maxLength:
            proteinsInFile[currentProt] = currentSeq
    elif len(currentSeq) >= minLength:
        if maxLength == -1:
            proteinsInFile[currentProt] = currentSeq
        elif len(currentSeq) <= maxLength:
            proteinsInFile[currentProt] = currentSeq

    if len(proteinsInFile.keys()) < 2:
        # There are to few protein sequences entered
        errorMessage = ("Not enough unique protein sequences have been entered." +
                        " This is possibly caused by not enough sequences of the required length being provided."
                        )
        return 3, errorMessage
    elif len(proteinsInFile.keys()) > 500:
        # There are to many protein sequences entered
        errorMessage = "Too many protein sequences have been entered. No more than 500 sequences can be entered. In order to cull more than 500 sequences, please download the standalone version of Leaf from http://www.bioinf.manchester.ac.uk/leaf/downloads/#SourceCode."
        return 3, errorMessage
    elif protDescription:
        # Return an indication that the FASTA file is correctly formatted
        outputString = ''
        for i in proteinsInFile.keys():
            outputString += i + '\n' + proteinsInFile[i] + '\n'
        return 0, outputString[:-1]
    else:
        # The file did not end with a protein sequence
        errorMessage = "Reached the end of the file, but no protein sequence found for the final protein."
        return 3, errorMessage