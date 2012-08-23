import sys
import os
import subprocess
import shutil

import checkfastaformat
import performBLAST
import adjlistcreation
import Leafcull

def help():
    """Prints out the usage information for the program."""
    
    print '\nUsage:'
    print '-h'
    print '\tDisplays usage information.'
    print '-i fastaFile'
    print '\t-i is used to provide the input file. fastaFile needs to be in the'
    print '\tfasta file format.'
    print '-p maxPercent'
    print '\t-p is used to provide the maximum allowable percentage sequence'
    print '\tidentity. The permissible range for maxPercent is 5 - 100.'
    print '-l minLength-maxLength'
    print '\tOPTIONAL: The permissible range of sequence lengths. The default'
    print '\tvalues are 40 for minLength and 10000 for maxLength.'
    print '-B exec'
    print '\tOPTIONAL: exec is the location of the folder which contains the'
    print '\tBLAST executables needed by the program. the default location'
    print '\tis the same location as the source files for the program.'
    print '-s'
    print '\tOPTIONAL: If this flag is set then SEG will be used prior to'
    print '\trunning BLAST. The default setting is to have SEG disabled.'
    print '-c cores'
    print '\tOPTIONAL: cores is the number of threads to create to run each BLAST'
    print '\tquery. The default number is 2. The value for cores must be an'
    print '\tinteger. If used the best result comes from setting this to the'
    print '\tnumber of cores your CPU has.'
    print '\n'
    print 'USAGE INFORMATION'
    print '================='
    print 'This program is deigned to cull a dataset of protein sequences so that no'
    print 'two sequences have a sequence identity greater than the specified threshold'
    print 'percentage. The method used is the Leaf heuristic, which is described in'
    print 'PAPER.'
    print 'The supplementary information of this paper contains further information'
    print 'about how the algorithm works.'


def main(inputFile, maxPercentage, minLength, maxLength, SEG, cullOperationID='CullResults'):
    """Runs the protein culling.
    
    @param arguments: The command line arguments minus the name of the program called.
    @type arguments: A Python list.
    
    """
    
    # Initialise the variables that will be potentially altered by the command line arguments.
    cores = 2
    
    # Get the default location for the BLAST executables.
    cwd = os.getcwd()
    srcLocation = os.path.abspath(__file__)
    srcLocation = '/'.join(srcLocation.split('/')[:-1])
    outputLocation = srcLocation + '/' + cullOperationID
    if os.path.isdir(outputLocation):
        shutil.rmtree(outputLocation)
    elif os.path.exists(outputLocation):
        os.remove(outputLocation)
    os.umask(00000)
    os.mkdir(outputLocation)
    
    if inputFile == '' or maxPercentage == 0:
        # If this occurs then the user has either not specified an input file or a percentage threshold.
        # Both of these arguments are mandatory.
        help()
        sys.exit()
    elif maxPercentage < 5 or maxPercentage > 100:
        print 'The valid range for the maximum allowable percentage sequence'
        print 'similarity is 5 - 100.'
        sys.exit()
    
    # Ensure that the fasta file input is appropriately formatted.
    fileToBLAST = outputLocation + '/InputCopy.fasta'
    inputFileToLoad = open(inputFile, 'r')
    inputFile = inputFileToLoad.read()
    inputFileToLoad.close()
    errorCode, message = checkfastaformat.main(inputFile)
    if errorCode != 0:
        print message
        sys.exit()
    writeOut = open(fileToBLAST, 'w')
    writeOut.write(message)
    writeOut.close()
    
    # Perform the BLASTing.
    similarities = performBLAST.main(fileToBLAST, fileToBLAST, cullOperationID + '/BLASTOutput', SEG, cores)
    
    # Create the sparsematrix of the protein similarity graph.
    adjList, protNames = adjlistcreation.user_seq_main(similarities, maxPercentage)
    
    # Choose which protein to remove from the similarity graph.
    if protNames == []:
        proteinsToCull = []
    else:
        proteinsToCull, proteinsToKeep = Leafcull.main(adjList, protNames)
    
    # Write out the proteins that were removed.
    writeOutRem = open(outputLocation + '/Removed.txt', 'w')
    for i in proteinsToCull:
        writeOutRem.write(i + '\n')
    writeOutRem.close()
    
    # Write out a fasta file of the proteins kept.
    writeOutKeepFasta = open(outputLocation + '/KeptFasta.fasta', 'w')
    writeOutKeepList = open(outputLocation + '/KeptList.txt', 'w')
    writeOutKeepList.write('IDs\tlength\n')
    readFasta = open(fileToBLAST, 'r')
    recording = False
    uniqueProteins = []  # Used to ensure no duplicates somehow get through.
    for line in readFasta:
        if line[0] == '>' and not line[1:-1] in proteinsToCull and not line in uniqueProteins:
            # If the line starts a new protein definition, and that protein is one of the ones to keep.
            recording = True
            uniqueProteins.append(line)
            writeOutKeepFasta.write(line)
            writeOutKeepList.write(line[1:-1])
        elif line[0] == '>':
            # If the line start a new protein definition, but the protein is not one of the ones to keep.
            recording = False
        else:
            # Otherwise the line is a protein sequence.
            if recording:
                # If we are currently working on a protein that is being kept.
                writeOutKeepFasta.write(line)
                writeOutKeepList.write('\t' + str(len(line[:-1])) + '\n')
    readFasta.close()
    writeOutKeepList.close()
    writeOutKeepFasta.close()
    
    return outputLocation