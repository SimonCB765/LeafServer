'''
Created on 7 Feb 2011

@author: Simon Bull
'''

import sys
import os
import subprocess
import shutil
import time

import processPSIoutput

def main(inputFile, databaseFile, blastOperationID, SEG=False, cores=2, minAlignLength=20, maxEValue=1):
    """
    
    """
    
    # Get the location of the BLAST executables.
    cwd = os.getcwd()
    srcLocation = os.path.abspath(__file__)
    srcLocation = '/'.join(srcLocation.split('/')[:-1])
    BLASTExecutables = srcLocation + '/BLASTExecutables'
    outputLocation = srcLocation + '/' + blastOperationID
    if os.path.exists(outputLocation):
        shutil.rmtree(outputLocation)
    os.mkdir(outputLocation)
    
    # Make a BLASTable database from the database file.
    databaseDir = outputLocation + '/TempDatabase'
    os.umask(00000)
    os.makedirs(databaseDir + '/TempDB')
    makeDBArgs = BLASTExecutables + '/makeblastdb -in ' + databaseFile + ' -out ' + databaseDir + '/TempDB -dbtype prot'
    subprocess.call(makeDBArgs.split(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    
    # Loop through the input file and create a FASTA format file for each individual protein.
    proteinDir = outputLocation + '/TempProteins'
    os.mkdir(proteinDir)
    fullFASTA = open(inputFile, 'r')
    protCount = 0
    for line in fullFASTA:
        if line[0] == '>':
            # If the line starts a new protein definition
            if protCount == 0:
                # If this is the first protein definition found
                proteinWrite = open(proteinDir + '/Prot' + str(protCount) + '.fasta', 'w')
                proteinWrite.write(line)
            else:
                # If this is not the first protein definition found
                proteinWrite.close()
                proteinWrite = open(proteinDir + '/Prot' + str(protCount) + '.fasta', 'w')
                proteinWrite.write(line)
            protCount += 1
        else:
            # Otherwise the line is a protein sequence
            proteinWrite.write(line)
    
    proteinWrite.close()
    fullFASTA.close()
    
    # Use each of the individual protein FASTA files just made, as a query to BLAST against the database created.
    processedBLAST = outputLocation + '/Processed.txt'
    proteinFiles = os.listdir(proteinDir)

    for file in proteinFiles:
        sequence_BLAST(processedBLAST, proteinDir + '/' + file, databaseDir + '/TempDB', BLASTExecutables + '/psiblast',
                       SEG, cores)
    
    similarities = {}
    readProcessedBLAST = open(processedBLAST, 'r')
    for line in readProcessedBLAST:
        chunks = line.split('\t')
        key = tuple(sorted([chunks[0], chunks[1]]))
        identity = float(chunks[2])
        alignLength = int(chunks[3])
        if alignLength < minAlignLength:
            continue
        evalue = float(chunks[4])
        if evalue > maxEValue:
            continue
        if similarities.has_key(key):
            oldSimilarity = similarities[key]['Identity']
            if identity > oldSimilarity:
                similarities[key] = {'Identity' : identity, 'Length' : alignLength, 'EValue' : evalue}
        else:
            similarities[key] = {'Identity' : identity, 'Length' : alignLength, 'EValue' : evalue}
    readProcessedBLAST.close()

    try:
        shutil.rmtree(outputLocation)
    except:
        time.sleep(60)
        shutil.rmtree(outputLocation)
    
    return similarities


def sequence_BLAST(processedBLAST, inputFile, database, BLASTLoc, SEG, cores):
    """Will perform the process of BLAST -> PROCESS OUTPUT -> CULL on inputFile.
    
    @param processedBLAST: The location at which to write the output of the processing of the BLAST output.
    @type processedBLAST: string
    @param inputFile: The FASTA file which needs to be submitted to PSI-BLAST.
    @type inputFile: string
    @param database: The database to BLAST the inputFile protein against
    @param database: string
    @param BLASTLoc: The location of the PSI-BLAST executable.
    @type BLASTLoc: string
    @param SEG: Set to True to use SEG to mask low complexity regions of the query.
    @type SEG: boolean
    @param cores: The number of threads to create to run BLAST with.
    @type cores: character
    
    """ 

    # Setup the parameters for the BLASTing.
    outputLoc = '.'.join(inputFile.split('.')[:-1]) + '.tmp'
    argsPSI = ['-query', inputFile, '-out', outputLoc, '-evalue', '1', '-inclusion_ethresh', '0.0001', '-num_iterations', '3', '-gap_trigger', '18', '-num_descriptions', '10000', '-num_alignments', '10000', '-dbsize', '0', '-db', database, '-outfmt', '7 qseqid sseqid pident length evalue']
    if SEG:
        argsPSI.extend(['-seg', 'yes'])
    else:
        argsPSI.extend(['-seg', 'no'])
    # Perform the BLASTing.
    argToBLAST = [BLASTLoc] + argsPSI
    subprocess.call(argToBLAST, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    # Process and cull the BLAST output.
    processPSIoutput.main(outputLoc, processedBLAST)