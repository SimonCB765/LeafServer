from django.core.files.base import ContentFile
from django.core.mail import send_mail

from Leaf.models import ProteinInformation
from Leaf.models import Similarity
from Leaf.models import Representative
from Leaf.models import DownloadableFiles

from threading import Thread
from time import gmtime, strftime, sleep
import shutil
import os

import traceback
import sys

import adjlistcreation
import Leafcull
import userseqcontroller

class RunCull(Thread):
    def __init__(self, request, requestType):
        Thread.__init__(self)
        os.umask(00000)
        logger = open('/srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp/ErrorLogs/ERROR.log', 'a')
        logger.write('=======================================================\n=======================================================\nEntered thread for request ' + str(request.id) + ', of type ' + str(requestType) + ', on ' + strftime('%Y/%m/%d/ at %H:%M:%S', gmtime()) + '.\n')
        logger.close()
        self.request = request
        self.requestType = requestType
        logger = open('/srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp/ErrorLogs/ERROR.log', 'a')
        logger.write('\tSTART: thread for request ' + str(self.request.id) + ', of type ' + str(self.requestType) + '.\n')
        logger.close()
        self.start()
    
    def run(self):
        if self.requestType == 'seq':
            try:
                outputLocation = userseqcontroller.main(self.request.userInput.path, self.request.sequenceIdentity, self.request.minLength, self.request.maxLength, self.request.SEG, 'UserSeq' + str(self.request.id))
                logger = open('/srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp/ErrorLogs/ERROR.log', 'a')
                logger.write('\tCULLED: thread for request ' + str(self.request.id) + ', of type ' + str(self.requestType) + '.\n')
                logger.close()
                readOutput = open(outputLocation + '/KeptList.txt', 'r')
                self.request.nonredNoSeq.save('', ContentFile(readOutput.read()))
                readOutput.close()
                readOutput = open(outputLocation + '/KeptFasta.fasta', 'r')
                self.request.nonredSeq.save('', ContentFile(readOutput.read()))
                readOutput.close()
                readOutput = open(outputLocation + '/Removed.txt', 'r')
                self.request.removed.save('', ContentFile(readOutput.read()))
                readOutput.close()
                self.request.completed = True
                self.request.save()
                logger = open('/srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp/ErrorLogs/ERROR.log', 'a')
                logger.write('\tFINISH: thread for request ' + str(self.request.id) + ', of type ' + str(self.requestType) + 'on ' + strftime('%Y/%m/%d/ at %H:%M:%S', gmtime()) + 'y.\n')
                logger.close()
            except:
                logger = open('/srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp/ErrorLogs/ERROR.log', 'a')
                logger.write('\tERROR: thread for request ' + str(self.request.id) + ', of type ' + str(self.requestType) + 'on ' + strftime('%Y/%m/%d/ at %H:%M:%S', gmtime()) + '.\n')
                excType, excValue, excTrace = sys.exc_info()
                logger.write('\t\tException type: ' + str(excType) + '\n')
                logger.write('\t\tException Value: ' + str(excValue) + '\n')
                errors = traceback.format_exception(excType, excValue, excTrace)
                for i in errors:
                    logger.write('\t\t' + i)
                logger.close()

            try:
                outputLocation == False
            except:
                pass
            else:
                try:
                    shutil.rmtree(outputLocation)
                except:
                    sleep(60)
                    shutil.rmtree(outputLocation)
        elif self.requestType == 'pdb':
            try:
                if self.request.skipNonXray:
                    skipNonXray = 'Yes'
                else:
                    skipNonXray = 'No'
                if self.request.skipAlphaCarbon:
                    skipAlphaCarbon = 'Yes'
                else:
                    skipAlphaCarbon = 'No'
                if self.request.cullByChain:
                    cullMethod = 'Chain'
                else:
                    cullMethod = 'Entry'
                if self.request.performIntraEntryCulling:
                    performIntraEntryCulling = True
                    intraEntrySequenceIdentity = self.request.intraEntrySequenceIdentity
                    intraEntryCull = 'Yes'
                    intraEntrySequenceIdentity = 'Within entry culling threshold: ' + str(self.request.intraEntrySequenceIdentity) + '\n'
                else:
                    performIntraEntryCulling = False
                    intraEntryCull = 'No'
                    intraEntrySequenceIdentity = ''
                if self.request.wholePDB:
                    wholePDB = True
                else:
                    wholePDB = False

                sequenceIdentity = self.request.sequenceIdentity
                minResolution = self.request.minResolution
                maxResolution = self.request.maxResolution
                maxRValue = self.request.maxRValue
                minLength = self.request.minLength
                maxLength = self.request.maxLength

                proteinData = DownloadableFiles.objects.filter(fileName__exact='ProteinInformation')[0].downloadFile.path
                similarityData = DownloadableFiles.objects.filter(fileName__exact='Similarity')[0].downloadFile.path
                representativeData = DownloadableFiles.objects.filter(fileName__exact='Representative')[0].downloadFile.path

                logger = open('/srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp/ErrorLogs/ERROR.log', 'a')
                logger.write('\tRequest parsed seqiden=' + str(sequenceIdentity) + ' res=' + str(minResolution) + '-' + str(maxResolution) + ' RVal=' + str(maxRValue) + '.\n')
                logger.close()
            
                if cullMethod == 'Entry':
                    # If the method of culling is 'by entry', record the entries and convert the entries to their corresponding chains.
                    if not wholePDB:
                        userInput = retVal.split('\n')
                        entriesUsed = userInput
                        potentialChains = []
                        chainsToCull = set([])
                        readProteinData = open(proteinData, 'r')
                        for i in readProteinData:
                            chunks = (i.strip()).split('\t')
                            chain = chunks[0]
                            entry = chunks[1]
                            if entry in userInput:
                                potentialChains.append(chain)
                                experimentType = chunks[2]
                                resolution = float(chunks[3])
                                rValueObs = float(chunks[4])
                                if chunks[6] == '0':
                                    alphaCarbonOnly = False
                                else:
                                    alphaCarbonOnly = True
                                sequence = chunks[11]
                                invalid = ((experimentType != 'XRAY' and skipNonXray) or
                                           (resolution < minResolution) or
                                           (resolution > maxResolution) or
                                           (rValueObs > maxRValue) or
                                           (alphaCarbonOnly and skipAlphaCarbon) or
                                           (minLength != -1 and len(sequence) < minLength) or
                                           (maxLength != -1 and len(sequence) > maxLength)
                                           )
                                if not invalid:
                                    chainsToCull.add(chain)
                    else:
                        entriesUsed = set([])
                        potentialChains = []
                        chainsToCull = set([])
                        readProteinData = open(proteinData, 'r')
                        for i in readProteinData:
                            chunks = (i.strip()).split('\t')
                            chain = chunks[0]
                            entry = chunks[1]
                            potentialChains.append(chain)
                            entriesUsed.add(entry)
                            experimentType = chunks[2]
                            resolution = float(chunks[3])
                            rValueObs = float(chunks[4])
                            if chunks[6] == '0':
                                alphaCarbonOnly = False
                            else:
                                alphaCarbonOnly = True
                            sequence = chunks[11]
                            invalid = ((experimentType != 'XRAY' and skipNonXray) or
                                       (resolution < minResolution) or
                                       (resolution > maxResolution) or
                                       (rValueObs > maxRValue) or
                                       (alphaCarbonOnly and skipAlphaCarbon) or
                                       (minLength != -1 and len(sequence) < minLength) or
                                       (maxLength != -1 and len(sequence) > maxLength)
                                       )
                            if not invalid:
                                chainsToCull.add(chain)
                        entriesUsed = list(entriesUsed)
                    entriesToCull = set([i[:4] for i in chainsToCull])
                elif cullMethod == 'Chain':
                    # If the method of culling is 'by chain', record the chains input by the user.
                    if not wholePDB:
                        userInput = retVal.split('\n')
                        potentialChains = []
                        chainsToCull = set([])
                        readProteinData = open(proteinData, 'r')
                        for i in readProteinData:
                            chunks = (i.strip()).split('\t')
                            chain = chunks[0]
                            if chain in userInput:
                                potentialChains.append(chain)
                                experimentType = chunks[2]
                                resolution = float(chunks[3])
                                rValueObs = float(chunks[4])
                                if chunks[6] == '0':
                                    alphaCarbonOnly = False
                                else:
                                    alphaCarbonOnly = True
                                sequence = chunks[11]
                                invalid = ((experimentType != 'XRAY' and skipNonXray) or
                                           (resolution < minResolution) or
                                           (resolution > maxResolution) or
                                           (rValueObs > maxRValue) or
                                           (alphaCarbonOnly and skipAlphaCarbon) or
                                           (minLength != -1 and len(sequence) < minLength) or
                                           (maxLength != -1 and len(sequence) > maxLength)
                                           )
                                if not invalid:
                                    chainsToCull.add(chain)
                        readProteinData.close()
                    else:
                        potentialChains = []
                        chainsToCull = set([])
                        readProteinData = open(proteinData, 'r')
                        for i in readProteinData:
                            chunks = (i.strip()).split('\t')
                            chain = chunks[0]
                            potentialChains.append(chain)
                            experimentType = chunks[2]
                            resolution = float(chunks[3])
                            rValueObs = float(chunks[4])
                            if chunks[6] == '0':
                                alphaCarbonOnly = False
                            else:
                                alphaCarbonOnly = True
                            sequence = chunks[11]
                            invalid = ((experimentType != 'XRAY' and skipNonXray) or
                                       (resolution < minResolution) or
                                       (resolution > maxResolution) or
                                       (rValueObs > maxRValue) or
                                       (alphaCarbonOnly and skipAlphaCarbon) or
                                       (minLength != -1 and len(sequence) < minLength) or
                                       (maxLength != -1 and len(sequence) > maxLength)
                                       )
                            if not invalid:
                                chainsToCull.add(chain)
                        readProteinData.close()

                logger = open('/srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp/ErrorLogs/ERROR.log', 'a')
                logger.write('\tInvalid chains removed. ' + str(len(chainsToCull)) + ' kept.\n')
                logger.close()

                # Determine representative chain information.
                # representatives records the non-representative to representative chain mapping for the non-representative
                # chains in the set of chains to cull.
                representatives = {}
                readRepresentativeData = open(representativeData, 'r')
                for i in readRepresentativeData:
                    chunks = (i.strip()).split('\t')
                    nonreprChain = chunks[0]
                    reprChain = chunks[1]
                    if nonreprChain in chainsToCull:
                        representatives[nonreprChain] = reprChain
                readRepresentativeData.close()
                # representativeChains records the set of representative chains that cover all the chains in the set of chains to cull.
                # This means that if a chain in chainsToCull is a representative itself then it is in representativeChains, and if
                # a chain in chainsToCull is not a representative, then its representative chain is in representativeChains.
                representativeChains = set([])
                for i in chainsToCull:
                    if representatives.has_key(i):
                        representativeChains.add(representatives[i])
                    else:
                        representativeChains.add(i)
                # representativesReverse records for each chain that represents at least one chain in chainsToCull, a set of the
                # non-representative chains in chainsToCull that it represents.
                # For example, if chainsToCull == [a, b, c, d], and a and b are non-representative chains represented by chain q,
                # then representativesReverse[q] = set([a, b]).
                representativesReverse = {}
                for i in representatives.keys():
                    reprChain = representatives[i]
                    if representativesReverse.has_key(reprChain):
                        representativesReverse[reprChain].add(i)
                    else:
                        representativesReverse[reprChain] = set([i])
                representativesReverseKeys = representativesReverse.keys()

                logger = open('/srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp/ErrorLogs/ERROR.log', 'a')
                logger.write('\t' + str(len(representativeChains)) + ' representatives found.\n')
                logger.close()

                if cullMethod == 'Entry':
                    removedInput = cull_main(similarityData, sequenceIdentity, representativeChains, 'entry', representativesReverse)
                    keptInput = set([i[:4] for i in entriesToCull if i[:4] not in removedInput])
                    if performIntraEntryCulling and intraEntrySequenceIdentity < 100:
                        entryToChain = {}
                        chainsOfInterest = set([])
                        readProteinData = open(proteinData, 'r')
                        for i in readProteinData:
                            chunks = (i.strip()).split('\t')
                            chain = chunks[0]
                            entry = chunks[1]
                            sequence = chunks[11]
                            invalid = (minLength != -1 and len(sequence) < minLength) or (maxLength != -1 and len(sequence) > maxLength)
                            if entry in keptInput and not invalid:
                                chainsOfInterest.add(chain)
                                if entryToChain.has_key(entry):
                                    entryToChain[entry].append(chain)
                                else:
                                    entryToChain[entry] = [chain]
                        readProteinData.close()
                        
                        representatives = {}
                        readRepresentativeData = open(representativeData, 'r')
                        for i in readRepresentativeData:
                            chunks = (i.strip()).split('\t')
                            nonreprChain = chunks[0]
                            reprChain = chunks[1]
                            if nonreprChain in chainsOfInterest:
                                representatives[nonreprChain] = reprChain
                        readRepresentativeData.close()
                        representativeChains = set([])
                        for i in chainsOfInterest:
                            if representatives.has_key(i):
                                representativeChains.add(representatives[i])
                            else:
                                representativeChains.add(i)
                        representativesReverse = {}
                        for i in representatives.keys():
                            reprChain = representatives[i]
                            if representativesReverse.has_key(reprChain):
                                representativesReverse[reprChain].add(i)
                            else:
                                representativesReverse[reprChain] = set([i])
                        
                        entryToRepChain = dict([(i, set([])) for i in entryToChain.keys()])
                        for i in chainsOfInterest:
                            entry = i[:4]
                            if representatives.has_key(i):
                                entryToRepChain[entry].add(representatives[i])
                            else:
                                entryToRepChain[entry].add(i)
                        
                        keptInputChains = set([])
                        
                        for i in keptInput:
                            if len(entryToRepChain[i]) == 1:
                                # If the entry's chains are all representated by one chain, then all the chains are identical. A random chain from the entry should be kept.
                                keptInputChains.add(entryToChain[i][0])
                                del entryToRepChain[i]
                        
                        adjList, namesList = adjlistcreation.intra_entry_main(similarityData, intraEntrySequenceIdentity, representativeChains, entryToRepChain)
                        
                        for i in range(len(adjList)):
                            chainsToCull = Leafcull.main(adjList[i], namesList[i])
                            keptReprChains = [j for j in namesList[i] if not j in chainsToCull]
                            for j in keptReprChains:
                                # Calculate the kept input chains.
                                if representativesReverse.has_key(i):
                                    # If the representative chain that was kept has non-representative chains in the input, then select one of them.
                                    keptInputChains.add(iter(representativesReverse[j]).next())
                                else:
                                    # If the representative chain that was kept has no non-representative chains in the input, then the representative chain was in the input. Keep it.
                                    keptInputChains.add(j)
                    else:
                        keptInputChains = set([i for i in chainsToCull if i[:4] in keptInput])
                else:
                    removedReprChains = set(cull_main(similarityData, sequenceIdentity, representativeChains, 'chain', {}))
                    keptReprChains = [i for i in representativeChains if i not in removedReprChains]
                    keptInputChains = set([])
                    for i in keptReprChains:
                        # Calculate the kept input chains.
                        if representativesReverse.has_key(i):
                            # If the representative chain that was kept has non-representative chains in the input, then select one of them.
                            keptInputChains.add(iter(representativesReverse[i]).next())
                        else:
                            # If the representative chain that was kept has no non-representative chains in the input, then the representative chain was in the input. Keep it.
                            keptInputChains.add(i)
                    removedInput = sorted([i for i in potentialChains if i not in keptInputChains])

                logger = open('/srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp/ErrorLogs/ERROR.log', 'a')
                logger.write('\tCULLED: thread for request ' + str(self.request.id) + ', of type ' + str(self.requestType) + ', on ' + strftime('%Y/%m/%d/ at %H:%M:%S', gmtime()) + '.\n')
                logger.close()
                
                keptInputOutput = 'IDs\tlength\tExptl.\tresolution\tR-factor\tFreeRvalue\n'
                fastaOutput = ''
                entryStats = {}
                readProteinData = open(proteinData, 'r')
                for i in readProteinData:
                    chunks = (i.strip()).split('\t')
                    chain = chunks[0]
                    entry = chunks[1]
                    experimentType = chunks[2]
                    resolution = chunks[3]
                    rValueObs = chunks[4]
                    rValueFree = chunks[5]
                    if chunks[6] == '0':
                        alphaCarbon = 'no'
                    else:
                        alphaCarbon = 'yes'
                    description = chunks[7]
                    dbName = chunks[8]
                    dbCode = chunks[9]
                    organism = chunks[10]
                    sequence = chunks[11]
                    if cullMethod == 'Entry':
                        if entry in keptInput:
                            entryStats[entry] = {'len' : str(len(sequence)), 'expt' : experimentType, 'res' : resolution, 'rval' : rValueObs, 'freeRval' : rValueFree}
                        if chain in keptInputChains:
                            fastaOutput += '>' + '\t'.join([chain, str(len(sequence)), experimentType, resolution, rValueObs, rValueFree, alphaCarbon, description, '<' + dbName + ' ' + dbCode + '>', '[' + organism + ']']) + '\n' + sequence + '\n'
                    else:
                        if chain in keptInputChains:
                            keptInputOutput += '\t'.join([chain, str(len(sequence)), experimentType, resolution, rValueObs, rValueFree]) + '\n'
                            fastaOutput += '>' + '\t'.join([chain, str(len(sequence)), experimentType, resolution, rValueObs, rValueFree, alphaCarbon, description, '<' + dbName + ' ' + dbCode + '>', '[' + organism + ']']) + '\n' + sequence + '\n'
                readProteinData.close()
                
                if cullMethod == 'Entry':
                    keptInputOutput += '\n'.join(['\t'.join([i, entryStats[i]['len'], entryStats[i]['expt'], entryStats[i]['res'], entryStats[i]['rval'], entryStats[i]['freeRval']]) for i in sorted(entryStats.keys())])

                self.request.removed.save('', ContentFile('\n'.join(removedInput)))
                self.request.nonredNoSeq.save('', ContentFile(keptInputOutput))
                self.request.nonredSeq.save('', ContentFile(fastaOutput))
                self.request.completed = True
                self.request.save()
            except:
                logger = open('/srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp/ErrorLogs/ERROR.log', 'a')
                logger.write('\tERROR: thread for request ' + str(self.request.id) + ', of type ' + str(self.requestType) + 'on ' + strftime('%Y/%m/%d/ at %H:%M:%S', gmtime()) + '.\n')
                excType, excValue, excTrace = sys.exc_info()
                logger.write('\t\tException type: ' + str(excType) + '\n')
                logger.write('\t\tException Value: ' + str(excValue) + '\n')
                errors = traceback.format_exception(excType, excValue, excTrace)
                for i in errors:
                    logger.write('\t\t' + i)
                logger.close()
        else:
            logger = open('/srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp/ErrorLogs/ERROR.log', 'a')
            logger.write('\tERROR: Request not for sequence or PDB culling!!!\n')
            logger.close()


def cull_main(similarities, thresholdPercentage, representativeChains, adjType, representativesReverse={}):
    """Perform the PDB redundancy removal.

    @param similarities: A record of the percentage sequence identity between the chains/entries up for culling.
    @type similarities : dictionary or file name (string)
    @param thresholdPercentage: The maximum permissible percentage sequence identity that any two chains/entries may possess.
    @type thresholdPercentage : float
    
    """

    # Create the sparsematrix of the protein similarity graph.
    if adjType == 'chain':
        adjacent, proteinNames = adjlistcreation.pdb_chain_main(similarities, thresholdPercentage, representativeChains)
    elif adjType == 'entry':
        adjacent, proteinNames = adjlistcreation.pdb_entry_main(similarities, thresholdPercentage, representativeChains, representativesReverse)
    
    # Choose which proteins to remove from the similarity graph.
    if proteinNames == []:
        # This is True if there are no similarities greater than the given percentage sequence identity. If there are no
        # chains that are too similar, then there is no need to cull any chains from the network.
        proteinsToCull = []
    else:
        # Choose which chains to remove from the similarity graph.
        proteinsToCull, proteinsToKeep = Leafcull.main(adjacent, proteinNames)

    return proteinsToCull