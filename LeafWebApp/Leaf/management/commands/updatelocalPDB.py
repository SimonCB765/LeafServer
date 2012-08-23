import sys
import os
import time
import gzip
import re

import traceback
import sys

from django.core.management.base import BaseCommand

from Leaf.models import ProteinInformation
from Leaf.models import ChainType
from Leaf.models import AllPDBEntries
from Leaf.models import Similarity
from Leaf.models import Representative
from Leaf.models import EntryRepresentative

srcLocation = os.path.abspath(__file__)
srcLocation = '/'.join(srcLocation.split('/')[:-3])
sys.path.append(srcLocation)
import cullinput.performBLAST

SECONDSINADAY = 86400*8

class Command(BaseCommand):
    
    def handle(self, *args, **options):
        os.umask(00000)
        # Perform the rsync update.
        
        # Determine what needs updating
        try:
            self.update(args[0])
            logger = open('/srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp/ErrorLogs/UPDATE.log', 'a')
            logger.write('\tSuccessfully updated on ' + time.strftime('%Y/%m/%d/ at %H:%M:%S', time.gmtime()) + '.\n')
            logger.close()
        except:
            logger = open('/srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp/ErrorLogs/UPDATE.log', 'a')
            logger.write('\tERROR encountered on ' + time.strftime('%Y/%m/%d/ at %H:%M:%S', time.gmtime()) + ':\n')
            excType, excValue, excTrace = sys.exc_info()
            logger.write('\t\tException type: ' + str(excType) + '\n')
            logger.write('\t\tException Value: ' + str(excValue) + '\n')
            errors = traceback.format_exception(excType, excValue, excTrace)
            for i in errors:
                logger.write('\t\t' + i)
            logger.close()

    
    def update(self, mmCIFLocalPDB):
        """mmCIFLocalPDB should be the folder where the mmCIF files of the PDB are stored"""

        logger = open('/srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp/ErrorLogs/UPDATE.log', 'a')
        logger.write('=======================================================\n=======================================================\nStarting update on ' + time.strftime('%Y/%m/%d/ at %H:%M:%S', time.gmtime()) + '.\n')
        logger.close()
    
        allPDBEntries = []  # Used to record all PDB entries. This enables any entries which have been deleted to be recognised.
    
        currentTime = time.time()
    
        mmCIFDirs = os.listdir(mmCIFLocalPDB)
    
        toParse = []
        for i in mmCIFDirs:
            mmCIFFolder = mmCIFLocalPDB + '/' + i
            mmCIFFiles = os.listdir(mmCIFFolder)
            for j in mmCIFFiles:
                allPDBEntries.append(j[:4].upper())
                currentFile = mmCIFFolder + '/' + j
                lastModified = os.path.getmtime(currentFile)
                creationChangeTime = os.path.getctime(currentFile)
                if currentTime - SECONDSINADAY < lastModified:# or currentTime - SECONDSINADAY < creationChangeTime:
                    # If the file has been created/modified within the last timeframe then it needs parsing.
                    toParse.append(currentFile)
    
        typeDict = {}
        proteinDict = {}
        inEntryRepresentativeDict = {}
        tempPDBRepresentativeDict = {}
        entriesChangedInLastWeek = set([])  # A list of the entries which have been added/altered within the time period. Only entries which have a 'type' (i.e. proteins, DNA, RNA) are recorded.
        alteredAndAddedChains = {}  # A mapping of entries that have been added/altered to their chains.

        logger = open('/srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp/ErrorLogs/UPDATE.log', 'a')
        logger.write('\tStarting parsing ' + str(len(toParse)) + ' entries.\n')
        logger.close()
    
        for i in toParse:
            # Parse mmCIF file.
            entryID, entityRecords, experimentalType, resolution, rFactorObs, rFactorFree = self.parse_file(i)
    
            # Chain type information.
            for j in entityRecords.keys():
                if entityRecords[j].has_key('type'):
                    entriesChangedInLastWeek.add(entryID[0])
                    chains = [entry + chain for entry in entryID for chain in entityRecords[j]['chains']]
                    type = entityRecords[j]['type'].strip()
                    if alteredAndAddedChains.has_key(entryID[0]):
                        alteredAndAddedChains[entryID[0]].extend(chains)
                    else:
                        alteredAndAddedChains[entryID[0]] = chains
    
                    if type == 'Protein':
                        if entityRecords[j].has_key('dbCode'):
                            dbCode = entityRecords[j]['dbCode'].strip()
                        else:
                            dbCode = ''
#                        dbCode = entityRecords[j]['dbCode'].strip() if entityRecords[j].has_key('dbCode') else ''
                        if dbCode in ['?', '.']:
                            dbCode = ''
                        if entityRecords[j].has_key('dbName'):
                            dbName = entityRecords[j]['dbName'].strip()
                        else:
                            dbName = ''
#                        dbName = entityRecords[j]['dbName'].strip() if entityRecords[j].has_key('dbName') else ''
                        if dbName in ['?', '.']:
                            dbName = ''
                        if entityRecords[j].has_key('description'):
                            description = entityRecords[j]['description'].strip()
                        else:
                            description = ''
#                        description = entityRecords[j]['description'].strip() if entityRecords[j].has_key('description') else ''
                        if description in ['?', '.']:
                            description = ''
                        onlyAlphaCarbon = entityRecords[j]['onlyAlphaCarbon']
                        if entityRecords[j].has_key('scientificName'):
                            scientificName = entityRecords[j]['scientificName'].strip()
                        else:
                            scientificName = ''
#                        scientificName = entityRecords[j]['scientificName'].strip() if entityRecords[j].has_key('scientificName') else ''
                        if scientificName in ['?', '.']:
                            scientificName = ''
                        sequence = entityRecords[j]['sequence'].upper()
                        # Determine if over 50% of the amino acids are X, in which case PISCES marks the protein as being a nonprotein.
                        if sequence.count('X') / float(len(sequence)) >= 0.5:
                            type = 'NonProtein'
                        else:
                            for k in chains:
                                proteinDict[k] = {'dbCode' : dbCode, 'dbName' : dbName, 'description' : description, 'entry' : entryID[0],
                                                 'experimentalType' : experimentalType, 'onlyAlphaCarbon' : onlyAlphaCarbon, 'resolution' : resolution,
                                                 'rFactorFree' : rFactorFree, 'rFactorObs' : rFactorObs, 'scientificName' : scientificName,
                                                 'sequence' : sequence}
                            # Set within entry representatives.
                            sortedChains = sorted(chains)
                            inEntryRepresentativeDict[sortedChains[0]] = sortedChains[1:]
                            # Update entire PDB representatives.
                            if not tempPDBRepresentativeDict.has_key(sequence):
                                tempPDBRepresentativeDict[sequence] = {'repr' : sortedChains[0], 'nonrepr' : sortedChains[1:]}
                            else:
                                currentRepr = tempPDBRepresentativeDict[sequence]['repr']
                                currentRes = proteinDict[currentRepr]['resolution']
                                currentRObs = proteinDict[currentRepr]['rFactorObs']
                                currentExpType = proteinDict[currentRepr]['experimentalType']
                                potentialRepr = sortedChains[0]
                                # X-ray structures take precedence.
                                if currentExpType != 'XRAY' and experimentalType == 'XRAY':
                                    tempPDBRepresentativeDict[sequence]['nonrepr'].append(currentRepr)
                                    tempPDBRepresentativeDict[sequence]['repr'] = potentialRepr
                                    tempPDBRepresentativeDict[sequence]['nonrepr'].extend(sortedChains[1:])
                                elif currentRes > resolution:
                                    tempPDBRepresentativeDict[sequence]['nonrepr'].append(currentRepr)
                                    tempPDBRepresentativeDict[sequence]['repr'] = potentialRepr
                                    tempPDBRepresentativeDict[sequence]['nonrepr'].extend(sortedChains[1:])
                                elif currentRObs > rFactorObs:
                                    tempPDBRepresentativeDict[sequence]['nonrepr'].append(currentRepr)
                                    tempPDBRepresentativeDict[sequence]['repr'] = potentialRepr
                                    tempPDBRepresentativeDict[sequence]['nonrepr'].extend(sortedChains[1:])
                                else:
                                    tempPDBRepresentativeDict[sequence]['nonrepr'].extend(sortedChains)
    
                    # Record the type of the entity.
                    for k in chains:
                        typeDict[k] = type

        # First get all the entries that are recorded in the database
        allRecordedEntries = [i.entry for i in AllPDBEntries.objects.all()]

        #=======================================================================
        # Delete information about any entries that no longer exist.
        #=======================================================================
        # Determine entries, and therefore chains, that need removing from the database.
        toRemove = [i for i in allRecordedEntries if not i in allPDBEntries]
        if toRemove != []:
            logger = open('/srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp/ErrorLogs/UPDATE.log', 'a')
            logger.write('\tDeleting ' + str(len(toRemove)) + ' entries that no longer exist.\n')
            logger.close()
            for i in toRemove:
                # Delete the entries.
                AllPDBEntries.objects.get(entry__exact=i).delete()
                # Delete the chain type.
                ChainType.objects.filter(chain__startswith=i).delete()
                # Delete the protein information.
                ProteinInformation.objects.filter(entry__exact=i).delete()
                # Delete the within entry representative information (as the entry no longer exists).
                EntryRepresentative.objects.filter(reprChain__startswith=i).delete()
                # Delete the representative chain information where the entry's chains are non-representative.
                Representative.objects.filter(nonreprChain__startswith=i).delete()
                # Determine if the entry contains a chain which is a representative one.
                chainsRepresented = Representative.objects.filter(reprChain__startswith=i)
                uniqueRepresentativeChains = set([j.reprChain for j in chainsRepresented])
                if chainsRepresented != []:
                    # If this is True then the entry contains a chain(s) which is the representative chain for at least
                    # one other chain.
                    for uniqueChain in uniqueRepresentativeChains:
                        # For every chain of the entry being removed which uniquely represents a set of chains.
                        uniqueChainsRepresentedSet = [j for j in chainsRepresented if j.reprChain == uniqueChain]
                        nonrepresentativeChains = [j.nonreprChain for j in uniqueChainsRepresentedSet]
                        proteins = ProteinInformation.objects.filter(chain__in=nonrepresentativeChains)
                        # Find the new representative chain out of the possibilities.
                        potentialReprChains = [j for j in proteins if j.experimentType == 'XRAY']
                        if len(potentialReprChains) == 1:
                            # Found the representative chain/entry as it is the only X-ray structure.
                            newReprChain = potentialReprChains[0].chain
                            newReprEntry = potentialReprChains[0].entry
                        elif len(potentialReprChains) > 1:
                            # Determine which of the many X-ray chains is the representative.
                            minResolution = min([j.resolution for j in potentialReprChains])
                            potentialReprChains = [j for j in potentialReprChains if j.resolution == minResolution]
                            if len(potentialReprChains) == 1:
                                # Found the representative chain/entry as it has the best resolution.
                                newReprChain = potentialReprChains[0].chain
                                newReprEntry = potentialReprChains[0].entry
                            else:
                                minRValue = min([j.rValueObs for j in potentialReprChains])
                                potentialReprChains = [j for j in potentialReprChains if j.rValueObs == minRValue]
                                if len(potentialReprChains) == 1:
                                    # Found the representative chain/entry as it has the best resolution, and a better Rvalue than the other
                                    # chains/entries with that resolution.
                                    newReprChain = potentialReprChains[0].chain
                                    newReprEntry = potentialReprChains[0].entry
                                else:
                                    potentialReprChains = sorted(potentialReprChains, key=lambda x : x.chain)
                                    newReprChain = potentialReprChains[0].chain
                                    newReprEntry = potentialReprChains[0].entry
                        else:
                            # There are no X-ray chains, so determine which of the non-X-ray chains is the representative.
                            minResolution = min([j.resolution for j in proteins])
                            potentialReprChains = [j for j in proteins if j.resolution == minResolution]
                            if len(potentialReprChains) == 1:
                                # Found the representative chain/entry as it has the best resolution.
                                newReprChain = potentialReprChains[0].chain
                                newReprEntry = potentialReprChains[0].entry
                            else:
                                minRValue = min([j.rValueObs for j in potentialReprChains])
                                potentialReprChains = [j for j in potentialReprChains if j.rValueObs == minRValue]
                                if len(potentialReprChains) == 1:
                                    # Found the representative chain/entry as it has the best resolution, and a better Rvalue than the other
                                    # chains/entries with that resolution.
                                    newReprChain = potentialReprChains[0].chain
                                    newReprEntry = potentialReprChains[0].entry
                                else:
                                    potentialReprChains = sorted(potentialReprChains, key=lambda x : x.chain)
                                    newReprChain = potentialReprChains[0].chain
                                    newReprEntry = potentialReprChains[0].entry
                        # Update the representative information to contain the new representative chain/entry.
                        updates = {'reprChain' : newReprChain}
                        Representative.objects.filter(reprChain__exact=uniqueChain).update(**updates)
                        # Now there should be one entry where both the reprChain and nonreprChain are the newReprChain,
                        # so remove this entry.
                        Representative.objects.filter(nonreprChain__exact=newReprChain).filter(reprChain__exact=newReprChain).delete()
                        # Update the similarity information related to the entry. Rather than deleting it, update the
                        # information so that it includes the new representative chain that is taking over.
                        updates = {'chainA' : newReprChain, 'entryA' : newReprEntry}
                        Similarity.objects.filter(chainA__exact=uniqueChain).update(**updates)
                        updates = {'chainB' : newReprChain, 'entryB' : newReprEntry}
                        Similarity.objects.filter(chainB__exact=uniqueChain).update(**updates)
                else:
                    # If this is True then the entry does not contain any representative chains.
                    # Delete the similarity information.
                    Similarity.objects.filter(entryA__exact=i).delete()
                    Similarity.objects.filter(entryB__exact=i).delete()

        logger = open('/srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp/ErrorLogs/UPDATE.log', 'a')
        logger.write('\tAdding/updating information for ' + str(len(entriesChangedInLastWeek)) + ' entries.\n')
        logger.close()

        # Add information about new entries, and alter information about old entries.
        sequencesNewlyAdded = {}  # The chains newly added or that had their sequence altered. Indexed by the sequence, each entry is a list of chains having the sequence.
        for i in entriesChangedInLastWeek:
            if i not in allRecordedEntries:
                # If this is True then the entry is a newly added one.
                # Add the entry in the existing entries table.
                allPDBRecord = AllPDBEntries(entry=i)
                allPDBRecord.save()
                # Determine the chains that correspond to the entry.
                chains = alteredAndAddedChains[i]
                for c in chains:
                    # Add the chain type.
                    chainRecord = ChainType(chain=c, chainType=typeDict[c])
                    chainRecord.save()
                    if proteinDict.has_key(c):
                        # Add the protein information.
                        proteinInfo = proteinDict[c]
                        proteinRecord = ProteinInformation(chain=c,
                                                           entry=i,
                                                           experimentType=proteinInfo['experimentalType'],
                                                           resolution=proteinInfo['resolution'],
                                                           rValueObs=proteinInfo['rFactorObs'],
                                                           rValueFree=proteinInfo['rFactorFree'],
                                                           alphaCarbonOnly=proteinInfo['onlyAlphaCarbon'],
                                                           description=proteinInfo['description'],
                                                           dbName=proteinInfo['dbName'],
                                                           dbCode=proteinInfo['dbCode'],
                                                           organism=proteinInfo['scientificName'],
                                                           sequence=proteinInfo['sequence'],
                                                           )
                        proteinRecord.save()
                        # Add the within entry representative information.
                        if inEntryRepresentativeDict.has_key(c):
                            for nonrepr in inEntryRepresentativeDict[c]:
                                reprRecord = EntryRepresentative(nonreprChain=nonrepr, reprChain=c)
                                reprRecord.save()
                        # Add the whole PDB representative information. First find if there are any other chains with the
                        # same sequence that are not from the same entry.
                        addedChainSeq = proteinInfo['sequence']
                        storedIdenticalSeqs = ProteinInformation.objects.filter(sequence__exact=addedChainSeq).exclude(entry__exact=i)

                        if bool(storedIdenticalSeqs):
                            # If this is True then there are chains already stored that have the same sequence as the
                            # chain being added. All of these already stored chains will have one representative.

                            # Check if the current chain is a within entry representative. If it is not then ignore this
                            # chain. Ignore it as only the within entry representative is needed to perform the update.
                            if inEntryRepresentativeDict.has_key(c):
                                # Find the representative of the stored chains.
                                chainForChecking = storedIdenticalSeqs[0]  # The chain which will be used to determine the representative chain for the sequence set.
                                determineRepresentative = Representative.objects.filter(nonreprChain__exact=chainForChecking.chain)
                                if bool(determineRepresentative):
                                    # If this is True then the chain selected for checking is a non-representative chain.
                                    # Extract the representative chain.
                                    representativeChain = determineRepresentative[0].reprChain
                                    representativeChainRecord = ProteinInformation.objects.filter(chain__exact=representativeChain)
                                    representativeChainRecord = representativeChainRecord[0]
                                else:
                                    # If this is True then the chain selected for checking is the representative chain for
                                    # the sequence set.
                                    representativeChainRecord = chainForChecking
                                    representativeChain = representativeChainRecord.chain
                            
                                # Now that the representative of the stored chains with sequence equivalent to addedChainSeq has
                                # been found, it must be determined whether this chain should represent the sequence set, or if
                                # the new chain c should represent it.
                                reprExpType = representativeChainRecord.experimentType
                                reprResolution = representativeChainRecord.resolution
                                reprRValue = representativeChainRecord.rValueObs
                                newExpType = proteinInfo['experimentalType']
                                newResolution = proteinInfo['resolution']
                                newRValue = proteinInfo['rFactorObs']
                                changeRepr = False
                                if (reprExpType != 'XRAY' and newExpType != 'XRAY') or (reprExpType == 'XRAY' and newExpType == 'XRAY'):
                                    # If this is True then the representative chain can't be determined based on experimental type.
                                    if reprResolution == newResolution:
                                        # If this is True then the representative chain can't be determined based on resolution.
                                        if reprRValue == newRValue:
                                            # If this is True then the representative chain can't be determined based on R value.
                                            if cmp(c, representativeChain):
                                                # If this is True then the new chain should become the representative purely based
                                                # on the alphanumeric chain name.
                                                changeRepr = True
                                        elif reprRValue > newRValue:
                                            # If this is True then the new chain should become the representative as its R value
                                            # is better than that of the current representative (everything else being equal).
                                            changeRepr = True
                                    elif reprResolution > newResolution:
                                        # If this is True then the new chain should become the representative as its resolution
                                        # is better than that of the current representative (the experiment type being equal).
                                        changeRepr = True
                                elif reprExpType != 'XRAY' and newExpType == 'XRAY':
                                    # If this is True then the new chain should become the representative as its experimental
                                    # type is 'better' (X-ray > anything else).
                                    changeRepr = True
    
                                if changeRepr:
                                    # If this is True then one of the criteria for changing the representative chain to the
                                    # chain being added has been met.
                                    # Update the representative table so that all the chains with representativeChain as their
                                    # representative chain now have c as it.
                                    Representative.objects.filter(reprChain__exact=representativeChain).update(reprChain=c)
                                    # Add in a new record recording the fact that reprChain is now represented by c.
                                    reprRecord = Representative(nonreprChain=representativeChain, reprChain=c)
                                    reprRecord.save()
                                    # For every chain within the entry that c comes from that is represented by c,
                                    # add a record of the fact that c represents that chain.
                                    for nonrepr in inEntryRepresentativeDict[c]:
                                        reprRecord = Representative(nonreprChain=nonrepr, reprChain=c)
                                        reprRecord.save()
                                    # Update the similarity information recorded about the old representative chain
                                    # (if there was any). Rather than deleting it, update the information so that it
                                    # includes the new representative chain that is taking over.
                                    updates = {'chainA' : c, 'entryA' : i}
                                    Similarity.objects.filter(chainA__exact=representativeChain).update(**updates)
                                    updates = {'chainB' : c, 'entryB' : i}
                                    Similarity.objects.filter(chainB__exact=representativeChain).update(**updates)
                                else:
                                    # If this is True then the criteria for changing the representative chain to the
                                    # chain being added has not been met. In this case simply create new whole
                                    # PDB representative entries where the chain being added, and the chains that it
                                    # represents, are representated by the current sequence set representative.
                                    reprRecord = Representative(nonreprChain=c, reprChain=representativeChain)
                                    reprRecord.save()
                                    for nonrepr in inEntryRepresentativeDict[c]:
                                        reprRecord = Representative(nonreprChain=nonrepr, reprChain=representativeChain)
                                        reprRecord.save()
                        else:
                            # If this is False then there are no stored chains with the same sequence as the chain
                            # being added. Therefore, just add the within entry representative information for the
                            # current chain to the whole PDB representative table. This is correct as the within
                            # entry representative information is also the whole PDB representative information due
                            # the fact that there are no chains in other entries with the same sequence.
                            if inEntryRepresentativeDict.has_key(c):
                                for nonrepr in inEntryRepresentativeDict[c]:
                                    reprRecord = Representative(nonreprChain=nonrepr, reprChain=c)
                                    reprRecord.save()
                            sequencesNewlyAdded[addedChainSeq] = c
            else:
                # If this is True then the entry already exists and has been altered.
                # No need to add the entry in the existing entries table since it's already there.
                chains = alteredAndAddedChains[i]
                for c in chains:
                    if proteinDict.has_key(c):
                        proteinInfo = proteinDict[c]
                        alteredExpType = proteinInfo['experimentalType']
                        alteredResolution = proteinInfo['resolution']
                        alteredRValue = proteinInfo['rFactorObs']
                        alteredSequence = proteinInfo['sequence']
                        # Check if c is already recorded in the database, or if it is a new chain.
                        if len(ProteinInformation.objects.filter(chain__exact=c)) == 0:
                            # If this is True then the chain is not already in the database.
                            newChain = ChainType(chain=c, chainType=typeDict[c])
                            newChain.save()
                            unchangedSequence = False
                        else:
                            # If this is True then the chain is already in the database, and the record needs updating.
                            ChainType.objects.filter(chain__exact=c).update(chainType=typeDict[c])  # Alter the chain type.
                            # Alter the protein information.
                            storedProteinInfo = ProteinInformation.objects.filter(chain__exact=c)
                            storedSequence = storedProteinInfo[0].sequence
                            unchangedSequence = storedSequence == alteredSequence

                        # Delete and the re-add the record for the chain. Simpler to do this (even if there is no change)
                        # than to check every field.
                        ProteinInformation.objects.filter(chain__exact=c).delete()
                        proteinRecord = ProteinInformation(chain=c,
                                                           entry=i,
                                                           experimentType=proteinInfo['experimentalType'],
                                                           resolution=proteinInfo['resolution'],
                                                           rValueObs=proteinInfo['rFactorObs'],
                                                           rValueFree=proteinInfo['rFactorFree'],
                                                           alphaCarbonOnly=proteinInfo['onlyAlphaCarbon'],
                                                           description=proteinInfo['description'],
                                                           dbName=proteinInfo['dbName'],
                                                           dbCode=proteinInfo['dbCode'],
                                                           organism=proteinInfo['scientificName'],
                                                           sequence=proteinInfo['sequence'],
                                                           )
                        proteinRecord.save()

                        # Alter the within entry representative chain information. There are three possibilities:
                        # 1) c is not the within entry representative in the altered or stored entry information.
                        #    In this case nothing will be deleted or added.
                        # 2) c is the within entry representative in the stored information, but not the altered.
                        #    In this case the stored representative relationships will be deleted, but no new ones added.
                        # 3) c is the within entry representative in the altered information, but not the stored.
                        #    In this case no information will be deleted, but the new representative relationships added.
                        EntryRepresentative.objects.filter(reprChain__exact=c).delete()
                        if inEntryRepresentativeDict.has_key(c):
                            withinEntryChainsRepresented = inEntryRepresentativeDict[c]
                            for nonrepr in withinEntryChainsRepresented:
                                reprRecord = EntryRepresentative(nonreprChain=nonrepr, reprChain=c)
                                reprRecord.save()
                            
                            if not unchangedSequence:
                                # If this is True then the sequence of c has changed. This means that either:
                                # 1) The old representative information where c was a non-representative must be deleted.
                                # OR
                                # 2) The old representative information where c was a representative must be updated.
                                # This also applies to all the chains within the entry that c represents. I.e. they must
                                # be removed from the old representative information.
                                for nonrepr in withinEntryChainsRepresented:
                                    # Delete all within entry chains that c represents from the whole PDB table.
                                    # This is done as c's sequence has changed, and therefore c and the chains it represents will be
                                    # moving to a new representative chain set.
                                    Representative.objects.filter(nonreprChain__exact=nonrepr).delete()
                                chainsRepresentedByC = Representative.objects.filter(reprChain__exact=c)
                                if bool(chainsRepresentedByC):
                                    # Condition 2).
                                    # If this is True, then c is the representative chain for a set of other chains.
                                    # Search for the new representative chain.
                                    chainsRepresentedByC = [j.nonreprChain for j in chainsRepresentedByC]
                                    proteins = ProteinInformation.objects.filter(chain__in=chainsRepresentedByC)
                                    # Find the new representative chain out of the possibilities.
                                    potentialReprChains = [j for j in proteins if j.experimentType == 'XRAY']
                                    if len(potentialReprChains) == 1:
                                        # Found the representative chain/entry as it is the only X-ray structure.
                                        newReprChain = potentialReprChains[0].chain
                                        newReprEntry = potentialReprChains[0].entry
                                    elif len(potentialReprChains) > 1:
                                        # Determine which of the many X-ray chains is the representative.
                                        minResolution = min([j.resolution for j in potentialReprChains])
                                        potentialReprChains = [j for j in potentialReprChains if j.resolution == minResolution]
                                        if len(potentialReprChains) == 1:
                                            # Found the representative chain/entry as it has the best resolution.
                                            newReprChain = potentialReprChains[0].chain
                                            newReprEntry = potentialReprChains[0].entry
                                        else:
                                            minRValue = min([j.rValueObs for j in potentialReprChains])
                                            potentialReprChains = [j for j in potentialReprChains if j.rValueObs == minRValue]
                                            if len(potentialReprChains) == 1:
                                                # Found the representative chain/entry as it has the best resolution, and a better Rvalue than the other
                                                # chains/entries with that resolution.
                                                newReprChain = potentialReprChains[0].chain
                                                newReprEntry = potentialReprChains[0].entry
                                            else:
                                                potentialReprChains = sorted(potentialReprChains, key=lambda x : x.chain)
                                                newReprChain = potentialReprChains[0].chain
                                                newReprEntry = potentialReprChains[0].entry
                                    else:
                                        # There are no X-ray chains, so determine which of the non-X-ray chains is the representative.
                                        minResolution = min([j.resolution for j in proteins])
                                        potentialReprChains = [j for j in proteins if j.resolution == minResolution]
                                        if len(potentialReprChains) == 1:
                                            # Found the representative chain/entry as it has the best resolution.
                                            newReprChain = potentialReprChains[0].chain
                                            newReprEntry = potentialReprChains[0].entry
                                        else:
                                            minRValue = min([j.rValueObs for j in potentialReprChains])
                                            potentialReprChains = [j for j in potentialReprChains if j.rValueObs == minRValue]
                                            if len(potentialReprChains) == 1:
                                                # Found the representative chain/entry as it has the best resolution, and a better Rvalue than the other
                                                # chains/entries with that resolution.
                                                newReprChain = potentialReprChains[0].chain
                                                newReprEntry = potentialReprChains[0].entry
                                            else:
                                                potentialReprChains = sorted(potentialReprChains, key=lambda x : x.chain)
                                                newReprChain = potentialReprChains[0].chain
                                                newReprEntry = potentialReprChains[0].entry
                                    # Update the representative information to contain the new representative chain/entry.
                                    updates = {'reprChain' : newReprChain}
                                    Representative.objects.filter(reprChain__exact=c).update(**updates)
                                    # Now there should be one entry where both the reprChain and nonreprChain are the newReprChain,
                                    # so remove this entry.
                                    Representative.objects.filter(nonreprChain__exact=newReprChain).filter(reprChain__exact=newReprChain).delete()
                                    # Update the similarity information related to the entry. Rather than deleting it, update the
                                    # information so that it includes the new representative chain that is taking over.
                                    updates = {'chainA' : newReprChain, 'entryA' : newReprEntry}
                                    Similarity.objects.filter(chainA__exact=c).update(**updates)
                                    updates = {'chainB' : newReprChain, 'entryB' : newReprEntry}
                                    Similarity.objects.filter(chainB__exact=c).update(**updates)
                                else:
                                    # Condition 1)
                                    # If this is not True, then c is a non-representative chain. Delete all representative
                                    # information stored about c (information about the chains c represents has already been deleted). This is done
                                    # because the sequence has changed, and therefore the chains will be moving to a
                                    # different representative sequence set.
                                    Representative.objects.filter(nonreprChain__exact=c).delete()

                            # Alter the whole PDB representative chain information. Only bother doing this when the chain
                            # c is a within entry representative.
                            # There are a few possibilities for this:
                            # 1) c is the whole PDB representative in the stored information, and should stay that way.
                            #    Make sure the within entry representative information is in the whole PDB stored representative information.
                            # 2) c is the whole PDB representative in the stored information, and now shouldn't be.
                            # 3) c is not the whole PDB representative in the stored information, and is still not.
                            #    Just need to make sure that all of the within entry representative information is
                            #    in the stored information, and if its not then it needs adding.
                            # 4) c is not the whole PDB representative in the stored information, and now is.
                            #    Need to update the stored information to reflect the fact that c is now the
                            #    representative, and also need to make sure all the within entry representative
                            #    information is stored.
                            storedIdenticalSeqs = ProteinInformation.objects.filter(sequence__exact=alteredSequence).exclude(entry__exact=i)
                            if bool(storedIdenticalSeqs):
                                # If this is True then there are chains already stored that have the same sequence as the
                                # chain being altered. All of these already stored chains will have one representative.

                                # Find the representative of the stored chains.
                                chainForChecking = storedIdenticalSeqs[0]  # The chain which will be used to determine the representative chain for the sequence set.
                                determineRepresentative = Representative.objects.filter(nonreprChain__exact=chainForChecking.chain)
                                if bool(determineRepresentative):
                                    # If this is True then the chain selected for checking is a non-representative chain.
                                    # Extract the representative chain.
                                    storedRepresentativeChain = determineRepresentative[0].reprChain
                                    storedRepresentativeChainRecord = ProteinInformation.objects.filter(chain__exact=storedRepresentativeChain)
                                    storedRepresentativeChainRecord = storedRepresentativeChainRecord[0]
                                else:
                                    # If this is True then the chain selected for checking is the representative chain for
                                    # the sequence set.
                                    storedRepresentativeChain = chainForChecking
                                    storedRepresentativeChainRecord = chainForChecking

                                # Determine the representative chain for the chains with the same sequence as c.
                                proteins = ProteinInformation.objects.filter(sequence__exact=alteredSequence)
                                # Find the new representative chain out of the possibilities.
                                potentialReprChains = [j for j in proteins if j.experimentType == 'XRAY']
                                if len(potentialReprChains) == 1:
                                    # Found the representative chain/entry as it is the only X-ray structure.
                                    newReprChain = potentialReprChains[0].chain
                                    newReprEntry = potentialReprChains[0].entry
                                elif len(potentialReprChains) > 1:
                                    # Determine which of the many X-ray chains is the representative.
                                    minResolution = min([j.resolution for j in potentialReprChains])
                                    potentialReprChains = [j for j in potentialReprChains if j.resolution == minResolution]
                                    if len(potentialReprChains) == 1:
                                        # Found the representative chain/entry as it has the best resolution.
                                        newReprChain = potentialReprChains[0].chain
                                        newReprEntry = potentialReprChains[0].entry
                                    else:
                                        minRValue = min([j.rValueObs for j in potentialReprChains])
                                        potentialReprChains = [j for j in potentialReprChains if j.rValueObs == minRValue]
                                        if len(potentialReprChains) == 1:
                                            # Found the representative chain/entry as it has the best resolution, and a better Rvalue than the other
                                            # chains/entries with that resolution.
                                            newReprChain = potentialReprChains[0].chain
                                            newReprEntry = potentialReprChains[0].entry
                                        else:
                                            potentialReprChains = sorted(potentialReprChains, key=lambda x : x.chain)
                                            newReprChain = potentialReprChains[0].chain
                                            newReprEntry = potentialReprChains[0].entry
                                else:
                                    # There are no X-ray chains, so determine which of the non-X-ray chains is the representative.
                                    minResolution = min([j.resolution for j in proteins])
                                    potentialReprChains = [j for j in proteins if j.resolution == minResolution]
                                    if len(potentialReprChains) == 1:
                                        # Found the representative chain/entry as it has the best resolution.
                                        newReprChain = potentialReprChains[0].chain
                                        newReprEntry = potentialReprChains[0].entry
                                    else:
                                        minRValue = min([j.rValueObs for j in potentialReprChains])
                                        potentialReprChains = [j for j in potentialReprChains if j.rValueObs == minRValue]
                                        if len(potentialReprChains) == 1:
                                            # Found the representative chain/entry as it has the best resolution, and a better Rvalue than the other
                                            # chains/entries with that resolution.
                                            newReprChain = potentialReprChains[0].chain
                                            newReprEntry = potentialReprChains[0].entry
                                        else:
                                            potentialReprChains = sorted(potentialReprChains, key=lambda x : x.chain)
                                            newReprChain = potentialReprChains[0].chain
                                            newReprEntry = potentialReprChains[0].entry

                                if newReprChain == storedRepresentativeChain:
                                    # The representative chain should not be changed.
                                    chainsRepresented = Representative.objects.filter(reprChain__exact=storedRepresentativeChain).values_list('nonreprChain', flat=True)
                                    for nonrepr in withinEntryChainsRepresented:
                                        # Make sure each within entry representative relationship is stored in the whole
                                        # PDB representative table.
                                        if not nonrepr in chainsRepresented:
                                            # If chain nonrepr is not recorded as being represented by c, then add it.
                                            reprRecord = Representative(nonreprChain=nonrepr, reprChain=storedRepresentativeChain)
                                            reprRecord.save()
                                    if newReprChain == c:
                                        # The representative chain remains as c.
                                        pass
                                    else:
                                        # The representative chain was not c, and remains not c.
                                        reprRecord = Representative(nonreprChain=c, reprChain=storedRepresentativeChain)
                                        reprRecord.save()
                                else:
                                    # The representative chain should change.
                                    # Update the representative table so that all the chains with storedRepresentativeChain as their
                                    # representative chain now have newReprChain as it.
                                    Representative.objects.filter(reprChain__exact=storedRepresentativeChain).update(reprChain=newReprChain)
                                    # Add in a new record recording the fact that storedRepresentativeChain is now represented by newReprChain.
                                    reprRecord = Representative(nonreprChain=storedRepresentativeChain, reprChain=newReprChain)
                                    reprRecord.save()
                                    # Remove the possible entry where c is recorded as both the non-representative and representative
                                    # chain (can only occur if newReprChain has not had it's sequence change).
                                    Representative.objects.filter(nonreprChain__exact=newReprChain, reprChain__exact=newReprChain).delete()
                                    # Update the similarity information recorded about the old representative chain
                                    # (if there was any). Rather than deleting it, update the information so that it
                                    # includes the new representative chain that is taking over.
                                    updates = {'chainA' : newReprChain, 'entryA' : newReprEntry}
                                    Similarity.objects.filter(chainA__exact=storedRepresentativeChain).update(**updates)
                                    updates = {'chainB' : newReprChain, 'entryB' : newReprEntry}
                                    Similarity.objects.filter(chainB__exact=storedRepresentativeChain).update(**updates)
                                    if newReprChain == c:
                                        # The representative chain changes to c.
                                        # Get a list of the stored chains that are represented by newReprChain.
                                        chainsRepresented = Representative.objects.filter(reprChain__exact=c).values_list('nonreprChain', flat=True)
                                        for nonrepr in withinEntryChainsRepresented:
                                            # Make sure each within entry representative relationship is stored in the whole
                                            # PDB representative table.
                                            if not nonrepr in chainsRepresented:
                                                # If chain nonrepr is not recorded as being represented by c, then add it.
                                                reprRecord = Representative(nonreprChain=nonrepr, reprChain=c)
                                                reprRecord.save()
                                    else:
                                        # The representative chain changes to a chain which isn't c.
                                        pass
                            else:
                                # If this is False then there are no stored chins with the same sequence as c.
                                # Therefore, just add the within entry representative information for the
                                # current chain to the whole PDB representative table. This is correct as the within
                                # entry representative information is also the whole PDB representative information due
                                # the fact that there are no chains in other entries with the same sequence.
                                for nonrepr in withinEntryChainsRepresented:
                                    reprRecord = Representative(nonreprChain=nonrepr, reprChain=c)
                                    reprRecord.save()
                                if not unchangedSequence:
                                    # If the sequence has been altered, and the new sequence is unique among the stored chains, then
                                    # record the sequence as being newly added.
                                    sequencesNewlyAdded[alteredSequence] = c

        if sequencesNewlyAdded != {}:
            # If there are sequences that have been newly added then it is necessary to perform BLAST on some sequences.
            # Determine the representative chains of the newly added sequences.
            logger = open('/srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp/ErrorLogs/UPDATE.log', 'a')
            logger.write('\tBLASTing ' + str(len(sequencesNewlyAdded.keys())) + ' entries with new sequences.\n')
            logger.close()
            newlyAddedChains = []
            for i in sequencesNewlyAdded.keys():
                chainForChecking = sequencesNewlyAdded[i]
                determineRepresentative = Representative.objects.filter(nonreprChain__exact=chainForChecking)
                if bool(determineRepresentative):
                    # If this is True then the chain selected for checking is a non-representative chain.
                    # Extract the representative chain.
                    sequencesNewlyAdded[i] = determineRepresentative[0].reprChain
                    newlyAddedChains.append(determineRepresentative[0].reprChain)
                else:
                    # If this is True then the chain selected for checking is the representative chain for
                    # the sequence set.
                    sequencesNewlyAdded[i] = chainForChecking
                    newlyAddedChains.append(chainForChecking)

            # Extract all the old representative sequences.
            storedReprSeqs = set(Representative.objects.all().values_list('reprChain', flat=True))
            tempSeqChainPairs = ProteinInformation.objects.filter(chain__in=storedReprSeqs).values('chain', 'sequence')
            seqChainPairs = [i for i in tempSeqChainPairs if not i['chain'] in newlyAddedChains]
            
            # Create the fasta files.
            srcLocation = os.path.abspath(__file__)
            srcLocation = '/'.join(srcLocation.split('/')[:-1])
            newlyAdded = srcLocation + '/NewlyAdded.fasta'
            storedSeqs = srcLocation + '/StoredSeqs.fasta'
            allSeqs = srcLocation + '/AllSeqs.fasta'
            
            writeNew = open(newlyAdded, 'w')
            writeStored = open(storedSeqs, 'w')
            writeAll = open(allSeqs, 'w')
            for i in sequencesNewlyAdded.keys():
                writeNew.write('>' + sequencesNewlyAdded[i] + '\n' + i + '\n')
                writeAll.write('>' + sequencesNewlyAdded[i] + '\n' + i + '\n')
            for i in seqChainPairs:
                writeStored.write('>' + i['chain'] + '\n' + i['sequence'] + '\n')
                writeAll.write('>' + i['chain'] + '\n' + i['sequence'] + '\n')
            writeAll.close()
            writeStored.close()
            writeNew.close()
            
            similarities = cullinput.performBLAST.main(newlyAdded, allSeqs, 'NewAgainstAll')
            if seqChainPairs != []:
                # If this is True then there are old sequences in the database (will be True every time except the first).
                storedAndNewSimilarities = cullinput.performBLAST.main(storedSeqs, newlyAdded, 'StoredAgainstNew')
                for i in storedAndNewSimilarities.keys():
                    if similarities.has_key(i):
                        if storedAndNewSimilarities[i]['Identity'] > similarities[i]['Identity']:
                            similarities[i] = storedAndNewSimilarities[i]
                    else:
                        similarities[i] = storedAndNewSimilarities[i]
            
            # Save the sequence identity information in the database.
            for i in similarities.keys():
                if i[0] == i[1] or similarities[i]['EValue'] >= 1 or similarities[i]['Length'] <= 20:
                    # Ignore any similarities where the two chains are the same, the EValue >= 1 or the alignment length <= 20.
                    continue
                similarityRecord = Similarity(chainA=i[0], entryA=i[0][:-1], chainB=i[1], entryB=i[1][:-1],
                                              similarity=similarities[i]['Identity'], matchLength=similarities[i]['Length'])
                similarityRecord.save()
            
            os.remove(newlyAdded)
            os.remove(storedSeqs)
            os.remove(allSeqs)

    
    def parse_file(self, fileToParse):
        parsedData, errorMessage = self.extract_info(fileToParse, ['_entry.id', '_entity', '_entity_poly', '_exptl.method', '_atom_site.label_atom_id', '_atom_site.label_entity_id', '_refine', '_reflns', '_struct_ref', '_entity_src_gen', '_entity_src_nat', '_pdbx_entity_src_syn'])
        entryID = parsedData['_entry']['id']
        entityRecords = {}
    
        # Entity ID and description information.
        for i in range(len(parsedData['_entity']['id'])):
            entityID = parsedData['_entity']['id'][i]
            if parsedData['_entity'].has_key('pdbx_description'):
                description = parsedData['_entity']['pdbx_description'][i]
            else:
                description = ''
            entityRecords[entityID] = {'description' : description}
    
        # Entity polymer information.
        if parsedData.has_key('_entity_poly'):
            for i in range(len(parsedData['_entity_poly']['entity_id'])):
                entityID = parsedData['_entity_poly']['entity_id'][i]
                chainIDs = (parsedData['_entity_poly']['pdbx_strand_id'][i].replace(' ', '')).split(',')
                entityRecords[entityID]['chains'] = chainIDs
                polyType = parsedData['_entity_poly']['type'][i]
                if polyType in ['polypeptide(L)', 'polypeptide(D)']:
                    entityRecords[entityID]['type'] = 'Protein'
                    entityRecords[entityID]['sequence'] = parsedData['_entity_poly']['pdbx_seq_one_letter_code_can'][i]
                else:
                    type = []
                    if 'polydeoxyribonucleotide' in polyType:
                        type.append('DNA')
                    if 'polyribonucleotide' in polyType:
                        type.append('RNA')
                    if 'polysaccharide' in polyType:
                        type.append('polysaccharide')
                    if type == []:
                        entityRecords[entityID]['type'] = 'other'
                    else:
                        entityRecords[entityID]['type'] = '/'.join(type)
    
        # Entity source information.
        if parsedData.has_key('_entity_src_nat'):
            for i in range(len(parsedData['_entity_src_nat']['entity_id'])):
                entityID = parsedData['_entity_src_nat']['entity_id'][i]
                scientificName = parsedData['_entity_src_nat']['pdbx_organism_scientific'][i]
                entityRecords[entityID]['scientificName'] = scientificName
        elif parsedData.has_key('_entity_src_gen'):
            for i in range(len(parsedData['_entity_src_gen']['entity_id'])):
                entityID = parsedData['_entity_src_gen']['entity_id'][i]
                scientificName = parsedData['_entity_src_gen']['pdbx_gene_src_scientific_name'][i]
                entityRecords[entityID]['scientificName'] = scientificName
        elif parsedData.has_key('_pdbx_entity_src_syn'):
            for i in range(len(parsedData['_pdbx_entity_src_syn']['entity_id'])):
                entityID = parsedData['_pdbx_entity_src_syn']['entity_id'][i]
                scientificName = parsedData['_pdbx_entity_src_syn']['organism_scientific'][i]
                entityRecords[entityID]['scientificName'] = scientificName
    
        # External reference (UniProt) information.
        if parsedData.has_key('_struct_ref'):
            for i in range(len(parsedData['_struct_ref']['entity_id'])):
                entityID = parsedData['_struct_ref']['entity_id'][i]
                if parsedData['_struct_ref'].has_key('db_name'):
                    dbName = parsedData['_struct_ref']['db_name'][i]
                else:
                    dbName = ''
                if dbName != '':
                    entityRecords[entityID]['dbName'] = dbName
                    if parsedData['_struct_ref'].has_key('db_code'):
                        dbCode = parsedData['_struct_ref']['db_code'][i]
                    else:
                        dbCode = ''
                    entityRecords[entityID]['dbCode'] = dbCode
                else:
                    entityRecords[entityID]['dbName'] = ''
                    entityRecords[entityID]['dbCode'] = ''
    
        # Experimental type.
        expTypes = {'NEUTRON DIFFRACTION' : 'NEUTRON', 'FIBER DIFFRACTION' : 'FIBER', 'X-RAY DIFFRACTION' : 'XRAY', 'ELECTRON MICROSCOPY' : 'EM',
                    'FLUORESCENCE TRANSFER' : 'NA', 'POWDER DIFFRACTION' : 'POWDER', 'SOLUTION NMR' : 'NMR', 'INFRARED SPECTROSCOPY' : 'FTIR',
                    'ELECTRON CRYSTALLOGRAPHY' : 'NA', 'SOLUTION SCATTERING' : 'NA', 'SOLID-STATE NMR' : 'NMR'}
        if parsedData.has_key('_exptl'):
            if expTypes.has_key(parsedData['_exptl']['method'][0]):
                experimentalType = expTypes[parsedData['_exptl']['method'][0]]
            else:
                experimentalType = 'NA'
        else:
            experimentalType = 'NA'
    
        # Resolution and R Factor Information
        if parsedData.has_key('_refine'):
            if parsedData['_refine'].has_key('ls_d_res_high'):
                try:
                    resolution = float(parsedData['_refine']['ls_d_res_high'][0])
                except:
                    resolution = 100
            else:
                resolution = 100
            if parsedData['_refine'].has_key('ls_R_factor_obs'):
                try:
                    rFactorObs = float(parsedData['_refine']['ls_R_factor_obs'][0])
                except:
                    rFactorObs = 1
            else:
                rFactorObs = 1
            if parsedData['_refine'].has_key('ls_R_factor_R_free'):
                try:
                    rFactorFree = float(parsedData['_refine']['ls_R_factor_R_free'][0])
                except:
                    rFactorFree = 1
            else:
                rFactorFree = 1
        elif parsedData.has_key('_reflns'):
            if parsedData['_reflns'].has_key('d_resolution_high'):
                try:
                    resolution = float(parsedData['_reflns']['d_resolution_high'][0])
                except:
                    resolution = 100
            else:
                resolution = 100
            rFactorObs = 1
            rFactorFree = 1
        else:
            resolution = 100
            rFactorObs = 1
            rFactorFree = 1
    
        # Determine structures with only alpha carbon atoms.
        onlyAlphaCarbon = dict([(i, True) for i in parsedData['_entity']['id']])
        for i in range(len(parsedData['_atom_site']['label_entity_id'])):
            if parsedData['_atom_site']['label_atom_id'][i] != 'CA':
                onlyAlphaCarbon[parsedData['_atom_site']['label_entity_id'][i]] = False
        for i in onlyAlphaCarbon.keys():
            if entityRecords[i].has_key('sequence'):
                entityRecords[i]['onlyAlphaCarbon'] = onlyAlphaCarbon[i]
    
        return entryID, entityRecords, experimentalType, resolution, rFactorObs, rFactorFree

    def extract_info(self, mmCIFFile, tokens='all'):
        """The parsing basically returns a dictionary with the same structure as the mmCIF file. all the top level block names like _entity and _entry
        have an entry in the returned token dictionary. The keys for each of these top level entries are the sub block names like id and type. Each data
        entry (so for example the data recorded at tokenDict['_entry']['id'] is a list. If there was a loop in the block then the list will have more
        than one entry, else it will have one. If there is more than one entry then each data item within the sub blocks is ordered in the same manner.
        For example, if tokenDict[i][j] == ['a', 'b', 'c'] and tokenDict[i][k] == [1, 2, 3], then one data record in the block had i.j == 'a' and i.k == 1.
        Similarly, i.j == 'b' and i.k == 2, and i.j == 'c' and i.k == 3.
        """
        assert(type(tokens) == list)
    
        tokenDict = self.parse_mmCIF(mmCIFFile)
    
        errorMessage = ''
        if tokens == 'all':
            return tokenDict, errorMessage
    
        tokensFound = tokenDict.keys()
        subDict = {}
        for i in tokens:
            if '.' in i:
                # If this is True, then the token is something like _entity.id not just _entity.
                chunks = i.split('.')
                mainToken = chunks[0]
                subToken = chunks[1]
                if mainToken in tokensFound:
                    if subToken in tokenDict[mainToken].keys():
                        if subDict.has_key(mainToken):
                            subDict[mainToken][subToken] = tokenDict[mainToken][subToken]
                        else:
                            subDict[mainToken] = {}
                            subDict[mainToken][subToken] = tokenDict[mainToken][subToken]
                    else:
                        errorMessage += 'Main token ' + mainToken + ' found, but sub-token ' + subToken + ' from requested token ' + i + ' not found.\n'
                else:
                    errorMessage += 'Main token ' + mainToken + ' from requested token ' + i + ' not found.\n'
            else:
                # If this is True, then the token is something like _entity not _entity.id or _entity.type.
                if i in tokensFound:
                    subDict[i] = tokenDict[i]
                else:
                    errorMessage += 'Main token ' + i + ' not found.\n'
        return subDict, errorMessage
    
    def parse_mmCIF(self, mmCIFFile):
        """Parse an mmCIF file.
        
        Only parses a single file. Parses files that record only one PDB entry in them.
        
        @param mmCIFFile: The location of the gzipped mmCIF file to parse.
        @type mmCIFFile: string
        
        """
    
        readFile = gzip.open(mmCIFFile, 'r')
        mmCIFContent = readFile.read()
        readFile.close()

        fileChunks = re.split('(?<=\n)# \n', mmCIFContent)[:-1]
    
        entryID = fileChunks[1].split()[1]
        dataDictionary = {}
        entityIDs = []

        for i in fileChunks[1:]:
            # As the lines all end with newline characters, the final list element is always ''. [:-1] removes this.
            # It is necessary to do it this way, rather than splitting (i.split()), as this method ensures that all lines that should end with a ' ' do so.
            blockChunks = i.split('\n')[:-1]
            if blockChunks[0] == 'loop_':
                # A loop has been found, and therefore the block contains more than one data record.
                blockDictionary = {}
                # Record the name of the block.
                blockName = re.match('_[a-zA-Z0-9_]+', blockChunks[1])
                blockName = blockName.group(0)
                subBlockNameOrder = []  # Stored the subBlockNames in the order that they are encountered.
                currentNameOrderIndex = 0
                subBlockData = ''
                lookingAtLongDataReocrd = False
                longDataRecordDelimiter = ''
                secondLongDataReocrdDelimiter = False
                lookingAtMultilineDataReocrd = False
                dictSetUp = False
                for j in blockChunks[1:]:
                    subBlockNameSearch = re.search('(?<=' + blockName + '.)[a-zA-Z0-9_]+', j)
                    if subBlockNameSearch:
                        # If this is True, then the line contains something like _entity.id or _cell.angle_gamma_esd.
                        subBlockNameOrder.append(subBlockNameSearch.group(0))
                    else:
                        # If this is True, then the line contains a data record.
                        if not dictSetUp:
                            blockDictionary = dict([(k, []) for k in subBlockNameOrder])
                            dictSetUp = True
                        previousCharacter = ''
                        for k in j:
                            # Go through every character on the line.
                            if lookingAtMultilineDataReocrd:
                                if previousCharacter == '' and k == ';':
                                    # If this is True, then the beginning of a new line has been found, and the first character on the new line is a ;.
                                    # This means that the end of a multiline data record has been found.
                                    blockDictionary[subBlockNameOrder[currentNameOrderIndex]].append(subBlockData)
                                    lookingAtMultilineDataReocrd = False
                                    currentNameOrderIndex = (currentNameOrderIndex + 1) % len(subBlockNameOrder)
                                    subBlockData = ''
                                else:
                                    # If this is True, then the current character is in the middle of a multiline data record.
                                    subBlockData += k
                                previousCharacter = k
                            elif lookingAtLongDataReocrd:
                                if k == longDataRecordDelimiter:
                                    # If this is True, then the end of a '' or "" delimited data record has been found.
                                    subBlockData += k
                                    secondLongDataReocrdDelimiter = True
#                                    lookingAtLongDataReocrd = False
#                                    longDataRecordDelimiter = ''
                                elif previousCharacter == longDataRecordDelimiter and k == ' ' and secondLongDataReocrdDelimiter:
                                    lookingAtLongDataReocrd = False
                                    longDataRecordDelimiter = ''
                                    blockDictionary[subBlockNameOrder[currentNameOrderIndex]].append(subBlockData[:-1])
                                    currentNameOrderIndex = (currentNameOrderIndex + 1) % len(subBlockNameOrder)
                                    subBlockData = ''
                                else:
                                    # If this is True, then the current character is in the middle of a long data record.
                                    subBlockData += k
                                previousCharacter = k
                            else:
                                if previousCharacter == '' and k == ';':
                                    # If this is True, then the beginning of a multiline data record has been found.
                                    lookingAtMultilineDataReocrd = True
                                    previousCharacter = k
                                elif k == '\'' or k == '"':
                                    # If this is True, then the beginning of a long data record has been found.
                                    lookingAtLongDataReocrd = True
                                    secondLongDataReocrdDelimiter = False
                                    longDataRecordDelimiter = k
                                    previousCharacter = k
                                elif k == ' ':
                                    if subBlockData != '':
                                        # If this is True, then the end of the current data record has been found.
                                        blockDictionary[subBlockNameOrder[currentNameOrderIndex]].append(subBlockData)
                                        currentNameOrderIndex = (currentNameOrderIndex + 1) % len(subBlockNameOrder)
                                        subBlockData = ''
                                    previousCharacter = k
                                else:
                                    # If this is True, then the current character is in the middle of a data record.
                                    subBlockData += k
                                    previousCharacter = k
            else:
                # There is no loop, and therefore the block only contains one data record.
                blockDictionary = {}
                # Record the name of the block.
                blockName = re.match('_[a-zA-Z0-9_]+', blockChunks[0])
                blockName = blockName.group(0)
                lookingAtMultilineDataReocrd = False
                for j in blockChunks:
                    subBlockNameSearch = re.search('(?<=' + blockName + '.)[a-zA-Z0-9_]+', j)
                    if subBlockNameSearch:
                        # If this is True, then the line contains something like _entity.id or _cell.angle_gamma_esd.
                        subBlockName = subBlockNameSearch.group(0)
                        subBlockChunks = j.split(None, 1)
                        if len(subBlockChunks) > 1:
                            # If this is True, then there is data on the line along with the subBlockName.
                            subBlockData = subBlockChunks[1].strip()
                            if subBlockData[-1] == '\'' or subBlockData[-1] == '"':
                                # Strip off leading and trailing '' or "".
                                blockDictionary[subBlockName] = [subBlockData[1:-1]]
                            else:
                                blockDictionary[subBlockName] = [subBlockData]
                        else:
                            # There is no data on the line with the subBlockName. This means that the data is spread over
                            # multiple lines, and is flanked by semi-colons.
                            pass
                    else:
                        # There is no subBlockName on the line. This means that the line contains some data to go along with
                        # the most recently found subBlockName.
                        if lookingAtMultilineDataReocrd:
                            if j.strip() == ';':
                                # If this is True, then the line ends a data record.
                                lookingAtMultilineDataReocrd = False
                                blockDictionary[subBlockName] = [subBlockData]
                            else:
                                subBlockData += j
                        elif (j[0] == '\'' and j[-2:] == '\' ') or (j[0] == '"' and j[-2:] == '" '):
                            # If this is true then the line is taken up with a long piece of data that is surrounded by '' or "".
                            blockDictionary[subBlockName] = [j.strip()[1:-1]]
                        elif j[0] == ';':
                            # If this is True, then the line starts a data record.
                            lookingAtMultilineDataReocrd = True
                            subBlockData = j[1:]
                        else:
                            # If this is True, then the line just contains a data record that is too long for one line, but does not need/use
                            # the '' or "" long data record delimiters.
                            blockDictionary[subBlockName] = [j.strip()]           

            dataDictionary[blockName] = blockDictionary
    
        return dataDictionary