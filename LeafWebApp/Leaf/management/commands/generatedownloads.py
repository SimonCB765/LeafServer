import os
import sys
import datetime
import time
import subprocess
from gzip import GzipFile
from cStringIO import StringIO
import tarfile

from django.core.management.base import BaseCommand
from django.core.files import File
from django.core.files.base import ContentFile
from django.core.files.uploadedfile import InMemoryUploadedFile
from django.core.files.temp import NamedTemporaryFile
from django.conf import settings

from Leaf.models import PDBCullRequest
from Leaf.models import ProteinInformation
from Leaf.models import Representative
from Leaf.models import EntryRepresentative
from Leaf.models import DownloadableFiles
from Leaf.models import ChainType
from Leaf.models import AllPDBEntries
from Leaf.models import Similarity

srcLocation = os.path.abspath(__file__)
srcLocation = '/'.join(srcLocation.split('/')[:-3])
sys.path.append(srcLocation)
import cullinput.performBLAST
import cullinput.controlthread

class Command(BaseCommand):

    def handle(self, *args, **options):
        os.umask(00000)
        logger = open('/srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp/ErrorLogs/GENERATION.log', 'a')
        logger.write('=======================================================\n=======================================================\nStarting generation on ' + time.strftime('%Y/%m/%d/ at %H:%M:%S', time.gmtime()) + '.\n')
        logger.close()

        # Turn all the protein chain information into a dictionary.
        pdbChainsData = ProteinInformation.objects.all()
        pdbChainsDict = dict([(i.chain, i) for i in pdbChainsData])
        del(pdbChainsData)
        
        # Turn all the within entry representative information into a dictionary.
        entryRepData = EntryRepresentative.objects.all()
        entryRepIgnore = []
        entryRepDict = {}
        for i in entryRepData:
            entryRepIgnore.append(i.nonreprChain)
            reprChain = i.reprChain
            if entryRepDict.has_key(reprChain):
                entryRepDict[reprChain].append(i.nonreprChain)
            else:
                entryRepDict[reprChain] = [i.nonreprChain]
        del(entryRepData)
        
        # Turn all the whole PDB representative information into a dictionary.
        wholeRepData = Representative.objects.all()
        wholeRepIgnore = []
        wholeRepDict = {}
        for i in wholeRepData:
            wholeRepIgnore.append(i.nonreprChain)
            reprChain = i.reprChain
            if wholeRepDict.has_key(reprChain):
                wholeRepDict[reprChain].append(i.nonreprChain)
            else:
                wholeRepDict[reprChain] = [i.nonreprChain]
        del(wholeRepData)

        # Generate the whole PDB fasta files.
        pdbaa = StringIO()
        gzf = GzipFile(fileobj=pdbaa, mode='wb')
        for i in sorted(pdbChainsDict.keys()):
            experimentType = str(pdbChainsDict[i].experimentType)
            resolution = str(pdbChainsDict[i].resolution)
            rValueObs = str(pdbChainsDict[i].rValueObs)
            rValueFree = str(pdbChainsDict[i].rValueFree)
            alphaCarbonOnly = ['no', 'yes'][pdbChainsDict[i].alphaCarbonOnly]
#            alphaCarbonOnly = 'yes' if pdbChainsDict[i].alphaCarbonOnly else 'no'
            description = str(pdbChainsDict[i].description)
            dbName = str(pdbChainsDict[i].dbName)
            dbCode = str(pdbChainsDict[i].dbCode)
            organism = str(pdbChainsDict[i].organism)
            sequence = str(pdbChainsDict[i].sequence)
            gzf.write('>' + str(i) + '\t' + str(len(sequence)) + '\t' + experimentType + '\t' + resolution + '\t' +
                      rValueObs + '\t' + rValueFree + '\t' + alphaCarbonOnly + '\t' + description + '\t<' + dbName +
                      ' ' + dbCode + '>\t[' + organism + ']\n')
            gzf.write(sequence + '\n')
        gzf.close()
        pdbaa.seek(0, 2)
        pdbaaFile = InMemoryUploadedFile(pdbaa, 'image', 'pdbaa.gz', None, pdbaa.tell(), None)
        pdbaaEntry = DownloadableFiles(fileName='pdbaa')
        pdbaaEntry.save()
        pdbaaEntry.downloadFile.save(pdbaaFile.name, pdbaaFile)
        pdbaa.close()
        
        pdbaaent = StringIO()
        gzf = GzipFile(fileobj=pdbaaent, mode='wb')
        for i in sorted([j for j in pdbChainsDict.keys() if not j in entryRepIgnore]):
            experimentType = str(pdbChainsDict[i].experimentType)
            resolution = str(pdbChainsDict[i].resolution)
            rValueObs = str(pdbChainsDict[i].rValueObs)
            rValueFree = str(pdbChainsDict[i].rValueFree)
            alphaCarbonOnly = ['no', 'yes'][pdbChainsDict[i].alphaCarbonOnly]
#            alphaCarbonOnly = 'yes' if pdbChainsDict[i].alphaCarbonOnly else 'no'
            description = str(pdbChainsDict[i].description)
            dbName = str(pdbChainsDict[i].dbName)
            dbCode = str(pdbChainsDict[i].dbCode)
            organism = str(pdbChainsDict[i].organism)
            sequence = str(pdbChainsDict[i].sequence)
            if entryRepDict.has_key(i):
                gzf.write('>' + str(i) + '\t' + str(len(sequence)) + '\t' + experimentType + '\t' + resolution + '\t' +
                          rValueObs + '\t' + rValueFree + '\t' + alphaCarbonOnly + '\t' + description + '\t<' + dbName +
                          ' ' + dbCode + '>\t[' + organism + ']\t||\t' + '\t'.join(entryRepDict[i]) + '\n')
            else:
                gzf.write('>' + str(i) + '\t' + str(len(sequence)) + '\t' + experimentType + '\t' + resolution + '\t' +
                          rValueObs + '\t' + rValueFree + '\t' + alphaCarbonOnly + '\t' + description + '\t<' + dbName +
                          ' ' + dbCode + '>\t[' + organism + ']\n')
            gzf.write(sequence + '\n')
        gzf.close()
        pdbaaent.seek(0, 2)
        pdbaaentFile = InMemoryUploadedFile(pdbaaent, 'image', 'pdbaaent.gz', None, pdbaaent.tell(), None)
        pdbaaentEntry = DownloadableFiles(fileName='pdbaaent')
        pdbaaentEntry.save()
        pdbaaentEntry.downloadFile.save(pdbaaentFile.name, pdbaaentFile)
        pdbaaent.close()
        
        pdbaanr = StringIO()
        gzf = GzipFile(fileobj=pdbaanr, mode='wb')
        for i in sorted([j for j in pdbChainsDict.keys() if not j in wholeRepIgnore]):
            sequence = str(pdbChainsDict[i].sequence)
            experimentType = str(pdbChainsDict[i].experimentType)
            resolution = str(pdbChainsDict[i].resolution)
            rValueObs = str(pdbChainsDict[i].rValueObs)
            rValueFree = str(pdbChainsDict[i].rValueFree)
            alphaCarbonOnly = ['no', 'yes'][pdbChainsDict[i].alphaCarbonOnly]
#            alphaCarbonOnly = 'yes' if pdbChainsDict[i].alphaCarbonOnly else 'no'
            description = str(pdbChainsDict[i].description)
            dbName = str(pdbChainsDict[i].dbName)
            dbCode = str(pdbChainsDict[i].dbCode)
            organism = str(pdbChainsDict[i].organism)
            if wholeRepDict.has_key(i):
                gzf.write('>' + str(i) + '\t' + str(len(sequence)) + '\t' + experimentType + '\t' + resolution + '\t' +
                          rValueObs + '\t' + rValueFree + '\t' + alphaCarbonOnly + '\t' + description + '\t<' + dbName +
                          ' ' + dbCode + '>\t[' + organism + ']\t||\t' + '\t'.join(wholeRepDict[i]) + '\n')
            else:
                gzf.write('>' + str(i) + '\t' + str(len(sequence)) + '\t' + experimentType + '\t' + resolution + '\t' +
                          rValueObs + '\t' + rValueFree + '\t' + alphaCarbonOnly + '\t' + description + '\t<' + dbName +
                          ' ' + dbCode + '>\t[' + organism + ']\n')
            gzf.write(sequence + '\n')
        gzf.close()
        pdbaanr.seek(0, 2)
        pdbaanrFile = InMemoryUploadedFile(pdbaanr, 'image', 'pdbaanr.gz', None, pdbaanr.tell(), None)
        pdbaanrEntry = DownloadableFiles(fileName='pdbaanr')
        pdbaanrEntry.save()
        pdbaanrEntry.downloadFile.save(pdbaanrFile.name, pdbaanrFile)
        pdbaanr.close()

        # Generate the files needed to run the culling locally (i.e. make files of the different database tables).
        proteinData = ProteinInformation.objects.all().values_list()
        proteinFile = NamedTemporaryFile()
#        sio = StringIO()
#        gzf = GzipFile(fileobj=sio, mode='wb')
        for i in proteinData:
#            gzf.write('\t'.join([str(j) for j in i]) + '\n')
            proteinFile.write('\t'.join([str(j) for j in i]) + '\n')
#        gzf.close()
#        sio.seek(0, 2)
#        proteinFile = InMemoryUploadedFile(sio, 'image', 'ProteinInformation.gz', None, sio.tell(), None)
        proteinFileEntry = DownloadableFiles(fileName='ProteinInformation', downloadFile=File(proteinFile))
#        proteinFileEntry = DownloadableFiles(fileName='ProteinInformation')
        proteinFileEntry.save()
#        proteinFileEntry.downloadFile.save(proteinFile.name, proteinFile)
#        sio.close()
        proteinFile.close()
        
        chainData = ChainType.objects.all().values_list()
        chainFile = NamedTemporaryFile()
#        sio = StringIO()
#        gzf = GzipFile(fileobj=sio, mode='wb')
        for i in chainData:
#            gzf.write('\t'.join([str(j) for j in i]) + '\n')
            chainFile.write('\t'.join([str(j) for j in i]) + '\n')
#        gzf.close()
#        sio.seek(0, 2)
#        chainFile = InMemoryUploadedFile(sio, 'image', 'ChainType.gz', None, sio.tell(), None)
        chainFileEntry = DownloadableFiles(fileName='ChainType', downloadFile=File(chainFile))
#        chainFileEntry = DownloadableFiles(fileName='ChainType')
        chainFileEntry.save()
#        chainFileEntry.downloadFile.save(chainFile.name, chainFile)
#        sio.close()
        chainFile.close()
        
        PDBEntriesData = AllPDBEntries.objects.all().values_list()
        entriesFile = NamedTemporaryFile()
#        sio = StringIO()
#        gzf = GzipFile(fileobj=sio, mode='wb')
        for i in PDBEntriesData:
#            gzf.write('\t'.join([str(j) for j in i]) + '\n')
            entriesFile.write('\t'.join([str(j) for j in i]) + '\n')
#        gzf.close()
#        sio.seek(0, 2)
#        entriesFile = InMemoryUploadedFile(sio, 'image', 'AllPDBEntries.gz', None, sio.tell(), None)
        entriesFileEntry = DownloadableFiles(fileName='AllPDBEntries', downloadFile=File(entriesFile))
#        entriesFileEntry = DownloadableFiles(fileName='AllPDBEntries')
        entriesFileEntry.save()
#        entriesFileEntry.downloadFile.save(entriesFile.name, entriesFile)
#        sio.close()
        entriesFile.close()

        similarityData = Similarity.objects.all().values_list()
        similarityFile = NamedTemporaryFile()
#        sio = StringIO()
#        gzf = GzipFile(fileobj=sio, mode='wb')
        for i in similarityData:
#            gzf.write('\t'.join([str(j) for j in i]) + '\n')
            similarityFile.write('\t'.join([str(j) for j in i]) + '\n')
#        gzf.close()
#        sio.seek(0, 2)
#        similarityFile = InMemoryUploadedFile(sio, 'image', 'Similarity.gz', None, sio.tell(), None)
        similarityFileEntry = DownloadableFiles(fileName='Similarity', downloadFile=File(similarityFile))
#        similarityFileEntry = DownloadableFiles(fileName='Similarity')
        similarityFileEntry.save()
#        similarityFileEntry.downloadFile.save(similarityFile.name, similarityFile)
#        sio.close()
        similarityFile.close()

        representativeData = Representative.objects.all().values_list()
        representativeFile = NamedTemporaryFile()
#        sio = StringIO()
#        gzf = GzipFile(fileobj=sio, mode='wb')
        for i in representativeData:
#            gzf.write('\t'.join([str(j) for j in i]) + '\n')
            representativeFile.write('\t'.join([str(j) for j in i]) + '\n')
#        gzf.close()
#        sio.seek(0, 2)
#        representativeFile = InMemoryUploadedFile(sio, 'image', 'Representative.gz', None, sio.tell(), None)
        representativeFileEntry = DownloadableFiles(fileName='Representative', downloadFile=File(representativeFile))
#        representativeFileEntry = DownloadableFiles(fileName='Representative')
        representativeFileEntry.save()
#        representativeFileEntry.downloadFile.save(representativeFile.name, representativeFile)
#        sio.close()
        representativeFile.close()
        
        destinationDir = settings.MEDIA_ROOT + 'TarData/PDBData'
#        tar = tarfile.open(settings.MEDIA_ROOT + 'TarData/PDBData.tar.gz', mode='w:gz')
        entriesInfo = DownloadableFiles.objects.filter(fileName__exact='AllPDBEntries')
        subprocess.call(['cp', str(entriesInfo[0].downloadFile.path), destinationDir + '/AllPDBEntries.txt'])
#        tar.add(entriesInfo[0].downloadFile.path, arcname='AllPDBEntries.txt')
        chainInfo = DownloadableFiles.objects.filter(fileName__exact='ChainType')
        subprocess.call(['cp', str(chainInfo[0].downloadFile.path), destinationDir + '/ChainType.txt'])
#        tar.add(chainInfo[0].downloadFile.path, arcname='ChainType.txt')
        proteinInfo = DownloadableFiles.objects.filter(fileName__exact='ProteinInformation')
        subprocess.call(['cp', str(proteinInfo[0].downloadFile.path), destinationDir + '/ProteinInformation.txt'])
#        tar.add(proteinInfo[0].downloadFile.path, arcname='ProteinInformation.txt')
        representativeInfo = DownloadableFiles.objects.filter(fileName__exact='Representative')
        subprocess.call(['cp', str(representativeInfo[0].downloadFile.path), destinationDir + '/Representative.txt'])
#        tar.add(representativeInfo[0].downloadFile.path, arcname='Representative.txt')
        similarityInfo = DownloadableFiles.objects.filter(fileName__exact='Similarity')
        subprocess.call(['cp', str(similarityInfo[0].downloadFile.path), destinationDir + '/Similarity.txt'])
#        tar.add(similarityInfo[0].downloadFile.path, arcname='Similarity.txt')
#        tar.close()
        subprocess.Popen(['tar', '-zcvf', 'PDBData.tar.gz', 'PDBData'], cwd=settings.MEDIA_ROOT + 'TarData', stdout=subprocess.PIPE, stderr=subprocess.STDOUT)


        proteinChains = []
        readChainTypeData = open(destinationDir + '/ChainType.txt', 'r')
        for i in readChainTypeData:
            chunks = (i.strip()).split('\t')
            if chunks[1] == 'Protein':
                proteinChains.append(chunks[0])
        readChainTypeData.close()
        proteinChains = '\n'.join(proteinChains)

        # Generate the culled PDB lists.
        resolutionsOfInterest = [1.6, 1.8, 2.0, 2.2, 2.5, 3.0]
        rValuesOfInterest = [0.25, 1.0]
        sequenceIdentities = [20, 25, 30, 40, 50, 60, 70, 80, 90]
        tuplesToDo = [(i, j, k) for i in resolutionsOfInterest for j in rValuesOfInterest for k in sequenceIdentities]
        tuplesToDo.extend([(100.0, 1.0, i) for i in sequenceIdentities])
        updatesPerformed = []
        while tuplesToDo != [] or updatesPerformed != []:
            # Continue looping until both lists are empty. If both are empty then all updates have been performed.
            updatesDone = []
            for i in updatesPerformed:
                # Determine if any of the culling operations that are currently running have finished.
                if i.completed:
                    if not i.skipNonXray:
                        if not i.skipAlphaCarbon:
                            xrayCAInfo = '_INCLNONXRAY_INCLCAONLY'
                        else:
                            xrayCAInfo = '_INCLNONXRAY'
                    else:
                        xrayCAInfo = ''

                    sio = StringIO()
                    gzf = GzipFile(fileobj=sio, mode='wb')
                    outputInfo = open(i.nonredSeq.path, 'r')
                    for line in outputInfo:
                        gzf.write(line)
                    outputInfo.close()
                    gzf.close()
                    sio.seek(0, 2)
                    culledFile = InMemoryUploadedFile(sio, 'image',
                                                      ('SeqIden_' + str(i.sequenceIdentity) + '_Res_' + str(i.maxResolution) +
                                                       '_RVal_' + str(i.maxRValue) + xrayCAInfo + '.gz'),
                                                       None, sio.tell(), None)
                    newDownload = DownloadableFiles(fileName='SeqIden_' + str(i.sequenceIdentity) + '_Res_' + str(i.maxResolution) +
                                                    '_RVal_' + str(i.maxRValue) + xrayCAInfo)
                    newDownload.save()
                    newDownload.downloadFile.save(culledFile.name, culledFile)
                    sio.close()
                    updatesDone.append(i)
                    i.delete()
            for i in updatesDone:
                updatesPerformed.remove(i)

            if tuplesToDo == [] or len(updatesPerformed) >= 1:
                # Only allow 1 culling operation to occur at once, to prevent overloading.
#                pass
                time.sleep(30)
            else:
                # Fill the buffer of culling operations up.
                while len(updatesPerformed) < 1 and tuplesToDo != []:
                    newCull = tuplesToDo.pop()
                    resolution = newCull[0]
                    rValue = newCull[1]
                    seqIdentity = newCull[2]
                    if resolution == 100.0:
                        skipNonXray = False
                        skipAlphaCarbon = False
                    else:
                        skipNonXray = True
                        skipAlphaCarbon = True
                    r = PDBCullRequest(
                                       wholePDB=False,
                                       sequenceIdentity=seqIdentity,
                                       minResolution=0.0,
                                       maxResolution=resolution,
                                       maxRValue=rValue,
                                       minLength = -1,
                                       maxLength = -1,
                                       skipNonXray=skipNonXray,
                                       skipAlphaCarbon=skipAlphaCarbon,
                                       cullByChain=True,
                                       performIntraEntryCulling=False,
                                       intraEntrySequenceIdentity=100.0,
                                       email='No Email',
                                       completed=False,
                                       requestDate=datetime.datetime.now()
                                       )
                    r.save()
                    r.userInput.save('', ContentFile(proteinChains))
                    cullinput.controlthread.RunCull(r, 'pdb')
                    updatesPerformed.append(r)

        logger = open('/srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp/ErrorLogs/GENERATION.log', 'a')
        logger.write('Finished generation on ' + time.strftime('%Y/%m/%d/ at %H:%M:%S', time.gmtime()) + '.\n')
        logger.close()