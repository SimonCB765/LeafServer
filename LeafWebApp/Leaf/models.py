from time import gmtime, strftime
import os

from django.db import models
from django.core.files.storage import default_storage

import settings

class UserCullRequest(models.Model):
    def upload_input_path(self, filename):
        uploadLoc = strftime('%Y/%m/%d/', gmtime()) + str(self.id) + '/'
        if not os.path.exists(settings.MEDIA_ROOT + uploadLoc):
            os.umask(00000)
            os.makedirs(settings.MEDIA_ROOT + uploadLoc, 0777)
        return uploadLoc + 'Input.txt'
    
    def upload_removed_path(self, filename):
        uploadLoc = strftime('%Y/%m/%d/', gmtime()) + str(self.id) + '/'
        if not os.path.exists(settings.MEDIA_ROOT + uploadLoc):
            os.umask(00000)
            os.makedirs(settings.MEDIA_ROOT + uploadLoc, 0777)
        return uploadLoc + 'Removed.txt'

    def upload_nonredNoSeq_path(self, filename):
        uploadLoc = strftime('%Y/%m/%d/', gmtime()) + str(self.id) + '/'
        if not os.path.exists(settings.MEDIA_ROOT + uploadLoc):
            os.umask(00000)
            os.makedirs(settings.MEDIA_ROOT + uploadLoc, 0777)
        return uploadLoc + 'NonRedundantList.txt'

    def upload_nonredSeq_path(self, filename):
        uploadLoc = strftime('%Y/%m/%d/', gmtime()) + str(self.id) + '/'
        if not os.path.exists(settings.MEDIA_ROOT + uploadLoc):
            os.umask(00000)
            os.makedirs(settings.MEDIA_ROOT + uploadLoc, 0777)
        return uploadLoc + 'NonRedundantFasta.fasta'
    
    id = models.AutoField(primary_key=True)
    userInput = models.FileField(upload_to=upload_input_path)
    sequenceIdentity = models.FloatField()
    minLength = models.IntegerField()
    maxLength = models.IntegerField()
    SEG = models.BooleanField()
    email = models.CharField(max_length=200)
    completed = models.BooleanField()
    nonredNoSeq = models.FileField(upload_to=upload_nonredNoSeq_path)
    nonredSeq = models.FileField(upload_to=upload_nonredSeq_path)
    removed = models.FileField(upload_to=upload_removed_path)
    requestDate = models.DateTimeField()
    
    def __unicode__(self):
        return str(self.email) + ' at ' + str(self.requestDate)

class PDBCullRequest(models.Model):
    def upload_input_path(self, filename):
        uploadLoc = strftime('%Y/%m/%d/', gmtime()) + str(self.id) + '/'
        if not os.path.exists(settings.MEDIA_ROOT + uploadLoc):
            os.umask(00000)
            os.makedirs(settings.MEDIA_ROOT + uploadLoc, 0777)
        return uploadLoc + 'Input.txt'
    
    def upload_removed_path(self, filename):
        uploadLoc = strftime('%Y/%m/%d/', gmtime()) + str(self.id) + '/'
        if not os.path.exists(settings.MEDIA_ROOT + uploadLoc):
            os.umask(00000)
            os.makedirs(settings.MEDIA_ROOT + uploadLoc, 0777)
        return uploadLoc + 'Removed.txt'

    def upload_nonredNoSeq_path(self, filename):
        uploadLoc = strftime('%Y/%m/%d/', gmtime()) + str(self.id) + '/'
        if not os.path.exists(settings.MEDIA_ROOT + uploadLoc):
            os.umask(00000)
            os.makedirs(settings.MEDIA_ROOT + uploadLoc, 0777)
        return uploadLoc + 'NonRedundantList.txt'

    def upload_nonredSeq_path(self, filename):
        uploadLoc = strftime('%Y/%m/%d/', gmtime()) + str(self.id) + '/'
        if not os.path.exists(settings.MEDIA_ROOT + uploadLoc):
            os.umask(00000)
            os.makedirs(settings.MEDIA_ROOT + uploadLoc, 0777)
        return uploadLoc + 'NonRedundantFasta.fasta'

    def upload_pairwiseScores_path(self, filename):
        uploadLoc = strftime('%Y/%m/%d/', gmtime()) + str(self.id) + '/'
        if not os.path.exists(settings.MEDIA_ROOT + uploadLoc):
            os.umask(00000)
            os.makedirs(settings.MEDIA_ROOT + uploadLoc, 0777)
        return uploadLoc + 'PairwiseScores.txt'
    
    id = models.AutoField(primary_key=True)
    wholePDB = models.BooleanField()
    userInput = models.FileField(upload_to=upload_input_path)
    sequenceIdentity = models.FloatField()
    minResolution = models.FloatField()
    maxResolution = models.FloatField()
    maxRValue = models.FloatField()
    minLength = models.IntegerField()
    maxLength = models.IntegerField()
    skipNonXray = models.BooleanField()
    skipAlphaCarbon = models.BooleanField()
    cullByChain = models.BooleanField()
    performIntraEntryCulling = models.BooleanField()
    intraEntrySequenceIdentity = models.FloatField()
    email = models.CharField(max_length=200)
    completed = models.BooleanField()
    nonredNoSeq = models.FileField(upload_to=upload_nonredNoSeq_path)
    nonredSeq = models.FileField(upload_to=upload_nonredSeq_path)
    pairwiseScores = models.FileField(upload_to=upload_pairwiseScores_path)
    removed = models.FileField(upload_to=upload_removed_path)
    requestDate = models.DateTimeField()
    
    def __unicode__(self):
        return str(self.email) + ' at ' + str(self.requestDate)

class ProteinInformation(models.Model):
    chain = models.CharField(max_length=5, primary_key=True)
    entry = models.CharField(max_length=4)
    experimentType = models.CharField(max_length=20)
    resolution = models.FloatField()
    rValueObs = models.FloatField()
    rValueFree = models.FloatField()
    alphaCarbonOnly = models.BooleanField()
    description = models.TextField()
    dbName = models.CharField(max_length=20)
    dbCode = models.CharField(max_length=20)
    organism = models.CharField(max_length=200)
    sequence = models.TextField()
    
    def __unicode__(self):
        return str(self.chain)

class ChainType(models.Model):
    chain = models.CharField(max_length=5, primary_key=True)
    chainType = models.CharField(max_length=40)
    
    def __unicode__(self):
        return str(self.chain) + ' ' + str(self.chainType)

class AllPDBEntries(models.Model):
    entry = models.CharField(max_length=4, primary_key=True)
    
    def __unicode__(self):
        return str(self.entry)

class Similarity(models.Model):
    id = models.AutoField(primary_key=True)
    chainA = models.CharField(max_length=5)
    entryA = models.CharField(max_length=4)
    chainB = models.CharField(max_length=5)
    entryB = models.CharField(max_length=4)
    similarity = models.FloatField()
    matchLength = models.IntegerField()
    
    def __unicode__(self):
        return str(self.chainA) + ' ' + str(self.chainB) + ' ' + str(self.similarity)

class Representative(models.Model):
    nonreprChain = models.CharField(max_length=5, primary_key=True)
    reprChain = models.CharField(max_length=5)

    def __unicode__(self):
        return str(self.nonreprChain) + ' ' + str(self.reprChain)

class EntryRepresentative(models.Model):
    nonreprChain = models.CharField(max_length=5, primary_key=True)
    reprChain = models.CharField(max_length=5)
    
    def __unicode__(self):
        return str(self.nonreprChain) + ' ' + str(self.reprChain)

class DownloadableFiles(models.Model):
    fileName = models.CharField(max_length=150, primary_key=True)
    downloadFile = models.FileField(upload_to='downloadable')
    
    def __unicode__(self):
        return str(self.fileName)
    
    def save(self, force_insert=False, force_update=False):
        try:
            old_obj = DownloadableFiles.objects.get(pk=self.pk)
            path = old_obj.downloadFile.path
            default_storage.delete(path)
        except:
            pass
        super(DownloadableFiles, self).save(force_insert, force_update)