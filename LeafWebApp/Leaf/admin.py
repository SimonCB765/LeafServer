from django.contrib import admin

from Leaf.models import UserCullRequest
from Leaf.models import PDBCullRequest
from Leaf.models import ProteinInformation
from Leaf.models import ChainType
from Leaf.models import AllPDBEntries
from Leaf.models import Similarity
from Leaf.models import Representative
from Leaf.models import EntryRepresentative
from Leaf.models import DownloadableFiles

class UserCullAdmin(admin.ModelAdmin):
    fieldsets = [
        ('Request Information', {'fields': ['userInput', 'sequenceIdentity', 'SEG', 'completed']}),
        ('Result Information',  {'fields': ['nonredNoSeq', 'nonredSeq']}),
        ('Contact Information', {'fields': ['email']}),
        ('Timestamp',           {'fields': ['requestDate']}),
    ]

class PDBCullAdmin(admin.ModelAdmin):
    fieldsets = [
        ('Request Information', {'fields': ['userInput', 'sequenceIdentity', 'minResolution', 'maxResolution', 'maxRValue',
                                            'skipNonXray', 'skipAlphaCarbon', 'cullByChain',
                                            'performIntraEntryCulling', 'intraEntrySequenceIdentity', 'completed']}),
        ('Contact Information', {'fields': ['email']}),
        ('Timestamp',           {'fields': ['requestDate']}),
    ]

class ProteinInformationAdmin(admin.ModelAdmin):
    fieldsets = [
        ('Protein Information', {'fields': ['chain', 'entry', 'experimentType', 'resolution', 'rValueObs', 'rValueFree',
                                            'alphaCarbonOnly', 'description', 'dbName', 'dbCode', 'organism', 'sequence']})
    ]
    
class ChainTypeAdmin(admin.ModelAdmin):
    fieldsets = [
        ('Chain Type Information', {'fields': ['chain', 'chainType']})
    ]

class AllPDBEntriesAdmin(admin.ModelAdmin):
    fieldsets = [
        ('Entries', {'fields': ['entry']})
    ]

class SimilarityAdmin(admin.ModelAdmin):
    fieldsets = [
        ('Similarity Information', {'fields': ['chainA', 'entryA', 'chainB', 'entryB', 'similarity', 'matchLength']})
    ]

class RepresentativeAdmin(admin.ModelAdmin):
    fieldsets = [
        ('Representative Information', {'fields': ['nonreprChain', 'reprChain']})
    ]

class EntryRepresentativeAdmin(admin.ModelAdmin):
    fieldsets = [
        ('Representative Information', {'fields': ['nonreprChain', 'reprChain']})
    ]

class DownloadableFilesAdmin(admin.ModelAdmin):
    fieldsets = [
        ('File Information', {'fields': ['fileName', 'downloadFile']})
    ]


admin.site.register(UserCullRequest, UserCullAdmin)
admin.site.register(PDBCullRequest, PDBCullAdmin)
admin.site.register(ProteinInformation, ProteinInformationAdmin)
admin.site.register(ChainType, ChainTypeAdmin)
admin.site.register(AllPDBEntries, AllPDBEntriesAdmin)
admin.site.register(Similarity, SimilarityAdmin)
admin.site.register(Representative, RepresentativeAdmin)
admin.site.register(EntryRepresentative, EntryRepresentativeAdmin)
admin.site.register(DownloadableFiles, DownloadableFilesAdmin)