from django.http import HttpResponseRedirect, HttpResponse, Http404
from django.core.files import File
from django.core.files.base import ContentFile
from django.core.servers.basehttp import FileWrapper
from django.core.urlresolvers import reverse
from django.shortcuts import render_to_response, get_object_or_404
from django.template import RequestContext
from django.conf import settings

import re
import datetime
import os
import shutil
import sys
import subprocess
import inspect
import tarfile
import time

import settings

from Leaf.models import ProteinInformation
from Leaf.models import UserCullRequest
from Leaf.models import PDBCullRequest
from Leaf.models import ChainType
from Leaf.models import AllPDBEntries
from Leaf.models import DownloadableFiles

import cullinput.checkPDBinput
import cullinput.controlthread
import cullinput.checkfastaformat

def help_page(request):
    return render_to_response('helppage.html',
                              {},
                              context_instance=RequestContext(request))

def contact(request):
    return render_to_response('contact.html',
                              {},
                              context_instance=RequestContext(request))

def culling(request):
    return render_to_response('culling.html',
                              {},
                              context_instance=RequestContext(request))

def cullingchoice(request):
    choiceMade = request.POST['choiceMade']
    if choiceMade == 'wholePDB':
        return HttpResponseRedirect(reverse('Leaf.views.whole_pdb_culling'))
    elif choiceMade == 'userPDB':
        return HttpResponseRedirect(reverse('Leaf.views.user_pdb_culling'))
    else:
        return HttpResponseRedirect(reverse('Leaf.views.user_culling'))

def user_culling(request):
    return render_to_response('user_culling.html',
                              {},
                              context_instance=RequestContext(request))

def user_pdb_culling(request):
    return render_to_response('user_pdb_culling.html',
                              {},
                              context_instance=RequestContext(request))

def whole_pdb_culling(request):
    return render_to_response('whole_pdb_culling.html',
                              {},
                              context_instance=RequestContext(request))

def downloads(request):
    fileObject = DownloadableFiles.objects.filter(fileName__exact='pdbaa')
    updateTime = time.gmtime(os.path.getmtime(fileObject[0].downloadFile.path))
    updateYear = updateTime[0]
    monthDict = {1 : 'January', 2 : 'February', 3 : 'March', 4 : 'April', 5 : 'May', 6 : 'June', 7 : 'July', 8 : 'August',
                 9 : 'September', 10 : 'October', 11 : 'November', 12 : 'December'}
    updateMonth = monthDict[updateTime[1]]
    daySuffix = {1 : 'st', 2 : 'nd', 3 : 'rd', 21 : 'st', 22 : 'nd', 23 : 'rd', 31 : 'st'}
    updateDay = updateTime[2]
    if daySuffix.has_key(updateDay):
        updateSuffix = daySuffix[updateDay]
    else:
        updateSuffix = 'th'
    updateDay = str(updateDay) + updateSuffix
    modifyDate = updateDay + ' of ' + updateMonth + ' ' + str(updateYear)
    return render_to_response('downloads.html',
                              {'currentPage' : request.path,
                               'modifyDate' : modifyDate,
                               'MEDIA_URL' : settings.MEDIA_URL}
                              )

def index(request):
    return render_to_response('index.html',
                              {})

def results(request, result_id, cullType):
    # Return a 404 if the ID is not in the database
    if cullType == 'user':
        req = UserCullRequest.objects.filter(id__exact=result_id)
    elif cullType == 'pdb':
        req = PDBCullRequest.objects.filter(id__exact=result_id)
    else:
        raise Http404
    
    if len(req) != 1:
        raise Http404
    
    # Determine if the processing is complete
    if not req[0].completed:
        return render_to_response('results.html',
                                  {'notComplete' : True})
    
    if cullType == 'user':
        pdbResult = False
        notWholePDBResult = True
    elif cullType == 'pdb':
        pdbResult = True
        if req[0].wholePDB:
            notWholePDBResult = False
        else:
            notWholePDBResult = True

    outputDict = {'currentPage' : request.path, 'pdb' : pdbResult, 'notWholePDB' : notWholePDBResult}

    outputDict['MEDIA_URL'] = settings.MEDIA_URL
    
    # Return a list of the links to the files otherwise
    return render_to_response('results.html',
                              outputDict)

def sent(request, result_id, cullType):
    server = 'http://www.bioinf.manchester.ac.uk/leaf/'
    if cullType == 'user':
        req = UserCullRequest.objects.filter(id__exact=result_id)
        req = req[0]
        if req.SEG:
            useSEG = ''
        else:
            useSEG = 'Not '
        date = req.requestDate
        resultsFolder = (server + 'results/' + str(date.year) + '/' + str(date.month) + '/' + str(date.day) + '/' +
                         str(result_id) + '/' + cullType)
        
        returnDict = {'percentage' : str(req.sequenceIdentity), 'seg' : useSEG + 'Used', 'resultsLink' : resultsFolder}
        
        if req.minLength == -1:
            returnDict['minlength'] = 'Not Enforced'
        else:
            returnDict['minlength'] = str(req.minLength)
        if req.maxLength == -1:
            returnDict['maxlength'] = 'Not Enforced'
        else:
            returnDict['maxlength'] = str(req.maxLength)
    elif cullType == 'pdb':
        req = PDBCullRequest.objects.filter(id__exact=result_id)
        req = req[0]
        if req.skipNonXray:
            skipNonXray = 'Yes'
        else:
            skipNonXray = 'No'
        if req.skipAlphaCarbon:
            skipAlphaCarbon = 'Yes'
        else:
            skipAlphaCarbon = 'No'
        if req.cullByChain:
            cullMethod = 'Chain'
        else:
            cullMethod = 'Entry'
        date = req.requestDate
        resultsFolder = (server + 'results/' + str(date.year) + '/' + str(date.month) + '/' + str(date.day) + '/' +
                         str(result_id) + '/' + cullType)
        
        returnDict = {'percentage' : str(req.sequenceIdentity), 'RVal' : str(req.maxRValue), 'xray' : skipNonXray,
                      'resolution' : str(req.minResolution) + '-' + str(req.maxResolution), 'alphacarbon' : skipAlphaCarbon,
                      'cullmethod' : cullMethod, 'resultsLink' : resultsFolder}
        
        if req.performIntraEntryCulling:
            returnDict['intraentrycull'] = 'Yes'
            returnDict['intraentrypc'] = str(req.intraEntrySequenceIdentity)
        else:
            returnDict['intraentrycull'] = 'No'
        
        if req.minLength == -1:
            returnDict['minlength'] = 'Not Enforced'
        else:
            returnDict['minlength'] = str(req.minLength)
        if req.maxLength == -1:
            returnDict['maxlength'] = 'Not Enforced'
        else:
            returnDict['maxlength'] = str(req.maxLength)

    returnDict['MEDIA_URL'] = settings.MEDIA_URL

    # Determine if the processing is complete
    if not req.completed:
        returnDict['notComplete'] = True
        return render_to_response('sent.html',
                                  returnDict)
    
    if cullType == 'user':
        pdbResult = False
        notWholePDBResult = True
    elif cullType == 'pdb':
        pdbResult = True
        if req.wholePDB:
            notWholePDBResult = False
        else:
            notWholePDBResult = True

    returnDict['currentPage'] = request.path
    returnDict['pdb'] = pdbResult
    returnDict['notWholePDB'] = notWholePDBResult
    
    # Return a list of the links to the files otherwise
    return render_to_response('sent.html',
                              returnDict)

def user_submit(request):

    errorDict = {}
    
    pastedFasta = request.POST['pastedInfo'].rstrip()
    pastedFasta = pastedFasta.lstrip()
    pastedFasta = pastedFasta.split('\n')
    try:
        uploaded = request.FILES['upload']
    except:
        uploaded = ""
    useSeg = request.POST['SEG']
    email = request.POST['email']
    
    incorrect = []
    
    try:
        percentIdentity = float(request.POST['pc'])
        if percentIdentity < 5 or percentIdentity >= 100:
            errorDict['errorPercent'] = True
            incorrect.append(u"The value entered for the percentage: " + str(percentIdentity) + " is not within the permissible range of values. The percentage must be greater than or equal to 5, and less than 100.")
    except:
        errorDict['errorPercent'] = True
        incorrect.append(u'The value for the sequence identity threshold must be greater than or equal to 5, and less than 100.')
    
    if request.POST['enforceMinLength'] == 'yes':
        try:
            minLength = int(request.POST['minLength'])
            if minLength < 0:
                errorDict['errorMinLength'] = True
                incorrect.append(u'The value for the minimum length must be greater than or equal to 0')
            if request.POST['enforceMaxLength'] == 'yes':
                try:
                    maxLength = int(request.POST['maxLength'])
                    if maxLength < 0:
                        errorDict['errorMaxLength'] = True
                        incorrect.append(u'The value for the maximum length must be greater than or equal to 0')
                    if minLength > maxLength:
                        errorDict['errorMinLength'] = True
                        incorrect.append(u"The value entered for the minimum length: " + str(minLength) + " is greater than the value for the maximum length: " + str(maxLength) + ".")
                except:
                    errorDict['errorMaxLength'] = True
                    incorrect.append(u'The value for the maximum sequence length must be an integer.')
            else:
                maxLength = -1
        except:
            errorDict['errorMinLength'] = True
            incorrect.append(u'The value for the minimum sequence length must be an integer.')
            if request.POST['enforceMaxLength'] == 'yes':
                try:
                    maxLength = int(request.POST['maxLength'])
                    if maxLength < 0:
                        errorDict['errorMaxLength'] = True
                        incorrect.append(u'The value for the maximum length must be greater than or equal to 0')
                    minLength = -1
                except:
                    errorDict['errorMaxLength'] = True
                    incorrect.append(u'The value for the maximum sequence length must be an integer.')
    elif request.POST['enforceMaxLength'] == 'yes':
        try:
            maxLength = int(request.POST['maxLength'])
            if maxLength < 0:
                errorDict['errorMaxLength'] = True
                incorrect.append(u'The value for the maximum length must be greater than or equal to 0')
            minLength = -1
        except:
            errorDict['errorMaxLength'] = True
            incorrect.append(u'The value for the maximum sequence length must be an integer.')
    else:
        minLength = -1
        maxLength = -1
    
    if not re.match("[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]+$", email):
        errorDict['errorEmail'] = True
        incorrect.append(u"The email address entered does not appear to be a valid email address.")

    # Determine whether the user input is correct.
    errorList = [pastedFasta == [''], uploaded == ""]
    if errorList.count(True) == 2:
        errorDict['errorPasted'] = True
        errorDict['errorUploaded'] = True
        incorrect.append(u"You have not provided a source of sequences to cull.")
    elif errorList.count(False) != 1:
        errorDict['errorPasted'] = True
        errorDict['errorUploaded'] = True
        incorrect.append(u'You have provided more than one source of sequences to cull.')
    else:
        # Determine if the sequence input is acceptable
        if pastedFasta != ['']:
            fileString = str(request.POST['pastedInfo'])
            retCode, retVal = cullinput.checkfastaformat.main(fileString, minLength, maxLength)
            if retCode != 0:
                errorDict['errorPasted'] = True
                incorrect.append('There was a problem with the sequences pasted in. ' + retVal)
        elif uploaded != '':
            fileString = ''
            for chunk in uploaded.chunks():
                fileString += chunk
            retCode, retVal = cullinput.checkfastaformat.main(fileString, minLength, maxLength)
            if retCode != 0:
                errorDict['errorUploaded'] = True
                incorrect.append('There was a problem with the uploaded file : ' + str(uploaded.name) + '. ' + retVal)

    if incorrect != []:
        if useSeg == 'yes':
            valueDict = {'errorMessage': incorrect, 'pastedSequences' : request.POST['pastedInfo'], 'upload' : uploaded,
                         'pc' : request.POST['pc'], 'SEG' : 'True', 'email' : request.POST['email']}
        else:
            valueDict = {'errorMessage': incorrect, 'pastedSequences' : request.POST['pastedInfo'], 'upload' : uploaded,
                         'pc' : request.POST['pc'], 'email' : request.POST['email']}
        if request.POST['enforceMinLength'] == 'yes':
            valueDict['enforceMinLength'] = 'yes'
            valueDict['minLength'] = request.POST['minLength']
        if request.POST['enforceMaxLength'] == 'yes':
            valueDict['enforceMaxLength'] = 'yes'
            valueDict['maxLength'] = request.POST['maxLength']
        for i in errorDict.keys():
            valueDict[i] = True
        return render_to_response('user_culling.html',
                                  valueDict,
                                  context_instance=RequestContext(request))

    
    # Setup the database entry
    if useSeg == 'yes':
        useSEG = True
    else:
        useSEG = False
    
#    if uploaded != "":
#        toSave = request.FILES['upload']
#        r = UserCullRequest(
#                    sequences=toSave,
#                    sequenceIdentity=percentIdentity,
#                    SEG=useSEG,
#                    email=str(email),
#                    completed=False,
#                    requestDate=datetime.datetime.now()
#                    )
#        r.save()
#    else:
    r = UserCullRequest(
                sequenceIdentity=percentIdentity,
                minLength = minLength,
                maxLength = maxLength,
                SEG=useSEG,
                email=str(email),
                completed=False,
                requestDate=datetime.datetime.now()
                )
    r.save()
    r.userInput.save('', ContentFile(retVal))
    
#    if r.userInput.size > 1048576:
#        # If the size of the uploaded file is greater than 2Mb then exit
#        if useSeg == 'yes':
#            incorrect = [u"There were errors in your input.",
#                         u"The size of the upload is greater than 2 megabytes."]
#            valueDict = {'errorMessage': incorrect, 'pastedSequences' : request.POST['pastedInfo'],
#                         'pc' : request.POST['pc'], 'SEG' : 'True', 'email' : request.POST['email']
#                         }
#        else:
#            valueDict = {'errorMessage': incorrect, 'pastedSequences' : request.POST['pastedInfo'],
#                         'pc' : request.POST['pc'], 'email' : request.POST['email']
#                         }
#        r.delete()
#        return render_to_response('user_culling.html',
#                                  valueDict,
#                                  context_instance=RequestContext(request))
    
    cullinput.controlthread.RunCull(r, 'seq')
    
#    return HttpResponseRedirect(reverse('Leaf.views.sent'))
    return HttpResponseRedirect(reverse('Leaf.views.sent', kwargs={'result_id': r.id, 'cullType' : 'user'}))

def user_pdb_submit(request):

    errorDict = {}

    pastedChains = request.POST['pastedInfo'].rstrip()
    pastedChains = pastedChains.lstrip()
    pastedChains = pastedChains.split('\n')
    speciesDropBox = request.POST['speciesDropBox']
    speciesTextBox = request.POST['speciesTextBox']
    try:
        uploaded = request.FILES['upload']
    except:
        uploaded = ""
    skipNonXray = request.POST['skipNonXray']
    skipAlphaCarbon = request.POST['skipAlphaCarbon']
    cullMethod = request.POST['cullMethod']
    intraEntryCull = request.POST['intraEntryCull']
    email = request.POST['email']
    
    incorrect = []
    
    try:
        percentIdentity = float(request.POST['pc'])
        if percentIdentity < 5 or percentIdentity >= 100:
            errorDict['errorPercent'] = True
            incorrect.append(u"The value entered for the percentage: " + str(percentIdentity) + " is not within the permissible range of values. The percentage must be greater than or equal to 5, and less than 100.")
    except:
        errorDict['errorPercent'] = True
        incorrect.append(u'The value for the sequence identity threshold must be greater than or equal to 5, and less than 100.')
    
    try:
        minResolution = float(request.POST['minRes'])
        if minResolution < 0:
            errorDict['errorMinRes'] = True
            incorrect.append(u"The value entered for the minimum resolution: " + str(minResolution) + " is less than 0.")
        if minResolution > 100:
            errorDict['errorMinRes'] = True
            incorrect.append(u"The value entered for the minimum resolution: " + str(minResolution) + " is greater than 100.")
        try:
            maxResolution = float(request.POST['maxRes'])
            if maxResolution < 0:
                errorDict['errorMaxRes'] = True
                incorrect.append(u"The value entered for the maximum resolution: " + str(maxResolution) + " is less than 0.")
            if maxResolution > 100:
                errorDict['errorMaxRes'] = True
                incorrect.append(u"The value entered for the maximum resolution: " + str(maxResolution) + " is greater than 100.")
            if minResolution > maxResolution:
                errorDict['errorMinRes'] = True
                incorrect.append(u"The value entered for the minimum resolution: " + str(minResolution) + " is greater than the value for the maximum resolution: " + str(maxResolution) + ".")
        except:
            errorDict['errorMaxRes'] = True
            incorrect.append(u'The value for the maximum resolution must be a number.')
    except:
        errorDict['errorMinRes'] = True
        incorrect.append(u'The value for the minimum resolution must be a number.')
        try:
            maxResolution = float(request.POST['maxRes'])
            if maxResolution < 0:
                errorDict['errorMaxRes'] = True
                incorrect.append(u"The value entered for the maximum resolution: " + str(maxResolution) + " is less than 0.")
            if maxResolution > 100:
                errorDict['errorMaxRes'] = True
                incorrect.append(u"The value entered for the maximum resolution: " + str(maxResolution) + " is greater than 100.")
        except:
            errorDict['errorMaxRes'] = True
            incorrect.append(u'The value for the maximum resolution must be a number.')
    
    try:
        maxRVal = float(request.POST['maxRVal'])
        if maxRVal < 0 or maxRVal > 1:
            errorDict['errorRVal'] = True
            incorrect.append(u"The value entered for the maximum R value: " + str(maxRVal) + " is not within the permissible range of values. The range allowed is: 0 to 1.")
    except:
        errorDict['errorRVal'] = True
        incorrect.append(u'The maximum R value must be a number in the range 0 - 1.')
    
    if request.POST['enforceMinLength'] == 'yes':
        try:
            minLength = int(request.POST['minLength'])
            if minLength < 0:
                errorDict['errorMinLength'] = True
                incorrect.append(u'The value for the minimum length must be greater than or equal to 0')
            if request.POST['enforceMaxLength'] == 'yes':
                try:
                    maxLength = int(request.POST['maxLength'])
                    if maxLength < 0:
                        errorDict['errorMaxLength'] = True
                        incorrect.append(u'The value for the maximum length must be greater than or equal to 0')
                    if minLength > maxLength:
                        errorDict['errorMinLength'] = True
                        incorrect.append(u"The value entered for the minimum length: " + str(minLength) + " is greater than the value for the maximum length: " + str(maxLength) + ".")
                except:
                    errorDict['errorMaxLength'] = True
                    incorrect.append(u'The value for the maximum sequence length must be an integer.')
            else:
                maxLength = -1
        except:
            errorDict['errorMinLength'] = True
            incorrect.append(u'The value for the minimum sequence length must be an integer.')
            if request.POST['enforceMaxLength'] == 'yes':
                try:
                    maxLength = int(request.POST['maxLength'])
                    if maxLength < 0:
                        errorDict['errorMaxLength'] = True
                        incorrect.append(u'The value for the maximum length must be greater than or equal to 0')
                    minLength = -1
                except:
                    errorDict['errorMaxLength'] = True
                    incorrect.append(u'The value for the maximum sequence length must be an integer.')
    elif request.POST['enforceMaxLength'] == 'yes':
        try:
            maxLength = int(request.POST['maxLength'])
            if maxLength < 0:
                errorDict['errorMaxLength'] = True
                incorrect.append(u'The value for the maximum length must be greater than or equal to 0')
            minLength = -1
        except:
            errorDict['errorMaxLength'] = True
            incorrect.append(u'The value for the maximum sequence length must be an integer.')
    else:
        minLength = -1
        maxLength = -1
    
    if not re.match("[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]+$", email):
        errorDict['errorEmail'] = True
        incorrect.append(u"The email address entered does not appear to be a valid email address.")

    # Determine whether the user input is correct.
    errorList = [pastedChains == [''], uploaded == "", speciesDropBox == 'Nothing', speciesTextBox == '']
    if errorList.count(True) == 4:
        errorDict['errorPasted'] = True
        errorDict['errorUploaded'] = True
        errorDict['errorSpeciesDrop'] = True
        errorDict['errorSpeciesText'] = True
        incorrect.append(u"You have not provided a source of chains/entries to cull.")
    elif errorList.count(False) != 1:
        if errorList[0] == False:
            errorDict['errorPasted'] = True
        if errorList[1] == False:
            errorDict['errorUploaded'] = True
        if errorList[2] == False:
            errorDict['errorSpeciesDrop'] = True
        if errorList[3] == False:
            errorDict['errorSpeciesText'] = True
        incorrect.append(u'You have provided more than one source of chains/entries to cull.')
    else:
        # Determine if the PDB chains/entries input are acceptable.
        if errorList[0] == False:
            fileString = str(request.POST['pastedInfo'])
            inputList = fileString.split()
            chainInfo = ChainType.objects.all()
            if cullMethod == 'chain':
                # If the user has supplied any entries in the list of input chains, then extract all chains corresponding to
                # the input entry.
                processedInputList = set([])
                for i in inputList:
                    if len(i) == 4:
                        # Possibly an entry.
                        chainsFromEntry = ProteinInformation.objects.filter(entry__exact=i)
                        if len(chainsFromEntry) == 0:
                            processedInputList.add(i)
                        else:
                            for j in chainsFromEntry:
                                processedInputList.add(str(j.chain))
                    else:
                        processedInputList.add(i)
                processedInputList = list(processedInputList)
                chainDict = {}
                for i in chainInfo:
                    chainDict[i.chain] = i.chainType
                retCode, retVal = cullinput.checkPDBinput.main(processedInputList, allChains=chainDict,  checkType='chain')
            else:
                allEntries = AllPDBEntries.objects.all().values_list('entry', flat=True)
                allProtEntries = set([i.chain[:4] for i in chainInfo if i.chainType == 'Protein'])
                retCode, retVal = cullinput.checkPDBinput.main(inputList, allEntries=allEntries, allProtEntries=allProtEntries, checkType='entry')
            if retCode != 0:
                errorDict['errorPasted'] = True
                incorrect.append(retVal)
        elif uploaded != '':
            fileString = ''
            for chunk in uploaded.chunks():
                fileString += chunk
            inputList = fileString.split()
            inputList = [i.strip() for i in inputList]
            chainInfo = ChainType.objects.all()
            if cullMethod == 'chain':
                # If the user has supplied any entries in the list of input chains, then extract all chains corresponding to
                # the input entry.
                processedInputList = set([])
                for i in inputList:
                    if len(i) == 4:
                        # Possibly an entry.
                        chainsFromEntry = ProteinInformation.objects.filter(entry__exact=i)
                        if len(chainsFromEntry) == 0:
                            processedInputList.add(i)
                        else:
                            for j in chainsFromEntry:
                                processedInputList.add(str(j.chain))
                    else:
                        processedInputList.add(i)
                processedInputList = list(processedInputList)
                chainDict = {}
                for i in chainInfo:
                    chainDict[i.chain] = i.chainType
                retCode, retVal = cullinput.checkPDBinput.main(processedInputList, allChains=chainDict,  checkType='chain')
            else:
                allEntries = AllPDBEntries.objects.all().values_list('entry', flat=True)
                allProtEntries = set([i.chain[:4] for i in chainInfo if i.chainType == 'Protein'])
                retCode, retVal = cullinput.checkPDBinput.main(inputList, allEntries=allEntries, allProtEntries=allProtEntries, checkType='entry')
            if retCode != 0:
                errorDict['errorUploaded'] = True
                incorrect.append(retVal)
        elif speciesDropBox != 'Nothing':
            allChains = ProteinInformation.objects.filter(organism__iexact=speciesDropBox)
            if cullMethod == 'chain':
                allChains = [str(i.chain) for i in allChains]
                if len(allChains) < 2:
                    errorDict['errorSpeciesDrop'] = True
                    incorrect.append(u'There are less than 2 chains in the PDB from the organism you selected.')
                retVal = '\n'.join(allChains)
            else:
                allEntries = list(set([str(i.entry) for i in allChains]))
                if len(allEntries) < 2:
                    errorDict['errorSpeciesDrop'] = True
                    incorrect.append(u'There are less than 2 entries in the PDB from the organism you selected.')
                retVal = '\n'.join(allEntries)
        elif speciesTextBox != '':
            allChains = ProteinInformation.objects.filter(organism__iexact=speciesTextBox)
            if cullMethod == 'chain':
                allChains = [i.chain for i in allChains]
                if len(allChains) < 2:
                    errorDict['errorSpeciesText'] = True
                    incorrect.append(u'There are less than 2 chains in the PDB from the organism you entered.' +
                                     ' This is possibly due to the species being spelled incorrectly.')
                retVal = '\n'.join(allChains)
            else:
                allEntries = list(set([i.entry for i in allChains]))
                if len(allEntries) < 2:
                    errorDict['errorSpeciesText'] = True
                    incorrect.append(u'There are less than 2 entries in the PDB from the organism you entered.' +
                                     ' This is possibly due to the species being spelled incorrectly.')
                retVal = '\n'.join(allEntries)

    if cullMethod == 'entry' and intraEntryCull == 'yes':
        try:
            intraEntryPC = float(request.POST['intraEntryPC'])
            if intraEntryPC < 5 or intraEntryPC >= 100:
                errorDict['errorIntraEntry'] = True
                incorrect.append(u"The value entered for the within entry sequence identity threshold: " + str(intraEntryPC) + " is not within the permissible range of values. The value must be a number greater than or equal to 5, and less than 100.")
        except:
            errorDict['errorIntraEntry'] = True
            incorrect.append(u'The value for the within entry sequence identity threshold must be a number greater than or equal to 5, and less than 100.')
    elif intraEntryCull == 'yes':
        incorrect.append(u'Intra entry culling is selected, but cull by entry is not.')
        

    if incorrect != []:
        valueDict = {'errorMessage': incorrect, 'pastedChains' : request.POST['pastedInfo'],
                     'speciesDropBox' : request.POST['speciesDropBox'], 'speciesTextBox' : request.POST['speciesTextBox'],
                     'pc' : request.POST['pc'], 'minRes' : request.POST['minRes'], 'maxRes' : request.POST['maxRes'],
                     'maxRVal' : request.POST['maxRVal'], 'intraEntryPC' : request.POST['intraEntryPC'],
                     'email' : request.POST['email']
                     }
        if skipNonXray == 'no':
            valueDict['skipNonXray'] = request.POST['skipNonXray']
        if skipAlphaCarbon == 'no':
            valueDict['skipAlphaCarbon'] = request.POST['skipAlphaCarbon']
        if cullMethod == 'entry':
            valueDict['cullMethod'] = request.POST['cullMethod']
        if intraEntryCull == 'yes':
            valueDict['intraEntryCull'] = request.POST['intraEntryCull']
        if request.POST['enforceMinLength'] == 'yes':
            valueDict['enforceMinLength'] = 'yes'
            valueDict['minLength'] = request.POST['minLength']
        if request.POST['enforceMaxLength'] == 'yes':
            valueDict['enforceMaxLength'] = 'yes'
            valueDict['maxLength'] = request.POST['maxLength']
        for i in errorDict.keys():
            valueDict[i] = True
        
        return render_to_response('user_pdb_culling.html',
                                  valueDict,
                                  context_instance=RequestContext(request))
    
    # Setup the database entry
    if skipNonXray == 'yes':
        skipNonXray = True
    else:
        skipNonXray = False
    if skipAlphaCarbon == 'yes':
        skipAlphaCarbon = True
    else:
        skipAlphaCarbon = False
    if cullMethod == 'chain':
        cullMethod = True
    else:
        cullMethod = False
    if intraEntryCull == 'yes':
        intraEntryCull = True
    else:
        intraEntryCull = False
    r = PDBCullRequest(
                       wholePDB=False,
                       sequenceIdentity=percentIdentity,
                       minResolution=minResolution,
                       maxResolution=maxResolution,
                       maxRValue=maxRVal,
                       minLength = minLength,
                       maxLength = maxLength,
                       skipNonXray=skipNonXray,
                       skipAlphaCarbon=skipAlphaCarbon,
                       cullByChain=cullMethod,
                       performIntraEntryCulling=intraEntryCull,
                       intraEntrySequenceIdentity=float(request.POST['intraEntryPC']),
                       email=str(email),
                       completed=False,
                       requestDate=datetime.datetime.now()
                       )
    r.save()
    r.userInput.save('', ContentFile(retVal))

    cullinput.controlthread.RunCull(r, 'pdb')

#    return HttpResponseRedirect(reverse('Leaf.views.sent'))
    return HttpResponseRedirect(reverse('Leaf.views.sent', kwargs={'result_id': r.id, 'cullType' : 'pdb'}))

def whole_pdb_submit(request):

    errorDict = {}

    skipNonXray = request.POST['skipNonXray']
    skipAlphaCarbon = request.POST['skipAlphaCarbon']
    cullMethod = request.POST['cullMethod']
    intraEntryCull = request.POST['intraEntryCull']
    email = request.POST['email']
    
    incorrect = []
    
    try:
        percentIdentity = float(request.POST['pc'])
        if percentIdentity < 5 or percentIdentity >= 100:
            errorDict['errorPercent'] = True
            incorrect.append(u"The value entered for the percentage: " + str(percentIdentity) + " is not within the permissible range of values. The percentage must be greater than or equal to 5, and less than 100.")
    except:
        errorDict['errorPercent'] = True
        incorrect.append(u'The value for the sequence identity threshold must be greater than or equal to 5, and less than 100')
    
    try:
        minResolution = float(request.POST['minRes'])
        if minResolution < 0:
            errorDict['errorMinRes'] = True
            incorrect.append(u"The value entered for the minimum resolution: " + str(minResolution) + " is less than 0.")
        if minResolution > 100:
            errorDict['errorMinRes'] = True
            incorrect.append(u"The value entered for the minimum resolution: " + str(minResolution) + " is greater than 100.")
        try:
            maxResolution = float(request.POST['maxRes'])
            if maxResolution < 0:
                errorDict['errorMaxRes'] = True
                incorrect.append(u"The value entered for the maximum resolution: " + str(maxResolution) + " is less than 0.")
            if maxResolution > 100:
                errorDict['errorMaxRes'] = True
                incorrect.append(u"The value entered for the maximum resolution: " + str(maxResolution) + " is greater than 100.")
            if minResolution > maxResolution:
                errorDict['errorMinRes'] = True
                incorrect.append(u"The value entered for the minimum resolution: " + str(minResolution) + " is greater than the value for the maximum resolution: " + str(maxResolution) + ".")
        except:
            errorDict['errorMaxRes'] = True
            incorrect.append(u'The value for the maximum resolution must be a number.')
    except:
        errorDict['errorMinRes'] = True
        incorrect.append(u'The value for the minimum resolution must be a number.')
        try:
            maxResolution = float(request.POST['maxRes'])
            if maxResolution < 0:
                errorDict['errorMaxRes'] = True
                incorrect.append(u"The value entered for the maximum resolution: " + str(maxResolution) + " is less than 0.")
            if maxResolution > 100:
                errorDict['errorMaxRes'] = True
                incorrect.append(u"The value entered for the maximum resolution: " + str(maxResolution) + " is greater than 100.")
        except:
            errorDict['errorMaxRes'] = True
            incorrect.append(u'The value for the maximum resolution must be a number.')
    
    try:
        maxRVal = float(request.POST['maxRVal'])
        if maxRVal < 0 or maxRVal > 1:
            errorDict['errorRVal'] = True
            incorrect.append(u"The value entered for the maximum R value: " + str(maxRVal) + " is not within the permissible range of values. The range allowed is: 0 to 1.")
    except:
        errorDict['errorRVal'] = True
        incorrect.append(u'The value for the R value must be a number in the range 0 - 1.')
    
    if request.POST['enforceMinLength'] == 'yes':
        try:
            minLength = int(request.POST['minLength'])
            if minLength < 0:
                errorDict['errorMinLength'] = True
                incorrect.append(u'The value for the minimum length must be greater than or equal to 0')
            if request.POST['enforceMaxLength'] == 'yes':
                try:
                    maxLength = int(request.POST['maxLength'])
                    if maxLength < 0:
                        errorDict['errorMaxLength'] = True
                        incorrect.append(u'The value for the maximum length must be greater than or equal to 0')
                    if minLength > maxLength:
                        errorDict['errorMinLength'] = True
                        incorrect.append(u"The value entered for the minimum length: " + str(minLength) + " is greater than the value for the maximum length: " + str(maxLength) + ".")
                except:
                    errorDict['errorMaxLength'] = True
                    incorrect.append(u'The value for the maximum sequence length must be an integer.')
            else:
                maxLength = -1
        except:
            errorDict['errorMinLength'] = True
            incorrect.append(u'The value for the minimum sequence length must be an integer.')
            if request.POST['enforceMaxLength'] == 'yes':
                try:
                    maxLength = int(request.POST['maxLength'])
                    if maxLength < 0:
                        errorDict['errorMaxLength'] = True
                        incorrect.append(u'The value for the maximum length must be greater than or equal to 0')
                    minLength = -1
                except:
                    errorDict['errorMaxLength'] = True
                    incorrect.append(u'The value for the maximum sequence length must be an integer.')
    elif request.POST['enforceMaxLength'] == 'yes':
        try:
            maxLength = int(request.POST['maxLength'])
            if maxLength < 0:
                errorDict['errorMaxLength'] = True
                incorrect.append(u'The value for the maximum length must be greater than or equal to 0')
            minLength = -1
        except:
            errorDict['errorMaxLength'] = True
            incorrect.append(u'The value for the maximum sequence length must be an integer.')
    else:
        minLength = -1
        maxLength = -1

    if not re.match("[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]+$", email):
        errorDict['errorEmail'] = True
        incorrect.append(u"The email address entered does not appear to be a valid email address.")

    if cullMethod == 'entry' and intraEntryCull == 'yes':
        try:
            intraEntryPC = float(request.POST['intraEntryPC'])
            if intraEntryPC < 5 or intraEntryPC >= 100:
                errorDict['errorIntraEntry'] = True
                incorrect.append(u"The value entered for the within entry sequence identity threshold: " + str(intraEntryPC) + " is not within the permissible range of values. The value must be a number greater than or equal to 5, and less than 100.")
        except:
            errorDict['errorIntraEntry'] = True
            incorrect.append(u'The value for the within entry sequence identity threshold must be a number greater than or equal to 5, and less than 100.')
    elif intraEntryCull == 'yes':
        incorrect.append(u'Intra entry culling is selected, but cull by entry is not.')

    if incorrect != []:
        valueDict = {'errorMessage': incorrect, 'pc' : request.POST['pc'], 'minRes' : request.POST['minRes'],
                     'maxRes' : request.POST['maxRes'], 'maxRVal' : request.POST['maxRVal'],
                     'intraEntryPC' : request.POST['intraEntryPC'], 'email' : request.POST['email']}
        if skipNonXray == 'no':
            valueDict['skipNonXray'] = request.POST['skipNonXray']
        if skipAlphaCarbon == 'no':
            valueDict['skipAlphaCarbon'] = request.POST['skipAlphaCarbon']
        if cullMethod == 'entry':
            valueDict['cullMethod'] = request.POST['cullMethod']
        if intraEntryCull == 'yes':
            valueDict['intraEntryCull'] = request.POST['intraEntryCull']
        if request.POST['enforceMinLength'] == 'yes':
            valueDict['enforceMinLength'] = 'yes'
            valueDict['minLength'] = request.POST['minLength']
        if request.POST['enforceMaxLength'] == 'yes':
            valueDict['enforceMaxLength'] = 'yes'
            valueDict['maxLength'] = request.POST['maxLength']
        for i in errorDict.keys():
            valueDict[i] = True
        
        return render_to_response('whole_pdb_culling.html',
                                  valueDict,
                                  context_instance=RequestContext(request))
    
    # Setup the database entry
    if skipNonXray == 'yes':
        skipNonXray = True
    else:
        skipNonXray = False
    if skipAlphaCarbon == 'yes':
        skipAlphaCarbon = True
    else:
        skipAlphaCarbon = False
    if cullMethod == 'chain':
        cullMethod = True
    else:
        cullMethod = False
    if intraEntryCull == 'yes':
        intraEntryCull = True
    else:
        intraEntryCull = False
    r = PDBCullRequest(
                       wholePDB=True,
                       sequenceIdentity=percentIdentity,
                       minResolution=minResolution,
                       maxResolution=maxResolution,
                       maxRValue=maxRVal,
                       minLength = minLength,
                       maxLength = maxLength,
                       skipNonXray=skipNonXray,
                       skipAlphaCarbon=skipAlphaCarbon,
                       cullByChain=cullMethod,
                       performIntraEntryCulling=intraEntryCull,
                       intraEntrySequenceIdentity=float(request.POST['intraEntryPC']),
                       email=str(email),
                       completed=False,
                       requestDate=datetime.datetime.now()
                       )
    r.save()
    
    cullinput.controlthread.RunCull(r, 'pdb')

#    return HttpResponseRedirect(reverse('Leaf.views.sent'))
    return HttpResponseRedirect(reverse('Leaf.views.sent', kwargs={'result_id': r.id, 'cullType' : 'pdb'}))

def txtlist(request, result_id, cullType, fileName):
    # Return a 404 if the ID is not in the database
    if cullType == 'user':
        req = UserCullRequest.objects.filter(id__exact=result_id)
    elif cullType == 'pdb':
        req = PDBCullRequest.objects.filter(id__exact=result_id)
    else:
        raise Http404
    
    # Determine if the processing is complete
    if len(req) != 1 or not req[0].completed:
        raise Http404
    
    req = req[0]
    
    if fileName == 'Input':
        file = req.userInput
    elif fileName == 'Removed':
        file = req.removed
    elif fileName == 'NonRedundantList':
        file = req.nonredNoSeq
    elif fileName == 'NonRedundantFasta':
        file = req.nonredSeq
#    elif fileName == 'Similarities':
#        file = req.pairwiseScores
    else:
        raise Http404
    
    # Display a text file with one line on each line
    response = HttpResponse(file, content_type='text/plain')
#    response['Content-Disposition'] = 'attachment; filename=' + fileName + '.txt'
    return response

def download_gzipped(request, fileName):
    if fileName == 'WindowsSourceCode':
        tempFile = open(settings.MEDIA_ROOT + 'TarData/WindowsLocalLeaf.tar.gz', 'r')
        outputFile = File(tempFile)
        response = HttpResponse(outputFile, content_type='application/gzip')
        response['Content-Disposition'] = 'attachment; filename=LocalLeaf.tar.gz'
        return response
    elif fileName == 'LinuxSourceCode':
        tempFile = open(settings.MEDIA_ROOT + 'TarData/LinuxLocalLeaf.tar.gz', 'r')
        outputFile = File(tempFile)
        response = HttpResponse(outputFile, content_type='application/gzip')
        response['Content-Disposition'] = 'attachment; filename=LocalLeaf.tar.gz'
        return response
    elif fileName == 'StandaloneFiles':
        tempFile = open(settings.MEDIA_ROOT + 'TarData/PDBData.tar.gz', 'r')
        outputFile = File(tempFile)
        response = HttpResponse(outputFile, content_type='application/gzip')
        response['Content-Disposition'] = 'attachment; filename=PDBData.tar.gz'
        return response
    else:
        fileObject = DownloadableFiles.objects.filter(fileName__exact=fileName)
        if fileObject:
            response = HttpResponse(fileObject[0].downloadFile, content_type='application/gzip')
            response['Content-Disposition'] = 'attachment; filename=' + fileName + '.gz'
            return response
        else:
            chunks = fileName.split('/')
            if len(chunks) == 2:
                tempFile = open(settings.MEDIA_ROOT + 'ModelOrganisms/' + fileName, 'r')
                outputFile = File(tempFile)
                if chunks[-1][-3:] == '.gz':
                    response = HttpResponse(outputFile, content_type='application/gzip')
                    response['Content-Disposition'] = 'attachment; filename=' + chunks[-1]
                else:
                    response = HttpResponse(outputFile, content_type='text/plain')
                    response['Content-Disposition'] = 'filename=' + chunks[0] + '_' + chunks[1]
                return response
            else:
                raise Http404