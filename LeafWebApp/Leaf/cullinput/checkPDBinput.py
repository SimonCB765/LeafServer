def main(listToCheck, allChains={}, allEntries=[], allProtEntries=[], checkType='chain'):

    listToCheck = set([i.strip().upper() for i in listToCheck])
    listToCheck -= set([''])  # Remove the record of any blank lines.

    if checkType == 'chain':
        nonExistantChains = [i for i in listToCheck if i not in allChains.keys()]
        nonProtChains = [i for i in listToCheck if not i in nonExistantChains and allChains[i] != 'Protein']
        if nonExistantChains != [] or nonProtChains != []:
            if nonExistantChains != []:
                errorMessage = u'The following chains were not found in the database:\n' + u', '.join(nonExistantChains) + u'\n'
            if nonProtChains != []:
                errorMessage = u'The following chains are not known to be protein chains:\n' + u', '.join(nonProtChains) + u'\n'
            return 1, errorMessage
        elif len(listToCheck) < 2:
            # If there aren't enough chains in the input file
            errorMessage = u'The input file must contain at least 2 chains.\n'
            return 1, errorMessage
        else:
            return 0, '\n'.join(listToCheck)
    elif checkType == 'entry':
        nonExistantEntries = [i for i in listToCheck if not i in allEntries]
        entriesWithNoProts = [i for i in listToCheck if not i in allProtEntries and not i in nonExistantEntries]

        if nonExistantEntries != [] or entriesWithNoProts != []:
            if nonExistantEntries != []:
                errorMessage = u'The following entries were not found in the database:\n' + u', '.join(nonExistantEntries) + u'\n'
            if entriesWithNoProts != []:
                errorMessage += (u'The following entries were found in the database, but have no protein chains recorded for them:\n' +
                                 u', '.join(entriesWithNoProts) + u'\n')
            return 1, errorMessage
        elif len(listToCheck) < 2:
            # If there aren't enough entries in the input file
            errorMessage = u'The input file must contain at least 2 entries.\n'
            return 1, errorMessage
        else:
            return 0, '\n'.join(listToCheck)
    else:
        return 2, u'The choice for paramter checkType was not one of \'chain\' or \'entry\'.'