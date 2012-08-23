from django.core.management.base import BaseCommand

import datetime
import sys
import traceback
import time

from Leaf.models import UserCullRequest
from Leaf.models import PDBCullRequest

class Command(BaseCommand):

    def handle(self, *args, **options):
        try:
            maximumDaysForStorage = 8

            allEntries = UserCullRequest.objects.all()
            for i in range(len(allEntries)):
                check = allEntries[i]
                elapsed = datetime.datetime.now() - check.requestDate
                daysElapsed = elapsed.days
                if daysElapsed >= maximumDaysForStorage:
                    try:
                        check.delete()
                    except:
                        raise
        
            allEntries = PDBCullRequest.objects.all()
            for i in range(len(allEntries)):
                check = allEntries[i]
                elapsed = datetime.datetime.now() - check.requestDate
                daysElapsed = elapsed.days
                if daysElapsed >= maximumDaysForStorage:
                    try:
                        check.delete()
                    except:
                        raise
        except:
            logger = open('/srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp/ErrorLogs/REMOVAL.log', 'a')
            logger.write('\tERROR encountered on ' + time.strftime('%Y/%m/%d/ at %H:%M:%S', time.gmtime()) + ':\n')
            excType, excValue, excTrace = sys.exc_info()
            logger.write('\t\tException type: ' + str(excType) + '\n')
            logger.write('\t\tException Value: ' + str(excValue) + '\n')
            errors = traceback.format_exception(excType, excValue, excTrace)
            for i in errors:
                logger.write('\t\t' + i)
            logger.close()