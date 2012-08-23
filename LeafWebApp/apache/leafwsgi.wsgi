import os
import sys

sys.path.append('/srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/')
sys.path.append('/srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp')

os.environ['DJANGO_SETTINGS_MODULE'] = 'LeafWebApp.settings'

import django.core.handlers.wsgi
application = django.core.handlers.wsgi.WSGIHandler()
