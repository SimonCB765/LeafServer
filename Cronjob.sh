#!/bin/sh
echo -n "Starting cron job on  " >> /srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp/ErrorLogs/CRON.log
date >> /srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp/ErrorLogs/CRON.log
rsync -rlpt -v -z --delete \
rsync.ebi.ac.uk::pub/databases/rcsb/pdb-remediated/data/structures/divided/mmCIF/ \
/fs/nas15/home/mqbpjsb2/LocalPDB/mmCIF
echo "Rsync finished" >> /srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp/ErrorLogs/CRON.log
python /srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp/manage.py updatelocalPDB /fs/nas15/home/mqbpjsb2/LocalPDB/mmCIF > /fs/nas15/home/mqbpjsb2/LeafCron.log 2>&1
echo "Updating database complete" >> /srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp/ErrorLogs/CRON.log
python /srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp/manage.py generatedownloads > /fs/nas15/home/mqbpjsb2/LeafCron.log 2>&1
echo "Downloads generated" >> /srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp/ErrorLogs/CRON.log
python /srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp/manage.py removeoldentries > /fs/nas15/home/mqbpjsb2/LeafCron.log 2>&1
echo "Old requests removed" >> /srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp/ErrorLogs/CRON.log
mv /fs/nas15/home/mqbpjsb2/LeafBackupCurrent.sql /fs/nas15/home/mqbpjsb2/LeafBackupOld.sql
mysqldump -uleaf -pxF_WT4aamHw3EU7GK7W6W leaf > /fs/nas15/home/mqbpjsb2/LeafBackupCurrent.sql
echo "Database backup created" >> /srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp/ErrorLogs/CRON.log
echo -n "Ending cron job on  " >> /srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp/ErrorLogs/CRON.log
date >> /srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp/ErrorLogs/CRON.log
echo -------------------------------------  >> /srv/www/vhosts.d/www.bioinf/html/doig/cgi-bin/django_projects/LeafWebApp/ErrorLogs/CRON.log
