#! /usr/bin/env python
'''This program is used to rsync FRATS files to a sigle locus node (locus013).

.. moduleauthor:: J. Emilio Enriquez <e.enriquez 'at' astro.ru.nl>

Revision History:
V1.0 created by E. Enriquez, Jan 2014

'''


from optparse import OptionParser
import os
import time
import subprocess
import socket
import smtplib
from email.mime.text import MIMEText


def converttime(sec):
    '''Convers times in seconds to D+date+T+time'''
    return time.strftime("D%Y%m%dT%H%M%S",time.gmtime(sec))

def filename2location(filename):
    ''' This function joins the filename with the location  '''
    name=os.path.split(filename)[-1]
    mainpath="/data/TBB/FRATS/tbb/data/"
    assert len(name)==47 # Incorrect filesize length
    newpath=mainpath+name[0:22]
    return os.path.join(newpath,name)

def rsyncfile(filename_from,filename_to=None,node_from="locus013",node_to="locus013"):
    """ This function rsyncs files from a list. It does this inteligently towards locus013."""
    if not filename_to:
        filename_to=filename2location(filename_from)
    path_to=os.path.split(filename_to)[0]

    if node_to != "locus013":
        raise NotImplementedError('We assume we move the data to locus013')

    node_local = socket.gethostname()

    if node_local != 'locus013':
        node_from = node_local
        proc=subprocess.call(["ssh",node_to,"mkdir","-pvm","775",path_to])
    else:
        proc = subprocess.call(["mkdir","-pvm","775",path_to])

    if node_from == 'locus013':
        source = filename_from
        destination = filename_to
        print "mv",source,destination
        proc = subprocess.call(["mv","-vi",source,destination])
    else:
        if node_local != 'locus013':
            source = filename_from
            destination =  node_to+":"+filename_to
        else:
            source = node_from+':'+filename_from
            destination = filename_to

        print "rsync","--bwlimit=2500","-rvvhp","--progress","--remove-source-files",source,destination
        proc=subprocess.Popen(["ionice","-c","3","rsync","--bwlimit=2500","-rvvhp","--progress","--remove-source-files",source,destination])


#------------------------------------------------------
#Command line options
#------------------------------------------------------
parser = OptionParser()
parser.add_option("-l","--locus",type="str",default='',help="Locus node from which you want to copy data.")
parser.add_option("-s","--size",type="str",default='G',help="File size to copy (M or G) for >Mb or >Gb sizes.")
parser.add_option("-t","--test",action="store_true",default=False,help="During test will not do rsync")
parser.add_option("-i","--id_obs",type="str",default=None,help="Give an the observation ID, with format D20120718T1215.")

(options, args) = parser.parse_args()

if not parser.get_prog_name()=="rsyncfratsfiles.py":
    #   Program was run from within python
    locus = ""
    size = 'G'
    test = False
    id_obs = None
else:
    locus = options.locus
    size = options.size
    test = options.test
    id_obs = options.id_obs

#------------------------------------------------------


if size != 'M':
    size_str = "'G 2'"
else:
    size_str = "'M 2'"

node_local = socket.gethostname()

if locus == '':
    locus = node_local

if locus != node_local and node_local != 'locus013':
    print '**Warning**: Using '+node_local+' as locus instead of '+locus+', since node_local is not locus013.'
    locus = node_local

if locus == 'superterp':
    locus = ['locus0'+str(i+11).zfill(2) for i in range(6)]
elif locus == 'core':
    locus = ['locus0'+str(i+10).zfill(2) for i in range(35)]
elif locus =='all':
    locus = ['locus0'+str(i+10).zfill(2) for i in range(72)]
else:
    locus = [locus]

#------------------------------------------------------
#File name information from triggers
file_std_loc= "/data/TBB/VHECRtest/"
triggers=open("/home/veen/logs/LORAreceived").readlines()

fratstriggers=[t for t in triggers if "FRATS" in t and "not allowed" not in t and "allowed" in t]
fratstimestamps=[t.split()[0] for t in fratstriggers]
fratstimestring=[converttime(int(t))[0:-2] for t in fratstimestamps]

if id_obs and 'D' in id_obs and 'T' in id_obs:
    fratstimestring = [id_obs]

msg = ''

for loc in locus:
    #Finding all Gb size files (To not list the VHECR files).Not that there are some FRATS files which are smaller, but must likely need to be deleted.
    pp = subprocess.Popen(["ssh",loc, "ls","-lh",file_std_loc,"|","grep",size_str], stdout=subprocess.PIPE)
    out, err = pp.communicate()
    out=out.split('\n')[:-1]
    files = [file_std_loc+o.split()[-1] for o in out]
    fratsfiles=[f for t in fratstimestring for f in files if t in f]

    if len(fratsfiles) < 1:
        print "**Note**: I did't find FRATS files in this location => " +loc+":"+file_std_loc
    else:
        if not msg:
            msg ='\n'
        msg += '\n **'+loc+'**'
        for filename in fratsfiles:
            msg += '\n  ' + filename

    for f in fratsfiles:
    #    print f,filename2location(f)
        rsyncfile(f,node_from=loc)

#------------------------------------------------------
#Mailing if there is data being transfer.
if msg:
    msg = MIMEText('Transfering the following file(s): '+ msg)
    msg['Subject'] = '[FRATS] Notice: rsyncing frats files from lofar locus nodes.'
    msg['From'] = 'Emilio_Enriquez'
    msg['To'] = 'e.enriquez@astro.ru.nl'
    send = smtplib.SMTP('localhost')
    send.sendmail('Emilio_Enriquez', ['e.enriquez@astro.ru.nl'], msg.as_string())
    send.quit()
