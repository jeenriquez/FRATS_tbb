#! /usr/bin/env python
'''This program is used to delete FRATS files in other locus nodes but locus013.'''
from optparse import OptionParser
import os
import time
import subprocess
import socket

# Need to have
'''
    - Need to check first which files are in locus013
     - Need to delete those files from other nodes.
        - need to check if those files are in locus013, by checking the name and the size, both need to match, otherwise raise an error.

'''
raise NotImplementedError



def converttime(sec):
    return time.strftime("D%Y%m%dT%H%M%S",time.gmtime(sec))

def filename2location(filename):
    name=os.path.split(filename)[-1]
    mainpath="/data/TBB/FRATS/tbb/data/"
    assert len(name)==47 # Incorrect filesize length
    newpath=mainpath+name[0:22]
    return os.path.join(newpath,name)

def rsyncfile(filename_from,filename_to=None,node_from="locus013",node_to="locus013"):
    """ """
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


#------------------------------------------------------
#Command line options
#------------------------------------------------------
parser = OptionParser()
parser.add_option("-l","--locus",type="str",default='',help="Locus node from which you want to copy data.")
(options, args) = parser.parse_args()

if not parser.get_prog_name()=="listfratsfiles.py":
    #   Program was run from within python
    locus = ""
else:
    locus = options.locus

node_local = socket.gethostname()

if locus == '':
    locus = node_local

if locus != node_local and node_local != 'locus013':
    print '**Warning**: Using '+node_local+' as locus instead of '+locus+', since node_local is not locus013.'
    locus = node_local

#------------------------------------------------------

#
file_std_loc= "/data/TBB/VHECRtest/"
#triggers=open("/home/veen/logs/LORAreceived").readlines()
#fratstriggers=[t for t in triggers if "FRATS" in t and "not allowed" not in t and "allowed" in t]
#fratstimestamps=[t.split()[0] for t in fratstriggers]
#fratstimestring=[converttime(int(t))[0:-2] for t in fratstimestamps]

#Finding all Gb size files (To not list the VHECR files).Not that there are some FRATS files which are smaller, but must likely need to be deleted.
if locus:
    pp = subprocess.Popen(["ssh",locus, "ls","-lh",file_std_loc,"|","grep",size_str], stdout=subprocess.PIPE)
    out, err = pp.communicate()
    out=out.split('\n')[:-1]
    files = [file_std_loc+o.split()[-1] for o in out]
else:
   files=[file_std_loc+f for f in os.listdir(file_std_loc) if os.path.getsize(file_std_loc+f) > size_flt]

fratsfiles=[f for t in fratstimestring for f in files if t in f]


for f in fratsfiles:
#    print f,filename2location(f)2
    rsyncfile(f,node_from=locus)
