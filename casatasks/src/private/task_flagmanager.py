from __future__ import absolute_import
import os
import time
import copy

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatools import agentflagger
    from casatasks import casalog
else:
    from taskinit import *

    agentflagger = casac.agentflagger


def flagmanager(
    vis=None,
    mode=None,
    versionname=None,
    oldname=None,
    comment=None,
    merge=None,
    ):

    casalog.origin('flagmanager')
    aflocal = agentflagger()

    try:
        if type(vis) == str and os.path.exists(vis):
            if mode != 'rename':
                aflocal.open(vis)
        else:
            raise ValueError('Visibility data set not found - please verify the name')
        if mode == 'list':
            flist = []
            flist = aflocal.getflagversionlist()
            
            # Get the name of the MS and properly add it to the dictionary
            if "\nMS : " in flist[0]:
                MS = flist.pop(0)
                MS = MS.strip("\nMS : ")
                
            flist.remove('main : working copy in main table')
            fdict = dict(enumerate(flist))
            fversionsdict = copy.deepcopy(fdict)
            for k in fdict:
                singleversion = {}
                # split each flagversion into versionname and comment
                # The below partitioning is a big fragile. If the string contains
                # other entries of the character ':', the spliting will fail
                (versionname, middle, comment) = fdict[k].partition(':')
                singleversion['name'] = versionname.rstrip()
                singleversion['comment'] = comment.lstrip()
                fversionsdict[k] = singleversion
            
            fversionsdict['MS'] = MS
            return fversionsdict
            
        elif mode == 'save':
            if versionname == '':
                raise ValueError("Illegal empty versionname: ''")
            
            newdir = vis+'.flagversions/flags.'+versionname
            if os.path.exists(newdir):
                tt = int(time.time())
                tmpname = versionname+'.old.'+str(tt)
                casalog.post("Version name \'%s\' already exist. Will rename it to %s"%(versionname,tmpname), 'WARN')
                
                tmpdir = newdir+'.old.'+str(tt)
                # Rename existing versionname to old name
                os.rename(newdir, tmpdir)
    
                # Edit entry in .flagversions/FLAG_VERSION_LIST
                # For realistic usecases, this file is short enough to keep in memory
                file = vis + '.flagversions/FLAG_VERSION_LIST'
                fd = open(file)
                lines = fd.readlines()
                fd.close()
    
                for i in range(len(lines)):
                    if (lines[i])[:len(versionname) + 3] == versionname + ' : ':
                        lines[i] = tmpname + ' : ' + comment + '\n'
                        break
    
                fd = open(file, 'w')
                fd.writelines(lines)
                fd.close()                                
            
            casalog.post('Save current flagversions to ' + versionname)
            aflocal.saveflagversion(versionname=versionname,
                                    comment=comment, merge=merge)
        elif mode == 'restore':
            if versionname == '':
                raise ValueError("Illegal versionname: ''")
            
            casalog.post('Restore flagversions ' + versionname)
            aflocal.restoreflagversion(versionname=versionname,
                    merge=merge)
        elif mode == 'delete':
            if versionname == '':
                raise ValueError("Illegal versionname: ''")
            
            aflocal.deleteflagversion(versionname=versionname)
        elif mode == 'rename':
            if versionname == '':
                raise ValueError("Illegal versionname: ''")

            if oldname == '':
                raise ValueError("Illegal oldname: ''")

            # The directory structure is unlikely to change
            olddir = vis + '.flagversions/flags.' + oldname
            newdir = vis + '.flagversions/flags.' + versionname
            if not os.path.isdir(olddir):
                raise ValueError('No such flagversions: ' + str(oldname))
            
            if os.path.exists(newdir):
                raise ValueError('Flagversions ' + str(versionname) + ' already exists!')

            casalog.post('Rename flagversions "%s" to "%s"' % (oldname,versionname))

            os.rename(olddir, newdir)

            # Edit entry in .flagversions/FLAG_VERSION_LIST
            # For realistic usecases, this file is short enough to keep in memory
            file = vis + '.flagversions/FLAG_VERSION_LIST'
            fd = open(file)
            lines = fd.readlines()
            fd.close()

            for i in range(len(lines)):
                if (lines[i])[:len(oldname) + 3] == oldname + ' : ':
                    lines[i] = versionname + ' : ' + comment + '\n'
                    break

            fd = open(file, 'w')
            fd.writelines(lines)
            fd.close()
            
        else:
            raise ValueError('Unknown mode' + str(mode))
        
    finally:
        aflocal.done()
