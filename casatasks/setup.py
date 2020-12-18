#!/usr/bin/env python3
# Copyright (C) 2018
# Associated Universities, Inc. Washington DC, USA.
#
# This library is free software; you can redistribute it and/or modify it
# under the terms of the GNU Library General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
# License for more details.
#
# You should have received a copy of the GNU Library General Public License
# along with this library; if not, write to the Free Software Foundation,
# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
#
# Correspondence concerning AIPS++ should be addressed as follows:
#        Internet email: aips2-request@nrao.edu.
#        Postal address: AIPS++ Project Office
#                        National Radio Astronomy Observatory
#                        520 Edgemont Road
#                        Charlottesville, VA 22903-2475 USA


"""CASAtasks Python Module

This is a standard python module that provides CASA tools and tasks
without regular CASA's bespoke CLI.
"""
from __future__ import division, print_function

classifiers = """\
Development Status :: 3 - Alpha
Intended Audience :: Developers
License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)
Programming Language :: Python :: 2.7
Programming Language :: C++
Topic :: Software Development
Topic :: Scientific/Engineering :: Astronomy
Topic :: Software Development :: Libraries :: Python Modules
Operating System :: MacOS :: MacOS X
Operating System :: POSIX
"""
##
##  Without this zlib import at the beginning, building the casataks
##  wheel on RHEL7 (w/ rh-python36-python-3.6.3-3.el7.x86_64) dies
##  with a segmentation violation at:
##
##    #0  0x00007f1d192ecf89 in fill_window () from /lib64/libz.so.1
##    #1  0x00007f1d192ed590 in deflate_slow () from /lib64/libz.so.1
##    #2  0x00007f1d192ee624 in deflate () from /lib64/libz.so.1
##    #3  0x00007f1d0851a195 in zlib_Compress_compress () from /opt/rh/rh-python36/root/usr/lib64/python3.6/lib-dynload/zlib.cpython-36m-x86_64-linux-gnu.so
##
import zlib
import sys
import os

try:
    import casatools
    from casatools.config import build as tools_config
except:
    print("cannot find CASAtools (https://open-bitbucket.nrao.edu/projects/CASA/repos/casatools/browse) in PYTHONPATH")
    os._exit(1)

from setuptools import setup, find_packages
from distutils.dir_util import copy_tree, remove_tree
from distutils.util import get_platform
from distutils.cmd import Command
from distutils.command.build import build
from subprocess import Popen, PIPE
from subprocess import call as Proc
from datetime import datetime
from textwrap import dedent
from shutil import copy2
import subprocess
import sysconfig
import pickle
import errno
import time
import re

from os import listdir
from os.path import isfile, join, islink
from itertools import chain
import argparse

parser=argparse.ArgumentParser()
parser.add_argument('--version', help='version')
parser.add_argument('bdist_wheel', help='bdist_wheel')
args=parser.parse_args()

print (args.version)

module_name = 'casatasks'

pyversion = float(sys.version_info[0]) + float(sys.version_info[1]) / 10.0

if pyversion < 3:
    str_encode = str
    str_decode = str
    def pipe_decode(output):
        return output
else:
    def str_encode(s):
        return bytes(s,sys.getdefaultencoding())
    def str_decode(bs):
        return bs.decode(sys.getdefaultencoding(),"strict")
    def pipe_decode(output):
        if isinstance(output,bytes) or isinstance(output,bytearray):
            return str_decode(output)
        elif isinstance(output,tuple):
            return (str_decode(output[0]),str_decode(output[1]))
        else:
            return ("","")

def compute_version( ):
    if (args.version != None ):
        print (args.version.split("."))
        (major, minor, patch, feature) = args.version.split(".")
        return(int(major), int(minor), int(patch), int(feature),"","","")
    else:
        proc = Popen( [ "./version" ], stdout=PIPE, stderr=PIPE )
        out,err = pipe_decode(proc.communicate( ))
        print(out)
        devbranchtag = out.split(" ")[0].strip()
        print(devbranchtag)
        releasetag = out.split(" ")[1].strip()
        dirty=""
        if (len(out.split(" ")) == 3):
            print("Latest commit doesn't have a tag. Adding -dirty flag to version string.")
            dirty="+" + out.split(" ")[2].strip() # "+" denotes local version identifier as described in PEP440
        print(releasetag)
        devbranchversion = ""
        devbranchrevision = ""
        if (devbranchtag != releasetag):
            devbranchrevision = devbranchtag.split("-")[-1]
            if (devbranchtag.startswith("CAS-")):
                devbranchversion=devbranchtag.split("-")[1]
            else:
                devbranchversion=100
            devbranchrevision = devbranchtag.split("-")[-1]
        else:
            isDevBranch = False
        (major, minor, patch, feature) = releasetag.split(".")
        return(int(major), int(minor), int(patch), int(feature), devbranchversion, devbranchrevision, dirty)


(casatasks_major,casatasks_minor,casatasks_patch,casatasks_feature, devbranchversion, devbranchrevision, dirty) = compute_version( )
casatasks_version = '%d.%d.%d.%d%s' % (casatasks_major,casatasks_minor,casatasks_patch,casatasks_feature, dirty)
if devbranchversion !="":
    casatasks_version = '%d.%d.%d.%da%s.dev%s%s' % (casatasks_major,casatasks_minor,casatasks_patch,casatasks_feature,devbranchversion,devbranchrevision,dirty)

public_scripts = [ 'src/scripts/config.py', 'src/scripts/LICENSE.txt', 'src/scripts/__main__.py' ]

private_scripts = [ 'src/scripts/userconfig.py',
                    'src/scripts/casa_transition.py',
                    'src/scripts/ialib.py',
                    'src/scripts/cvt.py',
                    'src/scripts/callibrary.py',
                    'src/scripts/flaghelper.py',
                    'src/scripts/partitionhelper.py',
                    'src/tasks/task_imhead.py',
                    'src/tasks/task_immoments.py',
                    'src/tasks/task_imhistory.py',
                    'src/tasks/task_applycal.py',
                    'src/tasks/task_bandpass.py',
                    'src/tasks/task_blcal.py',
                    'src/tasks/task_calstat.py',
                    'src/tasks/task_concat.py',
                    'src/scripts/concatephem.py',
                    'src/tasks/task_split.py',
                    'src/tasks/task_listobs.py',
                    'src/tasks/task_flagdata.py',
                    'src/tasks/task_flagcmd.py',
                    'src/tasks/task_setjy.py',
                    'src/scripts/setjy_helper.py',
                    'src/scripts/solar_system_setjy.py',
                    'src/scripts/mstools.py',
                    'src/scripts/update_spw.py',
                    'src/tasks/task_cvel.py',
                    'src/tasks/task_cvel2.py',
                    'src/tasks/task_importuvfits.py',
                    'src/tasks/task_importfits.py',
                    'src/tasks/task_exportfits.py',
                    'src/tasks/task_exportuvfits.py',
                    'src/tasks/task_partition.py',
                    'src/tasks/task_listpartition.py',
                    'src/tasks/task_flagmanager.py',
                    'src/tasks/task_mstransform.py',
                    'src/tasks/task_tclean.py',
                    'src/tasks/task_immath.py',
                    'src/tasks/task_vishead.py',
                    'src/scripts/vishead_util.py',
                    'src/tasks/task_uvsub.py',
                    'src/tasks/task_spxfit.py',
                    'src/tasks/task_splattotable.py',
                    'src/tasks/task_specsmooth.py',
                    'src/tasks/task_specflux.py',
                    'src/tasks/task_smoothcal.py',
                    'src/tasks/task_specfit.py',
                    'src/tasks/task_imstat.py',
                    'src/tasks/task_slsearch.py',
                    'src/tasks/task_delmod.py',
                    'src/tasks/task_imsubimage.py',
                    'src/tasks/task_accor.py',
                    'src/tasks/task_asdmsummary.py',
                    'src/tasks/task_clearcal.py',
                    'src/tasks/task_conjugatevis.py',
                    'src/tasks/task_exportasdm.py',
                    'src/tasks/task_importasdm.py',
                    'src/scripts/convertephem.py',
                    'src/tasks/task_clearstat.py',
                    'src/tasks/task_fixplanets.py',
                    'src/scripts/JPLephem_reader2.py',
                    'src/tasks/task_fixvis.py',
                    'src/tasks/task_fluxscale.py',
                    'src/tasks/task_ft.py',
                    'src/tasks/task_gaincal.py',
                    'src/tasks/task_gencal.py',
                    'src/scripts/correct_ant_posns.py',
                    'src/scripts/correct_ant_posns_alma.py',
                    'src/scripts/correct_ant_posns_evla.py',
                    'src/tasks/task_hanningsmooth.py',
                    'src/tasks/task_imcollapse.py',
                    'src/tasks/task_imcontsub.py',
                    'src/tasks/task_imdev.py',
                    'src/tasks/task_imfit.py',
                    'src/tasks/task_impbcor.py',
                    'src/tasks/task_importasap.py',
                    'src/tasks/task_importatca.py',
                    'src/tasks/task_importfitsidi.py',
                    'src/tasks/task_importgmrt.py',
                    'src/tasks/task_importnro.py',
                    'src/tasks/task_importvla.py',
                    'src/tasks/task_impv.py',
                    'src/tasks/task_imrebin.py',
                    'src/tasks/task_imreframe.py',
                    'src/tasks/task_imregrid.py',
                    'src/tasks/task_imsmooth.py',
                    'src/tasks/task_imtrans.py',
                    'src/tasks/task_imval.py',
                    'src/tasks/task_initweights.py',
                    'src/tasks/task_listcal.py',
                    'src/tasks/task_listfits.py',
                    'src/tasks/task_listhistory.py',
                    'src/tasks/task_listsdm.py',
                    'src/tasks/task_listvis.py',
                    'src/tasks/task_makemask.py',
                    'src/scripts/imtools.py',
                    'src/tasks/task_polcal.py',
                    'src/tasks/task_polfromgain.py',
                    'src/tasks/task_predictcomp.py',
                    'src/tasks/task_rerefant.py',
                    'src/scripts/restfreqtool.py',
                    'src/tasks/task_rmfit.py',
                    'src/tasks/task_rmtables.py',
                    'src/scripts/sdutil.py',
                    'src/tasks/task_sdbaseline.py',
                    'src/tasks/task_sdcal.py',
                    'src/tasks/task_sdfit.py',
                    'src/tasks/task_sdfixscan.py',
                    'src/tasks/task_sdgaincal.py',
                    'src/tasks/task_sdimaging.py',
                    'src/scripts/sdbeamutil.py',
                    'src/scripts/cleanhelper.py',
                    'src/tasks/task_sdimaging.py',
                    'src/tasks/task_sdsmooth.py',
                    'src/tasks/task_tsdimaging.py',
                    'src/tasks/task_nrobeamaverage.py',
                    'src/tasks/task_sdtimeaverage.py',
                    'src/tasks/task_sdpolaverage.py',
                    'src/tasks/task_sdsidebandsplit.py',
                    'src/tasks/task_plotprofilemap.py',
                    'src/scripts/simutil.py',
                    'src/tasks/task_simalma.py',
                    'src/tasks/task_simobserve.py',
                    'src/tasks/task_simanalyze.py',
                    'src/tasks/task_feather.py',
                    'src/tasks/task_statwt.py',
                    'src/tasks/task_virtualconcat.py',
                    'src/tasks/task_uvcontsub.py',
                    'src/tasks/task_uvmodelfit.py',
                    'src/tasks/task_visstat.py',
                    'src/tasks/task_widebandpbcor.py',
                    'src/tasks/task_importmiriad.py',
                    'src/tasks/task_plotweather.py',
                    'src/tasks/task_plotants.py',
                    'src/tasks/task_fringefit.py',
                    'src/tasks/task_plotbandpass.py',
                    'src/scripts/predictcomp_helper.py',
                    'src/scripts/task_logging.py',
                    'src/tasks/task_sdintimaging.py',
                    'src/scripts/sdint_helper.py',
                    'src/scripts/tec_maps.py',
]

private_modules = [ 'src/modules/parallel', 'src/modules/imagerhelpers' ]

xml_xlate = { 'casa-source/gcwrap/tasks/imhead.xml': 'xml/imhead.xml',
              'casa-source/gcwrap/tasks/immoments.xml': 'xml/immoments.xml',
              'casa-source/gcwrap/tasks/imhistory.xml': 'xml/imhistory.xml',
              'casa-source/gcwrap/tasks/applycal.xml': 'xml/applycal.xml',
              'casa-source/gcwrap/tasks/bandpass.xml': 'xml/bandpass.xml',
              'casa-source/gcwrap/tasks/blcal.xml': 'xml/blcal.xml',
              'casa-source/gcwrap/tasks/calstat.xml': 'xml/calstat.xml',
              'casa-source/gcwrap/tasks/concat.xml': 'xml/concat.xml',
              'casa-source/gcwrap/tasks/split.xml': 'xml/split.xml',
              'casa-source/gcwrap/tasks/listobs.xml': 'xml/listobs.xml',
              'casa-source/gcwrap/tasks/flagdata.xml': 'xml/flagdata.xml',
              'casa-source/gcwrap/tasks/flagcmd.xml': 'xml/flagcmd.xml',
              'casa-source/gcwrap/tasks/setjy.xml': 'xml/setjy.xml',
              'casa-source/gcwrap/tasks/cvel.xml': 'xml/cvel.xml',
              'casa-source/gcwrap/tasks/cvel2.xml': 'xml/cvel2.xml',
              'casa-source/gcwrap/tasks/importuvfits.xml': 'xml/importuvfits.xml',
              'casa-source/gcwrap/tasks/importfits.xml': 'xml/importfits.xml',
              'casa-source/gcwrap/tasks/exportfits.xml': 'xml/exportfits.xml',
              'casa-source/gcwrap/tasks/exportuvfits.xml': 'xml/exportuvfits.xml',
              'casa-source/gcwrap/tasks/partition.xml': 'xml/partition.xml',
              'casa-source/gcwrap/tasks/listpartition.xml': 'xml/listpartition.xml',
              'casa-source/gcwrap/tasks/flagmanager.xml': 'xml/flagmanager.xml',
              'casa-source/gcwrap/tasks/mstransform.xml': 'xml/mstransform.xml',
              'casa-source/gcwrap/tasks/tclean.xml': 'xml/tclean.xml',
              'casa-source/gcwrap/tasks/immath.xml': 'xml/immath.xml',
              'casa-source/gcwrap/tasks/vishead.xml': 'xml/vishead.xml',
              'casa-source/gcwrap/tasks/uvsub.xml': 'xml/uvsub.xml',
              'casa-source/gcwrap/tasks/spxfit.xml': 'xml/spxfit.xml',
              'casa-source/gcwrap/tasks/splattotable.xml': 'xml/splattotable.xml',
              'casa-source/gcwrap/tasks/specsmooth.xml': 'xml/specsmooth.xml',
              'casa-source/gcwrap/tasks/specflux.xml': 'xml/specflux.xml',
              'casa-source/gcwrap/tasks/smoothcal.xml': 'xml/smoothcal.xml',
              'casa-source/gcwrap/tasks/specfit.xml': 'xml/specfit.xml',
              'casa-source/gcwrap/tasks/imstat.xml': 'xml/imstat.xml',
              'casa-source/gcwrap/tasks/slsearch.xml': 'xml/slsearch.xml',
              'casa-source/gcwrap/tasks/delmod.xml': 'xml/delmod.xml',
              'casa-source/gcwrap/tasks/imsubimage.xml': 'xml/imsubimage.xml',
              'casa-source/gcwrap/tasks/accor.xml': 'xml/accor.xml',
              'casa-source/gcwrap/tasks/asdmsummary.xml': 'xml/asdmsummary.xml',
              'casa-source/gcwrap/tasks/clearcal.xml': 'xml/clearcal.xml',
              'casa-source/gcwrap/tasks/conjugatevis.xml': 'xml/conjugatevis.xml',
              'casa-source/gcwrap/tasks/exportasdm.xml': 'xml/exportasdm.xml',
              'casa-source/gcwrap/tasks/importasdm.xml': 'xml/importasdm.xml',
              'casa-source/gcwrap/tasks/clearstat.xml': 'xml/clearstat.xml',
              'casa-source/gcwrap/tasks/fixplanets.xml': 'xml/fixplanets.xml',
              'casa-source/gcwrap/tasks/fixvis.xml': 'xml/fixvis.xml',
              'casa-source/gcwrap/tasks/fluxscale.xml': 'xml/fluxscale.xml',
              'casa-source/gcwrap/tasks/ft.xml': 'xml/ft.xml',
              'casa-source/gcwrap/tasks/gaincal.xml': 'xml/gaincal.xml',
              'casa-source/gcwrap/tasks/gencal.xml': 'xml/gencal.xml',
              'casa-source/gcwrap/tasks/hanningsmooth.xml': 'xml/hanningsmooth.xml',
              'casa-source/gcwrap/tasks/imcollapse.xml': 'xml/imcollapse.xml',
              'casa-source/gcwrap/tasks/imcontsub.xml': 'xml/imcontsub.xml',
              'casa-source/gcwrap/tasks/imdev.xml': 'xml/imdev.xml',
              'casa-source/gcwrap/tasks/imfit.xml': 'xml/imfit.xml',
              'casa-source/gcwrap/tasks/impbcor.xml': 'xml/impbcor.xml',
              'casa-source/gcwrap/tasks/importasap.xml': 'xml/importasap.xml',
              'casa-source/gcwrap/tasks/importatca.xml': 'xml/importatca.xml',
              'casa-source/gcwrap/tasks/importfitsidi.xml': 'xml/importfitsidi.xml',
              'casa-source/gcwrap/tasks/importgmrt.xml': 'xml/importgmrt.xml',
              'casa-source/gcwrap/tasks/importnro.xml': 'xml/importnro.xml',
              'casa-source/gcwrap/tasks/importvla.xml': 'xml/importvla.xml',
              'casa-source/gcwrap/tasks/impv.xml': 'xml/impv.xml',
              'casa-source/gcwrap/tasks/imrebin.xml': 'xml/imrebin.xml',
              'casa-source/gcwrap/tasks/imreframe.xml': 'xml/imreframe.xml',
              'casa-source/gcwrap/tasks/imregrid.xml': 'xml/imregrid.xml',
              'casa-source/gcwrap/tasks/imsmooth.xml': 'xml/imsmooth.xml',
              'casa-source/gcwrap/tasks/imtrans.xml': 'xml/imtrans.xml',
              'casa-source/gcwrap/tasks/imval.xml': 'xml/imval.xml',
              'casa-source/gcwrap/tasks/initweights.xml': 'xml/initweights.xml',
              'casa-source/gcwrap/tasks/listcal.xml': 'xml/listcal.xml',
              'casa-source/gcwrap/tasks/listfits.xml': 'xml/listfits.xml',
              'casa-source/gcwrap/tasks/listhistory.xml': 'xml/listhistory.xml',
              'casa-source/gcwrap/tasks/listsdm.xml': 'xml/listsdm.xml',
              'casa-source/gcwrap/tasks/listvis.xml': 'xml/listvis.xml',
              'casa-source/gcwrap/tasks/makemask.xml': 'xml/makemask.xml',
              'casa-source/gcwrap/tasks/polcal.xml': 'xml/polcal.xml',
              'casa-source/gcwrap/tasks/polfromgain.xml': 'xml/polfromgain.xml',
              'casa-source/gcwrap/tasks/predictcomp.xml': 'xml/predictcomp.xml',
              'casa-source/gcwrap/tasks/rerefant.xml': 'xml/rerefant.xml',
              'casa-source/gcwrap/tasks/rmfit.xml': 'xml/rmfit.xml',
              'casa-source/gcwrap/tasks/rmtables.xml': 'xml/rmtables.xml',
              'casa-source/gcwrap/tasks/sdbaseline.xml': 'xml/sdbaseline.xml',
              'casa-source/gcwrap/tasks/sdcal.xml': 'xml/sdcal.xml',
              'casa-source/gcwrap/tasks/sdfit.xml': 'xml/sdfit.xml',
              'casa-source/gcwrap/tasks/sdfixscan.xml': 'xml/sdfixscan.xml',
              'casa-source/gcwrap/tasks/sdgaincal.xml': 'xml/sdgaincal.xml',
              'casa-source/gcwrap/tasks/sdimaging.xml': 'xml/sdimaging.xml',
              'casa-source/gcwrap/tasks/sdsmooth.xml': 'xml/sdsmooth.xml',
              'casa-source/gcwrap/tasks/tsdimaging.xml': 'xml/tsdimaging.xml',
              'casa-source/gcwrap/tasks/nrobeamaverage.xml': 'xml/nrobeamaverage.xml',
              'casa-source/gcwrap/tasks/sdtimeaverage.xml': 'xml/sdtimeaverage.xml',
              'casa-source/gcwrap/tasks/simalma.xml': 'xml/simalma.xml',
              'casa-source/gcwrap/tasks/simobserve.xml': 'xml/simobserve.xml',
              'casa-source/gcwrap/tasks/simanalyze.xml': 'xml/simanalyze.xml',
              'casa-source/gcwrap/tasks/feather.xml': 'xml/feather.xml',
              'casa-source/gcwrap/tasks/statwt.xml': 'xml/statwt.xml',
              'casa-source/gcwrap/tasks/virtualconcat.xml': 'xml/virtualconcat.xml',
              'casa-source/gcwrap/tasks/uvcontsub.xml': 'xml/uvcontsub.xml',
              'casa-source/gcwrap/tasks/uvmodelfit.xml': 'xml/uvmodelfit.xml',
              'casa-source/gcwrap/tasks/visstat.xml': 'xml/visstat.xml',
              'casa-source/gcwrap/tasks/widebandpbcor.xml': 'xml/widebandpbcor.xml',
              'casa-source/gcwrap/tasks/importmiriad.xml': 'xml/importmiriad.xml',
              'casa-source/gcwrap/tasks/plotweather.xml': 'xml/plotweather.xml',
              'casa-source/gcwrap/tasks/plotants.xml': 'xml/plotants.xml',
              'casa-source/gcwrap/tasks/fringefit.xml': 'xml/fringefit.xml',
              'casa-source/gcwrap/tasks/plotbandpass.xml': 'xml/plotbandpass.xml',
              'casa-source/gcwrap/tasks/sdintimaging.xml': 'xml/sdintimaging.xml',
              'casa-source/gcwrap/tasks/sdpolaverage.xml': 'xml/sdpolaverage.xml',
              'casa-source/gcwrap/tasks/sdsidebandsplit.xml': 'xml/sdsidebandsplit.xml',
              'casa-source/gcwrap/tasks/plotprofilemap.xml': 'xml/plotprofilemap.xml',
}

xml_files = [ 'xml/imhead.xml',
              'xml/immoments.xml',
              'xml/imhistory.xml',
              'xml/applycal.xml',
              'xml/bandpass.xml',
              'xml/blcal.xml',
              'xml/calstat.xml',
              'xml/concat.xml',
              'xml/split.xml',
              'xml/listobs.xml',
              'xml/flagdata.xml',
              'xml/flagcmd.xml',
              'xml/setjy.xml',
              'xml/cvel.xml',
              'xml/cvel2.xml',
              'xml/importuvfits.xml',
              'xml/importfits.xml',
              'xml/exportfits.xml',
              'xml/exportuvfits.xml',
              'xml/partition.xml',
              'xml/listpartition.xml',
              'xml/flagmanager.xml',
              'xml/mstransform.xml',
              'xml/tclean.xml',
              'xml/immath.xml',
              'xml/vishead.xml',
              'xml/uvsub.xml',
              'xml/spxfit.xml',
              'xml/splattotable.xml',
              'xml/specsmooth.xml',
              'xml/specflux.xml',
              'xml/smoothcal.xml',
              'xml/specfit.xml',
              'xml/imstat.xml',
              'xml/slsearch.xml',
              'xml/delmod.xml',
              'xml/imsubimage.xml',
              'xml/accor.xml',
              'xml/asdmsummary.xml',
              'xml/clearcal.xml',
              'xml/conjugatevis.xml',
              'xml/exportasdm.xml',
              'xml/importasdm.xml',
              'xml/clearstat.xml',
              'xml/fixplanets.xml',
              'xml/fixvis.xml',
              'xml/fluxscale.xml',
              'xml/ft.xml',
              'xml/gaincal.xml',
              'xml/gencal.xml',
              'xml/hanningsmooth.xml',
              'xml/imcollapse.xml',
              'xml/imcontsub.xml',
              'xml/imdev.xml',
              'xml/imfit.xml',
              'xml/impbcor.xml',
              'xml/importasap.xml',
              'xml/importatca.xml',
              'xml/importfitsidi.xml',
              'xml/importgmrt.xml',
              'xml/importnro.xml',
              'xml/importvla.xml',
              'xml/impv.xml',
              'xml/imrebin.xml',
              'xml/imreframe.xml',
              'xml/imregrid.xml',
              'xml/imsmooth.xml',
              'xml/imtrans.xml',
              'xml/imval.xml',
              'xml/initweights.xml',
              'xml/listcal.xml',
              'xml/listfits.xml',
              'xml/listhistory.xml',
              'xml/listsdm.xml',
              'xml/listvis.xml',
              'xml/makemask.xml',
              'xml/polcal.xml',
              'xml/polfromgain.xml',
              'xml/predictcomp.xml',
              'xml/rerefant.xml',
              'xml/rmfit.xml',
              'xml/rmtables.xml',
              'xml/sdbaseline.xml',
              'xml/sdcal.xml',
              'xml/sdfit.xml',
              'xml/sdfixscan.xml',
              'xml/sdgaincal.xml',
              'xml/sdimaging.xml',
              'xml/sdsmooth.xml',
              'xml/tsdimaging.xml',
              'xml/nrobeamaverage.xml',
              'xml/sdtimeaverage.xml',
              'xml/simalma.xml',
              'xml/simobserve.xml',
              'xml/simanalyze.xml',
              'xml/feather.xml',
              'xml/statwt.xml',
              'xml/virtualconcat.xml',
              'xml/uvcontsub.xml',
              'xml/uvmodelfit.xml',
              'xml/visstat.xml',
              'xml/widebandpbcor.xml',
              'xml/importmiriad.xml',
              'xml/plotweather.xml',
              'xml/plotants.xml',
              'xml/fringefit.xml',
              'xml/plotbandpass.xml',
              'xml/sdintimaging.xml',
              'xml/sdpolaverage.xml',
              'xml/sdsidebandsplit.xml',
              'xml/plotprofilemap.xml',
]

if pyversion < 3:
    str_encode = str
    str_decode = str
    def pipe_decode(output):
        return output
else:
    def str_encode(s):
        return bytes(s,sys.getdefaultencoding())
    def str_decode(bs):
        return bs.decode(sys.getdefaultencoding(),"strict")
    def pipe_decode(output):
        if isinstance(output,bytes) or isinstance(output,bytearray):
            return str_decode(output)
        elif isinstance(output,tuple):
            return ( None if output[0] is None else str_decode(output[0]), None if output[1] is None else str_decode(output[1]) )
        else:
            return ("","")

if sys.platform == 'darwin':
    def islib(l):
        return l.endswith('.dylib')
elif sys.platform == 'linux2' or sys.platform == 'linux':
    def islib(l):
        # 'lib_....so
        return l.find('.so') > 3
else:
    sys.exit("oops, had not planned on building on %s" % sys.platform)

def distutils_dir_name(dname):
    """Returns the name of a distutils build directory"""
    f = "{dirname}.{platform}-{version[0]}.{version[1]}"
    return f.format(dirname=dname,platform=get_platform(),version=sys.version_info)

def mkpath(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def upgrade_xml( conversions ):
    mkpath("xml")
    for k in conversions.keys( ):
        if not os.path.exists(conversions[k]) or os.path.getmtime(conversions[k]) < os.path.getmtime(k):

            print("upgrading %s" % k)

            proc = Popen( [tools_config['build.compiler.xml-casa'], "-upgrade", k],
                          stdout=subprocess.PIPE )

            (output, error) = pipe_decode(proc.communicate( ))

            exit_code = proc.wait( )
            if exit_code != 0:
                sys.exit('upgrading %s failed' % conversions[k])

            if k.endswith("/imval.xml"):
                print("fixing %s" % k)
                ## one test case test_imval.py is fundamentally broken, and fixing this through
                ## a commit causes the test case to proceed past the broken spot... so here
                ## we revert to patching at build time...
                output = output.replace('<output>','<bogus_output>')
                output = output.replace('</output>','</bogus_output>')

            xmlfd = open(conversions[k], 'w')
            xmlfd.write(output)
            xmlfd.close( )

def generate_pyinit(moduledir,tasks):
    """Generate __init__.py for the module
    """
    outfile = os.path.join(moduledir,'__init__.py')

    task_files_dict = list(map(lambda p: "'%s': %s" % p, zip(tasks,map(lambda x: "__os.path.join(__cwd,'%s')" % x.replace('xml','__xml__',1),xml_files))))

    with open(outfile, "w") as fd:
        fd.write("""###########################################################################\n""")
        fd.write("""########################## generated by setup.py ##########################\n""")
        fd.write("""###########################################################################\n""")
        fd.write("from __future__ import absolute_import\n")
        fd.write("from casatools import logsink as _logsink\n")
        fd.write("import os as __os\n")
        fd.write("__cwd = __os.path.dirname(__os.path.abspath(__file__))\n")
        fd.write("__name__ = '%s'\n" % module_name)
        fd.write("__all__ = [ \"casalog\", \"version\", \"version_string\",\n")
        for task in tasks:
            fd.write("            '%s',\n" % task)
        fd.write("          ]\n\n")
        fd.write("""from . import config\n""")
        fd.write("""casalog = _logsink( config.logfile )\n\n""")
        for task in tasks:
            fd.write("from .%s import %s\n" % (task,task))

        fd.write("\n")
        fd.write("def version( ): return [ %d, %d, %d, %d ]\n" % (casatasks_major,casatasks_minor,casatasks_patch,casatasks_feature))
        fd.write("def version_string( ): return \"%s\"\n" % casatasks_version)
        fd.write("casalog.setglobal(True)\n")
        fd.write("\n")
        fd.write("def xml_interface_defs( ): return { %s }\n" % ", ".join(task_files_dict))
        fd.write("\n")
        # Check if MPIEnvironment can be imported at all
        fd.write('mpi_env_found=False\n')
        fd.write('try:\n')
        fd.write('    from casampi.MPIEnvironment import MPIEnvironment\n')
        fd.write('    mpi_env_found=True\n')
        fd.write('except:\n')
        fd.write('    mpi_env_found=False\n')
        mpi_import_str = '\n'.join((
            "# When in MPI mode, this will put servers into their serve() loop.",
            "# From this point on user scripts can use tclean parallelization, Tier0 parallelization,",
            "# and MMS-parallel tasks",
            "try:",
            "    import importlib",
            "    _clith_spec = importlib.util.find_spec('casalith')",
            "    # Defer to later if in casalith",
            "    if _clith_spec is None:",
            "        import casampi.private.start_mpi",
            "except ImportError:",
            "    pass\n"))
        fd.write(mpi_import_str)
        fd.write("package_variant='modular'\n")
        fd.write("try:\n")
        fd.write("  if _clith_spec is not None:\n")
        fd.write("    package_variant='casalith'\n")
        fd.write("except:\n")
        fd.write("  pass\n")
        fd.write("try:\n")
        fd.write("  _pipe_spec = importlib.util.find_spec('pipeline')\n")
        fd.write("  if _pipe_spec is not None:\n")
        fd.write("    package_variant='pipeline'\n")
        fd.write("except:\n")
        fd.write("  pass\n")
        # Only the mpi "client" should write the version information (otherwise the logsink will crash)
        fd.write('if mpi_env_found and MPIEnvironment.is_mpi_enabled:\n')
        fd.write('    if MPIEnvironment.is_mpi_client:\n')
        fd.write('        try:\n')
        fd.write('            casalog.post("CASA Version " + package_variant.upper() + " %s")\n' %  casatasks_version)
        fd.write('            casalog.post("MPI Enabled")\n')
        fd.write('        except:\n')
        fd.write('            print("Error: the logfile is not writable")\n')
        fd.write('else:\n')
        fd.write('    try:\n')
        fd.write('        casalog.post("CASA Version " + package_variant.upper() + " %s")\n' % casatasks_version)
        fd.write('    except:\n')
        fd.write('        print("Error: the logfile is not writable")\n')  
        fd.write("\n")
        fd.write("from datetime import datetime as _time\n")
        fd.write("telemetry_starttime = str(_time.now())\n")
        fd.write("import platform\n")
        fd.write("import os\n")
        fd.write("serial_run = mpi_env_found and not MPIEnvironment.is_mpi_enabled\n")
        fd.write("mpi_run_client = mpi_env_found and MPIEnvironment.is_mpi_enabled and MPIEnvironment.is_mpi_client\n")
        fd.write("nompi_or_serial_or_client = not mpi_env_found or serial_run or mpi_run_client\n")
        fd.write("telemetry_available=False\n")
        fd.write("if nompi_or_serial_or_client:\n")
        fd.write("  try:\n")
        fd.write("    import casatelemetry\n")
        fd.write("    telemetry_available=True\n")
        fd.write("  except:\n")
        fd.write('    casalog.post("Can\'t import casatelemetry module.")\n')
        fd.write("if telemetry_available and config.telemetry_enabled and nompi_or_serial_or_client:\n")
        fd.write("  try:\n") 
        fd.write("    telemetrylogdirectory = None\n") 
        fd.write("    if config.rcdir != None:\n")
        fd.write("      telemetrylogdirectory = config.rcdir\n")
        fd.write("    if config.telemetry_log_directory != None:\n")
        fd.write("      telemetrylogdirectory = casatasks.config.telemetry_log_directory\n")
        fd.write("    telemetrylogger = casatelemetry.casatelemetry.telemetry(telemetrylogdirectory)\n")
        fd.write("  except:\n")       
        fd.write("    telemetrylogger = casatelemetry.casatelemetry.telemetry()\n")
        fd.write("  def logstop():\n")
        # Telemetry may be stopped during runtime so check if it is still enabled
        fd.write('    if telemetrylogger.telemetry_enabled:\n')
        fd.write("      telemetry_stoptime = str(_time.now())\n")            
        fd.write('      telemetrylogger.logger.info(telemetry_stoptime + " :: " + str(os.getpid()) + " :: CASAStop :: Stopping CASA at: " + telemetry_stoptime + \n') 
        fd.write('      " Version " + version_string() + " Platform: " + platform.platform() +  " Start time: " + telemetry_starttime + " Variant: " + package_variant)\n')
        #fd.write('   print("mpi_env " + str(mpi_env_found))\n')
        fd.write('  if not mpi_env_found or (mpi_env_found and MPIEnvironment.is_mpi_client):\n')
        fd.write('    telemetrylogger.submitStatistics()\n')
        fd.write('    if telemetrylogger.telemetry_enabled:\n')
        fd.write('      telemetrylogger.logger.info(telemetry_starttime + " :: " + str(os.getpid()) + " :: CASAStart :: Starting CASA at: " + telemetry_starttime + \n') 
        fd.write('      " Version " + version_string() + " Platform: " + platform.platform() + " Variant: " + package_variant)\n')
        fd.write("    import atexit\n")
        fd.write("    atexit.register(logstop)\n")
        fd.write('if telemetry_available and config.crashreporter_enabled:\n')
        fd.write("    casatelemetry.CrashReporter.init(config.logfile)\n")
        

class BuildCasa(build):
    description = "Description of the command"
    user_options = []

    # This method must be implemented
    def initialize_options(self):
        print("initializing options...")
        build.initialize_options(self)

    # This method must be implemented
    def finalize_options(self):
        #self.plat_name = get_platform( )
        print("finalizing options...")
        build.finalize_options(self)

    def run(self):
        upgrade_xml(xml_xlate)
        #scripts/xml-casa output-task=/tmp/debug/task -task taskxml/imhead.xml > /Users/drs/develop/casa/casatools/build/lib.macosx-10.12-x86_64-3.6/casatasks/imhead.py
        libdir = os.path.join("build",distutils_dir_name('lib'))
        moduledir = os.path.join(libdir,"casatasks")
        privatedir = os.path.join(moduledir,"private")
        xmldir = os.path.join(moduledir,"__xml__")

        mkpath(privatedir)
        print("generating task python files...")
        proc = Popen( [tools_config['build.compiler.xml-casa'], "output-task=%s" % moduledir, "-task"] + xml_files,
                      stdout=subprocess.PIPE )

        (output, error) = pipe_decode(proc.communicate( ))

        exit_code = proc.wait( )

        if exit_code != 0:
            sys.exit('python file generation failed')

        tasks = output.split( )
        generate_pyinit(moduledir,tasks)

        for f in public_scripts:
            copy2(f,moduledir)

        for f in private_scripts:
            copy2(f,privatedir)

        for m in private_modules:
            tgt = os.path.join(privatedir,os.path.basename(m))
            copy_tree(m,tgt)

        mkpath(xmldir)
        for x in xml_files:
            copy2(x,xmldir)

class TestCasa(Command):
    user_options = []

    def initialize_options(self):
        self.__test_dir = "build/%s" % distutils_dir_name('testing')
        self.__lib_dir = os.path.abspath("build/%s" % distutils_dir_name('lib'))
        self.__env = os.environ.copy( )
        if 'PYTHONPATH' in self.__env:
            existing_paths = ":".join(list(map(os.path.abspath,self.__env['PYTHONPATH'].split(':'))))
            self.__env['PYTHONPATH'] = "%s:%s" % (self.__lib_dir,existing_paths)
        else:
            self.__env['PYTHONPATH'] = self.__lib_dir
        self.__regression_dir = "build/%s" % distutils_dir_name('regression')
        #self.__regression_ref_dir = "tests/output-regression/reference-%d.%d" % (sys.version_info[0],sys.version_info[1])
        self.__regression_ref_dir = "tests/output-regression/reference"
        self.__regression_sample_gen = "scripts/output-snapshot"
        self.__regression_sample_dir = "%s/%s" % (self.__regression_dir,"output")


    def finalize_options(self):
        pass

    def __dump_output(self, working_dir, bname, out, err):
        stdout_path = "%s/%s-stdout.txt" % (working_dir,bname)
        stdout_fd = open(stdout_path, 'w')
        stdout_fd.write(out)
        stdout_fd.close( )
        stderr_path = "%s/%s-stderr.txt" % (working_dir,bname)
        stderr_fd = open(stderr_path, 'w')
        stderr_fd.write(err)
        stderr_fd.close( )
        return (stdout_path,stderr_path)

    def __run_test(self,tabwidth,test_path,working_dir):
        label = '.'.join(os.path.basename(test_path).split('.')[:-1])
        sys.stdout.write(label + '.' * (tabwidth - len(label)))
        sys.stdout.flush( )
        proc = Popen( [sys.executable,test_path], cwd=working_dir, env=self.__env,
                      stdout=subprocess.PIPE, stderr=subprocess.PIPE )
        (output, error) = pipe_decode(proc.communicate( ))
        exit_code = proc.wait( )
        (stdout_path,stderr_path) = self.__dump_output(working_dir,"log",output,error)
        print(" ok" if exit_code == 0 else " fail")
        return (exit_code, label, stdout_path, stderr_path)

    def __generate_sample(self):
        if os.path.exists(self.__regression_dir):
            remove_tree(self.__regression_dir)
        mkpath(self.__regression_sample_dir)
        proc = Popen( [self.__regression_sample_gen,"out=%s" % self.__regression_sample_dir],
                      env=self.__env, stdout=subprocess.PIPE, stderr=subprocess.PIPE )
        (output, error) = pipe_decode(proc.communicate( ))
        exit_code = proc.wait( )
        self.__dump_output(self.__regression_dir,"sample-generation",output,error)
        return exit_code

    def __compare_reg(self, tabwidth, label, refpath, samplepath):
        sys.stdout.write(label + '.' * (tabwidth - len(label)))
        sys.stdout.flush( )
        proc = Popen( ["/usr/bin/diff",refpath,samplepath], env=self.__env,
                      stdout=subprocess.PIPE, stderr=subprocess.PIPE )
        (output, error) = pipe_decode(proc.communicate( ))
        exit_code = proc.wait( )
        print(" ok" if exit_code == 0 else " fail")
        (op,ep) = self.__dump_output(os.path.dirname(samplepath),label,output,error)
        return (exit_code, label, op, ep)

    def __collect_tests(self, testdir):
        tests = [ ]
        for dir, subdirs, files in os.walk(testdir):
            for f in files:
                if f.endswith(".py") and f.startswith("test_"):
                    workingdir = "%s/%s" % (self.__test_dir,f[:-3])
                    mkpath(workingdir)
                    tests.append((os.path.abspath("%s/%s" % (dir,f)),workingdir))
        return tests

    def __collect_regression_files(self, regdir):
        regression = { }
        for dir, subdirs, files in os.walk(regdir):
            for f in files:
                if f != "log.txt":
                    regression[f] = "%s/%s" % (dir,f)
        return regression

    def __collect_regression(self):
        regression_ref = { }
        regression_sample = { }
        if os.path.isdir(self.__regression_ref_dir):
            if isexe(self.__regression_sample_gen):
                if self.__generate_sample( ) == 0:
                    regression_ref = self.__collect_regression_files(self.__regression_ref_dir)
                    regression_sample = self.__collect_regression_files(self.__regression_sample_dir)
                else:
                    print("warning, generation of regression sample failed; skipping regression test")
            else:
                print( "warning, regression sample generator (%s) does not exist; skipping regression test" %
                       self.__regression_sample_gen)
        else:
            print("warning, regression reference (%s) does not exist; skipping regression test" % self.__regression_ref_dir)
        return (regression_ref, regression_sample)

    def run(self):
        if os.path.exists(self.__test_dir):
            remove_tree(self.__test_dir)
        mkpath(self.__test_dir)
        tests = self.__collect_tests("tests/tasks")

        (regression_ref, regression_sample) = self.__collect_regression( )

        testwidth = 0 if len(tests) == 0 else max(map(lambda x: len(os.path.basename(x[0]))+3,tests))
        regressionwidth = 0 if len(regression_ref) == 0 else max(map(lambda x: len(x)+3,regression_ref.keys( )))
        tabwidth = max(testwidth,regressionwidth,45)

        testresults = list(map(lambda params: self.__run_test(tabwidth,*params),tests))

        len_message = "regression file count"
        print( len_message + '.' * (tabwidth - len(len_message)) +
               (" ok" if len(regression_ref) == len(regression_sample) else " fail") )

        regression_keys = filter(lambda k: k in regression_sample, regression_ref.keys( ))
        regressionresults = list(map(lambda k: self.__compare_reg(tabwidth,k,regression_ref[k],regression_sample[k]), regression_keys))

        results = testresults + regressionresults
        print('-' * (tabwidth + 8))
        passed = list(filter(lambda v: v[0] == 0,results))
        failed = list(filter(lambda v: v[0] != 0,results))
        print("test summary %d passed, %d failed" % (len(passed),len(failed)))
        print("OK" if len(failed) == 0 else "FAIL")
        sys.exit(0 if len(failed) == 0 else 1)

cmd_setup = { 'build': BuildCasa, 'test': TestCasa }
try:
    from wheel.bdist_wheel import bdist_wheel
    class casa_binary_wheel(bdist_wheel):
        user_options = bdist_wheel.user_options + [
            # The format is (long option, short option, description).
            ('debug', None, 'build a debugging version'),
            ('version=', None, 'Silence distutils when passing version for bdist_wheel from the command line'),
        ]

        def initialize_options(self):
            self.version = None
            bdist_wheel.initialize_options(self)
            if '--debug' in sys.argv:
                global _debug_build_
                _debug_build_ = True
                self.debug = 1
            else:
                self.debug = 0

        def run(self):
            global doing_wheel_build
            doing_wheel_build = True
            bdist_wheel.run(self)
            doing_wheel_build = False

    wheel_build = casa_binary_wheel

    cmd_setup['bdist_wheel'] = casa_binary_wheel
except ImportError:
    pass  # custom command not needed if wheel is not installed

setup( name=module_name,version=casatasks_version,
       maintainer="Darrell Schiebel",
       maintainer_email="drs@nrao.edu",
       author="CASA development team",
       author_email="aips2-request@nrao.edu",
       url="https://open-bitbucket.nrao.edu/projects/CASA/repos/casatools/browse",
       download_url="https://casa.nrao.edu/download/",
       license="GNU Library or Lesser General Public License (LGPL)",
       packages=[ module_name,
                  "%s.__xml__" % module_name,
                  "%s.private" % module_name,
                  "%s.private.parallel" % module_name,
                  "%s.private.imagerhelpers" % module_name ],
       classifiers=[ 'Programming Language :: Python :: %s' % pyversion ],
       description="the CASA tasks",
       long_description="The CASAtasks are a collection of (mostly) stateless functions for\nthe analysis of radio astronomy observations.",
       cmdclass=cmd_setup,
       package_dir={module_name: os.path.join('build',distutils_dir_name('lib'), module_name)},
       package_data={'': ['*.xml','*.txt']},
       install_requires=[ 'casatools==%s' % casatools.version_string( ), 'matplotlib', 'scipy' ]
)
