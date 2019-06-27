"""
casa6tools = [
    "agentflagger", "atcafiller", "atmosphere", "calanalysis", "calibrater", "coercetype", "componentlist", "config", "constants", "coordsys", "ctuser", "functional", "image", 
    "imagemetadata", "imagepol", "imager", "iterbotsink", "logsink", "measures", "miriadfiller", "ms", "msmetadata", "mstransformer", "platform", "quanta", "regionmanager", "sakura", 
    "sdm", "simulator", "singledishms", "spectralline", "synthesisdeconvolver", "synthesisimager", "synthesisimstore", "synthesisnormalizer", "synthesisutils", "table", "typecheck", "utils", 
    "vlafiller", "vpmanager"
    ]
"""
casa6tasks = set([
    "accor", "accum", "applycal", "asdmsummary", "bandpass", "blcal", "calstat","clearcal", "clearstat", "concat", "config", "conjugatevis", "cvel", "cvel2",
    "delmod", "exportasdm", "exportfits", "exportuvfits", "feather", "fixplanets","fixvis", "flagcmd", "flagdata", "flagmanager", "fluxscale", "ft", "gaincal",
    "gencal", "hanningsmooth", "imcollapse", "imcontsub", "imdev", "imfit", "imhead","imhistory", "immath", "immoments", "impbcor", "importasap", "importasdm",
    "importatca", "importfits", "importfitsidi", "importgmrt", "importmiriad","importnro", "importuvfits", "importvla", "impv", "imrebin", "imreframe",
    "imregrid", "imsmooth", "imstat", "imsubimage", "imtrans", "imval","initweights", "listcal", "listfits", "listhistory", "listobs", "listpartition",
    "listsdm", "listvis", "makemask", "mstransform", "partition", "polcal","predictcomp", "rerefant", "rmfit", "rmtables", "sdbaseline", "sdcal",
    "sdfit", "sdfixscan", "sdgaincal", "sdimaging", "sdsmooth", "setjy","simalma", "simanalyze", "simobserve", "slsearch", "smoothcal", "specfit",
    "specflux", "specsmooth", "splattotable", "split", "spxfit", "statwt","tclean", "uvcontsub", "uvmodelfit", "uvsub", "virtualconcat", "vishead", "visstat","widebandpbcor" ])

miscellaneous_tasks = set(['wvrgcal','plotms'])

import os, sys, time
from functools import wraps
import psutil
import fnmatch
import logging
import filecmp
import unittest
import tracemalloc
import pickle
import numpy
import math
import numbers



logging.basicConfig(level=logging.INFO,format='%(message)s')
#logging.basicConfig(level=logging.DEBUG,format='%(levelname)s-%(message)s')

"""

logging.debug('This is a debug message')
logging.info('This is an info message')
logging.warning('This is a warning message')
logging.error('This is an error message')
logging.critical('This is a critical message')
"""
casa5 = False
casa6 = False

try:
    # CASA 6
    logging.debug("Importing CASAtools")
    import casatools
    logging.debug("Importing CASAtasks")
    import casatasks

    tb  = casatools.table()
    tb2 = casatools.table()
    ms  = casatools.ms()
    ia  = casatools.image()

    from casatasks import casalog
    casa6 = True
except ImportError:
    # CASA 5
    from __main__ import default
    from taskinit import tbtool, mstool, iatool
    from taskinit import *
    from casa_stack_manip import stack_find, find_casa

    tb  = tbtool()
    tb2 = tbtool()
    ms  = mstool()
    ia  = iatool()

    casa = find_casa( )
    if casa.has_key('state') and casa['state'].has_key('init_version') and casa['state']['init_version'] > 0:
        casaglobals=True
        casac = stack_find("casac")
        casalog = stack_find("casalog")

    casa5 = True

############################################################################################
##################################       Classes        ##################################
############################################################################################

# Logger
class Logger:
    #TODO: This class needs work
    import sys
    import logging

    def verbose_logging_start():
        logger = logging.getLogger()
        logger.level = logging.DEBUG
        stream_handler = logging.StreamHandler(sys.stdout)
        logger.addHandler(stream_handler)

    def verbose_logging_stop():
        pass


# Weblog 
class Weblog:
    def __init__(self, taskname, localdict):
        self.localdict = localdict
        self.taskname = taskname

    def write_modal_style(self):


        html.write(' /* Style the Image Used to Trigger the Modal */' + '\n')
        html.write('.myImg {' + '\n') 
        html.write('border-radius: 5px; cursor: pointer; transition: 0.3s; }'+ '\n')

        html.write('.myImg:hover{' + '\n')
        html.write('opacity: 0.7;}'+ '\n')

        html.write('/* The Modal (background) */' + '\n')
        html.write('.modal {' + '\n')
        html.write('  display: none; /* Hidden by default */' + '\n')
        html.write('  position: fixed; /* Stay in place */' + '\n')
        html.write('  z-index: 1; /* Sit on top */' + '\n')
        html.write('  padding-top: 100px; /* Location of the box */' + '\n')
        html.write('  left: 0;' + '\n')
        html.write('  top: 0;' + '\n')
        html.write('  width: 100%; /* Full width */' + '\n')
        html.write('  height: 100%; /* Full height */' + '\n')
        html.write('  overflow: auto; /* Enable scroll if needed */' + '\n')
        html.write('  background-color: rgb(0,0,0); /* Fallback color */' + '\n')
        html.write('  background-color: rgba(0,0,0,0.9); /* Black w/ opacity */' + '\n')
        html.write('}' + '\n')
        html.write('/* Modal Content (Image) */' + '\n')
        html.write('.modal-content {' + '\n')
        html.write('  margin: auto;' + '\n')
        html.write('  display: block;' + '\n')
        html.write('  width: 80%;' + '\n')
        html.write('  max-width: 700px;' + '\n')
        html.write('}' + '\n')

        html.write('/* Caption of Modal Image (Image Text) - Same Width as the Image */' + '\n')
        html.write('#caption {' + '\n')
        html.write('  margin: auto;' + '\n')
        html.write('  display: block;' + '\n')
        html.write('  width: 80%;' + '\n')
        html.write('  max-width: 700px;' + '\n')
        html.write('  text-align: center;' + '\n')
        html.write('  color: #ccc;' + '\n')
        html.write('  padding: 10px 0;' + '\n')
        html.write('  height: 150px;' + '\n')
        html.write('}' + '\n')

        html.write('/* Add Animation - Zoom in the Modal */' + '\n')
        html.write('.modal-content, #caption {' + '\n')
        html.write('  animation-name: zoom;' + '\n')
        html.write('  animation-duration: 0.6s;' + '\n')
        html.write('}' + '\n')

        html.write('@keyframes zoom {' + '\n')
        html.write('  from {transform:scale(0)}' + '\n')
        html.write('  to {transform:scale(1)}' + '\n')
        html.write('}' + '\n')

        html.write('/* The Close Button */' + '\n')
        html.write('.close {' + '\n')
        html.write('  position: absolute;' + '\n')
        html.write('  top: 15px;' + '\n')
        html.write('  right: 35px;' + '\n')
        html.write('  color: #f1f1f1;' + '\n')
        html.write('  font-size: 40px;' + '\n')
        html.write('  font-weight: bold;' + '\n')
        html.write('  transition: 0.3s;' + '\n')
        html.write('}' + '\n')

        html.write('.close:hover,' + '\n')
        html.write('.close:focus {' + '\n')
        html.write('  color: #bbb;' + '\n')
        html.write('  text-decoration: none;' + '\n')
        html.write('  cursor: pointer;' + '\n')
        html.write('}' + '\n')

        html.write('/* 100% Image Width on Smaller Screens */' + '\n')
        html.write('@media only screen and (max-width: 700px){' + '\n')
        html.write('  .modal-content {' + '\n')
        html.write('    width: 100%;' + '\n')
        html.write('  }' + '\n')
        html.write('} ' + '\n')

    def generate_header(self, testname):
        html.write('<!doctype html>' + '\n')
        html.write('<html lang="en">' + '\n')
        html.write('<head>' + '\n')
        html.write('<meta charset="utf-8">' + '\n')
        html.write('<title>{}</title>'.format(testname) + '\n')
        html.write('<meta name="description" content="The HTML5 Herald">' + '\n')
        html.write('<meta name="author" content="SitePoint">' + '\n')
        html.write('<link rel="stylesheet" href="css/styles.css?v=1.0">' + '\n')
        html.write('</head>' + '\n')
        html.write('<body>' + '\n')
        html.write('<script src="js/scripts.js"></script>' + '\n')
        html.write('<h1>{}</h1>'.format(testname) + '\n')


    def generate_status_table_style(self,dictionary):
        html.write('<style type="text/css">' + '\n')
        html.write('.collapsible {background-color: #777;color: white;cursor: pointer;padding: 18px;width: 100%;border: none;text-align: left;outline: none;font-size: 15px;}' + '\n')
        html.write('.active, .collapsible:hover { background-color: #555;}' + '\n')
        html.write('.content {padding: 0 18px;display: none;overflow: hidden;background-color: #f1f1f1;}' + '\n')
        html.write(".boxed { border: 1px solid green ;padding: 0 5px 0 5px;margin: 50px}" + '\n')
        html.write('.tg  {border-collapse:collapse;border-spacing:0;}' + '\n')
        html.write('.tg td{font-family:Arial, sans-serif;font-size:14px;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:black;}' + '\n')
        html.write('.tg th{font-family:Arial, sans-serif;font-size:14px;font-weight:normal;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:black;}' + '\n')
        html.write('.tg .tg-0lax{text-align:left;vertical-align:top}' + '\n')
        html.write('.tg .tg-ck9b{background-color:#009901;color:#32cb00;text-align:left;vertical-align:top}' + '\n')
        html.write('.tg .tg-r50r{background-color:#cb0000;text-align:left;vertical-align:top}' + '\n')
        html.write('.tg-sort-header::-moz-selection{background:0 0}.tg-sort-header::selection{background:0 0}.tg-sort-header{cursor:pointer}.tg-sort-header:after{content:'';float:right;margin-top:7px;border-width:0 5px 5px;border-style:solid;border-color:#404040 transparent;visibility:hidden}.tg-sort-header:hover:after{visibility:visible}.tg-sort-asc:after,.tg-sort-asc:hover:after,.tg-sort-desc:after{visibility:visible;opacity:.4}.tg-sort-desc:after{border-bottom:none;border-width:5px 5px 0}' + '\n')

        Weblog(self.taskname, self.localdict).write_modal_style()

        html.write('</style>' + '\n')

    def generate_tail(self,dictionary):
        html.write('<script>' + '\n')
        html.write('var coll = document.getElementsByClassName("collapsible");' + '\n')
        html.write('var i;' + '\n')
        html.write('for (i = 0; i < coll.length; i++) {' + '\n')
        html.write('  coll[i].addEventListener("click", function() {' + '\n')
        html.write('    this.classList.toggle("active");' + '\n')
        html.write('    var content = this.nextElementSibling;' + '\n')
        html.write('    if (content.style.display === "block") {' + '\n')
        html.write('      content.style.display = "none";' + '\n')
        html.write('    } else {' + '\n')
        html.write('      content.style.display = "block";' + '\n')
        html.write('    }' + '\n')
        html.write('  });' + '\n')
        html.write('}' + '\n')

        html.write('// Get the modal' + '\n')
        html.write("var modal = document.getElementById('myModal');" + '\n')
        html.write('// Get the image and insert it inside the modal - use its "alt" text as a caption' + '\n')
        html.write("var img = $('.myImg');" + '\n')
        html.write('var modalImg = $("#img01");' + '\n')
        html.write('var captionText = document.getElementById("caption");' + '\n')
        html.write("$('.myImg').click(function(){" + '\n')
        html.write('    modal.style.display = "block";' + '\n')
        html.write('    var newSrc = this.src;' + '\n')
        html.write("    modalImg.attr('src', newSrc);" + '\n')
        html.write('    captionText.innerHTML = this.alt;' + '\n')
        html.write('});' + '\n')
        html.write('// Get the <span> element that closes the modal' + '\n')
        html.write('var span = document.getElementsByClassName("close")[0];' + '\n')
        html.write('// When the user clicks on <span> (x), close the modal' + '\n')
        html.write('span.onclick = function() {' + '\n')
        html.write('  modal.style.display = "none";' + '\n')
        html.write('}' + '\n')

        html.write('</script>' + '\n')
        html.write("</body>" + '\n')
        html.write("</html>" + '\n')

    def generate_table_row(self, test, description, runtime, status_color):
        
        html.write("<tr>" + "\n")
        html.write('<td class="tg-0lax">{}</td>'.format(test) + '\n')
        html.write('<td class="tg-0lax">{}</td>'.format(description) + '\n')
        html.write('<td class="tg-0lax">{}s</td>'.format(round(runtime,2)) + '\n')
        html.write('<td class={}></td>'.format(status_color) + '\n')
        html.write("</tr>" + "\n")

    def generate_status_table(self, dictionary):
        html.write('<table id="tg-cC48w" class="tg">' + '\n')
        html.write('<tr>' + '\n')
        html.write('<th class="tg-0lax">Test Name</th>' + '\n')
        html.write('<th class="tg-0lax">Description </th>' + '\n')
        html.write('<th class="tg-0lax">Run Time</th>' + '\n')
        html.write('<th class="tg-0lax">Status</th>' + '\n')
        html.write('</tr>' + '\n')
        for key, value in dictionary.items():
            Weblog(self.taskname, self.localdict).generate_table_row(str(key), dictionary[key]['description'],  dictionary[key]['runtime'], "tg-ck9b" if dictionary[key]['status'] == True else "tg-r50r" )
        html.write('</table>' + '\n')

    def generate_summary_box(self, dictionary):
        for key, value in dictionary.items():
            html.write('<button class="collapsible">{}</button>'.format(key) + '\n')
            html.write('<div class="content">'+ '\n')
            html.write('<div class="boxed">'+ '\n')
            html.write('<h3>{}</h3>'.format(key) + '\n')
            html.write('<i><sub>{}</sub></i>'.format(dictionary[key]['description'])+ '\n')
            html.write('<p><b>Elapsed Time:</b> {} Seconds</p>'.format(dictionary[key]['runtime'])+ '\n')
            html.write('<p><b>Status:</b> {}</p>'.format("PASS" if dictionary[key]['status'] == True else "FAIL" )+ '\n')
            html.write('<p><b>Task Executions:</b></p>'+ '\n')
            Weblog(self.taskname, self.localdict).write_inline_list( dictionary[key]['taskcall'] )
            Weblog(self.taskname, self.localdict).add_miscellaneous_info(dictionary[key])
            html.write('<p><b>Re-Run:</b> {}</p>'.format(dictionary[key]['rerun'])+ '\n')
            html.write('</div>' + '\n')
            html.write('</div>' + '\n') 

    def write_inline_list(self, array):
        html.write('<ul>' + '\n')
        for item in array:
            html.write('<li>{}</li>'.format(item) + '\n')
            if str(item).endswith(".png"):
                html.write('<script src="https://ajax.googleapis.com/ajax/libs/jquery/2.1.1/jquery.min.js"></script>'+ '\n')
                html.write('<img  class="myImg" src="{}" alt="{}" height="300" width="300">'.format(item, item) + '\n')
                html.write('<div id="myModal" class="modal">'+ '\n')
                html.write('  <span class="close" onclick="document.getElementById(\'myModal\').style.display=\'none\'">&times;</span>' + '\n')
                html.write('  <img class="modal-content" id="img01">' + '\n')
                html.write('  <div id="caption"></div>' + '\n')
                html.write('</div>' + '\n')
                html.write('</img>' + '\n')
        html.write('</ul>' + '\n')

    def write_inline_dict(self, dictionary):
        html.write('<ul>' + '\n')
        for key, value in sorted(dictionary.items()):
            html.write('<p><b>{}</b>: {}</p>'.format(key,value)+ '\n')
        html.write('</ul>' + '\n')

    def add_miscellaneous_info(self, subdictionary):
        default_keys = ['description','status','runtime','taskcall','rerun']
        for key, value in subdictionary.items():
            if key in default_keys:
                continue
            if type(subdictionary[key]) == list:
                html.write('<p><b>{}:</b></p>'.format(key)+ '\n')
                Weblog(self.taskname, self.localdict).write_inline_list(subdictionary[key])
            elif type(subdictionary[key]) == dict:
                html.write('<p><b>{}:</b></p>'.format(key)+ '\n')
                Weblog(self.taskname, self.localdict).write_inline_dict(subdictionary[key])
            else:
                html.write('<p><b>{}:</b> {}</p>'.format(key,subdictionary[key])+ '\n')

    def generate_weblog(self):
        #print("Generating Weblog: {}".format(self.taskname))
        Weblog(self.taskname, self.localdict).generate_header("Test {}".format(self.taskname))
        Weblog(self.taskname, self.localdict).generate_status_table_style(self.localdict)
        Weblog(self.taskname, self.localdict).generate_status_table(self.localdict)
        Weblog(self.taskname, self.localdict).generate_summary_box(self.localdict)
        Weblog(self.taskname, self.localdict).generate_tail(self.localdict)

############################################################################################
##################################       Functions        ##################################
############################################################################################
def compare_CASA_variable_cols(referencetab, testtab, varcol, tolerance=0.0):
    '''
    compare_CASA_variable_cols - Compare a variable column of two CASA tables.
       @param referencetab  --> a reference table
       @param testtab       --> a table to verify
       @param varcol        --> the name of a variable column (str)
       @param tolerance     --> Tolerance 

       @return: True if reference tab == test table else False
    '''
    logging.info("Comparing {} to {}".format(referencetab, testtab))
    logging.debug("Executing: compare_CASA_variable_cols(referencetab={},testtab={}, varcol={}, tolerance={})".format(referencetab, testtab, varcol, tolerance))
    retval = True

    tb.open(referencetab)
    cnames = tb.colnames()

    tb2.open(testtab)
    col = varcol
    if tb.isvarcol(col) and tb2.isvarcol(col):
        try:
            # First check
            if tb.nrows() != tb2.nrows():
                logging.error('Length of {} differ from {}, {} != {}'.format(referencetab,testtab,tb.nrows(),tb.nrows()))
                retval = False
            else:
                for therow in range(tb.nrows()):
                    rdata = tb.getcell(col,therow)
                    tdata = tb2.getcell(col,therow)

                    if not rdata.all()==tdata.all():
                        if (tolerance>0.0):
                            differs=False
                            for j in range(0,len(rdata)):

                                if ((isinstance(rdata[j],float)) or (isinstance(rdata[j],int))):
                                    if (abs(rdata[j]-tdata[j]) > tolerance*abs(rdata[j]+tdata[j])):
#                                        print('Column ', col,' differs in tables ', referencetab, ' and ', testtab)
#                                        print(therow, j)
#                                        print(rdata[j])
#                                        print(tdata[j])
                                        differs = True
                                elif (isinstance(rdata[j],list)) or (isinstance(rdata[j],numpy.ndarray)):
                                    for k in range(0,len(rdata[j])):
                                        if (abs(rdata[j][k]-tdata[j][k]) > tolerance*abs(rdata[j][k]+tdata[j][k])):
#                                            print('Column ', col,' differs in tables ', referencetab, ' and ', testtab)
#                                            print(therow, j, k)
#                                            print(rdata[j][k])
#                                            print(tdata[j][k])
                                            differs = True
                                if differs:
                                    print('ERROR: Column %s of %s and %s do not agree within tolerance %s'%(col,referencetab, testtab, tolerance))
                                    retval = False
                                    break
                        else:
                            print('ERROR: Column %s of %s and %s do not agree.'%(col,referencetab, testtab))
                            print('ERROR: First row to differ is row=%s'%therow)
                            retval = False
                            break
        finally:
            tb.close()
            tb2.close()
    
    else:
        logging.info('Columns are not varcolumns.')
        retval = False

    if retval:
        logging.info('Column {} of {} and {} agree'.format(col,referencetab, testtab))
        
    return retval


def compare_CASA_tables(referencetab, testtab, excludecols = [], tolerance=0.001, mode="percentage", startrow = 0, nrow = -1, rowincr = 1):
    '''
    compare_CASA_tables - compare two CASA tables
       @param referencetab - the table which is assumed to be correct
       @param testtab - the table which is to be compared to referencetab
       @param excludecols - list of column names which are to be ignored
       @param tolerance - permitted fractional difference (default 0.001 = 0.1 percent)
       @param mode - comparison is made as "percentage", "absolute", "phaseabsdeg" (for complex numbers = difference of the phases in degrees)

       @return: True if reference tab == test table else False
    '''
    logging.info("Comparing {} to {}".format(referencetab, testtab))
    logging.debug("Executing: compare_CASA_tables(referencetab = {}, testtab = {}, excludecols = {}, tolerance={}, mode={}, startrow = {}, nrow = {}, rowincr = {})".format(referencetab, testtab, excludecols, tolerance, mode, startrow, nrow , rowincr))

    if not isinstance(excludecols, list):
        logging.error("excludecols not in correct format")
        raise TypeError("excludecols must be a list")

    if referencetab.endswith("cal") or testtab.endswith("cal"):
        logging.warning("WARNING: Will compare caltables using compare_caltables")
        return compare_caltables(referencetab, testtab, cols= excludecols, rtol=8e-7, atol=1e-8)

    ##### Begin: Tempory Fix
    if len(excludecols) == 0:
        excludecols = ["FLAG_CATEGORY"]
    else:
        excludecols.append('FLAG_CATEGORY')
    """
    #TODO: Fix Error in checking FLAG_CATEGORY
        tb.getcol("FLAG_CATEGORY")
            SEVERE  getcol::FLAG_CATEGORY   Exception Reported: Table DataManager error: Invalid operation: TSM: no array in row 0 of column FLAG_CATEGORY in **
            RuntimeError: Table DataManager error: Invalid operation: TSM: no array in row 0 of column FLAG_CATEGORY in **
    """
    ##### End: Tempory Fix

    
    rval = True

    # Open reference table
    tb.open(referencetab)
    cnames = tb.colnames()

    # Open test table
    tb2.open(testtab)
    cnames2 = tb2.colnames()

    
    if sorted(cnames) != sorted(cnames2):
        logging.debug("Available columns in Reference Table {}: {}".format(referencetab,cnames))
        logging.debug("Available columns in Test Table{}: {}".format(testtab,cnames2))
        return False

    for excludecol in excludecols:
        if (excludecol not in cnames) and (excludecol not in cnames2):
            logging.warning("Column {} Not in {} or {}. Will Continue without Checking against this column".format(excludecol,referencetab,testtab))
            logging.debug("Available columns in Reference Table {}: {}".format(referencetab,cnames))
            logging.debug("Available columns in Test Table{}: {}".format(testtab,cnames2))

    try:
        for cname in cnames:
            if cname in excludecols:
                continue
            
            logging.info("\nTesting column: {}".format(cname))
            
            a = 0
            try:
                a = tb.getcol(cname,startrow=startrow,nrow=nrow,rowincr=rowincr)
            except:
                tb.getcol(cname)
                rval = False
                logging.critical('Error accessing column ', cname, ' in table ', referencetab)
                logging.critical(sys.exc_info()[0])
                break

            b = 0
            try:
                b = tb2.getcol(cname,startrow=startrow,nrow=nrow,rowincr=rowincr)
            except:
                rval = False
                logging.critical('Error accessing column ', cname, ' in table ', testtab)
                logging.critical(sys.exc_info()[0])
                break

            if not (len(a)==len(b)):
                logging.error('Column {} has different length in tables {} and {}'.format(cname, referencetab, testtab))
                logging.error(a)
                logging.error(b)
                rval = False
                break
            else:
                differs = False
                if not (a==b).all():
                    for i in range(0,len(a)):
                        if (isinstance(a[i],float)):
                            if ((mode=="percentage") and (abs(a[i]-b[i]) > tolerance*abs(a[i]))) or ((mode=="absolute") and (abs(a[i]-b[i]) > tolerance)):
                                print("Column " + cname + " differs")
                                print("Row=" + str(i))
                                print("Reference file value: " + str(a[i]))
                                print("Input file value: " + str(b[i]))
                                if (mode=="percentage"):
                                    print("Tolerance is {0}%; observed difference was {1} %".format (tolerance * 100, 100*abs(a[i]-b[i])/abs(a[i])))
                                else:
                                    print("Absolute tolerance is {0}; observed difference: {1}".format (tolerance, (abs(a[i]-b[i]))))
                                differs = True
                                rval = False
                                break
                        elif (isinstance(a[i],int) or isinstance(a[i],numpy.int32)):
                            if (abs(a[i]-b[i]) > 0):
                                print("Column " + cname + " differs")
                                print("Row=" + str(i))
                                print("Reference file value: " + str(a[i]))
                                print("Input file value: " + str(b[i]))
                                if (mode=="percentage"):
                                    print("tolerance in % should be " + str(100*abs(a[i]-b[i])/abs(a[i])))
                                else:
                                    print("absolute tolerance should be " + str(abs(a[i]-b[i])))
                                differs = True
                                rval = False
                                break
                        elif (isinstance(a[i],str) or isinstance(a[i],numpy.bool_)):
                            if not (a[i]==b[i]):
                                print("Column " + c + " differs")
                                print("Row=" + str(i))
                                print("Reference file value: " + str(a[i]))
                                print("Input file value: " + str(b[i]))
                                if (mode=="percentage"):   
                                    print("tolerance in % should be " + str(100*abs(a[i]-b[i])/abs(a[i])))
                                else:
                                    print("absolute tolerance should be " + str(abs(a[i]-b[i])))
                                differs = True
                                rval = False
                                break
                        elif (isinstance(a[i],list)) or (isinstance(a[i],numpy.ndarray)):
                            for j in range(0,len(a[i])):
                                if differs: break
                                if ((isinstance(a[i][j],float)) or (isinstance(a[i][j],int))):
                                    if ((mode=="percentage") and (abs(a[i][j]-b[i][j]) > tolerance*abs(a[i][j]))) or ((mode=="absolute") and (abs(a[i][j]-b[i][j]) > tolerance)):
                                        print("Column " + c + " differs")
                                        print("(Row,Element)=(" + str(j) + "," + str(i) + ")")
                                        print("Reference file value: " + str(a[i][j]))
                                        print("Input file value: " + str(b[i][j]))
                                        if (mode=="percentage"):
                                            print("Tolerance in % should be " + str(100*abs(a[i][j]-b[i][j])/abs(a[i][j])))
                                        else:
                                            print("Absolute tolerance should be " + str(abs(a[i][j]-b[i][j])))
                                        differs = True
                                        rval = False
                                        break
                                elif (isinstance(a[i][j],list)) or (isinstance(a[i][j],numpy.ndarray)):
                                    it = range(0,len(a[i][j]))
                                    if mode=="percentage":
                                        diff = numpy.abs(numpy.subtract(a[i][j], b[i][j])) > tolerance * numpy.abs(a[i][j])
                                        it = numpy.where(diff)[0]
                                    elif (mode=="absolute"):
                                        diff = numpy.abs(numpy.subtract(a[i][j], b[i][j])) > tolerance
                                        it = numpy.where(diff)[0]
                                    for k in it:
                                        if differs: break
                                        if ( ((mode=="percentage") and (abs(a[i][j][k]-b[i][j][k]) > tolerance*abs(a[i][j][k]))) \
                                                 or ((mode=="absolute") and (abs(a[i][j][k]-b[i][j][k]) > tolerance)) \
                                                 or ((mode=="phaseabsdeg") and (phasediffabsdeg(a[i][j][k],b[i][j][k])>tolerance)) \
                                                 ):
                                            print("Column " + c + " differs")
                                            print("(Row,Channel,Corr)=(" + str(k) + "," + str(j) + "," + str(i) + ")")
                                            print("Reference file value: " + str(a[i][j][k]))
                                            print("Input file value: " + str(b[i][j][k]))
                                            if (mode=="percentage"):
                                                print("Tolerance in % should be " + str(100*abs(a[i][j][k]-b[i][j][k])/abs(a[i][j][k])))
                                            elif (mode=="absolute"):
                                                print("Absolute tolerance should be " + str(abs(a[i][j][k]-b[i][j][k])))
                                            elif (mode=="phaseabsdeg"):
                                                print("Phase tolerance in degrees should be " + str(phasediffabsdeg(a[i][j][k],b[i][j][k])))
                                            else:
                                                print("Unknown comparison mode: ",mode)
                                            differs = True
                                            rval = False
                                            break

                        else:
                            print("Unknown data type: ",type(a[i]))
                            differs = True
                            rval = False
                            break
                
                if not differs: print("Column " + cname + " PASSED")
    finally:
        tb.close()
        tb2.close()

    logging.debug("compare_CASA_tables(referencetab = {}, testtab = {}): {}".format(referencetab,testtab, rval))
    return rval

def compare_files(file1, file2, shallow=False):
    '''
    compare_files - Compare two Files.
       @param file1       --> a reference file
       @param file2       --> a file to verify
       @param shallow     --> If shallow is true, files with identical os.stat() signatures are taken to be equal. Otherwise, the contents of the files are compared.

       @return: True if file1 & file2 seem equal, False otherwise
    '''
    logging.info("Comparing {} to {}".format(file1, file2))
    filecmp.clear_cache()
    return filecmp.cmp(file1, file2, shallow=shallow)

def compare_caltables(table1, table2, cols=[], rtol=8e-7, atol=1e-8):
    '''
    compare_caltables - Compare two caltables.
       @param table1       --> a reference table
       @param table2       --> a table to verify
       @param cols         --> the name of cols to compare (list). Leave Blank For All
       @param rtol         --> The relative tolerance parameter
       @param atol         --> The absolute tolerance parameter

       @return: True if table1 == table2 else False
    '''
    logging.info("Comparing {} to {}".format(table1, table2))
    logging.debug("Executing: compare_caltables(table1 = {}, table2 = {}, cols={}, rtol={}, atol={})".format(table1, table2, cols, rtol, atol))
    tableVal1 = {}
    tableVal2 = {}

    tb.open(table1)
    colname1 = tb.colnames()

    for col in colname1:
        try:
            tableVal1[col] = tb.getcol(col)
        except RuntimeError:
            pass
    tb.close()

    tb2.open(table2)
    colname2 = tb2.colnames()

    for col in colname2:
        try:
            tableVal2[col] = tb2.getcol(col)
        except RuntimeError:
            pass
    tb2.close()

    truthDict = {}

    for col in tableVal1.keys():
        try:
            if tableVal1[col].dtype == '<U16' and tableVal2[col].dtype == '<U16': #16-character string 
                truthDict[col] = numpy.array_equal(tableVal1[col],tableVal2[col])
            else:
                truthDict[col] = numpy.isclose(tableVal1[col], tableVal2[col], rtol=rtol, atol=atol)
        except TypeError:
            print(col, 'ERROR in finding truth value')
            casalog.post(message=col+': ERROR in determining the truth value')

    if len(cols) == 0:
        truths = [[x, numpy.all(truthDict[x] == True)] for x in truthDict.keys()]
    else:
        truths = [[x, numpy.all(truthDict[x] == True)] for x in cols]

    #Check that All Options are True
    for key in truthDict.keys():
        if isinstance(truthDict[key], bool):
            if not truthDict[key]:
                logging.info("{0} in caltables do not match".format(key))
                return False
        elif isinstance(truthDict[key], numpy.ndarray):
            if not numpy.all(truthDict[key]):
                return False
        else:
            logging.info('ERROR in finding truth value for Column: {}'.format(key))
            return False

    return True

def compare_dictionaries( dictionary1, dictionary2, skipkeys = [], rtol=8e-7, atol=1e-8):
    '''
    compare_dictionaries - compare two dictionaries
       @param dictionary1  --> the dictionary which is assumed to be correct
       @param dictionary2  --> the dictionary which is to be compared
       @param skipkeys     --> list of keys which are to be ignored
       @param rtol         --> The relative tolerance parameter
       @param atol         --> The absolute tolerance parameter

       @return: True if dictionary1 == dictionary2 else False
    '''
    if not isinstance(skipkeys, list):
        logging.error("skipkeys not in correct format")
        raise TypeError("skipkeys must be a list")

    key_list_1 = sorted(list(dictionary1.keys()))
    key_list_2 = sorted(list(dictionary2.keys()))

    #Checks if Keys are the same
    if key_list_1 != key_list_2:
        logging.debug("Keys Do Not Match")
        return False
    
    for key in key_list_1:
        if key in skipkeys:
            continue
        # Compare Numpy Arrays
        if isinstance(dictionary1[key], numpy.ndarray) and isinstance(dictionary2[key], numpy.ndarray):
            """
                For finite values, isclose uses the following equation to test whether two floating point values are equivalent.
                absolute(a - b) <= (atol + rtol * absolute(b))
            """
            if any( val == False for val in numpy.isclose(dictionary1[key], dictionary2[key], rtol=rtol, atol=atol, equal_nan=False)):
                logging.info("{0}:{1} != {0}:{2}".format(key,dictionary1[key],dictionary2[key]))
                return False

        # Compare Strings
        elif isinstance(dictionary1[key], str) and isinstance(dictionary2[key], str):

            if (dictionary1[key] == dictionary2[key]) and (dictionary1[key] is dictionary2[key]):
                pass
            else:
                logging.info("{0}:{1} != {0}:{2}".format(key,dictionary1[key],dictionary2[key]))
                return False

        # Compare lists
        elif isinstance(dictionary1[key], list) and isinstance(dictionary2[key], list):

            if dictionary1[key] != dictionary2[key]:
                logging.info("{0}:{1} != {0}:{2}".format(key,dictionary1[key],dictionary2[key]))
                return False

        # Compare Numerics
        elif isinstance(dictionary1[key], numbers.Number) and isinstance(dictionary2[key], numbers.Number):
            """
            rel_tol is the relative tolerance â€“ it is the maximum allowed difference between a and b, relative to the larger absolute value of a or b. 
            For example, to set a tolerance of 5%, pass rel_tol=0.05. The default tolerance is 1e-09, which assures that the two values are the same within about 9 decimal digits. rel_tol must be greater than zero.

            abs_tol is the minimum absolute tolerance â€“ useful for comparisons near zero. abs_tol must be at least zero.

            If no errors occur, the result will be: abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol).
            """

            if not math.isclose(dictionary1[key],dictionary2[key],rel_tol = rtol, abs_tol=atol):
                logging.info("{0}:{1} != {0}:{2}".format(key,dictionary1[key],dictionary2[key]))
                return False
        else:
            try:
                if dictionary1[key] != dictionary2[key]:
                    return False
            except:
                logging.error("Error in Comparing {0}:{1} != {0}:{2}".format(key,dictionary1[key],dictionary2[key]))
                return False
    return True

def compare_directories( directory1, directory2):
    '''
        Compare two directories recursively. Files in each directory are
        assumed to be equal if their names and contents are equal.

        @param directory1: First directory path
        @param directory2: Second directory path

        @return: True if the directory trees are the same and 
            there were no errors while accessing the directories or files, 
            False otherwise.
    '''
    dirs_cmp = filecmp.dircmp(directory1, directory2)
    if len(dirs_cmp.left_only)>0 or len(dirs_cmp.right_only)>0 or \
        len(dirs_cmp.funny_files)>0:
        return False
    (_, mismatch, errors) =  filecmp.cmpfiles(
        directory1, directory2, dirs_cmp.common_files, shallow=False)
    if len(mismatch)>0 or len(errors)>0:
        return False
    for common_dir in dirs_cmp.common_dirs:
        new_directory1 = os.path.join(directory1, common_dir)
        new_directory2 = os.path.join(directory2, common_dir)
        if not compare_directories(new_directory1, new_directory2):
            return False
    return True

def check_pixels(imagename='', loc=None, refval=None,  rtol=1e-05, atol=1e-08):
    '''
        Check pixels in an image to a specified reference value
       
        @param imagename: input image file
        @param loc: The index of the image to compare to the refval
        @param refval: The reference value to compare the selected pixel(s) to
        @param rtol: The relative tolerance used in the numpy.isclose function
        @param atol: The absolute tolerance used in the numpy.isclose function
       
        @return: True if the shape and value of the refval and selected pixel match.

    '''
    if not isinstance(loc,str):
        raise TypeError('Please give target location in string list format "20,30,2:4"')

    if os.path.exists(imagename):
        tb.open(imagename)
        image = tb.getcol('map')
        tb.close()
        
        if type(refval) != type(None):

            index = []
            to_slice = loc.split(',')
            
            for item in to_slice:
                if ':' not in item:
                    index.append(int(item))
                else:
                    item_split = item.split(':')
                    index.append(slice(int(item_split[0]),int(item_split[1])))
                    
            selected_slice = image[tuple(index)]
            
            if numpy.shape(selected_slice) != numpy.shape(refval):
                    logging.warning('Please check that the shape of the reference and selected slice are the same')
                    return False
            
            isequal = numpy.isclose(selected_slice, refval, rtol=rtol, atol=atol)

            logging.info("For pixel value check the obtained value was {}. The expected value was {} with a tolerance of {}. test success = {}.".format(selected_slice, refval, atol, isequal))
            return numpy.all(isequal == True)

        else:
            logging.warning('Please provide a refernce value to compare against')
        
    else:
        logging.warning('Not a valid Image name')
        
def compare_pixel_value(imagename=None, refimage=None, loc=None):
    '''
        Compare two images at a certain reference pixel
       
        @param imagename: Name of the image to be compared to a reference image
        @param refimage: Image to be compared against
        @param loc: The slice or pixel index to compare between the two images
       
        @return: True if the pixel values match at the provided index or slice. Returns False otherwise
    '''
    if imagename != None and refimage != None:
        
        if type(loc) == type(''):
            
            tb.open(imagename)
            image1 = tb.getcol('map')
            tb.close()
            
            tb.open(refimage)
            image2 = tb.getcol('map')
            tb.close()
            
            index = []
            to_slice = loc.split(',')
            # get index from the string array
            for item in to_slice:
                if ':' not in item:
                    index.append(int(item))
                else:
                    item_split = item.split(':')
                    index.append(slice(int(item_split[0]),int(item_split[1])))
                    
            selected_slice1 = image1[tuple(index)]
            selected_slice2 = image2[tuple(index)]
            
            isequal = numpy.isclose(selected_slice1, selected_slice2, rtol=1e-05, atol=1e-08)
            return numpy.all(isequal == True)
            
        else:
            logging.warning('Please give target location in string list format ("20,30,2:4")')
    else:
        logging.warning('Please provide both an image and reference image')
    
def compare_pixel_mask(maskname='', refmask=None, refval=None, loc=None):
    '''
        Compare to masks or mask values to a reference value
       
        @param maskname: The name of the maskfile to compare to either a reference mask file or value
        @param refmask: The reference mask image to be compared to
        @param refval: The reference value to compare the selected pixel(s) of the maskfile to
        @param loc: The index or slice of the mask image to compare to a refvalue.
       
        @return: True if the refmask and mask file are identical or if the selected slice of the mask file matches the refval
    '''
    
    if os.path.exists(maskname):
        if refmask == None and refval == None:
            logging.warning('Please select a mask or region to use for comparison')

        elif refmask != None and refval == None:
            # if comparing a refmask compare the values in the table
            if os.path.exists(refmask):
                tb.open(maskname)
                mask1 = tb.getcol('PagedArray')
                tb.close()

                tb.open(refmask)
                mask2 = tb.getcol('PagedArray')
                tb.close()
                
                return np.all(mask1 == mask2)
                
            else:
                logging.warning('Invalid refmask file name')
                
        elif refmask == None and refval != None:
            # If using a reference value compare the value/shape to the selected slice
            if type(loc) == type(''):
                
                tb.open(maskname)
                image = tb.getcol('PagedArray')
                tb.close()
                
                index = []
                to_slice = loc.split(',')
                # get index from the string array
                for item in to_slice:
                    if ':' not in item:
                        index.append(int(item))
                    else:
                        item_split = item.split(':')
                        index.append(slice(int(item_split[0]),int(item_split[1])))
                        
                selected_slice = image[tuple(index)]
                # return false if the shapes don't match up
                if numpy.shape(selected_slice) != numpy.shape(refval):
                    logging.warning('Please check that the shape of the reference and selected slice are the same')
                    return False
                
                isequal = numpy.all(selected_slice == refval)
                
                return isequal
                
            else:
                logging.warning('Please give target location in string list format ("20,30,2:4")')
        else:
            logging.warning('Please provide only a referance value or reference mask, not both')
    else:
        logging.warning('Invalid mask file name')


def add_to_dict(self, output=None, filename=None, dataset=None, **kwargs):
    '''
        This function adds key value pairs to a provided dictionary. Any additional keys and values can be added as keyword arguments to this function
       
        @param output: This is the dictionary that the key-value pairs will be appended to
        @param filename: This is the name of the test script file
        @param dataset: This is the name of the dataset used when executing this test case
       
        @return: Nothing is returned, the output dict is modified by this function
    '''
    testcase = unittest.TestCase.id(self)
    test_split = testcase.split('.')
    test_case = test_split[-1]
    taskname = test_split[1].split('_')[0]
    
    rerun = "python {} {}.{}".format(filename, test_split[1], test_split[2])
    
    current_case = None
    
    func_calls = []
    values = {key:kwargs[key] for key in kwargs}
    
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('def test_'):
                if line.split()[1][:-7].endswith(test_case):
                    current_case = test_case
                else:
                    current_case = None
                    
            if current_case == test_case:
                if line.startswith(taskname):
                    params = line.split(',')[1::]
                    call = "{}({},{})".format(taskname, dataset, ','.join(params))
                    func_calls.append(call)
                
    values['taskcall'] = func_calls
    values['rerun'] = rerun
    values['images'] = [ ]
    output[test_case] = values

def topickle(input_dict, picklefile):
    '''
        Add a new dictionary into the existing pickle file
       
        @param input_dict: The dictionary object to add to the pickle file
        @param picklefile: The picklefile containing a dictionary to be appended to
       
        @return: Nothing is returned by this function
    '''
    pickle_read = open(picklefile, 'rb')
    pickle_dict = pickle.load(pickle_read)
    # Make sure that the pickle file contains a dictionary
    if type(pickle_dict) != type({}):
        logging.warning('The pickle file is not a dictionary')
    # Add to the dictionary in the pickle file
    for item in list(input_dict.keys()):
        pickle_dict[item] = input_dict[item]
    # Re-write the pickle file with the new dictionary
    with open(picklefile, 'wb') as fout:
        pickle.dump(pickle_dict, fout)

def default_CASA_tasks():
    '''
    default_CASA_tasks - Default Casa Tasks
    Delete all *.last files and restore tasks to default 

       Returns 
    '''
    # Get a list of all files in directory
    for rootDir, subdirs, filenames in os.walk(os.getcwd()):
        # Find the files that matches the given patterm
        for filename in fnmatch.filter(filenames, '*.last'):
            try:
                os.remove(os.path.join(rootDir, filename))
            except OSError:
                logging.error("Error while deleting file")
    if casa5:
        for task in casa6tasks:
            default(task)
        for task in miscellaneous_tasks:
            default(task)
    return

def get_directory_size(directory):
    '''
    get_directory_size - Return the size of a directory in bytes
       directory  --> the directory which is to be summed

       Returns Return the size, in bytes, of directory
    '''
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(directory):
        for filename in filenames:
            fp = os.path.join(dirpath, filename)
            total_size += os.path.getsize(fp)
    return total_size

def get_table_column(table, colname):
    '''Return the requested variable column
       table    --> name of table or MS
       colname  --> column name
    Return the column as a dictionary
    '''

    col = {}
    tb.open(table)
    if tb.isvarcol(colname):
        col = tb.getvarcol(colname)
    else:
        logging.error("Error Returning Column {}".format(colname))
        return None

    tb.close()
    return col


def get_caltable_column(caltable, colname='CPARAM'):
    ''' Open a caltable and get the provided column
       caltable    --> name of cal table
       colname  --> column name
    Return the column as a dictionary
    '''
    tb.open(caltable)
    outtable = tb.getcol(colname)
    tb.close()
    return outtable

def get_column_shape(tab,col,start_row=0,nrow=1,row_inc=1):
    ''' 
    Get the shape of the given column.
    Keyword arguments:
        tab        --    input table or MS
        col        --    column to get the shape
        start_row  --    start row (default 0)
        nrow       --    number of rows to read (default 1)
        row_inc    --    increment of rows to read (default 1)
        
        Return a list of strings with the shape of each row in the column.
    
    '''

    col_shape = []
    try:
        try:
            tb.open(tab)
            col_shape = tb.getcolshapestring(col,start_row,nrow,row_inc)
        except:
            print('Cannot get shape of col %s from table %s '%(col,tab))

    finally:
        tb.close()
            
    return col_shape

def check_plotfile(plotfileName, min_size, max_size=None):
    ''' 
    Check if plotfile generated is cprrect size
        plotfileName --> Name of plotted Image
        min_size -- > Min Size of image
        max_size --> Max Size of image
    
        Return : True if image size > min_size ( and < max_size if max_size is provided )
    '''
    val = False
    if os.path.isfile(plotfileName):
        plotSize = os.path.getsize(plotfileName) # Return the size, in bytes, of path.
        logging.info( '{} file size is: {}'.format( plotfileName, plotSize))
        if plotSize > min_size:
            val = True
        if maxSize not None:
            if not plotSize < max_size:
                val = False
    else:
        logging.critical("Plot was not created")
    
    return val

def generate_weblog(task,dictionary):
    """Generate Test Summary Weblog

    Example: 
        generate_weblog("taskname", dictionary)
    """
    global html
    html = open("test_{}_weblog.html".format(task.lower()), 'w')
    Weblog(task, dictionary).generate_weblog()
    html.close()
############################################################################################
##################################       Decorators       ##################################
############################################################################################

#import casaTestHelper
#@casaTestHelper.skipIfMissingModule
def skipIfMissingModule(required_module,strict=False):
    '''
    Decorator: skip test if specified module is not avaliable 

    Example:
        @casaTestHelper.skipIfMissingModule('astropy')
        def test_test(self):
    '''
    import os
    try:
        __import__(required_module)
        flag = True
    except ImportError:
        flag = False
    def deco(function):
        if not CASA6:
            return deco
        
        def wrapper(self, *args, **kwargs):
            if not flag:
                # If there is a strict flag run the tests as normal
                print(sys.argv)
                if stri/ct:
                    function(self)
                    pass
                else:
                    # Module ImportError and no strict flag
                    self.skipTest("ModuleNotFoundError: No module named '{}'".format(required_module))
            else:
                function(self)
        return wrapper
    return deco

#import casaTestHelper
#@casaTestHelper.time_execution
def time_execution(out_dict):
    def time_decorator(function):
        '''
        Decorator: time execution of test 

        Example:
            @casaTestHelper.time_execution
            def test_test(self):
        '''
        @wraps(function)
        def function_timer(*args, **kwargs):
            t0 = time.time()
            result = function(*args, **kwargs)
            t1 = time.time()
            #print ("Total time running %s: %s seconds" % (function.__name__, str(t1-t0)))
            casalog.post("Total time running {}: {} seconds".format(function.__name__, str(t1-t0)))
            out_dict[function.__name__]['runtime'] = t1-t0
            return result
            
        return function_timer
    return time_decorator
    
def time_execution_alternative(out_dict):
    # TODO Ver if this is the better option
    def time_decorator(function):
        '''
        Decorator: time execution of test 

        Example:
            @casaTestHelper.time_execution
            def test_test(self):
        '''
        @wraps(function)
        def function_timer(*args, **kwargs):
            failed = False
            result = None
            t0 = time.time()
            try:
                result = function(*args, **kwargs)
            except:
                failed=True
                t1 = time.time()
                out_dict[function.__name__]['runtime'] = t1-t0
                casalog.post("Total time running {}: {} seconds".format(function.__name__, str(t1-t0)))
                out_dict[function.__name__]['status'] = False
                raise
                
            t1 = time.time()
            #print ("Total time running %s: %s seconds" % (function.__name__, str(t1-t0)))
            casalog.post("Total time running {}: {} seconds".format(function.__name__, str(t1-t0)))
            out_dict[function.__name__]['runtime'] = t1-t0
            out_dict[function.__name__]['status'] = True

            
            return result
            
        return function_timer
    return time_decorator
    
def cpu_usage(out_dict):
    def cpu_decorator(function):
        @wraps(function)
        def function_usage(*args, **kwargs):
            
            process = psutil.Process(os.getpid())
            snapshot1 = process.memory_info()
            open_files1 = process.open_files()
            num_file_descriptors1 = process.num_fds()

            #print ("Function: {}, {} MBs".format(function.__name__, megs1))
            #print ("Function: {}, Open Files: {}".format(function.__name__, open_files1))
            #print ("Function: {}, num_file_descriptors: {}".format(function.__name__, num_file_descriptors1))

            result = function(*args, **kwargs)

            process = psutil.Process(os.getpid())
            snapshot2 = process.memory_info()
            open_files2 = process.open_files()
            num_file_descriptors2 = process.num_fds()

            #print ("Function: {}, {} MBs".format(function.__name__, megs2))
            #print ("Function: {}, Open Files: {}".format(function.__name__, open_files2))
            #print ("Function: {}, num_file_descriptors: {}".format(function.__name__, num_file_descriptors2))
            #print('{:.2f} MB\n'.format(process.memory_info().rss / 1024 / 1024))

            #print ("Total Mem Info { }: {:.2f} MB".format(function.__name__,(process.memory_info().rss) / 1024 / 1024 ))
            out_dict[function.__name__]['cpu_usage'] = { "number of file descriptors opened" :  num_file_descriptors2 - num_file_descriptors1,
                                                      "Open files" : open_files2,
                                                      "Pre Memory Snapshot (bytes)" : snapshot1,
                                                      "Post Memory Snapshot (bytes)" : snapshot2
                                                    }
            return result
        return function_usage
    return cpu_decorator

def peakmem(out_dict):
    #TODO: https://pytracemalloc.readthedocs.io/examples.html
    # Only for python3.4+
    def mem_decorator(function):
        @wraps(function)
        def function_mem(*args, **kwargs):
            tracemalloc.clear_traces()
            tracemalloc.start()
            snapshot1 = tracemalloc.take_snapshot() # Snapshot of traces of memory blocks allocated by Python.

            result = function(*args, **kwargs)

            snapshot2 = tracemalloc.take_snapshot()
            peakmem = ("{} MiB".format(tracemalloc.get_traced_memory()[1] / 1024 /1024)) #Get the current size and peak size of memory blocks traced by the tracemalloc module as a tuple: (current: int, peak: int)
            tracemalloc.stop()
            top_stats = snapshot2.compare_to(snapshot1, 'lineno') # Compute the differences with an old snapshot.
            out_dict[function.__name__]['peakmem'] = peakmem 
            out_dict[function.__name__]['memleaks'] = top_stats[:10] #
            
            return result
        return function_mem
    return mem_decorator
    
    
def stats_dict(out_dict):
    def stats_decorator(function):
        @time_execution_alternative(out_dict)
        @cpu_usage(out_dict)
        @peakmem(out_dict)
        @wraps(function)
        def all_wrapped(*args, **kwargs):
            return function(*args, **kwargs)
        return all_wrapped
    return stats_decorator



