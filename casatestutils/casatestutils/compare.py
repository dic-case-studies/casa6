import logging
import sys
import filecmp
import os
import numpy
import numbers
import shutil
try:
    # CASA 6
    logging.debug("Importing CASAtools")
    import casatools
    _tb = casatools.table()
    _tb2 = casatools.table()
    # _casa6 = True

except ImportError:
    # CASA 5
    logging.debug("Import casa6 errors. Trying CASA5...")
    from taskinit import tbtool
    _tb = tbtool()
    _tb2 = tbtool()
    # _casa5 = True

#ignore_subversion = shutil.ignore_patterns('.svn')
################        ##################
################        ##################
################        ##################


#class TableCacheValidator(object):
#    def __init__(self):
#        self.original_cache = get_table_cache()
#
#    def validate(self):
#        cache = get_table_cache()
#        #print 'original {} current {}'.format(self.original_cache, cache)
#        return len(cache) == 0 or cache == self.original_cache

#class DictDiffer(object):
#    """
#    Calculate the difference between two dictionaries as:
#    (1) items added
#    (2) items removed
#    (3) keys same in both but changed values
#    (4) keys same in both and unchanged values
#    Example:
#            mydiff = DictDiffer(dict1, dict2)
#            mydiff.changed()  # to show what has changed
#    """
#    def __init__(self, current_dict, past_dict):
#        self.current_dict, self.past_dict = current_dict, past_dict
#        self.set_current, self.set_past = set(current_dict.keys()), set(past_dict.keys())
#        self.intersect = self.set_current.intersection(self.set_past)
#    def added(self):
#        return self.set_current - self.intersect 
#    def removed(self):
#        return self.set_past - self.intersect 
#    def changed(self):
#        return set(o for o in self.intersect if self.past_dict[o] != self.current_dict[o])            
#    def unchanged(self):
#        return set(o for o in self.intersect if self.past_dict[o] == self.current_dict[o])



################        ##################
################        ##################
################        ##################


def compare_CASA_var_col_tables(referencetab, testtab, varcol, tolerance=0.0):
    '''
    originally: compVarColTables(referencetab, testtab, varcol, tolerance=0.)
    compare_CASA_variable_cols - Compare a variable column of two CASA tables.
       @param referencetab  --> a reference table
       @param testtab       --> a table to verify
       @param varcol        --> the name of a variable column (str)
       @param tolerance     --> Tolerance
       @return: True if reference tab == test table else False
    '''
    logging.info("Comparing Column: {} within {} and {}".format(varcol,referencetab, testtab))
    logging.debug("Executing: compare_CASA_variable_cols(referencetab={},testtab={}, varcol={}, tolerance={})".format(referencetab, testtab, varcol, tolerance))
    retval = True

    _tb.open(referencetab)
    cnames = _tb.colnames()

    _tb2.open(testtab)
    col = varcol
    if _tb.isvarcol(col) and _tb2.isvarcol(col):
        try:
            # First check
            if _tb.nrows() != _tb2.nrows():
                #print('Length of %s differ from %s, %s!=%s'%(referencetab,testtab,len(rk),len(tk)))
                logging.error('Length of {} differ from {}, {} != {}'.format(referencetab,testtab,_tb.nrows(),_tb2.nrows()))
                retval = False
            else:
                for therow in range(_tb.nrows()):
                    rdata = _tb.getcell(col,therow)
                    tdata = _tb2.getcell(col,therow)
                    if not rdata.all()==tdata.all():
                        if (tolerance>0.):
                            differs=False
                            for j in range(0,len(rdata)):
                                if ((isinstance(rdata[j],float)) or (isinstance(rdata[j],int))):
                                    if (abs(rdata[j]-tdata[j]) > tolerance*abs(rdata[j]+tdata[j])):
#                                        print('Column ', col,' differs in tables ', referencetab, ' and ', testtab)
#                                        print(therow, j)
#                                        print(rdata[j])
#                                        print(tdata[j])
                                        differs = True
                                elif (isinstance(rdata[j],list)) or (isinstance(rdata[j],np.ndarray)):
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
            _tb.close()
            _tb2.close()

    else:
        print('Columns are not varcolumns.')
        retval = False

    if retval:
        print('Column %s of %s and %s agree'%(col,referencetab, testtab))
        
    return retval

def compare_CASA_tables(referencetab, testtab, excludecols, tolerance=0.001, mode="percentage", startrow = 0, nrow = -1, rowincr = 1):
    '''
    originally: testhelpers.compTables(referencetab, testtab, excludecols, tolerance=0.001, mode="percentage", startrow = 0, nrow = -1, rowincr = 1)
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

    if excludecols is None:
        excludecols = []
    if not isinstance(excludecols, list):
        logging.error("excludecols not in correct format")
        raise TypeError("excludecols must be a list")
    if referencetab.endswith(".cal") or testtab.endswith(".cal"):
        logging.warning("WARNING: Will compare caltables using compare_caltables")
        return compare_caltables(referencetab, testtab, cols= excludecols, rtol=8e-7, atol=1e-8)
    ##### Begin: Tempory Fix
    if len(excludecols) == 0:
        excludecols = ["FLAG_CATEGORY"]
    else:
        excludecols.append('FLAG_CATEGORY')
    """
    #TODO: Fix Error in checking FLAG_CATEGORY
        _tb.getcol("FLAG_CATEGORY")
            SEVERE  getcol::FLAG_CATEGORY   Exception Reported: Table DataManager error: Invalid operation: TSM: no array in row 0 of column FLAG_CATEGORY in **
            RuntimeError: Table DataManager error: Invalid operation: TSM: no array in row 0 of column FLAG_CATEGORY in **
    """
    ##### End: Tempory Fix
    rval = True
    _tb.open(referencetab)
    cnames = _tb.colnames()

    _tb2.open(testtab)
    cnames2 = _tb2.colnames()
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
        for c in cnames:
            if c in excludecols:
                continue

            print("\nTesting column {}".format(c))

            a = 0
            try:
                a = _tb.getcol(c,startrow=startrow,nrow=nrow,rowincr=rowincr)
            except:
                rval = False
                print('Error accessing column ', c, ' in table ', referencetab)
                print(sys.exc_info()[0])
                break

            b = 0
            try:
                b = _tb2.getcol(c,startrow=startrow,nrow=nrow,rowincr=rowincr)
            except:
                rval = False
                print('Error accessing column ', c, ' in table ', testtab)
                print(sys.exc_info()[0])
                break

            if not (len(a)==len(b)):
                print('Column ',c,' has different length in tables ', referencetab, ' and ', testtab)
                print(a)
                print(b)
                rval = False
                break
            else:
                differs = False
                if not (a==b).all():
                    for i in range(0,len(a)):
                        if (isinstance(a[i],float)):
                            if ((mode=="percentage") and (abs(a[i]-b[i]) > tolerance*abs(a[i]))) or ((mode=="absolute") and (abs(a[i]-b[i]) > tolerance)):
                                print("Column " + c + " differs")
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

                if not differs: print("Column " + c + " PASSED")
    finally:
        _tb.close()
        _tb2.close()

    logging.debug("compare_CASA_tables(referencetab = {}, testtab = {}): {}".format(referencetab,testtab, rval))
    return rval

def compare_files( file1, file2, shallow=False):
    '''
    compare_files - Compare two Files.
       @param file1       --> a reference file
       @param file2       --> a file to verify
       @param shallow     --> If shallow is true, files with identical os.stat() signatures are taken to be equal. Otherwise, the contents of the files are compared.
       @return: True if file1 & file2 seem equal, False otherwise
    '''
    logging.info("Comparing {} to {}".format(file1, file2))
    logging.debug("Executing: compare_files(file1 = {}, file2 = {}, shallow = {})".format(file1, file2, shallow))
    if sys.version_info > (3,0):
        filecmp.clear_cache()
    return filecmp.cmp(file1, file2, shallow=shallow)

def compare_caltables( table1, table2, cols=None, rtol=8e-7, atol=1e-8):
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
    if cols is None:
        cols = []
    tableVal1 = {}
    tableVal2 = {}
    _tb.open(table1)
    colname1 = _tb.colnames()
    for col in colname1:
        try:
            tableVal1[col] = _tb.getcol(col)
        except RuntimeError:
            pass
    _tb.close()
    _tb2.open(table2)
    colname2 = _tb2.colnames()
    for col in colname2:
        try:
            tableVal2[col] = _tb2.getcol(col)
        except RuntimeError:
            pass
    _tb2.close()
    truthDict = {}
    for col in tableVal1.keys():
        logging.debug("Column: {}, dtype: {}".format(col, tableVal1[col].dtype))
        try:
            if numpy.issubdtype(tableVal1[col].dtype, numpy.number):
                truthDict[col] = numpy.isclose(tableVal1[col], tableVal2[col], rtol=rtol, atol=atol)
            else:
                # Compare Non Numeric Types
                truthDict[col] = numpy.array_equal(tableVal1[col],tableVal2[col])
        except:
            print(col, 'ERROR in determining the truth value')
            #casalog.post(message=col+': ERROR in determining the truth value')
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

def compare_dictionaries(  dictionary1, dictionary2, skipkeys = None, rtol=8e-7, atol=1e-8):
    '''
    compare_dictionaries - compare two dictionaries
       Dictionaries will fail when 1st instance of a failure
       @param dictionary1  --> the dictionary which is assumed to be correct
       @param dictionary2  --> the dictionary which is to be compared
       @param skipkeys     --> list of keys which are to be ignored
       @param rtol         --> The relative tolerance parameter
       @param atol         --> The absolute tolerance parameter
       @return: True if dictionary1 == dictionary2 else False
    '''
    if skipkeys is None:
        skipkeys = []
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
            if numpy.issubdtype(dictionary1[key].dtype, numpy.number) and numpy.issubdtype(dictionary2[key].dtype, numpy.number):
                if any( val == False for val in numpy.isclose(dictionary1[key], dictionary2[key], rtol=rtol, atol=atol, equal_nan=False)):
                    logging.info("{0}:{1} != {0}:{2}".format(key,dictionary1[key],dictionary2[key]))
                    return False
            else:
                if any( val == False for val in numpy.array_equal(dictionary1[key], dictionary2[key])):
                    logging.info("{0}:{1} != {0}:{2}".format(key,dictionary1[key],dictionary2[key]))
                    return False
        # Compare Strings
        elif isinstance(dictionary1[key], str) and isinstance(dictionary2[key], str):
            if (dictionary1[key] == dictionary2[key]):
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
            rel_tol is the relative tolerance :  it is the maximum allowed difference between a and b, relative to the larger absolute value of a or b.
            For example, to set a tolerance of 5%, pass rel_tol=0.05. The default tolerance is 1e-09, which assures that the two values are the same within about 9 decimal digits. rel_tol must be greater than zero.
            abs_tol is the minimum absolute tolerance : useful for comparisons near zero. abs_tol must be at least zero.
            If no errors occur, the result will be: abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol).
            """
            if not numpy.isclose(dictionary1[key],dictionary2[key],rtol = rtol, atol=atol):
                logging.info("{0}:{1} != {0}:{2}".format(key,dictionary1[key],dictionary2[key]))
                return False
        elif isinstance(dictionary1[key], dict) and isinstance(dictionary2[key], dict):
            # Recursively compare dicts inside this dict
            res = compare_dictionaries(dictionary1[key], dictionary2[key], skipkeys,
                                       rtol, atol)
            if not res:
                return res
        else:
            try:
                if dictionary1[key] != dictionary2[key]:
                    return False
            except:
                logging.error("Error in Comparing {0}:{1} != {0}:{2}".format(key,dictionary1[key],dictionary2[key]))
                return False
    return True

def compare_directories(  directory1, directory2):
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

def compare_pixel_value( imagename=None, refimage=None, loc=None):
    '''
        Compare two images at a certain reference pixel
        @param imagename: Name of the image to be compared to a reference image
        @param refimage: Image to be compared against
        @param loc: The slice or pixel index to compare between the two images
        @return: True if the pixel values match at the provided index or slice. Returns False otherwise
    '''
    if imagename != None and refimage != None:
        if isinstance(loc, str):
            _tb.open(imagename)
            image1 = _tb.getcol('map')
            _tb.close()
            _tb.open(refimage)
            image2 = _tb.getcol('map')
            _tb.close()
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

def compare_pixel_mask( maskname='', refmask=None, refval=None, loc=None):
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
                _tb.open(maskname)
                mask1 = _tb.getcol('PagedArray')
                _tb.close()
                _tb.open(refmask)
                mask2 = _tb.getcol('PagedArray')
                _tb.close()
                return numpy.all(mask1 == mask2)
            else:
                logging.warning('Invalid refmask file name')
        elif refmask == None and refval != None:
            # If using a reference value compare the value/shape to the selected slice
            if isinstance(loc, str):
                _tb.open(maskname)
                image = _tb.getcol('PagedArray')
                _tb.close()
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
