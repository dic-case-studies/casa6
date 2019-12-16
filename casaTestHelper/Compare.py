import logging
import sys
import filecmp
import os
import numpy
import six
import numbers
try:
    # CASA 6
    logging.debug("Importing CASAtools")
    import casatools
    tb = casatools.table()
    tb2 = casatools.table()
    casa6 = True

except ImportError:
    # CASA 5
    logging.debug("Import casa6 errors. Trying CASA5...")
    from taskinit import tbtool
    tb = tbtool()
    tb2 = tbtool()
    casa5 = True

class Compare:
    def __init__(self):
        pass


    def compare_CASA_variable_cols(self, referencetab, testtab, varcol, tolerance=0.0):
        '''
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
        tb.open(referencetab)
        tb2.open(testtab)
        col = varcol
        if tb.isvarcol(col) and tb2.isvarcol(col):
            try:
                # First check
                if tb.nrows() != tb2.nrows():
                    logging.error('Length of {} differ from {}, {} != {}'.format(referencetab,testtab,tb.nrows(),tb2.nrows()))
                    retval = False
                else:
                    for therow in range(tb.nrows()):
                        rdata = tb.getcell(col,therow)
                        tdata = tb2.getcell(col,therow)
                        if not rdata.all()==tdata.all():
                            if (tolerance>0.0):
                                differs=False
                                for j in range(0,len(rdata)):
                                    if isinstance(rdata[j], (float, int)):
                                        if (abs(rdata[j]-tdata[j]) > tolerance*abs(rdata[j]+tdata[j])):
    #                                        print('Column ', col,' differs in tables ', referencetab, ' and ', testtab)
    #                                        print(therow, j)
    #                                        print(rdata[j])
    #                                        print(tdata[j])
                                            differs = True
                                    elif isinstance(rdata[j], (list, numpy.ndarray)):
                                        for k in range(0,len(rdata[j])):
                                            if (abs(rdata[j][k]-tdata[j][k]) > tolerance*abs(rdata[j][k]+tdata[j][k])):
    #                                            print('Column ', col,' differs in tables ', referencetab, ' and ', testtab)
    #                                            print(therow, j, k)
    #                                            print(rdata[j][k])
    #                                            print(tdata[j][k])
                                                differs = True
                                    if differs:
                                        print('ERROR: Column {} of {} and {} do not agree within tolerance {}'.format(col,referencetab, testtab, tolerance))
                                        retval = False
                                        break
                            else:
                                print('ERROR: Column {} of {} and {} do not agree.'.format(col,referencetab, testtab))
                                print('ERROR: First row to differ is row={}'.format(therow))
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

    def compare_CASA_tables(self, referencetab, testtab, excludecols = None, tolerance=0.001, mode="percentage", startrow = 0, nrow = -1, rowincr = 1):
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

        if excludecols is None:
            excludecols = []
        if not isinstance(excludecols, list):
            logging.error("excludecols not in correct format")
            raise TypeError("excludecols must be a list")
        if referencetab.endswith("cal") or testtab.endswith("cal"):
            logging.warning("WARNING: Will compare caltables using compare_caltables")
            return Compare().compare_caltables(referencetab, testtab, cols= excludecols, rtol=8e-7, atol=1e-8)
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
                                if ((mode == "percentage") and (abs(a[i]-b[i]) > tolerance*abs(a[i]))) or ((mode == "absolute") and (abs(a[i]-b[i]) > tolerance)):
                                    print("Column " + cname + " differs")
                                    print("Row=" + str(i))
                                    print("Reference file value: " + str(a[i]))
                                    print("Input file value: " + str(b[i]))
                                    if (mode == "percentage"):
                                        print("Tolerance is {0}%; observed difference was {1} %".format (tolerance * 100, 100*abs(a[i]-b[i])/abs(a[i])))
                                    else:
                                        print("Absolute tolerance is {0}; observed difference: {1}".format (tolerance, (abs(a[i]-b[i]))))
                                    differs = True
                                    rval = False
                                    break
                            elif isinstance(a[i], (int, numpy.int32)):
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
                            elif isinstance(a[i], (numpy.bool_, str)):
                                if not (a[i]==b[i]):
                                    #print("Column " + c + " differs")
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
                            elif isinstance(a[i], (list, numpy.ndarray)):
                                for j in range(0,len(a[i])):
                                    if differs: break
                                    if isinstance(a[i][j], (list, numpy.ndarray)):
                                        if ((mode=="percentage") and (abs(a[i][j]-b[i][j]) > tolerance*abs(a[i][j]))) or ((mode=="absolute") and (abs(a[i][j]-b[i][j]) > tolerance)):
                                            #print("Column " + c + " differs")
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
                                    elif isinstance(a[i][j], (list, numpy.ndarray)):
                                        it = range(0,len(a[i][j]))
                                        if mode == "percentage":
                                            diff = numpy.abs(numpy.subtract(a[i][j], b[i][j])) > tolerance * numpy.abs(a[i][j])
                                            it = numpy.where(diff)[0]
                                        elif (mode == "absolute"):
                                            diff = numpy.abs(numpy.subtract(a[i][j], b[i][j])) > tolerance
                                            it = numpy.where(diff)[0]
                                        for k in it:
                                            if differs: break
                                            if ( ((mode == "percentage") and (abs(a[i][j][k]-b[i][j][k]) > tolerance*abs(a[i][j][k]))) \
                                                     or ((mode == "absolute") and (abs(a[i][j][k]-b[i][j][k]) > tolerance)) \
                                                     #or ((mode == "phaseabsdeg") and (phasediffabsdeg(a[i][j][k],b[i][j][k])>tolerance)) \
                                                     ):
                                                #print("Column " + c + " differs")
                                                print("(Row,Channel,Corr)=(" + str(k) + "," + str(j) + "," + str(i) + ")")
                                                print("Reference file value: " + str(a[i][j][k]))
                                                print("Input file value: " + str(b[i][j][k]))
                                                if (mode=="percentage"):
                                                    print("Tolerance in % should be " + str(100*abs(a[i][j][k]-b[i][j][k])/abs(a[i][j][k])))
                                                elif (mode=="absolute"):
                                                    print("Absolute tolerance should be " + str(abs(a[i][j][k]-b[i][j][k])))
                                                elif (mode=="phaseabsdeg"):
                                                    #print("Phase tolerance in degrees should be " + str(phasediffabsdeg(a[i][j][k],b[i][j][k])))
                                                    print("N/A")
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

    def compare_files(self, file1, file2, shallow=False):
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

    def compare_caltables(self, table1, table2, cols=None, rtol=8e-7, atol=1e-8):
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

    def compare_dictionaries(self,  dictionary1, dictionary2, skipkeys = None, rtol=8e-7, atol=1e-8):
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
            elif isinstance(dictionary1[key], six.string_types) and isinstance(dictionary2[key], six.string_types):
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
            else:
                try:
                    if dictionary1[key] != dictionary2[key]:
                        return False
                except:
                    logging.error("Error in Comparing {0}:{1} != {0}:{2}".format(key,dictionary1[key],dictionary2[key]))
                    return False
        return True

    def compare_directories(self,  directory1, directory2):
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
            if not Compare().compare_directories(new_directory1, new_directory2):
                return False
        return True

    def compare_pixel_value(self, imagename=None, refimage=None, loc=None):
        '''
            Compare two images at a certain reference pixel
            @param imagename: Name of the image to be compared to a reference image
            @param refimage: Image to be compared against
            @param loc: The slice or pixel index to compare between the two images
            @return: True if the pixel values match at the provided index or slice. Returns False otherwise
        '''
        if imagename != None and refimage != None:
            if isinstance(loc, six.string_types):
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

    def compare_pixel_mask(self, maskname='', refmask=None, refval=None, loc=None):
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
                    return numpy.all(mask1 == mask2)
                else:
                    logging.warning('Invalid refmask file name')
            elif refmask == None and refval != None:
                # If using a reference value compare the value/shape to the selected slice
                if isinstance(loc, six.string_types):
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
