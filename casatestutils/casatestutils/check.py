import os
import numpy

import casatools
tb = casatools.table()

class Check(object):
    def __init__(self):
        self.val = False
        self.plotSize = ''

    def check_plotfile(self, plotfileName, min_size, max_size=None):
        '''
        Check if plotfile generated is cprrect size
            plotfileName --> Name of plotted Image
            min_size -- > Min Size of image
            max_size --> Max Size of image
            Return : True if image size > min_size ( and < max_size if max_size is provided )
        '''

        if os.path.isfile(plotfileName):
            self.plotSize = os.path.getsize(plotfileName) # Return the size, in bytes, of path.
            print('{} file size is: {}'.format(plotfileName, self.plotSize))
            if self.plotSize > min_size:
                self.val = True
            if max_size is not None:
                if not self.plotSize < max_size:
                    self.val = False

        return self.val

    def check_pixels(self, imagename='', loc=None, refval=None, rtol=1e-05, atol=1e-08):
        '''
            Check pixels in an image to a specified reference value
            @param imagename: input image file
            @param loc: The index of the image to compare to the refval
            @param refval: The reference value to compare the selected pixel(s) to
            @param rtol: The relative tolerance used in the numpy.isclose function
            @param atol: The absolute tolerance used in the numpy.isclose function
            @return: True if the shape and value of the refval and selected pixel match.
        '''
        if not isinstance(loc, str):
            raise TypeError('Please give target location in string list format "20,30,2:4"')
        if os.path.exists(imagename):
            tb.open(imagename)
            image = tb.getcol('map')
            tb.close()
            if type(refval) is not type(None):
                index = []
                to_slice = loc.split(',')
                for item in to_slice:
                    if ':' not in item:
                        index.append(int(item))
                    else:
                        item_split = item.split(':')
                        index.append(slice(int(item_split[0]), int(item_split[1])))
                selected_slice = image[tuple(index)]
                if numpy.shape(selected_slice) != numpy.shape(refval):
                    raise Exception('Please check that the shape of the reference and selected slice are the same')
                isequal = numpy.isclose(selected_slice, refval, rtol=rtol, atol=atol)
                print("For pixel value check the obtained value was {}. The expected value was {} with a tolerance of {}. test success = {}.".format(selected_slice, refval, atol, numpy.all(isequal == True)))
                return numpy.all(isequal == True)
            else:
                raise Exception('Please provide a refernce value to compare against')
        else:
            raise Exception('{} Does Not Exist'.format(imagename))
