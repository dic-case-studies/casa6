/******************** generated by xml-casa (v2) from imagemetadata.xml *************
********************* a6949a9b048f9c174974cdafb7d2b5c3 *****************************/
%module imagemetadata
%include <casa_typemaps.i>
%feature("kwargs");
%feature("autodoc", "0");

%feature("docstring", "

Summary:
    Add a key-value pair if possible.

Description:


Add a key-value pair if possible.


Input Parameters:
    key                       The name of the FITS or other keyword.
    value                     Associated value to add.

Example:

imd.open('myim.im')
# add a keyword 'test' with value 'first'
if add('test', 'first'):
print 'test=first has been added'
else:
print 'Unable to add key test'
imd.done()

--------------------------------------------------------------------------------
") add;

%feature("docstring", "

Summary:
    Close the image metadata tool. Synonym for done().

Description:



This function closes the image metadata tool.  This means that it detaches the
tool from its underlying metadata object. Methods cannot be run on it until it
is opened with another or the same image.

Example:

imd.open('myim.im')
# do stuff
imd.close()

--------------------------------------------------------------------------------
") close;

%feature("docstring", "

Summary:
    Close the image metadata tool. Synonym for close().

Description:



This function closes the image metadata tool.  This means that it detaches the
tool from its underlying metadata object. Methods cannot be run on it until it
is opened with another or the same image.



Example:

imd.open('myim.im')
# do stuff
imd.done()

--------------------------------------------------------------------------------
") done;

%feature("docstring", "

Summary:
    Get the value associated with the specified, case-insensitive FITS keyword.

Description:


Get the value associated with the specified, case-insensitive FITS keyword.


Input Parameters:
    key                       The name of the FITS or other keyword.

Example:

imd.open('myim.im')
imtype = imd.get('imtype')
imd.done()

--------------------------------------------------------------------------------
") get;

%feature("docstring", "

Summary:
    Get a dictionary of FITS-like header items.

Description:


Get a listing of traditional FITS-like 'header' items.


Input Parameters:
    verbose                   If true, print listing to logger

Example:

imd.open('myim.im')
mylist = imd.list(False)
imd.done()
crval1 = mylist{'crval1'}

--------------------------------------------------------------------------------
") list;

%feature("docstring", "

Summary:
    Open this image metadata tool providing access to an image's metadata.

Description:



This method creates access to the specified image's metadata.



Input Parameters:
    infile                    Image name. The image can be in any casa supported format.

Example:

immd.open('myim.im')
# do stuff with the tool and then close it.
immd.done()

--------------------------------------------------------------------------------
") open;

%feature("docstring", "

Summary:
    Remove or clear the value of a keyword if possible.

Description:


Remove or clear the value of a keyword if possible. If key='masks', a value specifying the mask
to remove may be specified. If no value is specified, all masks are removed.


Input Parameters:
    key                       The name of the FITS or other keyword.
    value                     Value to remove if the key is multi-valued. Only used in the case of key='masks'.

Example:

imd.open('myim.im')
# clear the brightness unit
if imd.remove('bunit'):
print 'bunit has been cleared'
else:
print 'Unable to clear bunit'
imd.done()

--------------------------------------------------------------------------------
") remove;

%feature("docstring", "

Summary:
    Set a keyword to the specified value if possible.

Description:


Set a key-value pair if possible.


Input Parameters:
    key                       The name of the FITS or other keyword.
    value                     Associated value to set.

Example:

Note that when setting the reference value of a polarizaiton axis, one must
provide an array of stokes/polarization strings (['I', 'Q', 'XX']) that is the
same length as the stokes axis. If the stokes axis is degenerate, one can alternatively
provide a string indicating the stokes value.

imd.open('myim.im')
# Set keyword 'telescope' with value 'Argus Array'
if imd.set('telescope', 'Argus Array'):
print 'telescope has been updated'
else:
print 'Unable to update telescope.'
imd.done()

# set polarizations for an image with three pixels on the stokes axis crval3
imd.open('myim.im')
if imd.set('crval3', [XY, LL, 'Q']):
print 'polarization values have been updated'
else:
print 'Unable to update polarization values.'
imd.done()

--------------------------------------------------------------------------------
") set;

%exception {
   try {
      $action
      } catch (const casacore::AipsError &ae) {
         PyErr_SetString(PyExc_RuntimeError, ae.what());
	 //PyErr_Print();
         return NULL;
      }
}
%include "imagemetadata_cmpt.h"

%{
#include <exception>
#include <imagemetadata_cmpt.h>
%}

