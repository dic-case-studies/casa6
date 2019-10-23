/******************** generated by xml-casa (v2) from functional.xml ****************
********************* 53398236e47ef47e219a969efc3c761f *****************************/
%module functional
%include <casa_typemaps.i>
%feature("kwargs");
%feature("autodoc", "0");

%feature("docstring", "

Summary:
    Calculate the value of the functional

Description:



Calculate the value of the functional.

Input Parameters:
    x                         real argument values

Example:

gfn = fn.gaussian1d(2, 0, 1)
#returns 0.125
gfn.f(1)
# returns array([  1.25000000e-01,   3.05175781e-05])
gfn.f([1,2])

--------------------------------------------------------------------------------
") f;

%feature("docstring", "

Summary:
    Get the number of dimensions

Description:


Return the number of dimensions.

Example:

a = fn.gaussian1d()
# nd is set to 1
nd = fn.ndim()

--------------------------------------------------------------------------------
") ndim;

%feature("docstring", "

Summary:
    Free resources of the functional

Description:



Free the functional's resources.

Example:

- a:=dfs.gaussian1d()
- a.setparameters([2,1,3])
T
- a.state()
[type=0, order=-1, ndim=1, npar=3, params=[2 1 3] , masks=[T T T] ]
- a.done()
T
- is_functional(a)
F
- a
F

--------------------------------------------------------------------------------
") done;

%feature("docstring", "

Summary:
    Create and return a new functional tool representing a 1D Gaussian function

Description:



Create a 1-dimensional Gaussian with the specified amplitude, fwhm, and
center.

Input Parameters:
    amplitude                 amplitude of Gaussian
    center                    center of Gaussian
    fwhm                      FWHM of Gaussian

Example:

# get the value and derivatives of a Gaussian with
# height=2; center at x=1; a width of 1 at x=[0,1,2]
gfn = fn.gaussian1d(2,1,1)

# returns array([ 0.125,  2.   ,  0.125])
vals = gfn.f([0, 1, 2])

--------------------------------------------------------------------------------
") gaussian1d;

%feature("docstring", "

Summary:
    Create a 2D Gaussian function

Description:



Create a 2-dimensional Gaussian with the specified amplitude, fwhms, and
center. The created functional has method {em f}  to
calculate the function value at a series of {em x, y} values, or the
value.

Input Parameters:
    amplitude                 Amplitude of Gaussian
    center                    Center (x,y) position. Must have exactly 2 elements.
    fwhm                      FWHM of the axes. Must have exactly 2 elements.
    pa                        The angle between the positive y axis and the major axis, measured counterclockwise.

Example:

# major axis along the y axis
g2d = fn.gaussian2d(1,[0,0],[3,2],'90deg')

# both these commands return 0.5
v = g2d([0, 1])
v = g2d([1.5, 0])

# returns array([ 0.5,  0.5])
v =  g2d.f([0,1,1.5,0])

--------------------------------------------------------------------------------
") gaussian2d;

%feature("docstring", "

Summary:
    Create and return a new functional tool representing a 1D polynomial function, y = c_0 + c_1*x + c_2*x**2 + ... + c_n*x**n

Description:



Create a 1-dimensional polynomial function with the specified coefficents.

Input Parameters:
    coefficients              Array of coefficients. Number of coefficients determines order of polynomial.

Example:

# get the value and derivatives of 3 + 2*x + 4*x*x
poly = fn.powerlogpoly(3, 2, 4)

# value at 3
vals = poly.f(3)

--------------------------------------------------------------------------------
") polynomial;

%feature("docstring", "

Summary:
    Create and return a new functional tool representing a 1D power log polynomial function, y = c_0 * x**( c_1 + c_2*ln(x) + c_3*ln(x)**2 + ... c_n*ln(x)**(n-1)

Description:



Create a 1-dimensional power log polynomial function with the specified coefficents.

Input Parameters:
    coefficients              Array of coefficients.

Example:

# get the value and derivatives of 2*x**(1+ln(x))
plp = fn.powerlogpoly(2,1,1)

# value at 3
vals = plp.f(3)

--------------------------------------------------------------------------------
") powerlogpoly;

%exception {
   try {
      $action
      } catch (const casacore::AipsError &ae) {
         PyErr_SetString(PyExc_RuntimeError, ae.what());
	 //PyErr_Print();
         return NULL;
      }
}
%include "functional_cmpt.h"

%{
#include <exception>
#include <functional_cmpt.h>
%}

