//# tImageConvolver.cc: 
//# Copyright (C) 1996,1997,1999,2000,2001,2002,2003,2004
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU Library General Public License as published by
//# the Free Software Foundation; either version 2 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
//# License for more details.
//#
//# You should have received a copy of the GNU Library General Public License
//# along with this library; if not, write to the Free Software Foundation,
//# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//# $Id: tImageConvolver.cc 20329 2008-06-06 07:59:22Z gervandiepen $
//
#include <casa/aips.h>
#include <casa/Arrays/Vector.h>
#include <casa/OS/EnvVar.h>
#include <casa/Exceptions/Error.h>
#include <casa/Inputs/Input.h>
#include <casa/Logging.h>
#include <tables/Tables/Table.h>
#include <casa/BasicSL/String.h>


#include <images/Images/PagedImage.h>
#include <casa/iostream.h>

using namespace std;
using namespace casacore;
// this header assumes these namespaces are available
#include <imageanalysis/ImageAnalysis/ImageConvolver.h>
using namespace casa;

int main (int argc, const char* argv[])
{

try {

   Input inputs(1);

// Get inputs
    String casapath = EnvironmentVariable::get("CASAPATH");
    if (casapath.empty()) {
        cerr << "CASAPATH env variable not defined. Can't find fixtures. Did you source the casainit.(c)sh file?" << endl;
        return 1;
    }

    String *parts = new String[2];
    split(casapath, parts, 2, String(" "));
    String datadir = parts[0] + "/data/regression/unittest/imageanalysis/ImageAnalysis/";
    delete [] parts;

   String name = datadir + "test_image.im";
   inputs.create("in", name, "Input file name");
   inputs.create("out", "tImageConvolver_tmp", "Output image name");
   inputs.readArguments(argc, argv);

   const String in = inputs.getString("in");
   const String out = inputs.getString("out");
//
   LogOrigin lor("tImageConvolver", "main()", WHERE);
   LogIO os(lor);
//
   if (in.empty()) {
      os << "You must give an input image name" << LogIO::EXCEPTION;
   }
   if (out.empty()) {
      os << "You must give an output image name" << LogIO::EXCEPTION;
   }
//
   DataType imageType = imagePixelType(in);
   if (imageType!=TpFloat) {
      os << "The image must be of type Float" << LogIO::EXCEPTION;
      return 1;
   }

// Construct image

   PagedImage<Float> inImage(in);

// Make kernel

   IPosition shape = inImage.shape();
   Array<Float> kernel(shape);
   kernel = 0.0;
   IPosition pos(shape);
   for (uInt i=0; i<pos.nelements(); i++) {
      pos(i) = shape(i) / 2;
   }
   kernel(pos) = 1.0;

// Convolve

   {
      ImageConvolver<Float> aic;
      PagedImage<Float> outImage(inImage.shape(), inImage.coordinates(), out);

// This function calls the other convolve function

      aic.convolve (os, outImage, inImage, kernel, ImageConvolver<Float>::AUTOSCALE, 
                    1.0, true);
   }


// Test some other things

   {
      ImageConvolver<Float> aic;
      ImageConvolver<Float> aic2(aic);
      aic2 = aic;
   }
//

   Table::deleteTable(out, true);
}

  catch (AipsError x) {
     cerr << "aipserror: error " << x.getMesg() << endl;
     return 1;
  } 

   return 0;
}

