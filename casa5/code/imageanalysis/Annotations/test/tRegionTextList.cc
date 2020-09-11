//# Copyright (C) 2000,2001,2003
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This program is free software; you can redistribute it and/or modify it
//# under the terms of the GNU General Public License as published by the Free
//# Software Foundation; either version 2 of the License, or (at your option)
//# any later version.
//#
//# This program is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
//# more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with this program; if not, write to the Free Software Foundation, Inc.,
//# 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#

#include <casa/aips.h>
#include <imageanalysis/Annotations/RegionTextList.h>
#include <imageanalysis/Annotations/AnnCircle.h>
#include <imageanalysis/Annotations/AnnEllipse.h>
#include <casacore/casa/OS/EnvVar.h>
#include <casacore/casa/BasicMath/Math.h>
#include <casacore/coordinates/Coordinates/CoordinateUtil.h>

#include <casacore/images/Images/ImageUtilities.h>

#include <casacore/casa/namespace.h>

#include <iomanip>

using namespace casa;

int main () {
	LogIO log;
	try {
		CoordinateSystem csys = CoordinateUtil::defaultCoords4D();
		Vector<Double> refVal = csys.referenceValue();
		cout << refVal << endl;
		Vector<String> units = csys.worldAxisUnits();
		units[0] = units[1] = "rad";
		csys.setWorldAxisUnits(units);
		refVal[0] = 4.296556;
		refVal[1] = 0.240673;
		csys.setReferenceValue(refVal);
		AnnotationBase::unitInit();
		String *parts = new String[2];
		split(EnvironmentVariable::get("CASAPATH"), parts, 2, String(" "));
		String goodFile = parts[0]
		    + "/data/regression/unittest/imageanalysis/Annotations/goodAsciiAnnotationsFile.txt";
		String imageName = parts[0]
		    + "/data/regression/unittest/imageanalysis/Annotations/test.im";
		String circleFile = parts[0]
		    + "/data/regression/unittest/imageanalysis/Annotations/circle.crtf";
		delete [] parts;

		RegionTextList list(
			goodFile, csys,
			IPosition(csys.nPixelAxes(), 2000, 2000, 4, 2000)
		);
		cout << std::setprecision(9) << list << endl;
		AlwaysAssert(list.nLines() == 33, AipsError);
		AsciiAnnotationFileLine line31 = list.lineAt(31);
		AlwaysAssert(
			line31.getType() == AsciiAnnotationFileLine::ANNOTATION,
			AipsError
		);
		AlwaysAssert(
			line31.getAnnotationBase()->getLineWidth() == 9,
			AipsError
		);
		AlwaysAssert(
			line31.getAnnotationBase()->isRegion(),
			AipsError
		);
		AlwaysAssert(
			line31.getAnnotationBase()->isAnnotationOnly(),
			AipsError
		);

		AsciiAnnotationFileLine line23 = list.lineAt(23);
		AlwaysAssert(
			line23.getType() == AsciiAnnotationFileLine::ANNOTATION,
			AipsError
		);
		AlwaysAssert(
			line23.getAnnotationBase()->getLineWidth() == 22,
			AipsError
		);
		AlwaysAssert(
			! line23.getAnnotationBase()->isRegion(),
			AipsError
		);
		Quantity k(4.1122334455667778, "km");
		cout << std::setprecision(10) << k << endl;
		{
			ImageInterface<Float> *image;
			ImageUtilities::openImage(image, imageName);
			csys = image->coordinates();
			RegionTextList list2(
				circleFile, csys,
				IPosition(csys.nPixelAxes(), 2000, 2000, 2000)
			);
			Vector<AsciiAnnotationFileLine>  lines2 = list2.getLines();
			for (uInt i=0; i < lines2.size(); ++i) {
				auto ann = lines2[i].getAnnotationBase();
				if (lines2[i].getType() != AsciiAnnotationFileLine::ANNOTATION) {
					continue;
				}
				auto type = ann->getType();
				if (i == 1) {
					AlwaysAssert(type == AnnotationBase::CIRCLE, AipsError);
					const AnnCircle *cir = dynamic_cast<const AnnCircle *>(ann.get());
					AlwaysAssert(cir->getRadius().getValue("arcsec") == 10, AipsError);
				}
				else if (i == 2) {
					cout << lines2[i] << endl;
					cout << AnnotationBase::typeToString(type);
					AlwaysAssert(type == AnnotationBase::ELLIPSE, AipsError);
					const AnnEllipse *ell = dynamic_cast<const AnnEllipse *>(ann.get());
					cout << "major" << ell->getSemiMajorAxis().getValue("arcsec") << endl;
					cout << "minor" << ell->getSemiMinorAxis().getValue("arcsec") << endl;
					cout << "pa " << ell->getPositionAngle() << endl;
					AlwaysAssert(near(ell->getSemiMajorAxis().getValue("arcsec"), 1e-8, 120.0), AipsError);
					AlwaysAssert(near(ell->getSemiMinorAxis().getValue("arcsec"), 1e-8, 80.0), AipsError);
					AlwaysAssert(near(ell->getPositionAngle().getValue("deg"), 0.0), AipsError);
					break;
				}

			}
		}
		{
			log << LogIO::NORMAL
				<< "Test RegionTextList ignores annotation regions outside image" 
				<< LogIO::POST;
			ImageInterface<Float> *image;
			ImageUtilities::openImage(image, imageName);
			csys = image->coordinates();
			RegionTextList list3(
				circleFile, csys,
				IPosition(csys.nPixelAxes(), 10, 10, 2000)
			);
			AlwaysAssert(list3.nLines() == 3, AipsError);
			AlwaysAssert(
				list3.lineAt(0).getType() == AsciiAnnotationFileLine::COMMENT, AipsError
			); // header
			AlwaysAssert(
				list3.lineAt(1).getType() == AsciiAnnotationFileLine::COMMENT, AipsError
			); // blank line
			AlwaysAssert(
				list3.lineAt(2).getType() == AsciiAnnotationFileLine::COMMENT, AipsError
			); // blank line

			log << LogIO::NORMAL
				<< "Test RegionTextList created regions outside image" 
				<< LogIO::POST;
			Bool requireRegion = false;
			RegionTextList list4(
				circleFile, csys,
				IPosition(csys.nPixelAxes(), 10, 10, 2000),
				"", "", "", RegionTextParser::CURRENT_VERSION, true, requireRegion
			);
			AlwaysAssert(list4.nLines() == 5, AipsError);
			AlwaysAssert(
				list4.lineAt(0).getType() == AsciiAnnotationFileLine::COMMENT, AipsError
			); // header
			AlwaysAssert(
				list4.lineAt(1).getType() == AsciiAnnotationFileLine::ANNOTATION, AipsError
			); // circle
			AlwaysAssert(
				list4.lineAt(1).getAnnotationBase()->getType() == AnnotationBase::CIRCLE, AipsError
			); // circle
			AlwaysAssert(
				list4.lineAt(2).getType() == AsciiAnnotationFileLine::ANNOTATION, AipsError
			); // ellipse
			AlwaysAssert(
				list4.lineAt(2).getAnnotationBase()->getType() == AnnotationBase::ELLIPSE, AipsError
			); // ellipse
		}
	}
	catch (const AipsError& x) {
		log << LogIO::SEVERE
			<< "Caught exception: " << x.getMesg()
			<< LogIO::POST;
		return 1;
	}

	cout << "OK" << endl;
	return 0;
}
