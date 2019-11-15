//# DisplayShape.h: Abstract base class for all shapes/annotations objects
//# Copyright (C) 1998,1999,2000,2001,2002
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
//# $Id:

#ifndef TRIALDISPLAY_DISPLAYSHAPE_H
#define TRIALDISPLAY_DISPLAYSHAPE_H

#include <casa/aips.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Containers/Record.h>

namespace casa { //# NAMESPACE CASA - BEGIN

	class DSShape;
	class DSClosed;
	class DisplayEnums;
	class PixelCanvas;
	class DParameterColorChoice;


// <summary>
// The abstract base class for all "DisplayShapes".
// </summary>
//
// <prerequisite>
// </prerequisite>
//
// <etymology>
// DisplayShape is a way of providing a consistant interface to a large number
// of different shapes.
// </etymology>
//
// <synopsis>
// DisplayShape provides a framework from which a large number of different
// shapes can be made, all with the same interface. Any new DisplayShape
// should inherit from this class, or higher level classes
// (see <linkto class="DSPoly">DSPoly</linkto> etc).
//
// There are generally two ways to make DisplayShape(s); To create them in
// "one hit" by providing arguments to the constructor, or by using the
// default constructor and then the "setOptions" method. A simple interface
// for all classes inheriting from the
// <linkto class="DisplayShape">DisplayShape</linkto> class is provided
// by <linkto class="DisplayShapeInterface">DisplayShapeInterface</linkto>.
// </synopsis>
//
// <motivation>
// A common interface to a large number of shapes was desired.
// </motivation>
//
// <example>
// </example>

	class DisplayShape {

	public:
		// Handle style
		enum HandleShape {Filled_Square=0, Open_Square, Filled_Circle,
		                  Open_Circle, Filled_Triangle, Open_Triangle
		                 };


		// Default constructor. Creates shape with default options set
		DisplayShape();
		// Copy constructor
		DisplayShape(const DisplayShape& other);
		// Destructor
		virtual ~DisplayShape();

		// These functions contol behaviour of handles during a call
		// to the display shape object.
		// (These calls should be propogated up through the class tree).
		// <group>
		virtual void draw(PixelCanvas* pc);
		virtual void rotateAbout(const casacore::Float& relAngle, const casacore::Float& aboutX,
		                         const casacore::Float& aboutY) ;
		virtual void move(const casacore::Float& dX, const casacore::Float& dY);
		// </group>


		// Rotate the supplied polygon (column 1 - x values, column 2 - y values)
		// about the supplied point by the supplied angle. NB Angle in radians
		virtual casacore::Matrix<casacore::Float> rotatePolygon(const casacore::Matrix<casacore::Float>& toRotate,
		                                    const casacore::Float& angle,
		                                    const casacore::Float& aboutX,
		                                    const casacore::Float& aboutY);

		// Rotates a point around the point specified. NB Angle in radians.
		virtual casacore::Vector<casacore::Float> rotatePoint(const casacore::Vector<casacore::Float>& toRotate,
		                                  const casacore::Float& angle,
		                                  const casacore::Float& aboutX, const casacore::Float& aboutY);

		// Translate an entire matrix by the specified dx / dy amounts.
		virtual casacore::Matrix<casacore::Float> translateMatrix(const casacore::Matrix<casacore::Float>& points,
		                                      const casacore::Float& dx, const casacore::Float& dy);

		// Is xPos, YPos inside the supplied points (column 1 - x values,
		// clolumn 2 - y values)
		virtual casacore::Bool inPolygon(const casacore::Matrix<casacore::Float>& points, const casacore::Float& xPos,
		                       const casacore::Float& yPos);

		// Determine the two vertices (firstVert, secondVert) which join the line
		// closest to the xPos, yPos point supplied. If closedPoly is left as
		// true, the points supplied are treated as a polygon, if not as a poly
		// line.
		virtual casacore::Bool closestLine(const casacore::Matrix<casacore::Float>& points, const casacore::Float& xPos,
		                         const casacore::Float& yPos,
		                         casacore::Int& firstVert, casacore::Int& secondVert,
		                         const casacore::Bool& closedPoly = true);

		// For a specified set of points, find the closest to xPos,YPos. out
		// relates the matrix index (row number) of the closest point.
		virtual casacore::Bool closestPoint(const casacore::Matrix<casacore::Float>& points,
		                          const casacore::Float& xPos, const casacore::Float& yPos,
		                          casacore::Int& out);

		// Find the closest two Points from a casacore::Matrix to the specified point.
		virtual casacore::Bool closestPoints(const casacore::Matrix<casacore::Float>& points,
		                           const casacore::Float& xPos, const casacore::Float& yPos,
		                           casacore::Int& outClosest, casacore::Int& outSecond);

		// Is the supplied point within the DisplayShape?
		virtual casacore::Bool inObject(const casacore::Float& xPos, const casacore::Float& yPos) = 0;

		// Convert degrees to radians
		virtual casacore::Float toRadians(const casacore::Float& degrees);

		// Conver radians to degree
		virtual casacore::Float toDegrees(const casacore::Float& radians);

		// Sets the center of the DisplayShape
		virtual void setCenter(const casacore::Float& xPos, const casacore::Float& yPos) = 0;

		// Returns the center of the DisplayShape (x,y).
		virtual casacore::Vector<casacore::Float> getCenter() = 0;

		// Changes the closest point to the supplied location to that location
		virtual void changePoint(const casacore::Vector<casacore::Float>& newPos) = 0;

		// Changes the nth point making up the DisplayShape ot the specified
		// location.
		virtual void changePoint(const casacore::Vector<casacore::Float>& newPoint,
		                         const casacore::Int nPoint) = 0;

		// If applicable, this function will add a point to the shape in the
		// most meaningful location.
		virtual void addPoint(const casacore::Vector<casacore::Float>& /*newPoint*/) { };

		// Rotate the shape about its center by a set angle (angle in degrees).
		virtual void rotate(const casacore::Float& angle) = 0;

		// Scale the shape about its center by the scaleFactor
		virtual void scale(const casacore::Float& scaleFactor) = 0;

		// Allow locking of other shapes onto this one. When a shape is locked,
		// if the current shape is moved, so to will the locked shape.
		virtual void addLocked(DisplayShape* toLock);

		// Removes a lock from the specified shape.
		virtual void removeLocked(DisplayShape* removeLock);

		// Handle management.
		// <group>
		virtual void buildHandles(const casacore::Matrix<casacore::Float>& startPoints);
		virtual casacore::Matrix<casacore::Float> getHandleLocations();
		virtual void setHandlePositions(const casacore::Matrix<casacore::Float>& newPoints);
		virtual DSClosed*  makeHandle(const casacore::Vector<casacore::Float>& newHandlePos);
		virtual void addHandle(const casacore::Vector<casacore::Float>& newHandlePos,
		                       const casacore::Bool& atEnd = true
		                               , const casacore::Int position = 0);
		virtual casacore::Bool removeHandle(const casacore::Vector<casacore::Float>& getRidOf);
		virtual casacore::Bool removeHandle(const casacore::Int nHandle);

		virtual casacore::Bool onHandles(const casacore::Float& xPos, const casacore::Float& yPos);
		virtual casacore::Bool whichHandle(const casacore::Float& xPos, const casacore::Float& yPos, casacore::Int& out);

		virtual void setDrawHandles(const casacore::Bool& shouldIDraw);
		virtual casacore::Bool drawingHandles() {
			return itsDrawHandles;
		}
		virtual void setHasHandles(const casacore::Bool& hasHandles);
		virtual void setHandleShape(const DisplayShape::HandleShape& shape);
		virtual void setHandleSize(const casacore::Int pixelSize);
		virtual void setHandleColor(const casacore::String& handleColor);
		virtual casacore::uInt nHandles();
		// </group>

		// Manage the color of object. (Does not include handles)
		// <group>
		virtual void setColor(const casacore::String& newColor);
		virtual casacore::String getColor();
		// </group>

		// Settings
		// <group>
		virtual casacore::Record getOptions();
		virtual casacore::Bool setOptions(const casacore::Record& settings);
		// </group>

		virtual void recalculateScreenPosition() {}

	private:
		// Set default options
		virtual void setDefaultOptions();


		// Object
		DParameterColorChoice* itsColor;

		// Handles
		casacore::PtrBlock<DSClosed*> itsHandles;

		// Locks
		casacore::PtrBlock<DisplayShape*> itsLocks;

		// Do I have handles / can a user resize me?
		// i.e. Do I *ever* want to draw handles (e.g. will be false for an object
		//      which IS a handle!)
		casacore::Bool itsHasHandles;

		// Should handles be shown if they exist
		casacore::Bool itsDrawHandles;

		// Have valid handles been made/supplied yet?
		casacore::Bool itsValidHandles;

		// Handle settings
		// <group>
		casacore::String itsHandleColor;
		DisplayShape::HandleShape itsHandleShape;
		casacore::Int itsHandleSize;
		// </group>

	};

} //# NAMESPACE CASA - END

#endif




