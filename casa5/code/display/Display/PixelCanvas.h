//# PixelCanvas.h: Base class defining interface to PixelCanvases
//# Copyright (C) 1995,1996,1997,1998,1999,2000,2001,2002,2003
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
//# $Id$

#ifndef DISPLAY_PIXELCANVAS_H
#define DISPLAY_PIXELCANVAS_H

#include <casa/aips.h>
#include <list>
#include <display/Display/DisplayEnums.h>
#include <display/Display/PixelCanvasColorTable.h>
#include <display/Display/PCVGBuffer.h>
#include <display/Display/DLFont.h>

namespace casacore{

	template <class T> class Vector;
	template <class T> class Matrix;
}

namespace casa { //# NAMESPACE CASA - BEGIN

	class Colormap;
	class PCMotionEH;
	class PCPositionEH;
	class PCRefreshEH;


// <summary>
// Base class defining interface to pixel-based output devices.
// </summary>
//
// <prerequisite>
// <li> <linkto class="Colormap">Colormap</linkto>
// <li> <linkto class="PixelCanvasColorTable">PixelCanvasColorTable</linkto>
// <li> <linkto class="PCMotionEH">PCMotionEH</linkto>
// <li> <linkto class="PCRefreshEH">PCRefreshEH</linkto>
// <li> <linkto class="PCPositionEH">PCPositionEH</linkto>
// </prerequisite>
//
// <etymology>
// PixelCanvas is the mechanism for drawing on the screen.
// </etymology>
//
// <synopsis>
// The philosophy of the PixelCanvas is to provide flexible fast interface to
// an underlying graphics system in terms of integer pixel positions and
// color values.  The interface should be as simple as possible and not demand
// complicated structures on its interface.
//
// The PixelCanvas performs minimal management, leaving up to the derived classes
// for most of the work required to interface to the underlying graphics library.
///
//
// To make it flexible, the fundamental interface accepts pointers to arrays of all
// the scalar AIPS++ types.  The casacore::Bool type is not acceptable at the PixelCanvas
// level because it cannot be used to represent a color index, as are the two
// complex number types.
//
// To make it fast, a caching mechanism is used to allow display lists to be
// created in the format native to the underlying graphics library.
// The caching works like OpenGL display lists.
//
// To create a display list:
//
// <ol>
// <li> Call the function newList() which will return a list id.  You need to
//      store the returned id somewhere so that it may be recalled later.
// <li> Perform some drawing commands.  These commands are output-only commands
//      that change the state of the PixelCanvas or draw some graphics or other
//      function that affects the canvas.
// <li> Call the function endList()
// </ol>
//
// To recall the drawing commands:
//
// <ol>
// <li> Call drawList(), passing the list id as the parameter
// </ol>
//
// To delete the list, call deleteList(), with the id
//
// The PixelCanvas maintains a translation stack which may be driven by
// calls to translate and calls to pushMatrix, popMatrix.
//
// The translation stack is an effective way to draw the same graphic in different
// places:
// <srcblock>
//
// casacore::uInt myGraphic = newList();
// ...
// endList();
// casacore::Matrix m(n,2);
// for (casacore::uInt i = 0; i < n; i++)
//   {
//     pc->pushMatrix();
//     pc->translate(m(i,0), m(i,1));
//     pc->drawList(myGraphic);
//     pc->popMatrix();
//   }
//
// </srcblock>
//
// Images are most correctly drawn through the following sequence of operations
//
// <ol>
// <li> Obtain the size of the current colormap with getColormapSize()
// <li> scale your data to fit in the range of [0,size-1];
// <li> call mapToColor() to get a proper color image (which you may consider
//      saving).
// <li> call drawImage()
// </ol>
//
// You may find that a class derived from
// <linkto class="WCDataScaleHandler">WCDataScaleHandler</linkto>, such as
// <linkto class="WCLinearScaleHandler">WCLinearScaleHandler</linkto>
// may be useful in step #2 above.
//
// mapToColor is also useful for transforming values that are associated with
// vector graphics as well (e.g., contour lines).
//
//
// The PixelCanvas layer is quite thin.  Most functionality is implemented
// in the derived classes.
//
// </synopsis>
//
// <motivation>
// Want a generic interface to possibly a variety of graphics hardware types.
// Want base class to maintain callback lists
// </motivation>
//
// <example>
// see the Display test directory
// </example>
//

	class PixelCanvas {
	public:
		virtual ~PixelCanvas();

		// add event handlers
		// <group>
		void addRefreshEventHandler(const PCRefreshEH &eh);
		void addMotionEventHandler(const PCMotionEH &eh);
		void addPositionEventHandler(const PCPositionEH &eh);
		// </group>

		// remove event handlers
		// <group>
		void removeRefreshEventHandler(const PCRefreshEH &eh);
		void removeMotionEventHandler(const PCMotionEH &eh);
		void removePositionEventHandler(const PCPositionEH &eh);
		// </group>

		// call event handlers
		// <group>
		void callRefreshEventHandlers(Display::RefreshReason reason);
		void callMotionEventHandlers(casacore::Int x, casacore::Int y, casacore::uInt state);
		void callPositionEventHandlers(Display::KeySym keysym, casacore::Bool keystate,
		                               casacore::Int x, casacore::Int y, casacore::uInt state);
		// </group>

		// enabling/disabling of event tracking
		// <group>
		virtual void enableMotionEvents() = 0;
		virtual void disableMotionEvents() = 0;
		virtual void enablePositionEvents() = 0;
		virtual void disablePositionEvents() = 0;
		// </group>

		// Does this canvas support cached display lists?  The user of the
		// canvas should always check this, because undefined behaviour can
		// result when an attempt is made to use a list on a PixelCanvas
		// which does not support lists.
		virtual casacore::Bool supportsLists() = 0;

		// begin caching display commands - return list ID
		virtual casacore::uInt newList() = 0;
		// end caching display commands
		virtual void endList() = 0;
		// (Cacheable) recall cached display commands
		virtual void drawList(casacore::uInt list) = 0;
		// translate all lists
		virtual void translateAllLists(casacore::Int xt, casacore::Int yt) = 0;
		// translate the list
		virtual void translateList(casacore::uInt list, casacore::Int xt, casacore::Int yt) = 0;
		// remove list from cache
		virtual void deleteList(casacore::uInt list) = 0;
		// flush all lists from the cache
		virtual void deleteLists() = 0;
		// return true if the list exists
		virtual casacore::Bool validList(casacore::uInt list) = 0;

		// (Cacheable) Set the font to the recognizable font name
		virtual casacore::Bool setFont(const casacore::String &fontName) = 0;

		// TODO: These should become abstract
		// Set the font via the DisplayLibrary Font class
		virtual casacore::Bool setFont(DLFont* /*font*/) {
			return false;
		}

		// Set the font to font name / size
		virtual casacore::Bool setFont(const casacore::String& /*fontName*/, const casacore::Int /*fontSize*/) {
			return false;
		}

		// (Cacheable) Draw text using that font aligned in some way to the
		// position
		virtual void drawText(casacore::Int x, casacore::Int y, const casacore::String &text,
		                      Display::TextAlign alignment = Display::AlignCenter) = 0;

		// TODO This should become abstract - NYI in GLPixelCanvas currently
		// Draw text at a specified angle.
		virtual void drawText(casacore::Int /*x*/, casacore::Int /*y*/, const casacore::String &/*text*/,
		                      const casacore::Float& /*angle*/,
		                      Display::TextAlign /*alignment*/ = Display::AlignCenter)
		{ }

		// TODO : This should become abstract
		// Determine the width / height of a string of text based on
		// current settings.
		// <group>
		virtual casacore::Int textWidth(const casacore::String& /*text*/) {
			return -1;
		}
		virtual casacore::Int textHeight(const casacore::String& /*text*/) {
			return -1;
		}
		// </group>

		// (Cacheable) Draw an array of 2D color data as a raster image for zoom = <1,1>
		// <group>
		virtual void drawImage(const casacore::Matrix<casacore::uInt> &data, casacore::Int x, casacore::Int y) = 0;
		virtual void drawImage(const casacore::Matrix<casacore::Int> &data, casacore::Int x, casacore::Int y) = 0;
		virtual void drawImage(const casacore::Matrix<casacore::uLong> &data, casacore::Int x, casacore::Int y) = 0;
		virtual void drawImage(const casacore::Matrix<casacore::Float> &data, casacore::Int x, casacore::Int y) = 0;
		virtual void drawImage(const casacore::Matrix<casacore::Double> &data, casacore::Int x, casacore::Int y) = 0;
		// </group>

		// (Cacheable) Draw an array of 2D color data as a raster image,
		// taking note of the <src>casacore::Bool</src> mask.
		// Set opaqueMask to true to draw masked pixels in the background color;
		// otherwise they will be transparent (letting whatever was drawn
		// previously at that point show through).
		// <group>
		virtual void drawImage(const casacore::Int &/*x*/, const casacore::Int &/*y*/,
		                       const casacore::Matrix<casacore::uInt> &/*data*/,
		                       const casacore::Matrix<casacore::Bool> &/*mask*/,
		                       casacore::Bool /*opaqueMask*/=false) {
			return;
		}
		// </group>

		// (Cacheable) Draw an array of 2D color data as a raster image for any positive integer zoom
		// <group>
		virtual void drawImage(const casacore::Matrix<casacore::uInt> &data, casacore::Int x, casacore::Int y,
		                       casacore::uInt xzoom, casacore::uInt yzoom) = 0;
		virtual void drawImage(const casacore::Matrix<casacore::Int> &data, casacore::Int x, casacore::Int y,
		                       casacore::uInt xzoom, casacore::uInt yzoom) = 0;
		virtual void drawImage(const casacore::Matrix<casacore::uLong> &data, casacore::Int x, casacore::Int y,
		                       casacore::uInt xzoom, casacore::uInt yzoom) = 0;
		virtual void drawImage(const casacore::Matrix<casacore::Float> &data, casacore::Int x, casacore::Int y,
		                       casacore::uInt xzoom, casacore::uInt yzoom) = 0;
		virtual void drawImage(const casacore::Matrix<casacore::Double> &data, casacore::Int x, casacore::Int y,
		                       casacore::uInt xzoom, casacore::uInt yzoom) = 0;
		// </group>

		// (Cacheable) Draw a component of a multi-channel image, storing it
		// in buffers until flushComponentImages() is called.
		virtual void drawImage(const casacore::Matrix<casacore::uInt> &data, const casacore::Int &x, const casacore::Int &y,
		                       const Display::ColorComponent &colorcomponent) = 0;

		// Fill one of the channel buffers.
		virtual void bufferComponent(const casacore::Matrix<casacore::uInt> &data,
		                             const casacore::Int &x, const casacore::Int &y,
		                             const Display::ColorComponent
		                             &colorcomponent) = 0;

		// (NOT CACHEABLE!) Flush the component buffers.
		virtual void flushComponentBuffers() = 0;

		// (Cacheable) Draw a single point using current color
		// <group>
		virtual void drawPoint(casacore::Int x1, casacore::Int y1) = 0;
		virtual void drawPoint(casacore::Float x1, casacore::Float y1) = 0;
		virtual void drawPoint(casacore::Double x1, casacore::Double y1) = 0;
		// </group>

		// (Cacheable) Draw N points specified as a Nx2 matrix
		// <group>
		virtual void drawPoints(const casacore::Matrix<casacore::Int> &verts) = 0;
		virtual void drawPoints(const casacore::Matrix<casacore::Float> &verts) = 0;
		virtual void drawPoints(const casacore::Matrix<casacore::Double> &verts) = 0;
		// </group>

		// (Cacheable) Draw a bunch of points using current color
		// <group>
		virtual void drawPoints(const casacore::Vector<casacore::Int> &x1, const casacore::Vector<casacore::Int> &y1) = 0;
		virtual void drawPoints(const casacore::Vector<casacore::Float> &x1,
		                        const casacore::Vector<casacore::Float> &y1) = 0;
		virtual void drawPoints(const casacore::Vector<casacore::Double> &x1,
		                        const casacore::Vector<casacore::Double> &y1) = 0;
		// </group>

		// (Cacheable) Draw a single line using current color
		// <group>
		virtual void drawLine(casacore::Int x1, casacore::Int y1, casacore::Int x2, casacore::Int y2) = 0;
		virtual void drawLine(casacore::Float x1, casacore::Float y1, casacore::Float x2, casacore::Float y2) = 0;
		virtual void drawLine(casacore::Double x1, casacore::Double y1, casacore::Double x2, casacore::Double y2) = 0;
		// </group>

		// (Cacheable) Draw N/2 lines from an Nx2 matrix
		// <group>
		virtual void drawLines(const casacore::Matrix<casacore::Int> &verts) = 0;
		virtual void drawLines(const casacore::Matrix<casacore::Float> &verts) = 0;
		virtual void drawLines(const casacore::Matrix<casacore::Double> &verts) = 0;
		// </group>

		// (Cacheable) Draw a bunch of unrelated lines using current color
		// <group>
		virtual void drawLines(const casacore::Vector<casacore::Int> &x1, const casacore::Vector<casacore::Int> &y1,
		                       const casacore::Vector<casacore::Int> &x2, const casacore::Vector<casacore::Int> &y2) = 0;
		virtual void drawLines(const casacore::Vector<casacore::Float> &x1, const casacore::Vector<casacore::Float> &y1,
		                       const casacore::Vector<casacore::Float> &x2,
		                       const casacore::Vector<casacore::Float> &y2) = 0;
		virtual void drawLines(const casacore::Vector<casacore::Double> &x1, const casacore::Vector<casacore::Double> &y1,
		                       const casacore::Vector<casacore::Double> &x2,
		                       const casacore::Vector<casacore::Double> &y2) = 0;
		// </group>

		// (Cacheable) Draw a single connected line between the points given
		// <group>
		virtual void drawPolyline(const casacore::Vector<casacore::Int> &x1,
		                          const casacore::Vector<casacore::Int> &y1) = 0;
		virtual void drawPolyline(const casacore::Vector<casacore::Float> &x1,
		                          const casacore::Vector<casacore::Float> &y1) = 0;
		virtual void drawPolyline(const casacore::Vector<casacore::Double> &x1,
		                          const casacore::Vector<casacore::Double> &y1) = 0;
		// </group>

		// (Cacheable) Draw N-1 connected lines from Nx2 matrix of vertices
		// <group>
		virtual void drawPolyline(const casacore::Matrix<casacore::Int> &verts) = 0;
		virtual void drawPolyline(const casacore::Matrix<casacore::Float> &verts) = 0;
		virtual void drawPolyline(const casacore::Matrix<casacore::Double> &verts) = 0;
		// </group>

		// Draw a "marker". See <linkto class="Display">Display</linkto>
		// for a list of available markers.
		// <group>
		virtual void drawMarker(const casacore::Int& x1, const casacore::Int& y1,
		                        const Display::Marker& marker,
		                        const casacore::Int& pixelHeight);
		virtual void drawMarker(const casacore::Float& x1, const casacore::Float& y1,
		                        const Display::Marker& marker,
		                        const casacore::Int& pixelHeight);
		virtual void drawMarker(const casacore::Double& x1, const casacore::Double& y1,
		                        const Display::Marker& marker,
		                        const casacore::Int& pixelHeight);
		// </group>

		// (Cacheable) Draw a closed polygon
		// <group>
		virtual void drawPolygon(const casacore::Vector<casacore::Int> &x1, const casacore::Vector<casacore::Int> &y1) = 0;
		virtual void drawPolygon(const casacore::Vector<casacore::Float> &x1,
		                         const casacore::Vector<casacore::Float> &y1) = 0;
		virtual void drawPolygon(const casacore::Vector<casacore::Double> &x1,
		                         const casacore::Vector<casacore::Double> &y1) = 0;
		// </group>

		// (Cacheable) Draw and fill a closed polygon
		// <group>
		virtual void drawFilledPolygon(const casacore::Vector<casacore::Int> &x1,
		                               const casacore::Vector<casacore::Int> &y1) = 0;
		virtual void drawFilledPolygon(const casacore::Vector<casacore::Float> &x1,
		                               const casacore::Vector<casacore::Float> &y1) = 0;
		virtual void drawFilledPolygon(const casacore::Vector<casacore::Double> &x1,
		                               const casacore::Vector<casacore::Double> &y1) = 0;
		// </group>

		// (Cacheable) Draw a closed N-sided polygon from Nx2 matrix of vertices
		// <group>
		virtual void drawPolygon(const casacore::Matrix<casacore::Int> &verts) = 0;
		virtual void drawPolygon(const casacore::Matrix<casacore::Float> &verts) = 0;
		virtual void drawPolygon(const casacore::Matrix<casacore::Double> &verts) = 0;
		// </group>

		// (Cacheable) Draw a rectangle
		// <group>
		virtual void drawRectangle(casacore::Int x1, casacore::Int y1, casacore::Int x2, casacore::Int y2) = 0;
		virtual void drawRectangle(casacore::Float x1, casacore::Float y1, casacore::Float x2, casacore::Float y2) = 0;
		virtual void drawRectangle(casacore::Double x1, casacore::Double y1, casacore::Double x2, casacore::Double y2) = 0;
		// </group>

		// (Cacheable) Draw a filled rectangle
		// <group>
		virtual void drawFilledRectangle(casacore::Int x1, casacore::Int y1, casacore::Int x2, casacore::Int y2) = 0;
		virtual void drawFilledRectangle(casacore::Float x1, casacore::Float y1, casacore::Float x2, casacore::Float y2) = 0;
		virtual void drawFilledRectangle(casacore::Double x1, casacore::Double y1, casacore::Double x2,
		                                 casacore::Double y2) = 0;
		// </group>

		// (Cacheable) Draw a set of points, specifying a color per point to be drawn.
		// <group>
		virtual void drawColoredPoints(const casacore::Vector<casacore::Int> &x1,
		                               const casacore::Vector<casacore::Int> &y1,
		                               const casacore::Vector<casacore::uInt> &colors) = 0;
		virtual void drawColoredPoints(const casacore::Vector<casacore::Float> &x1,
		                               const casacore::Vector<casacore::Float> &y1,
		                               const casacore::Vector<casacore::uInt> &colors) = 0;
		virtual void drawColoredPoints(const casacore::Vector<casacore::Double> &x1,
		                               const casacore::Vector<casacore::Double> &y1,
		                               const casacore::Vector<casacore::uInt> &colors) = 0;
		virtual void drawColoredPoints(const casacore::Matrix<casacore::Int> &xy,
		                               const casacore::Vector<casacore::uInt> &colors) {
			drawColoredPoints(xy.column(0), xy.column(1), colors);
		}
		virtual void drawColoredPoints(const casacore::Matrix<casacore::Float> &xy,
		                               const casacore::Vector<casacore::uInt> &colors) {
			drawColoredPoints(xy.column(0), xy.column(1), colors);
		}
		virtual void drawColoredPoints(const casacore::Matrix<casacore::Double> &xy,
		                               const casacore::Vector<casacore::uInt> &colors) {
			drawColoredPoints(xy.column(0), xy.column(1), colors);
		}
		// </group>

		// (Cacheable) Draw a set of lines, specifying a color per line to be drawn.
		// <group>
		virtual void drawColoredLines(const casacore::Vector<casacore::Int> &x1,
		                              const casacore::Vector<casacore::Int> &y1,
		                              const casacore::Vector<casacore::Int> &x2,
		                              const casacore::Vector<casacore::Int> &y2,
		                              const casacore::Vector<casacore::uInt> &colors) = 0;
		virtual void drawColoredLines(const casacore::Vector<casacore::Float> &x1,
		                              const casacore::Vector<casacore::Float> &y1,
		                              const casacore::Vector<casacore::Float> &x2,
		                              const casacore::Vector<casacore::Float> &y2,
		                              const casacore::Vector<casacore::uInt> &colors) = 0;
		virtual void drawColoredLines(const casacore::Vector<casacore::Double> &x1,
		                              const casacore::Vector<casacore::Double> &y1,
		                              const casacore::Vector<casacore::Double> &x2,
		                              const casacore::Vector<casacore::Double> &y2,
		                              const casacore::Vector<casacore::uInt> &colors) = 0;
		// </group>

		// Draw a single ellipse using the current pen (ie. color,
		// thickness, style).  The x and y location must be given, along
		// with the semi-major and -minor axis lengths, and the position
		// angle measured in degrees positive from the x axis in a
		// counter-clockwise direction. If outline is false, the
		// ellipse is solid filled, else it is just outlined.
		// xstretch, ystretch should be left defaulted to 1 in most cases;
		// see usage example in WorldCanvas::drawBeamEllipse(), where they are
		// used to stretch a beam ellipse when the display is also linearly
		// stretched away from the aspect ratio natural to the sky projection.
		// They multiply the relative x and y screen coordinates of ellipse
		// points, _before_ they are added to cx and cy.  xstretch, ystretch
		// can be negative (though smajor, sminor normally would not be).
		virtual void drawEllipse(const casacore::Float &cx, const casacore::Float &cy,
		                         const casacore::Float &smajor, const casacore::Float &sminor,
		                         const casacore::Float &pangle, casacore::Bool outline = true,
		                         casacore::Float xstretch = 1., casacore::Float ystretch = 1.);

		// Draw a set of colored ellipses, possibly with borders.  The x and
		// y locations must given, along with semi-major and semi-minor
		// axes, and position angle measured in degrees positive from the x
		// axis in a counter-clockwise direction.  The size of the ellipses
		// is globally scaled by the scale factor, and if <src>outline</src>
		// is <src>true</src>, then each ellipse will have an outline in the
		// current pen color.  <group>
		virtual void drawColoredEllipses(const casacore::Matrix<casacore::Float> &centres,
		                                 const casacore::Vector<casacore::Float> &smajor,
		                                 const casacore::Vector<casacore::Float> &sminor,
		                                 const casacore::Vector<casacore::Float> &pangle,
		                                 const casacore::Vector<casacore::uInt> &colors,
		                                 const casacore::Float &scale = 1.0,
		                                 const casacore::Bool &outline = true);
		// </group>

		// vector primitive buffering
		/*
		void bufferPoint(casacore::Int x, casacore::Int y)
		  { vgbuf_.accumPoint(x,y); }
		void bufferLine(casacore::Int x1, casacore::Int y1, casacore::Int x2, casacore::Int y2)
		  { vgbuf_.accumLine(x1,y1,x2,y2); }
		void bufferPolylinePoint(casacore::Int x, casacore::Int y)
		  { vgbuf_.accumPolylinePoint(x,y); }
		void bufferPolygonPoint(casacore::Int x, casacore::Int y)
		  { vgbuf_.accumPolygonPoint(x,y); }
		*/
		void bufferPoint(casacore::Float x, casacore::Float y) {
			vgbuf_.accumPoint(x,y);
		}
		void bufferLine(casacore::Float x1, casacore::Float y1, casacore::Float x2, casacore::Float y2) {
			vgbuf_.accumLine(x1,y1,x2,y2);
		}
		void bufferPolylinePoint(casacore::Float x, casacore::Float y) {
			vgbuf_.accumPolylinePoint(x,y);
		}
		void bufferPolygonPoint(casacore::Float x, casacore::Float y) {
			vgbuf_.accumPolygonPoint(x,y);
		}
		void flushBuffer() {
			vgbuf_.flush();
		}

		// Set Graphics Attributes
		// Options for functions with enum argument
		// listed in <linkto class=Display>DisplayEnums</linkto>
		// <group>
		virtual void setDrawFunction(Display::DrawFunction function) = 0;
		virtual void setForeground(casacore::uLong color) = 0;
		virtual void setBackground(casacore::uLong color) = 0;
		//virtual void setLineWidth(casacore::uInt width) = 0;
		virtual void setLineWidth(casacore::Float width) = 0;
		virtual void setLineStyle(Display::LineStyle style) = 0;
		virtual void setCapStyle(Display::CapStyle style) = 0;
		virtual void setJoinStyle(Display::JoinStyle style) = 0;
		virtual void setFillStyle(Display::FillStyle style) = 0;
		virtual void setFillRule(Display::FillRule rule) = 0;
		virtual void setArcMode(Display::ArcMode mode) = 0;
		// </group>

		// Get Graphics Attributes
		// <group>
		virtual Display::DrawFunction getDrawFunction() const = 0;
		virtual casacore::uLong                 getForeground()   const = 0;
		virtual casacore::uLong                 getBackground()   const = 0;
		//virtual casacore::uInt                  getLineWidth()    const = 0;
		virtual casacore::Float                 getLineWidth()    const = 0;
		virtual Display::LineStyle    getLineStyle()    const = 0;
		virtual Display::CapStyle     getCapStyle()     const = 0;
		virtual Display::JoinStyle    getJoinStyle()    const = 0;
		virtual Display::FillStyle    getFillStyle()    const = 0;
		virtual Display::FillRule     getFillRule()     const = 0;
		virtual Display::ArcMode      getArcMode()      const = 0;
		// </group>

		// (Cacheable) Option Control
		// Options listed in <linkto class=Display>DisplayEnums</linkto>
		// <group>
		virtual casacore::Bool enable(Display::Option option) = 0;
		virtual casacore::Bool disable(Display::Option option) = 0;
		// </group>

		// Control the image-caching strategy
		virtual void setImageCacheStrategy(Display::ImageCacheStrategy strategy) = 0;
		virtual Display::ImageCacheStrategy imageCacheStrategy() const = 0;

		// (Cacheable) Setup the clip window.  The clip window, when enabled, allows
		// a user to clip all graphics output to a rectangular region on
		// the screen.
		// <group>
		virtual void setClipWindow(casacore::Int x1, casacore::Int y1, casacore::Int x2, casacore::Int y2) = 0;
		virtual void getClipWindow(casacore::Int &x1, casacore::Int &y1, casacore::Int &x2, casacore::Int &y2) = 0;
		// </group>

		// (Not Cacheable) Redraw the window
		// <group>
		void redraw() {
			refresh();
		}
		virtual void refresh(const Display::RefreshReason &reason =
		                         Display::UserCommand,
		                     const casacore::Bool &explicitrequest = true) = 0;
		// </group>

		// Cause display to flush any graphics commands not yet drawn
		virtual void flush() = 0;

		// (Cacheable) Clear the window using the background color
		// <group>
		virtual void clear() = 0;
		virtual void clear(casacore::Int x1, casacore::Int y1, casacore::Int x2, casacore::Int y2) = 0;
		// </group>

		// (Cacheable) Set the color to use for clearing the display
		// <group>
		virtual void setClearColor(casacore::uInt colorIndex) = 0;
		virtual void setClearColor(const casacore::String &colorname) = 0;
		virtual void setClearColor(float r, float g, float b) = 0;
		// </group>

		// (Not Cacheable) Get the current color to use for clearing the display.
		virtual casacore::uInt clearColor() const = 0;
		virtual void getClearColor(float &r, float &g, float &b) const = 0;

		// Get/set the current foreground/background colors.  These colors
		// should be used when the special Strings "foreground" and "background"
		// are given for a color.
		// <group>
		virtual void setDeviceForegroundColor(const casacore::String colorname) = 0;
		virtual casacore::String deviceForegroundColor() const = 0;
		virtual void setDeviceBackgroundColor(const casacore::String colorname) = 0;
		virtual casacore::String deviceBackgroundColor() const = 0;
		// </group>

		// Return the width of the PixelCanvas in pixels
		virtual casacore::uInt width() const = 0;
		// Return the height of the PixelCanvas in pixels
		virtual casacore::uInt height() const = 0;
		// Return the depth of the PixelCanvas in bits
		virtual casacore::uInt depth() const = 0;
		// Get the pixel density (in dots per inch [dpi]) of the PixelCanvas
		virtual void pixelDensity(casacore::Float &xdpi, casacore::Float &ydpi) const = 0;

		// (Cacheable) Set current color (works in RGB or colormap mode)
		// <group>
		virtual void setColor(casacore::uInt colorIndex) = 0;
		virtual void setColor(const casacore::String &colorname) = 0;
		virtual void setRGBColor(float r, float g, float b) = 0;
		virtual void setHSVColor(float h, float s, float v);
		// </group>

		// Get color components in range 0 to 1 without actually
		// allocating the color.  This is needed to set up other
		// devices, for example PgPlot.
		virtual casacore::Bool getColorComponents(const casacore::String &colorname, casacore::Float &r,
		                                casacore::Float &g, casacore::Float &b) = 0;

		// (Not Cacheable) Returns the current color as a color index
		virtual casacore::uInt color() const = 0;

		// (Not Cacheable) Retuns the current color as an RGB triple
		virtual void getColor(float &r, float &g, float &b) const = 0;

		// (Not Cacheable) Get color index value (works in RGB or colormap mode)
		// <group>
		virtual casacore::Bool getColor(casacore::Int x, casacore::Int y, casacore::uInt &color) = 0;
		virtual casacore::Bool getRGBColor(casacore::Int x, casacore::Int y, float &r, float &g, float &b) = 0;
		virtual casacore::Bool getHSVColor(casacore::Int x, casacore::Int y, float &h, float &s, float &v);
		// </group>

		// (Not Cacheable) resize request.  returns true if window was resized.
		// Will refresh if doCallbacks is true.
		//#dk virtual casacore::Bool resize(casacore::uInt reqXSize, casacore::uInt reqYSize, casacore::Bool doCallbacks = true) = 0;
		virtual casacore::Bool resize(casacore::uInt /*reqXSize*/, casacore::uInt /*reqYSize*/, casacore::Bool /*doCallbacks*/ = true) {
			return false;
		}

		// (Not Cacheable) resize the colortable by requesting a new number of cells
		virtual casacore::Bool resizeColorTable(casacore::uInt newSize) = 0;

		// (Not Cacheable) resize the colortable by requesting a new RGB/HSV cube
		virtual casacore::Bool resizeColorTable(casacore::uInt nReds, casacore::uInt nGreens, casacore::uInt nBlues) = 0;

		// Need a mechanism to return the PixelCanvasColorTable so
		// drawing functions within classes can operate.
		virtual PixelCanvasColorTable * pcctbl() const = 0;

		virtual void setPcctbl(PixelCanvasColorTable * pcctbl) = 0;

		// (Not Cacheable) set/get the current colormap.  Undefined behavior results from
		// setting the colormap to a map that is not registered.  This cannot be cached
		// because it affects the result of getColormapSize().  Remember that once data
		// is mapped with mapToColor, it is no longer important to worry about which
		// colormap is active for purposes of drawing the image or vectors.
		// colormapRegistered() can be used to determine whether the currently set
		// colormap() is indeed registered/installed (i.e., color cells have been
		// allocated within the pcctbl).
		// <group>
		void setColormap(Colormap * map);
		Colormap * colormap() const {
			return colormap_;
		}
		casacore::Bool colormapRegistered() {
			return pcctbl()->member(colormap_);
		}
		// </group>

		// register a colormap to the pixelcanvas.  Registration counts are
		// remembered, so that map A is guaranteed to be available as long as
		// the number of register(A)'s exceed the number of unregister(A)'s.
		// Requests are forwarded to the PixelCanvas' PixelCanvasColorTable.
		void registerColormap(Colormap * dcmap, casacore::Float weight = 1.0);

		// Register the <src>cmap</src> Colormap on the PixelCanvas,
		// replacing the <src>cmapToReplace</src> Colormap if possible.
		void registerColormap(Colormap *cmap, Colormap *cmapToReplace);

		// unregister a colormap from a pixelcanvas.  Unregistering the
		// current map may result in undefined behavior.
		void unregisterColormap(Colormap * dcmap);

		// return the size of the current colormap
		casacore::uInt getColormapSize() const {
			return pcctbl()->getColormapSize(colormap_);
		}

		// map [0,N-1] into colorpixels, where N is the current colormap size
		// The values are returned as unsigned integers in their respective
		// array.
		// <note role="tip">The choice of what type to use should be guided by
		// the number of graphics bitplanes available.  For most systems with
		// 8-bit color, casacore::uChar is optimal.  Some systems with 12 bits per pixel
		// with an alpha channel may require using the uLong. </note>
		//
		// <note role="warning">casacore::uChar type may not have enough bits
		// to hold the pixel index on some high-end graphics systems </note>
		// <note role="warning">casacore::uShort type may not have enough bits
		// to hold the pixel index on some high-end graphics systems </note>
		// <group>

		void mapToColor(casacore::Array<casacore::uChar> &outArray,
		                const casacore::Array<casacore::uChar> &inArray, casacore::Bool rangeCheck = false) {
			pcctbl()->mapToColor(colormap_, outArray, inArray, rangeCheck);
		}
		void mapToColor(casacore::Array<casacore::uShort> &outArray,
		                const casacore::Array<casacore::uShort> &inArray, casacore::Bool rangeCheck = false) {
			pcctbl()->mapToColor(colormap_, outArray, inArray, rangeCheck);
		}
		void mapToColor(casacore::Array<casacore::uInt> &outArray,
		                const casacore::Array<casacore::uInt> &inArray, casacore::Bool rangeCheck = false) {
			pcctbl()->mapToColor(colormap_, outArray, inArray, rangeCheck);
		}
		void mapToColor(casacore::Array<casacore::uLong> &outArray,
		                const casacore::Array<casacore::uLong> &inArray, casacore::Bool rangeCheck = false) {
			pcctbl()->mapToColor(colormap_, outArray, inArray, rangeCheck);
		}
		// </group>

		// same as above except the matrix is operated on in place.  Only unsigned
		// values make sense here.
		// <group>
		void mapToColor(casacore::Array<casacore::uChar> &inOutArray, casacore::Bool rangeCheck = false) {
			pcctbl()->mapToColor(colormap_, inOutArray, rangeCheck);
		}
		void mapToColor(casacore::Array<casacore::uShort> &inOutArray, casacore::Bool rangeCheck = false) {
			pcctbl()->mapToColor(colormap_, inOutArray, rangeCheck);
		}
		void mapToColor(casacore::Array<casacore::uInt> &inOutArray, casacore::Bool rangeCheck = false) {
			pcctbl()->mapToColor(colormap_, inOutArray, rangeCheck);
		}
		void mapToColor(casacore::Array<casacore::uLong> &inOutArray, casacore::Bool rangeCheck = false) {
			pcctbl()->mapToColor(colormap_, inOutArray, rangeCheck);
		}
		// </group>



		void mapToColor(casacore::Array<casacore::uInt> &outArray,
				const casacore::Array<casacore::uInt> &inArrayRed, const casacore::Array<casacore::uInt> &inArrayGreen, const casacore::Array<casacore::uInt>& inArrayBlue) {
			pcctbl()->mapToColorRGB(colormap_, outArray, inArrayRed, inArrayGreen, inArrayBlue);
		}


		// Multi-Channel functions that combine separate array channels into
		// a single array of output colors for use with functions that take
		// color values.
		// These two functions accept color values from [0-1] for each channel.
		// The colorModel is compared with the model of the PixelCanvasColorTable
		// and is used to perform the correct mapping.
		// <group>
		void mapToColor3(casacore::Array<casacore::uLong> &out,
		                 const casacore::Array<casacore::Float> &chan1in,
		                 const casacore::Array<casacore::Float> &chan2in,
		                 const casacore::Array<casacore::Float> &chan3in);
		void mapToColor3(casacore::Array<casacore::uLong> &out,
		                 const casacore::Array<casacore::Double> &chan1in,
		                 const casacore::Array<casacore::Double> &chan2in,
		                 const casacore::Array<casacore::Double> &chan3in);
		// </group>

		// This one maps values between 0 and the integer maximum value for
		// each channel into a single output image suitable for
		// PixelCanvas::drawImage.
		// <group>
		virtual void mapToColor3(casacore::Array<casacore::uLong> &out,
		                         const casacore::Array<casacore::uInt> &chan1in,
		                         const casacore::Array<casacore::uInt> &chan2in,
		                         const casacore::Array<casacore::uInt> &chan3in);
		// </group>


		// save/restore the current translation.  This is called pushMatrix because
		// eventually we may want scaling or rotation to play a modest
		// role here.
		// <group>
		virtual void pushMatrix() = 0;
		virtual void popMatrix() = 0;
		// </group>
		// zero the current translation
		virtual void loadIdentity() = 0;

		// translation functions
		// translate applies a relative translation to the current matrix and
		// can be used to position graphics.  Together with pushMatrix and
		// popMatrix it can be used to build heirarchical scenes.
		// <group>
		virtual void translate(casacore::Int xt, casacore::Int yt) = 0;
		virtual void getTranslation(casacore::Int &xt, casacore::Int &yt) const = 0;
		virtual casacore::Int xTranslation() const = 0;
		virtual casacore::Int yTranslation() const = 0;
		// </group>

		// return the drawing buffer, the target destination for graphics
		Display::DrawBuffer drawBuffer() const {
			return drawBuffer_;
		}
		// (Not cacheable) set the draw buffer
		virtual void setDrawBuffer(Display::DrawBuffer buf) = 0;
		// buffer memory exchanges
		// (Not cacheable)
		// <group>
		virtual void copyBackBufferToFrontBuffer() = 0;
		virtual void copyFrontBufferToBackBuffer() = 0;
		virtual void swapBuffers() = 0;
		// </group>

		// partial buffer memory exchanges.  (x1,y1 are blc, x2,y2 are trc)
		// <group>
		virtual void copyBackBufferToFrontBuffer(casacore::Int x1, casacore::Int y1, casacore::Int x2, casacore::Int y2) = 0;
		virtual void copyFrontBufferToBackBuffer(casacore::Int x1, casacore::Int y1, casacore::Int x2, casacore::Int y2) = 0;
		virtual void swapBuffers(casacore::Int x1, casacore::Int y1, casacore::Int x2, casacore::Int y2) = 0;
		// </group>

		// return the drawmode (Compile or Draw)
		// Compile drawmode means that a display list is currently being built
		// Draw drawmode means that drawing commands are not cached but
		// instead are sent to the display.  This command is controlled by
		// the caching state using newList()
		Display::DrawMode drawMode() const {
			return drawMode_;
		}

		// return the colorModel used for multichannel color
		Display::ColorModel colorModel() const {
			return colorModel_;
		}
		// Set the input color model for multichannel color
		void setColorModel(Display::ColorModel colorModel);


		// return true if the refresh is active (Added for X11PixelCanvas DefaultBuffer)
		casacore::Bool refreshActive() const {
			return refreshActive_;
		}

		// return true if refresh is allowed right now...
		virtual casacore::Bool refreshAllowed() const {
			return true;
		}

		virtual casacore::Float pixelScaling() const {
			return 1.0;
		}

	protected:

		// Abstract base class idiom
		PixelCanvas();
		PixelCanvas( const PixelCanvas * );
		PixelCanvas(PixelCanvasColorTable * pcctbl);

		// Only allowed by derived classes
		void setDrawMode(Display::DrawMode mode) {
			drawMode_ = mode;
		}

		// Also allowed only by derived classes
		void setDrawBuffer_(Display::DrawBuffer buf) {
			drawBuffer_ = buf;
		}

		// true if the PixelCanvas was the one who
		// registered the existing colortable
		casacore::Bool defaultColormapActive_;

	private:
		casacore::Matrix<casacore::Float> getMarker(const Display::Marker& marker, const casacore::Float& pixelHeight);

		// This is the current colormap.
		Colormap * colormap_;


		// This is the current drawing mode
		Display::DrawMode drawMode_;

		// The current drawing buffer
		Display::DrawBuffer drawBuffer_;

		// The current color cube model to be used for all
		// multichannel color mapping.
		Display::ColorModel colorModel_;

		// true if refresh is active
		casacore::Bool refreshActive_;

		// Number of colormaps registered by external clients
		casacore::uInt nRegisteredColormaps_;

		// This is the list of registered refresh EventHandlers
		std::list<void *> refreshEHList_;
		// This is the list of registered motion EventHandlers
		std::list<void *> motionEHList_;
		// This is the list of registered position EventHandlers
		std::list<void *> positionEHList_;

		// The PCVGBuffer is used to accumulate lines and points until
		// flushed (or the buffer is full)
		PCVGBuffer vgbuf_;
	};


} //# NAMESPACE CASA - END

#endif


