//# PSPixelCanvasColorTable.cc: PostScript version of PixelCanvasColorTable
//# Copyright (C) 1999,2000
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

#ifndef TRIALDISPLAY_PSPIXELCANVASCOLORTABLE_H
#define TRIALDISPLAY_PSPIXELCANVASCOLORTABLE_H




#include <casa/aips.h>

#include <casa/Arrays/Vector.h>
#include <display/Display/DisplayEnums.h>
#include <display/Display/PixelCanvasColorTable.h>
#include <display/Display/PSDriver.h>


namespace casa { //# NAMESPACE CASA - BEGIN

// <summary>
// Implementation of PixelCanvasColorTable for PostScript device.
// </summary>
//
// <prerequisite>
// <li> <linkto class="PixelCanvasColorTable">PixelCanvasColorTable</linkto>
// <li> <linkto class="PSDriver">PSDriver</linkto>
// </prerequisite>
//
// <etymology>
// </etymology>
//
// <synopsis>
// To create a PSPixelCanvasColorTable, just pass the constructor a pointer
// to a
// <linkto class="PSDriver">PSDriver</linkto> and, optionally,
// supplying a color model. (The default is Index).
//<note role=tip> Unlike the <em>X11PixelCanvasColorTable</em>,
// PSPixelCanvasColorTable allows changing the color model on the fly.
//</note>
// PSPixelCanvasColorTable is not likely to be explicitly used by other
// than Display Library developers, particularly those creating
// "WorldCanvasApp"s. One exception is using PSPixelCanvasColorTable in non
// Index mode. Since
// <linkto class="PSWorldCanvasApp">PSWorldCanvasApp</linkto>
//  creates its PSPixelCanvasColorTable in
// Index mode, it will be necessary to get a pointer to the object and
// explicitly change to different color modes.
//<example>
//<srcblock>
// psapp = new PSWorldCanvasApp(psdriver);
// wCanvas = psapp->worldCanvas();
// pCanvas = (PSPixelCanvas *)wCanvas->pixelCanvas();
// PSPixelCanvasColorTable *psctbl = pCanvas->PSpcctbl();
// psctbl->setColorModel(Display::RGB);
//</srcblock>
// See Display/test/dMultichannelRaster.cc for an example.
//</example>
// </synopsis>
//
//<note role=tip> PostScript supports a 4096 entry color table for
// indexed color. PSPixelCanvasColorTable logically breaks this into two
// parts. One part is used for the changable colors. The other part is
// reserved for read only colors. (Those allocated by name). The number
// of read only colors is 512.
//</note>
//
// <motivation>
// </motivation>
//
// <todo>
// </todo>
//<use visibility=local>

	class PSPixelCanvasColorTable : public PixelCanvasColorTable {
	public:

		//# PostScript Level 2 definition allows 12 bits/component, we
		//# use 10 in RGB mode due to limitations in the DL spec.
		enum {	INDEXBPC = 12, RGBBPC=10, INDEXCOLORS = (1<<INDEXBPC),
		        RGBCOLORS = (1<<RGBBPC)
		     };
		// The last NUMROCOLORS of the table are reserved for RO values.
//#	NUMROCOLORS could be changed as long as enough RO slots remained.
		enum { NUMROCOLORS=512, NUMRWCOLORS=INDEXCOLORS-NUMROCOLORS};
		enum {	BMASK = 0xff, RGBMASK= (RGBCOLORS-1) ,
		        INDEXMASK = (INDEXCOLORS -1)
		     };
		// Amount to shift components when mapping RGB.
		enum { RSHIFT=RGBBPC*2, GSHIFT=RGBBPC, BSHIFT=0};
		// <group>
		PSPixelCanvasColorTable(PSDriver *ps,
		                        const Display::ColorModel = Display::Index);
		// </group>
		virtual ~PSPixelCanvasColorTable();

		// Resize the map if allowed.  Returns true if resize was accepted
		// <group>
		virtual casacore::Bool resize(casacore::uInt newSize);
		virtual casacore::Bool resize(casacore::uInt nReds, casacore::uInt nGreens, casacore::uInt nBlues);
		// </group>

		// Install colors into the color table.  Offset is zero-based.  Colors
		// are installed into the PixelCanvasColorTable until the Arrays run out
		// or until the end of the colortable is reached. Can be called in any
		// mode, but only affects graphics drawn in Index mode.
		// Values are clamped to [0.0,1.0].
//<thrown><li> casacore::AipsError </thrown>
		virtual casacore::Bool installRGBColors(const casacore::Vector<casacore::Float> & r,
		                              const casacore::Vector<casacore::Float> & g,
		                              const casacore::Vector<casacore::Float> & b,
		                              const casacore::Vector<casacore::Float> & alpha,
		                              casacore::uInt offset = 0);

		// Return the number of RW colors in the colormap.
		virtual casacore::uInt nColors() const;

		// Return the number of colors per component.
		// For RGB/HSV, returns the number of colors/component supported by PostScript.
		// For Index, returns the number of colors/component for the lookup table.
		//  (Limited by D.L. spec).
		virtual void nColors(casacore::uInt &n1, casacore::uInt &n2, casacore::uInt &n3) const;

		// Maximum number of colors in color table.
		casacore::uInt maxColors()const;
		// Return the depth in bits of the colors.
		virtual casacore::uInt depth() const;

		// Return the number of colors that are still unallocated.
		virtual casacore::uInt nSpareColors() const;

		// map [0,N-1] into colorpixels, where N is the current colormap size
		// The values are returned as unsigned integers in their respective
		// array.
		//
		// <note role="Warning">casacore::uChar type doesn't have enough bits
		// to hold the pixel index. </note>
		// <note role="Warning">casacore::uChar and casacore::uShort don't have enough bits to
		// hold RGB or HSV values.</note>
		// <group>
		virtual void mapToColor(const Colormap * map, casacore::Array<casacore::uChar> & outArray,
		                        const casacore::Array<casacore::uChar> & inArray, casacore::Bool rangeCheck = true) const;
		virtual void mapToColor(const Colormap * map, casacore::Array<casacore::uShort> & outArray,
		                        const casacore::Array<casacore::uShort> & inArray, casacore::Bool rangeCheck = true) const;
		virtual void mapToColor(const Colormap * map, casacore::Array<casacore::uInt> & outArray,
		                        const casacore::Array<casacore::uInt> & inArray, casacore::Bool rangeCheck = true) const;
		virtual void mapToColor(const Colormap * map, casacore::Array<casacore::uLong> & outArray,
		                        const casacore::Array<casacore::uLong> & inArray, casacore::Bool rangeCheck = true) const;
		// </group>

		// Same as above except the matrix is operated on in place. Only unsigned
		// values make sense here.
		// <group>
		virtual void mapToColor(const Colormap * map, casacore::Array<casacore::uChar> & inOutArray,
		                        casacore::Bool rangeCheck = true) const;
		virtual void mapToColor(const Colormap * map, casacore::Array<casacore::uShort> & inOutArray,
		                        casacore::Bool rangeCheck = true) const;
		virtual void mapToColor(const Colormap * map, casacore::Array<casacore::uInt> & inOutArray,
		                        casacore::Bool rangeCheck = true) const;
		virtual void mapToColor(const Colormap * map, casacore::Array<casacore::uLong> & inOutArray,
		                        casacore::Bool rangeCheck = true) const;
		// </group>

		// (Multichannel Color)
		// Merge separate channel data into an output image.
		// This function maps floating values between 0 and 1
		// into a output image suitable for PixelCanvas::drawImage().
		// <group>
		virtual void mapToColor3(casacore::Array<casacore::uLong> & out,
		                         const casacore::Array<casacore::Float> & chan1in,
		                         const casacore::Array<casacore::Float> & chan2in,
		                         const casacore::Array<casacore::Float> & chan3in);
		virtual void mapToColor3(casacore::Array<casacore::uLong> & out,
		                         const casacore::Array<casacore::Double> & chan1in,
		                         const casacore::Array<casacore::Double> & chan2in,
		                         const casacore::Array<casacore::Double> & chan3in);
		// </group>

		// This one maps values between 0 and the integer
		// maximum value for each channel into a single
		// output image suitable for PixelCanvas::drawImage().
		// <group>
		virtual void mapToColor3(casacore::Array<casacore::uLong> & out,
		                         const casacore::Array<casacore::uShort> & chan1in,
		                         const casacore::Array<casacore::uShort> & chan2in,
		                         const casacore::Array<casacore::uShort> & chan3in);
		virtual void mapToColor3(casacore::Array<casacore::uLong> & out,
		                         const casacore::Array<casacore::uInt> & chan1in,
		                         const casacore::Array<casacore::uInt> & chan2in,
		                         const casacore::Array<casacore::uInt> & chan3in);
		// </group>

		// Convert from a packed array of RGB triples to an array of color values.
		// The output array needs to be 3 times as long as the input array.
		// Used by PSPixelCanvas to convert from D.L. RGB format to an array
		// the PostScript driver can use.
		void mapFromColor3(const casacore::Array<casacore::uLong> & inArray,
		                   casacore::Array<casacore::uShort> & outArray) const;

		// (Multichannel Color)
		// Transform arrays from the passed color model into
		// the colormodel of the PSPCCT.
		// Does nothing if colorModel is Display::Index.
		// It is assumed that input arrays are in the range of [0,1]
		virtual casacore::Bool colorSpaceMap(Display::ColorModel,
		                           const casacore::Array<casacore::Float> & chan1in,
		                           const casacore::Array<casacore::Float> & chan2in,
		                           const casacore::Array<casacore::Float> & chan3in,
		                           casacore::Array<casacore::Float> & chan1out,
		                           casacore::Array<casacore::Float> & chan2out,
		                           casacore::Array<casacore::Float> & chan3out);


		// Return the color model for multichannel color
		Display::ColorModel colorModel() const;
		// Returns the current # of color components (1 for Indexed, 3 for RGB/HSV).
		casacore::uInt numComponents()const;
		// Changeable at any time.
		virtual void setColorModel(const Display::ColorModel);

		PSDriver *getPSDriver()const {
			return ps;
		}

//#  ////////////////////////////////////////////////////////////////
//			X11 emulation routines
//#  ////////////////////////////////////////////////////////////////
//<srcblock>
		// Convert a colorname to a color triple. Returns true for success,
		// false if name can't be found. The color spec can also be in the form:
		//  "#xxxxxx"		 A '#' character followed by exactly 6 hex digits.
		// 			(This form is considered obsolete and is not
		//			 completely implemented here).
		//  "rgb:<red>/<green>/<blue>"	Where <red>, <green> and <blue> are
		//			Each 1 to 4 hex digits. Each value is divided
		//			by 1.0/(2^n -1) where n is the # of hex chars in
		//			the term. The result is 3 floating point numbers in
		//			the range 0..1.
		// "rgbi:<red>/<green>/<blue>"	Where <red>, <green> and <blue> are
		//					floating point #s in the range 0..1.
//</srcblock>
		// See <em>XParseColor</em> for more information.

//<group>
		static casacore::Bool parseColor(const char *name,
		                       float &red, float &green, float &blue);
		static casacore::Bool parseColor(const casacore::String &name,
		                       float &red, float &green, float &blue);
//</group>

		// Return contents of colormap at the given index. Returns false if
		// the index is out of range. The valid range of index is 0..4095.
		casacore::Bool queryColor(const int index, float &r, float &g, float &b);
		// Sets the contents of colormap at the given index. Returns false if
		// the index is out of range. The valid range of index is 0..nColors().
		// ( Can't change read only values).
		casacore::Bool storeColor(const int index,
		                const float r, const float g, const float b);
		// Allocate the color value in the color table. index is set to the
		// index allocated. Returns true for success, else false.
		casacore::Bool allocColor(const float r, const float g, const float b, int &index);
		casacore::Bool allocColor(const casacore::String &name, int &index);
		casacore::Bool allocColor(const char *name, int &index);

		// Whether to put tracing comments in the output. This may be helpful
		// when trying to decipher the PostScript file.
		casacore::Bool annotate()const {
			return annotate_;
		}
		void annotate(const casacore::Bool a) {
			annotate_ = a;
		}

		// print details of class to ostream
		friend casacore::ostream & operator << (casacore::ostream & os,
		                              const PSPixelCanvasColorTable & pcc);

		// Convert a packed pixel (from mapToColor3) to three color components.
		static inline void pixelToComponents( const casacore::uLong pixel,
		                                      casacore::uShort &r, casacore::uShort &g, casacore::uShort &b) {
			r =  (pixel >> RSHIFT) & RGBMASK;
			g = (pixel >> GSHIFT) & RGBMASK;
			b = (pixel >> BSHIFT) & RGBMASK;
		}
		// Pack RGB or HSV color components into a single unsigned long.
		static inline void componentsToPixel(	const casacore::uShort r, const casacore::uShort g,
		                                        casacore::uShort &b, casacore::uLong &pixel) {
			pixel = ((r & RGBMASK) << RSHIFT) | ((g & RGBMASK) << GSHIFT)
			        | ((b & RGBMASK) << BSHIFT);
		}

	private:
		PSPixelCanvasColorTable();
		void pspcinit(PSDriver *ps, const Display::ColorModel);
		// Finds the index of a color triple coming 'close' to the RGB args.
		// Returns: true if a match is found, else false.
		// If lookROColor finds a match with a deallocated cell, it reallocates it.
		casacore::Bool lookupROColor(const float r, const float g, const float b, int &index);
		casacore::Bool lookupRWColor(const float r, const float g, const float b, int &index);
		int allocColor_(const float r, const float g, const float b, int &index);

		// Mark a RO color as unallocated.
		void deallocate(casacore::uLong index);

		// (Valid Always) number of total colors available. Changed by resize.
		casacore::uInt nColors_;
		// Number of bits/per component. Determined by colorModel_.
		casacore::uInt bpc_;

		// (Valid Always)
		// The colormodel that this PSPixelCanvasColorTable is currently
		// configured as.
		Display::ColorModel colorModel_;

		// PS
		PSDriver	*ps;
		// Copies of the color table.
		float	red_[INDEXCOLORS],
		        blue_[INDEXCOLORS],
		        green_[INDEXCOLORS];
		// true if index has been allocated.
		casacore::Bool	allocated_[NUMROCOLORS];
		casacore::Bool	annotate_;
	};


} //# NAMESPACE CASA - END

#endif
