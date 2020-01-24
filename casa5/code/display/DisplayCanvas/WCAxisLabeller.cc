//# WCAxisLabeller.cc: base class for labelling axes on the WorldCanvas
//# Copyright (C) 1999,2000,2001,2002
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

//# aips includes:
#include <casa/aips.h>
#include <casa/Arrays/Vector.h>
#include <casa/Containers/Record.h>
#include <casa/Utilities/Regex.h>
#include <casa/System/Aipsrc.h>
#include <casa/System/AipsrcValue.h>

//# trial includes:

//# display library includes:

//# this include:
#include <display/DisplayCanvas/WCAxisLabeller.h>

using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN

	const String WCAxisLabeller::LABEL_CHAR_SIZE = "labelcharsize";
	const String WCAxisLabeller::PLOT_TITLE="titletext";

	WCAxisLabeller::WCAxisLabeller() {
		// set defaults from .aipsrc, if they exist there.
		// Otherwise, use the hard-coded values below as defaults.
		String onoff;
		Aipsrc::find(onoff,"display.axislabels","on");
		itsDefaultSwitch = ! onoff.matches(Regex(
		                                       "[ \t]*(([nN]o)|([oO]ff)|([fF](alse)*))[ \t\n]*"   ));

		AipsrcValue<Float>::find(itsDefaultCharSize,"display.axislabels.charsize",
		                         1.2f);
		titleChanged = false;
	}

	WCAxisLabeller::~WCAxisLabeller()
	{}

	void WCAxisLabeller::setDefaultOptions() {
		itsOptionsAxisLabelSwitch = itsDefaultSwitch;
		itsOptionsCharSize = itsDefaultCharSize;

		itsOptionsTitleText = String("");
		itsOptionsTitleTextColor = String("foreground");
		itsOptionsXAxisText = String("X axis");
		//itsOptionsXAxisTextUnset = false;
		itsOptionsXAxisTextUnset = true;
		itsOptionsXAxisTextColor = String("foreground");
		itsOptionsYAxisText = String("Y axis");
		//itsOptionsYAxisTextUnset = false;
		itsOptionsYAxisTextUnset = true;
		itsOptionsYAxisTextColor = String("foreground");
		itsOptionsXGridType = String("Tick marks");
		itsOptionsXGridColor = String("foreground");
		itsOptionsYGridType = String("Tick marks");
		itsOptionsYGridColor = String("foreground");
		itsOptionsTickLength = Float(4.0);
		itsOptionsLabelPos = String("Auto");
		itsOptionsPlotOutline = true;
		itsOptionsPlotOutlineColor = String("foreground");
		itsOptionsCharFont = String("normal");
		AipsrcValue<Float>::find(itsOptionsLineWidth,"display.axislabels.linewidth",
		                         1.4f);
	}

	Bool WCAxisLabeller::setOptions(const Record &rec, Record &) {
		Bool ret = false;
		Bool localchange = false;

		Bool error;
		localchange = (readOptionRecord(itsOptionsAxisLabelSwitch, error,
		                                rec, "axislabelswitch")
		               || localchange);
		localchange = (readOptionRecord(itsOptionsTitleText, error,
		                                rec, PLOT_TITLE) || localchange);
		localchange = (readOptionRecord(itsOptionsTitleTextColor, error,
		                                rec, "titletextcolor") || localchange);
		localchange = (readOptionRecord(itsOptionsXAxisText,
		                                itsOptionsXAxisTextUnset,
		                                error,
		                                rec, "xaxistext") || localchange);
		localchange = (readOptionRecord(itsOptionsXAxisTextColor, error,
		                                rec, "xaxistextcolor") || localchange);
		localchange = (readOptionRecord(itsOptionsYAxisText,
		                                itsOptionsYAxisTextUnset,
		                                error,
		                                rec, "yaxistext") || localchange);
		localchange = (readOptionRecord(itsOptionsYAxisTextColor, error,
		                                rec, "yaxistextcolor") || localchange);
		localchange = (readOptionRecord(itsOptionsXGridType, error,
		                                rec, "xgridtype") || localchange);
		localchange = (readOptionRecord(itsOptionsXGridColor, error,
		                                rec, "xgridcolor") || localchange);
		localchange = (readOptionRecord(itsOptionsYGridType, error,
		                                rec, "ygridtype") || localchange);
		localchange = (readOptionRecord(itsOptionsYGridColor, error,
		                                rec, "ygridcolor") || localchange);
		localchange = (readOptionRecord(itsOptionsTickLength, error,
		                                rec, "ticklength") || localchange);
		localchange = (readOptionRecord(itsOptionsPlotOutline, error,
		                                rec, "plotoutline") || localchange);
		localchange = (readOptionRecord(itsOptionsPlotOutlineColor, error,
		                                rec, "plotoutlinecolor") ||
		               localchange);
		localchange = (readOptionRecord(itsOptionsLabelPos, error,
		                                rec, "labelposition") || localchange);
		localchange = (readOptionRecord(itsOptionsCharSize, error,
		                                rec, LABEL_CHAR_SIZE) || localchange);
		localchange = (readOptionRecord(itsOptionsCharFont, error,
		                                rec, "labelcharfont") || localchange);
		localchange = (readOptionRecord(itsOptionsLineWidth, error,
		                                rec, "labellinewidth") || localchange);

		if (localchange) {
			// invalidate existing draw lists, etc.
			invalidate();
		}

		ret = (ret || localchange);
		return ret;
	}

	Record WCAxisLabeller::getOptions() const {
		Record rec;

		Record axislabelswitch;
		axislabelswitch.define("dlformat", "axislabelswitch");
		axislabelswitch.define("listname", "axis labelling & annotation?");
		axislabelswitch.define("ptype", "boolean");
		axislabelswitch.define("default", itsDefaultSwitch);
		axislabelswitch.define("value", itsOptionsAxisLabelSwitch);
		axislabelswitch.define("allowunset", false);
		axislabelswitch.define("context", "axis_labels");
		rec.defineRecord("axislabelswitch", axislabelswitch);

		Record titletext;
		titletext.define("dlformat", "titletext");
		titletext.define("listname", "title");
		titletext.define("ptype", "string");
		titletext.define("default", String(""));
		titletext.define("value", itsOptionsTitleText);
		titletext.define("allowunset", false);
		titletext.define("context", "axis_labels");
		rec.defineRecord(PLOT_TITLE, titletext);

		Record titletextcolor;
		titletextcolor.define("dlformat", "titletextcolor");
		titletextcolor.define("listname", "title color");
		titletextcolor.define("ptype", "userchoice");
		Vector<String> vcolor(11);
		vcolor(0) = "foreground";
		vcolor(1) = "background";
		vcolor(2) = "black";
		vcolor(3) = "white";
		vcolor(4) = "red";
		vcolor(5) = "green";
		vcolor(6) = "blue";
		vcolor(7) = "cyan";
		vcolor(8) = "magenta";
		vcolor(9) = "yellow";
		vcolor(10) = "gray";
		titletextcolor.define("popt", vcolor);
		titletextcolor.define("default", "foreground");
		titletextcolor.define("value", itsOptionsTitleTextColor);
		titletextcolor.define("allowunset", false);
		titletextcolor.define("context", "axis_label_properties");
		rec.defineRecord("titletextcolor", titletextcolor);

		Record xaxistext;
		xaxistext.define("dlformat", "xaxistext");
		xaxistext.define("listname", "x axis label");
		xaxistext.define("ptype", "string");
		//xaxistext.define("default", String("X axis"));
		//Record unset;
		//unset.define("i_am_unset", "i_am_unset");
		xaxistext.defineRecord("default", unset());
		if (itsOptionsXAxisTextUnset) {
			xaxistext.defineRecord("value", unset());
		} else {
			xaxistext.define("value", itsOptionsXAxisText);
		}
		xaxistext.define("allowunset", true);
		xaxistext.define("context", "axis_labels");
		rec.defineRecord("xaxistext", xaxistext);

		Record xaxistextcolor;
		xaxistextcolor.define("dlformat", "xaxistextcolor");
		xaxistextcolor.define("listname", "x axis label color");
		xaxistextcolor.define("ptype", "userchoice");
		xaxistextcolor.define("popt", vcolor);
		xaxistextcolor.define("default", "foreground");
		xaxistextcolor.define("value", itsOptionsXAxisTextColor);
		xaxistextcolor.define("allowunset", false);
		xaxistextcolor.define("context", "axis_label_properties");
		rec.defineRecord("xaxistextcolor", xaxistextcolor);

		Record yaxistext;
		yaxistext.define("dlformat", "yaxistext");
		yaxistext.define("listname", "y axis label");
		yaxistext.define("ptype", "string");
		//yaxistext.define("default", String("Y axis"));
		yaxistext.defineRecord("default", unset());
		if (itsOptionsYAxisTextUnset) {
			yaxistext.defineRecord("value", unset());
		} else {
			yaxistext.define("value", itsOptionsYAxisText);
		}
		yaxistext.define("allowunset", true);
		yaxistext.define("context", "axis_labels");
		rec.defineRecord("yaxistext", yaxistext);

		Record yaxistextcolor;
		yaxistextcolor.define("dlformat", "yaxistextcolor");
		yaxistextcolor.define("listname", "y axis label color");
		yaxistextcolor.define("ptype", "userchoice");
		yaxistextcolor.define("popt", vcolor);
		yaxistextcolor.define("default", "foreground");
		yaxistextcolor.define("value", itsOptionsYAxisTextColor);
		yaxistextcolor.define("allowunset", false);
		yaxistextcolor.define("context", "axis_label_properties");
		rec.defineRecord("yaxistextcolor", yaxistextcolor);

		Record xgridtype;
		xgridtype.define("dlformat", "xgridtype");
		xgridtype.define("listname", "x grid type");
		xgridtype.define("ptype", "choice");
		Vector<String> vgridtype(3);
		vgridtype(0) = "None";
		vgridtype(1) = "Tick marks";
		vgridtype(2) = "Full grid";
		xgridtype.define("popt", vgridtype);
		xgridtype.define("default", "Tick marks");
		xgridtype.define("value", itsOptionsXGridType);
		xgridtype.define("allowunset", false);
		xgridtype.define("context", "axis_labels");
		rec.defineRecord("xgridtype", xgridtype);

		Record xgridcolor;
		xgridcolor.define("dlformat", "xgridcolor");
		xgridcolor.define("listname", "x grid/tick color");
		xgridcolor.define("ptype", "userchoice");
		xgridcolor.define("popt", vcolor);
		xgridcolor.define("default", "foreground");
		xgridcolor.define("value", itsOptionsXGridColor);
		xgridcolor.define("allowunset", false);
		xgridcolor.define("context", "axis_label_properties");
		rec.defineRecord("xgridcolor", xgridcolor);

		Record ygridtype;
		ygridtype.define("dlformat", "ygridtype");
		ygridtype.define("listname", "y grid type");
		ygridtype.define("ptype", "choice");
		ygridtype.define("popt", vgridtype);
		ygridtype.define("default", "Tick marks");
		ygridtype.define("value", itsOptionsYGridType);
		ygridtype.define("allowunset", false);
		ygridtype.define("context", "axis_labels");
		rec.defineRecord("ygridtype", ygridtype);

		Record ygridcolor;
		ygridcolor.define("dlformat", "ygridcolor");
		ygridcolor.define("listname", "y grid/tick color");
		ygridcolor.define("ptype", "userchoice");
		ygridcolor.define("popt", vcolor);
		ygridcolor.define("default", "foreground");
		ygridcolor.define("value", itsOptionsYGridColor);
		ygridcolor.define("allowunset", false);
		ygridcolor.define("context", "axis_label_properties");
		rec.defineRecord("ygridcolor", ygridcolor);

		Record ticklength;
		ticklength.define("dlformat", "ticklength");
		ticklength.define("listname", "tick mark length");
		ticklength.define("ptype", "floatrange");
		ticklength.define("pmin", Float(0.0));
		ticklength.define("pmax", Float(20.0));
		ticklength.define("presolution", Float(0.1));
		ticklength.define("default", Float(4.0));
		ticklength.define("value", itsOptionsTickLength);
		ticklength.define("allowunset", false);
		ticklength.define("context", "axis_labels");
		rec.defineRecord("ticklength", ticklength);

		Record plotoutline;
		plotoutline.define("dlformat", "plotoutline");
		plotoutline.define("listname", "plot border?");
		plotoutline.define("ptype", "boolean");
		plotoutline.define("default", Bool(true));
		plotoutline.define("value", itsOptionsPlotOutline);
		plotoutline.define("allowunset", false);
		plotoutline.define("context", "axis_labels");
		rec.defineRecord("plotoutline", plotoutline);

		Record plotoutlinecolor;
		plotoutlinecolor.define("dlformat", "plotoutlinecolor");
		plotoutlinecolor.define("listname", "plot border color");
		plotoutlinecolor.define("ptype", "userchoice");
		plotoutlinecolor.define("popt", vcolor);
		plotoutlinecolor.define("default", "foreground");
		plotoutlinecolor.define("value", itsOptionsPlotOutlineColor);
		plotoutlinecolor.define("allowunset", false);
		plotoutlinecolor.define("context", "axis_label_properties");
		rec.defineRecord("plotoutlinecolor", plotoutlinecolor);

		Record labelposition;
		labelposition.define("dlformat", "labelposition");
		labelposition.define("listname", "label position");
		labelposition.define("ptype", "choice");
		Vector<String> lblencod(5);
		lblencod(0) = "Auto";
		lblencod(1) = "bottom-left";
		lblencod(2) = "bottom-right";
		lblencod(3) = "top-left";
		lblencod(4) = "top-right";
		labelposition.define("popt", lblencod);
		labelposition.define("default", "Auto");
		labelposition.define("value", itsOptionsLabelPos);
		labelposition.define("allowunset", false);
		labelposition.define("context", "axis_label_properties");
		rec.defineRecord("labelposition", labelposition);

		Record labelcharsize;
		labelcharsize.define("dlformat", LABEL_CHAR_SIZE );
		labelcharsize.define("listname", "character size");
		labelcharsize.define("ptype", "floatrange");
		labelcharsize.define("pmin", Float(0.2));
		labelcharsize.define("pmax", Float(4.0));
		labelcharsize.define("presolution", Float(0.05));
		labelcharsize.define("default", itsDefaultCharSize);
		labelcharsize.define("value", itsOptionsCharSize);
		labelcharsize.define("allowunset", false);
		labelcharsize.define("context", "axis_label_properties");
		rec.defineRecord(LABEL_CHAR_SIZE, labelcharsize);

		Record labelcharfont;
		labelcharfont.define("dlformat", "labelcharfont");
		labelcharfont.define("listname", "character font");
		labelcharfont.define("ptype", "choice");
		Vector<String> vlabelcharfont(4);
		vlabelcharfont(0) = "normal";
		vlabelcharfont(1) = "roman";
		vlabelcharfont(2) = "italic";
		vlabelcharfont(3) = "script";
		labelcharfont.define("popt", vlabelcharfont);
		labelcharfont.define("default", String("normal"));
		labelcharfont.define("value", itsOptionsCharFont);
		labelcharfont.define("allowunset", false);
		labelcharfont.define("context", "axis_label_properties");
		rec.defineRecord("labelcharfont", labelcharfont);

		Record labellinewidth;
		labellinewidth.define("dlformat", "labellinewidth");
		labellinewidth.define("listname", "line width");
		labellinewidth.define("ptype", "floatrange");
		labellinewidth.define("pmin", Float(0.0));
		labellinewidth.define("pmax", Float(5.0));
		labellinewidth.define("presolution", Float(0.1));
		labellinewidth.define("default", Float(0.5));
		labellinewidth.define("value", itsOptionsLineWidth);
		labellinewidth.define("allowunset", false);
		labellinewidth.define("context", "axis_label_properties");
		rec.defineRecord("labellinewidth", labellinewidth);

		return rec;
	}

	Bool WCAxisLabeller::setAxisLabelSwitch(const Bool labelswitch) {
		Bool ret = (itsOptionsAxisLabelSwitch != labelswitch);
		itsOptionsAxisLabelSwitch = labelswitch;
		return ret;
	}

	/*Bool WCAxisLabeller::setTitleText(const String text) {
		Bool ret = (itsOptionsTitleText != text);
		itsOptionsTitleText = text;
		return ret;
	}*/
	void WCAxisLabeller::setSubstituteTitleText( const String substituteImageName){
		if ( substituteTitleText != substituteImageName ){
			substituteTitleText = substituteImageName;
			titleChanged = true;
		}
	}

	Bool WCAxisLabeller::setTitleTextColor(const String color) {
		Bool ret = (itsOptionsTitleTextColor != color);
		itsOptionsTitleTextColor = color;
		return ret;
	}

	Bool WCAxisLabeller::setXAxisText(const String text) {
		Bool ret = (itsOptionsXAxisText != text);
		itsOptionsXAxisText = text;
		itsOptionsXAxisTextUnset = false;
		return ret;
	}

	Bool WCAxisLabeller::setYAxisText(const String text) {
		Bool ret = (itsOptionsYAxisText != text);
		itsOptionsYAxisText = text;
		itsOptionsYAxisTextUnset = false;
		return ret;
	}

	Bool WCAxisLabeller::unsetXAxisText() {
		Bool ret = (!itsOptionsXAxisTextUnset);
		itsOptionsXAxisTextUnset = true;
		return ret;
	}

	Bool WCAxisLabeller::unsetYAxisText() {
		Bool ret = (!itsOptionsYAxisTextUnset);
		itsOptionsYAxisTextUnset = true;
		return ret;
	}

	String WCAxisLabeller::xAxisText() const {
		if (!itsOptionsXAxisTextUnset) {
			return itsOptionsXAxisText;
		} else {
			return String("");
		}
	}

	String WCAxisLabeller::yAxisText() const {
		if (!itsOptionsYAxisTextUnset) {
			return itsOptionsYAxisText;
		} else {
			return String("");
		}
	}

	Bool WCAxisLabeller::setXAxisTextColor(const String color) {
		Bool ret = (itsOptionsXAxisTextColor != color);
		itsOptionsXAxisTextColor = color;
		return ret;
	}

	Bool WCAxisLabeller::setYAxisTextColor(const String color) {
		Bool ret = (itsOptionsYAxisTextColor != color);
		itsOptionsYAxisTextColor = color;
		return ret;
	}

	Bool WCAxisLabeller::setXGridType(const String type) {
		Bool ret = (itsOptionsXGridType != type);
		itsOptionsXGridType = type;
		return ret;
	}

	Bool WCAxisLabeller::setYGridType(const String type) {
		Bool ret = (itsOptionsYGridType != type);
		itsOptionsYGridType = type;
		return ret;
	}

	Bool WCAxisLabeller::setXGridColor(const String color) {
		Bool ret = (itsOptionsXGridColor != color);
		itsOptionsXGridColor = color;
		return ret;
	}

	Bool WCAxisLabeller::setYGridColor(const String color) {
		Bool ret = (itsOptionsYGridColor != color);
		itsOptionsYGridColor = color;
		return ret;
	}

	Bool WCAxisLabeller::setTickLength(const Float length) {
		Bool ret = (itsOptionsTickLength != length);
		itsOptionsTickLength = length;
		return ret;
	}

	Bool WCAxisLabeller::setLabelPosition(const String position) {
		Bool ret = (itsOptionsLabelPos != position);
		itsOptionsLabelPos = position;
		return ret;
	}

	Bool WCAxisLabeller::setPlotOutline(const Bool outline) {
		Bool ret = (itsOptionsPlotOutline != outline);
		itsOptionsPlotOutline = outline;
		return ret;
	}

	Bool WCAxisLabeller::setPlotOutlineColor(const String color) {
		Bool ret = (itsOptionsPlotOutlineColor != color);
		itsOptionsPlotOutlineColor = color;
		return ret;
	}

	Bool WCAxisLabeller::setCharSize(const Float size) {
		Bool ret = (itsOptionsCharSize != size);
		itsOptionsCharSize = size;
		return ret;
	}

	Bool WCAxisLabeller::setCharFont(const String font) {
		Bool ret = (itsOptionsCharFont != font);
		itsOptionsCharFont = font;
		return ret;
	}

	Bool WCAxisLabeller::setLineWidth(const Float width) {
		Bool ret = (itsOptionsLineWidth != width);
		itsOptionsLineWidth = width;
		return ret;
	}

} //# NAMESPACE CASA - END

