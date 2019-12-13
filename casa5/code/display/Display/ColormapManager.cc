//# ColormapManager.cc: management of dynamic color allocation
//# Copyright (C) 1994,1995,1996,1997,1998,1999,2000,2002
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

#include <display/Display/ColormapManager.h>
#include <display/Display/ColormapInfo.h>
#include <display/Display/PixelCanvasColorTable.h>
#include <casa/iostream.h>

using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN

	ColormapManager::ColormapManager(PixelCanvasColorTable * pcctbl) :
		itsPCColorTable(pcctbl) {
	}

	ColormapManager::~ColormapManager() {
	}

	void ColormapManager::registerColormap(Colormap *cmap, Float weight) {
		if ( cmap ){
			if (itsInfoMap.find(cmap) != itsInfoMap.end( )) {
				// already known
				ColormapInfo * mi = itsInfoMap[cmap];
				mi->ref();
				// compare weights -- if different call redistribute
				if (mi->weight() != weight)	{
					mi->setWeight(weight);
					redistributeColormaps();
				}
			}
			else {
				// new colormap
				ColormapInfo *mi = new ColormapInfo(cmap, weight, 0, 0);
				mi->ref();
				itsInfoMap[cmap] = mi;
				redistributeColormaps();
			}
			// now tell the cmap that it is used by the PixelCanvasColorTable
			// that this is a ColormapManager for...
			cmap->registerPCColorTable(itsPCColorTable);
		}
	}

	void ColormapManager::registerColormap(Colormap *cmap,
	                                       Colormap *cmapToReplace) {
		if (cmap == cmapToReplace) {
			return;
		}
		if (itsInfoMap.find(cmap) != itsInfoMap.end( )) {
			// already defined, so add to ref count, and decrement ref count
			// of cmapToReplace by attempting to unregister it
			ColormapInfo *mi = itsInfoMap[cmap];
			mi->ref();
			cmap->registerPCColorTable(itsPCColorTable);
			if (itsInfoMap.find(cmapToReplace) != itsInfoMap.end( )) {
				unregisterColormap(cmapToReplace);
			}
			return;
		} else if (itsInfoMap.find(cmapToReplace) == itsInfoMap.end( )) {
			// but cmapToReplace is not defined, so just register in the usual
			// way!
			registerColormap(cmap);
			return;
		} else {
			// ok, let's see if we can replace cmapToReplace:
            auto iter = itsInfoMap.find(cmapToReplace);
            if ( iter != itsInfoMap.end( ) ) {
                ColormapInfo *mi = iter->second;
                itsInfoMap.erase(iter);
                mi->unref(); // decrement ref of cmapToReplace no matter what.
                cmapToReplace->unregisterPCColorTable(itsPCColorTable);
                if (mi->refCount() != 0) {
                    // we cannot, so just register in the usual way.
                    registerColormap(cmap);
                    return;
                }
                // okilidokile, on with the substitution
                ColormapInfo *minew = new ColormapInfo(cmap, 1.0, 0, 0);
                minew->ref();
                itsInfoMap[cmap] = minew;
                minew->setWeight(mi->weight());
                minew->setOffset(mi->offset());
                minew->setSize(mi->size());
                cmapToReplace->unregisterPCColorTable(itsPCColorTable);
                cmap->registerPCColorTable(itsPCColorTable);
                reinstallColormaps();
                delete mi;
            }
		}
		return;
	}

	Bool ColormapManager::unregisterColormap(Colormap *cmap) {
        auto miptr = itsInfoMap.find(cmap);
		if (miptr != itsInfoMap.end( )) {
			ColormapInfo * mi = miptr->second;
			mi->unref();
			if (mi->refCount() == 0) {
				itsInfoMap.erase(miptr);
				AlwaysAssert(itsInfoMap.find(cmap) == itsInfoMap.end( ), AipsError);
				delete mi;
				redistributeColormaps();
			}
			// now tell the cmap that it is no longer used by the
			// PixelCanvasColorTable for which this is a ColormapManager ...
			cmap->unregisterPCColorTable(itsPCColorTable);
			return true;
		} else {
			cerr << "Colormap: " << *cmap << endl;
			throw(AipsError("unregisterColormap passed unknown colormap"));
		}
		return false;
	}

	uInt ColormapManager::getColormapSize(const Colormap *cmap) const {
		uInt retval = 0;
		auto miptr = itsInfoMap.find(cmap);
		if (miptr != itsInfoMap.end( )) {
			ColormapInfo * mi = miptr->second;
			retval = mi->size();
		} else {
			cerr << "Colormap: " << *cmap << endl;
			throw(AipsError("getColormapSize passed unknown colormap"));
		}
		return retval;
	}

	uInt ColormapManager::getColormapOffset(const Colormap *cmap) const {
		uInt retval = 0;
		auto miptr = itsInfoMap.find(cmap);
		if (miptr != itsInfoMap.end( )) {
			ColormapInfo * mi = miptr->second;
			retval = mi->offset();
		} else {
			cerr << "Colormap: " << *cmap << endl;
			throw(AipsError("getColormapOffset passed unknown colormap"));
		}
		return retval;
	}

	Bool ColormapManager::member(const Colormap *cmap) const {
		// needs to return base_pixel plus dcmap's offset
		return (itsInfoMap.find(cmap) != itsInfoMap.end( )) ? true : false;
	}

	const Colormap *ColormapManager::getMap(const uInt mapnum) const {
		if (mapnum < nMaps()) {
			uInt count = 0;
			for ( auto iter = itsInfoMap.begin( ); iter != itsInfoMap.end( ); ++iter, ++count )
				if ( count == mapnum ) return iter->first;
		}
		return (Colormap *)NULL;
	}

	void ColormapManager::redistributeColormaps() {
		// add up total frozen space, subtract from size
		// if negative throw exception and bail.
		// if positive, distribute remaining color cells according
		// to weight.

		uInt nRigid = 0;
		uInt nFlexible = 0;
		uInt n = itsInfoMap.size( );
		uInt size = itsPCColorTable->nColors();
		uInt p = 0;
		Float totalWeight = 0;

		if (n == 0) {
			// use default colormap
			return;
		}

		// compute statistics for controlled colormaps
		for (auto iter = itsInfoMap.begin( ); iter != itsInfoMap.end( ); ++iter) {
			ColormapInfo * mi = iter->second;
			if (mi->colormap()->rigid()) {
				nRigid++;
				p += mi->colormap()->rigidSize();
			} else {
				nFlexible++;
				totalWeight += mi->weight();
			}
		}

		// Bail if we can't fix it
		if (p + 2*nFlexible > size) {
			cerr << "# fixed maps: " << nRigid << ", total colors = " << p << endl;
			cerr << "# var.  maps: " << nFlexible << ", total req colors = "
			     << 2*nFlexible << endl;
			cerr << "        size: " << size << endl;
			cerr.flush();
			throw(AipsError("ColormapManager unable to distribute maps"));
			return;
		}

		// set distribution for for fixed maps
		p = 0;
		for (auto iter = itsInfoMap.begin( ); iter != itsInfoMap.end( ); ++iter) {
			ColormapInfo * mi = iter->second;
			if (mi->colormap()->rigid()) {
				mi->setOffset(p);
				mi->setSize(mi->colormap()->rigidSize());
				p += mi->colormap()->rigidSize();
			}
		}
		uInt spaceLeft = size - p;

		// Compute distribution for remaining flexible data colormaps
		if (nFlexible == 1) {
			for (auto iter = itsInfoMap.begin( ); iter != itsInfoMap.end( ); ++iter) {
				ColormapInfo * mi = iter->second;
				if (!mi->colormap()->rigid()) {
					mi->setOffset(p);
					mi->setSize(spaceLeft);
				}
			}
		} else {
			uInt count = 0;
            for (auto iter = itsInfoMap.begin( ); iter != itsInfoMap.end( ); ++iter) {
				ColormapInfo * mi = iter->second;
				if (!mi->colormap()->rigid()) {
					count++;
					uInt share;
					if (count < nFlexible) {
						share = (uInt) (0.5 + spaceLeft * mi->weight() / totalWeight);
						if (share < 2 && mi->size() > 1) {
							share = 2;
						}
					} else {
						share = (uInt) size - p;
					}
					mi->setSize(share);
					mi->setOffset(p);
					p += share;
				}
			}
		}
		itsPCColorTable->doResizeCallbacks(Display::ClearPriorToColorChange);
		reinstallColormaps();
		itsPCColorTable->doResizeCallbacks();
	}

	void ColormapManager::reinstallColormaps() {
		Vector<Float> redMap(1u,0.0);
		Vector<Float> greenMap(1u,0.0);
		Vector<Float> blueMap(1u,0.0);
		Vector<Float> alphaMap(1u,0.0);

        for (auto iter = itsInfoMap.begin( ); iter != itsInfoMap.end( ); ++iter) {
			ColormapInfo * mi = iter->second;
			redMap.resize(mi->size());
			greenMap.resize(mi->size());
			blueMap.resize(mi->size());
			mi->colormap()->calcRGBMaps(mi->size(), redMap, greenMap, blueMap, alphaMap );
			itsPCColorTable->installRGBColors(redMap, greenMap, blueMap, alphaMap, mi->offset());
		}
		itsPCColorTable->doResizeCallbacks(Display::ColormapChange);
	}

	ostream & operator << (ostream & os, const ColormapManager & cm) {
		cout << "-------------------- Colormap Manager ----------------------\n";
		uInt nMaps = cm.itsInfoMap.size();
		cout << "    Range of values: 0 to " << cm.itsPCColorTable->nColors()-1
		     << endl;
		cout << "Number of Colormaps: " << nMaps << endl;
        uInt count = 0;
		for (auto iter = cm.itsInfoMap.begin( ); iter != cm.itsInfoMap.end( ); ++iter, ++count) {
			ColormapInfo * mi = iter->second;
			cout << "Map # " << count+1 << " :" << *mi->colormap() << " Occupies ["
			     << mi->offset() << ","
			     << mi->size() + mi->offset() - 1 << "]";
			cout << "." << endl;
		}
		cout << "------------------- END ColormapManager --------------------\n";
		return os;
	}


} //# NAMESPACE CASA - END

