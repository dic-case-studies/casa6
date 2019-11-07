/*
 * SimpleSIImageStore is an implementation without much checking the caller must make sure pointers exist before using them
 * Copyright (C) 2019  Associated Universities, Inc. Washington DC, USA.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * 
 * Queries concerning CASA should be submitted at
 *        https://help.nrao.edu
 *
 *         Postal address: CASA Project Manager 
 *         National Radio Astronomy Observatory
 *         520 Edgemont Road
 *         Charlottesville, VA 22903-2475 USA
 */
#include <casa/Exceptions/Error.h>
#include <casa/iostream.h>
#include <synthesis/ImagerObjects/SimpleSIImageStore.h>

using namespace std;

using namespace casacore;

namespace casa { //# NAMESPACE CASA - BEGIN
	
	SimpleSIImageStore::SimpleSIImageStore  (const shared_ptr<ImageInterface<Float> > &modelim,
	       const shared_ptr<ImageInterface<Float> > &residim, const shared_ptr<ImageInterface<Float> > &psfim,
	       const shared_ptr<ImageInterface<Float> > &weightim,const shared_ptr<ImageInterface<Float> > &restoredim,
	       const shared_ptr<ImageInterface<Float> > &maskim,const shared_ptr<ImageInterface<Float> > &sumwtim,
	       const shared_ptr<ImageInterface<Float> > &gridwtim, const shared_ptr<ImageInterface<Float> > &pbim,
	       const shared_ptr<ImageInterface<Float> > &restoredpbcorim,const Bool useweightimage) : SIImageStore() {
			if(!psfim && !residim)  {
				throw(AipsError("SimpleSIImagestore has to have a valid residual or psf image"));
			}
			else{
				shared_ptr<ImageInterface<Float> > theim= psfim ? psfim : residim;
				itsCoordSys=theim->coordinates();
				itsImageShape=theim->shape();
			}
			if(useweightimage && !weightim)
				throw(AipsError("SimpleSIImagestore has to have a valid weightimage for this kind of weighting scheme"));
			itsPsf = psfim;
			itsModel=modelim;
			itsResidual=residim;
			itsWeight=weightim;
			itsImage=restoredim;
			itsSumWt=sumwtim;
			itsMask=maskim;
			itsImagePBcor=restoredpbcorim;
			
			itsUseWeight=useweightimage;
			
		
		
		
	}
	shared_ptr<ImageInterface<Float> > SimpleSIImageStore::psf(uInt){
		if(!itsPsf)
			throw(AipsError("Programmer's error: calling for psf without setting it"));
		return itsPsf;
		
	}
	shared_ptr<ImageInterface<Float> > SimpleSIImageStore::residual(uInt){
		if(!itsResidual)
			throw(AipsError("Programmer's error: calling for residual without setting it"));
		return itsResidual;
		
	}
	shared_ptr<ImageInterface<Float> > SimpleSIImageStore::weight(uInt){
		if(!itsWeight)
			throw(AipsError("Programmer's error: calling for weight without setting it"));
		return itsWeight;
		
	}
	shared_ptr<ImageInterface<Float> > SimpleSIImageStore::model(uInt){
		if(!itsModel)
			throw(AipsError("Programmer's error: calling for model without setting it"));
		return itsModel;
		
	}
	shared_ptr<ImageInterface<Float> > SimpleSIImageStore::sumwt(uInt){
		if(!itsSumWt)
			throw(AipsError("Programmer's error: calling for sumweight without setting it"));
		return itsSumWt;
		
	}
	
	shared_ptr<ImageInterface<Complex> > SimpleSIImageStore::forwardGrid(uInt){
		if ( !itsForwardGrid ) {
			Vector<Int> whichStokes ( 0 );
			IPosition cimageShape;
			cimageShape=itsImageShape;
			MFrequency::Types freqframe = itsCoordSys.spectralCoordinate ( itsCoordSys.findCoordinate ( Coordinate::SPECTRAL ) ).frequencySystem ( True );
			// No need to set a conversion layer if image is in LSRK already or it is 'Undefined'
			if ( freqframe != MFrequency::LSRK && freqframe!=MFrequency::Undefined && freqframe!=MFrequency::REST ) {
				itsCoordSys.setSpectralConversion ( "LSRK" );
			}
			CoordinateSystem cimageCoord = StokesImageUtil::CStokesCoord ( itsCoordSys,
                                       whichStokes, itsDataPolRep );
			cimageShape ( 2 ) =whichStokes.nelements();

			//cout << "Making forward grid of shape : " << cimageShape << " for imshape : " << itsImageShape << endl;
			itsForwardGrid.reset ( new TempImage<Complex> ( TiledShape ( cimageShape, tileShape() ), cimageCoord, memoryBeforeLattice() ) );

		}
		return itsForwardGrid;
		
	}
	shared_ptr<ImageInterface<Complex> > SimpleSIImageStore::backwardGrid(uInt){
		if(!itsBackwardGrid){
			Vector<Int> whichStokes(0);
			IPosition cimageShape;
			cimageShape=itsImageShape;
			MFrequency::Types freqframe = itsCoordSys.spectralCoordinate(itsCoordSys.findCoordinate(Coordinate::SPECTRAL)).frequencySystem(True);
			// No need to set a conversion layer if image is in LSRK already or it is 'Undefined'
			if(freqframe != MFrequency::LSRK && freqframe!=MFrequency::Undefined && freqframe!=MFrequency::REST ) 
			{ itsCoordSys.setSpectralConversion("LSRK"); }
			CoordinateSystem cimageCoord = StokesImageUtil::CStokesCoord( itsCoordSys,
								  whichStokes, itsDataPolRep);
			cimageShape(2)=whichStokes.nelements();
			//cout << "Making backward grid of shape : " << cimageShape << " for imshape : " << itsImageShape << endl;
			itsBackwardGrid.reset( new TempImage<Complex>(TiledShape(cimageShape, tileShape()), cimageCoord, memoryBeforeLattice()) );
			
		}
		return itsBackwardGrid;
		
	}
	
	
	
	
} //# NAMESPACE CASA - END
