//# Copyright (C) 1998,1999,2000,2001,2003
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
//# $Id: $

#include <imageanalysis/ImageAnalysis/ImageConcatenator.h>

#include <casacore/casa/BasicSL/STLIO.h>
#include <casacore/casa/Exceptions/Error.h>
#include <casacore/images/Images/ImageConcat.h>
#include <imageanalysis/ImageAnalysis/ImageFactory.h>

namespace casa {

template <class T>
const casacore::String ImageConcatenator<T>::_class = "ImageConcatenator";
/*
template <class T>
ImageConcatenator<T>::ImageConcatenator(
	SPIIT image, const casacore::String& outname,
	casacore::Bool overwrite
) : ImageTask<T>(
		image, "", 0, "", "", "",
		"", outname, overwrite
	), _axis(-1), _tempClose(false), _relax(false), _reorder(false) {
	this->_construct();
}
*/

template <class T>
ImageConcatenator<T>::ImageConcatenator(
    std::vector<casacore::String>& imageNames, const casacore::String& outname,
	casacore::Bool overwrite
) : _imageNames(imageNames), _outname(_outname), _overwrite(overwrite) {
	// this->_construct();
    ThrowIf(
        _imageNames.size() < 2,
        "You must give at least two extant images to concatenate"
    );
	// SPCIIT myImage = this->_getImage();
    // vector<casacore::String> imageList;
	// imageList.push_back(myImage->name());
	// imageList.insert(
	//     imageList.end(), imageNames.begin(), imageNames.end()
	// );
    casacore::File p(_outname);
    ThrowIf(
        p.exists() && ! _overwrite,
        _outname + " exists and overwrite is false"
    );
    LogIO log;
    log << casacore::LogOrigin(_class, __func__, WHERE);
    log << casacore::LogIO::NORMAL << "Number of images to concatenate = "
        << (_imageNames.size()) << casacore::LogIO::POST;
}

template <class T>
ImageConcatenator<T>::~ImageConcatenator() {
    if (_mode == MOVEVIRTUAL) {
        ImageFactory::remove(this->_getImage(), false);
    }
}

template <class T>
void ImageConcatenator<T>::setAxis(int axis) {
	ThrowIf(
		axis >= (casacore::Int)this->_getImage()->ndim(),
		"Specified zero-based value of axis exceeds "
		"number of dimensions in image"
	);
	if (axis < 0) {
		ThrowIf(
			! this->_getImage()->coordinates().hasSpectralAxis(),
			"This image has no spectral axis"
		);
		_axis = this->_getImage()->coordinates().spectralAxisNumber(false);
	}
	else {
		_axis = axis;
	}
}

template <class T>
SPIIT ImageConcatenator<T>::concatenate() {
    ThrowIf(
        _outname.empty() && _mode != PAGED,
        "An empty outname can be used only if mode == PAGED"
    );
    if (_mode != PAGED) {
        casacore::Path t(_outname);
        for (const auto& name: _imageNames) {
            casacore::Path p(name);
            ThrowIf(
                p.absoluteName() == t.absoluteName(),
                "Cannot have one of the input images also be "
                "the output image is mode is not PAGED"
            );
        }
    }
    auto myImage = ImageFactory::fromFile(_imageNames[0], T(0), false);
	const auto& csys = myImage->coordinates();
	casacore::Int whichCoordinate, axisInCoordinate;
	csys.findPixelAxis(whichCoordinate, axisInCoordinate, _axis);
	const auto ctype = csys.coordinate(whichCoordinate).type();
	const auto pix = csys.referencePixel();
   	auto isIncreasing = false;
	vector<casacore::Double> minVals, maxVals;
	casacore::uInt n = 0;
	if (! _relax || _reorder) {
        // The image held by ImageTask is compared to
		n = _imageNames.size();
		minVals.resize(n);
		maxVals.resize(n);
		isIncreasing = _minMaxAxisValues(
			minVals[0], maxVals[0], myImage->ndim(),
			csys, myImage->shape()
		);
	}
    auto dataType = myImage->dataType();
    uInt i = 0;
    for(casacore::String name: _imageNames) {
        if (i == 0) {
            // The imageNames[0] is in myImage, so don't open it again
            continue;
        }
        auto myIm = ImageFactory::fromFile(name, T(0), false);
        auto oDType = myIm->dataType();
   		ThrowIf(
			oDType != dataType,
			"Concatenation of images of different data types is not supported"
		);
        if (! _relax || _reorder) {
            auto shape = myIm->shape();
            const auto& oCsys = myIm->coordinates();
			ThrowIf(
				_minMaxAxisValues(
					minVals[i], maxVals[i], shape.size(), oCsys, shape
				) != isIncreasing,
				"Coordinate axes in different images with opposing increment "
				"signs is not permitted if relax=false or reorder=true"
			);
		}
		if (! _relax) {
            const auto& oCsys = myIm->coordinates();
			oCsys.findPixelAxis(
				whichCoordinate, axisInCoordinate, _axis
			);
			ThrowIf(
				oCsys.coordinate(whichCoordinate).type() != ctype,
				"Cannot concatenate different coordinates in different images "
				"if relax=false"
			);
		}
		++i;
    }
	if (_reorder) {
		casacore::Sort sorter;
		sorter.sortKey(
			minVals.data(), TpDouble, 0,
			isIncreasing ? casacore::Sort::Ascending
			    : casacore::Sort::Descending
		);
		casacore::Vector<casacore::uInt> indices;
		sorter.sort(indices, n);
		auto tmp = _imageNames;
		auto iter = tmp.begin();
		auto end = tmp.end();
		auto index = indices.begin();
		while (iter != end) {
			*iter++ = _imageNames[*index++];
		}
		_imageNames = tmp;
        LogIO log;
		log << LogOrigin("ImageConcatenator", __func__)
            << casacore::LogIO::NORMAL
			<< "Images will be concatenated in the order "
			<< _imageNames << " and the coordinate system of "
			<< _imageNames[0] << " will be used as the reference"
			<< casacore::LogIO::POST;
	}
	std::shared_ptr<casacore::ImageConcat<T>> pConcat(
	    new casacore::ImageConcat<T>(_axis, _tempClose)
	);
	ThrowIf(
		! pConcat.get(), "Failed to create ImageConcat object"
	);
    // auto copyNames = imageList;
    if (_mode == COPYVIRTUAL || _mode == MOVEVIRTUAL) {
        casacore::File p(_outname);
        casacore::Directory eldir(p);
        eldir.create(_overwrite);
        auto rootdir = eldir.path().absoluteName();
        /*
        uInt k = 0;
        for (auto& name: copyNames) {
            name = rootdir + Path(imageNames[k]).baseName();
            ++k;
        }
        */
        // if mode == MOVE the first image gets removed during object destruction.
        /*
        auto inname = myImage->name();
        Directory myFirst(inname);
        ThrowIf(
            myFirst.copy(outname),
            "Error copying image " + inname
        );
        */
        casacore::uInt k = 0;
        auto copyNames = _imageNames;
        // FIXME probably needs to be imageList since imageList can get screwed around above
        for (auto imname: _imageNames) {
            casacore::Directory elim(imname);
            copyNames[k] = rootdir + "/" + elim.path().baseName();
            if (_mode == MOVEVIRTUAL) {
                elim.move(copyNames[k]);
            }
            else if (_mode == COPYVIRTUAL) {
                elim.copy(copyNames[k]);
            }
            else {
                ThrowCc("Logic Error");
            }
            ++k;
        }
        //	    imageList.insert(imageList.begin(), copyNames.begin(), copyNames.end());
        _imageNames = copyNames;
    }
	if (_axis < 0) {
        setAxis(-1);
	}
    /*
	*this->_getLog() << casacore::LogOrigin(_class, __func__, WHERE);
	*this->_getLog() << casacore::LogIO::NORMAL << "Number of images to concatenate = "
		<< (imageNames.size() + 1) << casacore::LogIO::POST;
	const auto& csys = myImage->coordinates();
	casacore::Int whichCoordinate, axisInCoordinate;
	csys.findPixelAxis(whichCoordinate, axisInCoordinate, _axis);
	casacore::Coordinate::Type ctype = csys.coordinate(whichCoordinate).type();
	casacore::Vector<casacore::Double> pix = csys.referencePixel();
    vector<casacore::String> imageList;
	imageList.push_back(myImage->name());
	imageList.insert(
	    imageList.end(), copyNames.begin(), copyNames.end()
	);
	casacore::Bool isIncreasing = false;
	vector<casacore::Double> minVals, maxVals;
	casacore::uInt n = 0;
	if (! _relax || _reorder) {
		n = imageList.size();
		minVals.resize(n);
		maxVals.resize(n);
		isIncreasing = _minMaxAxisValues(
			minVals[0], maxVals[0], myImage->ndim(),
			csys, myImage->shape()
		);
	}
    */
	//casacore::uInt i = 1;
	//for(casacore::String name: copyNames) {
		//auto myIm = ImageFactory::fromFile(name, T(0), false);
        //oDType = myIm->dataType();
        // oCsys = myIm->coordinates();
        // shape = myIm->shape();
        /*
        SPCIIF imageF;
        SPCIIC imageC;
        SPCIID imageD;
        SPCIIDC imageDC;
        std::tie(imageF, imageC, imageD, imageDC) = imagePtrs;
		casacore::DataType oDType;
		casacore::CoordinateSystem oCsys;
		casacore::IPosition shape;
		if (imageF) {
		    oDType = imageF->dataType();
		    oCsys = imageF->coordinates();
		    shape = imageF->shape();
		}
		else if (imageC) {
		    oDType = imageC->dataType();
		    oCsys = imageC->coordinates();
		    shape = imageC->shape();
		}
		else if (imageD) {
		    oDType = imageD->dataType();
		    oCsys = imageD->coordinates();
		    shape = imageD->shape();
		}
		else if (imageDC) {
		    oDType = imageDC->dataType();
		    oCsys = imageDC->coordinates();
		    shape = imageDC->shape();
		}
		else {
		    ThrowCc("Logic error");
		}
		ThrowIf(
			oDType != dataType,
			"Concatenation of images of different data types is not supported"
		);
        if (! _relax || _reorder) {
			ThrowIf(
				_minMaxAxisValues(
					minVals[i], maxVals[i], shape.size(), oCsys, shape
				) != isIncreasing,
				"Coordinate axes in different images with opposing increment "
				"signs is not permitted if relax=false or reorder=true"
			);
		}
		if (! _relax) {
			oCsys.findPixelAxis(
				whichCoordinate, axisInCoordinate, _axis
			);
			ThrowIf(
				oCsys.coordinate(whichCoordinate).type() != ctype,
				"Cannot concatenate different coordinates in different images "
				"if relax=false"
			);
		}
		++i;
        */
	// }
    auto first = true;
	for(const auto& name: _imageNames) {
		_addImage(pConcat, name, first);
        first = false;
	}
    if (_mode == PAGED) {
        // return this->_prepareOutputImage(*pConcat);
        static const casacore::Record empty;
        static const casacore::String emptyString;
        return SubImageFactory<T>::createImage(
            *pConcat, _outname, empty, emptyString, false,
            _overwrite, true, false, false
        );  
    }
    else {
        return SPIIT(ImageFactory::fromFile(_outname), T(0), false);
    }
    /*
    switch (_mode) i{
        case PAGED:
	    return this->_prepareOutputImage(*pConcat);
      default: {
          cout << __FILE__ << " " << __LINE__ << endl;
          auto outname = this->_getOutname();
          cout << __FILE__ << " " << __LINE__ << endl;
          if (outname.empty()) {
          cout << __FILE__ << " " << __LINE__ << endl;
              return pConcat;
          }
          else {
          cout << __FILE__ << " " << __LINE__ << endl;
            pConcat->save(outname);
          cout << __FILE__ << " " << __LINE__ << endl;
              return ImageFactory::fromFile(outname, T(0), False);
          }
        }
    }
    */
}

template <class T> void ImageConcatenator<T>::setMode(const casacore::String& mymode) {
    auto m = mymode;
    m.downcase();
    if (m.startsWith("m")) {
        _mode = MOVEVIRTUAL;
    }
    else if (m.startsWith("c")) {
        _mode = COPYVIRTUAL;
    }
    else if (m.startsWith("n")) {
        _mode = NOMOVEVIRTUAL;
    }
    else if (m.startsWith("p")) {
        _mode = PAGED;
    }
    else {
        ThrowCc("Unsupported mode " + mymode);
    }
}

template <class T> void ImageConcatenator<T>::_addImage(
	std::shared_ptr<casacore::ImageConcat<T>> pConcat, const casacore::String& name,
    casacore::Bool first
) const {
    cout << this->_getImage()->name() << " " << name << endl;
    if (first) {
        /*
        auto outname = _mode == PAGED ? "" : name;
        SPIIT mycopy = SubImageFactory<T>::createImage(
			*this->_getImage(), outname, casacore::Record(), "", false, false, false, false
		);
        cout << "added temporary image copy" << endl;
		pConcat->setImage(*mycopy, _relax);
        */
        pConcat->setImage(*ImageFactory::fromFile(name, T(0), false), _relax);
		return;
	}
	casacore::Bool doneOpen = false;
	try {
		SPIIT im2 = casacore::ImageUtilities::openImage<T>(name);
        cout << im2->imageType() << endl;
        doneOpen = true;
		pConcat->setImage(*im2, _relax);
	}
	catch (const casacore::AipsError& x) {
		ThrowIf(doneOpen, x.getMesg());
		ThrowCc(
			"Failed to open file " + name
			+ "This may mean you have too many files open simultaneously. "
			"Try using tempclose=T in the imageconcat constructor. "
			"Exception message " + x.getMesg()
		);
	}
}

template <class T> casacore::Bool ImageConcatenator<T>::_minMaxAxisValues(
	casacore::Double& mymin, casacore::Double& mymax, casacore::uInt ndim, const casacore::CoordinateSystem& csys,
	const casacore::IPosition& shape
) const {
	ThrowIf(
		ndim != this->_getImage()->ndim(),
		"All images must have the same number of dimensions"
	);
	casacore::Vector<casacore::Double> pix = csys.referencePixel();
	pix[_axis] = 0;
	mymin = csys.toWorld(pix)[_axis];
	if (shape[_axis] == 1) {
		mymax = mymin;
		return this->_getImage()->coordinates().increment()[_axis] > 0;
	}
	pix[_axis] = shape[_axis] - 1;
	mymax = csys.toWorld(pix)[_axis];
	casacore::Bool isIncreasing = mymax > mymin;
	if (! isIncreasing) {
        std::swap(mymin, mymax);
	}
	return isIncreasing;
}

template <class T>
String ImageConcatenator<T>::getClass() const {
	return _class;
}

}


