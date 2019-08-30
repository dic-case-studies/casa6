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

#ifndef ANNOTATIONS_ANNROTBOX_H
#define ANNOTATIONS_ANNROTBOX_H

#include <casa/aips.h>
#include <imageanalysis/Annotations/AnnPolygon.h>

namespace casa {

// <summary>
// This class represents an annotation for rectangular (in position coordinates) region specified
// in an ascii region file as proposed in CAS-2285. It is specified by its center position
// and side widths and a position angle.
// </summary>
// <author>Dave Mehringer</author>
// <use visibility=export>
// <reviewed reviewer="" date="yyyy/mm/dd" tests="" demos="">
// </reviewed>
// <prerequisite>

// </prerequisite>

// <etymology>
// Holds the specification of an annotation for a rectangular region as specified in ASCII format.
// Specified by center position and widths of sides and a position angle
// </etymology>

// <synopsis>
// This class represents an annotation for a rectangular region in coordinate space specified by
// center and widths of sides and a position angle.
// </synopsis>


class AnnRotBox: public AnnPolygon {

public:

	// <src>positionAngle</src> is measured in the usual astronomical
	// way; starting at north through east (counterclockwise)
	AnnRotBox(
		const casacore::Quantity& xcenter,
		const casacore::Quantity& ycenter,
		const casacore::Quantity& xwidth,
		const casacore::Quantity& ywidth, const casacore::Quantity& positionAngle,
		const casacore::String& dirRefFrameString,
		const casacore::CoordinateSystem& csys,
		const casacore::IPosition& imShape,
		const casacore::Quantity& beginFreq,
		const casacore::Quantity& endFreq,
		const casacore::String& freqRefFrameString,
		const casacore::String& dopplerString,
		const casacore::Quantity& restfreq,
		const casacore::Vector<casacore::Stokes::StokesTypes> stokes,
		const casacore::Bool annotationOnly,
		const casacore::Bool requireImageRegion=true
	);

	// Simplified constructor.
	// all frequencies are used (these can be set after construction).
	// xcenter and ycenter
	// must be in the same frame as the csys direction coordinate.
	// is a region (not just an annotation), although this value can be changed after
	// construction.
	AnnRotBox(
		const casacore::Quantity& xcenter,
		const casacore::Quantity& ycenter,
		const casacore::Quantity& xwidth,
		const casacore::Quantity& ywidth, const casacore::Quantity& positionAngle,
		const casacore::CoordinateSystem& csys,
		const casacore::IPosition& imShape,
		const casacore::Vector<casacore::Stokes::StokesTypes>& stokes,
		const casacore::Bool requireImageRegion=true
	);

	// implicit copy constructor and destructor are fine

	AnnRotBox& operator=(const AnnRotBox& other);

	virtual std::ostream& print(std::ostream &os) const;

private:
	AnnotationBase::Direction _inputCenter;
	casacore::Vector<casacore::Quantity> _inputWidths;
	casacore::Quantity _positionAngle;

};

}

#endif
