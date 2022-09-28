//# tPowerLogPoly.cc: tests the PowerLogPoly class
//# Copyright (C) 1998,1999,2000,2001,2003
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

#include <casacore/casa/Arrays/ArrayLogical.h>
#include <casa/aips.h>
#include <components/ComponentModels/ComponentType.h>
#include <components/ComponentModels/PowerLogPoly.h>
#include <casa/Arrays/Vector.h>
#include <casa/Containers/Record.h>
#include <casa/Containers/RecordFieldId.h>
#include <casa/Exceptions/Error.h>
#include <casa/BasicMath/Math.h>
#include <measures/Measures/MFrequency.h>
#include <measures/Measures/MeasureHolder.h>
#include <measures/Measures/MeasRef.h>
#include <casa/Quanta/MVFrequency.h>
#include <casa/Quanta/Quantum.h>
#include <casa/Utilities/Assert.h>
#include <casa/BasicSL/String.h>
#include <casa/iostream.h>

#include <casa/namespace.h>

#include <memory>

using namespace casa;
int main() {
    try {
        std::unique_ptr<SpectralModel> siPtr;
        const MFrequency f1(Quantity(1.0, "GHz"), MFrequency::LSRK);
        const MFrequency f2(Quantity(2.0, "GHz"), MFrequency::LSRK);
        const MFrequency f4(Quantity(4.0, "GHz"), MFrequency::LSRK);
        {
            PowerLogPoly plp;
            AlwaysAssert(plp.ok(), AipsError);
            AlwaysAssert(
                plp.type() == ComponentType::PLP, AipsError
            );
            AlwaysAssert(
                plp.nParameters() == 1 && plp.parameters()[0] == 0,
                AipsError
            );
            AlwaysAssert(
                plp.refFrequency().get("GHz").getValue() == 1,
                AipsError
            );
            AlwaysAssert(
                plp.refFrequency().get("MHz").getValue() == 1000,
                AipsError
            );
            AlwaysAssert(near(plp.sample(f1), 1.0), AipsError);
            AlwaysAssert(near(plp.sample(f2), 1.0), AipsError);
            plp.setRefFrequency(f2);
            siPtr.reset(plp.clone());
            cout << "Passed the default constructor test" << endl;
        }
        {
            Vector<Double> indices(2);
            indices[0] = 2;
            indices[1] = 1;
            siPtr->setParameters(indices);
            // show that its a real copy not a reference
            indices[0] = 1000;
            indices[1] = 1000;
            auto p = siPtr->parameters();
            AlwaysAssert((p[0] == 2 && p[1] == 1), AipsError);
            siPtr->setRefFrequency(f2);
            AlwaysAssert(
                siPtr->refFrequency().get("GHz").getValue() == 2.0, 
		        AipsError
            );
            AlwaysAssert(siPtr->nParameters() == 2, AipsError);
            PowerLogPoly copy(*((PowerLogPoly*)(siPtr.get())));
            PowerLogPoly assigned;
            assigned = copy;
            indices.resize(3);
            indices[0] = 3;
            indices[1] = 4;
            indices[2] = 5;
            copy.setParameters(indices);
            copy.setRefFrequency(f1);
            assigned.setRefFrequency(f4);
            indices.resize(2);
            indices[0] = 10;
            indices[1] = 11;
            assigned.setParameters(indices);
            AlwaysAssert(copy.nParameters() == 3, AipsError);
            auto k = copy.parameters();
            AlwaysAssert(
                k.size() == 3 && k[0] == 3 && k[1] == 4 && k[2] == 5,
                AipsError
            );
            AlwaysAssert(
                siPtr->refFrequency().get("GHz").getValue() == 2, 
		        AipsError
            );
            AlwaysAssert(
                copy.refFrequency().get("GHz").getValue() == 1, 
		        AipsError
            );
            AlwaysAssert(
                assigned.refFrequency().get("GHz").getValue() == 4.0, 
		        AipsError
            );
            auto x = assigned.parameters();
            AlwaysAssert(assigned.nParameters() == 2, AipsError);
            AlwaysAssert(
                x.size() == 2 && x[0] == 10 && x[1] == 11, AipsError
            );
            cout << "Passed the copy semantics test" << endl;
        }
        {

            Vector<Double> indices(4);
            indices[0] = 2;
            indices[1] = 1.5;
            indices[2] = -0.2;
            indices[3] = -1.8;
            PowerLogPoly plp(f2, indices);
            AlwaysAssert(
                plp.refFrequency().get("GHz").getValue() == 2.0, 
		        AipsError
            );
            auto r1 = plp.sample(f1);
            auto r2 = plp.sample(f2);
            auto r4 = plp.sample(f4);
            AlwaysAssert(near(r1, 0.362578931147, 1e-12), AipsError);
            AlwaysAssert(r2 == 1, AipsError);
            AlwaysAssert(near(r4, 5.077737427813, 1e-12), AipsError);
            Vector<MVFrequency> freqs(3);
            freqs(0) = f1.getValue(); 
            freqs(1) = f2.getValue();
            freqs(2) = f4.getValue();
            MFrequency::Ref ref(MFrequency::LSRK);
            Vector<Double> results(3);
            plp.sample(results, freqs, ref);
            AlwaysAssert(near(results(0), r1, 1e-10), AipsError);
            AlwaysAssert(near(results(1), 1.0), AipsError);
            AlwaysAssert(near(results(2), r4, 1e-10), AipsError);
            Vector<Double> stokes(4, 1);
            stokes[0] = 5;
            auto expec = stokes.copy();
            expec[0] = stokes[0]*r4;
            plp.sampleStokes(f4, stokes);
            AlwaysAssert(allNear(stokes, expec, 1e-10), AipsError);
            Matrix<Double> vstokes(3,4);
            vstokes(0, 0) = 5;
            vstokes(1, 0) = 6;
            vstokes(2, 0) = 7;
            Matrix<Double> vexp = vstokes.copy();
            vexp(0, 0) = r1 * vstokes(0, 0);
            vexp(1, 0) = r2 * vstokes(1, 0);
            vexp(2, 0) = r4 * vstokes(2, 0);
            plp.sampleStokes(vstokes, freqs, ref);
            cout << vstokes << endl;
            cout << vexp << endl;
            AlwaysAssert(allNear(vstokes, vexp, 1e-10), AipsError);
            cout << "Passed the multi-sample test" << endl;
            Record rec;
            String errMsg;
            AlwaysAssert(plp.toRecord(errMsg, rec), AipsError); 
            AlwaysAssert(errMsg == "", AipsError);
            AlwaysAssert(rec.isDefined("type"), AipsError);
            AlwaysAssert(rec.isDefined("frequency"), AipsError);
            AlwaysAssert(rec.isDefined("coeffs"), AipsError);
            AlwaysAssert(rec.isDefined("error"), AipsError);
            String type;
            rec.get(RecordFieldId("type"), type);
            AlwaysAssert(type == "Power Logarithmic Polynomial", AipsError);
            Vector<Double> coeffs;
            rec.get(RecordFieldId("coeffs"), coeffs);
            AlwaysAssert(allNear(coeffs, indices, 1e-10), AipsError); 
            Record freqRec = rec.asRecord(RecordFieldId("frequency"));
            MeasureHolder mh;
            mh.fromRecord(errMsg, freqRec);
            AlwaysAssert(errMsg.length() == 0, AipsError);
            AlwaysAssert(mh.isMFrequency(), AipsError);
            mh = f1;
            Record newRec;
            newRec.define(RecordFieldId("type"), "poWer LoGarithmic poLynomial");
            newRec.define(RecordFieldId("coeffs"), Vector<Double>(2, 1));
            Record newFreq;
            AlwaysAssert(mh.toRecord(errMsg, newFreq), AipsError);
            AlwaysAssert(errMsg.length() == 0, AipsError);
            newRec.defineRecord(RecordFieldId("frequency"), newFreq);
            AlwaysAssert(plp.fromRecord(errMsg, newRec), AipsError);
            AlwaysAssert(
                allNear(plp.parameters(), Vector<Double>(2, 1), 1e-10),
                AipsError
            );
            AlwaysAssert(
                near(plp.refFrequency().get("GHz").getValue(), 1.0), 
		        AipsError
            );
            Record emptyRec;
            AlwaysAssert(plp.convertUnit(errMsg, emptyRec), AipsError);
            emptyRec.define(RecordFieldId("index"), "");
            AlwaysAssert(plp.convertUnit(errMsg, emptyRec), AipsError);
            emptyRec.define(RecordFieldId("index"), "deg");
            AlwaysAssert(plp.convertUnit(errMsg, emptyRec), AipsError);
            PowerLogPoly plp2(f4, indices);
            plp2.toRecord(errMsg, rec);
            PowerLogPoly plp3;
            plp3.fromRecord(errMsg, rec);
            AlwaysAssert(
                allNear(plp3.parameters(), plp2.parameters(), 1e-10),
                AipsError
            );
            cout << "Passed the record handling test" << endl;
        }
    }
    catch (const AipsError& x) {
        cerr << x.getMesg() << endl;
        cout << "FAIL" << endl;
        return 1;
    }
    catch (...) {
        cerr << "Exception not derived from AipsError" << endl;
        cout << "FAIL" << endl;
        return 2;
    }
    cout << "OK" << endl;
    return 0;
}
// compile-command: "gmake OPTLIB=1 tSpectralIndex"
