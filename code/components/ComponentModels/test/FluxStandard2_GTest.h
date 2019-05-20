//# FluxStandard2_GTest.h: definition of FluxStandard2 (tests newer standards with coefs stored external to the code )  google test            
//#                                                                            
//# Copyright (C) 2015                                                         
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
//# 51 Franklin Street, Fifth FloorBoston, MA 02110-1335, USA                  
//#                                                                            

#ifndef COMPONENTS_COMPONENTMODELS_TEST_FLUXSTANDARD2_GTEST_H
#define COMPONENTS_COMPONENTMODELS_TEST_FLUXSTANDARD2_GTEST_H

#include <memory>
#include <gtest/gtest.h>
#include <components/ComponentModels/test/FluxStandardTest.h>

#include <casa/aips.h>
#include <casa/aipstype.h>
#include <casa/BasicSL/String.h>
#include <components/ComponentModels/FluxStandard.h>
#include <measures/Measures/MEpoch.h>
#include <measures/Measures/MFrequency.h>

namespace test 
{

class NewFluxStandardTest: public FluxStandardTest 
{

public:

    NewFluxStandardTest();
    virtual ~NewFluxStandardTest();
    
protected:
    
    virtual void SetUp();
    virtual void TearDown();
    casacore::Bool modelExists(casacore::String modeltablename);

    casacore::String foundModelPath;
    casacore::String flxStdName;          
    casacore::String srcName;
    casacore::Double freq;
    casacore::MFrequency mfreq;
    casacore::MEpoch mtime;
    casacore::MDirection srcDir;
    casacore::String interpMethod;
    casa::FluxStandard::FluxScale flxStdEnum;
    casacore::Vector<casacore::Double> fluxUsed;
    casa::Flux<casacore::Double> returnFlux, returnFluxErr;
    casacore::Bool foundStd;
    casacore::String coeffsTbName;
    // expected values
    casa::FluxStandard::FluxScale expFlxStdEnum;
    casacore::Double expFlxVal;
};

class MatchStandardTest: public NewFluxStandardTest,
       public ::testing::WithParamInterface<std::tr1::tuple<casacore::String, casa::FluxStandard::FluxScale>>
{

public:

    MatchStandardTest();
    virtual ~MatchStandardTest();
    
protected:
    
    virtual void SetUp();
    virtual void TearDown();
    casacore::Bool matchedStandard;
    casacore::String flxStdDesc;

};//matchStandardTest


class AltSrcNameTest: public NewFluxStandardTest,
       public ::testing::WithParamInterface<std::tr1::tuple<casacore::String, casacore::String, casacore::Double, casa::FluxStandard::FluxScale>>
{

public:

    AltSrcNameTest();
    virtual ~AltSrcNameTest();
    
protected:
    
    virtual void SetUp();
    virtual void TearDown();
};

class AltSrcNamePB2017Test: public NewFluxStandardTest,
       public ::testing::WithParamInterface<std::tr1::tuple<casacore::String, casacore::String, casacore::Double, casa::FluxStandard::FluxScale>>
{

public:

    AltSrcNamePB2017Test();
    virtual ~AltSrcNamePB2017Test();
    
protected:
    
    virtual void SetUp();
    virtual void TearDown();
};

class SetInterpMethodTest: public NewFluxStandardTest,
       public ::testing::WithParamInterface<std::tr1::tuple<casacore::String, casacore::Double>>
{

public:

    SetInterpMethodTest();
    ~SetInterpMethodTest();

protected:
     
    virtual void SetUp();
    virtual void TearDown();
};

class FluxValueTest: public NewFluxStandardTest, 
       public ::testing::WithParamInterface<std::tr1::tuple<casacore::String, casacore::String, casacore::Double, casa::FluxStandard::FluxScale, casacore::Double>>
{

public:

    FluxValueTest();
    virtual ~FluxValueTest();
    
protected:
    
    virtual void SetUp();
    virtual void TearDown();
    //casacore::Double expFlxVal;
};

class PB2017FluxValueTest: public NewFluxStandardTest, 
       public ::testing::WithParamInterface<std::tr1::tuple<casacore::String, casacore::String, casacore::Double, casa::FluxStandard::FluxScale, casacore::Double>>
{

public:

    PB2017FluxValueTest();
    virtual ~PB2017FluxValueTest();
    
protected:
    
    virtual void SetUp();
    virtual void TearDown();
    //casacore::Double expFlxVal;
};

}// end namespace test

#endif /* COMPONENTS_COMPONENTMODELS_TEST_FLUXSTANDARD_GTEST_H */
