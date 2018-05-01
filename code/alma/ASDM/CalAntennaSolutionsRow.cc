
/*
 * ALMA - Atacama Large Millimeter Array
 * (c) European Southern Observatory, 2002
 * (c) Associated Universities Inc., 2002
 * Copyright by ESO (in the framework of the ALMA collaboration),
 * Copyright by AUI (in the framework of the ALMA collaboration),
 * All rights reserved.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY, without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA 02111-1307  USA
 *
 * Warning!
 *  -------------------------------------------------------------------- 
 * | This is generated code!  Do not modify this file.                  |
 * | If you do, all changes will be lost when the file is re-generated. |
 *  --------------------------------------------------------------------
 *
 * File CalAntennaSolutionsRow.cpp
 */
 
#include <vector>
using std::vector;

#include <set>
using std::set;

#include <ASDM.h>
#include <CalAntennaSolutionsRow.h>
#include <CalAntennaSolutionsTable.h>

#include <CalDataTable.h>
#include <CalDataRow.h>

#include <CalReductionTable.h>
#include <CalReductionRow.h>
	

using asdm::ASDM;
using asdm::CalAntennaSolutionsRow;
using asdm::CalAntennaSolutionsTable;

using asdm::CalDataTable;
using asdm::CalDataRow;

using asdm::CalReductionTable;
using asdm::CalReductionRow;


#include <Parser.h>
using asdm::Parser;

#include <EnumerationParser.h>
#include <ASDMValuesParser.h>
 
#include <InvalidArgumentException.h>
using asdm::InvalidArgumentException;

namespace asdm {
	CalAntennaSolutionsRow::~CalAntennaSolutionsRow() {
	}

	/**
	 * Return the table to which this row belongs.
	 */
	CalAntennaSolutionsTable &CalAntennaSolutionsRow::getTable() const {
		return table;
	}

	bool CalAntennaSolutionsRow::isAdded() const {
		return hasBeenAdded;
	}	

	void CalAntennaSolutionsRow::isAdded(bool added) {
		hasBeenAdded = added;
	}
	
#ifndef WITHOUT_ACS
	using asdmIDL::CalAntennaSolutionsRowIDL;
#endif
	
#ifndef WITHOUT_ACS
	/**
	 * Return this row in the form of an IDL struct.
	 * @return The values of this row as a CalAntennaSolutionsRowIDL struct.
	 */
	CalAntennaSolutionsRowIDL *CalAntennaSolutionsRow::toIDL() const {
		CalAntennaSolutionsRowIDL *x = new CalAntennaSolutionsRowIDL ();
		
		// Fill the IDL structure.
	
		
	
  		
		
		
			
				
		x->antennaName = CORBA::string_dup(antennaName.c_str());
				
 			
		
	

	
  		
		
		
			
				
		x->atmPhaseCorrection = atmPhaseCorrection;
 				
 			
		
	

	
  		
		
		
			
				
		x->basebandName = basebandName;
 				
 			
		
	

	
  		
		
		
			
				
		x->receiverBand = receiverBand;
 				
 			
		
	

	
  		
		
		
			
		x->startValidTime = startValidTime.toIDLArrayTime();
			
		
	

	
  		
		
		
			
		x->endValidTime = endValidTime.toIDLArrayTime();
			
		
	

	
  		
		
		
			
				
		x->numReceptor = numReceptor;
 				
 			
		
	

	
  		
		
		
			
				
		x->refAntennaName = CORBA::string_dup(refAntennaName.c_str());
				
 			
		
	

	
  		
		
		
			
		x->direction.length(direction.size());
		for (unsigned int i = 0; i < direction.size(); ++i) {
			
			x->direction[i] = direction.at(i).toIDLAngle();
			
	 	}
			
		
	

	
  		
		
		
			
		x->frequencyRange.length(frequencyRange.size());
		for (unsigned int i = 0; i < frequencyRange.size(); ++i) {
			
			x->frequencyRange[i] = frequencyRange.at(i).toIDLFrequency();
			
	 	}
			
		
	

	
  		
		
		
			
		x->integrationTime = integrationTime.toIDLInterval();
			
		
	

	
  		
		
		
			
		x->polarizationTypes.length(polarizationTypes.size());
		for (unsigned int i = 0; i < polarizationTypes.size(); ++i) {
			
				
			x->polarizationTypes[i] = polarizationTypes.at(i);
	 			
	 		
	 	}
			
		
	

	
  		
		
		
			
				
		x->correctionValidity = correctionValidity;
 				
 			
		
	

	
  		
		
		
			
		x->phaseAnt.length(phaseAnt.size());
		for (unsigned int i = 0; i < phaseAnt.size(); ++i) {
			
				
			x->phaseAnt[i] = phaseAnt.at(i);
	 			
	 		
	 	}
			
		
	

	
  		
		
		
			
		x->phaseAntRMS.length(phaseAntRMS.size());
		for (unsigned int i = 0; i < phaseAntRMS.size(); ++i) {
			
				
			x->phaseAntRMS[i] = phaseAntRMS.at(i);
	 			
	 		
	 	}
			
		
	

	
  		
		
		
			
		x->amplitudeAnt.length(amplitudeAnt.size());
		for (unsigned int i = 0; i < amplitudeAnt.size(); ++i) {
			
				
			x->amplitudeAnt[i] = amplitudeAnt.at(i);
	 			
	 		
	 	}
			
		
	

	
  		
		
		
			
		x->amplitudeAntRMS.length(amplitudeAntRMS.size());
		for (unsigned int i = 0; i < amplitudeAntRMS.size(); ++i) {
			
				
			x->amplitudeAntRMS[i] = amplitudeAntRMS.at(i);
	 			
	 		
	 	}
			
		
	

	
	
		
	
  	
 		
		
	 	
			
		x->calDataId = calDataId.toIDLTag();
			
	 	 		
  	

	
  	
 		
		
	 	
			
		x->calReductionId = calReductionId.toIDLTag();
			
	 	 		
  	

	
		
	

	

		
		return x;
	
	}
	
	void CalAntennaSolutionsRow::toIDL(asdmIDL::CalAntennaSolutionsRowIDL& x) const {
		// Set the x's fields.
	
		
	
  		
		
		
			
				
		x.antennaName = CORBA::string_dup(antennaName.c_str());
				
 			
		
	

	
  		
		
		
			
				
		x.atmPhaseCorrection = atmPhaseCorrection;
 				
 			
		
	

	
  		
		
		
			
				
		x.basebandName = basebandName;
 				
 			
		
	

	
  		
		
		
			
				
		x.receiverBand = receiverBand;
 				
 			
		
	

	
  		
		
		
			
		x.startValidTime = startValidTime.toIDLArrayTime();
			
		
	

	
  		
		
		
			
		x.endValidTime = endValidTime.toIDLArrayTime();
			
		
	

	
  		
		
		
			
				
		x.numReceptor = numReceptor;
 				
 			
		
	

	
  		
		
		
			
				
		x.refAntennaName = CORBA::string_dup(refAntennaName.c_str());
				
 			
		
	

	
  		
		
		
			
		x.direction.length(direction.size());
		for (unsigned int i = 0; i < direction.size(); ++i) {
			
			x.direction[i] = direction.at(i).toIDLAngle();
			
	 	}
			
		
	

	
  		
		
		
			
		x.frequencyRange.length(frequencyRange.size());
		for (unsigned int i = 0; i < frequencyRange.size(); ++i) {
			
			x.frequencyRange[i] = frequencyRange.at(i).toIDLFrequency();
			
	 	}
			
		
	

	
  		
		
		
			
		x.integrationTime = integrationTime.toIDLInterval();
			
		
	

	
  		
		
		
			
		x.polarizationTypes.length(polarizationTypes.size());
		for (unsigned int i = 0; i < polarizationTypes.size(); ++i) {
			
				
			x.polarizationTypes[i] = polarizationTypes.at(i);
	 			
	 		
	 	}
			
		
	

	
  		
		
		
			
				
		x.correctionValidity = correctionValidity;
 				
 			
		
	

	
  		
		
		
			
		x.phaseAnt.length(phaseAnt.size());
		for (unsigned int i = 0; i < phaseAnt.size(); ++i) {
			
				
			x.phaseAnt[i] = phaseAnt.at(i);
	 			
	 		
	 	}
			
		
	

	
  		
		
		
			
		x.phaseAntRMS.length(phaseAntRMS.size());
		for (unsigned int i = 0; i < phaseAntRMS.size(); ++i) {
			
				
			x.phaseAntRMS[i] = phaseAntRMS.at(i);
	 			
	 		
	 	}
			
		
	

	
  		
		
		
			
		x.amplitudeAnt.length(amplitudeAnt.size());
		for (unsigned int i = 0; i < amplitudeAnt.size(); ++i) {
			
				
			x.amplitudeAnt[i] = amplitudeAnt.at(i);
	 			
	 		
	 	}
			
		
	

	
  		
		
		
			
		x.amplitudeAntRMS.length(amplitudeAntRMS.size());
		for (unsigned int i = 0; i < amplitudeAntRMS.size(); ++i) {
			
				
			x.amplitudeAntRMS[i] = amplitudeAntRMS.at(i);
	 			
	 		
	 	}
			
		
	

	
	
		
	
  	
 		
		
	 	
			
		x.calDataId = calDataId.toIDLTag();
			
	 	 		
  	

	
  	
 		
		
	 	
			
		x.calReductionId = calReductionId.toIDLTag();
			
	 	 		
  	

	
		
	

	

	
	}
#endif
	

#ifndef WITHOUT_ACS
	/**
	 * Fill the values of this row from the IDL struct CalAntennaSolutionsRowIDL.
	 * @param x The IDL struct containing the values used to fill this row.
	 */
	void CalAntennaSolutionsRow::setFromIDL (CalAntennaSolutionsRowIDL x){
		try {
		// Fill the values from x.
	
		
	
		
		
			
		setAntennaName(string (x.antennaName));
			
 		
		
	

	
		
		
			
		setAtmPhaseCorrection(x.atmPhaseCorrection);
  			
 		
		
	

	
		
		
			
		setBasebandName(x.basebandName);
  			
 		
		
	

	
		
		
			
		setReceiverBand(x.receiverBand);
  			
 		
		
	

	
		
		
			
		setStartValidTime(ArrayTime (x.startValidTime));
			
 		
		
	

	
		
		
			
		setEndValidTime(ArrayTime (x.endValidTime));
			
 		
		
	

	
		
		
			
		setNumReceptor(x.numReceptor);
  			
 		
		
	

	
		
		
			
		setRefAntennaName(string (x.refAntennaName));
			
 		
		
	

	
		
		
			
		direction .clear();
		for (unsigned int i = 0; i <x.direction.length(); ++i) {
			
			direction.push_back(Angle (x.direction[i]));
			
		}
			
  		
		
	

	
		
		
			
		frequencyRange .clear();
		for (unsigned int i = 0; i <x.frequencyRange.length(); ++i) {
			
			frequencyRange.push_back(Frequency (x.frequencyRange[i]));
			
		}
			
  		
		
	

	
		
		
			
		setIntegrationTime(Interval (x.integrationTime));
			
 		
		
	

	
		
		
			
		polarizationTypes .clear();
		for (unsigned int i = 0; i <x.polarizationTypes.length(); ++i) {
			
			polarizationTypes.push_back(x.polarizationTypes[i]);
  			
		}
			
  		
		
	

	
		
		
			
		setCorrectionValidity(x.correctionValidity);
  			
 		
		
	

	
		
		
			
		phaseAnt .clear();
		for (unsigned int i = 0; i <x.phaseAnt.length(); ++i) {
			
			phaseAnt.push_back(x.phaseAnt[i]);
  			
		}
			
  		
		
	

	
		
		
			
		phaseAntRMS .clear();
		for (unsigned int i = 0; i <x.phaseAntRMS.length(); ++i) {
			
			phaseAntRMS.push_back(x.phaseAntRMS[i]);
  			
		}
			
  		
		
	

	
		
		
			
		amplitudeAnt .clear();
		for (unsigned int i = 0; i <x.amplitudeAnt.length(); ++i) {
			
			amplitudeAnt.push_back(x.amplitudeAnt[i]);
  			
		}
			
  		
		
	

	
		
		
			
		amplitudeAntRMS .clear();
		for (unsigned int i = 0; i <x.amplitudeAntRMS.length(); ++i) {
			
			amplitudeAntRMS.push_back(x.amplitudeAntRMS[i]);
  			
		}
			
  		
		
	

	
	
		
	
		
		
			
		setCalDataId(Tag (x.calDataId));
			
 		
		
	

	
		
		
			
		setCalReductionId(Tag (x.calReductionId));
			
 		
		
	

	
		
	

	

		} catch (IllegalAccessException err) {
			throw ConversionException (err.getMessage(),"CalAntennaSolutions");
		}
	}
#endif
	
	/**
	 * Return this row in the form of an XML string.
	 * @return The values of this row as an XML string.
	 */
	string CalAntennaSolutionsRow::toXML() const {
		string buf;
		buf.append("<row> \n");
		
	
		
  	
 		
		
		Parser::toXML(antennaName, "antennaName", buf);
		
		
	

  	
 		
		
			buf.append(EnumerationParser::toXML("atmPhaseCorrection", atmPhaseCorrection));
		
		
	

  	
 		
		
			buf.append(EnumerationParser::toXML("basebandName", basebandName));
		
		
	

  	
 		
		
			buf.append(EnumerationParser::toXML("receiverBand", receiverBand));
		
		
	

  	
 		
		
		Parser::toXML(startValidTime, "startValidTime", buf);
		
		
	

  	
 		
		
		Parser::toXML(endValidTime, "endValidTime", buf);
		
		
	

  	
 		
		
		Parser::toXML(numReceptor, "numReceptor", buf);
		
		
	

  	
 		
		
		Parser::toXML(refAntennaName, "refAntennaName", buf);
		
		
	

  	
 		
		
		Parser::toXML(direction, "direction", buf);
		
		
	

  	
 		
		
		Parser::toXML(frequencyRange, "frequencyRange", buf);
		
		
	

  	
 		
		
		Parser::toXML(integrationTime, "integrationTime", buf);
		
		
	

  	
 		
		
			buf.append(EnumerationParser::toXML("polarizationTypes", polarizationTypes));
		
		
	

  	
 		
		
		Parser::toXML(correctionValidity, "correctionValidity", buf);
		
		
	

  	
 		
		
		Parser::toXML(phaseAnt, "phaseAnt", buf);
		
		
	

  	
 		
		
		Parser::toXML(phaseAntRMS, "phaseAntRMS", buf);
		
		
	

  	
 		
		
		Parser::toXML(amplitudeAnt, "amplitudeAnt", buf);
		
		
	

  	
 		
		
		Parser::toXML(amplitudeAntRMS, "amplitudeAntRMS", buf);
		
		
	

	
	
		
  	
 		
		
		Parser::toXML(calDataId, "calDataId", buf);
		
		
	

  	
 		
		
		Parser::toXML(calReductionId, "calReductionId", buf);
		
		
	

	
		
	

	

		
		buf.append("</row>\n");
		return buf;
	}

	/**
	 * Fill the values of this row from an XML string 
	 * that was produced by the toXML() method.
	 * @param x The XML string being used to set the values of this row.
	 */
	void CalAntennaSolutionsRow::setFromXML (string rowDoc) {
		Parser row(rowDoc);
		string s = "";
		try {
	
		
	
  		
			
	  	setAntennaName(Parser::getString("antennaName","CalAntennaSolutions",rowDoc));
			
		
	

	
		
		
		
		atmPhaseCorrection = EnumerationParser::getAtmPhaseCorrection("atmPhaseCorrection","CalAntennaSolutions",rowDoc);
		
		
		
	

	
		
		
		
		basebandName = EnumerationParser::getBasebandName("basebandName","CalAntennaSolutions",rowDoc);
		
		
		
	

	
		
		
		
		receiverBand = EnumerationParser::getReceiverBand("receiverBand","CalAntennaSolutions",rowDoc);
		
		
		
	

	
  		
			
	  	setStartValidTime(Parser::getArrayTime("startValidTime","CalAntennaSolutions",rowDoc));
			
		
	

	
  		
			
	  	setEndValidTime(Parser::getArrayTime("endValidTime","CalAntennaSolutions",rowDoc));
			
		
	

	
  		
			
	  	setNumReceptor(Parser::getInteger("numReceptor","CalAntennaSolutions",rowDoc));
			
		
	

	
  		
			
	  	setRefAntennaName(Parser::getString("refAntennaName","CalAntennaSolutions",rowDoc));
			
		
	

	
  		
			
					
	  	setDirection(Parser::get1DAngle("direction","CalAntennaSolutions",rowDoc));
	  			
	  		
		
	

	
  		
			
					
	  	setFrequencyRange(Parser::get1DFrequency("frequencyRange","CalAntennaSolutions",rowDoc));
	  			
	  		
		
	

	
  		
			
	  	setIntegrationTime(Parser::getInterval("integrationTime","CalAntennaSolutions",rowDoc));
			
		
	

	
		
		
		
		polarizationTypes = EnumerationParser::getPolarizationType1D("polarizationTypes","CalAntennaSolutions",rowDoc);			
		
		
		
	

	
  		
			
	  	setCorrectionValidity(Parser::getBoolean("correctionValidity","CalAntennaSolutions",rowDoc));
			
		
	

	
  		
			
					
	  	setPhaseAnt(Parser::get1DFloat("phaseAnt","CalAntennaSolutions",rowDoc));
	  			
	  		
		
	

	
  		
			
					
	  	setPhaseAntRMS(Parser::get1DFloat("phaseAntRMS","CalAntennaSolutions",rowDoc));
	  			
	  		
		
	

	
  		
			
					
	  	setAmplitudeAnt(Parser::get1DFloat("amplitudeAnt","CalAntennaSolutions",rowDoc));
	  			
	  		
		
	

	
  		
			
					
	  	setAmplitudeAntRMS(Parser::get1DFloat("amplitudeAntRMS","CalAntennaSolutions",rowDoc));
	  			
	  		
		
	

	
	
		
	
  		
			
	  	setCalDataId(Parser::getTag("calDataId","CalData",rowDoc));
			
		
	

	
  		
			
	  	setCalReductionId(Parser::getTag("calReductionId","CalReduction",rowDoc));
			
		
	

	
		
	

	

		} catch (IllegalAccessException err) {
			throw ConversionException (err.getMessage(),"CalAntennaSolutions");
		}
	}
	
	void CalAntennaSolutionsRow::toBin(EndianOSStream& eoss) {
	
	
	
	
		
						
			eoss.writeString(antennaName);
				
		
	

	
	
		
					
			eoss.writeString(CAtmPhaseCorrection::name(atmPhaseCorrection));
			/* eoss.writeInt(atmPhaseCorrection); */
				
		
	

	
	
		
					
			eoss.writeString(CReceiverBand::name(receiverBand));
			/* eoss.writeInt(receiverBand); */
				
		
	

	
	
		
					
			eoss.writeString(CBasebandName::name(basebandName));
			/* eoss.writeInt(basebandName); */
				
		
	

	
	
		
	calDataId.toBin(eoss);
		
	

	
	
		
	calReductionId.toBin(eoss);
		
	

	
	
		
	startValidTime.toBin(eoss);
		
	

	
	
		
	endValidTime.toBin(eoss);
		
	

	
	
		
						
			eoss.writeInt(numReceptor);
				
		
	

	
	
		
						
			eoss.writeString(refAntennaName);
				
		
	

	
	
		
	Angle::toBin(direction, eoss);
		
	

	
	
		
	Frequency::toBin(frequencyRange, eoss);
		
	

	
	
		
	integrationTime.toBin(eoss);
		
	

	
	
		
		
			
		eoss.writeInt((int) polarizationTypes.size());
		for (unsigned int i = 0; i < polarizationTypes.size(); i++)
				
			eoss.writeString(CPolarizationType::name(polarizationTypes.at(i)));
			/* eoss.writeInt(polarizationTypes.at(i)); */
				
				
						
		
	

	
	
		
						
			eoss.writeBoolean(correctionValidity);
				
		
	

	
	
		
		
			
		eoss.writeInt((int) phaseAnt.size());
		for (unsigned int i = 0; i < phaseAnt.size(); i++)
				
			eoss.writeFloat(phaseAnt.at(i));
				
				
						
		
	

	
	
		
		
			
		eoss.writeInt((int) phaseAntRMS.size());
		for (unsigned int i = 0; i < phaseAntRMS.size(); i++)
				
			eoss.writeFloat(phaseAntRMS.at(i));
				
				
						
		
	

	
	
		
		
			
		eoss.writeInt((int) amplitudeAnt.size());
		for (unsigned int i = 0; i < amplitudeAnt.size(); i++)
				
			eoss.writeFloat(amplitudeAnt.at(i));
				
				
						
		
	

	
	
		
		
			
		eoss.writeInt((int) amplitudeAntRMS.size());
		for (unsigned int i = 0; i < amplitudeAntRMS.size(); i++)
				
			eoss.writeFloat(amplitudeAntRMS.at(i));
				
				
						
		
	


	
	
	}
	
void CalAntennaSolutionsRow::antennaNameFromBin(EndianIStream& eis) {
		
	
	
		
			
		antennaName =  eis.readString();
			
		
	
	
}
void CalAntennaSolutionsRow::atmPhaseCorrectionFromBin(EndianIStream& eis) {
		
	
	
		
			
		atmPhaseCorrection = CAtmPhaseCorrection::literal(eis.readString());
			
		
	
	
}
void CalAntennaSolutionsRow::receiverBandFromBin(EndianIStream& eis) {
		
	
	
		
			
		receiverBand = CReceiverBand::literal(eis.readString());
			
		
	
	
}
void CalAntennaSolutionsRow::basebandNameFromBin(EndianIStream& eis) {
		
	
	
		
			
		basebandName = CBasebandName::literal(eis.readString());
			
		
	
	
}
void CalAntennaSolutionsRow::calDataIdFromBin(EndianIStream& eis) {
		
	
		
		
		calDataId =  Tag::fromBin(eis);
		
	
	
}
void CalAntennaSolutionsRow::calReductionIdFromBin(EndianIStream& eis) {
		
	
		
		
		calReductionId =  Tag::fromBin(eis);
		
	
	
}
void CalAntennaSolutionsRow::startValidTimeFromBin(EndianIStream& eis) {
		
	
		
		
		startValidTime =  ArrayTime::fromBin(eis);
		
	
	
}
void CalAntennaSolutionsRow::endValidTimeFromBin(EndianIStream& eis) {
		
	
		
		
		endValidTime =  ArrayTime::fromBin(eis);
		
	
	
}
void CalAntennaSolutionsRow::numReceptorFromBin(EndianIStream& eis) {
		
	
	
		
			
		numReceptor =  eis.readInt();
			
		
	
	
}
void CalAntennaSolutionsRow::refAntennaNameFromBin(EndianIStream& eis) {
		
	
	
		
			
		refAntennaName =  eis.readString();
			
		
	
	
}
void CalAntennaSolutionsRow::directionFromBin(EndianIStream& eis) {
		
	
		
		
			
	
	direction = Angle::from1DBin(eis);	
	

		
	
	
}
void CalAntennaSolutionsRow::frequencyRangeFromBin(EndianIStream& eis) {
		
	
		
		
			
	
	frequencyRange = Frequency::from1DBin(eis);	
	

		
	
	
}
void CalAntennaSolutionsRow::integrationTimeFromBin(EndianIStream& eis) {
		
	
		
		
		integrationTime =  Interval::fromBin(eis);
		
	
	
}
void CalAntennaSolutionsRow::polarizationTypesFromBin(EndianIStream& eis) {
		
	
	
		
			
	
		polarizationTypes.clear();
		
		unsigned int polarizationTypesDim1 = eis.readInt();
		for (unsigned int  i = 0 ; i < polarizationTypesDim1; i++)
			
			polarizationTypes.push_back(CPolarizationType::literal(eis.readString()));
			
	

		
	
	
}
void CalAntennaSolutionsRow::correctionValidityFromBin(EndianIStream& eis) {
		
	
	
		
			
		correctionValidity =  eis.readBoolean();
			
		
	
	
}
void CalAntennaSolutionsRow::phaseAntFromBin(EndianIStream& eis) {
		
	
	
		
			
	
		phaseAnt.clear();
		
		unsigned int phaseAntDim1 = eis.readInt();
		for (unsigned int  i = 0 ; i < phaseAntDim1; i++)
			
			phaseAnt.push_back(eis.readFloat());
			
	

		
	
	
}
void CalAntennaSolutionsRow::phaseAntRMSFromBin(EndianIStream& eis) {
		
	
	
		
			
	
		phaseAntRMS.clear();
		
		unsigned int phaseAntRMSDim1 = eis.readInt();
		for (unsigned int  i = 0 ; i < phaseAntRMSDim1; i++)
			
			phaseAntRMS.push_back(eis.readFloat());
			
	

		
	
	
}
void CalAntennaSolutionsRow::amplitudeAntFromBin(EndianIStream& eis) {
		
	
	
		
			
	
		amplitudeAnt.clear();
		
		unsigned int amplitudeAntDim1 = eis.readInt();
		for (unsigned int  i = 0 ; i < amplitudeAntDim1; i++)
			
			amplitudeAnt.push_back(eis.readFloat());
			
	

		
	
	
}
void CalAntennaSolutionsRow::amplitudeAntRMSFromBin(EndianIStream& eis) {
		
	
	
		
			
	
		amplitudeAntRMS.clear();
		
		unsigned int amplitudeAntRMSDim1 = eis.readInt();
		for (unsigned int  i = 0 ; i < amplitudeAntRMSDim1; i++)
			
			amplitudeAntRMS.push_back(eis.readFloat());
			
	

		
	
	
}

		
	
	CalAntennaSolutionsRow* CalAntennaSolutionsRow::fromBin(EndianIStream& eis, CalAntennaSolutionsTable& table, const vector<string>& attributesSeq) {
		CalAntennaSolutionsRow* row = new  CalAntennaSolutionsRow(table);
		
		map<string, CalAntennaSolutionsAttributeFromBin>::iterator iter ;
		for (unsigned int i = 0; i < attributesSeq.size(); i++) {
			iter = row->fromBinMethods.find(attributesSeq.at(i));
			if (iter != row->fromBinMethods.end()) {
				(row->*(row->fromBinMethods[ attributesSeq.at(i) ] ))(eis);			
			}
			else {
				BinaryAttributeReaderFunctor* functorP = table.getUnknownAttributeBinaryReader(attributesSeq.at(i));
				if (functorP)
					(*functorP)(eis);
				else
					throw ConversionException("There is not method to read an attribute '"+attributesSeq.at(i)+"'.", "CalAntennaSolutionsTable");
			}
				
		}				
		return row;
	}

	//
	// A collection of methods to set the value of the attributes from their textual value in the XML representation
	// of one row.
	//
	
	// Convert a string into an String 
	void CalAntennaSolutionsRow::antennaNameFromText(const string & s) {
		 
		antennaName = ASDMValuesParser::parse<string>(s);
		
	}
	
	
	// Convert a string into an AtmPhaseCorrection 
	void CalAntennaSolutionsRow::atmPhaseCorrectionFromText(const string & s) {
		 
		atmPhaseCorrection = ASDMValuesParser::parse<AtmPhaseCorrection>(s);
		
	}
	
	
	// Convert a string into an ReceiverBand 
	void CalAntennaSolutionsRow::receiverBandFromText(const string & s) {
		 
		receiverBand = ASDMValuesParser::parse<ReceiverBand>(s);
		
	}
	
	
	// Convert a string into an BasebandName 
	void CalAntennaSolutionsRow::basebandNameFromText(const string & s) {
		 
		basebandName = ASDMValuesParser::parse<BasebandName>(s);
		
	}
	
	
	// Convert a string into an Tag 
	void CalAntennaSolutionsRow::calDataIdFromText(const string & s) {
		 
		calDataId = ASDMValuesParser::parse<Tag>(s);
		
	}
	
	
	// Convert a string into an Tag 
	void CalAntennaSolutionsRow::calReductionIdFromText(const string & s) {
		 
		calReductionId = ASDMValuesParser::parse<Tag>(s);
		
	}
	
	
	// Convert a string into an ArrayTime 
	void CalAntennaSolutionsRow::startValidTimeFromText(const string & s) {
		 
		startValidTime = ASDMValuesParser::parse<ArrayTime>(s);
		
	}
	
	
	// Convert a string into an ArrayTime 
	void CalAntennaSolutionsRow::endValidTimeFromText(const string & s) {
		 
		endValidTime = ASDMValuesParser::parse<ArrayTime>(s);
		
	}
	
	
	// Convert a string into an int 
	void CalAntennaSolutionsRow::numReceptorFromText(const string & s) {
		 
		numReceptor = ASDMValuesParser::parse<int>(s);
		
	}
	
	
	// Convert a string into an String 
	void CalAntennaSolutionsRow::refAntennaNameFromText(const string & s) {
		 
		refAntennaName = ASDMValuesParser::parse<string>(s);
		
	}
	
	
	// Convert a string into an Angle 
	void CalAntennaSolutionsRow::directionFromText(const string & s) {
		 
		direction = ASDMValuesParser::parse1D<Angle>(s);
		
	}
	
	
	// Convert a string into an Frequency 
	void CalAntennaSolutionsRow::frequencyRangeFromText(const string & s) {
		 
		frequencyRange = ASDMValuesParser::parse1D<Frequency>(s);
		
	}
	
	
	// Convert a string into an Interval 
	void CalAntennaSolutionsRow::integrationTimeFromText(const string & s) {
		 
		integrationTime = ASDMValuesParser::parse<Interval>(s);
		
	}
	
	
	// Convert a string into an PolarizationType 
	void CalAntennaSolutionsRow::polarizationTypesFromText(const string & s) {
		 
		polarizationTypes = ASDMValuesParser::parse1D<PolarizationType>(s);
		
	}
	
	
	// Convert a string into an boolean 
	void CalAntennaSolutionsRow::correctionValidityFromText(const string & s) {
		 
		correctionValidity = ASDMValuesParser::parse<bool>(s);
		
	}
	
	
	// Convert a string into an float 
	void CalAntennaSolutionsRow::phaseAntFromText(const string & s) {
		 
		phaseAnt = ASDMValuesParser::parse1D<float>(s);
		
	}
	
	
	// Convert a string into an float 
	void CalAntennaSolutionsRow::phaseAntRMSFromText(const string & s) {
		 
		phaseAntRMS = ASDMValuesParser::parse1D<float>(s);
		
	}
	
	
	// Convert a string into an float 
	void CalAntennaSolutionsRow::amplitudeAntFromText(const string & s) {
		 
		amplitudeAnt = ASDMValuesParser::parse1D<float>(s);
		
	}
	
	
	// Convert a string into an float 
	void CalAntennaSolutionsRow::amplitudeAntRMSFromText(const string & s) {
		 
		amplitudeAntRMS = ASDMValuesParser::parse1D<float>(s);
		
	}
	

		
	
	void CalAntennaSolutionsRow::fromText(const std::string& attributeName, const std::string&  t) {
		map<string, CalAntennaSolutionsAttributeFromText>::iterator iter;
		if ((iter = fromTextMethods.find(attributeName)) == fromTextMethods.end())
			throw ConversionException("I do not know what to do with '"+attributeName+"' and its content '"+t+"' (while parsing an XML document)", "CalAntennaSolutionsTable");
		(this->*(iter->second))(t);
	}
			
	////////////////////////////////////////////////
	// Intrinsic Table Attributes getters/setters //
	////////////////////////////////////////////////
	
	

	
 	/**
 	 * Get antennaName.
 	 * @return antennaName as string
 	 */
 	string CalAntennaSolutionsRow::getAntennaName() const {
	
  		return antennaName;
 	}

 	/**
 	 * Set antennaName with the specified string.
 	 * @param antennaName The string value to which antennaName is to be set.
 	 
 	
 		
 	 * @throw IllegalAccessException If an attempt is made to change this field after is has been added to the table.
 	 	
 	 */
 	void CalAntennaSolutionsRow::setAntennaName (string antennaName)  {
  	
  	
  		if (hasBeenAdded) {
 		
			throw IllegalAccessException("antennaName", "CalAntennaSolutions");
		
  		}
  	
 		this->antennaName = antennaName;
	
 	}
	
	

	

	
 	/**
 	 * Get atmPhaseCorrection.
 	 * @return atmPhaseCorrection as AtmPhaseCorrectionMod::AtmPhaseCorrection
 	 */
 	AtmPhaseCorrectionMod::AtmPhaseCorrection CalAntennaSolutionsRow::getAtmPhaseCorrection() const {
	
  		return atmPhaseCorrection;
 	}

 	/**
 	 * Set atmPhaseCorrection with the specified AtmPhaseCorrectionMod::AtmPhaseCorrection.
 	 * @param atmPhaseCorrection The AtmPhaseCorrectionMod::AtmPhaseCorrection value to which atmPhaseCorrection is to be set.
 	 
 	
 		
 	 * @throw IllegalAccessException If an attempt is made to change this field after is has been added to the table.
 	 	
 	 */
 	void CalAntennaSolutionsRow::setAtmPhaseCorrection (AtmPhaseCorrectionMod::AtmPhaseCorrection atmPhaseCorrection)  {
  	
  	
  		if (hasBeenAdded) {
 		
			throw IllegalAccessException("atmPhaseCorrection", "CalAntennaSolutions");
		
  		}
  	
 		this->atmPhaseCorrection = atmPhaseCorrection;
	
 	}
	
	

	

	
 	/**
 	 * Get basebandName.
 	 * @return basebandName as BasebandNameMod::BasebandName
 	 */
 	BasebandNameMod::BasebandName CalAntennaSolutionsRow::getBasebandName() const {
	
  		return basebandName;
 	}

 	/**
 	 * Set basebandName with the specified BasebandNameMod::BasebandName.
 	 * @param basebandName The BasebandNameMod::BasebandName value to which basebandName is to be set.
 	 
 	
 		
 	 * @throw IllegalAccessException If an attempt is made to change this field after is has been added to the table.
 	 	
 	 */
 	void CalAntennaSolutionsRow::setBasebandName (BasebandNameMod::BasebandName basebandName)  {
  	
  	
  		if (hasBeenAdded) {
 		
			throw IllegalAccessException("basebandName", "CalAntennaSolutions");
		
  		}
  	
 		this->basebandName = basebandName;
	
 	}
	
	

	

	
 	/**
 	 * Get receiverBand.
 	 * @return receiverBand as ReceiverBandMod::ReceiverBand
 	 */
 	ReceiverBandMod::ReceiverBand CalAntennaSolutionsRow::getReceiverBand() const {
	
  		return receiverBand;
 	}

 	/**
 	 * Set receiverBand with the specified ReceiverBandMod::ReceiverBand.
 	 * @param receiverBand The ReceiverBandMod::ReceiverBand value to which receiverBand is to be set.
 	 
 	
 		
 	 * @throw IllegalAccessException If an attempt is made to change this field after is has been added to the table.
 	 	
 	 */
 	void CalAntennaSolutionsRow::setReceiverBand (ReceiverBandMod::ReceiverBand receiverBand)  {
  	
  	
  		if (hasBeenAdded) {
 		
			throw IllegalAccessException("receiverBand", "CalAntennaSolutions");
		
  		}
  	
 		this->receiverBand = receiverBand;
	
 	}
	
	

	

	
 	/**
 	 * Get startValidTime.
 	 * @return startValidTime as ArrayTime
 	 */
 	ArrayTime CalAntennaSolutionsRow::getStartValidTime() const {
	
  		return startValidTime;
 	}

 	/**
 	 * Set startValidTime with the specified ArrayTime.
 	 * @param startValidTime The ArrayTime value to which startValidTime is to be set.
 	 
 	
 		
 	 */
 	void CalAntennaSolutionsRow::setStartValidTime (ArrayTime startValidTime)  {
  	
  	
  		if (hasBeenAdded) {
 		
  		}
  	
 		this->startValidTime = startValidTime;
	
 	}
	
	

	

	
 	/**
 	 * Get endValidTime.
 	 * @return endValidTime as ArrayTime
 	 */
 	ArrayTime CalAntennaSolutionsRow::getEndValidTime() const {
	
  		return endValidTime;
 	}

 	/**
 	 * Set endValidTime with the specified ArrayTime.
 	 * @param endValidTime The ArrayTime value to which endValidTime is to be set.
 	 
 	
 		
 	 */
 	void CalAntennaSolutionsRow::setEndValidTime (ArrayTime endValidTime)  {
  	
  	
  		if (hasBeenAdded) {
 		
  		}
  	
 		this->endValidTime = endValidTime;
	
 	}
	
	

	

	
 	/**
 	 * Get numReceptor.
 	 * @return numReceptor as int
 	 */
 	int CalAntennaSolutionsRow::getNumReceptor() const {
	
  		return numReceptor;
 	}

 	/**
 	 * Set numReceptor with the specified int.
 	 * @param numReceptor The int value to which numReceptor is to be set.
 	 
 	
 		
 	 */
 	void CalAntennaSolutionsRow::setNumReceptor (int numReceptor)  {
  	
  	
  		if (hasBeenAdded) {
 		
  		}
  	
 		this->numReceptor = numReceptor;
	
 	}
	
	

	

	
 	/**
 	 * Get refAntennaName.
 	 * @return refAntennaName as string
 	 */
 	string CalAntennaSolutionsRow::getRefAntennaName() const {
	
  		return refAntennaName;
 	}

 	/**
 	 * Set refAntennaName with the specified string.
 	 * @param refAntennaName The string value to which refAntennaName is to be set.
 	 
 	
 		
 	 */
 	void CalAntennaSolutionsRow::setRefAntennaName (string refAntennaName)  {
  	
  	
  		if (hasBeenAdded) {
 		
  		}
  	
 		this->refAntennaName = refAntennaName;
	
 	}
	
	

	

	
 	/**
 	 * Get direction.
 	 * @return direction as vector<Angle >
 	 */
 	vector<Angle > CalAntennaSolutionsRow::getDirection() const {
	
  		return direction;
 	}

 	/**
 	 * Set direction with the specified vector<Angle >.
 	 * @param direction The vector<Angle > value to which direction is to be set.
 	 
 	
 		
 	 */
 	void CalAntennaSolutionsRow::setDirection (vector<Angle > direction)  {
  	
  	
  		if (hasBeenAdded) {
 		
  		}
  	
 		this->direction = direction;
	
 	}
	
	

	

	
 	/**
 	 * Get frequencyRange.
 	 * @return frequencyRange as vector<Frequency >
 	 */
 	vector<Frequency > CalAntennaSolutionsRow::getFrequencyRange() const {
	
  		return frequencyRange;
 	}

 	/**
 	 * Set frequencyRange with the specified vector<Frequency >.
 	 * @param frequencyRange The vector<Frequency > value to which frequencyRange is to be set.
 	 
 	
 		
 	 */
 	void CalAntennaSolutionsRow::setFrequencyRange (vector<Frequency > frequencyRange)  {
  	
  	
  		if (hasBeenAdded) {
 		
  		}
  	
 		this->frequencyRange = frequencyRange;
	
 	}
	
	

	

	
 	/**
 	 * Get integrationTime.
 	 * @return integrationTime as Interval
 	 */
 	Interval CalAntennaSolutionsRow::getIntegrationTime() const {
	
  		return integrationTime;
 	}

 	/**
 	 * Set integrationTime with the specified Interval.
 	 * @param integrationTime The Interval value to which integrationTime is to be set.
 	 
 	
 		
 	 */
 	void CalAntennaSolutionsRow::setIntegrationTime (Interval integrationTime)  {
  	
  	
  		if (hasBeenAdded) {
 		
  		}
  	
 		this->integrationTime = integrationTime;
	
 	}
	
	

	

	
 	/**
 	 * Get polarizationTypes.
 	 * @return polarizationTypes as vector<PolarizationTypeMod::PolarizationType >
 	 */
 	vector<PolarizationTypeMod::PolarizationType > CalAntennaSolutionsRow::getPolarizationTypes() const {
	
  		return polarizationTypes;
 	}

 	/**
 	 * Set polarizationTypes with the specified vector<PolarizationTypeMod::PolarizationType >.
 	 * @param polarizationTypes The vector<PolarizationTypeMod::PolarizationType > value to which polarizationTypes is to be set.
 	 
 	
 		
 	 */
 	void CalAntennaSolutionsRow::setPolarizationTypes (vector<PolarizationTypeMod::PolarizationType > polarizationTypes)  {
  	
  	
  		if (hasBeenAdded) {
 		
  		}
  	
 		this->polarizationTypes = polarizationTypes;
	
 	}
	
	

	

	
 	/**
 	 * Get correctionValidity.
 	 * @return correctionValidity as bool
 	 */
 	bool CalAntennaSolutionsRow::getCorrectionValidity() const {
	
  		return correctionValidity;
 	}

 	/**
 	 * Set correctionValidity with the specified bool.
 	 * @param correctionValidity The bool value to which correctionValidity is to be set.
 	 
 	
 		
 	 */
 	void CalAntennaSolutionsRow::setCorrectionValidity (bool correctionValidity)  {
  	
  	
  		if (hasBeenAdded) {
 		
  		}
  	
 		this->correctionValidity = correctionValidity;
	
 	}
	
	

	

	
 	/**
 	 * Get phaseAnt.
 	 * @return phaseAnt as vector<float >
 	 */
 	vector<float > CalAntennaSolutionsRow::getPhaseAnt() const {
	
  		return phaseAnt;
 	}

 	/**
 	 * Set phaseAnt with the specified vector<float >.
 	 * @param phaseAnt The vector<float > value to which phaseAnt is to be set.
 	 
 	
 		
 	 */
 	void CalAntennaSolutionsRow::setPhaseAnt (vector<float > phaseAnt)  {
  	
  	
  		if (hasBeenAdded) {
 		
  		}
  	
 		this->phaseAnt = phaseAnt;
	
 	}
	
	

	

	
 	/**
 	 * Get phaseAntRMS.
 	 * @return phaseAntRMS as vector<float >
 	 */
 	vector<float > CalAntennaSolutionsRow::getPhaseAntRMS() const {
	
  		return phaseAntRMS;
 	}

 	/**
 	 * Set phaseAntRMS with the specified vector<float >.
 	 * @param phaseAntRMS The vector<float > value to which phaseAntRMS is to be set.
 	 
 	
 		
 	 */
 	void CalAntennaSolutionsRow::setPhaseAntRMS (vector<float > phaseAntRMS)  {
  	
  	
  		if (hasBeenAdded) {
 		
  		}
  	
 		this->phaseAntRMS = phaseAntRMS;
	
 	}
	
	

	

	
 	/**
 	 * Get amplitudeAnt.
 	 * @return amplitudeAnt as vector<float >
 	 */
 	vector<float > CalAntennaSolutionsRow::getAmplitudeAnt() const {
	
  		return amplitudeAnt;
 	}

 	/**
 	 * Set amplitudeAnt with the specified vector<float >.
 	 * @param amplitudeAnt The vector<float > value to which amplitudeAnt is to be set.
 	 
 	
 		
 	 */
 	void CalAntennaSolutionsRow::setAmplitudeAnt (vector<float > amplitudeAnt)  {
  	
  	
  		if (hasBeenAdded) {
 		
  		}
  	
 		this->amplitudeAnt = amplitudeAnt;
	
 	}
	
	

	

	
 	/**
 	 * Get amplitudeAntRMS.
 	 * @return amplitudeAntRMS as vector<float >
 	 */
 	vector<float > CalAntennaSolutionsRow::getAmplitudeAntRMS() const {
	
  		return amplitudeAntRMS;
 	}

 	/**
 	 * Set amplitudeAntRMS with the specified vector<float >.
 	 * @param amplitudeAntRMS The vector<float > value to which amplitudeAntRMS is to be set.
 	 
 	
 		
 	 */
 	void CalAntennaSolutionsRow::setAmplitudeAntRMS (vector<float > amplitudeAntRMS)  {
  	
  	
  		if (hasBeenAdded) {
 		
  		}
  	
 		this->amplitudeAntRMS = amplitudeAntRMS;
	
 	}
	
	

	
	///////////////////////////////////////////////
	// Extrinsic Table Attributes getters/setters//
	///////////////////////////////////////////////
	
	

	
 	/**
 	 * Get calDataId.
 	 * @return calDataId as Tag
 	 */
 	Tag CalAntennaSolutionsRow::getCalDataId() const {
	
  		return calDataId;
 	}

 	/**
 	 * Set calDataId with the specified Tag.
 	 * @param calDataId The Tag value to which calDataId is to be set.
 	 
 	
 		
 	 * @throw IllegalAccessException If an attempt is made to change this field after is has been added to the table.
 	 	
 	 */
 	void CalAntennaSolutionsRow::setCalDataId (Tag calDataId)  {
  	
  	
  		if (hasBeenAdded) {
 		
			throw IllegalAccessException("calDataId", "CalAntennaSolutions");
		
  		}
  	
 		this->calDataId = calDataId;
	
 	}
	
	

	

	
 	/**
 	 * Get calReductionId.
 	 * @return calReductionId as Tag
 	 */
 	Tag CalAntennaSolutionsRow::getCalReductionId() const {
	
  		return calReductionId;
 	}

 	/**
 	 * Set calReductionId with the specified Tag.
 	 * @param calReductionId The Tag value to which calReductionId is to be set.
 	 
 	
 		
 	 * @throw IllegalAccessException If an attempt is made to change this field after is has been added to the table.
 	 	
 	 */
 	void CalAntennaSolutionsRow::setCalReductionId (Tag calReductionId)  {
  	
  	
  		if (hasBeenAdded) {
 		
			throw IllegalAccessException("calReductionId", "CalAntennaSolutions");
		
  		}
  	
 		this->calReductionId = calReductionId;
	
 	}
	
	


	//////////////////////////////////////
	// Links Attributes getters/setters //
	//////////////////////////////////////
	
	
	
	
		

	/**
	 * Returns the pointer to the row in the CalData table having CalData.calDataId == calDataId
	 * @return a CalDataRow*
	 * 
	 
	 */
	 CalDataRow* CalAntennaSolutionsRow::getCalDataUsingCalDataId() {
	 
	 	return table.getContainer().getCalData().getRowByKey(calDataId);
	 }
	 

	

	
	
	
		

	/**
	 * Returns the pointer to the row in the CalReduction table having CalReduction.calReductionId == calReductionId
	 * @return a CalReductionRow*
	 * 
	 
	 */
	 CalReductionRow* CalAntennaSolutionsRow::getCalReductionUsingCalReductionId() {
	 
	 	return table.getContainer().getCalReduction().getRowByKey(calReductionId);
	 }
	 

	

	
	/**
	 * Create a CalAntennaSolutionsRow.
	 * <p>
	 * This constructor is private because only the
	 * table can create rows.  All rows know the table
	 * to which they belong.
	 * @param table The table to which this row belongs.
	 */ 
	CalAntennaSolutionsRow::CalAntennaSolutionsRow (CalAntennaSolutionsTable &t) : table(t) {
		hasBeenAdded = false;
		
	
	

	

	

	

	

	

	

	

	

	

	

	

	

	

	

	

	

	
	

	

	
	
	
	

	
// This attribute is scalar and has an enumeration type. Let's initialize it to some valid value (the 1st of the enumeration).		
atmPhaseCorrection = CAtmPhaseCorrection::from_int(0);
	

	
// This attribute is scalar and has an enumeration type. Let's initialize it to some valid value (the 1st of the enumeration).		
basebandName = CBasebandName::from_int(0);
	

	
// This attribute is scalar and has an enumeration type. Let's initialize it to some valid value (the 1st of the enumeration).		
receiverBand = CReceiverBand::from_int(0);
	

	

	

	

	

	

	

	

	

	

	

	

	

	

	
	
	 fromBinMethods["antennaName"] = &CalAntennaSolutionsRow::antennaNameFromBin; 
	 fromBinMethods["atmPhaseCorrection"] = &CalAntennaSolutionsRow::atmPhaseCorrectionFromBin; 
	 fromBinMethods["receiverBand"] = &CalAntennaSolutionsRow::receiverBandFromBin; 
	 fromBinMethods["basebandName"] = &CalAntennaSolutionsRow::basebandNameFromBin; 
	 fromBinMethods["calDataId"] = &CalAntennaSolutionsRow::calDataIdFromBin; 
	 fromBinMethods["calReductionId"] = &CalAntennaSolutionsRow::calReductionIdFromBin; 
	 fromBinMethods["startValidTime"] = &CalAntennaSolutionsRow::startValidTimeFromBin; 
	 fromBinMethods["endValidTime"] = &CalAntennaSolutionsRow::endValidTimeFromBin; 
	 fromBinMethods["numReceptor"] = &CalAntennaSolutionsRow::numReceptorFromBin; 
	 fromBinMethods["refAntennaName"] = &CalAntennaSolutionsRow::refAntennaNameFromBin; 
	 fromBinMethods["direction"] = &CalAntennaSolutionsRow::directionFromBin; 
	 fromBinMethods["frequencyRange"] = &CalAntennaSolutionsRow::frequencyRangeFromBin; 
	 fromBinMethods["integrationTime"] = &CalAntennaSolutionsRow::integrationTimeFromBin; 
	 fromBinMethods["polarizationTypes"] = &CalAntennaSolutionsRow::polarizationTypesFromBin; 
	 fromBinMethods["correctionValidity"] = &CalAntennaSolutionsRow::correctionValidityFromBin; 
	 fromBinMethods["phaseAnt"] = &CalAntennaSolutionsRow::phaseAntFromBin; 
	 fromBinMethods["phaseAntRMS"] = &CalAntennaSolutionsRow::phaseAntRMSFromBin; 
	 fromBinMethods["amplitudeAnt"] = &CalAntennaSolutionsRow::amplitudeAntFromBin; 
	 fromBinMethods["amplitudeAntRMS"] = &CalAntennaSolutionsRow::amplitudeAntRMSFromBin; 
		
	
	
	
	
	
				 
	fromTextMethods["antennaName"] = &CalAntennaSolutionsRow::antennaNameFromText;
		 
	
				 
	fromTextMethods["atmPhaseCorrection"] = &CalAntennaSolutionsRow::atmPhaseCorrectionFromText;
		 
	
				 
	fromTextMethods["receiverBand"] = &CalAntennaSolutionsRow::receiverBandFromText;
		 
	
				 
	fromTextMethods["basebandName"] = &CalAntennaSolutionsRow::basebandNameFromText;
		 
	
				 
	fromTextMethods["calDataId"] = &CalAntennaSolutionsRow::calDataIdFromText;
		 
	
				 
	fromTextMethods["calReductionId"] = &CalAntennaSolutionsRow::calReductionIdFromText;
		 
	
				 
	fromTextMethods["startValidTime"] = &CalAntennaSolutionsRow::startValidTimeFromText;
		 
	
				 
	fromTextMethods["endValidTime"] = &CalAntennaSolutionsRow::endValidTimeFromText;
		 
	
				 
	fromTextMethods["numReceptor"] = &CalAntennaSolutionsRow::numReceptorFromText;
		 
	
				 
	fromTextMethods["refAntennaName"] = &CalAntennaSolutionsRow::refAntennaNameFromText;
		 
	
				 
	fromTextMethods["direction"] = &CalAntennaSolutionsRow::directionFromText;
		 
	
				 
	fromTextMethods["frequencyRange"] = &CalAntennaSolutionsRow::frequencyRangeFromText;
		 
	
				 
	fromTextMethods["integrationTime"] = &CalAntennaSolutionsRow::integrationTimeFromText;
		 
	
				 
	fromTextMethods["polarizationTypes"] = &CalAntennaSolutionsRow::polarizationTypesFromText;
		 
	
				 
	fromTextMethods["correctionValidity"] = &CalAntennaSolutionsRow::correctionValidityFromText;
		 
	
				 
	fromTextMethods["phaseAnt"] = &CalAntennaSolutionsRow::phaseAntFromText;
		 
	
				 
	fromTextMethods["phaseAntRMS"] = &CalAntennaSolutionsRow::phaseAntRMSFromText;
		 
	
				 
	fromTextMethods["amplitudeAnt"] = &CalAntennaSolutionsRow::amplitudeAntFromText;
		 
	
				 
	fromTextMethods["amplitudeAntRMS"] = &CalAntennaSolutionsRow::amplitudeAntRMSFromText;
		 
	

		
	}
	
	CalAntennaSolutionsRow::CalAntennaSolutionsRow (CalAntennaSolutionsTable &t, CalAntennaSolutionsRow &row) : table(t) {
		hasBeenAdded = false;
		
		if (&row == 0) {
	
	
	

	

	

	

	

	

	

	

	

	

	

	

	

	

	

	

	

	
	

	
		
		}
		else {
	
		
			antennaName = row.antennaName;
		
			atmPhaseCorrection = row.atmPhaseCorrection;
		
			receiverBand = row.receiverBand;
		
			basebandName = row.basebandName;
		
			calDataId = row.calDataId;
		
			calReductionId = row.calReductionId;
		
		
		
		
			startValidTime = row.startValidTime;
		
			endValidTime = row.endValidTime;
		
			numReceptor = row.numReceptor;
		
			refAntennaName = row.refAntennaName;
		
			direction = row.direction;
		
			frequencyRange = row.frequencyRange;
		
			integrationTime = row.integrationTime;
		
			polarizationTypes = row.polarizationTypes;
		
			correctionValidity = row.correctionValidity;
		
			phaseAnt = row.phaseAnt;
		
			phaseAntRMS = row.phaseAntRMS;
		
			amplitudeAnt = row.amplitudeAnt;
		
			amplitudeAntRMS = row.amplitudeAntRMS;
		
		
		
		
		}
		
		 fromBinMethods["antennaName"] = &CalAntennaSolutionsRow::antennaNameFromBin; 
		 fromBinMethods["atmPhaseCorrection"] = &CalAntennaSolutionsRow::atmPhaseCorrectionFromBin; 
		 fromBinMethods["receiverBand"] = &CalAntennaSolutionsRow::receiverBandFromBin; 
		 fromBinMethods["basebandName"] = &CalAntennaSolutionsRow::basebandNameFromBin; 
		 fromBinMethods["calDataId"] = &CalAntennaSolutionsRow::calDataIdFromBin; 
		 fromBinMethods["calReductionId"] = &CalAntennaSolutionsRow::calReductionIdFromBin; 
		 fromBinMethods["startValidTime"] = &CalAntennaSolutionsRow::startValidTimeFromBin; 
		 fromBinMethods["endValidTime"] = &CalAntennaSolutionsRow::endValidTimeFromBin; 
		 fromBinMethods["numReceptor"] = &CalAntennaSolutionsRow::numReceptorFromBin; 
		 fromBinMethods["refAntennaName"] = &CalAntennaSolutionsRow::refAntennaNameFromBin; 
		 fromBinMethods["direction"] = &CalAntennaSolutionsRow::directionFromBin; 
		 fromBinMethods["frequencyRange"] = &CalAntennaSolutionsRow::frequencyRangeFromBin; 
		 fromBinMethods["integrationTime"] = &CalAntennaSolutionsRow::integrationTimeFromBin; 
		 fromBinMethods["polarizationTypes"] = &CalAntennaSolutionsRow::polarizationTypesFromBin; 
		 fromBinMethods["correctionValidity"] = &CalAntennaSolutionsRow::correctionValidityFromBin; 
		 fromBinMethods["phaseAnt"] = &CalAntennaSolutionsRow::phaseAntFromBin; 
		 fromBinMethods["phaseAntRMS"] = &CalAntennaSolutionsRow::phaseAntRMSFromBin; 
		 fromBinMethods["amplitudeAnt"] = &CalAntennaSolutionsRow::amplitudeAntFromBin; 
		 fromBinMethods["amplitudeAntRMS"] = &CalAntennaSolutionsRow::amplitudeAntRMSFromBin; 
			
	
			
	}

	
	bool CalAntennaSolutionsRow::compareNoAutoInc(string antennaName, AtmPhaseCorrectionMod::AtmPhaseCorrection atmPhaseCorrection, ReceiverBandMod::ReceiverBand receiverBand, BasebandNameMod::BasebandName basebandName, Tag calDataId, Tag calReductionId, ArrayTime startValidTime, ArrayTime endValidTime, int numReceptor, string refAntennaName, vector<Angle > direction, vector<Frequency > frequencyRange, Interval integrationTime, vector<PolarizationTypeMod::PolarizationType > polarizationTypes, bool correctionValidity, vector<float > phaseAnt, vector<float > phaseAntRMS, vector<float > amplitudeAnt, vector<float > amplitudeAntRMS) {
		bool result;
		result = true;
		
	
		
		result = result && (this->antennaName == antennaName);
		
		if (!result) return false;
	

	
		
		result = result && (this->atmPhaseCorrection == atmPhaseCorrection);
		
		if (!result) return false;
	

	
		
		result = result && (this->receiverBand == receiverBand);
		
		if (!result) return false;
	

	
		
		result = result && (this->basebandName == basebandName);
		
		if (!result) return false;
	

	
		
		result = result && (this->calDataId == calDataId);
		
		if (!result) return false;
	

	
		
		result = result && (this->calReductionId == calReductionId);
		
		if (!result) return false;
	

	
		
		result = result && (this->startValidTime == startValidTime);
		
		if (!result) return false;
	

	
		
		result = result && (this->endValidTime == endValidTime);
		
		if (!result) return false;
	

	
		
		result = result && (this->numReceptor == numReceptor);
		
		if (!result) return false;
	

	
		
		result = result && (this->refAntennaName == refAntennaName);
		
		if (!result) return false;
	

	
		
		result = result && (this->direction == direction);
		
		if (!result) return false;
	

	
		
		result = result && (this->frequencyRange == frequencyRange);
		
		if (!result) return false;
	

	
		
		result = result && (this->integrationTime == integrationTime);
		
		if (!result) return false;
	

	
		
		result = result && (this->polarizationTypes == polarizationTypes);
		
		if (!result) return false;
	

	
		
		result = result && (this->correctionValidity == correctionValidity);
		
		if (!result) return false;
	

	
		
		result = result && (this->phaseAnt == phaseAnt);
		
		if (!result) return false;
	

	
		
		result = result && (this->phaseAntRMS == phaseAntRMS);
		
		if (!result) return false;
	

	
		
		result = result && (this->amplitudeAnt == amplitudeAnt);
		
		if (!result) return false;
	

	
		
		result = result && (this->amplitudeAntRMS == amplitudeAntRMS);
		
		if (!result) return false;
	

		return result;
	}	
	
	
	
	bool CalAntennaSolutionsRow::compareRequiredValue(ArrayTime startValidTime, ArrayTime endValidTime, int numReceptor, string refAntennaName, vector<Angle > direction, vector<Frequency > frequencyRange, Interval integrationTime, vector<PolarizationTypeMod::PolarizationType > polarizationTypes, bool correctionValidity, vector<float > phaseAnt, vector<float > phaseAntRMS, vector<float > amplitudeAnt, vector<float > amplitudeAntRMS) {
		bool result;
		result = true;
		
	
		if (!(this->startValidTime == startValidTime)) return false;
	

	
		if (!(this->endValidTime == endValidTime)) return false;
	

	
		if (!(this->numReceptor == numReceptor)) return false;
	

	
		if (!(this->refAntennaName == refAntennaName)) return false;
	

	
		if (!(this->direction == direction)) return false;
	

	
		if (!(this->frequencyRange == frequencyRange)) return false;
	

	
		if (!(this->integrationTime == integrationTime)) return false;
	

	
		if (!(this->polarizationTypes == polarizationTypes)) return false;
	

	
		if (!(this->correctionValidity == correctionValidity)) return false;
	

	
		if (!(this->phaseAnt == phaseAnt)) return false;
	

	
		if (!(this->phaseAntRMS == phaseAntRMS)) return false;
	

	
		if (!(this->amplitudeAnt == amplitudeAnt)) return false;
	

	
		if (!(this->amplitudeAntRMS == amplitudeAntRMS)) return false;
	

		return result;
	}
	
	
	/**
	 * Return true if all required attributes of the value part are equal to their homologues
	 * in x and false otherwise.
	 *

	 * @param x a pointer on the CalAntennaSolutionsRow whose required attributes of the value part 

	 * will be compared with those of this.
	 * @return a boolean.
	 */
	bool CalAntennaSolutionsRow::equalByRequiredValue(CalAntennaSolutionsRow*  x ) {
		
			
		if (this->startValidTime != x->startValidTime) return false;
			
		if (this->endValidTime != x->endValidTime) return false;
			
		if (this->numReceptor != x->numReceptor) return false;
			
		if (this->refAntennaName != x->refAntennaName) return false;
			
		if (this->direction != x->direction) return false;
			
		if (this->frequencyRange != x->frequencyRange) return false;
			
		if (this->integrationTime != x->integrationTime) return false;
			
		if (this->polarizationTypes != x->polarizationTypes) return false;
			
		if (this->correctionValidity != x->correctionValidity) return false;
			
		if (this->phaseAnt != x->phaseAnt) return false;
			
		if (this->phaseAntRMS != x->phaseAntRMS) return false;
			
		if (this->amplitudeAnt != x->amplitudeAnt) return false;
			
		if (this->amplitudeAntRMS != x->amplitudeAntRMS) return false;
			
		
		return true;
	}	
	
/*
	 map<string, CalAntennaSolutionsAttributeFromBin> CalAntennaSolutionsRow::initFromBinMethods() {
		map<string, CalAntennaSolutionsAttributeFromBin> result;
		
		result["antennaName"] = &CalAntennaSolutionsRow::antennaNameFromBin;
		result["atmPhaseCorrection"] = &CalAntennaSolutionsRow::atmPhaseCorrectionFromBin;
		result["receiverBand"] = &CalAntennaSolutionsRow::receiverBandFromBin;
		result["basebandName"] = &CalAntennaSolutionsRow::basebandNameFromBin;
		result["calDataId"] = &CalAntennaSolutionsRow::calDataIdFromBin;
		result["calReductionId"] = &CalAntennaSolutionsRow::calReductionIdFromBin;
		result["startValidTime"] = &CalAntennaSolutionsRow::startValidTimeFromBin;
		result["endValidTime"] = &CalAntennaSolutionsRow::endValidTimeFromBin;
		result["numReceptor"] = &CalAntennaSolutionsRow::numReceptorFromBin;
		result["refAntennaName"] = &CalAntennaSolutionsRow::refAntennaNameFromBin;
		result["direction"] = &CalAntennaSolutionsRow::directionFromBin;
		result["frequencyRange"] = &CalAntennaSolutionsRow::frequencyRangeFromBin;
		result["integrationTime"] = &CalAntennaSolutionsRow::integrationTimeFromBin;
		result["polarizationTypes"] = &CalAntennaSolutionsRow::polarizationTypesFromBin;
		result["correctionValidity"] = &CalAntennaSolutionsRow::correctionValidityFromBin;
		result["phaseAnt"] = &CalAntennaSolutionsRow::phaseAntFromBin;
		result["phaseAntRMS"] = &CalAntennaSolutionsRow::phaseAntRMSFromBin;
		result["amplitudeAnt"] = &CalAntennaSolutionsRow::amplitudeAntFromBin;
		result["amplitudeAntRMS"] = &CalAntennaSolutionsRow::amplitudeAntRMSFromBin;
		
		
			
		
		return result;	
	}
*/	
} // End namespace asdm
 
