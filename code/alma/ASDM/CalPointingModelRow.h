
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
 * File CalPointingModelRow.h
 */
 
#ifndef CalPointingModelRow_CLASS
#define CalPointingModelRow_CLASS

#include <vector>
#include <string>
#include <set>

#ifndef WITHOUT_ACS
#include <asdmIDLC.h>
#endif






	 
#include <Angle.h>
	

	 
#include <ArrayTime.h>
	

	 
#include <Tag.h>
	




	

	
#include "CReceiverBand.h"
	

	

	

	
#include "CAntennaMake.h"
	

	
#include "CPointingModelMode.h"
	

	
#include "CPolarizationType.h"
	

	

	

	

	

	

	

	

	

	

	

	



#include <ConversionException.h>
#include <NoSuchRow.h>
#include <IllegalAccessException.h>

#include <RowTransformer.h>
//#include <TableStreamReader.h>

/*\file CalPointingModel.h
    \brief Generated from model's revision "-1", branch ""
*/

namespace asdm {

//class asdm::CalPointingModelTable;


// class asdm::CalDataRow;
class CalDataRow;

// class asdm::CalReductionRow;
class CalReductionRow;
	

class CalPointingModelRow;
typedef void (CalPointingModelRow::*CalPointingModelAttributeFromBin) (EndianIStream& eis);
typedef void (CalPointingModelRow::*CalPointingModelAttributeFromText) (const string& s);

/**
 * The CalPointingModelRow class is a row of a CalPointingModelTable.
 * 
 * Generated from model's revision "-1", branch ""
 *
 */
class CalPointingModelRow {
friend class asdm::CalPointingModelTable;
friend class asdm::RowTransformer<CalPointingModelRow>;
//friend class asdm::TableStreamReader<CalPointingModelTable, CalPointingModelRow>;

public:

	virtual ~CalPointingModelRow();

	/**
	 * Return the table to which this row belongs.
	 */
	CalPointingModelTable &getTable() const;
	
	/**
	 * Has this row been added to its table ?
	 * @return true if and only if it has been added.
	 */
	bool isAdded() const;
		
	////////////////////////////////
	// Intrinsic Table Attributes //
	////////////////////////////////
	
	
	// ===> Attribute antennaName
	
	
	

	
 	/**
 	 * Get antennaName.
 	 * @return antennaName as string
 	 */
 	string getAntennaName() const;
	
 
 	
 	
 	/**
 	 * Set antennaName with the specified string.
 	 * @param antennaName The string value to which antennaName is to be set.
 	 
 		
 			
 	 * @throw IllegalAccessException If an attempt is made to change this field after is has been added to the table.
 	 		
 	 */
 	void setAntennaName (string antennaName);
  		
	
	
	


	
	// ===> Attribute receiverBand
	
	
	

	
 	/**
 	 * Get receiverBand.
 	 * @return receiverBand as ReceiverBandMod::ReceiverBand
 	 */
 	ReceiverBandMod::ReceiverBand getReceiverBand() const;
	
 
 	
 	
 	/**
 	 * Set receiverBand with the specified ReceiverBandMod::ReceiverBand.
 	 * @param receiverBand The ReceiverBandMod::ReceiverBand value to which receiverBand is to be set.
 	 
 		
 			
 	 * @throw IllegalAccessException If an attempt is made to change this field after is has been added to the table.
 	 		
 	 */
 	void setReceiverBand (ReceiverBandMod::ReceiverBand receiverBand);
  		
	
	
	


	
	// ===> Attribute startValidTime
	
	
	

	
 	/**
 	 * Get startValidTime.
 	 * @return startValidTime as ArrayTime
 	 */
 	ArrayTime getStartValidTime() const;
	
 
 	
 	
 	/**
 	 * Set startValidTime with the specified ArrayTime.
 	 * @param startValidTime The ArrayTime value to which startValidTime is to be set.
 	 
 		
 			
 	 */
 	void setStartValidTime (ArrayTime startValidTime);
  		
	
	
	


	
	// ===> Attribute endValidTime
	
	
	

	
 	/**
 	 * Get endValidTime.
 	 * @return endValidTime as ArrayTime
 	 */
 	ArrayTime getEndValidTime() const;
	
 
 	
 	
 	/**
 	 * Set endValidTime with the specified ArrayTime.
 	 * @param endValidTime The ArrayTime value to which endValidTime is to be set.
 	 
 		
 			
 	 */
 	void setEndValidTime (ArrayTime endValidTime);
  		
	
	
	


	
	// ===> Attribute antennaMake
	
	
	

	
 	/**
 	 * Get antennaMake.
 	 * @return antennaMake as AntennaMakeMod::AntennaMake
 	 */
 	AntennaMakeMod::AntennaMake getAntennaMake() const;
	
 
 	
 	
 	/**
 	 * Set antennaMake with the specified AntennaMakeMod::AntennaMake.
 	 * @param antennaMake The AntennaMakeMod::AntennaMake value to which antennaMake is to be set.
 	 
 		
 			
 	 */
 	void setAntennaMake (AntennaMakeMod::AntennaMake antennaMake);
  		
	
	
	


	
	// ===> Attribute pointingModelMode
	
	
	

	
 	/**
 	 * Get pointingModelMode.
 	 * @return pointingModelMode as PointingModelModeMod::PointingModelMode
 	 */
 	PointingModelModeMod::PointingModelMode getPointingModelMode() const;
	
 
 	
 	
 	/**
 	 * Set pointingModelMode with the specified PointingModelModeMod::PointingModelMode.
 	 * @param pointingModelMode The PointingModelModeMod::PointingModelMode value to which pointingModelMode is to be set.
 	 
 		
 			
 	 */
 	void setPointingModelMode (PointingModelModeMod::PointingModelMode pointingModelMode);
  		
	
	
	


	
	// ===> Attribute polarizationType
	
	
	

	
 	/**
 	 * Get polarizationType.
 	 * @return polarizationType as PolarizationTypeMod::PolarizationType
 	 */
 	PolarizationTypeMod::PolarizationType getPolarizationType() const;
	
 
 	
 	
 	/**
 	 * Set polarizationType with the specified PolarizationTypeMod::PolarizationType.
 	 * @param polarizationType The PolarizationTypeMod::PolarizationType value to which polarizationType is to be set.
 	 
 		
 			
 	 */
 	void setPolarizationType (PolarizationTypeMod::PolarizationType polarizationType);
  		
	
	
	


	
	// ===> Attribute numCoeff
	
	
	

	
 	/**
 	 * Get numCoeff.
 	 * @return numCoeff as int
 	 */
 	int getNumCoeff() const;
	
 
 	
 	
 	/**
 	 * Set numCoeff with the specified int.
 	 * @param numCoeff The int value to which numCoeff is to be set.
 	 
 		
 			
 	 */
 	void setNumCoeff (int numCoeff);
  		
	
	
	


	
	// ===> Attribute coeffName
	
	
	

	
 	/**
 	 * Get coeffName.
 	 * @return coeffName as vector<string >
 	 */
 	vector<string > getCoeffName() const;
	
 
 	
 	
 	/**
 	 * Set coeffName with the specified vector<string >.
 	 * @param coeffName The vector<string > value to which coeffName is to be set.
 	 
 		
 			
 	 */
 	void setCoeffName (vector<string > coeffName);
  		
	
	
	


	
	// ===> Attribute coeffVal
	
	
	

	
 	/**
 	 * Get coeffVal.
 	 * @return coeffVal as vector<float >
 	 */
 	vector<float > getCoeffVal() const;
	
 
 	
 	
 	/**
 	 * Set coeffVal with the specified vector<float >.
 	 * @param coeffVal The vector<float > value to which coeffVal is to be set.
 	 
 		
 			
 	 */
 	void setCoeffVal (vector<float > coeffVal);
  		
	
	
	


	
	// ===> Attribute coeffError
	
	
	

	
 	/**
 	 * Get coeffError.
 	 * @return coeffError as vector<float >
 	 */
 	vector<float > getCoeffError() const;
	
 
 	
 	
 	/**
 	 * Set coeffError with the specified vector<float >.
 	 * @param coeffError The vector<float > value to which coeffError is to be set.
 	 
 		
 			
 	 */
 	void setCoeffError (vector<float > coeffError);
  		
	
	
	


	
	// ===> Attribute coeffFixed
	
	
	

	
 	/**
 	 * Get coeffFixed.
 	 * @return coeffFixed as vector<bool >
 	 */
 	vector<bool > getCoeffFixed() const;
	
 
 	
 	
 	/**
 	 * Set coeffFixed with the specified vector<bool >.
 	 * @param coeffFixed The vector<bool > value to which coeffFixed is to be set.
 	 
 		
 			
 	 */
 	void setCoeffFixed (vector<bool > coeffFixed);
  		
	
	
	


	
	// ===> Attribute azimuthRMS
	
	
	

	
 	/**
 	 * Get azimuthRMS.
 	 * @return azimuthRMS as Angle
 	 */
 	Angle getAzimuthRMS() const;
	
 
 	
 	
 	/**
 	 * Set azimuthRMS with the specified Angle.
 	 * @param azimuthRMS The Angle value to which azimuthRMS is to be set.
 	 
 		
 			
 	 */
 	void setAzimuthRMS (Angle azimuthRMS);
  		
	
	
	


	
	// ===> Attribute elevationRms
	
	
	

	
 	/**
 	 * Get elevationRms.
 	 * @return elevationRms as Angle
 	 */
 	Angle getElevationRms() const;
	
 
 	
 	
 	/**
 	 * Set elevationRms with the specified Angle.
 	 * @param elevationRms The Angle value to which elevationRms is to be set.
 	 
 		
 			
 	 */
 	void setElevationRms (Angle elevationRms);
  		
	
	
	


	
	// ===> Attribute skyRMS
	
	
	

	
 	/**
 	 * Get skyRMS.
 	 * @return skyRMS as Angle
 	 */
 	Angle getSkyRMS() const;
	
 
 	
 	
 	/**
 	 * Set skyRMS with the specified Angle.
 	 * @param skyRMS The Angle value to which skyRMS is to be set.
 	 
 		
 			
 	 */
 	void setSkyRMS (Angle skyRMS);
  		
	
	
	


	
	// ===> Attribute reducedChiSquared
	
	
	

	
 	/**
 	 * Get reducedChiSquared.
 	 * @return reducedChiSquared as double
 	 */
 	double getReducedChiSquared() const;
	
 
 	
 	
 	/**
 	 * Set reducedChiSquared with the specified double.
 	 * @param reducedChiSquared The double value to which reducedChiSquared is to be set.
 	 
 		
 			
 	 */
 	void setReducedChiSquared (double reducedChiSquared);
  		
	
	
	


	
	// ===> Attribute numObs, which is optional
	
	
	
	/**
	 * The attribute numObs is optional. Return true if this attribute exists.
	 * @return true if and only if the numObs attribute exists. 
	 */
	bool isNumObsExists() const;
	

	
 	/**
 	 * Get numObs, which is optional.
 	 * @return numObs as int
 	 * @throws IllegalAccessException If numObs does not exist.
 	 */
 	int getNumObs() const;
	
 
 	
 	
 	/**
 	 * Set numObs with the specified int.
 	 * @param numObs The int value to which numObs is to be set.
 	 
 		
 	 */
 	void setNumObs (int numObs);
		
	
	
	
	/**
	 * Mark numObs, which is an optional field, as non-existent.
	 */
	void clearNumObs ();
	


	
	// ===> Attribute coeffFormula, which is optional
	
	
	
	/**
	 * The attribute coeffFormula is optional. Return true if this attribute exists.
	 * @return true if and only if the coeffFormula attribute exists. 
	 */
	bool isCoeffFormulaExists() const;
	

	
 	/**
 	 * Get coeffFormula, which is optional.
 	 * @return coeffFormula as vector<string >
 	 * @throws IllegalAccessException If coeffFormula does not exist.
 	 */
 	vector<string > getCoeffFormula() const;
	
 
 	
 	
 	/**
 	 * Set coeffFormula with the specified vector<string >.
 	 * @param coeffFormula The vector<string > value to which coeffFormula is to be set.
 	 
 		
 	 */
 	void setCoeffFormula (vector<string > coeffFormula);
		
	
	
	
	/**
	 * Mark coeffFormula, which is an optional field, as non-existent.
	 */
	void clearCoeffFormula ();
	


	////////////////////////////////
	// Extrinsic Table Attributes //
	////////////////////////////////
	
	
	// ===> Attribute calDataId
	
	
	

	
 	/**
 	 * Get calDataId.
 	 * @return calDataId as Tag
 	 */
 	Tag getCalDataId() const;
	
 
 	
 	
 	/**
 	 * Set calDataId with the specified Tag.
 	 * @param calDataId The Tag value to which calDataId is to be set.
 	 
 		
 			
 	 * @throw IllegalAccessException If an attempt is made to change this field after is has been added to the table.
 	 		
 	 */
 	void setCalDataId (Tag calDataId);
  		
	
	
	


	
	// ===> Attribute calReductionId
	
	
	

	
 	/**
 	 * Get calReductionId.
 	 * @return calReductionId as Tag
 	 */
 	Tag getCalReductionId() const;
	
 
 	
 	
 	/**
 	 * Set calReductionId with the specified Tag.
 	 * @param calReductionId The Tag value to which calReductionId is to be set.
 	 
 		
 			
 	 * @throw IllegalAccessException If an attempt is made to change this field after is has been added to the table.
 	 		
 	 */
 	void setCalReductionId (Tag calReductionId);
  		
	
	
	


	///////////
	// Links //
	///////////
	
	

	
		
	/**
	 * calDataId pointer to the row in the CalData table having CalData.calDataId == calDataId
	 * @return a CalDataRow*
	 * 
	 
	 */
	 CalDataRow* getCalDataUsingCalDataId();
	 

	

	

	
		
	/**
	 * calReductionId pointer to the row in the CalReduction table having CalReduction.calReductionId == calReductionId
	 * @return a CalReductionRow*
	 * 
	 
	 */
	 CalReductionRow* getCalReductionUsingCalReductionId();
	 

	

	
	
	
	/**
	 * Compare each mandatory attribute except the autoincrementable one of this CalPointingModelRow with 
	 * the corresponding parameters and return true if there is a match and false otherwise.
	 	
	 * @param antennaName
	    
	 * @param receiverBand
	    
	 * @param calDataId
	    
	 * @param calReductionId
	    
	 * @param startValidTime
	    
	 * @param endValidTime
	    
	 * @param antennaMake
	    
	 * @param pointingModelMode
	    
	 * @param polarizationType
	    
	 * @param numCoeff
	    
	 * @param coeffName
	    
	 * @param coeffVal
	    
	 * @param coeffError
	    
	 * @param coeffFixed
	    
	 * @param azimuthRMS
	    
	 * @param elevationRms
	    
	 * @param skyRMS
	    
	 * @param reducedChiSquared
	    
	 */ 
	bool compareNoAutoInc(string antennaName, ReceiverBandMod::ReceiverBand receiverBand, Tag calDataId, Tag calReductionId, ArrayTime startValidTime, ArrayTime endValidTime, AntennaMakeMod::AntennaMake antennaMake, PointingModelModeMod::PointingModelMode pointingModelMode, PolarizationTypeMod::PolarizationType polarizationType, int numCoeff, vector<string > coeffName, vector<float > coeffVal, vector<float > coeffError, vector<bool > coeffFixed, Angle azimuthRMS, Angle elevationRms, Angle skyRMS, double reducedChiSquared);
	
	

	
	/**
	 * Compare each mandatory value (i.e. not in the key) attribute  with 
	 * the corresponding parameters and return true if there is a match and false otherwise.
	 	
	 * @param startValidTime
	    
	 * @param endValidTime
	    
	 * @param antennaMake
	    
	 * @param pointingModelMode
	    
	 * @param polarizationType
	    
	 * @param numCoeff
	    
	 * @param coeffName
	    
	 * @param coeffVal
	    
	 * @param coeffError
	    
	 * @param coeffFixed
	    
	 * @param azimuthRMS
	    
	 * @param elevationRms
	    
	 * @param skyRMS
	    
	 * @param reducedChiSquared
	    
	 */ 
	bool compareRequiredValue(ArrayTime startValidTime, ArrayTime endValidTime, AntennaMakeMod::AntennaMake antennaMake, PointingModelModeMod::PointingModelMode pointingModelMode, PolarizationTypeMod::PolarizationType polarizationType, int numCoeff, vector<string > coeffName, vector<float > coeffVal, vector<float > coeffError, vector<bool > coeffFixed, Angle azimuthRMS, Angle elevationRms, Angle skyRMS, double reducedChiSquared); 
		 
	
	/**
	 * Return true if all required attributes of the value part are equal to their homologues
	 * in x and false otherwise.
	 *
	 * @param x a pointer on the CalPointingModelRow whose required attributes of the value part 
	 * will be compared with those of this.
	 * @return a boolean.
	 */
	bool equalByRequiredValue(CalPointingModelRow* x) ;
	
#ifndef WITHOUT_ACS
	/**
	 * Return this row in the form of an IDL struct.
	 * @return The values of this row as a CalPointingModelRowIDL struct.
	 */
	asdmIDL::CalPointingModelRowIDL *toIDL() const;
	
	/**
	 * Define the content of a CalPointingModelRowIDL struct from the values
	 * found in this row.
	 *
	 * @param x a reference to the CalPointingModelRowIDL struct to be set.
	 *
	 */
	 void toIDL(asdmIDL::CalPointingModelRowIDL& x) const;
#endif
	
#ifndef WITHOUT_ACS
	/**
	 * Fill the values of this row from the IDL struct CalPointingModelRowIDL.
	 * @param x The IDL struct containing the values used to fill this row.
	 * @throws ConversionException
	 */
	void setFromIDL (asdmIDL::CalPointingModelRowIDL x) ;
#endif
	
	/**
	 * Return this row in the form of an XML string.
	 * @return The values of this row as an XML string.
	 */
	std::string toXML() const;

	/**
	 * Fill the values of this row from an XML string 
	 * that was produced by the toXML() method.
	 * @param rowDoc the XML string being used to set the values of this row.
	 * @throws ConversionException
	 */
	void setFromXML (std::string rowDoc) ;

	/// @cond DISPLAY_PRIVATE	
	////////////////////////////////////////////////////////////
	// binary-deserialization material from an EndianIStream  //
	////////////////////////////////////////////////////////////

	std::map<std::string, CalPointingModelAttributeFromBin> fromBinMethods;
void antennaNameFromBin( EndianIStream& eis);
void receiverBandFromBin( EndianIStream& eis);
void calDataIdFromBin( EndianIStream& eis);
void calReductionIdFromBin( EndianIStream& eis);
void startValidTimeFromBin( EndianIStream& eis);
void endValidTimeFromBin( EndianIStream& eis);
void antennaMakeFromBin( EndianIStream& eis);
void pointingModelModeFromBin( EndianIStream& eis);
void polarizationTypeFromBin( EndianIStream& eis);
void numCoeffFromBin( EndianIStream& eis);
void coeffNameFromBin( EndianIStream& eis);
void coeffValFromBin( EndianIStream& eis);
void coeffErrorFromBin( EndianIStream& eis);
void coeffFixedFromBin( EndianIStream& eis);
void azimuthRMSFromBin( EndianIStream& eis);
void elevationRmsFromBin( EndianIStream& eis);
void skyRMSFromBin( EndianIStream& eis);
void reducedChiSquaredFromBin( EndianIStream& eis);

void numObsFromBin( EndianIStream& eis);
void coeffFormulaFromBin( EndianIStream& eis);


	 /**
	  * Deserialize a stream of bytes read from an EndianIStream to build a PointingRow.
	  * @param eiss the EndianIStream to be read.
	  * @param table the CalPointingModelTable to which the row built by deserialization will be parented.
	  * @param attributesSeq a vector containing the names of the attributes . The elements order defines the order 
	  * in which the attributes are written in the binary serialization.
	  */
	 static CalPointingModelRow* fromBin(EndianIStream& eis, CalPointingModelTable& table, const std::vector<std::string>& attributesSeq);	 
 
 	 /**
 	  * Parses a string t and assign the result of the parsing to the attribute of name attributeName.
 	  *
 	  * @param attributeName the name of the attribute whose value is going to be defined.
 	  * @param t the string to be parsed into a value given to the attribute of name attributeName.
 	  */
 	 void fromText(const std::string& attributeName, const std::string&  t);
     /// @endcond			

private:
	/**
	 * The table to which this row belongs.
	 */
	CalPointingModelTable &table;
	/**
	 * Whether this row has been added to the table or not.
	 */
	bool hasBeenAdded;

	// This method is used by the Table class when this row is added to the table.
	void isAdded(bool added);


	/**
	 * Create a CalPointingModelRow.
	 * <p>
	 * This constructor is private because only the
	 * table can create rows.  All rows know the table
	 * to which they belong.
	 * @param table The table to which this row belongs.
	 */ 
	CalPointingModelRow (CalPointingModelTable &table);

	/**
	 * Create a CalPointingModelRow using a copy constructor mechanism.
	 * <p>
	 * Given a CalPointingModelRow row and a CalPointingModelTable table, the method creates a new
	 * CalPointingModelRow owned by table. Each attribute of the created row is a copy (deep)
	 * of the corresponding attribute of row. The method does not add the created
	 * row to its table, its simply parents it to table, a call to the add method
	 * has to be done in order to get the row added (very likely after having modified
	 * some of its attributes).
	 * If row is null then the method returns a row with default values for its attributes. 
	 *
	 * This constructor is private because only the
	 * table can create rows.  All rows know the table
	 * to which they belong.
	 * @param table The table to which this row belongs.
	 * @param row  The row which is to be copied.
	 */
	 CalPointingModelRow (CalPointingModelTable &table, CalPointingModelRow &row);
	 	
	////////////////////////////////
	// Intrinsic Table Attributes //
	////////////////////////////////
	
	
	// ===> Attribute antennaName
	
	

	string antennaName;

	
	
 	

	
	// ===> Attribute receiverBand
	
	

	ReceiverBandMod::ReceiverBand receiverBand;

	
	
 	

	
	// ===> Attribute startValidTime
	
	

	ArrayTime startValidTime;

	
	
 	

	
	// ===> Attribute endValidTime
	
	

	ArrayTime endValidTime;

	
	
 	

	
	// ===> Attribute antennaMake
	
	

	AntennaMakeMod::AntennaMake antennaMake;

	
	
 	

	
	// ===> Attribute pointingModelMode
	
	

	PointingModelModeMod::PointingModelMode pointingModelMode;

	
	
 	

	
	// ===> Attribute polarizationType
	
	

	PolarizationTypeMod::PolarizationType polarizationType;

	
	
 	

	
	// ===> Attribute numCoeff
	
	

	int numCoeff;

	
	
 	

	
	// ===> Attribute coeffName
	
	

	vector<string > coeffName;

	
	
 	

	
	// ===> Attribute coeffVal
	
	

	vector<float > coeffVal;

	
	
 	

	
	// ===> Attribute coeffError
	
	

	vector<float > coeffError;

	
	
 	

	
	// ===> Attribute coeffFixed
	
	

	vector<bool > coeffFixed;

	
	
 	

	
	// ===> Attribute azimuthRMS
	
	

	Angle azimuthRMS;

	
	
 	

	
	// ===> Attribute elevationRms
	
	

	Angle elevationRms;

	
	
 	

	
	// ===> Attribute skyRMS
	
	

	Angle skyRMS;

	
	
 	

	
	// ===> Attribute reducedChiSquared
	
	

	double reducedChiSquared;

	
	
 	

	
	// ===> Attribute numObs, which is optional
	
	
	bool numObsExists;
	

	int numObs;

	
	
 	

	
	// ===> Attribute coeffFormula, which is optional
	
	
	bool coeffFormulaExists;
	

	vector<string > coeffFormula;

	
	
 	

	////////////////////////////////
	// Extrinsic Table Attributes //
	////////////////////////////////
	
	
	// ===> Attribute calDataId
	
	

	Tag calDataId;

	
	
 	

	
	// ===> Attribute calReductionId
	
	

	Tag calReductionId;

	
	
 	

	///////////
	// Links //
	///////////
	
	
		

	 

	

	
		

	 

	

	
/*
	////////////////////////////////////////////////////////////
	// binary-deserialization material from an EndianIStream  //
	////////////////////////////////////////////////////////////
	std::map<std::string, CalPointingModelAttributeFromBin> fromBinMethods;
void antennaNameFromBin( EndianIStream& eis);
void receiverBandFromBin( EndianIStream& eis);
void calDataIdFromBin( EndianIStream& eis);
void calReductionIdFromBin( EndianIStream& eis);
void startValidTimeFromBin( EndianIStream& eis);
void endValidTimeFromBin( EndianIStream& eis);
void antennaMakeFromBin( EndianIStream& eis);
void pointingModelModeFromBin( EndianIStream& eis);
void polarizationTypeFromBin( EndianIStream& eis);
void numCoeffFromBin( EndianIStream& eis);
void coeffNameFromBin( EndianIStream& eis);
void coeffValFromBin( EndianIStream& eis);
void coeffErrorFromBin( EndianIStream& eis);
void coeffFixedFromBin( EndianIStream& eis);
void azimuthRMSFromBin( EndianIStream& eis);
void elevationRmsFromBin( EndianIStream& eis);
void skyRMSFromBin( EndianIStream& eis);
void reducedChiSquaredFromBin( EndianIStream& eis);

void numObsFromBin( EndianIStream& eis);
void coeffFormulaFromBin( EndianIStream& eis);

*/
	
	///////////////////////////////////
	// text-deserialization material //
	///////////////////////////////////
	std::map<std::string, CalPointingModelAttributeFromText> fromTextMethods;
	
void antennaNameFromText (const string & s);
	
	
void receiverBandFromText (const string & s);
	
	
void calDataIdFromText (const string & s);
	
	
void calReductionIdFromText (const string & s);
	
	
void startValidTimeFromText (const string & s);
	
	
void endValidTimeFromText (const string & s);
	
	
void antennaMakeFromText (const string & s);
	
	
void pointingModelModeFromText (const string & s);
	
	
void polarizationTypeFromText (const string & s);
	
	
void numCoeffFromText (const string & s);
	
	
void coeffNameFromText (const string & s);
	
	
void coeffValFromText (const string & s);
	
	
void coeffErrorFromText (const string & s);
	
	
void coeffFixedFromText (const string & s);
	
	
void azimuthRMSFromText (const string & s);
	
	
void elevationRmsFromText (const string & s);
	
	
void skyRMSFromText (const string & s);
	
	
void reducedChiSquaredFromText (const string & s);
	

	
void numObsFromText (const string & s);
	
	
void coeffFormulaFromText (const string & s);
	
	
	
	/**
	 * Serialize this into a stream of bytes written to an EndianOSStream.
	 * @param eoss the EndianOSStream to be written to
	 */
	 void toBin(EndianOSStream& eoss);
	 	 
	 /**
	  * Deserialize a stream of bytes read from an EndianIStream to build a PointingRow.
	  * @param eiss the EndianIStream to be read.
	  * @param table the CalPointingModelTable to which the row built by deserialization will be parented.
	  * @param attributesSeq a vector containing the names of the attributes . The elements order defines the order 
	  * in which the attributes are written in the binary serialization.

	 static CalPointingModelRow* fromBin(EndianIStream& eis, CalPointingModelTable& table, const std::vector<std::string>& attributesSeq);	 
		*/
};

} // End namespace asdm

#endif /* CalPointingModel_CLASS */
