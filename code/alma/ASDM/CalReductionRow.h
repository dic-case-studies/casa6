
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
 * File CalReductionRow.h
 */
 
#ifndef CalReductionRow_CLASS
#define CalReductionRow_CLASS

#include <vector>
#include <string>
#include <set>

#ifndef WITHOUT_ACS
#include <asdmIDLC.h>
#endif






	 
#include <ArrayTime.h>
	

	 
#include <Tag.h>
	




	

	

	

	

	

	

	
#include "CInvalidatingCondition.h"
	

	

	

	

	



#include <ConversionException.h>
#include <NoSuchRow.h>
#include <IllegalAccessException.h>

#include <RowTransformer.h>
//#include <TableStreamReader.h>

/*\file CalReduction.h
    \brief Generated from model's revision "-1", branch ""
*/

namespace asdm {

//class asdm::CalReductionTable;

	

class CalReductionRow;
typedef void (CalReductionRow::*CalReductionAttributeFromBin) (EndianIStream& eis);
typedef void (CalReductionRow::*CalReductionAttributeFromText) (const string& s);

/**
 * The CalReductionRow class is a row of a CalReductionTable.
 * 
 * Generated from model's revision "-1", branch ""
 *
 */
class CalReductionRow {
friend class asdm::CalReductionTable;
friend class asdm::RowTransformer<CalReductionRow>;
//friend class asdm::TableStreamReader<CalReductionTable, CalReductionRow>;

public:

	virtual ~CalReductionRow();

	/**
	 * Return the table to which this row belongs.
	 */
	CalReductionTable &getTable() const;
	
	/**
	 * Has this row been added to its table ?
	 * @return true if and only if it has been added.
	 */
	bool isAdded() const;
		
	////////////////////////////////
	// Intrinsic Table Attributes //
	////////////////////////////////
	
	
	// ===> Attribute calReductionId
	
	
	

	
 	/**
 	 * Get calReductionId.
 	 * @return calReductionId as Tag
 	 */
 	Tag getCalReductionId() const;
	
 
 	
 	
	
	


	
	// ===> Attribute numApplied
	
	
	

	
 	/**
 	 * Get numApplied.
 	 * @return numApplied as int
 	 */
 	int getNumApplied() const;
	
 
 	
 	
 	/**
 	 * Set numApplied with the specified int.
 	 * @param numApplied The int value to which numApplied is to be set.
 	 
 		
 			
 	 */
 	void setNumApplied (int numApplied);
  		
	
	
	


	
	// ===> Attribute appliedCalibrations
	
	
	

	
 	/**
 	 * Get appliedCalibrations.
 	 * @return appliedCalibrations as vector<string >
 	 */
 	vector<string > getAppliedCalibrations() const;
	
 
 	
 	
 	/**
 	 * Set appliedCalibrations with the specified vector<string >.
 	 * @param appliedCalibrations The vector<string > value to which appliedCalibrations is to be set.
 	 
 		
 			
 	 */
 	void setAppliedCalibrations (vector<string > appliedCalibrations);
  		
	
	
	


	
	// ===> Attribute numParam
	
	
	

	
 	/**
 	 * Get numParam.
 	 * @return numParam as int
 	 */
 	int getNumParam() const;
	
 
 	
 	
 	/**
 	 * Set numParam with the specified int.
 	 * @param numParam The int value to which numParam is to be set.
 	 
 		
 			
 	 */
 	void setNumParam (int numParam);
  		
	
	
	


	
	// ===> Attribute paramSet
	
	
	

	
 	/**
 	 * Get paramSet.
 	 * @return paramSet as vector<string >
 	 */
 	vector<string > getParamSet() const;
	
 
 	
 	
 	/**
 	 * Set paramSet with the specified vector<string >.
 	 * @param paramSet The vector<string > value to which paramSet is to be set.
 	 
 		
 			
 	 */
 	void setParamSet (vector<string > paramSet);
  		
	
	
	


	
	// ===> Attribute numInvalidConditions
	
	
	

	
 	/**
 	 * Get numInvalidConditions.
 	 * @return numInvalidConditions as int
 	 */
 	int getNumInvalidConditions() const;
	
 
 	
 	
 	/**
 	 * Set numInvalidConditions with the specified int.
 	 * @param numInvalidConditions The int value to which numInvalidConditions is to be set.
 	 
 		
 			
 	 */
 	void setNumInvalidConditions (int numInvalidConditions);
  		
	
	
	


	
	// ===> Attribute invalidConditions
	
	
	

	
 	/**
 	 * Get invalidConditions.
 	 * @return invalidConditions as vector<InvalidatingConditionMod::InvalidatingCondition >
 	 */
 	vector<InvalidatingConditionMod::InvalidatingCondition > getInvalidConditions() const;
	
 
 	
 	
 	/**
 	 * Set invalidConditions with the specified vector<InvalidatingConditionMod::InvalidatingCondition >.
 	 * @param invalidConditions The vector<InvalidatingConditionMod::InvalidatingCondition > value to which invalidConditions is to be set.
 	 
 		
 			
 	 */
 	void setInvalidConditions (vector<InvalidatingConditionMod::InvalidatingCondition > invalidConditions);
  		
	
	
	


	
	// ===> Attribute timeReduced
	
	
	

	
 	/**
 	 * Get timeReduced.
 	 * @return timeReduced as ArrayTime
 	 */
 	ArrayTime getTimeReduced() const;
	
 
 	
 	
 	/**
 	 * Set timeReduced with the specified ArrayTime.
 	 * @param timeReduced The ArrayTime value to which timeReduced is to be set.
 	 
 		
 			
 	 */
 	void setTimeReduced (ArrayTime timeReduced);
  		
	
	
	


	
	// ===> Attribute messages
	
	
	

	
 	/**
 	 * Get messages.
 	 * @return messages as string
 	 */
 	string getMessages() const;
	
 
 	
 	
 	/**
 	 * Set messages with the specified string.
 	 * @param messages The string value to which messages is to be set.
 	 
 		
 			
 	 */
 	void setMessages (string messages);
  		
	
	
	


	
	// ===> Attribute software
	
	
	

	
 	/**
 	 * Get software.
 	 * @return software as string
 	 */
 	string getSoftware() const;
	
 
 	
 	
 	/**
 	 * Set software with the specified string.
 	 * @param software The string value to which software is to be set.
 	 
 		
 			
 	 */
 	void setSoftware (string software);
  		
	
	
	


	
	// ===> Attribute softwareVersion
	
	
	

	
 	/**
 	 * Get softwareVersion.
 	 * @return softwareVersion as string
 	 */
 	string getSoftwareVersion() const;
	
 
 	
 	
 	/**
 	 * Set softwareVersion with the specified string.
 	 * @param softwareVersion The string value to which softwareVersion is to be set.
 	 
 		
 			
 	 */
 	void setSoftwareVersion (string softwareVersion);
  		
	
	
	


	////////////////////////////////
	// Extrinsic Table Attributes //
	////////////////////////////////
	
	///////////
	// Links //
	///////////
	
	
	
	
	/**
	 * Compare each mandatory attribute except the autoincrementable one of this CalReductionRow with 
	 * the corresponding parameters and return true if there is a match and false otherwise.
	 	
	 * @param numApplied
	    
	 * @param appliedCalibrations
	    
	 * @param numParam
	    
	 * @param paramSet
	    
	 * @param numInvalidConditions
	    
	 * @param invalidConditions
	    
	 * @param timeReduced
	    
	 * @param messages
	    
	 * @param software
	    
	 * @param softwareVersion
	    
	 */ 
	bool compareNoAutoInc(int numApplied, vector<string > appliedCalibrations, int numParam, vector<string > paramSet, int numInvalidConditions, vector<InvalidatingConditionMod::InvalidatingCondition > invalidConditions, ArrayTime timeReduced, string messages, string software, string softwareVersion);
	
	

	
	/**
	 * Compare each mandatory value (i.e. not in the key) attribute  with 
	 * the corresponding parameters and return true if there is a match and false otherwise.
	 	
	 * @param numApplied
	    
	 * @param appliedCalibrations
	    
	 * @param numParam
	    
	 * @param paramSet
	    
	 * @param numInvalidConditions
	    
	 * @param invalidConditions
	    
	 * @param timeReduced
	    
	 * @param messages
	    
	 * @param software
	    
	 * @param softwareVersion
	    
	 */ 
	bool compareRequiredValue(int numApplied, vector<string > appliedCalibrations, int numParam, vector<string > paramSet, int numInvalidConditions, vector<InvalidatingConditionMod::InvalidatingCondition > invalidConditions, ArrayTime timeReduced, string messages, string software, string softwareVersion); 
		 
	
	/**
	 * Return true if all required attributes of the value part are equal to their homologues
	 * in x and false otherwise.
	 *
	 * @param x a pointer on the CalReductionRow whose required attributes of the value part 
	 * will be compared with those of this.
	 * @return a boolean.
	 */
	bool equalByRequiredValue(CalReductionRow* x) ;
	
#ifndef WITHOUT_ACS
	/**
	 * Return this row in the form of an IDL struct.
	 * @return The values of this row as a CalReductionRowIDL struct.
	 */
	asdmIDL::CalReductionRowIDL *toIDL() const;
	
	/**
	 * Define the content of a CalReductionRowIDL struct from the values
	 * found in this row.
	 *
	 * @param x a reference to the CalReductionRowIDL struct to be set.
	 *
	 */
	 void toIDL(asdmIDL::CalReductionRowIDL& x) const;
#endif
	
#ifndef WITHOUT_ACS
	/**
	 * Fill the values of this row from the IDL struct CalReductionRowIDL.
	 * @param x The IDL struct containing the values used to fill this row.
	 * @throws ConversionException
	 */
	void setFromIDL (asdmIDL::CalReductionRowIDL x) ;
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

	std::map<std::string, CalReductionAttributeFromBin> fromBinMethods;
void calReductionIdFromBin( EndianIStream& eis);
void numAppliedFromBin( EndianIStream& eis);
void appliedCalibrationsFromBin( EndianIStream& eis);
void numParamFromBin( EndianIStream& eis);
void paramSetFromBin( EndianIStream& eis);
void numInvalidConditionsFromBin( EndianIStream& eis);
void invalidConditionsFromBin( EndianIStream& eis);
void timeReducedFromBin( EndianIStream& eis);
void messagesFromBin( EndianIStream& eis);
void softwareFromBin( EndianIStream& eis);
void softwareVersionFromBin( EndianIStream& eis);

	

	 /**
	  * Deserialize a stream of bytes read from an EndianIStream to build a PointingRow.
	  * @param eiss the EndianIStream to be read.
	  * @param table the CalReductionTable to which the row built by deserialization will be parented.
	  * @param attributesSeq a vector containing the names of the attributes . The elements order defines the order 
	  * in which the attributes are written in the binary serialization.
	  */
	 static CalReductionRow* fromBin(EndianIStream& eis, CalReductionTable& table, const std::vector<std::string>& attributesSeq);	 
 
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
	CalReductionTable &table;
	/**
	 * Whether this row has been added to the table or not.
	 */
	bool hasBeenAdded;

	// This method is used by the Table class when this row is added to the table.
	void isAdded(bool added);


	/**
	 * Create a CalReductionRow.
	 * <p>
	 * This constructor is private because only the
	 * table can create rows.  All rows know the table
	 * to which they belong.
	 * @param table The table to which this row belongs.
	 */ 
	CalReductionRow (CalReductionTable &table);

	/**
	 * Create a CalReductionRow using a copy constructor mechanism.
	 * <p>
	 * Given a CalReductionRow row and a CalReductionTable table, the method creates a new
	 * CalReductionRow owned by table. Each attribute of the created row is a copy (deep)
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
	 CalReductionRow (CalReductionTable &table, CalReductionRow &row);
	 	
	////////////////////////////////
	// Intrinsic Table Attributes //
	////////////////////////////////
	
	
	// ===> Attribute calReductionId
	
	

	Tag calReductionId;

	
	
 	
 	/**
 	 * Set calReductionId with the specified Tag value.
 	 * @param calReductionId The Tag value to which calReductionId is to be set.
		
 		
			
 	 * @throw IllegalAccessException If an attempt is made to change this field after is has been added to the table.
 	 		
 	 */
 	void setCalReductionId (Tag calReductionId);
  		
	

	
	// ===> Attribute numApplied
	
	

	int numApplied;

	
	
 	

	
	// ===> Attribute appliedCalibrations
	
	

	vector<string > appliedCalibrations;

	
	
 	

	
	// ===> Attribute numParam
	
	

	int numParam;

	
	
 	

	
	// ===> Attribute paramSet
	
	

	vector<string > paramSet;

	
	
 	

	
	// ===> Attribute numInvalidConditions
	
	

	int numInvalidConditions;

	
	
 	

	
	// ===> Attribute invalidConditions
	
	

	vector<InvalidatingConditionMod::InvalidatingCondition > invalidConditions;

	
	
 	

	
	// ===> Attribute timeReduced
	
	

	ArrayTime timeReduced;

	
	
 	

	
	// ===> Attribute messages
	
	

	string messages;

	
	
 	

	
	// ===> Attribute software
	
	

	string software;

	
	
 	

	
	// ===> Attribute softwareVersion
	
	

	string softwareVersion;

	
	
 	

	////////////////////////////////
	// Extrinsic Table Attributes //
	////////////////////////////////
	
	///////////
	// Links //
	///////////
	
	
/*
	////////////////////////////////////////////////////////////
	// binary-deserialization material from an EndianIStream  //
	////////////////////////////////////////////////////////////
	std::map<std::string, CalReductionAttributeFromBin> fromBinMethods;
void calReductionIdFromBin( EndianIStream& eis);
void numAppliedFromBin( EndianIStream& eis);
void appliedCalibrationsFromBin( EndianIStream& eis);
void numParamFromBin( EndianIStream& eis);
void paramSetFromBin( EndianIStream& eis);
void numInvalidConditionsFromBin( EndianIStream& eis);
void invalidConditionsFromBin( EndianIStream& eis);
void timeReducedFromBin( EndianIStream& eis);
void messagesFromBin( EndianIStream& eis);
void softwareFromBin( EndianIStream& eis);
void softwareVersionFromBin( EndianIStream& eis);

	
*/
	
	///////////////////////////////////
	// text-deserialization material //
	///////////////////////////////////
	std::map<std::string, CalReductionAttributeFromText> fromTextMethods;
	
void calReductionIdFromText (const string & s);
	
	
void numAppliedFromText (const string & s);
	
	
void appliedCalibrationsFromText (const string & s);
	
	
void numParamFromText (const string & s);
	
	
void paramSetFromText (const string & s);
	
	
void numInvalidConditionsFromText (const string & s);
	
	
void invalidConditionsFromText (const string & s);
	
	
void timeReducedFromText (const string & s);
	
	
void messagesFromText (const string & s);
	
	
void softwareFromText (const string & s);
	
	
void softwareVersionFromText (const string & s);
	

		
	
	/**
	 * Serialize this into a stream of bytes written to an EndianOSStream.
	 * @param eoss the EndianOSStream to be written to
	 */
	 void toBin(EndianOSStream& eoss);
	 	 
	 /**
	  * Deserialize a stream of bytes read from an EndianIStream to build a PointingRow.
	  * @param eiss the EndianIStream to be read.
	  * @param table the CalReductionTable to which the row built by deserialization will be parented.
	  * @param attributesSeq a vector containing the names of the attributes . The elements order defines the order 
	  * in which the attributes are written in the binary serialization.

	 static CalReductionRow* fromBin(EndianIStream& eis, CalReductionTable& table, const std::vector<std::string>& attributesSeq);	 
		*/
};

} // End namespace asdm

#endif /* CalReduction_CLASS */
