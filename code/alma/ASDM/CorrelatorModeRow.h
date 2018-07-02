
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
 * File CorrelatorModeRow.h
 */
 
#ifndef CorrelatorModeRow_CLASS
#define CorrelatorModeRow_CLASS

#include <vector>
#include <string>
#include <set>

#ifndef WITHOUT_ACS
#include <asdmIDLC.h>
#endif






	 
#include <Tag.h>
	




	

	

	
#include "CBasebandName.h"
	

	

	
#include "CAccumMode.h"
	

	

	

	
#include "CAxisName.h"
	

	
#include "CFilterMode.h"
	

	
#include "CCorrelatorName.h"
	



#include <ConversionException.h>
#include <NoSuchRow.h>
#include <IllegalAccessException.h>

#include <RowTransformer.h>
//#include <TableStreamReader.h>

/*\file CorrelatorMode.h
    \brief Generated from model's revision "-1", branch ""
*/

namespace asdm {

//class asdm::CorrelatorModeTable;

	

class CorrelatorModeRow;
typedef void (CorrelatorModeRow::*CorrelatorModeAttributeFromBin) (EndianIStream& eis);
typedef void (CorrelatorModeRow::*CorrelatorModeAttributeFromText) (const string& s);

/**
 * The CorrelatorModeRow class is a row of a CorrelatorModeTable.
 * 
 * Generated from model's revision "-1", branch ""
 *
 */
class CorrelatorModeRow {
friend class asdm::CorrelatorModeTable;
friend class asdm::RowTransformer<CorrelatorModeRow>;
//friend class asdm::TableStreamReader<CorrelatorModeTable, CorrelatorModeRow>;

public:

	virtual ~CorrelatorModeRow();

	/**
	 * Return the table to which this row belongs.
	 */
	CorrelatorModeTable &getTable() const;
	
	/**
	 * Has this row been added to its table ?
	 * @return true if and only if it has been added.
	 */
	bool isAdded() const;
		
	////////////////////////////////
	// Intrinsic Table Attributes //
	////////////////////////////////
	
	
	// ===> Attribute correlatorModeId
	
	
	

	
 	/**
 	 * Get correlatorModeId.
 	 * @return correlatorModeId as Tag
 	 */
 	Tag getCorrelatorModeId() const;
	
 
 	
 	
	
	


	
	// ===> Attribute numBaseband
	
	
	

	
 	/**
 	 * Get numBaseband.
 	 * @return numBaseband as int
 	 */
 	int getNumBaseband() const;
	
 
 	
 	
 	/**
 	 * Set numBaseband with the specified int.
 	 * @param numBaseband The int value to which numBaseband is to be set.
 	 
 		
 			
 	 */
 	void setNumBaseband (int numBaseband);
  		
	
	
	


	
	// ===> Attribute basebandNames
	
	
	

	
 	/**
 	 * Get basebandNames.
 	 * @return basebandNames as vector<BasebandNameMod::BasebandName >
 	 */
 	vector<BasebandNameMod::BasebandName > getBasebandNames() const;
	
 
 	
 	
 	/**
 	 * Set basebandNames with the specified vector<BasebandNameMod::BasebandName >.
 	 * @param basebandNames The vector<BasebandNameMod::BasebandName > value to which basebandNames is to be set.
 	 
 		
 			
 	 */
 	void setBasebandNames (vector<BasebandNameMod::BasebandName > basebandNames);
  		
	
	
	


	
	// ===> Attribute basebandConfig
	
	
	

	
 	/**
 	 * Get basebandConfig.
 	 * @return basebandConfig as vector<int >
 	 */
 	vector<int > getBasebandConfig() const;
	
 
 	
 	
 	/**
 	 * Set basebandConfig with the specified vector<int >.
 	 * @param basebandConfig The vector<int > value to which basebandConfig is to be set.
 	 
 		
 			
 	 */
 	void setBasebandConfig (vector<int > basebandConfig);
  		
	
	
	


	
	// ===> Attribute accumMode
	
	
	

	
 	/**
 	 * Get accumMode.
 	 * @return accumMode as AccumModeMod::AccumMode
 	 */
 	AccumModeMod::AccumMode getAccumMode() const;
	
 
 	
 	
 	/**
 	 * Set accumMode with the specified AccumModeMod::AccumMode.
 	 * @param accumMode The AccumModeMod::AccumMode value to which accumMode is to be set.
 	 
 		
 			
 	 */
 	void setAccumMode (AccumModeMod::AccumMode accumMode);
  		
	
	
	


	
	// ===> Attribute binMode
	
	
	

	
 	/**
 	 * Get binMode.
 	 * @return binMode as int
 	 */
 	int getBinMode() const;
	
 
 	
 	
 	/**
 	 * Set binMode with the specified int.
 	 * @param binMode The int value to which binMode is to be set.
 	 
 		
 			
 	 */
 	void setBinMode (int binMode);
  		
	
	
	


	
	// ===> Attribute numAxes
	
	
	

	
 	/**
 	 * Get numAxes.
 	 * @return numAxes as int
 	 */
 	int getNumAxes() const;
	
 
 	
 	
 	/**
 	 * Set numAxes with the specified int.
 	 * @param numAxes The int value to which numAxes is to be set.
 	 
 		
 			
 	 */
 	void setNumAxes (int numAxes);
  		
	
	
	


	
	// ===> Attribute axesOrderArray
	
	
	

	
 	/**
 	 * Get axesOrderArray.
 	 * @return axesOrderArray as vector<AxisNameMod::AxisName >
 	 */
 	vector<AxisNameMod::AxisName > getAxesOrderArray() const;
	
 
 	
 	
 	/**
 	 * Set axesOrderArray with the specified vector<AxisNameMod::AxisName >.
 	 * @param axesOrderArray The vector<AxisNameMod::AxisName > value to which axesOrderArray is to be set.
 	 
 		
 			
 	 */
 	void setAxesOrderArray (vector<AxisNameMod::AxisName > axesOrderArray);
  		
	
	
	


	
	// ===> Attribute filterMode
	
	
	

	
 	/**
 	 * Get filterMode.
 	 * @return filterMode as vector<FilterModeMod::FilterMode >
 	 */
 	vector<FilterModeMod::FilterMode > getFilterMode() const;
	
 
 	
 	
 	/**
 	 * Set filterMode with the specified vector<FilterModeMod::FilterMode >.
 	 * @param filterMode The vector<FilterModeMod::FilterMode > value to which filterMode is to be set.
 	 
 		
 			
 	 */
 	void setFilterMode (vector<FilterModeMod::FilterMode > filterMode);
  		
	
	
	


	
	// ===> Attribute correlatorName
	
	
	

	
 	/**
 	 * Get correlatorName.
 	 * @return correlatorName as CorrelatorNameMod::CorrelatorName
 	 */
 	CorrelatorNameMod::CorrelatorName getCorrelatorName() const;
	
 
 	
 	
 	/**
 	 * Set correlatorName with the specified CorrelatorNameMod::CorrelatorName.
 	 * @param correlatorName The CorrelatorNameMod::CorrelatorName value to which correlatorName is to be set.
 	 
 		
 			
 	 */
 	void setCorrelatorName (CorrelatorNameMod::CorrelatorName correlatorName);
  		
	
	
	


	////////////////////////////////
	// Extrinsic Table Attributes //
	////////////////////////////////
	
	///////////
	// Links //
	///////////
	
	
	
	
	/**
	 * Compare each mandatory attribute except the autoincrementable one of this CorrelatorModeRow with 
	 * the corresponding parameters and return true if there is a match and false otherwise.
	 	
	 * @param numBaseband
	    
	 * @param basebandNames
	    
	 * @param basebandConfig
	    
	 * @param accumMode
	    
	 * @param binMode
	    
	 * @param numAxes
	    
	 * @param axesOrderArray
	    
	 * @param filterMode
	    
	 * @param correlatorName
	    
	 */ 
	bool compareNoAutoInc(int numBaseband, vector<BasebandNameMod::BasebandName > basebandNames, vector<int > basebandConfig, AccumModeMod::AccumMode accumMode, int binMode, int numAxes, vector<AxisNameMod::AxisName > axesOrderArray, vector<FilterModeMod::FilterMode > filterMode, CorrelatorNameMod::CorrelatorName correlatorName);
	
	

	
	/**
	 * Compare each mandatory value (i.e. not in the key) attribute  with 
	 * the corresponding parameters and return true if there is a match and false otherwise.
	 	
	 * @param numBaseband
	    
	 * @param basebandNames
	    
	 * @param basebandConfig
	    
	 * @param accumMode
	    
	 * @param binMode
	    
	 * @param numAxes
	    
	 * @param axesOrderArray
	    
	 * @param filterMode
	    
	 * @param correlatorName
	    
	 */ 
	bool compareRequiredValue(int numBaseband, vector<BasebandNameMod::BasebandName > basebandNames, vector<int > basebandConfig, AccumModeMod::AccumMode accumMode, int binMode, int numAxes, vector<AxisNameMod::AxisName > axesOrderArray, vector<FilterModeMod::FilterMode > filterMode, CorrelatorNameMod::CorrelatorName correlatorName); 
		 
	
	/**
	 * Return true if all required attributes of the value part are equal to their homologues
	 * in x and false otherwise.
	 *
	 * @param x a pointer on the CorrelatorModeRow whose required attributes of the value part 
	 * will be compared with those of this.
	 * @return a boolean.
	 */
	bool equalByRequiredValue(CorrelatorModeRow* x) ;
	
#ifndef WITHOUT_ACS
	/**
	 * Return this row in the form of an IDL struct.
	 * @return The values of this row as a CorrelatorModeRowIDL struct.
	 */
	asdmIDL::CorrelatorModeRowIDL *toIDL() const;
	
	/**
	 * Define the content of a CorrelatorModeRowIDL struct from the values
	 * found in this row.
	 *
	 * @param x a reference to the CorrelatorModeRowIDL struct to be set.
	 *
	 */
	 void toIDL(asdmIDL::CorrelatorModeRowIDL& x) const;
#endif
	
#ifndef WITHOUT_ACS
	/**
	 * Fill the values of this row from the IDL struct CorrelatorModeRowIDL.
	 * @param x The IDL struct containing the values used to fill this row.
	 * @throws ConversionException
	 */
	void setFromIDL (asdmIDL::CorrelatorModeRowIDL x) ;
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

	std::map<std::string, CorrelatorModeAttributeFromBin> fromBinMethods;
void correlatorModeIdFromBin( EndianIStream& eis);
void numBasebandFromBin( EndianIStream& eis);
void basebandNamesFromBin( EndianIStream& eis);
void basebandConfigFromBin( EndianIStream& eis);
void accumModeFromBin( EndianIStream& eis);
void binModeFromBin( EndianIStream& eis);
void numAxesFromBin( EndianIStream& eis);
void axesOrderArrayFromBin( EndianIStream& eis);
void filterModeFromBin( EndianIStream& eis);
void correlatorNameFromBin( EndianIStream& eis);

	

	 /**
	  * Deserialize a stream of bytes read from an EndianIStream to build a PointingRow.
	  * @param eiss the EndianIStream to be read.
	  * @param table the CorrelatorModeTable to which the row built by deserialization will be parented.
	  * @param attributesSeq a vector containing the names of the attributes . The elements order defines the order 
	  * in which the attributes are written in the binary serialization.
	  */
	 static CorrelatorModeRow* fromBin(EndianIStream& eis, CorrelatorModeTable& table, const std::vector<std::string>& attributesSeq);	 
 
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
	CorrelatorModeTable &table;
	/**
	 * Whether this row has been added to the table or not.
	 */
	bool hasBeenAdded;

	// This method is used by the Table class when this row is added to the table.
	void isAdded(bool added);


	/**
	 * Create a CorrelatorModeRow.
	 * <p>
	 * This constructor is private because only the
	 * table can create rows.  All rows know the table
	 * to which they belong.
	 * @param table The table to which this row belongs.
	 */ 
	CorrelatorModeRow (CorrelatorModeTable &table);

	/**
	 * Create a CorrelatorModeRow using a copy constructor mechanism.
	 * <p>
	 * Given a CorrelatorModeRow row and a CorrelatorModeTable table, the method creates a new
	 * CorrelatorModeRow owned by table. Each attribute of the created row is a copy (deep)
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
	 CorrelatorModeRow (CorrelatorModeTable &table, CorrelatorModeRow &row);
	 	
	////////////////////////////////
	// Intrinsic Table Attributes //
	////////////////////////////////
	
	
	// ===> Attribute correlatorModeId
	
	

	Tag correlatorModeId;

	
	
 	
 	/**
 	 * Set correlatorModeId with the specified Tag value.
 	 * @param correlatorModeId The Tag value to which correlatorModeId is to be set.
		
 		
			
 	 * @throw IllegalAccessException If an attempt is made to change this field after is has been added to the table.
 	 		
 	 */
 	void setCorrelatorModeId (Tag correlatorModeId);
  		
	

	
	// ===> Attribute numBaseband
	
	

	int numBaseband;

	
	
 	

	
	// ===> Attribute basebandNames
	
	

	vector<BasebandNameMod::BasebandName > basebandNames;

	
	
 	

	
	// ===> Attribute basebandConfig
	
	

	vector<int > basebandConfig;

	
	
 	

	
	// ===> Attribute accumMode
	
	

	AccumModeMod::AccumMode accumMode;

	
	
 	

	
	// ===> Attribute binMode
	
	

	int binMode;

	
	
 	

	
	// ===> Attribute numAxes
	
	

	int numAxes;

	
	
 	

	
	// ===> Attribute axesOrderArray
	
	

	vector<AxisNameMod::AxisName > axesOrderArray;

	
	
 	

	
	// ===> Attribute filterMode
	
	

	vector<FilterModeMod::FilterMode > filterMode;

	
	
 	

	
	// ===> Attribute correlatorName
	
	

	CorrelatorNameMod::CorrelatorName correlatorName;

	
	
 	

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
	std::map<std::string, CorrelatorModeAttributeFromBin> fromBinMethods;
void correlatorModeIdFromBin( EndianIStream& eis);
void numBasebandFromBin( EndianIStream& eis);
void basebandNamesFromBin( EndianIStream& eis);
void basebandConfigFromBin( EndianIStream& eis);
void accumModeFromBin( EndianIStream& eis);
void binModeFromBin( EndianIStream& eis);
void numAxesFromBin( EndianIStream& eis);
void axesOrderArrayFromBin( EndianIStream& eis);
void filterModeFromBin( EndianIStream& eis);
void correlatorNameFromBin( EndianIStream& eis);

	
*/
	
	///////////////////////////////////
	// text-deserialization material //
	///////////////////////////////////
	std::map<std::string, CorrelatorModeAttributeFromText> fromTextMethods;
	
void correlatorModeIdFromText (const string & s);
	
	
void numBasebandFromText (const string & s);
	
	
void basebandNamesFromText (const string & s);
	
	
void basebandConfigFromText (const string & s);
	
	
void accumModeFromText (const string & s);
	
	
void binModeFromText (const string & s);
	
	
void numAxesFromText (const string & s);
	
	
void axesOrderArrayFromText (const string & s);
	
	
void filterModeFromText (const string & s);
	
	
void correlatorNameFromText (const string & s);
	

		
	
	/**
	 * Serialize this into a stream of bytes written to an EndianOSStream.
	 * @param eoss the EndianOSStream to be written to
	 */
	 void toBin(EndianOSStream& eoss);
	 	 
	 /**
	  * Deserialize a stream of bytes read from an EndianIStream to build a PointingRow.
	  * @param eiss the EndianIStream to be read.
	  * @param table the CorrelatorModeTable to which the row built by deserialization will be parented.
	  * @param attributesSeq a vector containing the names of the attributes . The elements order defines the order 
	  * in which the attributes are written in the binary serialization.

	 static CorrelatorModeRow* fromBin(EndianIStream& eis, CorrelatorModeTable& table, const std::vector<std::string>& attributesSeq);	 
		*/
};

} // End namespace asdm

#endif /* CorrelatorMode_CLASS */
