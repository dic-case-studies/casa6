
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
 * File AnnotationRow.h
 */
 
#ifndef AnnotationRow_CLASS
#define AnnotationRow_CLASS

#include <vector>
#include <string>
#include <set>

#ifndef WITHOUT_ACS
#include <asdmIDLC.h>
#endif



#include <stdint.h>




	 
#include <ArrayTime.h>
	

	 
#include <Tag.h>
	

	 
#include <Interval.h>
	




	

	

	

	

	

	
#include "CBasebandName.h"
	

	

	

	

	

	

	

	

	



#include <ConversionException.h>
#include <NoSuchRow.h>
#include <IllegalAccessException.h>

#include <RowTransformer.h>
//#include <TableStreamReader.h>

/*\file Annotation.h
    \brief Generated from model's revision "-1", branch ""
*/

namespace asdm {

//class asdm::AnnotationTable;


// class asdm::AntennaRow;
class AntennaRow;
	

class AnnotationRow;
typedef void (AnnotationRow::*AnnotationAttributeFromBin) (EndianIStream& eis);
typedef void (AnnotationRow::*AnnotationAttributeFromText) (const string& s);

/**
 * The AnnotationRow class is a row of a AnnotationTable.
 * 
 * Generated from model's revision "-1", branch ""
 *
 */
class AnnotationRow {
friend class asdm::AnnotationTable;
friend class asdm::RowTransformer<AnnotationRow>;
//friend class asdm::TableStreamReader<AnnotationTable, AnnotationRow>;

public:

	virtual ~AnnotationRow();

	/**
	 * Return the table to which this row belongs.
	 */
	AnnotationTable &getTable() const;
	
	/**
	 * Has this row been added to its table ?
	 * @return true if and only if it has been added.
	 */
	bool isAdded() const;
		
	////////////////////////////////
	// Intrinsic Table Attributes //
	////////////////////////////////
	
	
	// ===> Attribute annotationId
	
	
	

	
 	/**
 	 * Get annotationId.
 	 * @return annotationId as Tag
 	 */
 	Tag getAnnotationId() const;
	
 
 	
 	
	
	


	
	// ===> Attribute time
	
	
	

	
 	/**
 	 * Get time.
 	 * @return time as ArrayTime
 	 */
 	ArrayTime getTime() const;
	
 
 	
 	
 	/**
 	 * Set time with the specified ArrayTime.
 	 * @param time The ArrayTime value to which time is to be set.
 	 
 		
 			
 	 */
 	void setTime (ArrayTime time);
  		
	
	
	


	
	// ===> Attribute issue
	
	
	

	
 	/**
 	 * Get issue.
 	 * @return issue as string
 	 */
 	string getIssue() const;
	
 
 	
 	
 	/**
 	 * Set issue with the specified string.
 	 * @param issue The string value to which issue is to be set.
 	 
 		
 			
 	 */
 	void setIssue (string issue);
  		
	
	
	


	
	// ===> Attribute details
	
	
	

	
 	/**
 	 * Get details.
 	 * @return details as string
 	 */
 	string getDetails() const;
	
 
 	
 	
 	/**
 	 * Set details with the specified string.
 	 * @param details The string value to which details is to be set.
 	 
 		
 			
 	 */
 	void setDetails (string details);
  		
	
	
	


	
	// ===> Attribute numAntenna, which is optional
	
	
	
	/**
	 * The attribute numAntenna is optional. Return true if this attribute exists.
	 * @return true if and only if the numAntenna attribute exists. 
	 */
	bool isNumAntennaExists() const;
	

	
 	/**
 	 * Get numAntenna, which is optional.
 	 * @return numAntenna as int
 	 * @throws IllegalAccessException If numAntenna does not exist.
 	 */
 	int getNumAntenna() const;
	
 
 	
 	
 	/**
 	 * Set numAntenna with the specified int.
 	 * @param numAntenna The int value to which numAntenna is to be set.
 	 
 		
 	 */
 	void setNumAntenna (int numAntenna);
		
	
	
	
	/**
	 * Mark numAntenna, which is an optional field, as non-existent.
	 */
	void clearNumAntenna ();
	


	
	// ===> Attribute basebandName, which is optional
	
	
	
	/**
	 * The attribute basebandName is optional. Return true if this attribute exists.
	 * @return true if and only if the basebandName attribute exists. 
	 */
	bool isBasebandNameExists() const;
	

	
 	/**
 	 * Get basebandName, which is optional.
 	 * @return basebandName as vector<BasebandNameMod::BasebandName >
 	 * @throws IllegalAccessException If basebandName does not exist.
 	 */
 	std::vector<BasebandNameMod::BasebandName > getBasebandName() const;
	
 
 	
 	
 	/**
 	 * Set basebandName with the specified vector<BasebandNameMod::BasebandName >.
 	 * @param basebandName The vector<BasebandNameMod::BasebandName > value to which basebandName is to be set.
 	 
 		
 	 */
 	void setBasebandName (std::vector<BasebandNameMod::BasebandName > basebandName);
		
	
	
	
	/**
	 * Mark basebandName, which is an optional field, as non-existent.
	 */
	void clearBasebandName ();
	


	
	// ===> Attribute numBaseband, which is optional
	
	
	
	/**
	 * The attribute numBaseband is optional. Return true if this attribute exists.
	 * @return true if and only if the numBaseband attribute exists. 
	 */
	bool isNumBasebandExists() const;
	

	
 	/**
 	 * Get numBaseband, which is optional.
 	 * @return numBaseband as int
 	 * @throws IllegalAccessException If numBaseband does not exist.
 	 */
 	int getNumBaseband() const;
	
 
 	
 	
 	/**
 	 * Set numBaseband with the specified int.
 	 * @param numBaseband The int value to which numBaseband is to be set.
 	 
 		
 	 */
 	void setNumBaseband (int numBaseband);
		
	
	
	
	/**
	 * Mark numBaseband, which is an optional field, as non-existent.
	 */
	void clearNumBaseband ();
	


	
	// ===> Attribute interval, which is optional
	
	
	
	/**
	 * The attribute interval is optional. Return true if this attribute exists.
	 * @return true if and only if the interval attribute exists. 
	 */
	bool isIntervalExists() const;
	

	
 	/**
 	 * Get interval, which is optional.
 	 * @return interval as Interval
 	 * @throws IllegalAccessException If interval does not exist.
 	 */
 	Interval getInterval() const;
	
 
 	
 	
 	/**
 	 * Set interval with the specified Interval.
 	 * @param interval The Interval value to which interval is to be set.
 	 
 		
 	 */
 	void setInterval (Interval interval);
		
	
	
	
	/**
	 * Mark interval, which is an optional field, as non-existent.
	 */
	void clearInterval ();
	


	
	// ===> Attribute dValue, which is optional
	
	
	
	/**
	 * The attribute dValue is optional. Return true if this attribute exists.
	 * @return true if and only if the dValue attribute exists. 
	 */
	bool isDValueExists() const;
	

	
 	/**
 	 * Get dValue, which is optional.
 	 * @return dValue as double
 	 * @throws IllegalAccessException If dValue does not exist.
 	 */
 	double getDValue() const;
	
 
 	
 	
 	/**
 	 * Set dValue with the specified double.
 	 * @param dValue The double value to which dValue is to be set.
 	 
 		
 	 */
 	void setDValue (double dValue);
		
	
	
	
	/**
	 * Mark dValue, which is an optional field, as non-existent.
	 */
	void clearDValue ();
	


	
	// ===> Attribute vdValue, which is optional
	
	
	
	/**
	 * The attribute vdValue is optional. Return true if this attribute exists.
	 * @return true if and only if the vdValue attribute exists. 
	 */
	bool isVdValueExists() const;
	

	
 	/**
 	 * Get vdValue, which is optional.
 	 * @return vdValue as vector<double >
 	 * @throws IllegalAccessException If vdValue does not exist.
 	 */
 	std::vector<double > getVdValue() const;
	
 
 	
 	
 	/**
 	 * Set vdValue with the specified vector<double >.
 	 * @param vdValue The vector<double > value to which vdValue is to be set.
 	 
 		
 	 */
 	void setVdValue (std::vector<double > vdValue);
		
	
	
	
	/**
	 * Mark vdValue, which is an optional field, as non-existent.
	 */
	void clearVdValue ();
	


	
	// ===> Attribute vvdValues, which is optional
	
	
	
	/**
	 * The attribute vvdValues is optional. Return true if this attribute exists.
	 * @return true if and only if the vvdValues attribute exists. 
	 */
	bool isVvdValuesExists() const;
	

	
 	/**
 	 * Get vvdValues, which is optional.
 	 * @return vvdValues as vector<vector<double > >
 	 * @throws IllegalAccessException If vvdValues does not exist.
 	 */
 	std::vector<std::vector<double > > getVvdValues() const;
	
 
 	
 	
 	/**
 	 * Set vvdValues with the specified vector<vector<double > >.
 	 * @param vvdValues The vector<vector<double > > value to which vvdValues is to be set.
 	 
 		
 	 */
 	void setVvdValues (std::vector<std::vector<double > > vvdValues);
		
	
	
	
	/**
	 * Mark vvdValues, which is an optional field, as non-existent.
	 */
	void clearVvdValues ();
	


	
	// ===> Attribute llValue, which is optional
	
	
	
	/**
	 * The attribute llValue is optional. Return true if this attribute exists.
	 * @return true if and only if the llValue attribute exists. 
	 */
	bool isLlValueExists() const;
	

	
 	/**
 	 * Get llValue, which is optional.
 	 * @return llValue as int64_t
 	 * @throws IllegalAccessException If llValue does not exist.
 	 */
 	int64_t getLlValue() const;
	
 
 	
 	
 	/**
 	 * Set llValue with the specified int64_t.
 	 * @param llValue The int64_t value to which llValue is to be set.
 	 
 		
 	 */
 	void setLlValue (int64_t llValue);
		
	
	
	
	/**
	 * Mark llValue, which is an optional field, as non-existent.
	 */
	void clearLlValue ();
	


	
	// ===> Attribute vllValue, which is optional
	
	
	
	/**
	 * The attribute vllValue is optional. Return true if this attribute exists.
	 * @return true if and only if the vllValue attribute exists. 
	 */
	bool isVllValueExists() const;
	

	
 	/**
 	 * Get vllValue, which is optional.
 	 * @return vllValue as vector<int64_t >
 	 * @throws IllegalAccessException If vllValue does not exist.
 	 */
 	std::vector<int64_t > getVllValue() const;
	
 
 	
 	
 	/**
 	 * Set vllValue with the specified vector<int64_t >.
 	 * @param vllValue The vector<int64_t > value to which vllValue is to be set.
 	 
 		
 	 */
 	void setVllValue (std::vector<int64_t > vllValue);
		
	
	
	
	/**
	 * Mark vllValue, which is an optional field, as non-existent.
	 */
	void clearVllValue ();
	


	
	// ===> Attribute vvllValue, which is optional
	
	
	
	/**
	 * The attribute vvllValue is optional. Return true if this attribute exists.
	 * @return true if and only if the vvllValue attribute exists. 
	 */
	bool isVvllValueExists() const;
	

	
 	/**
 	 * Get vvllValue, which is optional.
 	 * @return vvllValue as vector<vector<int64_t > >
 	 * @throws IllegalAccessException If vvllValue does not exist.
 	 */
 	std::vector<std::vector<int64_t > > getVvllValue() const;
	
 
 	
 	
 	/**
 	 * Set vvllValue with the specified vector<vector<int64_t > >.
 	 * @param vvllValue The vector<vector<int64_t > > value to which vvllValue is to be set.
 	 
 		
 	 */
 	void setVvllValue (std::vector<std::vector<int64_t > > vvllValue);
		
	
	
	
	/**
	 * Mark vvllValue, which is an optional field, as non-existent.
	 */
	void clearVvllValue ();
	


	////////////////////////////////
	// Extrinsic Table Attributes //
	////////////////////////////////
	
	
	// ===> Attribute antennaId, which is optional
	
	
	
	/**
	 * The attribute antennaId is optional. Return true if this attribute exists.
	 * @return true if and only if the antennaId attribute exists. 
	 */
	bool isAntennaIdExists() const;
	

	
 	/**
 	 * Get antennaId, which is optional.
 	 * @return antennaId as vector<Tag> 
 	 * @throws IllegalAccessException If antennaId does not exist.
 	 */
 	std::vector<Tag>  getAntennaId() const;
	
 
 	
 	
 	/**
 	 * Set antennaId with the specified vector<Tag> .
 	 * @param antennaId The vector<Tag>  value to which antennaId is to be set.
 	 
 		
 	 */
 	void setAntennaId (std::vector<Tag>  antennaId);
		
	
	
	
	/**
	 * Mark antennaId, which is an optional field, as non-existent.
	 */
	void clearAntennaId ();
	


	///////////
	// Links //
	///////////
	
	
 		
 	/**
 	 * Set antennaId[i] with the specified Tag.
 	 * @param i The index in antennaId where to set the Tag value.
 	 * @param antennaId The Tag value to which antennaId[i] is to be set. 
 	 * @throws OutOfBoundsException
  	 */
  	void setAntennaId (int i, Tag antennaId)  ;
 			
	

	
		 
/**
 * Append a Tag to antennaId.
 * @param id the Tag to be appended to antennaId
 */
 void addAntennaId(Tag id); 

/**
 * Append a vector of Tag to antennaId.
 * @param id an array of Tag to be appended to antennaId
 */
 void addAntennaId(const std::vector<Tag> & id); 
 

 /**
  * Returns the Tag stored in antennaId at position i.
  * @param i the position in antennaId where the Tag is retrieved.
  * @return the Tag stored at position i in antennaId.
  */
 const Tag getAntennaId(int i);
 
 /**
  * Returns the AntennaRow linked to this row via the tag stored in antennaId
  * at position i.
  * @param i the position in antennaId.
  * @return a pointer on a AntennaRow whose key (a Tag) is equal to the Tag stored at position
  * i in the antennaId. 
  */
 AntennaRow* getAntennaUsingAntennaId(int i); 
 
 /**
  * Returns the vector of AntennaRow* linked to this row via the Tags stored in antennaId
  * @return an array of pointers on AntennaRow.
  */
 std::vector<AntennaRow *> getAntennasUsingAntennaId(); 
  

	

	
	
	
	/**
	 * Compare each mandatory attribute except the autoincrementable one of this AnnotationRow with 
	 * the corresponding parameters and return true if there is a match and false otherwise.
	 	
	 * @param time
	    
	 * @param issue
	    
	 * @param details
	    
	 */ 
	bool compareNoAutoInc(ArrayTime time, string issue, string details);
	
	

	
	/**
	 * Compare each mandatory value (i.e. not in the key) attribute  with 
	 * the corresponding parameters and return true if there is a match and false otherwise.
	 	
	 * @param time
	    
	 * @param issue
	    
	 * @param details
	    
	 */ 
	bool compareRequiredValue(ArrayTime time, string issue, string details); 
		 
	
	/**
	 * Return true if all required attributes of the value part are equal to their homologues
	 * in x and false otherwise.
	 *
	 * @param x a pointer on the AnnotationRow whose required attributes of the value part 
	 * will be compared with those of this.
	 * @return a boolean.
	 */
	bool equalByRequiredValue(AnnotationRow* x) ;
	
#ifndef WITHOUT_ACS
	/**
	 * Return this row in the form of an IDL struct.
	 * @return The values of this row as a AnnotationRowIDL struct.
	 */
	asdmIDL::AnnotationRowIDL *toIDL() const;
	
	/**
	 * Define the content of a AnnotationRowIDL struct from the values
	 * found in this row.
	 *
	 * @param x a reference to the AnnotationRowIDL struct to be set.
	 *
	 */
	 void toIDL(asdmIDL::AnnotationRowIDL& x) const;
#endif
	
#ifndef WITHOUT_ACS
	/**
	 * Fill the values of this row from the IDL struct AnnotationRowIDL.
	 * @param x The IDL struct containing the values used to fill this row.
	 * @throws ConversionException
	 */
	void setFromIDL (asdmIDL::AnnotationRowIDL x) ;
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

	std::map<std::string, AnnotationAttributeFromBin> fromBinMethods;
void annotationIdFromBin( EndianIStream& eis);
void timeFromBin( EndianIStream& eis);
void issueFromBin( EndianIStream& eis);
void detailsFromBin( EndianIStream& eis);

void numAntennaFromBin( EndianIStream& eis);
void basebandNameFromBin( EndianIStream& eis);
void numBasebandFromBin( EndianIStream& eis);
void intervalFromBin( EndianIStream& eis);
void dValueFromBin( EndianIStream& eis);
void vdValueFromBin( EndianIStream& eis);
void vvdValuesFromBin( EndianIStream& eis);
void llValueFromBin( EndianIStream& eis);
void vllValueFromBin( EndianIStream& eis);
void vvllValueFromBin( EndianIStream& eis);
void antennaIdFromBin( EndianIStream& eis);


	 /**
	  * Deserialize a stream of bytes read from an EndianIStream to build a PointingRow.
	  * @param eiss the EndianIStream to be read.
	  * @param table the AnnotationTable to which the row built by deserialization will be parented.
	  * @param attributesSeq a vector containing the names of the attributes . The elements order defines the order 
	  * in which the attributes are written in the binary serialization.
	  */
	 static AnnotationRow* fromBin(EndianIStream& eis, AnnotationTable& table, const std::vector<std::string>& attributesSeq);	 
 
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
	AnnotationTable &table;
	/**
	 * Whether this row has been added to the table or not.
	 */
	bool hasBeenAdded;

	// This method is used by the Table class when this row is added to the table.
	void isAdded(bool added);


	/**
	 * Create a AnnotationRow.
	 * <p>
	 * This constructor is private because only the
	 * table can create rows.  All rows know the table
	 * to which they belong.
	 * @param table The table to which this row belongs.
	 */ 
	AnnotationRow (AnnotationTable &table);

	/**
	 * Create a AnnotationRow using a copy constructor mechanism.
	 * <p>
	 * Given a AnnotationRow row and a AnnotationTable table, the method creates a new
	 * AnnotationRow owned by table. Each attribute of the created row is a copy (deep)
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
	 AnnotationRow (AnnotationTable &table, AnnotationRow &row);
	 	
	////////////////////////////////
	// Intrinsic Table Attributes //
	////////////////////////////////
	
	
	// ===> Attribute annotationId
	
	

	Tag annotationId;

	
	
 	
 	/**
 	 * Set annotationId with the specified Tag value.
 	 * @param annotationId The Tag value to which annotationId is to be set.
		
 		
			
 	 * @throw IllegalAccessException If an attempt is made to change this field after is has been added to the table.
 	 		
 	 */
 	void setAnnotationId (Tag annotationId);
  		
	

	
	// ===> Attribute time
	
	

	ArrayTime time;

	
	
 	

	
	// ===> Attribute issue
	
	

	string issue;

	
	
 	

	
	// ===> Attribute details
	
	

	string details;

	
	
 	

	
	// ===> Attribute numAntenna, which is optional
	
	
	bool numAntennaExists;
	

	int numAntenna;

	
	
 	

	
	// ===> Attribute basebandName, which is optional
	
	
	bool basebandNameExists;
	

	std::vector<BasebandNameMod::BasebandName > basebandName;

	
	
 	

	
	// ===> Attribute numBaseband, which is optional
	
	
	bool numBasebandExists;
	

	int numBaseband;

	
	
 	

	
	// ===> Attribute interval, which is optional
	
	
	bool intervalExists;
	

	Interval interval;

	
	
 	

	
	// ===> Attribute dValue, which is optional
	
	
	bool dValueExists;
	

	double dValue;

	
	
 	

	
	// ===> Attribute vdValue, which is optional
	
	
	bool vdValueExists;
	

	std::vector<double > vdValue;

	
	
 	

	
	// ===> Attribute vvdValues, which is optional
	
	
	bool vvdValuesExists;
	

	std::vector<std::vector<double > > vvdValues;

	
	
 	

	
	// ===> Attribute llValue, which is optional
	
	
	bool llValueExists;
	

	int64_t llValue;

	
	
 	

	
	// ===> Attribute vllValue, which is optional
	
	
	bool vllValueExists;
	

	std::vector<int64_t > vllValue;

	
	
 	

	
	// ===> Attribute vvllValue, which is optional
	
	
	bool vvllValueExists;
	

	std::vector<std::vector<int64_t > > vvllValue;

	
	
 	

	////////////////////////////////
	// Extrinsic Table Attributes //
	////////////////////////////////
	
	
	// ===> Attribute antennaId, which is optional
	
	
	bool antennaIdExists;
	

	std::vector<Tag>  antennaId;

	
	
 	

	///////////
	// Links //
	///////////
	
	
		


	

	
/*
	////////////////////////////////////////////////////////////
	// binary-deserialization material from an EndianIStream  //
	////////////////////////////////////////////////////////////
	std::map<std::string, AnnotationAttributeFromBin> fromBinMethods;
void annotationIdFromBin( EndianIStream& eis);
void timeFromBin( EndianIStream& eis);
void issueFromBin( EndianIStream& eis);
void detailsFromBin( EndianIStream& eis);

void numAntennaFromBin( EndianIStream& eis);
void basebandNameFromBin( EndianIStream& eis);
void numBasebandFromBin( EndianIStream& eis);
void intervalFromBin( EndianIStream& eis);
void dValueFromBin( EndianIStream& eis);
void vdValueFromBin( EndianIStream& eis);
void vvdValuesFromBin( EndianIStream& eis);
void llValueFromBin( EndianIStream& eis);
void vllValueFromBin( EndianIStream& eis);
void vvllValueFromBin( EndianIStream& eis);
void antennaIdFromBin( EndianIStream& eis);

*/
	
	///////////////////////////////////
	// text-deserialization material //
	///////////////////////////////////
	std::map<std::string, AnnotationAttributeFromText> fromTextMethods;
	
void annotationIdFromText (const string & s);
	
	
void timeFromText (const string & s);
	
	
void issueFromText (const string & s);
	
	
void detailsFromText (const string & s);
	

	
void numAntennaFromText (const string & s);
	
	
void basebandNameFromText (const string & s);
	
	
void numBasebandFromText (const string & s);
	
	
void intervalFromText (const string & s);
	
	
void dValueFromText (const string & s);
	
	
void vdValueFromText (const string & s);
	
	
void vvdValuesFromText (const string & s);
	
	
void llValueFromText (const string & s);
	
	
void vllValueFromText (const string & s);
	
	
void vvllValueFromText (const string & s);
	
	
void antennaIdFromText (const string & s);
	
	
	
	/**
	 * Serialize this into a stream of bytes written to an EndianOSStream.
	 * @param eoss the EndianOSStream to be written to
	 */
	 void toBin(EndianOSStream& eoss);
	 	 
	 /**
	  * Deserialize a stream of bytes read from an EndianIStream to build a PointingRow.
	  * @param eiss the EndianIStream to be read.
	  * @param table the AnnotationTable to which the row built by deserialization will be parented.
	  * @param attributesSeq a vector containing the names of the attributes . The elements order defines the order 
	  * in which the attributes are written in the binary serialization.

	 static AnnotationRow* fromBin(EndianIStream& eis, AnnotationTable& table, const std::vector<std::string>& attributesSeq);	 
		*/
};

} // End namespace asdm

#endif /* Annotation_CLASS */
