
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
 * File DelayModelFixedParametersTable.h
 */
 
#ifndef DelayModelFixedParametersTable_CLASS
#define DelayModelFixedParametersTable_CLASS

#include <string>
#include <vector>
#include <map>



	
#include <alma/ASDM/Speed.h>
	

	
#include <alma/ASDM/AngularRate.h>
	

	
#include <alma/ASDM/Length.h>
	

	
#include <alma/ASDM/Tag.h>
	




	

	

	

	

	

	

	

	

	

	

	

	

	

	

	

	

	

	

	



#include <alma/ASDM/ConversionException.h>
#include <alma/ASDM/DuplicateKey.h>
#include <alma/ASDM/UniquenessViolationException.h>
#include <alma/ASDM/NoSuchRow.h>
#include <alma/ASDM/DuplicateKey.h>


#ifndef WITHOUT_ACS
#include <asdmIDLC.h>
#endif

#include <alma/ASDM/Representable.h>

#include <pthread.h>

namespace asdm {

//class asdm::ASDM;
//class asdm::DelayModelFixedParametersRow;

class ASDM;
class DelayModelFixedParametersRow;
/**
 * The DelayModelFixedParametersTable class is an Alma table.
 * <BR>
 * 
 * \par Role
 * 
 * <BR>
 
 * Generated from model's revision "-1", branch ""
 *
 * <TABLE BORDER="1">
 * <CAPTION> Attributes of DelayModelFixedParameters </CAPTION>
 * <TR BGCOLOR="#AAAAAA"> <TH> Name </TH> <TH> Type </TH> <TH> Expected shape  </TH> <TH> Comment </TH></TR>
 
 * <TR> <TH BGCOLOR="#CCCCCC" colspan="4" align="center"> Key </TD></TR>
	
 * <TR>
 		
 * <TD><I> delayModelFixedParametersId </I></TD>
 		 
 * <TD> Tag</TD>
 * <TD> &nbsp; </TD>
 * <TD> &nbsp;identifies a unique row in the table. </TD>
 * </TR>
	


 * <TR> <TH BGCOLOR="#CCCCCC"  colspan="4" valign="center"> Value <br> (Mandatory) </TH></TR>
	
 * <TR>
 * <TD> delayModelVersion </TD> 
 * <TD> std::string </TD>
 * <TD>  &nbsp;  </TD> 
 * <TD> &nbsp; should include the name of the software and its version.  Something like  "CALC v11" or "VDT v1.0" or "MODEST v2.1".  </TD>
 * </TR>
	
 * <TR>
 * <TD> execBlockId </TD> 
 * <TD> Tag </TD>
 * <TD>  &nbsp;  </TD> 
 * <TD> &nbsp;refers to a unique row of the ExecBlock table. </TD>
 * </TR>
	


 * <TR> <TH BGCOLOR="#CCCCCC"  colspan="4" valign="center"> Value <br> (Optional) </TH></TR>
	
 * <TR>
 * <TD> gaussConstant</TD> 
 * <TD> AngularRate </TD>
 * <TD>  &nbsp; </TD>
 * <TD>&nbsp; the Gauss gravitational constant (should be of order \f$ 1.720209895.10^{-2} rad/d \f$ but in SI units of \f$ rad s^{-1} \f$).  </TD>
 * </TR>
	
 * <TR>
 * <TD> newtonianConstant</TD> 
 * <TD> double </TD>
 * <TD>  &nbsp; </TD>
 * <TD>&nbsp; the newtonian constant of gravitation (should be of order \f$ 6.67259.10^{-11} m^3  kg^{-1}  s^2 \f$).  </TD>
 * </TR>
	
 * <TR>
 * <TD> gravity</TD> 
 * <TD> double </TD>
 * <TD>  &nbsp; </TD>
 * <TD>&nbsp; the gravity acceleration in \f$ m s^{-2} \f$. </TD>
 * </TR>
	
 * <TR>
 * <TD> earthFlattening</TD> 
 * <TD> double </TD>
 * <TD>  &nbsp; </TD>
 * <TD>&nbsp; the ratio of equatorial to polar radii. </TD>
 * </TR>
	
 * <TR>
 * <TD> earthRadius</TD> 
 * <TD> Length </TD>
 * <TD>  &nbsp; </TD>
 * <TD>&nbsp; the earth equatorial radius in \f$ m \f$. </TD>
 * </TR>
	
 * <TR>
 * <TD> moonEarthMassRatio</TD> 
 * <TD> double </TD>
 * <TD>  &nbsp; </TD>
 * <TD>&nbsp;  </TD>
 * </TR>
	
 * <TR>
 * <TD> ephemerisEpoch</TD> 
 * <TD> std::string </TD>
 * <TD>  &nbsp; </TD>
 * <TD>&nbsp; should always be 'J2000'. </TD>
 * </TR>
	
 * <TR>
 * <TD> earthTideLag</TD> 
 * <TD> double </TD>
 * <TD>  &nbsp; </TD>
 * <TD>&nbsp;  </TD>
 * </TR>
	
 * <TR>
 * <TD> earthGM</TD> 
 * <TD> double </TD>
 * <TD>  &nbsp; </TD>
 * <TD>&nbsp; the earth gravitation constant in \f$ m^3 s^{-2} \f$. </TD>
 * </TR>
	
 * <TR>
 * <TD> moonGM</TD> 
 * <TD> double </TD>
 * <TD>  &nbsp; </TD>
 * <TD>&nbsp; the moon gravitation constant in \f$ m^3 s^{-2} \f$. </TD>
 * </TR>
	
 * <TR>
 * <TD> sunGM</TD> 
 * <TD> double </TD>
 * <TD>  &nbsp; </TD>
 * <TD>&nbsp; the sun gravitation constant in \f$ m^3 s^{-2} \f$. </TD>
 * </TR>
	
 * <TR>
 * <TD> loveNumberH</TD> 
 * <TD> double </TD>
 * <TD>  &nbsp; </TD>
 * <TD>&nbsp; the earth global Love number H. </TD>
 * </TR>
	
 * <TR>
 * <TD> loveNumberL</TD> 
 * <TD> double </TD>
 * <TD>  &nbsp; </TD>
 * <TD>&nbsp; the earth global Love number L. </TD>
 * </TR>
	
 * <TR>
 * <TD> precessionConstant</TD> 
 * <TD> AngularRate </TD>
 * <TD>  &nbsp; </TD>
 * <TD>&nbsp; the general precession constant in \f$ arcsec \  s^{-1} \f$. </TD>
 * </TR>
	
 * <TR>
 * <TD> lightTime1AU</TD> 
 * <TD> double </TD>
 * <TD>  &nbsp; </TD>
 * <TD>&nbsp; the light time for 1 AU in seconds. </TD>
 * </TR>
	
 * <TR>
 * <TD> speedOfLight</TD> 
 * <TD> Speed </TD>
 * <TD>  &nbsp; </TD>
 * <TD>&nbsp; the speed of light in \f$ m s^{-1} \f$.
 </TD>
 * </TR>
	
 * <TR>
 * <TD> delayModelFlags</TD> 
 * <TD> std::string </TD>
 * <TD>  &nbsp; </TD>
 * <TD>&nbsp; the delay model switches. </TD>
 * </TR>
	

 * </TABLE>
 */
class DelayModelFixedParametersTable : public Representable {
	friend class ASDM;

public:


	/**
	 * Return the list of field names that make up key key
	 * as an array of strings.
	 * @return a vector of string.
	 */	
	static const std::vector<std::string>& getKeyName();


	virtual ~DelayModelFixedParametersTable();
	
	/**
	 * Return the container to which this table belongs.
	 *
	 * @return the ASDM containing this table.
	 */
	ASDM &getContainer() const;
	
	/**
	 * Return the number of rows in the table.
	 *
	 * @return the number of rows in an unsigned int.
	 */
	unsigned int size() const;
	
	/**
	 * Return the name of this table.
	 *
	 * This is a instance method of the class.
	 *
	 * @return the name of this table in a string.
	 */
	std::string getName() const;
	
	/**
	 * Return the name of this table.
	 *
	 * This is a static method of the class.
	 *
	 * @return the name of this table in a string.
	 */
	static std::string name() ;	
	
	/**
	 * Return the version information about this table.
	 *
	 */
	 std::string getVersion() const ;
	
	/**
	 * Return the names of the attributes of this table.
	 *
	 * @return a vector of string
	 */
	 static const std::vector<std::string>& getAttributesNames();

	/**
	 * Return the default sorted list of attributes names in the binary representation of the table.
	 *
	 * @return a const reference to a vector of string
	 */
	 static const std::vector<std::string>& defaultAttributesNamesInBin();
	 
	/**
	 * Return this table's Entity.
	 */
	Entity getEntity() const;

	/**
	 * Set this table's Entity.
	 * @param e An entity. 
	 */
	void setEntity(Entity e);
		
	/**
	 * Produces an XML representation conform
	 * to the schema defined for DelayModelFixedParameters (DelayModelFixedParametersTable.xsd).
	 *
	 * @returns a string containing the XML representation.
	 * @throws ConversionException
	 */
	std::string toXML()  ;

#ifndef WITHOUT_ACS
	// Conversion Methods
	/**
	 * Convert this table into a DelayModelFixedParametersTableIDL CORBA structure.
	 *
	 * @return a pointer to a DelayModelFixedParametersTableIDL
	 */
	asdmIDL::DelayModelFixedParametersTableIDL *toIDL() ;
	
	/**
	 * Fills the CORBA data structure passed in parameter
	 * with the content of this table.
	 *
	 * @param x a reference to the asdmIDL::DelayModelFixedParametersTableIDL to be populated
	 * with the content of this.
	 */
	 void toIDL(asdmIDL::DelayModelFixedParametersTableIDL& x) const;
	 
#endif

#ifndef WITHOUT_ACS
	/**
	 * Populate this table from the content of a DelayModelFixedParametersTableIDL Corba structure.
	 *
	 * @throws DuplicateKey Thrown if the method tries to add a row having a key that is already in the table.
	 * @throws ConversionException
	 */	
	void fromIDL(asdmIDL::DelayModelFixedParametersTableIDL x) ;
#endif
	
	//
	// ====> Row creation.
	//
	
	/**
	 * Create a new row with default values.
	 * @return a pointer on a DelayModelFixedParametersRow
	 */
	DelayModelFixedParametersRow *newRow();
	
	
	/**
	 * Create a new row initialized to the specified values.
	 * @return a pointer on the created and initialized row.
	
 	 * @param delayModelVersion
	
 	 * @param execBlockId
	
     */
	DelayModelFixedParametersRow *newRow(std::string delayModelVersion, Tag execBlockId);
	


	/**
	 * Create a new row using a copy constructor mechanism.
	 * 
	 * The method creates a new DelayModelFixedParametersRow owned by this. Each attribute of the created row 
	 * is a (deep) copy of the corresponding attribute of row. The method does not add 
	 * the created row to this, its simply parents it to this, a call to the add method
	 * has to be done in order to get the row added (very likely after having modified
	 * some of its attributes).
	 * If row is null then the method returns a new DelayModelFixedParametersRow with default values for its attributes. 
	 *
	 * @param row the row which is to be copied.
	 */
	 DelayModelFixedParametersRow *newRow(DelayModelFixedParametersRow *row); 

	//
	// ====> Append a row to its table.
	//

	
	
	
	/** 
	 * Add a row.
	 * If there table contains a row whose key's fields are equal
	 * to x's ones then return a pointer on this row (i.e. no actual insertion is performed) 
	 * otherwise add x to the table and return x.
	 * @param x . A pointer on the row to be added.
 	 * @returns a pointer to a DelayModelFixedParametersRow.	 
	 */	 
	 
 	 DelayModelFixedParametersRow* add(DelayModelFixedParametersRow* x) ;



	//
	// ====> Methods returning rows.
	//
		
	/**
	 * Get a collection of pointers on the rows of the table.
	 * @return Alls rows in a vector of pointers of DelayModelFixedParametersRow. The elements of this vector are stored in the order 
	 * in which they have been added to the DelayModelFixedParametersTable.
	 */
	std::vector<DelayModelFixedParametersRow *> get() ;
	
	/**
	 * Get a const reference on the collection of rows pointers internally hold by the table.
	 * @return A const reference of a vector of pointers of DelayModelFixedParametersRow. The elements of this vector are stored in the order 
	 * in which they have been added to the DelayModelFixedParametersTable.
	 *
	 */
	 const std::vector<DelayModelFixedParametersRow *>& get() const ;
	


 
	
	/**
 	 * Returns a DelayModelFixedParametersRow* given a key.
 	 * @return a pointer to the row having the key whose values are passed as parameters, or 0 if
 	 * no row exists for that key.
	
	 * @param delayModelFixedParametersId
	
 	 *
	 */
 	DelayModelFixedParametersRow* getRowByKey(Tag delayModelFixedParametersId);

 	 	



	/**
 	 * Look up the table for a row whose all attributes  except the autoincrementable one 
 	 * are equal to the corresponding parameters of the method.
 	 * @return a pointer on this row if any, null otherwise.
 	 *
			
 	 * @param delayModelVersion
 	 		
 	 * @param execBlockId
 	 		 
 	 */
	DelayModelFixedParametersRow* lookup(std::string delayModelVersion, Tag execBlockId); 


	void setUnknownAttributeBinaryReader(const std::string& attributeName, BinaryAttributeReaderFunctor* barFctr);
	BinaryAttributeReaderFunctor* getUnknownAttributeBinaryReader(const std::string& attributeName) const;

private:

	/**
	 * Create a DelayModelFixedParametersTable.
	 * <p>
	 * This constructor is private because only the
	 * container can create tables.  All tables must know the container
	 * to which they belong.
	 * @param container The container to which this table belongs.
	 */ 
	DelayModelFixedParametersTable (ASDM & container);

	ASDM & container;
	
	bool archiveAsBin; // If true archive binary else archive XML
	bool fileAsBin ; // If true file binary else file XML	
	
	std::string version ; 
	
	Entity entity;
	

	
	

	// A map for the autoincrementation algorithm
	std::map<std::string,int>  noAutoIncIds;
	void autoIncrement(std::string key, DelayModelFixedParametersRow* x);


	/**
	 * If this table has an autoincrementable attribute then check if *x verifies the rule of uniqueness and throw exception if not.
	 * Check if *x verifies the key uniqueness rule and throw an exception if not.
	 * Append x to its table.
	 * @throws DuplicateKey
	 
	 * @throws UniquenessViolationException
	 
	 */
	DelayModelFixedParametersRow* checkAndAdd(DelayModelFixedParametersRow* x, bool skipCheckUniqueness=false) ;
	
	/**
	 * Brutally append an DelayModelFixedParametersRow x to the collection of rows already stored in this table. No uniqueness check is done !
	 *
	 * @param DelayModelFixedParametersRow* x a pointer onto the DelayModelFixedParametersRow to be appended.
	 */
	 void append(DelayModelFixedParametersRow* x) ;
	 
	/**
	 * Brutally append an DelayModelFixedParametersRow x to the collection of rows already stored in this table. No uniqueness check is done !
	 *
	 * @param DelayModelFixedParametersRow* x a pointer onto the DelayModelFixedParametersRow to be appended.
	 */
	 void addWithoutCheckingUnique(DelayModelFixedParametersRow* x) ;
	 
	 



// A data structure to store the pointers on the table's rows.

// In all cases we maintain a private vector of DelayModelFixedParametersRow s.
   std::vector<DelayModelFixedParametersRow * > privateRows;
   

			
	std::vector<DelayModelFixedParametersRow *> row;

	
	void error() ; //throw(ConversionException);

	
	/**
	 * Populate this table from the content of a XML document that is required to
	 * be conform to the XML schema defined for a DelayModelFixedParameters (DelayModelFixedParametersTable.xsd).
	 * @throws ConversionException
	 * 
	 */
	void fromXML(std::string& xmlDoc) ;
		
	std::map<std::string, BinaryAttributeReaderFunctor *> unknownAttributes2Functors;

	/**
	  * Private methods involved during the build of this table out of the content
	  * of file(s) containing an external representation of a DelayModelFixedParameters table.
	  */
	void setFromMIMEFile(const std::string& directory);
	/*
	void openMIMEFile(const std::string& directory);
	*/
	void setFromXMLFile(const std::string& directory);
	
		 /**
	 * Serialize this into a stream of bytes and encapsulates that stream into a MIME message.
	 * @returns a string containing the MIME message.
	 *
	 * @param byteOrder a const pointer to a static instance of the class ByteOrder.
	 * 
	 */
	std::string toMIME(const asdm::ByteOrder* byteOrder=asdm::ByteOrder::Machine_Endianity);
  
	
   /** 
     * Extracts the binary part of a MIME message and deserialize its content
	 * to fill this with the result of the deserialization. 
	 * @param mimeMsg the string containing the MIME message.
	 * @throws ConversionException
	 */
	 void setFromMIME(const std::string & mimeMsg);
	
	/**
	  * Private methods involved during the export of this table into disk file(s).
	  */
	std::string MIMEXMLPart(const asdm::ByteOrder* byteOrder=asdm::ByteOrder::Machine_Endianity);
	
	/**
	  * Stores a representation (binary or XML) of this table into a file.
	  *
	  * Depending on the boolean value of its private field fileAsBin a binary serialization  of this (fileAsBin==true)  
	  * will be saved in a file "DelayModelFixedParameters.bin" or an XML representation (fileAsBin==false) will be saved in a file "DelayModelFixedParameters.xml".
	  * The file is always written in a directory whose name is passed as a parameter.
	 * @param directory The name of directory  where the file containing the table's representation will be saved.
	  * 
	  */
	  void toFile(std::string directory);
	  
	  /**
	   * Load the table in memory if necessary.
	   */
	  bool loadInProgress;
	  void checkPresenceInMemory() {
		if (!presentInMemory && !loadInProgress) {
			loadInProgress = true;
			setFromFile(getContainer().getDirectory());
			presentInMemory = true;
			loadInProgress = false;
	  	}
	  }
	/**
	 * Reads and parses a file containing a representation of a DelayModelFixedParametersTable as those produced  by the toFile method.
	 * This table is populated with the result of the parsing.
	 * @param directory The name of the directory containing the file te be read and parsed.
	 * @throws ConversionException If any error occurs while reading the 
	 * files in the directory or parsing them.
	 *
	 */
	 void setFromFile(const std::string& directory);	
 
};

} // End namespace asdm

#endif /* DelayModelFixedParametersTable_CLASS */
