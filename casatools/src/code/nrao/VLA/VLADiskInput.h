//# VLADiskInput.h: This class reads VLA archive files from a Disk
//# Copyright (C) 1999,2000
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
//# $Id: VLADiskInput.h,v 19.3.14.1 2006/02/12 05:56:12 wyoung Exp $

#ifndef NRAO_VLADISKINPUT_H
#define NRAO_VLADISKINPUT_H

#include <casacore/casa/aips.h>
#include <casacore/casa/IO/RegularFileIO.h>
#include <nrao/VLA/VLAArchiveInput.h>

#include <casacore/casa/namespace.h>
namespace casacore{

class Path;
}

namespace casa { //# NAMESPACE CASA - BEGIN
} //# NAMESPACE CASA - END


// <summary>This class reads VLA archive records from a Disk file</summary>

// <reviewed reviewer="" date="" tests="" demos="">

// <prerequisite>
// <ol>
//   <li> The IO Module
// </ol>
// </prerequisite>
//
// <etymology>
// This class is designed to reads VLA archive records from a Disk file
// </etymology>
//
// <synopsis> 

// This class is designed to read VLA archive data.  The data may be read from
// a disk, tape drive or any other data source supported by the IO module.  A
// call to the operator++() function assembles the next reconstructed VLA
// archive data record from the input. A reference to this data can be
// obtained using the logicalRecord function.
//
// Refer to the "VLA Archive casacore::Data casacore::Format", VLA Computer Memorandum 186
// by G.C. Hunt, K.P. Sowinski, and T.J. Bottomly; June 1993.
// (This is also available as AIPS++ note 159)
//
// The VLA archive records are always a multiple of 2048 bytes.  The
// record sizes were designed for use with magnetic tapes, so there is
// a maximum physical record size of 13*2048=26624 bytes.
//
// The low level class (blockio), that actually does the I/O, allows
// for a record (hereinafter chunk) size and for a input block size of
// a multiple of the chunk size.  The low level read operation tests
// for the number of bytes actually read from the device.
//
// The helper classes VlaDiskInput, VlaTapeInput, and VlaStdInput are
// designed to deal with the low level input from the devices in an
// analogous fashion to the ones used for casacore::FITS input.
//
// Since a read may be issued for an arbitrary number of bytes from a
// disk, the chunk multiple is arbitrary and may be used to tune the
// speed of operation.  There is an obvious trade-off between the
// block size created in the blockio class and the number of read
// operations.
//
// The story is quite different for tape input.  A read request for at
// least the maximum physical record size must be made to avoid loss of
// data.  Since a single tape record will be read with a single read
// operation, there is no point is having it any larger.  The chunk
// multiple must be exactly 13 so that the block size is 26624.
//
// The reconstitution algorithm is as follows:
//
// 1. Read a 2048 chunk from the input.
//
// The first two 16-bit integers should contain the values 1 and n,
// where n is the number of "physical records" in the current "logical  
// record."  (If the first value is not 1, then the chunk is rejected
// and a new one read until the first 16-bit value is 1.)  These two
// values are not part of the reconstituted "logical record."
//
// 2. The next 32-bit integer contains the length of the "logical
// record" in 16-bit words.  The buffer is resized so that it can
// contain the whole reconstituted "logical record."
//
// 3. The remainder of the chunk is copied to the buffer.
//
// 4. Successive chunks are read from the input.
//
// The chunks are copied into the buffer until the "logical record"
// is complete.  For "logical records" longer than 26620 byte, this is
// not the whole story.  Each "physical record" contains maximum of 13
// chunks.  When the first "physical record" of 13 chunks has been read,
// the next chunk will be the first of the next "physical record."  The
// first two 16-bit integers will now be 2 and n, to indicate that this
// is the second "physical record" of the sequence.  These 4 bytes are
// decoded and the rest of this chunk is copied to the buffer.  And so
// on...
//
// An end-of-file condition on the input will cause record processing
// to be declared complete.  
// </synopsis> 
//
// <example>
// To open and read a VLA archive data file
// <code>
//    VLAArchiveInput *in;
//    casacore::Block <casacore::Char> *buff;
//    casacore::String fileName = " ";
//    casacore::String fileType = "tape";
//
//    if (fileType == casacore::String("tape")) {
//        in = new VLAArchiveInput(fileName.chars(), VLAArchiveInput::Tape);
//    } else {
//        in = new VLAArchiveInput(fileName.chars(), VLAArchiveInput::Disk);
//    }
//
//    casacore::uInt record = 0;
//    for (buff=&(in->next()); in->atEnd()==false; buff=&(in->next()), record++) {
//        cout << "casacore::Record" << record << endl;
//        // process record pointed to by buff
//    }
// </code>
//
// </example>
//
// <motivation>
// </motivation>
//
// <todo asof="">
// <ol>
//   <li> Bulletproofing - check for realistic buffer size (<1e6)
//   <li> Bulletproofing - check newn and newm on each read
//   <li> What happens after a single end-of-file on a tape drive?
//   <li> Add record skipping
//   <li> Should it work with stdin?  This is in place but not debugged.
// </ol>
// </todo>

class VLADiskInput: public VLAArchiveInput
{
public: 
  // Create an object that reads its data from the specified file on disk.
  // An exception is thrown if the file does not exist or cannot be read.
  VLADiskInput(const casacore::Path& fileName);

  // Use this method to attach the VLA online data repository
  // onlineFlag will typically be "online", but could be "online:-1"
  // for previous day's file.  The VLA online data repository typically
  // holds 14 days
  VLADiskInput(const casacore::String &onlineFlag);
    
  // The destructor closes the file.
  virtual ~VLADiskInput();

  // Reads the next logical record from specified IO source. Returns false if
  // there was a problem assembling the next record ie., it returns the value
  // of the hasData() member function.
  virtual casacore::Bool read();

  // Returns the number of bytes read so far.
  virtual casacore::uInt bytesRead();

  // Returns the size of the file in bytes.
  virtual casacore::uInt totalBytes();
  
private: 
  //# The default constructor is private and undefined
  VLADiskInput();
    
  //# The copy constructor is private and undefined
  VLADiskInput(const VLADiskInput& other);

  //# The assignment operator is private and undefined
  VLADiskInput& operator=(const VLADiskInput& other);

  //# Reads through a VLA archive looking for the first physical record in a
  //# logical record. Returns false if the first record could not be found. If
  //# It was found then this function also returns the number of physical
  //# records in this logical record.
  casacore::Bool findFirstRecord(casacore::Short& m);

  //# This object that provides the data input.
  casacore::RegularFileIO *itsFile;

  // The following privates are used for accessing the VLA online repository
  // casacore::Directory containing the visibility data
  casacore::Bool   onlineFill;
  casacore::String visDir; 
  casacore::Int oldDay;
  casacore::Int previousDay;

  casacore::String getTodaysFile(int relDay=0);
  casacore::Int whatsToday();
  void reattachCurrent();
};

#endif
