
/***
 * Framework independent implementation file for tableindex...
 *
 * Implement the tableindex component here.
 * 
 * // TODO: WRITE YOUR DESCRIPTION HERE! 
 *
 * @author
 * @version 
 ***/

#include <iostream>
#include <tableindex_cmpt.h>
#include <casa/aips.h>

using namespace std;
using namespace casacore;
using namespace casa;

using namespace casacore;
namespace casac {

tableindex::tableindex()
{

}

tableindex::~tableindex()
{

}

bool
tableindex::set(const ::casac::record& /*tab*/, const std::vector<std::string>& /*columns*/, const bool /*sort*/)
{

    // TODO : IMPLEMENT ME HERE !
  return false;
}

bool
tableindex::isunique()
{

    // TODO : IMPLEMENT ME HERE !
  return false;
}

bool
tableindex::setchanged(const std::vector<std::string>& /*columns*/)
{

    // TODO : IMPLEMENT ME HERE !
  return false;
}

long
tableindex::rownr(const ::casac::record& /*key*/)
{

    long rstat(0);
    return rstat;
    // TODO : IMPLEMENT ME HERE !
}

std::vector<long>
tableindex::rownrs(const ::casac::record& /*key*/, const ::casac::record& /*upperkey*/, const bool /*lowerincl*/, const bool /*upperincl*/)
{

    // TODO : IMPLEMENT ME HERE !
   return std::vector<long>(0);
}

bool
tableindex::close()
{

    // TODO : IMPLEMENT ME HERE !
  return false;
}

bool
tableindex::done()
{

    // TODO : IMPLEMENT ME HERE !
  return false;
}

} // casac namespace

