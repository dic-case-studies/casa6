
/***
 * Framework independent implementation file for tableiterator...
 *
 * Implement the tableiterator component here.
 * 
 * // TODO: WRITE YOUR DESCRIPTION HERE! 
 *
 * @author
 * @version 
 ***/

#include <iostream>
#include <tableiterator_cmpt.h>
#include <casa/aips.h>

using namespace std;
using namespace casacore;
using namespace casa;

using namespace casacore;
namespace casac {

tableiterator::tableiterator()
{

}

tableiterator::~tableiterator()
{

}

::casac::record*
tableiterator::table()
{

    // TODO : IMPLEMENT ME HERE !
    return 0;
}

bool
tableiterator::reset()
{

    // TODO : IMPLEMENT ME HERE !
  return false;
}

bool
tableiterator::next()
{

    // TODO : IMPLEMENT ME HERE !
  return false;
}

bool
tableiterator::terminate()
{

    // TODO : IMPLEMENT ME HERE !
  return false;
}

bool
tableiterator::done()
{

    // TODO : IMPLEMENT ME HERE !
  return false;
}

} // casac namespace

