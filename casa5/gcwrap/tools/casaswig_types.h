#ifndef CASASWIG_TYPES_H
#define CASASWIG_TYPES_H 1

#include <stdarg.h>
#include <string>
#include <vector>
#include <stdcasa/Quantity.h>
#include <stdcasa/record.h>
#include <stdcasa/variant.h>

#define USING_NUMPY_ARRAYS 1

namespace casac {

struct complex {
	complex() {}
	complex(double arg0, double arg1): re(arg0), im(arg1){}
	double re;
	double im;
};

template<class T>
std::vector<T> initialize_vector(int count, T v1, ...) {
   va_list ap;
   va_start(ap, v1);
   std::vector<T> result(count);
   result[0] = v1;
   for ( int i=1; i < count; ++i ) {
       T val = va_arg(ap,T);
       result[i] = val;
   }
   return result;
}

//***
//*** In file included from /opt/casa/02/include/python2.7/unicodeobject.h:4:0,
//*** from /opt/casa/02/include/python2.7/Python.h:85,
//***    from coordsysPYTHON_wrap.cxx:174:
//*** ../../../../src/gcwrap/tools/casaswig_types.h: In function ‘std::vector<T> casac::initialize_vector(int, T, ...) [with T = bool]’:
//*** ../../../../src/gcwrap/tools/casaswig_types.h:29:26: warning: ‘bool’ is promoted to ‘int’ when passed through ‘...’
//***    T val = va_arg(ap,T);
//***                      ^
//*** ../../../../src/gcwrap/tools/casaswig_types.h:29:26: note: (so you should pass ‘int’ not ‘bool’ to ‘va_arg’)
//*** ../../../../src/gcwrap/tools/casaswig_types.h:29:26: note: if this code is reached, the program will abort
//***
template<> std::vector<bool> initialize_vector<bool>(int count, bool v1, ...);

struct BoolAry
{
  BoolAry( ) { }
  BoolAry(std::vector<bool> arg0, std::vector<int> arg1) : value(arg0), shape(arg1) { }
  std::vector<bool> value;
  std::vector<int> shape;

};

struct ComplexAry
{
  ComplexAry( ) { }
  ComplexAry(std::vector<casac::complex> arg0, std::vector<int> arg1) : value(arg0), shape(arg1) { }
  std::vector<casac::complex> value;
  std::vector<int> shape;

};

struct DoubleAry
{
  DoubleAry( ) { }
  DoubleAry(std::vector<double> arg0, std::vector<int> arg1) : value(arg0), shape(arg1) { }
  std::vector<double> value;
  std::vector<int> shape;

};

struct IntAry
{
  IntAry( ) { }
  IntAry(std::vector<long> arg0, std::vector<int> arg1) : value(arg0), shape(arg1) { }
  std::vector<long> value;
  std::vector<int> shape;

};

struct uIntAry
{
  uIntAry( ) { }
  uIntAry(std::vector<unsigned long> arg0, std::vector<int> arg1) : value(arg0), shape(arg1) { }
  std::vector<unsigned long> value;
  std::vector<int> shape;

};


struct StringAry
{
  StringAry( ) { }
  StringAry(std::vector<std::string> arg0, std::vector<int> arg1) : value(arg0), shape(arg1) { }
  std::vector<std::string> value;
  std::vector<int> shape;

};
typedef  std::vector<record> RecordVec;
typedef StringAry StringVec;
typedef IntAry IntVec;
typedef uIntAry uIntVec;
typedef DoubleAry DoubleVec;
typedef BoolAry BoolVec;
typedef ComplexAry ComplexVec;
/*
typedef  std::vector<std::string> StringVec;
typedef  std::vector<double> DoubleVec;
typedef  std::vector<long> IntVec;
typedef  std::vector<bool> BoolVec;
typedef  std::vector<std::complex<double> > ComplexVec;
*/
typedef std::string MDirection;
typedef std::string MRadialVelocity;
typedef std::string MPosition;
typedef std::string Region;

} // casac namespace
#endif
