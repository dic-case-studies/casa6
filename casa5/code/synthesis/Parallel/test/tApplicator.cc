#include <synthesis/Parallel/Applicator.h>
#include <synthesis/Parallel/Algorithm.h>
#include <synthesis/Parallel/PTransport.h>
#include <casa/Arrays/Vector.h>
#include <casa/BasicSL/Complex.h>
#include <casa/BasicSL/String.h>

#include <casa/namespace.h>
extern casa::Applicator casa::applicator;

const float VALINT = 3;
const float VALFLOAT = 3.141593;
const double VALDOUBLE = 3.141592654;

class TestAlgorithm : public casa::Algorithm {
   public :
    TestAlgorithm(bool multiproc): myName(String("Test Algorithm")),
                                   m_multiproc(multiproc){}
     ~TestAlgorithm(){};

      void get();
      void put();
      String &name(){return myName;}
   private :
      Int      one;
      Float    two;
      Double   three;
      Complex  four;
      DComplex five;
      String   six;
      Bool     seven;
      
      Vector<Int>      aOne;
      Vector<Float>    aTwo;
      Vector<Double>   aThree;
      Vector<Complex>  aFour;
      Vector<DComplex> aFive;
      Record           rec;
      void             task();
      String           myName;
    Int m_multiproc;
};

void TestAlgorithm::get(){

      cout << "In TestAlgorithm::get\n";

      // This sequence of gets must match the sequence of puts in the "applicator" below
      casa::applicator.get(one);
      cout << "got one " << one << '\n';
      casa::applicator.get(two);
      cout << "got two " << two << '\n';
      casa::applicator.get(three);
      cout << "got three " << three << '\n';
      casa::applicator.get(four);
      cout << "got four " << four << '\n';
      casa::applicator.get(five);
      cout << "got five " << five << '\n';
      casa::applicator.get(six);
      cout << "got six " << six << '\n';
      casa::applicator.get(seven);
      cout << "got seven " << seven << '\n';

      casa::applicator.get(aOne);
      cout << "got aOne " << '\n';
      casa::applicator.get(aTwo);
      cout << "got aTwo " << '\n';
      casa::applicator.get(aThree);
      cout << "got aThree " << '\n';
      casa::applicator.get(aFour);
      cout << "got aFour " << '\n';
      casa::applicator.get(aFive);
      cout << "got aFive " << '\n';

      casa::applicator.get(rec);
      cout << "got record: " << rec << '\n';


      DComplex val;
      unsigned int ncx = 100;
      for(unsigned int idx=0; idx<ncx; idx++) {
          casa::applicator.get(val);
      }
      cout << "got " << ncx << " complex, last val: " << val << '\n';

      Array<Float> farray;
      casa::applicator.get(farray);
      cout << "got casacore multidimensional array (float) of " << farray.ndim()
           << " dimensions. Sizes: " << farray.shape() << '\n';
      IPosition fpos(2,0,0);
      cout << " value of first element: " << farray(fpos) << ", last element: "
           << farray(farray.endPosition()) << '\n';

      Array<Double> darray;
      casa::applicator.get(darray);
      cout << "got casacore multidimensional array (Double) of " << darray.ndim()
           << " dimensions. Sizes: " << darray.shape() << '\n';
      IPosition dpos(4,0,0,0,0);
      cout << " value of first element: " << darray(dpos) << '\n';

      Array<Int> iarray;
      casa::applicator.get(iarray);
      cout << "got casacore multidimensional array (Int) of " << iarray.ndim()
           << " dimensions. Sizes: " << iarray.shape() << '\n';
      IPosition ipos(5,0,0,0,0,0);
      cout << " value of first element: " << iarray(ipos) << '\n';

      Array<Complex> carray;
      casa::applicator.get(carray);
      cout << "got casacore multidimensional array (Complex) of " << carray.ndim()
           << " dimensions. Sizes: " << carray.shape() << '\n';
      IPosition cpos(3,0,0,0);
      cout << " value of first element: " << carray(cpos) << '\n';

      Array<DComplex> dcarray;
      casa::applicator.get(dcarray);
      cout << "got casacore multidimensional array (DComplex) of " << dcarray.ndim()
           << " dimensions. Sizes: " << dcarray.shape() << '\n';
      IPosition dcpos(3,0,0,0);
      cout << " value of first element: " << dcarray(dcpos) << '\n';
}

void TestAlgorithm::put(){
      casa::applicator.put(true);
}

void TestAlgorithm::task(){
      cout << "Do work now!" << '\n';
}

// This serves at the moment as a basic test on the put/get methods of
// the "transports" defined in code/synthesis/Parallel/
// Not really clear what we want from the Algorithm and Applicator interfaces defined
// in code/synthesis/Parallel/ (May 2019).
//
// To run the test in MPI mode, use mpirun or mpicasa
// make unit_test_tApplicator
// mpirun ./tApplicator
// mpirun -n 32 ./tApplicater
// etc.
int main(Int argc, Char *argv[]){

   TestAlgorithm *testMe = new TestAlgorithm(casa::applicator.numProcs() > 1);
   casa::applicator.defineAlgorithm(testMe);
   casa::applicator.init(argc, argv);

   cout << "Number of procs: " << casa::applicator.numProcs() << '\n';

   if(casa::applicator.isController()){
      Int rank(1);
      casa::applicator.nextAvailProcess(testMe, rank);
      // nextAvailProcess does a put(int)!

      // This sequence of puts must match the sequence of gets in the "algorithm" above
      Int      one(1);
      Float    two(2.0f);
      Double   three(3.0);
      Complex  four(4.0,4.0);
      DComplex five(5.0,5.0);
      String   six("Six");
      Bool     seven(true);
      cout << "one " << one << '\n';
      cout << "two " << two << '\n';
      cout << "three " << three << '\n';
      cout << "four " << four << '\n';
      cout << "five " << five << '\n';
      cout << "six " << six << '\n';
      cout << "seven " << seven << '\n';

      casa::applicator.put(one);
      casa::applicator.put(two);
      casa::applicator.put(three);
      casa::applicator.put(four);
      casa::applicator.put(five);
      casa::applicator.put(six);
      casa::applicator.put(seven);

      Vector<Int>      aOne(3,one);
      Vector<Float>    aTwo(4,two);
      Vector<Double>   aThree(5,three);
      Vector<Complex>  aFour(6,four);
      Vector<DComplex> aFive(7,five);

      casa::applicator.put(aOne);
      casa::applicator.put(aTwo);
      casa::applicator.put(aThree);
      casa::applicator.put(aFour);
      casa::applicator.put(aFive);

      Record rec;
      rec.define("field_a", "example");
      rec.define("field_b", 3.45);
      rec.define("field_c", aFour);
      Record subrec;
      subrec.define("fs1", true);
      subrec.define("fs2", "bla1 bla2");
      subrec.define("fs3", aOne);
      rec.defineRecord("field_d", subrec);
      casa::applicator.put(rec);

      // send a bunch of items
      DComplex val(3.4, -5.6);
      for(int idx=0; idx<100; idx++) {
          casa::applicator.put(val);
      }

      // send casacore arrays of different types
      IPosition fdims(2,2,10);
      Array<Float> farray2d(fdims);
      farray2d = VALFLOAT;
      casa::applicator.put(farray2d);

      IPosition ddims(4,2,4,5,3);
      Array<Double> darray4d(ddims);
      darray4d = VALDOUBLE;
      casa::applicator.put(darray4d);

      IPosition idims(5,2,2,2,3,2);
      Array<Int> iarray5d(idims);
      iarray5d = VALINT;
      casa::applicator.put(iarray5d);

      IPosition cdims(3,3,4,2);
      Array<Complex> carray3d(cdims);
      carray3d = Complex(VALDOUBLE, -VALDOUBLE);
      casa::applicator.put(carray3d);

      IPosition dcdims(3,1,2,3);
      Array<DComplex> dcarray3d(dcdims);
      dcarray3d = DComplex(VALDOUBLE, -VALDOUBLE);
      casa::applicator.put(dcarray3d);

      casa::applicator.apply(testMe);
      Bool status;
      casa::applicator.get(status);
/*
      Bool allDone;
      applicator.nextProcessDone(testMe, allDone);
*/

      return 0;

   }
   return 0;
}
