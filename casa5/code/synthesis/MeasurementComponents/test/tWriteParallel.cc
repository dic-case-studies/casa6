
#include <casa/Arrays/Array.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/ArrayMath.h>
#include <synthesis/MeasurementComponents/WriteMSAlgorithm.h>
#include <msvis/MSVis/VisSetUtil.h>
#include <ms/MeasurementSets/MeasurementSet.h>


#include <casa/namespace.h>
using namespace casacore;
using namespace casa;
using namespace std;
extern casa::Applicator casa::applicator;

int main(Int argc, Char *argv[])
{
WriteMSAlgorithm *pw=new WriteMSAlgorithm();
casa::applicator.defineAlgorithm(pw);
  casa::applicator.init(argc, argv);
if(casa::applicator.isController()){
  for (int lala=0; lala < 2 ; ++lala){
    cerr << "doing loop  "<< lala << endl;
  WriteMSAlgorithm pwrite;
  if (argc<2) {
    cout <<"Usage: tWriteParallel ms-table-name  "<<endl;
    exit(0);
  }
  String msname(argv[1]);
  if(!Table::isWritable(msname))
    cout << msname << " is not a writable Table"<< endl;
  {
    //Make sure it has a model column
    MeasurementSet thems(msname, TableLock(TableLock::AutoNoReadLocking),Table::Update);
    VisSetUtil::addScrCols(thems, true, false, true, false);

  }
  
  
  cerr << "Number of procs: " << casa::applicator.numProcs() << endl;
  Int rank(0);
  Bool assigned; //(casa::casa::applicator.nextAvailProcess(pwrite, rank));
  Bool allDone(false);
  
    for (int k=0; k < 4; ++k){
      assigned=casa::applicator.nextAvailProcess(pwrite, rank);
      cerr << "assigned "<< assigned << endl;
      while(!assigned){
        rank = casa::applicator.nextProcessDone(pwrite, allDone);
        cerr << "while rank " << rank << endl;
        Int status;
        casa::applicator.get(status);
        cerr <<"STATUS " << status << endl;
        if(status != 0)
          cerr << k << " rank " << rank << " successful " << endl;
        else
          cerr << k << " rank " << rank << " failed " << endl;
        assigned = casa::applicator.nextAvailProcess(pwrite, rank);
       
      }
      
      ///send process to write the next 100 rows in ms
      Vector<Int> rowids(100);
      indgen(rowids, k*100);
      casa::applicator.put(msname);
      casa::applicator.put(rowids);
      casa::applicator.put(Complex(k));
      casa::applicator.apply(pwrite);

    }
    // Wait for all outstanding processes to return
    rank = casa::applicator.nextProcessDone(pwrite, allDone);
    while (!allDone) {
      Int status;
      if(casa::applicator.isSerial())
        casa::applicator.get(status);// get that extra put
      casa::applicator.get(status);
      cerr << "STATUSend " << status << endl;
      if(status != 0)
        cerr << "remainder rank " << rank << " successful " << endl;
      else
        cerr << "remainder rank " << rank << " failed " << endl;
     
      rank = casa::applicator.nextProcessDone(pwrite, allDone);
      if(casa::applicator.isSerial())
        allDone=true;
    }

  }
 }
  return(1);
}


