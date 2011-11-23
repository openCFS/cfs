#include "MatFile.hh"

// This sample program shows how to  use the MatFile class to write Matlab 7.3
// .mat files, which are HDF5 based.
// Source: http://www.mathworks.com/matlabcentral/fileexchange/27350-c-class-to-write-hdf5-mat-files

// The  .mat files  are normal  HDF5 files,  which have  a prepended  512 byte
// header.  To view  them  in hdfview  one  can remove  the  header using  the
// following command:

// dd if=test.mat bs=512 skip=1 of=out.h5

int main(int argc, char** argv) 
{
  // Create new MatFile object.
  CoupledField::MatFile mf("test.mat");

  // Prepare some data
  double data[6];
  for(unsigned i=0; i<6; i++) 
  {
    data[i] = i+1;
  }  
  
  // Write a matrix and a vector to the data file.
  mf.WriteMatrix("A", 3, 2, data);
  mf.WriteVector("b", 3, data+2);
  
  // Close .mat file and automatically prepend a header.
  mf.Close();
}
