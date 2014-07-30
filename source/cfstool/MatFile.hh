// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef MATFILE_H
#define MATFILE_H

#include "hdf5.h"
#include "hdf5_hl.h"

#include <string>

namespace CoupledField {

  /** Class for writing Matlab 7.3 .mat files.
   *
   * The  .mat files  are normal  HDF5 files,  which have  a simple (and pointless)
   * prepended  512 byte header.  They can be viewed  in hdfview.
   * Currently only simple matrix and vector variables are supported.
   *
   * TODO: cell arrays would be nice.
   *
   * Source: http://www.mathworks.com/matlabcentral/fileexchange/27350-c-class-to-write-hdf5-mat-files */
  class MatFile
  {
  public:
    // Create an unopened mat file. Use Open after this.
    MatFile();
    // Create and open a mat file.
    explicit MatFile(const std::string& fn);
    
    // Close the mat file.
    ~MatFile();
    
    // Open a new mat file. If the file already exists it will be overwritten.
    bool Open(const std::string& fn);
    
    // Close the current file.
    bool Close();

    // Return true if the mat file is open and generally ok.
    bool IsOpen() const;
    
    // Write a matrix with the given name to the mat file. data is in row major format (i.e. idx = y*width + x).
    bool WriteMatrix(const std::string& varname, int nx, int ny, const double* data);
    
    // Write a vector.
    bool WriteVector(const std::string& varname, int nx, const double* data);
    
    // Convenience function.
    bool WriteValue(const std::string& varname, double v)
    {
      return WriteVector(varname, 1, &v);
    }
    
    
  private:
    
    // No copy.
    MatFile(const MatFile&);
    const MatFile& operator=(const MatFile&);
    
    bool PrependHeader(const char* filename);
    
    std::string filename;
    hid_t file; // -ve is invalid.
  };
  
} // end namespace CoupledField


#endif // MATFILE_H
