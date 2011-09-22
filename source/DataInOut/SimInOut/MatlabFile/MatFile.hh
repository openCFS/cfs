#ifndef MATFILE_H
#define MATFILE_H

#include "hdf5.h"
#include "hdf5_hl.h"

#include <string>

// Class for writing matlab's .mat files. These are HDF5 files with a simple (and pointless) header.
// Currently only simple matrix and vector variables are supported. TODO: cell arrays would be nice.
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




#endif // MATFILE_H
