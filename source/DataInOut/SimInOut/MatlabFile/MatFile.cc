#include "MatFile.hh"

using namespace std;

#include <cstdio>
#include <cstring>

#include "General/exception.hh"

namespace CoupledField {

MatFile::MatFile() : file(-1)
{
}
MatFile::MatFile(const string& fn) : file(-1)
{
	Open(fn);
}

MatFile::~MatFile()
{
	Close();
}

bool MatFile::Open(const string& fn)
{
  herr_t status;
  hid_t fcpl;
  
  if (file >= 0)
    Close();
  
  filename = fn;
  
  // Create a copy  or instance of the File  Creation property list with
  // H5Pcreate.
  fcpl = H5Pcreate (H5P_FILE_CREATE);
  
  //Modify the property list with one of the File Creation property list
  //functions. For example, modify  the User Block with H5Pset_userblock
  //:
  status = H5Pset_userblock(fcpl, 512);
  if(status < 0) return false;  
 
  //Call H5Fcreate, passing in the  identifier of the property list that
  //was just modified. For example:
  
  file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, fcpl, H5P_DEFAULT);
  
  // Lastly, close the property list with H5Pclose:
  status = H5Pclose (fcpl);
  
  return file >= 0;
}

bool MatFile::Close()
{
	if (file >= 0)
	{
		H5Fclose(file);
		file = -1;
		return PrependHeader(filename.c_str());
	}
	return false;
}

bool MatFile::IsOpen() const
{
	return file >= 0;
}

bool MatFile::WriteMatrix(const string& varname, int nx, int ny, const double* data)
{
	if (file < 0)
		return false;

	if (nx <= 0 || ny <= 0)
		return true;

	hsize_t dimsf[2];
	dimsf[0] = nx;
	dimsf[1] = ny;
	return H5LTmake_dataset_double(file, varname.c_str(), 2, dimsf, data) >= 0;
}

bool MatFile::WriteVector(const string& varname, int nx, const double* data)
{
	if (file < 0)
		return false;

	if (nx <= 0)
		return true;

	hsize_t dimsf[1];
	dimsf[0] = nx;
	return H5LTmake_dataset_double(file, varname.c_str(), 1, dimsf, data) >= 0;
}



bool MatFile::PrependHeader(const char* filename)
{
	// Prepend the header.
	FILE* fd = fopen(filename, "r+b");
	if (!fd)
        {
		EXCEPTION("Couldn't open file for writing header.");
		return false;
	}

	char header[512];
	memset(header, 0, sizeof(header));
	sprintf(header, "MATLAB 7.3 format file written by MatFile class, by Tim Hutt"); // Do not make this longer than 116 chars (I think)
	header[124] = 0;
	header[125] = 2;
	header[126] = 'I';
	header[127] = 'M';

#if 0
	// Get file length.
	fseek(fd, 0, SEEK_END);
	size_t length = ftell(fd);
	fseek(fd, 0, SEEK_SET);

	// TODO: Do this properly without reading entire file into memory.
        size_t one_gb = 1024L*1024L*1024L*1;
	if (length > one_gb) // 1 GB.
	{
                fclose(fd);
		EXCEPTION("File too big to write header, sorry.");
		return false;
	}

	unsigned char* buffer = new unsigned char[length];
	if (fread(buffer, 1, length, fd) != length)
	{
                fclose(fd);
		delete[] buffer;
		EXCEPTION("Couldn't read file, sorry.");
		return false;
	}

	fseek(fd, 0, SEEK_SET);
#endif

	fwrite(header, 1, sizeof(header), fd);
        //	fwrite(buffer, 1, length, fd);
	fclose(fd);

        //	delete[] buffer;

	return true;
}

} // end namespace CoupledField
