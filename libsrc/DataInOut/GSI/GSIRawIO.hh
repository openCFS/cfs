// Definition of the GSIRawIO class

#ifndef GSIRawIO_class
#define GSIRawIO_class

#include <cstdio>

#include "GSIBaseIO.hh"
#include "GSIIOException.hh"

class GSIRawIO : public GSIBaseIO
{
 public:
  GSIRawIO(FILE *file_read, FILE *file_write);

  // the Send Message function
  virtual void writeMsg(const std::string& s) throw(GSIIOException);

  //!the function to recieve the Message; Returns a pointer to the string
  virtual void readMsg(std::string& s) throw(GSIIOException);

  //!the function to recieve integers
  virtual int readInt() throw(GSIIOException);
  
  virtual void readIntVector(std::vector<int>& vec) throw(GSIIOException); 

  virtual u_int readUInt() throw(GSIIOException);
  
  virtual void readUIntVector(std::vector<u_int>& vec) throw(GSIIOException); 

  virtual float readFloat() throw(GSIIOException);
  
  virtual void readFloatVector(std::vector<float>& vec) throw(GSIIOException); 

  virtual double readDouble() throw(GSIIOException);
  
  virtual void readDoubleVector(std::vector<double>& vec) throw(GSIIOException); 

  //!the function to recieve the vector;
  //  virtual data readVector();

  //!the function to Send integers
  virtual void  writeInt(const int buff) throw(GSIIOException);
  
  virtual void writeIntArray(const int *array, const int size) throw(GSIIOException);

  virtual void writeIntVector(const std::vector<int>& vec) throw(GSIIOException);

  virtual void  writeUInt(const u_int buff) throw(GSIIOException);
  
  virtual void writeUIntArray(const u_int *array, const int size) throw(GSIIOException);

  virtual void writeUIntVector(const std::vector<u_int>& vec) throw(GSIIOException);

  virtual void  writeFloat(const float buff) throw(GSIIOException);
  
  virtual void writeFloatArray(const float *array, const int size) throw(GSIIOException);

  virtual void writeFloatVector(const std::vector<float>& vec) throw(GSIIOException);

  virtual void  writeDouble(const double buff) throw(GSIIOException);
  
  virtual void writeDoubleArray(const double *array, const int size) throw(GSIIOException);
  
  virtual void writeDoubleVector(const std::vector<double>& vec) throw(GSIIOException);
  //!the function to Send the vector;
  //  virtual void writeVector(data outgoing_xdr_struct);
};


#endif
