// Definition of the CFS_BaseIO class

#ifndef CFS_BaseIO_class
#define CFS_BaseIO_class

#include <sys/types.h>
#include <sys/time.h>

#ifdef __sgi
#include <stdio.h>
#else
#include <cstdio>
#endif

#include <string>
#include <vector>

#include "GSIIOException.hh"
#include "GSISocket.hh"

class GSIBaseIO
{
 public:

  //! Constructor for file IO
  GSIBaseIO(FILE* fr, FILE* fw) : file_read( fr ), file_write( fw ), sock_(NULL), timeout_(-1) {};
  
  //! Constructor for socket IO
  GSIBaseIO(GSISocket *sock, int timeout = -1) : file_read( NULL ), file_write( NULL ), sock_(sock), timeout_(timeout) {};

  //!
  virtual ~GSIBaseIO() {};

  // the Send Message function
  virtual void writeMsg(const std::string& s) throw(GSIIOException) = 0;

  //!the function to recieve the Message; Returns a pointer to the string
  virtual void readMsg(std::string& s) throw(GSIIOException) = 0;

  //!the function to recieve integers
  virtual int readInt() throw(GSIIOException) = 0;
  
  virtual void readIntVector(std::vector<int>& vec) throw(GSIIOException) = 0; 

  virtual u_int readUInt() throw(GSIIOException) = 0;
  
  virtual void readUIntVector(std::vector<u_int>& vec) throw(GSIIOException) = 0; 

  virtual float readFloat() throw(GSIIOException) = 0;
  
  virtual void readFloatVector(std::vector<float>& vec) throw(GSIIOException) = 0; 

  virtual double readDouble() throw(GSIIOException) = 0;
  
  virtual void readDoubleVector(std::vector<double>& vec) throw(GSIIOException) = 0; 

  //!the function to recieve the vector;
  //  virtual data readVector() = 0;

  //!the function to Send integers
  virtual void  writeInt(const int buff) throw(GSIIOException) = 0;
  
  virtual void writeIntArray(const int *array, const int size) throw(GSIIOException) = 0;

  virtual void writeIntVector(const std::vector<int>& vec) throw(GSIIOException) = 0;

  virtual void  writeUInt(const u_int buff) throw(GSIIOException) = 0;
  
  virtual void writeUIntArray(const u_int *array, const int size) throw(GSIIOException) = 0;

  virtual void writeUIntVector(const std::vector<u_int>& vec) throw(GSIIOException) = 0;


  virtual void  writeFloat(const float buff) throw(GSIIOException) = 0;
  
  virtual void writeFloatArray(const float *array, const int size) throw(GSIIOException) = 0;

  virtual void writeFloatVector(const std::vector<float>& vec) throw(GSIIOException) = 0;

  virtual void  writeDouble(const double buff) throw(GSIIOException) = 0;
  
  virtual void writeDoubleArray(const double *array, const int size) throw(GSIIOException) = 0;
  
  virtual void writeDoubleVector(const std::vector<double>& vec) throw(GSIIOException) = 0;

  virtual int ioOK() {return io_ok;};
  //!the function to Send the vector;
  //  virtual void writeVector(data outgoing_xdr_struct) = 0;

 protected:

  FILE* file_read;
  FILE* file_write;

  int desc_read;
  int desc_write;

  int timeout_;
  int io_ok;

  GSISocket* sock_;

  //  virtual bool mayRead() throw(GSIIOException) = 0;
  //  virtual bool mayWrite() throw(GSIIOException) = 0;

  void is_valid() throw(GSIIOException)
  {
    if(sock_ == NULL)
    {    
      if((file_read == NULL) && (file_write == NULL))
        throw GSIIOException("Read/Write Streams not valid!");
    }
    else
    {
      if(!sock_->is_valid())
        throw GSIIOException("CFSSocket not valid!");
    }
    
  }
  
};

GSIBaseIO& operator << (GSIBaseIO& io, const std::string& s ) throw(GSIIOException);
GSIBaseIO& operator >> (GSIBaseIO& io, std::string& s) throw(GSIIOException);

GSIBaseIO& operator << (GSIBaseIO& io, const int &i ) throw(GSIIOException);
GSIBaseIO& operator >> (GSIBaseIO& io, int& i) throw(GSIIOException);

GSIBaseIO& operator << (GSIBaseIO& io, const float &f ) throw(GSIIOException);
GSIBaseIO& operator >> (GSIBaseIO& io, float& f) throw(GSIIOException);

GSIBaseIO& operator << (GSIBaseIO& io, const double &d ) throw(GSIIOException);
GSIBaseIO& operator >> (GSIBaseIO& io, double& d) throw(GSIIOException);

GSIBaseIO& operator << (GSIBaseIO& io, const std::vector<int> &iv ) throw(GSIIOException);
GSIBaseIO& operator >> (GSIBaseIO& io, std::vector<int>& iv) throw(GSIIOException);

GSIBaseIO& operator << (GSIBaseIO& io, const std::vector<u_int> &uiv ) throw(GSIIOException);
GSIBaseIO& operator >> (GSIBaseIO& io, std::vector<u_int>& uiv) throw(GSIIOException);

GSIBaseIO& operator << (GSIBaseIO& io, const std::vector<float> &fv ) throw(GSIIOException);
GSIBaseIO& operator >> (GSIBaseIO& io, std::vector<float>& fv) throw(GSIIOException);

GSIBaseIO& operator << (GSIBaseIO& io, const std::vector<double> &dv ) throw(GSIIOException);
GSIBaseIO& operator >> (GSIBaseIO& io, std::vector<double>& dv) throw(GSIIOException);


#endif
