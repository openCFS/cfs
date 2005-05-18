// Definition of the GSI::BaseIO class

#ifndef GSI_BASEIO
#define GSI_BASEIO

#include <sys/types.h>
#include <sys/time.h>

#ifdef __sgi
#include <stdio.h>
#else
#include <cstdio>
#endif

#include <string>
#include <vector>

#include "GSITypeDefs.hh"
#include "GSIIOException.hh"
#include "GSISocket.hh"

namespace GridlibSocketInterface
{

class BaseIO
{
 public:

  //! Constructor for file IO
  BaseIO(FILE* fr, FILE* fw) : file_read( fr ), file_write( fw ), sock_(NULL), timeout_(-1) {};
  
  //! Constructor for socket IO
  BaseIO(Socket *sock, int32 timeout = -1) : file_read( NULL ), file_write( NULL ), sock_(sock), timeout_(timeout) {};

  //!
  virtual ~BaseIO() {};

  virtual int32 EOFRead() throw(IOException) = 0;

  // the Send Message function
  virtual void writeMsg(const std::string& s) throw(IOException) = 0;

  //!the function to recieve the Message; Returns a pointer to the string
  virtual void readMsg(std::string& s) throw(IOException) = 0;

  //!the function to recieve integers
  virtual int32 readInt() throw(IOException) = 0;
  
  virtual void readIntVector(std::vector<int32>& vec) throw(IOException) = 0;

  virtual uint32 readUInt() throw(IOException) = 0;
  
  virtual void readUIntVector(std::vector<uint32>& vec) throw(IOException) = 0;

  virtual float32 readFloat() throw(IOException) = 0;
  
  virtual void readFloatVector(std::vector<float32>& vec) throw(IOException) = 0;

  virtual float64 readDouble() throw(IOException) = 0;

  virtual void readDoubleVector(std::vector<float64>& vec) throw(IOException) = 0;

  //!the function to recieve the vector;
  //  virtual data readVector() = 0;

  //!the function to Send integers
  virtual void  writeInt(const int32 buff) throw(IOException) = 0;

  virtual void writeIntArray(const int32 *array, const int32 size) throw(IOException) = 0;

  virtual void writeIntVector(const std::vector<int32>& vec) throw(IOException) = 0;

  virtual void  writeUInt(const uint32 buff) throw(IOException) = 0;

  virtual void writeUIntArray(const uint32 *array, const int32 size) throw(IOException) = 0;

  virtual void writeUIntVector(const std::vector<uint32>& vec) throw(IOException) = 0;


  virtual void  writeFloat(const float32 buff) throw(IOException) = 0;

  virtual void writeFloatArray(const float32 *array, const int32 size) throw(IOException) = 0;

  virtual void writeFloatVector(const std::vector<float32>& vec) throw(IOException) = 0;

  virtual void  writeDouble(const float64 buff) throw(IOException) = 0;
  
  virtual void writeDoubleArray(const float64 *array, const int32 size) throw(IOException) = 0;
  
  virtual void writeDoubleVector(const std::vector<float64>& vec) throw(IOException) = 0;

  virtual int32 ioOK() {return io_ok;};
  //!the function to Send the vector;
  //  virtual void writeVector(data outgoing_xdr_struct) = 0;
  virtual void Write(const void* data, int32 size) throw(IOException) = 0;
  virtual void Read(void* data, int32 size) throw(IOException) = 0;

 protected:

  FILE* file_read;
  FILE* file_write;

  int32 desc_read;
  int32 desc_write;

  int32 timeout_;
  int32 io_ok;

  Socket* sock_;

  //  virtual bool mayRead() throw(IOException) = 0;
  //  virtual bool mayWrite() throw(IOException) = 0;

  void is_valid() throw(IOException)
  {
    if(sock_ == NULL)
    {
      if((file_read == NULL) && (file_write == NULL))
        throw IOException("Read/Write Streams not valid!");
    }
    else
    {
      if(!sock_->is_valid())
        throw IOException("Socket not valid!");
    }

  }

};

BaseIO& operator << (BaseIO& io, const std::string& s ) throw(IOException);
BaseIO& operator >> (BaseIO& io, std::string& s) throw(IOException);

BaseIO& operator << (BaseIO& io, const int32 &i ) throw(IOException);
BaseIO& operator >> (BaseIO& io, int32& i) throw(IOException);

BaseIO& operator << (BaseIO& io, const float32 &f ) throw(IOException);
BaseIO& operator >> (BaseIO& io, float32& f) throw(IOException);

BaseIO& operator << (BaseIO& io, const float64 &d ) throw(IOException);
BaseIO& operator >> (BaseIO& io, float64& d) throw(IOException);

BaseIO& operator << (BaseIO& io, const std::vector<int32> &iv ) throw(IOException);
BaseIO& operator >> (BaseIO& io, std::vector<int32>& iv) throw(IOException);

BaseIO& operator << (BaseIO& io, const std::vector<uint32> &uiv ) throw(IOException);
BaseIO& operator >> (BaseIO& io, std::vector<uint32>& uiv) throw(IOException);

BaseIO& operator << (BaseIO& io, const std::vector<float32> &fv ) throw(IOException);
BaseIO& operator >> (BaseIO& io, std::vector<float32>& fv) throw(IOException);

BaseIO& operator << (BaseIO& io, const std::vector<float64> &dv ) throw(IOException);
BaseIO& operator >> (BaseIO& io, std::vector<float64>& dv) throw(IOException);

}


#endif //GSI_BASEIO
