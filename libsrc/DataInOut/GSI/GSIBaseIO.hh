/*----------------------------------------------------------------------
  |
  | $Id$
  |
  +---------------------------------------------------------------------*/

#ifndef GSI_BASEIO
#define GSI_BASEIO

#include <sys/types.h>
#ifndef WIN32
#include <sys/time.h>
#endif
#include <cstdio>
#include <stdio.h>
#include <string>
#include <vector>

#include "GSITypeDefs.hh"
#include "GSIIOException.hh"
#ifndef WIN32
#include "GSISocket.hh"
#endif

namespace GridlibSocketInterface
{

  class BaseIO
  {
  public:
    
    //! Constructor for file IO
    BaseIO(FILE* fr, FILE* fw) : 
      file_read( fr ), file_write( fw ), timeout_(-1), sock_(NULL)  {};

#ifndef WIN32
    //! Constructor for socket IO
    BaseIO(Socket *sock, int32 timeout = -1) : 
      file_read( NULL ), file_write( NULL ), timeout_(timeout),sock_(sock) {};
#endif

    virtual ~BaseIO() {};

    enum SeekPos {SET, CUR, END};

    // tell the current stream position
    virtual int32 TellRead() throw(IOException) = 0;
    virtual int32 TellWrite() throw(IOException) = 0;
    
    // seek the position
    virtual int32 SeekRead(int32 offset, SeekPos whence) throw(IOException) = 0;
    virtual int32 SeekWrite(int32 offset, SeekPos whence) throw(IOException) = 0;
    
    virtual int32 EOFRead() throw(IOException) = 0;
    
    virtual void   ReadMsg(std::string& s) throw(IOException) = 0;
    
    virtual int16  ReadShort() throw(IOException) = 0;
    virtual void   ReadShortVector(std::vector<int16>& vec) throw(IOException) = 0;
    virtual uint16 ReadUShort() throw(IOException) = 0;
    virtual void   ReadUShortVector(std::vector<uint16>& vec) throw(IOException) = 0;
    
    virtual int32  ReadInt() throw(IOException) = 0;
    virtual void   ReadIntVector(std::vector<int32>& vec) throw(IOException) = 0;
    virtual uint32 ReadUInt() throw(IOException) = 0;
    virtual void   ReadUIntVector(std::vector<uint32>& vec) throw(IOException) = 0;
    
    virtual real32 ReadFloat() throw(IOException) = 0;
    virtual void   ReadFloatVector(std::vector<real32>& vec) throw(IOException) = 0;
    
    virtual real64 ReadDouble() throw(IOException) = 0;
    virtual void   ReadDoubleVector(std::vector<real64>& vec) throw(IOException) = 0;
    

    virtual void WriteMsg(const std::string& s) throw(IOException) = 0;
    
    virtual void WriteShort(const int16 buff) throw(IOException) = 0;
    virtual void WriteShortArray(const int16 *array, const int32 size) throw(IOException) = 0;
    virtual void WriteShortVector(const std::vector<int16>& vec) throw(IOException) = 0;
    virtual void WriteUShort(const uint16 buff) throw(IOException) = 0;
    virtual void WriteUShortArray(const uint16 *array, const int32 size) throw(IOException) = 0;
    virtual void WriteUShortVector(const std::vector<uint16>& vec) throw(IOException) = 0;

    virtual void WriteInt(const int32 buff) throw(IOException) = 0;
    virtual void WriteIntArray(const int32 *array, const int32 size) throw(IOException) = 0;
    virtual void WriteIntVector(const std::vector<int32>& vec) throw(IOException) = 0;
    virtual void WriteUInt(const uint32 buff) throw(IOException) = 0;
    virtual void WriteUIntArray(const uint32 *array, const int32 size) throw(IOException) = 0;
    virtual void WriteUIntVector(const std::vector<uint32>& vec) throw(IOException) = 0;
    
    virtual void WriteFloat(const real32 buff) throw(IOException) = 0;
    virtual void WriteFloatArray(const real32 *array, const int32 size) throw(IOException) = 0;
    virtual void WriteFloatVector(const std::vector<real32>& vec) throw(IOException) = 0;
    
    virtual void WriteDouble(const real64 buff) throw(IOException) = 0;
    virtual void WriteDoubleArray(const real64 *array, const int32 size) throw(IOException) = 0;
    virtual void WriteDoubleVector(const std::vector<real64>& vec) throw(IOException) = 0;
    
    virtual int32 ioOK() {return io_ok;};

    virtual void Write(const void* data, int32 size) throw(IOException) = 0;
    virtual void Read(void* data, int32 size) throw(IOException) = 0;
    
  protected:
    
    FILE* file_read;
    FILE* file_write;
    
    int32 desc_read;
    int32 desc_write;
    
    int32 timeout_;
    int32 io_ok;
#ifndef WIN32
    Socket* sock_;
#else
    int *sock_;
#endif
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
#ifndef WIN32
          if(!sock_->is_valid())
            throw IOException("Socket not valid!");
#else
          return;
#endif
        }
      
    }
  };

  BaseIO& operator << (BaseIO& io, const std::string& s ) throw(IOException);
  BaseIO& operator >> (BaseIO& io, std::string& s) throw(IOException);
  
  BaseIO& operator << (BaseIO& io, const int16 &i ) throw(IOException);
  BaseIO& operator >> (BaseIO& io, int16& i) throw(IOException);
  BaseIO& operator << (BaseIO& io, const uint16 &i ) throw(IOException);
  BaseIO& operator >> (BaseIO& io, uint16& i) throw(IOException);
  
  BaseIO& operator << (BaseIO& io, const int32 &i ) throw(IOException);
  BaseIO& operator >> (BaseIO& io, int32& i) throw(IOException);
  BaseIO& operator << (BaseIO& io, const uint32 &i ) throw(IOException);
  BaseIO& operator >> (BaseIO& io, uint32& i) throw(IOException);

  BaseIO& operator << (BaseIO& io, const real32 &f ) throw(IOException);
  BaseIO& operator >> (BaseIO& io, real32& f) throw(IOException);
  
  BaseIO& operator << (BaseIO& io, const real64 &d ) throw(IOException);
  BaseIO& operator >> (BaseIO& io, real64& d) throw(IOException);
  
  BaseIO& operator << (BaseIO& io, const std::vector<int16> &iv ) throw(IOException);
  BaseIO& operator >> (BaseIO& io, std::vector<int16>& iv) throw(IOException);
  BaseIO& operator << (BaseIO& io, const std::vector<uint16> &uiv ) throw(IOException);
  BaseIO& operator >> (BaseIO& io, std::vector<uint16>& uiv) throw(IOException);
  
  BaseIO& operator << (BaseIO& io, const std::vector<int32> &iv ) throw(IOException);
  BaseIO& operator >> (BaseIO& io, std::vector<int32>& iv) throw(IOException);
  BaseIO& operator << (BaseIO& io, const std::vector<uint32> &uiv ) throw(IOException);
  BaseIO& operator >> (BaseIO& io, std::vector<uint32>& uiv) throw(IOException);
  
  
  BaseIO& operator << (BaseIO& io, const std::vector<real32> &fv ) throw(IOException);
  BaseIO& operator >> (BaseIO& io, std::vector<real32>& fv) throw(IOException);
  
  BaseIO& operator << (BaseIO& io, const std::vector<real64> &dv ) throw(IOException);
  BaseIO& operator >> (BaseIO& io, std::vector<real64>& dv) throw(IOException);

}


#endif //GSI_BASEIO
