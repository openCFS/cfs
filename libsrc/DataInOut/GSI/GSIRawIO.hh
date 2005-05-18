/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/

#ifndef GSI_RAWIO
#define GSI_RAWIO

#include <cstdio>

#include "GSIBaseIO.hh"
#include "GSIIOException.hh"
#include "GSITypeDefs.hh"

namespace GridlibSocketInterface
{

  class RawIO : public BaseIO
  {
  public:
    RawIO(FILE *file_read, FILE *file_write, bool bigendian = true);
#ifndef WIN32
    RawIO(Socket *sock, int32 timeout = -1, bool bigendian = true);
#endif

    // tell the current stream position
    virtual int32 TellRead() throw(IOException);
    virtual int32 TellWrite() throw(IOException);
    
    // seek the position
    virtual int32 SeekRead(int32 offset, SeekPos whence) throw(IOException);
    virtual int32 SeekWrite(int32 offset, SeekPos whence) throw(IOException);
    
    virtual int32 EOFRead() throw(IOException);
    //  int32 EOFWrite();
    
    virtual void   ReadMsg(std::string& s) throw(IOException);
    
    virtual int16  ReadShort() throw(IOException);
    virtual void   ReadShortVector(std::vector<int16>& vec) throw(IOException);
    virtual uint16 ReadUShort() throw(IOException);
    virtual void   ReadUShortVector(std::vector<uint16>& vec) throw(IOException);
    
    virtual int32  ReadInt() throw(IOException);
    virtual void   ReadIntVector(std::vector<int32>& vec) throw(IOException);
    virtual uint32 ReadUInt() throw(IOException);
    virtual void   ReadUIntVector(std::vector<uint32>& vec) throw(IOException);
    
    virtual real32 ReadFloat() throw(IOException);
    virtual void   ReadFloatVector(std::vector<real32>& vec) throw(IOException);
    
    virtual real64 ReadDouble() throw(IOException);
    virtual void   ReadDoubleVector(std::vector<real64>& vec) throw(IOException);
    

    virtual void WriteMsg(const std::string& s) throw(IOException);
    
    virtual void WriteShort(const int16 buff) throw(IOException);
    virtual void WriteShortArray(const int16 *array, const int32 size) throw(IOException);
    virtual void WriteShortVector(const std::vector<int16>& vec) throw(IOException);
    virtual void WriteUShort(const uint16 buff) throw(IOException);
    virtual void WriteUShortArray(const uint16 *array, const int32 size) throw(IOException);
    virtual void WriteUShortVector(const std::vector<uint16>& vec) throw(IOException);

    virtual void WriteInt(const int32 buff) throw(IOException);
    virtual void WriteIntArray(const int32 *array, const int32 size) throw(IOException);
    virtual void WriteIntVector(const std::vector<int32>& vec) throw(IOException);
    virtual void WriteUInt(const uint32 buff) throw(IOException);
    virtual void WriteUIntArray(const uint32 *array, const int32 size) throw(IOException);
    virtual void WriteUIntVector(const std::vector<uint32>& vec) throw(IOException);
    
    virtual void WriteFloat(const real32 buff) throw(IOException);
    virtual void WriteFloatArray(const real32 *array, const int32 size) throw(IOException);
    virtual void WriteFloatVector(const std::vector<real32>& vec) throw(IOException);
    
    virtual void WriteDouble(const real64 buff) throw(IOException);
    virtual void WriteDoubleArray(const real64 *array, const int32 size) throw(IOException);
    virtual void WriteDoubleVector(const std::vector<real64>& vec) throw(IOException);

    virtual void Write(const void* data, int32 size) throw(IOException);
    virtual void Read(void* data, int32 size) throw(IOException);
    
    void SetBigEndianOn() {bigendian_ = true;}
    void SetLittleEndianOn() {bigendian_ = false;}
  protected:
    
    bool bigendian_;
  };

}

#endif //GSI_RAWIO
