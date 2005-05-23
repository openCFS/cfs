/*----------------------------------------------------------------------
  |
  | $Id$
  |
  +---------------------------------------------------------------------*/

#ifndef GSI_XDRIO
#define GSI_XDRIO

#include "GSITypeDefs.hh"
#include "GSIBaseIO.hh"
#include "GSIRawIO.hh"
#include "GSIIOException.hh"

namespace GridlibSocketInterface
{

  class XDRIO : public RawIO
  {
  public:
    XDRIO(FILE *file_read, FILE *file_write, bool bigendian = true);
#ifndef WIN32
    XDRIO(Socket *sock, int32 timeout = -1, bool bigendian = true);
#endif

    virtual void ReadMsg(std::string& s) throw(IOException);

    virtual void ReadShortVector(std::vector<int16>& vec) throw(IOException);
    virtual void ReadUShortVector(std::vector<uint16>& vec) throw(IOException);

    virtual void ReadIntVector(std::vector<int32>& vec) throw(IOException);
    virtual void ReadUIntVector(std::vector<uint32>& vec) throw(IOException);

    virtual void ReadFloatVector(std::vector<real32>& vec) throw(IOException);

    virtual void ReadDoubleVector(std::vector<real64>& vec) throw(IOException);



    virtual void WriteMsg(const std::string& s) throw(IOException);

    virtual void WriteShortArray(const int16 *array, const int32 size) throw(IOException);
    virtual void WriteShortVector(const std::vector<int16>& vec) throw(IOException);
    virtual void WriteUShortArray(const uint16 *array, const int32 size) throw(IOException);
    virtual void WriteUShortVector(const std::vector<uint16>& vec) throw(IOException);

    virtual void WriteIntArray(const int32 *array, const int32 size) throw(IOException);
    virtual void WriteIntVector(const std::vector<int32>& vec) throw(IOException);
    virtual void WriteUIntArray(const uint32 *array, const int32 size) throw(IOException);
    virtual void WriteUIntVector(const std::vector<uint32>& vec) throw(IOException);

    virtual void WriteFloatArray(const real32 *array, const int32 size) throw(IOException);
    virtual void WriteFloatVector(const std::vector<real32>& vec) throw(IOException);
    virtual void WriteDoubleArray(const real64 *array, const int32 size) throw(IOException);
    virtual void WriteDoubleVector(const std::vector<real64>& vec) throw(IOException);

  private:
    /*
      int32 genericRead(void* buf, int32 size);
      int32 genericWrite(void* buf, int32 size);
      int32 genericReadInt(int32& val);
      int32 genericWriteInt(int32 val);
    */
  };

}

#endif //GSI_XDRIO
