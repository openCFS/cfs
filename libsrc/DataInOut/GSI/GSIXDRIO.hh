// Definition of the GSI::XDRIO class

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
  XDRIO(Socket *sock, int32 timeout = -1, bool bigendian = true);

  // the Send Message function
  virtual void writeMsg(const std::string& s) throw(IOException);

  //!the function to recieve the Message; Returns a pointer to the string
  virtual void readMsg(std::string& s) throw(IOException);

  virtual void readIntVector(std::vector<int32>& vec) throw(IOException);
  virtual void readUIntVector(std::vector<uint32>& vec) throw(IOException);
  virtual void readFloatVector(std::vector<float32>& vec) throw(IOException);
  virtual void readDoubleVector(std::vector<float64>& vec) throw(IOException);

  virtual void writeIntArray(const int32 *array, const int32 size) throw(IOException);
  virtual void writeIntVector(const std::vector<int32>& vec) throw(IOException);
  virtual void writeUIntArray(const uint32 *array, const int32 size) throw(IOException);
  virtual void writeUIntVector(const std::vector<uint32>& vec) throw(IOException);
  virtual void writeFloatArray(const float32 *array, const int32 size) throw(IOException);
  virtual void writeFloatVector(const std::vector<float32>& vec) throw(IOException);
  virtual void writeDoubleArray(const float64 *array, const int32 size) throw(IOException);
  virtual void writeDoubleVector(const std::vector<float64>& vec) throw(IOException);

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
