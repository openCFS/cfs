// Definition of the GSIRawIO class

#ifndef GSI_RAWIO
#define GSI_RAWIO

#ifdef __sgi
#include <stdio.h>
#else
#include <cstdio>
#endif

#include "GSIBaseIO.hh"
#include "GSIIOException.hh"
#include "GSITypeDefs.hh"

namespace GridlibSocketInterface
{

class RawIO : public BaseIO
{
 public:
  RawIO(FILE *file_read, FILE *file_write, bool bigendian = true);
  RawIO(Socket *sock, int32 timeout = -1, bool bigendian = true);

  enum SeekPos {SET, CUR, END};

  // tell the current stream position
  virtual int32 TellRead() throw(IOException);
  virtual int32 TellWrite() throw(IOException);

  // seek the position
  virtual int32 SeekRead(int32 offset, SeekPos whence) throw(IOException);
  virtual int32 SeekWrite(int32 offset, SeekPos whence) throw(IOException);

  virtual int32 EOFRead() throw(IOException);
//  int32 EOFWrite();

  // the Send Message function
  virtual void writeMsg(const std::string& s) throw(IOException);

  //!the function to recieve the Message; Returns a pointer to the string
  virtual void readMsg(std::string& s) throw(IOException);

  //!the function to recieve integers
  virtual int32 readInt() throw(IOException);

  virtual void readIntVector(std::vector<int32>& vec) throw(IOException);

  virtual uint32 readUInt() throw(IOException);

  virtual void readUIntVector(std::vector<uint32>& vec) throw(IOException);

  virtual float32 readFloat() throw(IOException);

  virtual void readFloatVector(std::vector<float32>& vec) throw(IOException);

  virtual float64 readDouble() throw(IOException);

  virtual void readDoubleVector(std::vector<float64>& vec) throw(IOException);

  //!the function to recieve the vector;
  //  virtual data readVector();

  //!the function to Send integers
  virtual void  writeInt(const int32 buff) throw(IOException);

  virtual void writeIntArray(const int32 *array, const int32 size) throw(IOException);

  virtual void writeIntVector(const std::vector<int32>& vec) throw(IOException);

  virtual void  writeUInt(const uint32 buff) throw(IOException);

  virtual void writeUIntArray(const uint32 *array, const int size) throw(IOException);

  virtual void writeUIntVector(const std::vector<uint32>& vec) throw(IOException);

  virtual void  writeFloat(const float32 buff) throw(IOException);

  virtual void writeFloatArray(const float32 *array, const int32 size) throw(IOException);

  virtual void writeFloatVector(const std::vector<float32>& vec) throw(IOException);

  virtual void  writeDouble(const float64 buff) throw(IOException);

  virtual void writeDoubleArray(const float64 *array, const int32 size) throw(IOException);

  virtual void writeDoubleVector(const std::vector<float64>& vec) throw(IOException);
  //!the function to Send the vector;
  //  virtual void writeVector(data outgoing_xdr_struct);

  virtual void Write(const void* data, int32 size) throw(IOException);
  virtual void Read(void* data, int32 size) throw(IOException);

  void SetBigEndianOn() {bigendian_ = true;}
  void SetLittleEndianOn() {bigendian_ = false;}
  protected:

  bool bigendian_;
};

}

#endif //GSI_RAWIO
