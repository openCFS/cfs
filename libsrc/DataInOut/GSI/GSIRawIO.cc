#ifdef USE_RCSID
static const char RCSid_GSIRawIO[] = "$Id$";
#endif

/*----------------------------------------------------------------------
  |
  |
  | $Log$
  | Revision 1.4  2005/05/23 22:15:31  ahauck
  | Corrected indentation of all cc/hh-files according to xemacs' 'gnu' style
  | and according to our coding rules, so all TABs were replaced by spaces.
  |
  | Revision 1.3  2005/05/18 19:26:03  strieben
  | Upgraded GSI library to newest available version.
  |
  | Revision 1.2  2004/09/01 15:24:58  simon
  | Added support for writing int16 and uint16
  |
  | Revision 1.1.1.1  2004/08/31 15:53:00  simon
  | Initial GSI import
  |
  |
  +---------------------------------------------------------------------*/

#include <cstdlib>
#include <cstring>
#include <iostream>

#include "GSIRawIO.hh"
#include "GSIEndianConverter.hh"

namespace GridlibSocketInterface
{

#define CHECK_READ() if(file_read == NULL) throw IOException("Attempt to access NULL file_read")
#define CHECK_WRITE() if(file_write == NULL) throw IOException("Attempt to access NULL file_write")

  RawIO::RawIO(FILE *file_read, FILE *file_write, bool bigendian) :
    BaseIO(file_read, file_write),
    bigendian_(bigendian)
  {
    InitEndian();
  }
  
#ifndef WIN32
  RawIO::RawIO(Socket *sock, int32 timeout, bool bigendian) :
    BaseIO(sock, timeout),
    bigendian_(bigendian)
  {
    InitEndian();
  }
#endif

  int32 RawIO :: TellRead() throw(IOException)
  {
    CHECK_READ();
    int32 ret;

    ret = ftell(file_read);
    if(ret == -1)
      {
        IOException ex("Error when calling ftell(file_read)");
        ex.SetErrno(errno);
        throw ex;
      }
    return ret;
  }

  int32 RawIO :: TellWrite() throw(IOException)
  {
    CHECK_WRITE();

    int32 ret;
    ret = ftell(file_write);
    if(ret == -1)
      {
        IOException ex("Error when calling ftell(file_write)");
        ex.SetErrno(errno);
        throw ex;
      }
    return ret;
  }

  int32 RawIO :: SeekRead(int32 offset, SeekPos whence) throw(IOException)
  {
    CHECK_READ();
    int32 w;
    int32 ret;

    switch(whence)
      {
      case SET:
        w = SEEK_SET;
        break;
      case CUR:
        w = SEEK_CUR;
        break;
      case END:
        w = SEEK_END;
        break;
      default:
        w = SEEK_SET;
        break;
      }

    ret = fseek(file_read, offset, w);
    if(ret == -1)
      {
        IOException ex("Error when calling fseek(file_read)");
        ex.SetErrno(errno);
        throw ex;
      }
    return ret;
  }

  int32 RawIO :: SeekWrite(int32 offset, SeekPos whence) throw(IOException)
  {
    CHECK_WRITE();
    int32 w;
    int32 ret;

    switch(whence)
      {
      case SET:
        w = SEEK_SET;
        break;
      case CUR:
        w = SEEK_CUR;
        break;
      case END:
        w = SEEK_END;
        break;
      default:
        w = SEEK_SET;
        break;
      }

    ret = fseek(file_write, offset, w);
    if(ret == -1)
      {
        IOException ex("Error when calling fseek(file_write)");
        ex.SetErrno(errno);
        throw ex;
      }
    return ret;
  }

  int32 RawIO :: EOFRead() throw(IOException)
  {
    CHECK_READ();

    int32 ret;
    ret = feof(file_read);
    if(ret == -1)
      {
        IOException ex("Error when calling feof(file_read)");
        ex.SetErrno(errno);
        throw ex;
      }
    return ret;
  }

  void RawIO :: ReadMsg(std::string& s) throw(IOException)
  {
    CHECK_READ();
    int32 length = ReadInt();
    int32 pad = 4 - (length % 4);
    int32 bufsize = (pad == 4) ? length : length + pad;

    char* msg= new char[bufsize+1];
    memset(msg, 0, bufsize+1);

    try{
      Read(msg, bufsize);
    }
    catch(IOException& ex) {
      delete[] msg;
      ex.SetDescription("An Error occured while reading the actual message");
      throw ex;
    }

    msg[length] = 0;
    s = msg;
    delete[] msg;
  }

  int16 RawIO :: ReadShort() throw(IOException)
  {
    CHECK_READ();
    // I use int32 because I want to make sure that that the files
    // have always a size that is a multiple of 4 bytes
    int32 buff;  

    try{
      Read(&buff, sizeof(int32));
    }
    catch(IOException& ex) {
      ex.SetDescription("An Error occured while reading a Short");
      throw ex;
    }

    buff = bigendian_ ? BigINT32(buff) : LittleINT32(buff);

    return((int16) buff);
  }


  void RawIO :: ReadShortVector(std::vector<int16>& vec) throw(IOException)
  {
    CHECK_READ();
    int32 size = ReadInt();
    int16 *buff= new int16[size];
    int16 dummy;

    try{
      Read(buff, sizeof(int16)*size);
      // make sure the 4 byte padding is right
      if((size % 2) == 1)
        Read(&dummy, sizeof(int16));
    }
    catch(IOException& ex) {
      delete[] buff;
      ex.SetDescription("An Error occured while reading a short array");
      throw ex;
    }

    vec.resize(size);
    for(int32 i=0; i<size; i++) {
      vec[i] = bigendian_ ? BigINT16(buff[i]) : LittleINT16(buff[i]);
    }

    delete[] buff;
  }

  uint16 RawIO :: ReadUShort() throw(IOException)
  {
    CHECK_READ();
    uint32 buff;

    try{
      Read(&buff, sizeof(uint32));
    }
    catch(IOException& ex) {
      ex.SetDescription("An Error occured while reading an unsigned short");
      throw ex;
    }

    buff = bigendian_ ? BigUINT32(buff) : LittleUINT32(buff);

    return((uint16) buff);
  }


  void RawIO :: ReadUShortVector(std::vector<uint16>& vec) throw(IOException)
  {
    CHECK_READ();
    int32 size = ReadInt();
    uint16 *buff= new uint16[size];
    uint16 dummy;

    try{
      Read(buff, sizeof(uint16)*size);
      // make sure the 4 byte padding is right
      if((size % 2) == 1)
        Read(&dummy, sizeof(uint16));
    }
    catch(IOException& ex) {
      delete[] buff;
      ex.SetDescription("An Error occured while reading an unsigned integer vector");
      throw ex;
    }

    vec.resize(size);
    for(int32 i=0; i<size; i++) {
      vec[i] = bigendian_ ? BigINT16(buff[i]) : LittleINT16(buff[i]);
    }

    delete[] buff;
  }

  int32 RawIO :: ReadInt() throw(IOException)
  {
    CHECK_READ();
    int32 buff;

    try{
      Read(&buff, sizeof(int32));
    }
    catch(IOException& ex) {
      ex.SetDescription("An Error occured while reading an Integer");
      throw ex;
    }

    buff = bigendian_ ? BigINT32(buff) : LittleINT32(buff);

    return(buff);
  }


  void RawIO :: ReadIntVector(std::vector<int32>& vec) throw(IOException)
  {
    CHECK_READ();
    int32 size = ReadInt();
    int32 *buff= new int32[size];

    try{
      Read(buff, sizeof(int32)*size);
    }
    catch(IOException& ex) {
      delete[] buff;
      ex.SetDescription("An Error occured while reading an integer vector");
      throw ex;
    }

    vec.resize(size);
    for(int32 i=0; i<size; i++) {
      vec[i] = bigendian_ ? BigINT32(buff[i]) : LittleINT32(buff[i]);
    }

    delete[] buff;
  }

  uint32 RawIO :: ReadUInt() throw(IOException)
  {
    CHECK_READ();
    uint32 buff;

    try{
      Read(&buff, sizeof(uint32));
    }
    catch(IOException& ex) {
      ex.SetDescription("An Error occured while reading an Integer");
      throw ex;
    }

    buff = bigendian_ ? BigUINT32(buff) : LittleUINT32(buff);

    return(buff);
  }


  void RawIO :: ReadUIntVector(std::vector<uint32>& vec) throw(IOException)
  {
    CHECK_READ();
    int32 size = ReadInt();
    uint32 *buff= new uint32[size];

    try{
      Read(buff, sizeof(uint32)*size);
    }
    catch(IOException& ex) {
      delete[] buff;
      ex.SetDescription("An Error occured while reading an integer vector");
      throw ex;
    }

    vec.resize(size);
    for(int32 i=0; i<size; i++) {
      vec[i] = bigendian_ ? BigINT32(buff[i]) : LittleINT32(buff[i]);
    }

    delete[] buff;
  }

  real32 RawIO :: ReadFloat() throw(IOException)
  {
    CHECK_READ();
    real32 buff;

    try{
      Read(&buff, sizeof(real32));
    }
    catch(IOException& ex) {
      ex.SetDescription("An Error occured while reading a Float");
      throw ex;
    }

    buff = bigendian_ ? BigREAL32(buff) : LittleREAL32(buff);

    return(buff);
  }


  void RawIO :: ReadFloatVector(std::vector<real32>& vec) throw(IOException)
  {
    CHECK_READ();
    int32 size = ReadInt();
    real32 *buff= new real32[size];

    try{
      Read(buff, sizeof(real32)*size);
    }
    catch(IOException& ex) {
      delete[] buff;
      ex.SetDescription("An Error occured while reading a float array");
      throw ex;
    }

    vec.resize(size);
    for(int32 i=0; i<size; i++) {
      vec[i] = bigendian_ ? BigREAL32(buff[i]) : LittleREAL32(buff[i]);
    }

    delete[] buff;
  }

  real64 RawIO :: ReadDouble() throw(IOException)
  {
    CHECK_READ();
    real64 buff;

    try{
      Read(&buff, sizeof(real64));
    }
    catch(IOException& ex) {
      ex.SetDescription("An Error occured while reading a Double");
      throw ex;
    }

    buff = bigendian_ ? BigREAL64(buff) : LittleREAL64(buff);

    return(buff);
  }


  void RawIO :: ReadDoubleVector(std::vector<real64>& vec) throw(IOException)
  {
    CHECK_READ();
    int32 size = ReadInt();
    real64 *buff= new real64[size];

    try{
      Read(buff, sizeof(real64)*size);
    }
    catch(IOException& ex) {
      delete[] buff;
      ex.SetDescription("An Error occured while reading a double array");
      throw ex;
    }

    vec.resize(size);
    for(int32 i=0; i<size; i++) {
      vec[i] = bigendian_ ? BigREAL64(buff[i]) : LittleREAL64(buff[i]);
    }

    delete[] buff;
  }

  void RawIO :: WriteMsg(const std::string& s) throw(IOException)
  {
    CHECK_WRITE();
    int32 length = s.length();
    int32 l = bigendian_ ? BigINT32(length) : LittleINT32(length);
    int32 pad = 4 - (length % 4);
    int32 bufsize = (pad == 4) ? length : length + pad;

    try{
      Write(&l, sizeof(int32));
    }
    catch(IOException& ex) {
      ex.SetDescription("An Error occured while writing the message length");
      throw ex;
    }

    char* msg= new char[bufsize];
    memset(msg, 0, bufsize);
    memcpy(msg, s.c_str(), length);

    /*
      std::cout << "length: " << length << "\n";
      std::cout << "pad: " << pad << "\n";
      std::cout << "buf_size: " << bufsize << "\n";
    */

    try{
      Write(msg, bufsize);
    }
    catch(IOException& ex) {
      delete[] msg;
      ex.SetDescription("An Error occured while writing the actual message");
      throw ex;
    }
    delete[] msg;
  }

  void RawIO :: WriteShort(const int16 buff) throw(IOException)
  {
    CHECK_WRITE();
    int32 out = bigendian_ ? BigINT32((int32)buff) : LittleINT32((int32)buff);
    try{
      Write(&out, sizeof(int32));
    }
    catch(IOException& ex) {
      ex.SetDescription("An Error occured while writing a short");
      throw ex;
    }
  }


  void RawIO :: WriteShortArray(const int16 *buff , const int32 size) throw(IOException)
  {
    CHECK_WRITE();
    int32 s = bigendian_ ? BigINT32(size) : LittleINT32(size);
    int16 dummy = 0;

    try{
      Write(&s, sizeof(int32));
    }
    catch(IOException& ex) {
      ex.SetDescription("Error while writing the length of short array");
      throw ex;
    }

    int16 *array = new int16[size];
    for(int32 i=0; i<size; i++)
      array[i] = bigendian_ ? BigINT16(buff[i]) : LittleINT16(buff[i]);

    try{
      Write(array, sizeof(int16)*size);
      if((size % 2) == 1)
        Write(&dummy, sizeof(int16));
    }
    catch(IOException& ex) {
      delete[] array;
      ex.SetDescription("An Error occured while writing an integer array");
      throw ex;
    }

    delete[] array;
  }

  void RawIO :: WriteShortVector(const std::vector<int16>& vec) throw(IOException)
  {
    CHECK_WRITE();
    int32 size = vec.size();
    int32 s = bigendian_ ? BigINT32(size) : LittleINT32(size);
    int16 dummy = 0;

    try{
      Write(&s, sizeof(int32));
    }
    catch(IOException& ex) {
      ex.SetDescription("Error while writing the length of short vector");
      throw ex;
    }

    int16 *array = new int16[size];
    for(int32 i=0; i<size; i++)
      array[i] = bigendian_ ? BigINT16(vec[i]) : LittleINT16(vec[i]);

    try{
      Write(array, sizeof(int16)*size);
      if((size % 2) == 1)
        Write(&dummy, sizeof(int16));
    }
    catch(IOException& ex) {
      delete[] array;
      ex.SetDescription("An Error occured while writing a short vector");
      throw ex;
    }

    delete[] array;
  }

  void RawIO :: WriteUShort(const uint16 buff) throw(IOException)
  {
    CHECK_WRITE();
    uint32 out = bigendian_ ? BigUINT32((int32)buff) : LittleUINT16((int32)buff);

    try{
      Write(&out, sizeof(uint32));
    }
    catch(IOException& ex) {
      ex.SetDescription("An Error occured while writing an unsigned short");
      throw ex;
    }
  }


  void RawIO :: WriteUShortArray(const uint16 *buff , const int32 size) throw(IOException)
  {
    CHECK_WRITE();
    int32 s = bigendian_ ? BigINT32(size) : LittleINT32(size);
    uint16 dummy = 0;

    try{
      Write(&s, sizeof(int32));
    }
    catch(IOException& ex) {
      ex.SetDescription("Error while writing the length of unsigned short array");
      throw ex;
    }

    uint16 *array = new uint16[size];
    for(int32 i=0; i<size; i++)
      array[i] = bigendian_ ? BigUINT16(buff[i]) : LittleUINT16(buff[i]);

    try{
      Write(array, sizeof(uint16)*size);
      if((size % 2) == 1)
        Write(&dummy, sizeof(uint16));
    }
    catch(IOException& ex) {
      delete[] array;
      ex.SetDescription("An Error occured while writing an unsigned short array");
      throw ex;
    }

    delete[] array;
  }

  void RawIO :: WriteUShortVector(const std::vector<uint16>& vec) throw(IOException)
  {
    CHECK_WRITE();
    int32 size = vec.size();
    int32 s = bigendian_ ? BigINT32(size) : LittleINT32(size);
    uint16 dummy = 0;

    try{
      Write(&s, sizeof(int32));
    }
    catch(IOException& ex) {
      ex.SetDescription("Error while writing the length of unsigned short vector");
      throw ex;
    }

    uint16 *array = new uint16[size];
    for(int32 i=0; i<size; i++)
      array[i] = bigendian_ ? BigUINT16(vec[i]) : LittleUINT16(vec[i]);

    try{
      Write(array, sizeof(uint16)*size);
      if((size % 2) == 1)
        Write(&dummy, sizeof(uint16));
    }
    catch(IOException& ex) {
      delete[] array;
      ex.SetDescription("An Error occured while writing an unsigned short vector");
      throw ex;
    }

    delete[] array;
  }

  void RawIO :: WriteInt(const int32 buff) throw(IOException)
  {
    CHECK_WRITE();
    int32 out = bigendian_ ? BigINT32(buff) : LittleINT32(buff);
    try{
      Write(&out, sizeof(int32));
    }
    catch(IOException& ex) {
      ex.SetDescription("An Error occured while writing an integer");
      throw ex;
    }
  }


  void RawIO :: WriteIntArray(const int32 *buff , const int32 size) throw(IOException)
  {
    CHECK_WRITE();
    int32 s = bigendian_ ? BigINT32(size) : LittleINT32(size);

    try{
      Write(&s, sizeof(int32));
    }
    catch(IOException& ex) {
      ex.SetDescription("Error while writing the length of integer array");
      throw ex;
    }

    int32 *array = new int32[size];
    for(int32 i=0; i<size; i++)
      array[i] = bigendian_ ? BigINT32(buff[i]) : LittleINT32(buff[i]);

    try{
      Write(array, sizeof(int32)*size);
    }
    catch(IOException& ex) {
      delete[] array;
      ex.SetDescription("An Error occured while writing an integer array");
      throw ex;
    }

    delete[] array;
  }

  void RawIO :: WriteIntVector(const std::vector<int32>& vec) throw(IOException)
  {
    CHECK_WRITE();
    int32 size = vec.size();
    int32 s = bigendian_ ? BigINT32(size) : LittleINT32(size);

    try{
      Write(&s, sizeof(int32));
    }
    catch(IOException& ex) {
      ex.SetDescription("Error while writing the length of integer vector");
      throw ex;
    }

    int32 *array = new int32[size];
    for(int32 i=0; i<size; i++)
      array[i] = bigendian_ ? BigINT32(vec[i]) : LittleINT32(vec[i]);

    try{
      Write(array, sizeof(int32)*size);
    }
    catch(IOException& ex) {
      delete[] array;
      ex.SetDescription("An Error occured while writing an integer vector");
      throw ex;
    }

    delete[] array;
  }

  void RawIO :: WriteUInt(const uint32 buff) throw(IOException)
  {
    CHECK_WRITE();
    uint32 out = bigendian_ ? BigUINT32(buff) : LittleUINT32(buff);

    try{
      Write(&out, sizeof(uint32));
    }
    catch(IOException& ex) {
      ex.SetDescription("An Error occured while writing an uint");
      throw ex;
    }
  }


  void RawIO :: WriteUIntArray(const uint32 *buff , const int32 size) throw(IOException)
  {
    CHECK_WRITE();
    int32 s = bigendian_ ? BigINT32(size) : LittleINT32(size);

    try{
      Write(&s, sizeof(int32));
    }
    catch(IOException& ex) {
      ex.SetDescription("Error while writing the length of uint array");
      throw ex;
    }

    uint32 *array = new uint32[size];
    for(int32 i=0; i<size; i++)
      array[i] = bigendian_ ? BigUINT32(buff[i]) : LittleUINT32(buff[i]);

    try{
      Write(array, sizeof(uint32)*size);
    }
    catch(IOException& ex) {
      delete[] array;
      ex.SetDescription("An Error occured while writing an uint array");
      throw ex;
    }

    delete[] array;
  }

  void RawIO :: WriteUIntVector(const std::vector<uint32>& vec) throw(IOException)
  {
    CHECK_WRITE();
    int32 size = vec.size();
    int32 s = bigendian_ ? BigINT32(size) : LittleINT32(size);

    try{
      Write(&s, sizeof(int32));
    }
    catch(IOException& ex) {
      ex.SetDescription("Error while writing the length of uint vector");
      throw ex;
    }

    uint32 *array = new uint32[size];
    for(int32 i=0; i<size; i++)
      array[i] = bigendian_ ? BigUINT32(vec[i]) : LittleUINT32(vec[i]);

    try{
      Write(array, sizeof(uint32)*size);
    }
    catch(IOException& ex) {
      delete[] array;
      ex.SetDescription("An Error occured while writing an uint vector");
      throw ex;
    }

    delete[] array;
  }

  void RawIO :: WriteFloat(const real32 buff) throw(IOException)
  {
    CHECK_WRITE();
    real32 out = bigendian_ ? BigREAL32(buff) : LittleREAL32(buff);

    try{
      Write(&out, sizeof(real32));
    }
    catch(IOException& ex) {
      ex.SetDescription("An Error occured while writing a float");
      throw ex;
    }
  }


  void RawIO :: WriteFloatArray(const real32 *buff , const int32 size) throw(IOException)
  {
    CHECK_WRITE();
    int32 s = bigendian_ ? BigINT32(size) : LittleINT32(size);

    try{
      Write(&s, sizeof(int32));
    }
    catch(IOException& ex) {
      ex.SetDescription("Error while writing the length of float array");
      throw ex;
    }

    real32 *array = new real32[size];
    for(int32 i=0; i<size; i++)
      array[i] = bigendian_ ? BigREAL32(buff[i]) : LittleREAL32(buff[i]);

    try{
      Write(array, sizeof(real32)*size);
    }
    catch(IOException& ex) {
      delete[] array;
      ex.SetDescription("An Error occured while writing a float array");
      throw ex;
    }

    delete[] array;
  }

  void RawIO :: WriteFloatVector(const std::vector<real32>& vec) throw(IOException)
  {
    CHECK_WRITE();
    int32 size = vec.size();
    int32 s = bigendian_ ? BigINT32(size) : LittleINT32(size);

    try{
      Write(&s, sizeof(int32));
    }
    catch(IOException& ex) {
      ex.SetDescription("Error while writing the length of float vector");
      throw ex;
    }

    real32 *array = new real32[size];
    for(int32 i=0; i<size; i++)
      array[i] = bigendian_ ? BigREAL32(vec[i]) : LittleREAL32(vec[i]);

    try{
      Write(array, sizeof(real32)*size);
    }
    catch(IOException& ex) {
      delete[] array;
      ex.SetDescription("An Error occured while writing a float vector");
      throw ex;
    }

    delete[] array;
  }

  void RawIO :: WriteDouble(const real64 buff) throw(IOException)
  {
    CHECK_WRITE();
    real64 out = bigendian_ ? BigREAL64(buff) : LittleREAL64(buff);

    try{
      Write(&out, sizeof(real64));
    }
    catch(IOException& ex) {
      ex.SetDescription("An Error occured while writing a double");
      throw ex;
    }
  }


  void RawIO :: WriteDoubleArray(const real64 *buff , const int32 size) throw(IOException)
  {
    CHECK_WRITE();
    int32 s = bigendian_ ? BigINT32(size) : LittleINT32(size);

    try{
      Write(&s, sizeof(int32));
    }
    catch(IOException& ex) {
      ex.SetDescription("Error while writing the length of float array");
      throw ex;
    }

    real64 *array = new real64[size];
    for(int32 i=0; i<size; i++)
      array[i] = bigendian_ ? BigREAL64(buff[i]) : LittleREAL64(buff[i]);

    try{
      Write(array, sizeof(real64)*size);
    }
    catch(IOException& ex) {
      delete[] array;
      ex.SetDescription("An Error occured while writing a double array");
      throw ex;
    }

    delete[] array;
  }

  void RawIO :: WriteDoubleVector(const std::vector<real64>& vec) throw(IOException)
  {
    CHECK_WRITE();
    int32 size = vec.size();
    int32 s = bigendian_ ? BigINT32(size) : LittleINT32(size);

    try{
      Write(&s, sizeof(int32));
    }
    catch(IOException& ex) {
      ex.SetDescription("Error while writing the length of double vector");
      throw ex;
    }

    real64 *array = new real64[size];
    for(int32 i=0; i<size; i++)
      array[i] = bigendian_ ? BigREAL64(vec[i]) : LittleREAL64(vec[i]);

    try{
      Write(array, sizeof(real64)*size);
    }
    catch(IOException& ex) {
      delete[] array;
      ex.SetDescription("An Error occured while writing a double vector");
      throw ex;
    }

    delete[] array;
  }

  void RawIO :: Write(const void* data, int32 size) throw(IOException)
  {
    CHECK_WRITE();
    fwrite(data, size, 1, file_write);
    if(ferror(file_write)) {
      throw IOException ( "Error while writing to file_write" );
    }
    if(fflush(file_write) != 0) {
      throw IOException ( "Error while flushing file_write" );
    }
  }

  void RawIO :: Read(void* data, int32 size) throw(IOException)
  {
    CHECK_READ();
    fread(data, size, 1, file_read);
    if(ferror(file_read) || feof(file_read)) {
      throw IOException ( "Error while reading from file_read" );
    }
  }

  /*
    int32 main(int32 argc, char** argv)
    {
    FILE* out = fopen("test.out", "wb");
    FILE* in = fopen("test.in", "rb");
    RawIO *io = new RawIO(in, out, false);
    (*io) << "Test message";
    delete io;
    fclose(in);
    fclose(out);

    return 0;
    }
  */

}
