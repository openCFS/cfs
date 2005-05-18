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

RawIO::RawIO(Socket *sock, int32 timeout, bool bigendian) :
BaseIO(sock, timeout),
bigendian_(bigendian)
{
  InitEndian();
}


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
}

void RawIO :: readMsg(std::string& s) throw(IOException)
{
  CHECK_READ();
  int32 length = readInt();
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

int32 RawIO :: readInt() throw(IOException)
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


void RawIO :: readIntVector(std::vector<int32>& vec) throw(IOException)
{
  CHECK_READ();
  int32 size = readInt();
  int32 *buff= new int32[size];

  try{
    Read(buff, sizeof(int32)*size);
  }
  catch(IOException& ex) {
    delete[] buff;
    ex.SetDescription("An Error occured while reading an integer array");
    throw ex;
  }

  vec.resize(size);
  for(int32 i=0; i<size; i++) {
    vec[i] = bigendian_ ? BigINT32(buff[i]) : LittleINT32(buff[i]);
  }

  delete[] buff;
}

uint32 RawIO :: readUInt() throw(IOException)
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


void RawIO :: readUIntVector(std::vector<uint32>& vec) throw(IOException)
{
  CHECK_READ();
  int32 size = readInt();
  uint32 *buff= new uint32[size];

  try{
    Read(buff, sizeof(uint32)*size);
  }
  catch(IOException& ex) {
    delete[] buff;
    ex.SetDescription("An Error occured while reading an integer array");
    throw ex;
  }

  vec.resize(size);
  for(int32 i=0; i<size; i++) {
    vec[i] = bigendian_ ? BigINT32(buff[i]) : LittleINT32(buff[i]);
  }

  delete[] buff;
}


float32 RawIO :: readFloat() throw(IOException)
{
  CHECK_READ();
  float32 buff;

  try{
    Read(&buff, sizeof(float32));
  }
  catch(IOException& ex) {
    ex.SetDescription("An Error occured while reading a Float");
    throw ex;
  }

  buff = bigendian_ ? BigFLOAT32(buff) : LittleFLOAT32(buff);

  return(buff);
}


void RawIO :: readFloatVector(std::vector<float32>& vec) throw(IOException)
{
  CHECK_READ();
  int32 size = readInt();
  float32 *buff= new float32[size];

  try{
    Read(buff, sizeof(float32)*size);
  }
  catch(IOException& ex) {
    delete[] buff;
    ex.SetDescription("An Error occured while reading a float array");
    throw ex;
  }

  vec.resize(size);
  for(int32 i=0; i<size; i++) {
    vec[i] = bigendian_ ? BigFLOAT32(buff[i]) : LittleFLOAT32(buff[i]);
  }

  delete[] buff;
}

float64 RawIO :: readDouble() throw(IOException)
{
  CHECK_READ();
  float64 buff;

  try{
    Read(&buff, sizeof(float64));
  }
  catch(IOException& ex) {
    ex.SetDescription("An Error occured while reading a Double");
    throw ex;
  }

  buff = bigendian_ ? BigFLOAT64(buff) : LittleFLOAT64(buff);

  return(buff);
}


void RawIO :: readDoubleVector(std::vector<float64>& vec) throw(IOException)
{
  CHECK_READ();
  int32 size = readInt();
  float64 *buff= new float64[size];

  try{
    Read(buff, sizeof(float64)*size);
  }
  catch(IOException& ex) {
    delete[] buff;
    ex.SetDescription("An Error occured while reading a double array");
    throw ex;
  }

  vec.resize(size);
  for(int32 i=0; i<size; i++) {
    vec[i] = bigendian_ ? BigFLOAT64(buff[i]) : LittleFLOAT64(buff[i]);
  }

  delete[] buff;
}

void RawIO :: writeMsg(const std::string& s) throw(IOException)
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

  std::cout << "length: " << length << "\n";
  std::cout << "pad: " << pad << "\n";
  std::cout << "buf_size: " << bufsize << "\n";

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

void RawIO :: writeInt(const int32 buff) throw(IOException)
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


void RawIO :: writeIntArray(const int32 *buff , const int32 size) throw(IOException)
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

void RawIO :: writeIntVector(const std::vector<int32>& vec) throw(IOException)
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

void RawIO :: writeUInt(const uint32 buff) throw(IOException)
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


void RawIO :: writeUIntArray(const uint32 *buff , const int32 size) throw(IOException)
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

void RawIO :: writeUIntVector(const std::vector<uint32>& vec) throw(IOException)
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


void RawIO :: writeFloat(const float32 buff) throw(IOException)
{
  CHECK_WRITE();
  float32 out = bigendian_ ? BigFLOAT32(buff) : LittleFLOAT32(buff);

  try{
    Write(&out, sizeof(float32));
  }
  catch(IOException& ex) {
    ex.SetDescription("An Error occured while writing a float");
    throw ex;
  }
}


void RawIO :: writeFloatArray(const float32 *buff , const int32 size) throw(IOException)
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

  float32 *array = new float32[size];
  for(int32 i=0; i<size; i++)
     array[i] = bigendian_ ? BigFLOAT32(buff[i]) : LittleFLOAT32(buff[i]);

  try{
    Write(array, sizeof(float32)*size);
  }
  catch(IOException& ex) {
    delete[] array;
    ex.SetDescription("An Error occured while writing a float array");
    throw ex;
  }

  delete[] array;
}

void RawIO :: writeFloatVector(const std::vector<float32>& vec) throw(IOException)
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

  float32 *array = new float32[size];
  for(int32 i=0; i<size; i++)
     array[i] = bigendian_ ? BigFLOAT32(vec[i]) : LittleFLOAT32(vec[i]);

  try{
    Write(array, sizeof(float32)*size);
  }
  catch(IOException& ex) {
    delete[] array;
    ex.SetDescription("An Error occured while writing a float vector");
    throw ex;
  }

  delete[] array;
}

void RawIO :: writeDouble(const float64 buff) throw(IOException)
{
  CHECK_WRITE();
  float64 out = bigendian_ ? BigFLOAT64(buff) : LittleFLOAT64(buff);

  try{
    Write(&out, sizeof(float64));
  }
  catch(IOException& ex) {
    ex.SetDescription("An Error occured while writing a double");
    throw ex;
  }
}


void RawIO :: writeDoubleArray(const float64 *buff , const int32 size) throw(IOException)
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

  float64 *array = new float64[size];
  for(int32 i=0; i<size; i++)
     array[i] = bigendian_ ? BigFLOAT64(buff[i]) : LittleFLOAT64(buff[i]);

  try{
    Write(array, sizeof(float64)*size);
  }
  catch(IOException& ex) {
    delete[] array;
    ex.SetDescription("An Error occured while writing a double array");
    throw ex;
  }

  delete[] array;
}

void RawIO :: writeDoubleVector(const std::vector<float64>& vec) throw(IOException)
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

  float64 *array = new float64[size];
  for(int32 i=0; i<size; i++)
     array[i] = bigendian_ ? BigFLOAT64(vec[i]) : LittleFLOAT64(vec[i]);

  try{
    Write(array, sizeof(float64)*size);
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
