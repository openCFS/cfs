#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <memory>

#include <sys/signal.h>
#include <sys/ioctl.h>

#include "GSIXDRIO.hh"
#include "GSIIOException.hh"

using namespace std;

GSIXDRIO::GSIXDRIO(FILE *file_read, FILE *file_write) : GSIBaseIO(file_read, file_write)
{
  // If we get STREAMS as input we always do blocking IO

  if(file_read)
  {
      xdrstdio_create( &xdr_in, file_read, XDR_DECODE );
      desc_read = fileno(file_read);
  }
  
  if(file_write)
  {
      xdrstdio_create( &xdr_out, file_write, XDR_ENCODE );
      desc_write = fileno(file_write);
  }
  
  // If we get STREAMS as input we always do blocking IO
  io_ok = 1;

  /*
  if(file_read!=NULL) {
     xdrstdio_create( &xdr_in, file_read, XDR_DECODE );
  }
  else {
    throw GSIIOException ( "Could not open XDR input stream." );
  }

  if(file_write!=NULL) {
     xdrstdio_create( &xdr_out, file_write, XDR_ENCODE );
  }
  else {
    throw GSIIOException ( "Could not open XDR output stream." );
  }
  */
}

GSIXDRIO::GSIXDRIO(GSISocket *sock, int timeout) : GSIBaseIO(sock, timeout)
{

  //  file_read = sock_->getReadHandle();
  //  file_write = sock_->getWriteHandle();

  //  xdrstdio_create( &xdr_in, file_read, XDR_DECODE );
  //  xdrstdio_create( &xdr_out, file_write, XDR_ENCODE );

  //  desc_read = fileno(file_read);
  //  desc_write = fileno(file_write);

  // Here we have to differentiate between blocking an timed out
  if(timeout_ < 0) 
  {
      io_ok = 1;
  }
  else 
  {
      io_ok = 0;
  }

  std::cout << "XDRIO timeout " << timeout_ << std::endl;
}

int GSIXDRIO::genericRead(void* buf, int size)
{
  int ret;

  //  std::cout << "XDR timeout: " << timeout_ << std::endl;

  if(sock_)
    ret = sock_->read(buf, size, timeout_);
  else
  {
      sigset_t oldm, newm;
      
      if(sigemptyset(&newm) < 0)
        return -1;
      
      if(sigaddset(&newm, SIGPIPE) < 0)
        return -1;

      if(sigprocmask(SIG_BLOCK, &newm, &oldm) < 0)
        return -1;

      do 
      {
          ret = fread(buf, size, 1, file_read);
      } while ( ret != 1 && errno == EINTR);

      if(ferror(file_read))
        return -1;

      if(sigprocmask(SIG_SETMASK, &oldm, NULL) < 0)
        return -1;

      if(ret != 1)
        ret = -1;
  }

  return ret;
}

int GSIXDRIO::genericWrite(void* buf, int size)
{
  int ret = 0;
  
  if(sock_)
    ret = sock_->write(buf, size, timeout_);
  else
  {
      //      while(ret != size) 
      //      {
      //          std::cout << "fw ret: " << ret << " size " << size << std::endl;

      sigset_t oldm, newm;
      
      if(sigemptyset(&newm) < 0)
        return -1;
      
      if(sigaddset(&newm, SIGPIPE) < 0)
        return -1;

      if(sigprocmask(SIG_BLOCK, &newm, &oldm) < 0)
        return -1;

      do 
      {
          ret = fwrite(buf, size, 1, file_write);
      } while ( ret != 1 && errno == EINTR);

          //      }
      

      //      std::cout << "fwrite ret: " << ret << std::endl;
      
      //      int err = ferror(file_write);

      //      std::cout << "ferror ret: " << err << std::endl;      

      if(ferror(file_write))
        return -1;

      //      perror("generic Write");

      if(ret != 1)
        return -1;

      do 
      {
          ret = fflush(file_write);
      } while ( ret != 0 && errno == EINTR);
      
      if(sigprocmask(SIG_SETMASK, &oldm, NULL) < 0)
        return -1;
      
      if(ret != 0)
        ret = -1;
      else
        ret = 1;
  }

  //  std::cout << "genericWrite ret " << ret<< std::endl;

  return ret;
}

// reads an int and won't throw any exceptions
int GSIXDRIO::genericReadInt(int& val)
{
  char* buf = new char[4];
  int i, ret;
  
  ret = genericRead(buf, 4);
  
  if(ret <= 0) {
      delete[] buf;
      return ret;
  }

  XDR memhandle;
  
  xdrmem_create(&memhandle, buf, 4, XDR_DECODE);
  
  ret = xdr_int(&memhandle, &i);
  
  xdr_destroy(&memhandle);
  
  delete[] buf;

  if(ret == 0) {
      return -1;
  }

  val = i;
  
  return 1;
}

int GSIXDRIO::genericWriteInt(int val)
{
  char* buf = new char[4];
  int ret;

  XDR memhandle;
  
  xdrmem_create(&memhandle, buf, 4, XDR_ENCODE);
  
  ret = xdr_int(&memhandle, &val);
  
  if(ret == 0) {
      xdr_destroy(&memhandle);
      delete[] buf;
      return -1;
  }
 
  ret = genericWrite(buf, xdr_getpos(&memhandle));
  
  xdr_destroy(&memhandle);
  delete[] buf;

  return ret;
}


void GSIXDRIO :: readMsg(std::string& s) throw(GSIIOException)
{
  int ret;
  int size;

  is_valid();

  // read the message size
  ret = genericReadInt(size);

  if(ret < 0) 
  {
      throw GSIIOException ( "Error while reading the message size" );
  }
  
  // a timeout occured
  if(ret == 0)
  {
      io_ok = 0;
      return;
  }
  
  char* buf = new char[size];
  //  std::cout << "buffersize: "<< size << std::endl;
  

  ret = genericRead(buf, size);
  
  if(ret < 0) {
      delete[] buf;
      throw GSIIOException ( "Error while reading the message" );
  }
  
  if( ret == 0 ) 
  {
      delete[] buf;
      io_ok = 0;
      return;
  }

  XDR memhandle;

  xdrmem_create(&memhandle, buf, size, XDR_DECODE);
  
  char* str = new char[size+1];
  ret = xdr_string(&memhandle, &str, size);
  
  xdr_destroy(&memhandle);
  
  
  if( ret == 0 )
  {
      delete[] str;
      delete[] buf;
      throw GSIIOException ( "An Error occured while reading the message");
  }
  
          
  str[size] = 0;
  s = str;
  delete[] str;
  delete[] buf;
  io_ok = 1;
}

int GSIXDRIO :: readInt() throw(GSIIOException)
{
  int ret;
  int i;
  
  is_valid();

  ret = genericReadInt(i);
  
  if(ret < 0) {
      throw GSIIOException ( "An Error occured while reading an integer" );
  }

  if( ret == 0 ) 
  {
      io_ok = 0;
      return 0;
  }

  io_ok = 1;
  return i;
}


void GSIXDRIO :: readIntVector(std::vector<int>& vec) throw(GSIIOException)
{
  int size, length;
  int ret;
  int val;

  is_valid();

  // read the message length
  ret = genericReadInt(size); // size in bytes

  //  std::cout << "intvec size " << size << std::endl;
  

  if(ret < 0) 
  {
      throw GSIIOException ( "Error while reading the length of IntVecData" );
  }
  
  // a timeout occured
  if(ret == 0)
  {
      io_ok = 0;
      return;
  }
  
  // read the pure message
  char* buf = new char[size];
  
  ret = genericRead(buf, size);
  
  if(ret < 0) {
      delete[] buf;
      throw GSIIOException ( "Error while reading the IntVec" );
  }
  
  if( ret == 0 ) 
  {
      delete[] buf;
      io_ok = 0;
      return;
  }

  // create a XDR memory handle from the buffer
  XDR memhandle;
  
  xdrmem_create(&memhandle, buf, size, XDR_DECODE);
  
  ret = xdr_int(&memhandle, &length);

  //  std::cout << "intvec length " << length << std::endl;

  if( ret == 0 )
  {
      xdr_destroy(&memhandle);
      delete[] buf;
      throw GSIIOException ( "An Error occured while decoding IntVecLength");
  }  
  
  int* intbuf = new int[length];

  for(int l=0; l < length; l++) 
  {
      ret = xdr_int(&memhandle, &val);

      if( ret == 0 )
      {
          delete[] intbuf;
          xdr_destroy(&memhandle);
          delete[] buf;
          throw GSIIOException ( "An Error occured while decoding the IntVec");
      }
      intbuf[l] = val;
  }

  vec.clear();
 
  for(int l=0; l < length; l++) 
  {
      vec.push_back(intbuf[l]);
  }

  xdr_destroy(&memhandle);
  delete[] intbuf;
  delete[] buf;
  io_ok = 1;
}


u_int GSIXDRIO :: readUInt() throw(GSIIOException)
{
  u_int buff;
  int ret;

  is_valid();
  
  char* buf = new char[4];

  ret = genericRead(buf, 4);
  
  if(ret < 0) {
      delete[] buf;
      throw GSIIOException ( "Error while reading a UInt" );
  }

  if(ret == 0)
  {
      delete[] buf;
      io_ok = 0;
      return 0;
  }

  XDR memhandle;
  
  xdrmem_create(&memhandle, buf, 4, XDR_DECODE);
  
  ret = xdr_u_int(&memhandle, &buff);
  
  xdr_destroy(&memhandle);
  
  delete[] buf;

  if(ret == 0) {
      throw GSIIOException ( "Error while decoding a UInt" );
  }

  io_ok = 1;
  return buff;
}


void GSIXDRIO :: readUIntVector(std::vector<u_int>& vec) throw(GSIIOException)
{
  int size, length;
  int ret;
  u_int val;

  is_valid();

  // read the message length
  ret = genericReadInt(size); // size in bytes

  if(ret < 0) 
  {
      throw GSIIOException ( "Error while reading the length of UIntVecData" );
  }
  
  // a timeout occured
  if(ret == 0)
  {
      io_ok = 0;
      return;
  }
  
  // read the pure message
  char* buf = new char[size];
  
  ret = genericRead(buf, size);
  
  if(ret < 0) {
      delete[] buf;
      throw GSIIOException ( "Error while reading the UIntVec" );
  }
  
  if( ret == 0 ) 
  {
      delete[] buf;
      io_ok = 0;
      return;
  }

  // create a XDR memory handle from the buffer
  XDR memhandle;
  
  xdrmem_create(&memhandle, buf, size, XDR_DECODE);
  
  ret = xdr_int(&memhandle, &length);

  if( ret == 0 )
  {
      xdr_destroy(&memhandle);
      delete[] buf;
      throw GSIIOException ( "An Error occured while decoding UIntVecLength");
  }  
  
  u_int* uintbuf = new u_int[length];

  for(int l=0; l < length; l++) 
  {
      ret = xdr_u_int(&memhandle, &val);

      if( ret == 0 )
      {
          delete[] uintbuf;
          xdr_destroy(&memhandle);
          delete[] buf;
          throw GSIIOException ( "An Error occured while decoding the UIntVec");
      }
      uintbuf[l] = val;
  }
  
  vec.clear();
  
  for(int l=0; l < length; l++) 
  {
      vec.push_back(uintbuf[l]);
  }

  xdr_destroy(&memhandle);
  delete[] uintbuf;
  delete[] buf;
  io_ok = 1;
}


float GSIXDRIO :: readFloat() throw(GSIIOException)
{
  float buff;
  int ret;

  is_valid();
  
  char* buf = new char[4];

  ret = genericRead(buf, 4);
  
  if(ret < 0) {
      delete[] buf;
      throw GSIIOException ( "Error while reading a Float" );
  }

  if(ret == 0)
  {
      delete[] buf;
      io_ok = 0;
      return 0.0f;
  }

  XDR memhandle;
  
  xdrmem_create(&memhandle, buf, 4, XDR_DECODE);
  
  ret = xdr_float(&memhandle, &buff);
  
  xdr_destroy(&memhandle);
  
  delete[] buf;

  if(ret == 0) {
      throw GSIIOException ( "Error while decoding a Float" );
  }

  io_ok = 1;
  return buff;
}



void GSIXDRIO :: readFloatVector(std::vector<float>& vec) throw(GSIIOException)
{
  int size, length;
  int ret;
  float val;

  is_valid();

  // read the message length
  ret = genericReadInt(size); // size in bytes

  if(ret < 0) 
  {
      throw GSIIOException ( "Error while reading the length of FloatVecData" );
  }
  
  // a timeout occured
  if(ret == 0)
  {
      io_ok = 0;
      return;
  }
  
  // read the pure message
  char* buf = new char[size];
  
  ret = genericRead(buf, size);
  
  if(ret < 0) {
      delete[] buf;
      throw GSIIOException ( "Error while reading the FloatVec" );
  }
  
  if( ret == 0 ) 
  {
      delete[] buf;
      io_ok = 0;
      return;
  }

  // create a XDR memory handle from the buffer
  XDR memhandle;
  
  xdrmem_create(&memhandle, buf, size, XDR_DECODE);
  
  ret = xdr_int(&memhandle, &length);

  if( ret == 0 )
  {
      xdr_destroy(&memhandle);
      delete[] buf;
      throw GSIIOException ( "An Error occured while decoding FloatVecLength");
  }  
  
  float* floatbuf = new float[length];

  for(int l=0; l < length; l++) 
  {
      ret = xdr_float(&memhandle, &val);

      if( ret == 0 )
      {
          delete[] floatbuf;
          xdr_destroy(&memhandle);
          delete[] buf;
          throw GSIIOException ( "An Error occured while decoding the FloatVec");
      }
      floatbuf[l] = val;
  }
  
  vec.clear();
  
  for(int l=0; l < length; l++) 
  {
      vec.push_back(floatbuf[l]);
  }

  xdr_destroy(&memhandle);
  delete[] floatbuf;
  delete[] buf;
  io_ok = 1;
}

double GSIXDRIO :: readDouble() throw(GSIIOException)
{
  double buff;
  int ret;

  is_valid();
  
  char* buf = new char[8];

  ret = genericRead(buf, 8);
  
  if(ret < 0) {
      delete[] buf;
      throw GSIIOException ( "Error while reading a Double" );
  }

  if(ret == 0)
  {
      delete[] buf;
      io_ok = 0;
      return 0.0;
  }

  XDR memhandle;
  
  xdrmem_create(&memhandle, buf, 8, XDR_DECODE);
  
  ret = xdr_double(&memhandle, &buff);
  
  xdr_destroy(&memhandle);
  
  delete[] buf;

  if(ret == 0) {
      throw GSIIOException ( "Error while decoding a Double" );
  }

  io_ok = 1;
  return buff;
}


void GSIXDRIO :: readDoubleVector(std::vector<double>& vec) throw(GSIIOException)
{
  int size, length;
  int ret;
  double val;

  is_valid();

  // read the message length
  ret = genericReadInt(size); // size in bytes

  if(ret < 0) 
  {
      throw GSIIOException ( "Error while reading the length of FloatVecData" );
  }
  
  // a timeout occured
  if(ret == 0)
  {
      io_ok = 0;
      return;
  }
  
  // read the pure message
  char* buf = new char[size];
  
  ret = genericRead(buf, size);
  
  if(ret < 0) {
      delete[] buf;
      throw GSIIOException ( "Error while reading the FloatVec" );
  }
  
  if( ret == 0 ) 
  {
      delete[] buf;
      io_ok = 0;
      return;
  }

  // create a XDR memory handle from the buffer
  XDR memhandle;
  
  xdrmem_create(&memhandle, buf, size, XDR_DECODE);
  
  ret = xdr_int(&memhandle, &length);

  if( ret == 0 )
  {
      xdr_destroy(&memhandle);
      delete[] buf;
      throw GSIIOException ( "An Error occured while decoding FloatVecLength");
  }  
  
  double* doublebuf = new double[length];

  for(int l=0; l < length; l++) 
  {
      ret = xdr_double(&memhandle, &val);

      if( ret == 0 )
      {
          delete[] doublebuf;
          xdr_destroy(&memhandle);
          delete[] buf;
          throw GSIIOException ( "An Error occured while decoding the FloatVec");
      }
      doublebuf[l] = val;
  }
  
  vec.clear();
  
  for(int l=0; l < length; l++) 
  {
      vec.push_back(doublebuf[l]);
  }

  xdr_destroy(&memhandle);
  delete[] doublebuf;
  delete[] buf;
  io_ok = 1;
}

void GSIXDRIO :: writeMsg(const std::string& s) throw(GSIIOException)
{
  int ret;
  int size;

  is_valid();

  int bufsize = s.size()+128;
  char* buf = new char[bufsize];

  XDR memhandle;

  xdrmem_create(&memhandle, buf, bufsize, XDR_ENCODE);

  char* str = strdup(s.c_str());
  
  ret = xdr_string(&memhandle, &str, strlen(str));

  free(str);

  if( ret == 0 )
  {
      xdr_destroy(&memhandle);
      delete[] buf;
      throw GSIIOException ( "An Error occured while encoding the message");
  }

  int pos = xdr_getpos(&memhandle);

  xdr_destroy(&memhandle);


  //  std::cout << "writeMsg xdr_getpos: " << pos << " " << timeout_ << std::endl;

  ret = genericWriteInt(pos);

  if(ret < 0) 
  {
      delete[] buf;
      throw GSIIOException ( "Error while writing the message size" );
  }
  
  // a timeout occured
  if(ret == 0)
  {
      delete[] buf;
      io_ok = 0;
      return;
  }
  
  ret = genericWrite(buf, pos);

  if(ret < 0) 
  {
      delete[] buf;
      throw GSIIOException ( "Error while writing the message" );
  }
  
  // a timeout occured
  if(ret == 0)
  {
      delete[] buf;
      io_ok = 0;
      return;
  }

  //  std::cout << "writeMsg ret " << ret << std::endl;

  delete[] buf;
  io_ok = 1;
  //  std::cout << "io_ok " << ioOK() << std::endl;
}

void GSIXDRIO :: writeInt(const int buff) throw(GSIIOException)
{
  int ret;
  int i = buff;
  
  is_valid();

  ret = genericWriteInt(i);
  
  if(ret < 0) {
      throw GSIIOException ( "An Error occured while writing an integer" );
  }

  if( ret == 0 ) 
  {
      io_ok = 0;
      return;
  }

  io_ok = 1;
}


void GSIXDRIO :: writeIntArray(const int *buff , const int size) throw(GSIIOException)
{
  int ret;
  int bufsize;
  
  is_valid();

  bufsize = size*4+128;
  char* buf = new char[size*4+128];

  XDR memhandle;

  xdrmem_create(&memhandle, buf, bufsize, XDR_ENCODE);

  int sz = size;
  ret = xdr_int(&memhandle, &sz);
  
  if( ret == 0 )
  {
      xdr_destroy(&memhandle);
      delete[] buf;
      throw GSIIOException ( "An Error occured while encoding the IntArrayLength");
  }

  for(int i=0; i<size; i++) 
  {
      int val = buff[i];

      ret = xdr_int(&memhandle, &val);

      if( ret == 0 )
      {
          xdr_destroy(&memhandle);
          delete[] buf;
          throw GSIIOException ( "An Error occured while encoding the IntArray");
      }
  }

  int pos = xdr_getpos(&memhandle);

  xdr_destroy(&memhandle);


  //  std::cout<< "writeIntArray xdr_getpos: " << pos << std::endl;

  ret = genericWriteInt(pos);

  if(ret < 0) 
  {
      delete[] buf;
      throw GSIIOException ( "Error while writing the IA message size" );
  }
  
  // a timeout occured
  if(ret == 0)
  {
      delete[] buf;
      io_ok = 0;
      return;
  }
  
  ret = genericWrite(buf, pos);

  if(ret < 0) 
  {
      delete[] buf;
      throw GSIIOException ( "Error while writing the IntArray" );
  }
  
  // a timeout occured
  if(ret == 0)
  {
      delete[] buf;
      io_ok = 0;
      return;
  }

  delete[] buf;
  io_ok = 1;
}

void GSIXDRIO :: writeIntVector(const std::vector<int>& vec) throw(GSIIOException)
{
  int ret;
  int size = vec.size();
  int bufsize;

  is_valid();

  bufsize = size*4+128;
  char* buf = new char[bufsize];

  XDR memhandle;

  xdrmem_create(&memhandle, buf, bufsize, XDR_ENCODE);
  
  ret = xdr_int(&memhandle, &size);
  
  if( ret == 0 )
  {
      xdr_destroy(&memhandle);
      delete[] buf;
      throw GSIIOException ( "An Error occured while encoding the IntVecLength");
  }

  for(int i=0; i<vec.size(); i++) 
  {
      int val = vec[i];

      //      std::cout << val << std::endl;

      ret = xdr_int(&memhandle, &val);

      if( ret == 0 )
      {
          xdr_destroy(&memhandle);
          delete[] buf;
          throw GSIIOException ( "An Error occured while encoding the IntVector");
      }
  }

  int pos = xdr_getpos(&memhandle);

  xdr_destroy(&memhandle);


  //  std::cout << "writeIntArray xdr_getpos: " << pos << std::endl;

  ret = genericWriteInt(pos);

  if(ret < 0) 
  {
      delete[] buf;
      throw GSIIOException ( "Error while writing the IV message size" );
  }
  
  // a timeout occured
  if(ret == 0)
  {
      delete[] buf;
      io_ok = 0;
      return;
  }
  
  ret = genericWrite(buf, pos);

  //  std::cout << "genericWrite IntVec ret: " << ret << std::endl;

  if(ret < 0) 
  {
      delete[] buf;
      throw GSIIOException ( "Error while writing the IntVec" );
  }
  
  // a timeout occured
  if(ret == 0)
  {
      delete[] buf;
      io_ok = 0;
      return;
  }

  io_ok = 1;
  delete[] buf;
}

void GSIXDRIO :: writeUInt(const u_int val) throw(GSIIOException)
{
  u_int buff=val;
  int ret;

  is_valid();
  
  char* buf = new char[4];

  XDR memhandle;
  
  xdrmem_create(&memhandle, buf, 4, XDR_ENCODE);
  
  ret = xdr_u_int(&memhandle, &buff);
  
  int pos = xdr_getpos(&memhandle);

  xdr_destroy(&memhandle);
  
  if(ret == 0) {
      delete[] buf;
      throw GSIIOException ( "Error while encoding a UInt" );
  }

  ret = genericWrite(buf, pos);
  
  if(ret < 0) {
      delete[] buf;
      throw GSIIOException ( "Error while reading a UInt" );
  }

  if(ret == 0)
  {
      delete[] buf;
      io_ok = 0;
      return;
  }

  delete[] buf;
  io_ok = 1;
}



void GSIXDRIO :: writeUIntArray(const u_int *buff , const int size) throw(GSIIOException)
{
  int ret;
  int bufsize;

  is_valid();

  bufsize = size*4+128;
  char* buf = new char[bufsize];

  XDR memhandle;

  xdrmem_create(&memhandle, buf, bufsize, XDR_ENCODE);

  int sz = size;
  ret = xdr_int(&memhandle, &sz);
  
  if( ret == 0 )
  {
      xdr_destroy(&memhandle);
      delete[] buf;
      throw GSIIOException ( "An Error occured while encoding the UIntArrayLength");
  }

  for(int i=0; i<size; i++) 
  {
      u_int val = buff[i];

      //      std::cout << val << std::endl;

      ret = xdr_u_int(&memhandle, &val);

      if( ret == 0 )
      {
          xdr_destroy(&memhandle);
          delete[] buf;
          throw GSIIOException ( "An Error occured while encoding the UIntArray");
      }
  }

  int pos = xdr_getpos(&memhandle);

  xdr_destroy(&memhandle);


  //  std::cout << "writeIntArray xdr_getpos: " << pos << std::endl;

  ret = genericWriteInt(pos);

  if(ret < 0) 
  {
      delete[] buf;
      throw GSIIOException ( "Error while writing the UIA message size" );
  }
  
  // a timeout occured
  if(ret == 0)
  {
      delete[] buf;
      io_ok = 0;
      return;
  }
  
  ret = genericWrite(buf, pos);

  if(ret < 0) 
  {
      delete[] buf;
      throw GSIIOException ( "Error while writing the UIntArray" );
  }
  
  // a timeout occured
  if(ret == 0)
  {
      delete[] buf;
      io_ok = 0;
      return;
  }

  io_ok = 1;
  delete[] buf;
}

void GSIXDRIO :: writeUIntVector(const std::vector<u_int>& vec) throw(GSIIOException)
{
  int ret;
  int size = vec.size();
  int bufsize;

  is_valid();

  bufsize = size*4+128;
  char* buf = new char[bufsize];

  XDR memhandle;

  xdrmem_create(&memhandle, buf, bufsize, XDR_ENCODE);
  
  ret = xdr_int(&memhandle, &size);
  
  if( ret == 0 )
  {
      xdr_destroy(&memhandle);
      delete[] buf;
      throw GSIIOException ( "An Error occured while encoding the UIntVecLength");
  }

  for(int i=0; i<vec.size(); i++) 
  {
      u_int val = vec[i];

      //      std::cout << val << std::endl;

      ret = xdr_u_int(&memhandle, &val);

      if( ret == 0 )
      {
          xdr_destroy(&memhandle);
          delete[] buf;
          throw GSIIOException ( "An Error occured while encoding the UIntVector");
      }
  }

  int pos = xdr_getpos(&memhandle);

  xdr_destroy(&memhandle);


  //  std::cout << "writeIntArray xdr_getpos: " << pos << std::endl;

  ret = genericWriteInt(pos);

  if(ret < 0) 
  {
      delete[] buf;
      throw GSIIOException ( "Error while writing the UIV message size" );
  }
  
  // a timeout occured
  if(ret == 0)
  {
      delete[] buf;
      io_ok = 0;
      return;
  }
  
  ret = genericWrite(buf, pos);

  if(ret < 0) 
  {
      delete[] buf;
      throw GSIIOException ( "Error while writing the UIntVec" );
  }
  
  // a timeout occured
  if(ret == 0)
  {
      delete[] buf;
      io_ok = 0;
      return;
  }

  io_ok = 1;
  delete[] buf;
}

void GSIXDRIO :: writeFloat(const float val) throw(GSIIOException)
{
  float buff=val;
  int ret;

  is_valid();
  
  char* buf = new char[4];

  XDR memhandle;
  
  xdrmem_create(&memhandle, buf, 4, XDR_ENCODE);
  
  ret = xdr_float(&memhandle, &buff);
  
  int pos = xdr_getpos(&memhandle);

  xdr_destroy(&memhandle);
  
  if(ret == 0) {
      delete[] buf;
      throw GSIIOException ( "Error while encoding a Float" );
  }

  ret = genericWrite(buf, pos);
  
  if(ret < 0) {
      delete[] buf;
      throw GSIIOException ( "Error while writing a Float" );
  }

  if(ret == 0)
  {
      delete[] buf;
      io_ok = 0;
      return;
  }

  delete[] buf;
  io_ok = 1;
}


void GSIXDRIO :: writeFloatArray(const float *buff , const int size) throw(GSIIOException)
{
  int ret;
  int bufsize;

  is_valid();

  bufsize = size*4+128;
  char* buf = new char[bufsize];

  XDR memhandle;

  xdrmem_create(&memhandle, buf, bufsize, XDR_ENCODE);

  int sz = size;
  ret = xdr_int(&memhandle, &sz);
  
  if( ret == 0 )
  {
      xdr_destroy(&memhandle);
      delete[] buf;
      throw GSIIOException ( "An Error occured while encoding the FloatArrayLength");
  }

  for(int i=0; i<size; i++) 
  {
      float val = buff[i];

      //      std::cout << val << std::endl;

      ret = xdr_float(&memhandle, &val);

      if( ret == 0 )
      {
          xdr_destroy(&memhandle);
          delete[] buf;
          throw GSIIOException ( "An Error occured while encoding the FLoatArray");
      }
  }

  int pos = xdr_getpos(&memhandle);

  xdr_destroy(&memhandle);


  ret = genericWriteInt(pos);

  if(ret < 0) 
  {
      delete[] buf;
      throw GSIIOException ( "Error while writing the FA message size" );
  }
  
  // a timeout occured
  if(ret == 0)
  {
      delete[] buf;
      io_ok = 0;
      return;
  }
  
  ret = genericWrite(buf, pos);

  if(ret < 0) 
  {
      delete[] buf;
      throw GSIIOException ( "Error while writing the FloatArray" );
  }
  
  // a timeout occured
  if(ret == 0)
  {
      delete[] buf;
      io_ok = 0;
      return;
  }

  io_ok = 1;
  delete[] buf;
}

void GSIXDRIO :: writeFloatVector(const std::vector<float>& vec) throw(GSIIOException)
{
  int ret;
  int size = vec.size();
  int bufsize;

  is_valid();

  bufsize = size*4+128;
  char* buf = new char[bufsize];

  XDR memhandle;

  xdrmem_create(&memhandle, buf, bufsize, XDR_ENCODE);
  
  ret = xdr_int(&memhandle, &size);
  
  if( ret == 0 )
  {
      xdr_destroy(&memhandle);
      delete[] buf;
      throw GSIIOException ( "An Error occured while encoding the FloatVecLength");
  }

  for(int i=0; i<vec.size(); i++) 
  {
      float val = vec[i];

      //      std::cout << val << std::endl;

      ret = xdr_float(&memhandle, &val);

      if( ret == 0 )
      {
          xdr_destroy(&memhandle);
          delete[] buf;
          throw GSIIOException ( "An Error occured while encoding the FloatVector");
      }
  }

  int pos = xdr_getpos(&memhandle);

  xdr_destroy(&memhandle);

  ret = genericWriteInt(pos);

  if(ret < 0) 
  {
      delete[] buf;
      throw GSIIOException ( "Error while writing the FLOATV message size" );
  }
  
  // a timeout occured
  if(ret == 0)
  {
      delete[] buf;
      io_ok = 0;
      return;
  }
  
  ret = genericWrite(buf, pos);

  if(ret < 0) 
  {
      delete[] buf;
      throw GSIIOException ( "Error while writing the FloatVec" );
  }
  
  // a timeout occured
  if(ret == 0)
  {
      delete[] buf;
      io_ok = 0;
      return;
  }

  io_ok = 1;
  delete[] buf;
}

void GSIXDRIO :: writeDouble(const double val) throw(GSIIOException)
{
  double buff=val;
  int ret;

  is_valid();
  
  char* buf = new char[8];

  XDR memhandle;
  
  xdrmem_create(&memhandle, buf, 8, XDR_ENCODE);
  
  ret = xdr_double(&memhandle, &buff);
  
  int pos = xdr_getpos(&memhandle);

  xdr_destroy(&memhandle);
  
  if(ret == 0) {
      delete[] buf;
      throw GSIIOException ( "Error while encoding a Double" );
  }

  ret = genericWrite(buf, pos);
  
  if(ret < 0) {
      delete[] buf;
      throw GSIIOException ( "Error while writing a Double" );
  }

  if(ret == 0)
  {
      delete[] buf;
      io_ok = 0;
      return;
  }

  delete[] buf;
  io_ok = 1;
}


void GSIXDRIO :: writeDoubleArray(const double *buff , const int size) throw(GSIIOException)
{
  int ret;
  int bufsize;

  is_valid();

  bufsize = size*8+128;
  char* buf = new char[bufsize];

  XDR memhandle;

  xdrmem_create(&memhandle, buf, bufsize, XDR_ENCODE);

  int sz = size;
  ret = xdr_int(&memhandle, &sz);
  
  if( ret == 0 )
  {
      xdr_destroy(&memhandle);
      delete[] buf;
      throw GSIIOException ( "An Error occured while encoding the DoubleArrayLength");
  }

  for(int i=0; i<size; i++) 
  {
      double val = buff[i];

      //      std::cout << val << std::endl;

      ret = xdr_double(&memhandle, &val);

      if( ret == 0 )
      {
          xdr_destroy(&memhandle);
          delete[] buf;
          throw GSIIOException ( "An Error occured while encoding the DoubleArray");
      }
  }

  int pos = xdr_getpos(&memhandle);

  xdr_destroy(&memhandle);


  ret = genericWriteInt(pos);

  if(ret < 0) 
  {
      delete[] buf;
      throw GSIIOException ( "Error while writing the DA message size" );
  }
  
  // a timeout occured
  if(ret == 0)
  {
      delete[] buf;
      io_ok = 0;
      return;
  }
  
  ret = genericWrite(buf, pos);

  if(ret < 0) 
  {
      delete[] buf;
      throw GSIIOException ( "Error while writing the DoubleArray" );
  }
  
  // a timeout occured
  if(ret == 0)
  {
      delete[] buf;
      io_ok = 0;
      return;
  }

  io_ok = 1;
  delete[] buf;
}

void GSIXDRIO :: writeDoubleVector(const std::vector<double>& vec) throw(GSIIOException)
{
  int ret;
  int size = vec.size();
  int bufsize;

  is_valid();

  bufsize = size*8+128;
  char* buf = new char[bufsize];

  XDR memhandle;

  xdrmem_create(&memhandle, buf, bufsize, XDR_ENCODE);
  
  ret = xdr_int(&memhandle, &size);
  
  if( ret == 0 )
  {
      xdr_destroy(&memhandle);
      delete[] buf;
      throw GSIIOException ( "An Error occured while encoding the DoubleVecLength");
  }

  for(int i=0; i<vec.size(); i++) 
  {
      double val = vec[i];

      //      std::cout << val << std::endl;

      ret = xdr_double(&memhandle, &val);

      if( ret == 0 )
      {
          xdr_destroy(&memhandle);
          delete[] buf;
          throw GSIIOException ( "An Error occured while encoding the DoubleVector");
      }
  }

  int pos = xdr_getpos(&memhandle);

  xdr_destroy(&memhandle);

  ret = genericWriteInt(pos);

  if(ret < 0) 
  {
      delete[] buf;
      throw GSIIOException ( "Error while writing the DoubleV message size" );
  }
  
  // a timeout occured
  if(ret == 0)
  {
      delete[] buf;
      io_ok = 0;
      return;
  }
  
  ret = genericWrite(buf, pos);

  if(ret < 0) 
  {
      delete[] buf;
      throw GSIIOException ( "Error while writing the DoubleVec" );
  }
  
  // a timeout occured
  if(ret == 0)
  {
      delete[] buf;
      io_ok = 0;
      return;
  }

  io_ok = 1;
  delete[] buf;  
}
