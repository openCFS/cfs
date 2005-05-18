#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>

#include "GSIXDRIO.hh"

namespace GridlibSocketInterface
{

XDRIO::XDRIO(FILE *file_read, FILE *file_write, bool bigendian) :
RawIO(file_read, file_write, bigendian)
{
  // If we get STREAMS as input we always do blocking IO
/*
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
*/
  /*
  if(file_read!=NULL) {
     xdrstdio_create( &xdr_in, file_read, XDR_DECODE );
  }
  else {
    throw IOException ( "Could not open XDR input stream." );
  }

  if(file_write!=NULL) {
     xdrstdio_create( &xdr_out, file_write, XDR_ENCODE );
  }
  else {
    throw IOException ( "Could not open XDR output stream." );
  }
  */
}

XDRIO::XDRIO(Socket *sock, int timeout, bool bigendian) :
RawIO(sock, timeout, bigendian)
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

/*
int XDRIO::genericRead(void* buf, int size)
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

int XDRIO::genericWrite(void* buf, int size)
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
int XDRIO::genericReadInt(int& val)
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

int XDRIO::genericWriteInt(int val)
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
*/

void XDRIO :: readMsg(std::string& s) throw(IOException)
{
  RawIO::readInt();
  RawIO::readMsg(s);
}

void XDRIO :: readIntVector(std::vector<int32>& vec) throw(IOException)
{
  RawIO::readInt();
  RawIO::readIntVector(vec);
}

void XDRIO :: readUIntVector(std::vector<uint32>& vec) throw(IOException)
{
  RawIO::readInt();
  RawIO::readUIntVector(vec);
}


void XDRIO :: readFloatVector(std::vector<float32>& vec) throw(IOException)
{
  RawIO::readInt();
  RawIO::readFloatVector(vec);
}


void XDRIO :: readDoubleVector(std::vector<float64>& vec) throw(IOException)
{
  RawIO::readInt();
  RawIO::readDoubleVector(vec);
}

void XDRIO :: writeMsg(const std::string& s) throw(IOException)
{
  int32 length = s.length();
  int32 pad = 4 - (length % 4);
  int32 bufsize = (pad == 4) ? length : length + pad;

  RawIO::writeInt(bufsize+4);
  RawIO::writeMsg(s);
}

void XDRIO :: writeIntArray(const int32 *buff , const int32 size) throw(IOException)
{
  RawIO::writeInt(size*sizeof(int32)+4);
  RawIO::writeIntArray(buff, size);
}

void XDRIO :: writeIntVector(const std::vector<int32>& vec) throw(IOException)
{
  RawIO::writeInt(vec.size()*sizeof(int32)+4);
  RawIO::writeIntVector(vec);
}

void XDRIO :: writeUIntArray(const uint32 *buff , const int32 size) throw(IOException)
{
  RawIO::writeInt(size*sizeof(uint32)+4);
  RawIO::writeUIntArray(buff, size);
}

void XDRIO :: writeUIntVector(const std::vector<uint32>& vec) throw(IOException)
{
  RawIO::writeInt(vec.size()*sizeof(uint32)+4);
  RawIO::writeUIntVector(vec);
}


void XDRIO :: writeFloatArray(const float32 *buff , const int32 size) throw(IOException)
{
  RawIO::writeInt(size*sizeof(float32)+4);
  RawIO::writeFloatArray(buff, size);
}

void XDRIO :: writeFloatVector(const std::vector<float32>& vec) throw(IOException)
{
  RawIO::writeInt(vec.size()*sizeof(float32)+4);
  RawIO::writeFloatVector(vec);
}

void XDRIO :: writeDoubleArray(const float64 *buff , const int32 size) throw(IOException)
{
  RawIO::writeInt(size*sizeof(float64)+4);
  RawIO::writeDoubleArray(buff, size);
}

void XDRIO :: writeDoubleVector(const std::vector<float64>& vec) throw(IOException)
{
  RawIO::writeInt(vec.size()*sizeof(float64)+4);
  RawIO::writeDoubleVector(vec);
}

}
