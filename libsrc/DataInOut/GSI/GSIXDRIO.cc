#ifdef USE_RCSID
static const char RCSid_GSIXDRIO[] = "$Id$";
#endif

/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.4  2005/05/18 19:26:03  strieben
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

#ifndef WIN32
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
#endif

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

void XDRIO :: ReadMsg(std::string& s) throw(IOException)
{
  RawIO::ReadInt();
  RawIO::ReadMsg(s);
}

void XDRIO :: ReadShortVector(std::vector<int16>& vec) throw(IOException)
{
  RawIO::ReadInt();
  RawIO::ReadShortVector(vec);
}

void XDRIO :: ReadUShortVector(std::vector<uint16>& vec) throw(IOException)
{
  RawIO::ReadInt();
  RawIO::ReadUShortVector(vec);
}

void XDRIO :: ReadIntVector(std::vector<int32>& vec) throw(IOException)
{
  RawIO::ReadInt();
  RawIO::ReadIntVector(vec);
}

void XDRIO :: ReadUIntVector(std::vector<uint32>& vec) throw(IOException)
{
  RawIO::ReadInt();
  RawIO::ReadUIntVector(vec);
}

void XDRIO :: ReadFloatVector(std::vector<real32>& vec) throw(IOException)
{
  RawIO::ReadInt();
  RawIO::ReadFloatVector(vec);
}


void XDRIO :: ReadDoubleVector(std::vector<real64>& vec) throw(IOException)
{
  RawIO::ReadInt();
  RawIO::ReadDoubleVector(vec);
}

void XDRIO :: WriteMsg(const std::string& s) throw(IOException)
{
  int32 length = s.length();
  int32 pad = 4 - (length % 4);
  int32 bufsize = (pad == 4) ? length : length + pad;

  RawIO::WriteInt(bufsize+4);
  RawIO::WriteMsg(s);
}

void XDRIO :: WriteShortArray(const int16 *buff , const int32 size) throw(IOException)
{
  int32 numitems = size;
  numitems += (size % 2) == 1 ? 1 : 0;

  RawIO::WriteInt(numitems*sizeof(int16)+4);
  RawIO::WriteShortArray(buff, size);
}

void XDRIO :: WriteShortVector(const std::vector<int16>& vec) throw(IOException)
{
  int32 numitems = vec.size();
  numitems += (vec.size() % 2) == 1 ? 1 : 0;

  RawIO::WriteInt(numitems*sizeof(int16)+4);
  RawIO::WriteShortVector(vec);
}

void XDRIO :: WriteUShortArray(const uint16 *buff , const int32 size) throw(IOException)
{
  int32 numitems = size;
  numitems += (size % 2) == 1 ? 1 : 0;

  RawIO::WriteInt(numitems*sizeof(uint16)+4);
  RawIO::WriteUShortArray(buff, size);
}

void XDRIO :: WriteUShortVector(const std::vector<uint16>& vec) throw(IOException)
{
  int32 numitems = vec.size();
  numitems += (vec.size() % 2) == 1 ? 1 : 0;

  RawIO::WriteInt(numitems*sizeof(int16)+4);
  RawIO::WriteUShortVector(vec);
}

void XDRIO :: WriteIntArray(const int32 *buff , const int32 size) throw(IOException)
{
  RawIO::WriteInt(size*sizeof(int32)+4);
  RawIO::WriteIntArray(buff, size);
}

void XDRIO :: WriteIntVector(const std::vector<int32>& vec) throw(IOException)
{
  RawIO::WriteInt(vec.size()*sizeof(int32)+4);
  RawIO::WriteIntVector(vec);
}

void XDRIO :: WriteUIntArray(const uint32 *buff , const int32 size) throw(IOException)
{
  RawIO::WriteInt(size*sizeof(uint32)+4);
  RawIO::WriteUIntArray(buff, size);
}

void XDRIO :: WriteUIntVector(const std::vector<uint32>& vec) throw(IOException)
{
  RawIO::WriteInt(vec.size()*sizeof(uint32)+4);
  RawIO::WriteUIntVector(vec);
}

void XDRIO :: WriteFloatArray(const real32 *buff , const int32 size) throw(IOException)
{
  RawIO::WriteInt(size*sizeof(real32)+4);
  RawIO::WriteFloatArray(buff, size);
}

void XDRIO :: WriteFloatVector(const std::vector<real32>& vec) throw(IOException)
{
  RawIO::WriteInt(vec.size()*sizeof(real32)+4);
  RawIO::WriteFloatVector(vec);
}

void XDRIO :: WriteDoubleArray(const real64 *buff , const int32 size) throw(IOException)
{
  RawIO::WriteInt(size*sizeof(real64)+4);
  RawIO::WriteDoubleArray(buff, size);
}

void XDRIO :: WriteDoubleVector(const std::vector<real64>& vec) throw(IOException)
{
  RawIO::WriteInt(vec.size()*sizeof(real64)+4);
  RawIO::WriteDoubleVector(vec);
}

}
