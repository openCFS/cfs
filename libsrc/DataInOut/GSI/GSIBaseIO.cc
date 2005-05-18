#ifdef USE_RCSID
static const char RCSid_GSIBaseIO[] = "$Id$";
#endif

/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.3  2005/05/18 19:26:02  strieben
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

#include "GSIBaseIO.hh"

namespace GridlibSocketInterface
{


BaseIO& operator << (BaseIO& io, const std::string& s ) throw(IOException)
{
  io.WriteMsg(s);

  return io;
}


BaseIO& operator >> (BaseIO& io, std::string& s ) throw(IOException)
{
  io.ReadMsg(s);

  return io;
}

BaseIO& operator << (BaseIO& io, const int16 &i ) throw(IOException)
{
  io.WriteShort(i);

  return io;
}

BaseIO& operator >> (BaseIO& io, int16& i) throw(IOException)
{
  i = io.ReadShort();

  return io;
}

BaseIO& operator << (BaseIO& io, const uint16 &i ) throw(IOException)
{
  io.WriteUShort(i);

  return io;
}

BaseIO& operator >> (BaseIO& io, uint16& i) throw(IOException)
{
  i = io.ReadUShort();

  return io;
}

BaseIO& operator << (BaseIO& io, const int32 &i ) throw(IOException)
{
  io.WriteInt(i);

  return io;
}

BaseIO& operator >> (BaseIO& io, int32& i) throw(IOException)
{
  i = io.ReadInt();

  return io;
}

BaseIO& operator << (BaseIO& io, const uint32 &i ) throw(IOException)
{
  io.WriteUInt(i);

  return io;
}

BaseIO& operator >> (BaseIO& io, uint32& i) throw(IOException)
{
  i = io.ReadUInt();

  return io;
}

BaseIO& operator << (BaseIO& io, const real32 &f ) throw(IOException)
{
  io.WriteFloat(f);

  return io;
}

BaseIO& operator >> (BaseIO& io, real32& f) throw(IOException)
{
  f = io.ReadFloat();

  return io;
}

BaseIO& operator << (BaseIO& io, const real64 &d ) throw(IOException)
{
  io.WriteDouble(d);

  return io;
}

BaseIO& operator >> (BaseIO& io, real64& d) throw(IOException)
{
  d = io.ReadDouble();

  return io;
}

BaseIO& operator << (BaseIO& io, const std::vector<int16> &iv ) throw(IOException)
{
  io.WriteShortVector(iv);

  return io;
}

BaseIO& operator >> (BaseIO& io, std::vector<int16>& iv) throw(IOException)
{
  io.ReadShortVector(iv);

  return io;
}

BaseIO& operator << (BaseIO& io, const std::vector<uint16> &iv ) throw(IOException)
{
  io.WriteUShortVector(iv);

  return io;
}

BaseIO& operator >> (BaseIO& io, std::vector<uint16>& iv) throw(IOException)
{
  io.ReadUShortVector(iv);

  return io;
}

BaseIO& operator << (BaseIO& io, const std::vector<int32> &iv ) throw(IOException)
{
  io.WriteIntVector(iv);

  return io;
}

BaseIO& operator >> (BaseIO& io, std::vector<int32>& iv) throw(IOException)
{
  io.ReadIntVector(iv);

  return io;
}

BaseIO& operator << (BaseIO& io, const std::vector<uint32> &iv ) throw(IOException)
{
  io.WriteUIntVector(iv);

  return io;
}

BaseIO& operator >> (BaseIO& io, std::vector<uint32>& iv) throw(IOException)
{
  io.ReadUIntVector(iv);

  return io;
}

BaseIO& operator << (BaseIO& io, const std::vector<real32> &fv ) throw(IOException)
{
  io.WriteFloatVector(fv);

  return io;
}

BaseIO& operator >> (BaseIO& io, std::vector<real32>& fv) throw(IOException)
{
  io.ReadFloatVector(fv);

  return io;
}

BaseIO& operator << (BaseIO& io, const std::vector<real64> &dv ) throw(IOException)
{
  io.WriteDoubleVector(dv);

  return io;
}

BaseIO& operator >> (BaseIO& io, std::vector<real64>& dv) throw(IOException)
{
  io.ReadDoubleVector(dv);

  return io;
}

 
}
