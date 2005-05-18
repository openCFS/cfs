#include "GSIBaseIO.hh"

namespace GridlibSocketInterface
{


BaseIO& operator << (BaseIO& io, const std::string& s ) throw(IOException)
{
  io.writeMsg(s);

  return io;
}


BaseIO& operator >> (BaseIO& io, std::string& s ) throw(IOException)
{
  io.readMsg(s);

  return io;
}

BaseIO& operator << (BaseIO& io, const int &i ) throw(IOException)
{
  io.writeInt(i);

  return io;
}

BaseIO& operator >> (BaseIO& io, int& i) throw(IOException)
{
  i = io.readInt();

  return io;
}

BaseIO& operator << (BaseIO& io, const float &f ) throw(IOException)
{
  io.writeFloat(f);

  return io;
}

BaseIO& operator >> (BaseIO& io, float& f) throw(IOException)
{
  f = io.readFloat();

  return io;
}

BaseIO& operator << (BaseIO& io, const double &d ) throw(IOException)
{
  io.writeDouble(d);

  return io;
}

BaseIO& operator >> (BaseIO& io, double& d) throw(IOException)
{
  d = io.readDouble();

  return io;
}

BaseIO& operator << (BaseIO& io, const std::vector<int> &iv ) throw(IOException)
{
  io.writeIntVector(iv);

  return io;
}

BaseIO& operator >> (BaseIO& io, std::vector<int>& iv) throw(IOException)
{
  io.readIntVector(iv);

  return io;
}

BaseIO& operator << (BaseIO& io, const std::vector<u_int> &iv ) throw(IOException)
{
  io.writeUIntVector(iv);

  return io;
}

BaseIO& operator >> (BaseIO& io, std::vector<u_int>& iv) throw(IOException)
{
  io.readUIntVector(iv);

  return io;
}

BaseIO& operator << (BaseIO& io, const std::vector<float> &fv ) throw(IOException)
{
  io.writeFloatVector(fv);

  return io;
}

BaseIO& operator >> (BaseIO& io, std::vector<float>& fv) throw(IOException)
{
  io.readFloatVector(fv);

  return io;
}

BaseIO& operator << (BaseIO& io, const std::vector<double> &dv ) throw(IOException)
{
  io.writeDoubleVector(dv);

  return io;
}

BaseIO& operator >> (BaseIO& io, std::vector<double>& dv) throw(IOException)
{
  io.readDoubleVector(dv);

  return io;
}

 
}
