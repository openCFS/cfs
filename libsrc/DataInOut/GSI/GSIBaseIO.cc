#include "GSIBaseIO.hh"

GSIBaseIO& operator << (GSIBaseIO& io, const std::string& s ) throw(GSIIOException)
{
  io.writeMsg(s);

  return io;
}


GSIBaseIO& operator >> (GSIBaseIO& io, std::string& s ) throw(GSIIOException)
{
  io.readMsg(s);

  return io;
}

GSIBaseIO& operator << (GSIBaseIO& io, const int &i ) throw(GSIIOException)
{
  io.writeInt(i);

  return io;
}

GSIBaseIO& operator >> (GSIBaseIO& io, int& i) throw(GSIIOException)
{
  i = io.readInt();

  return io;
}

GSIBaseIO& operator << (GSIBaseIO& io, const float &f ) throw(GSIIOException)
{
  io.writeFloat(f);

  return io;
}

GSIBaseIO& operator >> (GSIBaseIO& io, float& f) throw(GSIIOException)
{
  f = io.readFloat();

  return io;
}

GSIBaseIO& operator << (GSIBaseIO& io, const double &d ) throw(GSIIOException)
{
  io.writeDouble(d);

  return io;
}

GSIBaseIO& operator >> (GSIBaseIO& io, double& d) throw(GSIIOException)
{
  d = io.readDouble();

  return io;
}

GSIBaseIO& operator << (GSIBaseIO& io, const std::vector<int> &iv ) throw(GSIIOException)
{
  io.writeIntVector(iv);

  return io;
}

GSIBaseIO& operator >> (GSIBaseIO& io, std::vector<int>& iv) throw(GSIIOException)
{
  io.readIntVector(iv);

  return io;
}

GSIBaseIO& operator << (GSIBaseIO& io, const std::vector<u_int> &iv ) throw(GSIIOException)
{
  io.writeUIntVector(iv);

  return io;
}

GSIBaseIO& operator >> (GSIBaseIO& io, std::vector<u_int>& iv) throw(GSIIOException)
{
  io.readUIntVector(iv);

  return io;
}

GSIBaseIO& operator << (GSIBaseIO& io, const std::vector<float> &fv ) throw(GSIIOException)
{
  io.writeFloatVector(fv);

  return io;
}

GSIBaseIO& operator >> (GSIBaseIO& io, std::vector<float>& fv) throw(GSIIOException)
{
  io.readFloatVector(fv);

  return io;
}

GSIBaseIO& operator << (GSIBaseIO& io, const std::vector<double> &dv ) throw(GSIIOException)
{
  io.writeDoubleVector(dv);

  return io;
}

GSIBaseIO& operator >> (GSIBaseIO& io, std::vector<double>& dv) throw(GSIIOException)
{
  io.readDoubleVector(dv);

  return io;
}
