// cvXDRIO.H

#ifndef _CVXDRIO_H_
#define _CVXDRIO_H_

// requires the SUN RPC stuff for platform independend storage
#include <sys/types.h>
#include <rpc/rpc.h>
#include <rpc/xdr.h>

/** A Class for reading data from a file in platform independent format.
 */

class cvXDR
{
public:
  cvXDRInput(const char *filename);
  virtual ~cvXDRInput();

  bool  allOk();

  bool  getInt(int &val);
  bool  getUInt(u_int &val);
  bool  getLong(long &val);
  bool  getULong(u_long &val);
  bool  getShort(short &val);
  bool  getUShort(u_short &val);
  bool  getChar(char &val);
  bool  getUChar(u_char &val);

  bool  getFloat(float &val);
  bool  getDouble(double &val);

  bool  getBool(bool &val);
  bool  getString(char *val,u_int maxlen);

protected:
  FILE  *dataFile;
  XDR    xdrs;
};

/** A Class for writing data to a file in platform independent format
 */

class cvXDROutput
{
public:
  cvXDROutput(const char *filename);
  virtual ~cvXDROutput();

  bool  allOk();

  bool  putInt(int val);
  bool  putUInt(u_int val);
  bool  putLong(long val);
  bool  putULong(u_long val);
  bool  putShort(short val);
  bool  putUShort(u_short val);
  bool  putChar(char val);
  bool  putUChar(u_char val);

  bool  putFloat(float val);
  bool  putDouble(double val);

  bool  putBool(bool val);
  bool  putString(char *val);

protected:
  FILE  *dataFile;
  XDR    xdrs;
};


#endif
