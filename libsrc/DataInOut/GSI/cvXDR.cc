// cvXDRIO.C

#include <string.h>
#include "cvXDRIO.H"

// ----- Input -----

cvXDRInput::cvXDRInput(const char *filename)
{
  dataFile = fopen(filename,"r");
  if(dataFile!=NULL) {
     xdrstdio_create( &xdrs, dataFile, XDR_DECODE );
  }
}

cvXDRInput::~cvXDRInput()
{
  if(dataFile!=NULL)
    fclose(dataFile);
}

bool cvXDRInput::allOk()
{
  return (dataFile!=NULL) && (ferror(dataFile)==0) && (feof(dataFile)==0);
}

bool cvXDRInput::getInt(int &val)
{
  return xdr_int(&xdrs,&val);
}

bool cvXDRInput::getUInt(u_int &val)
{
  return xdr_u_int(&xdrs,&val);
}

bool cvXDRInput::getLong(long &val)
{
  return xdr_long(&xdrs,&val);
}

bool cvXDRInput::getULong(u_long &val)
{
  return xdr_u_long(&xdrs,&val);
}

bool cvXDRInput::getShort(short &val)
{
  return xdr_short(&xdrs,&val);
}

bool cvXDRInput::getUShort(u_short &val)
{
  return xdr_u_short(&xdrs,&val);
}

bool cvXDRInput::getChar(char &val)
{
  return xdr_char(&xdrs,&val);
}

bool cvXDRInput::getUChar(u_char &val)
{
  return xdr_u_char(&xdrs,&val);
}

bool cvXDRInput::getFloat(float &val)
{
  return xdr_float(&xdrs,&val);
}

bool cvXDRInput::getDouble(double &val)
{
  return xdr_double(&xdrs,&val);
}

bool cvXDRInput::getBool(bool &val)
{
  bool_t mybool;
  bool result;

  result = xdr_bool(&xdrs,&mybool);
  val = mybool;
  return result;
}

bool cvXDRInput::getString(char *val,u_int maxlen)
{
  return xdr_string(&xdrs,&val,maxlen);
}

// ----- Output -----

cvXDROutput::cvXDROutput(const char *filename)
{
  dataFile = fopen(filename,"w");
  if(dataFile!=NULL) {
     xdrstdio_create( &xdrs, dataFile, XDR_ENCODE );
  }
}

cvXDROutput::~cvXDROutput()
{
  if(dataFile!=NULL)
    fclose(dataFile);
}

bool cvXDROutput::allOk()
{
  return (dataFile!=NULL) && (ferror(dataFile)==0) && (feof(dataFile)==0);
}

bool cvXDROutput::putInt(int val)
{
  return xdr_int(&xdrs,&val);
}

bool cvXDROutput::putUInt(u_int val)
{
  return xdr_u_int(&xdrs,&val);
}

bool cvXDROutput::putLong(long val)
{
  return xdr_long(&xdrs,&val);
}

bool cvXDROutput::putULong(u_long val)
{
  return xdr_u_long(&xdrs,&val);
}

bool cvXDROutput::putShort(short val)
{
  return xdr_short(&xdrs,&val);
}

bool cvXDROutput::putUShort(u_short val)
{
  return xdr_u_short(&xdrs,&val);
}

bool cvXDROutput::putChar(char val)
{
  return xdr_char(&xdrs,&val);
}

bool cvXDROutput::putUChar(u_char val)
{
  return xdr_u_char(&xdrs,&val);
}

bool cvXDROutput::putFloat(float val)
{
  return xdr_float(&xdrs,&val);
}

bool cvXDROutput::putDouble(double val)
{
  return xdr_double(&xdrs,&val);
}

bool cvXDROutput::putBool(bool val)
{
  bool_t mybool = val;

  return xdr_bool(&xdrs,&mybool);
}

bool cvXDROutput::putString(char *val)
{
  return xdr_string(&xdrs,&val,strlen(val)+1);
}
