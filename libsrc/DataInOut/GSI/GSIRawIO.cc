#include <cstdlib>
#include <cstring>
#include <iostream>

#include "GSIRawIO.hh"
#include "GSIIOException.hh"

GSIRawIO::GSIRawIO(FILE *file_read, FILE *file_write) : GSIBaseIO(file_read, file_write)
{
}

void GSIRawIO :: readMsg(std::string& s) throw(GSIIOException)
{
  int length;
  fread(&length, sizeof(int), 1, file_read);

  if(ferror(file_read)) {
    throw GSIIOException ( "An Error occured while reading the message length" );
  }

  char* Msg=(char*)malloc(length+1);
  memset( Msg, 0, length+1);
	     	      
  fread(Msg, length, 1, file_read);
  if(ferror(file_read)) {
    throw GSIIOException ( "An Error occured while reading the actual message" );
  }

  s = Msg;
  free(Msg);
}

int GSIRawIO :: readInt() throw(GSIIOException)
{
  int buff;

  int n = fread(&buff, sizeof(int), 1, file_read);
  //  std::cout << "int #: " << n << "\n";

  if(ferror(file_read)) {
    throw GSIIOException ( "An Error occured while reading an Integer" );
  }

  return(buff);
}


void GSIRawIO :: readIntVector(std::vector<int>& vec) throw(GSIIOException)
{
  int size = readInt();
  int *buff= new int[size];
  
  vec.resize(size);

  fread(buff, sizeof(int)*size, 1, file_read);

  if(ferror(file_read)) {
    throw GSIIOException ( "An Error occured while reading an integer array" );
  }

  for(int i=0; i<size; i++) {
    vec[i] = buff[i];
  }

  delete[] buff;
}

u_int GSIRawIO :: readUInt() throw(GSIIOException)
{
  u_int buff;

  fread(&buff, sizeof(u_int), 1, file_read);

  if(ferror(file_read)) {
    throw GSIIOException ( "An Error occured while reading an Integer" );
  }

  return(buff);
}


void GSIRawIO :: readUIntVector(std::vector<u_int>& vec) throw(GSIIOException)
{
  int size = readInt();
  u_int *buff= new u_int[size];
  
  vec.resize(size);

  fread(buff, sizeof(u_int)*size, 1, file_read);

  if(ferror(file_read)) {
    throw GSIIOException ( "An Error occured while reading an integer array" );
  }

  for(int i=0; i<size; i++) {
    vec[i] = buff[i];
  }

  delete[] buff;
}


float GSIRawIO :: readFloat() throw(GSIIOException)
{
  float buff;

  fread(&buff, sizeof(float), 1, file_read);

  if(ferror(file_read)) {
    throw GSIIOException ( "An Error occured while reading a Float" );
  }

  return(buff);
}


void GSIRawIO :: readFloatVector(std::vector<float>& vec) throw(GSIIOException)
{
  int size = readInt();
  float *buff= new float[size];
  
  vec.resize(size);

  fread(buff, sizeof(float)*size, 1, file_read);

  if(ferror(file_read)) {
    throw GSIIOException ( "An Error occured while reading a float array" );
  }

  for(int i=0; i<size; i++) {
    vec[i] = buff[i];
  }

  delete[] buff;
}

double GSIRawIO :: readDouble() throw(GSIIOException)
{
  double buff;

  fread(&buff, sizeof(double), 1, file_read);

  if(ferror(file_read)) {
    throw GSIIOException ( "An Error occured while reading a Double" );
  }

  return(buff);
}


void GSIRawIO :: readDoubleVector(std::vector<double>& vec) throw(GSIIOException)
{
  int size = readInt();
  double *buff= new double[size];
  
  vec.resize(size);

  fread(buff, sizeof(double)*size, 1, file_read);

  if(ferror(file_read)) {
    throw GSIIOException ( "An Error occured while reading a double array" );
  }

  for(int i=0; i<size; i++) {
    vec[i] = buff[i];
  }

  delete[] buff;
}

void GSIRawIO :: writeMsg(const std::string& s) throw(GSIIOException)
{
  int length = s.size();
  fwrite(&length, sizeof(int), 1, file_write);
  if(ferror(file_read)) {
    throw GSIIOException ( "An Error occured while writing the message length" );
  }

  fwrite(s.c_str(), length, 1, file_write);
  if(ferror(file_read)) {
    throw GSIIOException ( "An Error occured while writing the actual message" );
  }
  if(fflush(file_write) != 0) {
    throw GSIIOException ( "An Error occured while flushing the message" );
  }
}

void GSIRawIO :: writeInt(const int buff) throw(GSIIOException)
{
  fwrite(&buff, sizeof(int), 1, file_write);
  if(ferror(file_read)) {
    throw GSIIOException ( "An Error occured while writing an integer" );
  }
  if(fflush(file_write) != 0) {
    throw GSIIOException ( "An Error occured while flushing an integer" );
  }
}


void GSIRawIO :: writeIntArray(const int *buff , const int size) throw(GSIIOException)
{
  int s = size;
  fwrite(&s, sizeof(int), 1, file_write);
  if(ferror(file_read)) {
    throw GSIIOException ( "Error while writing the length of integer array" );
  }  
  fwrite(buff, sizeof(int)*size, 1, file_write);
  if(ferror(file_read)) {
    throw GSIIOException ( "An Error occured while writing an integer array" );
  }
  if(fflush(file_write) != 0) {
    throw GSIIOException ( "An Error occured while flushing an integer array" );
  }
}

void GSIRawIO :: writeIntVector(const std::vector<int>& vec) throw(GSIIOException)
{
  int size = vec.size();
  int *buff = new int[size];

  for(int i=0; i<size; i++) {
    buff[i] = vec[i];
  }
  
  writeIntArray(buff, size);

  delete[] buff;
}

void GSIRawIO :: writeUInt(const u_int buff) throw(GSIIOException)
{
  fwrite(&buff, sizeof(u_int), 1, file_write);
  if(ferror(file_read)) {
    throw GSIIOException ( "An Error occured while writing an integer" );
  }
  if(fflush(file_write) != 0) {
    throw GSIIOException ( "An Error occured while flushing an integer" );
  }
}


void GSIRawIO :: writeUIntArray(const u_int *buff , const int size) throw(GSIIOException)
{
  int s = size;
  fwrite(&s, sizeof(int), 1, file_write);
  if(ferror(file_read)) {
    throw GSIIOException ( "Error while writing the length of u_int array" );
  }  
  fwrite(buff, sizeof(u_int)*size, 1, file_write);
  if(ferror(file_read)) {
    throw GSIIOException ( "An Error occured while writing an u_int array" );
  }
  if(fflush(file_write) != 0) {
    throw GSIIOException ( "An Error occured while flushing an u_int array" );
  }
}

void GSIRawIO :: writeUIntVector(const std::vector<u_int>& vec) throw(GSIIOException)
{
  int size = vec.size();
  u_int *buff = new u_int[size];

  for(int i=0; i<size; i++) {
    buff[i] = vec[i];
  }
  
  writeUIntArray(buff, size);

  delete[] buff;
}


void GSIRawIO :: writeFloat(const float buff) throw(GSIIOException)
{
  fwrite(&buff, sizeof(float), 1, file_write);
  if(ferror(file_read)) {
    throw GSIIOException ( "An Error occured while writing a Float" );
  }
  if(fflush(file_write) != 0) {
    throw GSIIOException ( "An Error occured while flushing a Float" );
  }
}

void GSIRawIO :: writeFloatVector(const std::vector<float>& vec) throw(GSIIOException)
{
  int size = vec.size();
  float *buff = new float[size];

  for(int i=0; i<size; i++) {
    buff[i] = vec[i];
  }
  
  writeFloatArray(buff, size);

  delete[] buff;
}

void GSIRawIO :: writeFloatArray(const float *buff , const int size) throw(GSIIOException)
{
  int s = size;
  fwrite(&s, sizeof(int), 1, file_write);
  if(ferror(file_read)) {
    throw GSIIOException ( "Error while writing the length of Float array" );
  }
  fwrite(buff, sizeof(float)*size, 1, file_write);
  if(ferror(file_read)) {
    throw GSIIOException ( "An Error occured while writing a Float array" );
  }
  if(fflush(file_write) != 0) {
    throw GSIIOException ( "An Error occured while flushing a Float array" );
  }
}

void GSIRawIO :: writeDouble(const double buff) throw(GSIIOException)
{
  fwrite(&buff, sizeof(double), 1, file_write);
  if(ferror(file_read)) {
    throw GSIIOException ( "An Error occured while writing a Double" );
  }
  if(fflush(file_write) != 0) {
    throw GSIIOException ( "An Error occured while flushing a double" );
  }
}


void GSIRawIO :: writeDoubleArray(const double *buff , const int size) throw(GSIIOException)
{
  int s = size;
  fwrite(&s, sizeof(int), 1, file_write);
  if(ferror(file_read)) {
    throw GSIIOException ( "Error while writing the length of integer array" );
  }
  fwrite(buff, sizeof(double)*size, 1, file_write);
  if(ferror(file_read)) {
    throw GSIIOException ( "An Error occured while writing a double array" );
  }
  if(fflush(file_write) != 0) {
    throw GSIIOException ( "An Error occured while flushing a double array" );
  }
}

void GSIRawIO :: writeDoubleVector(const std::vector<double>& vec) throw(GSIIOException)
{
  int size = vec.size();
  double *buff = new double[size];

  for(int i=0; i<size; i++) {
    buff[i] = vec[i];
  }
  
  writeDoubleArray(buff, size);

  delete[] buff;
}
