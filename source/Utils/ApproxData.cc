#include <iostream>
#include <fstream>

#include "Utils/tools.hh"
//#include "DataInOut/simInput.hh"
#include "ApproxData.hh"
#include "DataInOut/Logging/LogConfigurator.hh"  

namespace CoupledField
{ 

DEFINE_LOG(approxdata, "approxdata")

ApproxData::ApproxData(std::string nlFileName, MaterialType matType, UInt numIndep)
{
  numIndepend_ = numIndep;
  matType_ = matType;
  nlFileName_ = nlFileName;
  factor_ = 1.;
  if (numIndepend_ == 1) {
    ReadNlinFunc(nlFileName);
    PerformChecksOnInputData(nlFileName);
  } else if (numIndepend_ == 2) {
    ReadNlinFuncTwoIndep(nlFileName);
  } else if (numIndepend_ == 3) {
    ReadNlinFuncThreeIndep(nlFileName);
  }
  else
    EXCEPTION("Cannot handle number of independent variables: " << numIndepend_);

}


void ApproxData::ReadNlinFunc(std::string fncName)  {

  std::ifstream datafile;

  datafile.open(fncName.c_str());
  if ( !datafile ) {
    std::string str = "Failed to open file '" + fncName +
        "' suspected to contain nonlinear data";
    EXCEPTION( str );
  }

  datafile.clear(); // clear flags

  // we don't trust .eof() =)
  datafile.seekg(0,std::ios::end);
  std::string::size_type pos = 0, line_start_pos = 0,
      pos_end = datafile.tellg();

  datafile.seekg(0,std::ios::beg); // start from the beginning
  std::string     buf;
  double          xval,yval;
  std::vector<double> xx, yy;

  while(pos <= pos_end)
  {
    line_start_pos = datafile.tellg();
    std::getline(datafile,buf);

    // big choice of signs for comment's
    if ((buf[0] != '#' || buf[0] != '%' || buf[0] != '!') && buf.size() > 0)
    {
      datafile.seekg(line_start_pos); // rewind

      datafile >> xval >> yval;
      yy.push_back(yval);
      xx.push_back(xval);
      datafile.ignore(100,'\n');

    }

    pos = datafile.tellg();  // and, where we are ?
  }

  datafile.close();
  numMeas_ = xx.size();

  x_.Resize(numMeas_);
  y_.Resize(numMeas_);

  for ( UInt k=0; k<numMeas_; k++ ) {
    x_[k] = xx[k];
    y_[k] = yy[k];
  }
}

void ApproxData::ReadNlinFuncThreeIndep(std::string fncName)  {
  // general idea:
  // filename given is a list of <value, filename> pairs
  std::ifstream datafile;

  datafile.open(fncName.c_str());
  if ( !datafile ) {
    std::string str = "Failed to open file '" + fncName +
        "' suspected to contain nonlinear data";
    EXCEPTION( str );
  }

  datafile.clear(); // clear flags

  // we don't trust .eof() =)
  datafile.seekg(0,std::ios::end);
  std::string::size_type pos = 0, line_start_pos = 0,
      pos_end = datafile.tellg();

  datafile.seekg(0,std::ios::beg); // start from the beginning
  std::string     buf;
  std::vector<double> x;
  std::vector<std::string> filenames;
  double xt;
  std::string iFilename;
  while(pos <= pos_end)
  {
    line_start_pos = datafile.tellg();
    std::getline(datafile,buf);

    if ((buf[0] != '#' || buf[0] != '%' || buf[0] != '!') && buf.size() > 0)
    {
      datafile.seekg(line_start_pos); // rewind
      datafile >> xt >> iFilename;
      x.push_back(xt);
      filenames.push_back(iFilename);
      datafile.ignore(100,'\n');
    }
    pos = datafile.tellg();  // and, where we are ?
  }

  datafile.close();
  numMeas_ = x.size();
  if ( x.size() != filenames.size() || x.size() < 2) {
    EXCEPTION("Error in slices data file: You need to provide a list of <double, string> tuples, each providing the data for a 2d slice. You need to provide at least two slices.");
  }


  LOG_DBG(approxdata) << "Read 3d slices file '" << fncName << "' with " << numMeas_ << " slice files ";

  x_.Resize(numMeas_);
  slicesFiles_.Resize(numMeas_);
  for (UInt k=0; k < numMeas_; k++) {

    x_[k] = x[k];
    slicesFiles_[k] = filenames[k];
  }

  for (UInt l=1; l<numMeas_; l++) {
    if ( x_[l-1] >= x_[l]){
      EXCEPTION("coordinate x must increase monotonically");
    }
  }
  LOG_DBG(approxdata) << "x Vector: " << x_.ToString();
  LOG_DBG(approxdata) << "slices filenames: \n" << slicesFiles_.ToString();


}

void ApproxData::findBracketIndices(const double &x, const Vector<double> & axis, UInt & klo, UInt & khi, double &diff) const {
  // instead of findBracketIndices one can use
  // khi = std::upper_bound(axis.GetPointer(), axis.GetPointer() + axis.GetSize() - 1, x) - axis.GetPointer();

  const UInt kend = axis.GetSize() - 1;
  UInt k;
  klo=0;
  khi=kend;
  // We will find the right place in the table by means of bisection.
  //  klo and khi bracket the input value of xEntry
  while (khi-klo > 1) {
    k=(khi+klo) >> 1; // binary right shift
    if (axis[k] > x)
      khi=k;
    else
      klo=k;
  }
  // size of x interval
  diff = axis[khi] - axis[klo];

  // The x-values must be distinct!
  if (diff == 0.0) {
    EXCEPTION("You cannot have two equal x values!" );
  }

}


void ApproxData::ReadNlinFuncTwoIndep(std::string fncName)  {

  std::ifstream datafile;

  datafile.open(fncName.c_str());
  if ( !datafile ) {
    std::string str = "Failed to open file '" + fncName +
        "' suspected to contain nonlinear data";
    EXCEPTION( str );
  }

  datafile.clear(); // clear flags

  // we don't trust .eof() =)
  datafile.seekg(0,std::ios::end);
  std::string::size_type pos = 0, line_start_pos = 0,
      pos_end = datafile.tellg();

  datafile.seekg(0,std::ios::beg); // start from the beginning
  std::string     buf;
  double          x0val,x1val, yval;
  std::vector<double> xx0, xx1;
  std::vector<double> yy;
  std::vector<UInt> blockSizes;
  UInt nblocks;
  UInt nMeas = 0;

  while(pos <= pos_end)
  {
    line_start_pos = datafile.tellg();
    std::getline(datafile,buf);
    // Data file must be in gnuplot syntax:
    // x0 y0 z0
    // x0 y1 z1
    // x0 y2 z3
    //
    // x1 y0 z4
    // x1 y1 z5
    // x1 y2 z6
    //
    // ...
    // big choice of signs for comment's
    if ((buf[0] != '#' || buf[0] != '%' || buf[0] != '!') && buf.size() > 0)
    {
      datafile.seekg(line_start_pos); // rewind

      datafile >> x0val >> x1val >> yval;
      xx0.push_back(x0val);
      xx1.push_back(x1val);
      yy.push_back(yval);
      datafile.ignore(100,'\n');
      nMeas += 1;
    }
    else if (buf.size() == 0) {
      // advance block
      if (nMeas>0)
        blockSizes.push_back(nMeas);
      nMeas=0;
    }

    pos = datafile.tellg();  // and, where we are ?
  }

  datafile.close();

  UInt firstBlockSize = blockSizes[0];
  for (UInt k=1; k< blockSizes.size(); k ++ ){
    // blocks have to be the same size
    if (blockSizes[k] != firstBlockSize)
      EXCEPTION("each block size must be the same in the data file: "
          << fncName);
    // the x1 coordinates must be the same for each block
    for (UInt l=0; l< firstBlockSize; l++) {
      if (xx1[l] != xx1[k*firstBlockSize + l])
        EXCEPTION("coordinate x1 must be the same in each block of the data file: "
            << fncName );
    }
  }
  for (UInt k=0; k<blockSizes.size(); k++) {
    for (UInt l=0; l<firstBlockSize; l++) {
      if (xx0[k*firstBlockSize] != xx0[k*firstBlockSize + l])
        EXCEPTION("coordinate x0 must be the same throughout a block of the data file: "
            << fncName);
    }
  }

  numMeas_ = xx0.size(); // FIXME: do I need it somewhere?
  nblocks = blockSizes.size();
  LOG_DBG(approxdata) << "Read 2d file '" << fncName << "' with " << nblocks <<
      " blocks with " << firstBlockSize << " measurements each.";

  x_.Resize(nblocks);
  z_.Resize(nblocks);
  x1_.Resize(firstBlockSize);

  for (UInt l=0; l<firstBlockSize; l++) {
    x1_[l] = xx1[l];
  }
  for (UInt l=1; l<firstBlockSize; l++) {
    if ( x1_[l-1] >= x1_[l]){
      LOG_DBG(approxdata) << "x1 Vector: " << x1_.ToString();
      EXCEPTION("coordinate x1 must increase monotonically. In file : "
          << fncName);
    }
  }

  for ( UInt k=0; k<nblocks; k++ ) {
    z_[k].Resize(firstBlockSize);
    x_[k] = xx0[k*firstBlockSize];
    for (UInt l=0; l<firstBlockSize; l++) {
      z_[k][l] = yy[k*firstBlockSize + l];
    }
  }
  for (UInt l=1; l<nblocks; l++) {
    if ( x_[l-1] >= x_[l]){
      LOG_DBG(approxdata) << "x0 Vector: " << x_.ToString();
      EXCEPTION("coordinate x0 must increase monotonically. In file: "
          << fncName);
    }

  }
  LOG_DBG(approxdata) << "x0 Vector: " << x_.ToString();
  LOG_DBG(approxdata) << "x1 Vector: " << x1_.ToString();
  LOG_DBG(approxdata) << "z   Array: \n" << z_.ToString();


}


void ApproxData::PerformChecksOnInputData( std::string nlFileName ) {

  if  ( matType_ == MAG_PERMEABILITY_SCALAR ) {
    //all x=H have to be greater than zero
    for ( UInt k=0; k<numMeas_; k++ ) {
      if ( x_[k] < 0.0 ) {
        std::string str = "There are values for H in file '" + nlFileName +
            "', that are smaller than zero!";
        EXCEPTION( str );
      }
    }
    //all y=B have to be greater than zero
    for ( UInt k=0; k<numMeas_; k++ ) {
      if ( y_[k] < 0.0 ) {
        std::string str = "There are values for B in file '" + nlFileName +
            "', that are smaller than zero!";
        EXCEPTION( str );
      }
    }

    // starting point should be (0,0)
    if ( x_[0] != 0.0 ) {
      WARN("First measured point in BH curve read from file '"
          << nlFileName << "' is not zero. We will add such a point!");

      Vector<double> xx = x_;
      Vector<double> yy = y_;
      numMeas_ +=1;

      x_.Resize(numMeas_);
      y_.Resize(numMeas_);
      x_[0] = 0.0;
      y_[0] = 0.0;

      for ( UInt k=1; k<numMeas_; k++ ) {
        x_[k] = xx[k-1];
        y_[k] = yy[k-1];
      }
    }

    //now check, if we have monoton increasing x- and y-values
    for ( UInt k=0; k<numMeas_-1; k++ ) {
      double epsX = 1e-2;
      double epsY = 1e-3;
      if ( x_[k+1] - x_[k] < epsX ) {
        std::string str = "The H-values in file '" + nlFileName +
            "', are not monoton increasing (epsH = 1e-2)! ";
        EXCEPTION( str );
      }
      if ( y_[k+1] - y_[k] < epsY ) {        
        std::string str = "The B-values in file '" + nlFileName +
            "', are not monoton increasing (epsB = 1e-3)!";
        EXCEPTION( str );
      }
    }
  }
}

}
