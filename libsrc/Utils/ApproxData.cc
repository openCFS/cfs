#include <iostream>
#include <fstream>

#include "Utils/tools.hh"
#include "ApproxData.hh"

namespace CoupledField
{ 

ApproxData :: ApproxData(std::string nlFileName)
{
  ENTER_FCN("ApproxData::ApproxData" );
  ReadNlinFunc(nlFileName);
}

ApproxData :: ~ApproxData()
{
  delete [] x;
  delete [] y;
}


void ApproxData ::ReadNlinFunc(std::string fncName)
{
  ENTER_FCN("ApproxData::ReadNlinFunc" );
  
  std::ifstream datafile;
  
  datafile.open(fncName.c_str());
  if (!datafile)  
    Error("Can't open file with nonlinear data");
  
  datafile.clear(); // clear flags
  
  // we don't trust .eof() =)
  datafile.seekg(0,std::ios::end);
  std::string::size_type pos = 0, pos_end = datafile.tellg();
  
  datafile.seekg(0,std::ios::beg); // start from the beginning
  std::string     buf;
  Double          xval,yval;
  std::vector<Double> xx, yy;

  while(pos <= pos_end)
    {	  
      std::getline(datafile,buf);
      
      // big choice of signs for comment's
      if ((buf[0] != '#' || buf[0] != '%' || buf[0] != '!') && buf.size() > 0) 
	{
	  datafile.seekg(- buf.size() - 1,std::ios::cur); // rewind
	  
	  datafile >> xval >> yval;		     
	  yy.push_back(yval);
	  xx.push_back(xval);
	  datafile.ignore(100,'\n');
	  
	}
      
      pos = datafile.tellg();  // and, where we are ?    
    }
  
  datafile.close();
  nummeas = xx.size();

  x = new double[nummeas];
  y = new double[nummeas];

  for (Integer k=0; k<nummeas; k++) {
    x[k] = xx[k];
    y[k] = yy[k];
  }

}

}
