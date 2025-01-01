#include <fstream>
#include <forward_list>
#include <boost/tokenizer.hpp>
#include <iostream>

#include "Domain/CoefFunction/CoefFunctionFileDataMeas.hh"
#include "Driver/TransientDriver.hh"
#include "DataInOut/ResultHandler.hh"
#include "DataInOut/Logging/LogConfigurator.hh"


using std::string;

DEFINE_LOG(cfdm, "coefFunctionFileDataMeas");

namespace CoupledField
{

CoefFunctionFileDataMeas::CoefFunctionFileDataMeas(Domain* ptDomain, const std::string& pdename,
                                                   PtrParamNode pn, int dim) : CoefFunction()
{
  domain_ = ptDomain;
  pdeName_ = pdename;
  Init(dim);
  filename_ = pn->Get("file")->As<string>();
  missing_  = missing.Parse(pn->Get("missing")->As<string>());

  std::ifstream file(filename_);
  if(!file.good())
    throw Exception("Cannot open fileData file '" + filename_ + "'");

  ReadData(file, dim);
}

void CoefFunctionFileDataMeas::GetVector(Vector<Double>& vec, const LocPointMapped& lpm) {
  //get the actual time step, time
  BaseDriver* driver = domain_->GetDriver();
  TransientDriver * tDriver = dynamic_cast<TransientDriver*>(driver);
  UInt timeStep = 1;
  Double dt = 1.0;
  Double actTime = 1.0;
  if (tDriver != NULL) {
    timeStep = tDriver->GetActStep(pdeName_);
    dt       = tDriver->GetDeltaT();
    actTime  = timeStep*dt;
  }
  //std::cout << "actStep: " << timeStep << "  actTime: " << actTime << std::endl;

  //get correct index
  std::list<StdVector<Double> >::iterator iterData = dataH_.begin();
  std::list<Double>::iterator iterTime = time_.begin();

  for(UInt i=0; i< timeStep-1; i++) {
        ++iterData;
        ++iterTime;
  }

  StdVector<Double> dataVec = *iterData;
  Double time = *iterTime;
  //std::cout << "time: " << time << std::endl;

  //check time
  if ( (std::abs(actTime - time)) / std::abs(time) > 1e-10 )
    EXCEPTION("In CoefFunctionFileDataMeas::GetVector: time from file and transient analysis does not fit!")

  vec.Resize(dataVec.size());
  for ( UInt k=0; k<dataVec.size();k++)
     vec[k] = dataVec[k];

  //std::cout << "DataVec: " << vec << std::endl;
 }


void CoefFunctionFileDataMeas::Init(int dim)
{
  isAnalytic_ = false;
  isComplex_ = false;
  supportDerivative_ = false;
  dimType_ = VECTOR;
  dependType_ = SPACE;
  dim_ = dim;

  // no need to have this static
  missing.SetName("CoefFunctionFileDataMeas");
  missing.Add(THROW_EXCEPTION, "exception");
  missing.Add(WARNING_ZERO, "warning");
  missing.Add(ZERO, "zero");
}

void CoefFunctionFileDataMeas::ReadData(std::istream& input, int dim)
{
  // read file line by line
  string line;
  // we allow almost anything as separation character in any combination and repetitions. 
  boost::char_separator<char> sep(", \t;");

  int line_cnt = 0;
  while (std::getline(input, line))
  {
    line_cnt++;
    // skip empty and comment lines
    if(line.length() < 3 || line.at(0) == '#')
      continue;

    boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
    // tokenizer has no size() :(
    auto iter = tokens.begin();
    StdVector<Double> vec(dim);  //time and 2D / 3D of magnetic intensity)

    try
    {
      // the element number
      if(iter == tokens.end())
        EXCEPTION("cannot parse number token in line " << line_cnt << " in dataFile " << filename_);
      time_.push_back(std::stod(*iter));
      iter++;

      // the data elements
      for(int i = 0; i < dim; i++)
      {
        if(iter == tokens.end())
          EXCEPTION("cannot parse " << dim << " data tokens in line " << line_cnt << " in dataFile " << filename_);
        vec[i] = std::stod(*iter);
        iter++;
      }
    }
    catch(const std::invalid_argument&) {
      EXCEPTION("fail to parse fileData " << filename_ << " line:" << line);
    }
    catch(const std::out_of_range&) {
       EXCEPTION("fail to parse fileData " << filename_ << " line:" << line);
    }

    if(iter != tokens.end())
      EXCEPTION("too much data tokens in line " << line_cnt << " in dataFile " << filename_);
    LOG_DBG2(cfdm) << "RD: read vec=" << ToStringCont(vec);
    dataH_.push_back(vec);
  }

}


} // end of namespace
