#include <fstream>
#include <forward_list>
#include <boost/tokenizer.hpp>
#include <iostream>

#include "Domain/CoefFunction/CoefFunctionFileData.hh"
#include "Utils/tools.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

using std::string;

DEFINE_LOG(cfd, "coefFileData");

namespace CoupledField
{

CoefFunctionFileData::CoefFunctionFileData(PtrParamNode pn, int dim) : CoefFunction()
{
  Init(dim);
  filename_ = pn->Get("file")->As<string>();
  missing_  = missing.Parse(pn->Get("missing")->As<string>());

  std::ifstream file(filename_);
  if(!file.good())
    throw Exception("Cannot open fileData file '" + filename_ + "'");

  ReadData(file, dim);
}

CoefFunctionFileData::CoefFunctionFileData(std::istream& input, int dim) : CoefFunction()
{
  Init(dim);
  filename_ = "std::istream";

  ReadData(input, dim);
}


void CoefFunctionFileData::Init(int dim)
{
  isAnalytic_ = false;
  isComplex_ = false;
  supportDerivative_ = false;
  dimType_ = dim == 1 ? SCALAR : VECTOR;
  dependType_ = SPACE;
  node_guess_.Set(0);

  // no need to have this static
  missing.SetName("CoefFunctionFileData");
  missing.Add(THROW_EXCEPTION, "exception");
  missing.Add(WARNING_ZERO, "warning");
  missing.Add(ZERO, "zero");
}

int CoefFunctionFileData::GetIndex(const LocPointMapped& lpm)
{
  assert(lpm.lp.number != LocPoint::NOT_SET);

  unsigned int number = lpm.lp.number;
  int idx = node_.Find(number, node_guess_.Mine());

  LOG_DBG2(cfd) << "GI: number=" << number << " idx=" << idx << " rows=" << data_.GetNumRows();

  if(idx == -1)
  {
    if(missing_ == THROW_EXCEPTION)
      throw Exception("fileData: element number " + std::to_string(number) + " not given in file '" + filename_ + "' and 'missing' set to 'exception'");
    if(missing_ == WARNING_ZERO && !warning_issued_) {
      PtrParamNode pn = domain->GetInfoRoot()->Get(ParamNode::PROCESS);
      pn->SetWarning("fileData: at least element number " + std::to_string(number) + " not given in file '" + filename_ + "' and set to zero");
      warning_issued_ = true;
    }
  }

  assert(idx < (int) data_.GetNumRows());
  return idx;
}


void CoefFunctionFileData::ReadData(std::istream& input, int dim)
{

  // we don't know the size of the data file in advance. Hence read to lists and copy later for compact storage to StdVector and Matrix
  std::list<int> num;
  std::list<std::vector<double> > dat;

  // read file line by line
  string line;
  // we allow almost anything as separation character in any combination and repetitions. see also testbed.cc StringParse
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
    std::vector<double> vec(dim);

    try
    {
      // the element number
      if(iter == tokens.end())
        EXCEPTION("cannot parse number token in line " << line_cnt << " in dataFile " << filename_);
      num.push_back(std::stoi(*iter));
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
    LOG_DBG2(cfd) << "RD: read vec=" << ToStringCont(vec);
    dat.push_back(vec);
  }

  // copy list stuff to own containers
  node_.Import(num.begin(), num.end(), num.size());
  data_.Resize(num.size(),dim);
  int idx = 0;
  for(const auto& vec : dat) {
    for(unsigned j = 0; j < (unsigned int) dim; j++) {
      LOG_DBG2(cfd) << "RD: idx=" << idx << " j=" << j << " v=" << vec.at(j);
      data_.SetEntry(idx, j, vec.at(j));
    }
    idx++;
  }


  assert(idx == (int) node_.GetSize());
}


} // end of namespace
