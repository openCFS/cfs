#include "SimOutputStreaming.hh"

#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include <boost/asio.hpp>
#include <boost/bind.hpp>

using boost::asio::ip::tcp;

using namespace CoupledField;
using std::string;

SimOutputStreaming::SimOutputStreaming(PtrParamNode outputNode) :
  SimOutput("", outputNode),
  silent_(false)
{
  formatName_ = "streaming";
  capabilities_.insert(MESH_RESULTS);

  host_ = outputNode->Get("host")->As<string>();
  port_ = outputNode->Get("port")->As<string>();
  path_ = outputNode->Get("path")->As<string>();
  send_mesh_ = outputNode->Get("sendMesh")->As<bool>();
  compressed_ = outputNode->Get("compressed")->As<bool>();

  if(outputNode->Has("silent"))
    silent_ = outputNode->Get("silent")->As<bool>();

  info_root = info->Get("streaming")->Get(ParamNode::PROCESS); // TODO!
}

SimOutputStreaming::~SimOutputStreaming()
{
}

void SimOutputStreaming::Init(Grid * ptGrid, bool printGridOnly)
{
}

//! Begin multisequence step
void SimOutputStreaming::BeginMultiSequenceStep(UInt step, BasePDE::AnalysisType type, UInt numSteps)
{
  actMSStep_ = step;
}

//! Register result (within one multisequence step)
void SimOutputStreaming::RegisterResult(shared_ptr<BaseResult> br, UInt saveBegin, UInt saveInc, UInt saveEnd, bool isHistory)
{
   shared_ptr<ResultInfo> ri = br->GetResultInfo();
   PtrParamNode in = info_root->GetByVal("sequence", "step", boost::lexical_cast<string>(actMSStep_));
   // create new
   in = in->Get("result", ParamNode::APPEND);
   br->SetInfoNode(in);
   in->Get("data")->SetValue(ri->resultName);
   in->Get("location")->SetValue(br->GetEntityList()->GetName());
   string loc_type;
   ResultInfo::Enum2String(ri->definedOn, loc_type);
   in->Get("definedOn")->SetValue(loc_type);
}

void SimOutputStreaming::BeginStep( UInt stepNum, Double stepVal)
{
  actStep_ = stepNum;
  actStepVal_ = stepVal;
}

void SimOutputStreaming::AddResult( shared_ptr<BaseResult> sol)
{
  results_.Push_back(sol);
}

void SimOutputStreaming::FinishStep()
{
  boost::asio::io_service io_service;
  Client client(io_service, host_, port_, path_, this);
  io_service.run();
  results_.Clear();
}

UInt SimOutputStreaming::GetContentLength()
{
  UInt v = 40;
  for(UInt r = 0, s = results_.GetSize(); r < s; ++r){
    if(results_[r]->GetResultInfo()->resultName != "physicalPseudoDensity")
      continue;
    Vector<Double>& resultVec = dynamic_cast<Result<Double>&>(*results_[r]).GetVector();
    v += resultVec.GetSize() * 7; 
  }
  return(v + 2);
}


void SimOutputStreaming::Transmit(std::ostream& out)
{
  // create a new Param Node
  // ParamNode content = PtrParamNode(new ParamNode(ParamNode::INSERT, ParamNode::ELEMENT ));
  ParamNode content(ParamNode::INSERT);
  content.SetName("cfsStreaming");

  // add the complete info.xml treee
  content.Get("info")->SetValue(info, true); // use own name and make sure we are not seed as stupid boost::any

  // add mesh
  if(send_mesh_) domain->GetGrid()->ExportGrid(content.Get("mesh"));

  // now the actual results
  PtrParamNode results = content.Get("results");
  for(unsigned int r = 0; r < results_.GetSize(); r++)
  {
    shared_ptr<BaseResult> sol = results_[r];
    shared_ptr<ResultInfo> resInfo = sol->GetResultInfo();

    PtrParamNode rpn = results->Get("result", ParamNode::APPEND);
    std::string tmp;
    resInfo->Enum2String(resInfo->definedOn, tmp); // how ugly! :(
    rpn->Get("solution")->SetValue(tmp);

    rpn->Get("name")->SetValue(resInfo->resultName);
    rpn->Get("region")->SetValue(sol->GetEntityList()->GetName());

    unsigned int dofs = resInfo->dofNames.GetSize();
    rpn->Get("dofs")->SetValue(dofs);
    for(unsigned int d = 0; d < dofs; d++)
      rpn->Get("dof_" + boost::lexical_cast<std::string>(d))->SetValue(resInfo->dofNames[d]);

    rpn->Get("unit")->SetValue(resInfo->unit);

    // TODO now only real valued!
    // the stuff goes into the data element and we set the elements directly
    PtrParamNode data = rpn->Get("data");
    ParamNodeList& list = data->GetChildren();
    // the result vector is a set of dofs
    Vector<double>& resultVec = dynamic_cast<Result<double>&>(*sol).GetVector();
    // the actual item size
    unsigned int items = resultVec.GetSize() / dofs;

    list.Resize(items);
    for(unsigned int e = 0; e < items; e++)
    {
      PtrParamNode el = data->SetNewChild("item", e);

      for(unsigned int d = 0; d < dofs; d++)
        el->Get("v_" + boost::lexical_cast<std::string>(d))->SetValue(resultVec[dofs * e + d]);

      // add the optional index only at the end to allow save attribute retrieval by index
      if(!compressed_)
        el->Get("idx")->SetValue(e);
    }
  }


  // write the stuff
  content.ToXML(out, compressed_ ? -99 : 0, true); // adjust element type!

}

/*
void SimOutputStreaming::Transmit(std::ostream& out)
{
  out << "Iteration: "; 
  out.width(8);
  out << actStep_ << std::endl;
  for(UInt r = 0; r < results_.GetSize(); r++)
  {
    shared_ptr<BaseResult> sol = results_[r];

    //std::string regionName = sol->GetEntityList()->GetName();
    shared_ptr<ResultInfo> resInfo = sol->GetResultInfo();
    std::string resultName = resInfo->resultName;
    //std::cout << "result: " << resultName << std::endl;
    // UInt numDOFs = resInfo->dofNames.GetSize();

    if(resultName != "physicalPseudoDensity") continue;
    Vector<Double>& resultVec = dynamic_cast<Result<Double>&>(*sol).GetVector();
    out << "Densities: "; 
    out.width(8);
    out << resultVec.GetSize() << std::endl;

    out.setf(std::ios::fixed);
    out.precision(4); // we need to be below 16K for the current reader or reading might break
    for(UInt d = 0; d < resultVec.GetSize(); d++)
    {
      out << resultVec[d] << std::endl;
    }
  }
  out << std::endl << std::endl;
}
*/

SimOutputStreaming::Client::Client(boost::asio::io_service& io_service,
      const std::string& server, const std::string& port, const std::string& path, SimOutputStreaming* base)
    : resolver_(io_service),
      socket_(io_service)
{
    // it is more robust known the content length a priori, therefore we use a memory stream
    std::stringstream mem;
    // fill with all our data
    base->Transmit(mem);

    std::ostream request_stream(&request_);
    request_stream << "POST " << path << " HTTP/1.1\r\n";
    request_stream << "Host: " << server << "\r\n";
    request_stream << "Accept: */*\r\n";
    request_stream << "Connection: close\r\n";
    request_stream << "Content-Length: " << mem.tellp() << "\r\n";
    request_stream << "Content-Type: text/plain\r\n\r\n";

    // copy data
    request_stream << mem.rdbuf() << "\r\n\r\n";

    if(!base->silent_)
      std::cout << "try to connect " << server << " port " << port << " to transmit " << mem.tellp() << " bytes" << std::endl;

    // Start an asynchronous resolve to translate the server and service names
    // into a list of endpoints.
    tcp::resolver::query query(server, port); // we use http anyway
    resolver_.async_resolve(query,
        boost::bind(&Client::handle_resolve, this,
          boost::asio::placeholders::error,
          boost::asio::placeholders::iterator));
}

void SimOutputStreaming::Client::handle_resolve(const boost::system::error_code& err,
      tcp::resolver::iterator endpoint_iterator)
{
    if (!err)
    {
      // Attempt a connection to the first endpoint in the list. Each endpoint
      // will be tried until we successfully establish a connection.
      tcp::endpoint endpoint = *endpoint_iterator;
      socket_.async_connect(endpoint,
          boost::bind(&Client::handle_connect, this,
            boost::asio::placeholders::error, ++endpoint_iterator));
    }
    else
    {
      std::cout << "Error: " << err.message() << "\n";
    }
}

void SimOutputStreaming::Client::handle_connect(const boost::system::error_code& err,
      tcp::resolver::iterator endpoint_iterator)
{
    if (!err)
    {
      // The connection was successful. Send the request.
      boost::asio::async_write(socket_, request_,
          boost::bind(&Client::handle_write_request, this,
            boost::asio::placeholders::error));
    }
    else if (endpoint_iterator != tcp::resolver::iterator())
    {
      // The connection failed. Try the next endpoint in the list.
      socket_.close();
      tcp::endpoint endpoint = *endpoint_iterator;
      socket_.async_connect(endpoint,
          boost::bind(&Client::handle_connect, this,
            boost::asio::placeholders::error, ++endpoint_iterator));
    }
    else
    {
      std::cout << "Error: " << err.message() << "\n";
    }
}

void SimOutputStreaming::Client::handle_write_request(const boost::system::error_code& err)
{
    if (!err)
    {
      // Read the response status line.
      boost::asio::async_read_until(socket_, response_, "\r\n",
          boost::bind(&Client::handle_read_status_line, this,
            boost::asio::placeholders::error));
    }
    else
    {
      std::cout << "Error: " << err.message() << "\n";
    }
}

void SimOutputStreaming::Client::handle_read_status_line(const boost::system::error_code& err)
{
    if (!err)
    {
      // Check that response is OK.
      std::istream response_stream(&response_);
      std::string http_version;
      response_stream >> http_version;
      unsigned int status_code;
      response_stream >> status_code;
      std::string status_message;
      std::getline(response_stream, status_message);
      if (!response_stream || http_version.substr(0, 5) != "HTTP/")
      {
        std::cout << "Invalid response\n";
        return;
      }
      if (status_code != 200)
      {
        std::cout << "Response returned with status code ";
        std::cout << status_code << "\n";
        return;
      }

      // Read the response headers, which are terminated by a blank line.
      boost::asio::async_read_until(socket_, response_, "\r\n\r\n",
          boost::bind(&Client::handle_read_headers, this,
            boost::asio::placeholders::error));
    }
    else
    {
      std::cout << "Error: " << err << "\n";
    }
}

void SimOutputStreaming::Client::handle_read_headers(const boost::system::error_code& err)
{
    if (!err)
    {
      // Process the response headers.
      std::istream response_stream(&response_);
      std::string header;
      while (std::getline(response_stream, header) && header != "\r");
        // std::cout << header << "\n";
      //std::cout << "\n";

      // Write whatever content we already have to output.
      if (response_.size() > 0)
        std::cout << &response_;

      // Start reading remaining data until EOF.
      boost::asio::async_read(socket_, response_,
          boost::asio::transfer_at_least(1),
          boost::bind(&Client::handle_read_content, this,
            boost::asio::placeholders::error));
    }
    else
    {
      std::cout << "Error: " << err << "\n";
    }
}

void SimOutputStreaming::Client::handle_read_content(const boost::system::error_code& err)
{
    if (!err)
    {
      // Write all of the data that has been read so far.
      std::cout << &response_;

      // Continue reading remaining data until EOF.
      boost::asio::async_read(socket_, response_,
          boost::asio::transfer_at_least(1),
          boost::bind(&Client::handle_read_content, this,
            boost::asio::placeholders::error));
    }
    else if (err != boost::asio::error::eof)
    {
      std::cout << "Error: " << err << "\n";
    }
}

