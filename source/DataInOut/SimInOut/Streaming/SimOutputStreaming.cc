#include <fstream>
#include <iostream>
#include <new>
#include <set>
#include <sstream>

#include "Domain/Domain.hh"
#include "Domain/Mesh/Grid.hh"
#include "Domain/Results/ResultInfo.hh"
#include "SimOutputStreaming.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "boost/bind/bind.hpp"

using boost::asio::ip::tcp;

using namespace CoupledField;

DEFINE_LOG(SOS, "streaming")


SimOutputStreaming::SimOutputStreaming(PtrParamNode outputNode, PtrParamNode infoNode, bool isRestart) :
  SimOutput("", outputNode, infoNode, isRestart),
  io_service_thread(SimOutputStreaming::io_service_runner_wrapper, this), // initialize the sending thread
  io_context_work(boost::asio::make_work_guard(io_context))
{
  formatName_ = "streaming";
  capabilities_.insert(MESH_RESULTS);

  http_ = outputNode->Get("protocol")->As<string>() == "http";
  host_ = outputNode->Get("host")->As<string>();
  port_ = outputNode->Get("port")->As<string>();
  path_ = outputNode->Get("path")->As<string>();
  send_mesh_ = outputNode->Get("sendMesh")->As<bool>();
  silent_ = outputNode->Has("silent") ? outputNode->Get("silent")->As<bool>() : false;
  content_ = PtrParamNode(new ParamNode(ParamNode::INSERT));
  content_->SetName("cfsStreaming");

  infoNode->Get("streaming/settings")->SetValue(outputNode, false);
  if(!http_ && path_ == "/")
    infoNode->Get("streaming")->SetWarning("with streaming file protocol, give an output filename in path");

  //wait_var = new std::condition_variable();
  current_client = NULL;

  // start the sending thread
  io_service_thread.detach();
}

SimOutputStreaming::~SimOutputStreaming()
{
  // send last .info.xml to receive status="finished" and memory data
  // set force to true
  TransmitData(true);

  io_context_work.reset(); // tell the dummy work to stop

  io_context.run(); // joining the C++ thread is buggy for some reason
  // this will do it as well
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
   // augment the info.xml by a streaming tag to inform the receiving unit about the results to be processed
   PtrParamNode info_root = domain->GetInfoRoot()->Get("streaming")->Get(ParamNode::PROCESS);
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
  results_.Clear();
  actStep_ = stepNum;
  actStepVal_ = stepVal;
}

void SimOutputStreaming::AddResult( shared_ptr<BaseResult> sol)
{
  results_.Push_back(sol);
}

void SimOutputStreaming::FinishStep() {
  TransmitData(false);
}

void SimOutputStreaming::TransmitData(bool force) {
  if (http_) {
    /** we do not want to send more than one thing at a time,
     * otherwise client suicide will mess everything up **/
    if (force) {
      int i=0;
      while(current_client != NULL) {
        // kindof busy sleep because this will only be executed at the very end
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        i++;
        if (i>200) {
          // if this takes longer than 2 seconds, orphan the running client:
          current_client = NULL; // because we still want to try sending data
          // (to avoid deadlock)
          break;
        }
      }
    }
    if ( current_client == NULL ) {
      // this client will destroy itself upon finishing
      // is the requests fails, it will lock up and skip sending data
      // (except for the final data which will still be tried)
      current_client = new Client(io_context, this);
      current_client->Send(host_, port_, path_);
    }
  } else {
    std::ofstream out(path_.c_str());
    Transmit(out);
    out.close();
  }
}

UInt SimOutputStreaming::GetContentLength() {
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
  Grid* grid = domain->GetGrid();

  // add the complete info.xml treee
  // use own name and make sure we are not seed as stupid boost::any
  // also don't write warnings againt to cerr
  content_->Get("cfsInfo")->SetValue(domain->GetInfoRoot(), false, false);

  // add mesh if it is new
  if(send_mesh_ && !content_->Has("grid"))
    domain->GetGrid()->ExportGrid(content_->Get("grid"));

  // now the actual results
  PtrParamNode results = content_->Get("results");
  // delete all potential children -- shared pointers!
  results->GetChildren().Resize(0);

  for(unsigned int r = 0; r < results_.GetSize(); r++)
  {
    shared_ptr<BaseResult> sol = results_[r];
    bool cplx = sol->GetEntryType() == BaseMatrix::COMPLEX;
    shared_ptr<ResultInfo> resInfo = sol->GetResultInfo();

    PtrParamNode rpn = results->Get("result", ParamNode::APPEND);
    std::string tmp;
    resInfo->Enum2String(resInfo->definedOn, tmp); // how ugly! :(
    rpn->Get("solution")->SetValue(tmp);

    rpn->Get("name")->SetValue(resInfo->resultName);
    rpn->Get("region")->SetValue(sol->GetEntityList()->GetName());
    RegionIdType regid = grid->GetRegion().Parse(sol->GetEntityList()->GetName());

    unsigned int dofs = resInfo->dofNames.GetSize();
    rpn->Get("dofs")->SetValue(dofs);
    for(unsigned int d = 0; d < dofs; d++)
      rpn->Get("dof_" + boost::lexical_cast<std::string>(d))->SetValue(resInfo->dofNames[d]);

    rpn->Get("unit")->SetValue(resInfo->unit);
    if(cplx)
      rpn->Get("desiredComplexFormat")->SetValue(resInfo->complexFormat == REAL_IMAG ? "realImag" : "amplPhase");

    // we do the fast bulk block stuff as it saves a lot of time!
    StdVector<std::string>& block = rpn->GetFastBulkBlock();

    SingleVector* vec = sol->GetSingleVector();
    assert(vec != NULL && vec->GetSize() >= 0);

    // the result vector is a set of dofs
    unsigned int items = vec->GetSize() / dofs;

    StdVector<unsigned int> nodes; // nodes per region
    StdVector<Elem*> elems;   // // elements per region
    // TODO handle only simple cases here!
    assert(resInfo->definedOn == ResultInfo::NODE || resInfo->definedOn == ResultInfo::ELEMENT);
    if(resInfo->definedOn == ResultInfo::NODE)
      grid->GetNodesByRegion(nodes, regid);
    else
      grid->GetElems(elems, regid);

    block.Resize(items);

    for(unsigned int e = 0; e < items; e++)
    {
      std::stringstream ss;

      ss << "<item";

      for(unsigned int d = 0; d < dofs; d++)
      {
        ss << " v_" << d << "=\"";
        if(cplx)
          ss << vec->GetComplexEntry(dofs * e + d);
        else
          ss << vec->GetDoubleEntry(dofs * e + d);
        ss << "\"";
      }

      ss << " id=\"";
      // mandatory for nodes
      if(resInfo->definedOn == ResultInfo::NODE)
        ss << nodes[e];
      else
        ss << elems[e]->elemNum;
      ss << "\"/>";
      block[e] = ss.str();
    }
  }
  // write the stuff
  content_->ToXML(out, -99, true); // adjust element type!
}

void SimOutputStreaming::io_service_runner(void) {
  LOG_DBG(SOS) << "++ starting io_context thread ..." << std::endl;
  io_context.run();
  LOG_DBG(SOS) << "++ io_context thread ended" << std::endl;
}

void SimOutputStreaming::io_service_runner_wrapper(SimOutputStreaming* this_) {
  this_->io_service_runner();
}

SimOutputStreaming::Client::Client(boost::asio::io_context& io_context, SimOutputStreaming* base)
: resolver_(io_context),
  socket_(io_context) {
  base_ = base;}

void SimOutputStreaming::Client::Send(const std::string& server, const std::string& port, const std::string& path)
{
    // it is more robust known the content length a priori, therefore we use a memory stream
    std::stringstream mem;
    // fill with all our data
    base_->Transmit(mem);


    std::ostream request_stream(&request_);
    request_stream << "POST " << path << " HTTP/1.1\r\n";
    request_stream << "Host: " << server << "\r\n";
    request_stream << "Accept: */*\r\n";
    request_stream << "Connection: close\r\n";
    request_stream << "Content-Length: " << mem.tellp() << "\r\n";
    request_stream << "Content-Type: text/plain\r\n\r\n";

    // copy data
    request_stream << mem.rdbuf() << "\r\n\r\n";

    LOG_DBG(SOS) << "try to connect " << server << " port " << port << " to transmit " << mem.tellp() << " bytes" << std::endl;

    // Start an asynchronous resolve to translate the server and service names
    // into a list of endpoints.
    resolver_.async_resolve(server, port,
        boost::bind(&Client::handle_resolve, this,
          boost::asio::placeholders::error,
          boost::asio::placeholders::results));
}

void SimOutputStreaming::Client::handle_resolve(const boost::system::error_code& err,
      tcp::resolver::results_type endpoints)
{
    if (!err)
    {
      // Attempt a connection to the endpoints. Each endpoint
      // will be tried until we successfully establish a connection.
      boost::asio::async_connect(socket_, endpoints,
          boost::bind(&Client::handle_connect, this,
            boost::asio::placeholders::error,
            boost::asio::placeholders::endpoint));
    }
    else
    {
      std::cout << "Error: " << err.message() << "\n";
    }
}

void SimOutputStreaming::Client::handle_connect(const boost::system::error_code& err,
      const tcp::endpoint& endpoint)
{
    if (!err)
    {
      // The connection was successful. Send the request.
      boost::asio::async_write(socket_, request_,
          boost::bind(&Client::handle_write_request, this,
            boost::asio::placeholders::error));
      LOG_DBG(SOS) << "Connected to " << endpoint << std::endl;
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

  string msg = "Streaming to " + base_->host_ + ":" + base_->port_ + " results in ";

  if(!err)
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
      std::cout << msg + "invalid response\n";
      return;
    }
    if (status_code != 200)
    {
      std::cout << msg + "status code " << status_code << "\n";
      return;
    }

    // Read the response headers, which are terminated by a blank line.
    boost::asio::async_read_until(socket_, response_, "\r\n\r\n",
        boost::bind(&Client::handle_read_headers, this,
            boost::asio::placeholders::error));
  }
  else
  {
    std::cout << msg + "error: " << err << "\n";
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
      if (response_.size() > 0) {
        LOG_DBG(SOS) << "S:C:hrh: " << &response_;
      }

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
      LOG_DBG(SOS) << "S:C:hrc: " << &response_;

      // Continue reading remaining data until EOF.
      boost::asio::async_read(socket_, response_,
          boost::asio::transfer_at_least(1),
          boost::bind(&Client::handle_read_content, this,
            boost::asio::placeholders::error));
    }
    else {
      if (err != boost::asio::error::eof) {
        std::cout << "Error: " << err << "\n";
      } else {
        SimOutputStreaming* base_copy = base_;

        // commit suicide here
        // DO NOT ACCESS ANY MEMBER VARIABLES AFTER DELETE!
        // delete base_->current_client;
        delete this; // we are using this since the client might have been overridden by forcefully sending

        if (this == base_copy->current_client) // we are just accessing the pointer so this comparison is fine
          base_copy->current_client = NULL; // tell streaming that the client comitted suicide
      }
    }
}

