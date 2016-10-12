#ifndef SIMOUTSTREAMING_HH_
#define SIMOUTSTREAMING_HH_

#include <iosfwd>
#include <string>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/SimOutput.hh"
#include "Utils/StdVector.hh"
#include <boost/asio.hpp>
#include <boost/bind.hpp>

namespace CoupledField {
class BaseResult;
class Grid;
}  // namespace CoupledField

using boost::asio::ip::tcp;

namespace CoupledField
{

  /** This class collects the information, the basic text writer writes to files
   * in to info.xml output.
   * This writer is intialized by default */
  class SimOutputStreaming : public SimOutput
  {
  public:
    SimOutputStreaming(PtrParamNode outputNode, PtrParamNode infoNode, bool isRestart);

    virtual ~SimOutputStreaming();

    //! Initialize class
    void Init(Grid * ptGrid, bool printGridOnly );

    //! Register result (within one multisequence step)
    void RegisterResult(shared_ptr<BaseResult> sol, UInt saveBegin, UInt saveInc,
                                 UInt saveEnd, bool isHistory );

    //! Begin multisequence step
    void BeginMultiSequenceStep(UInt step, BasePDE::AnalysisType type, UInt numSteps);


    //! Begin single analysis step
    void BeginStep( UInt stepNum, Double stepVal );

    /** This transmits the data */
    void FinishStep();

    /** Mesh results are node and cell results */
    void AddResult( shared_ptr<BaseResult> sol );

    /** Used for calculating size of Http Content Length */
    UInt GetContentLength();

    /** Called by the Http Client.
     * @param out here the data has to be written, the wrapping is done by Client */
    void Transmit(std::ostream& out);

  private:

    /** accumulated results by AddMeshResult() to be written and released by FinishStep() */
    StdVector<shared_ptr<BaseResult> > results_;

    /** copied from boost
     * http://www.boost.org/doc/libs/1_37_0/doc/html/boost_asio/example/http/client/async_client.cpp */
    class Client
    {
    public:
      Client(boost::asio::io_service& io_service,
          const std::string& server, const std::string& port, const std::string& path, SimOutputStreaming* base);

    private:
      void handle_resolve(const boost::system::error_code& err,
          tcp::resolver::iterator endpoint_iterator);

      void handle_connect(const boost::system::error_code& err,
          tcp::resolver::iterator endpoint_iterator);

      void handle_write_request(const boost::system::error_code& err);

      void handle_read_status_line(const boost::system::error_code& err);

      void handle_read_headers(const boost::system::error_code& err);

      void handle_read_content(const boost::system::error_code& err);

      tcp::resolver resolver_;
      tcp::socket socket_;
      boost::asio::streambuf request_;
      boost::asio::streambuf response_;
    };

    /** do http streaming or file? */
    bool http_;

    /** by default localhost */
    std::string host_;

    /** this is a number or http which is equivalent to 80 */
    std::string port_;
    
    /** this is the path on the server the POST is sent to */
    std::string path_;

    /** add the mesh ? */
    bool send_mesh_;

    /** compression is not zipped but no identation and no 'nr' in result/data/item */
    bool compressed_;

    /** should we output more information to the command line? */
    bool silent_;

    /** this is the output data. overwritten to save multiple grid writing */
    PtrParamNode content_;
  };
}



#endif /* SIMOUTSTREAMING_HH_ */
