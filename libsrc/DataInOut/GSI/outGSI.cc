#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

#include "outGSI.hh"
#include "DataInOut/ParamHandling/ConfFile.hh"
#include "DataInOut/socketopts.hh"

#include "CFSSocketComm/CFS_ServerSocket.h"
#include "CFSSocketComm/CFS_SocketException.h"
#include "CFSSocketComm/CFS_BaseIO.h"
#include "CFSSocketComm/CFS_RawIO.h"
#include "CFSSocketComm/CFS_XDRIO.h"
#include "CFSSocketComm/CFS_IOException.h"

namespace CoupledField
{

WriteResultsGSI :: WriteResultsGSI(const Char * const filename, Boolean withHistory)
: WriteResults(filename,withHistory)
{
#ifdef TRACE
  (*trace) << "entering WriteResultsGSI :: WriteResultsGSI" << std::endl;
#endif

 ptgrid=NULL;
 currstep_=0;

}

int WriteResultsGSI::waitForConnection() {
#ifdef TRACE
  (*trace) << "entering WriteResultsGSI :: waitForConnection" << std::endl;
#endif

  std::string okmsg = "++++ CFS++ SERVER READY! ++++";
  std::string failedmsg = "++++ SOCKET CREATION FAILED ++++";
  std::string connfailmsg = "++++ ERROR ESTABLISHING CONNECTION ++++";

  std::cout << "++++ sockopts port: " << sockopts->getPort() << " ++++" << std::endl;
  CFS_ServerSocket ss2_(sockopts->getPort());
  if(!ss2_.is_valid()) 
  {
      std::cout << failedmsg << std::endl;
      return 1;
  }
  

  sock_ = new CFS_ServerSocket();
  std::cout << okmsg << std::endl;

  try
  {
      ss2_.accept ( *sock_ );
  }
  catch ( CFS_SocketException& e )
  {
      std::cerr << "Exception was caught:" << e.description() << std::endl;
      std::cout << connfailmsg << std::endl;
      return 2;
  }

  if(!sock_->is_valid()) {
      std::cout << connfailmsg << std::endl;      
      return 2;
  }
  
  io_ = new CFS_XDRIO(sock_->getReadHandle(), sock_->getWriteHandle());
  return 0;
}

WriteResultsGSI ::~WriteResultsGSI()
{
#ifdef TRACE
  (*trace) << "entering WriteResultsGSI::~ WriteResultsGSI" << std::endl;
#endif

 // write keyword
  //  (*io_) << "endgmv  ";

 // if (ascii_) (*output) << std::endl;

 // delete [] namedir_;

 // delete output;
 delete io_;
 delete sock_;
}

void WriteResultsGSI::WriteGrid(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering WriteResultsGSI :: WriteGrid" << std::endl;
#endif
  //  std::cout << "Writing Grid...\n";
  Dataset666(level);
  Dataset781(level);
  Dataset780(level);
}

void WriteResultsGSI::WriteSolution(const Vector<Double> & sol, const Integer step, const Double time, const std::string title)
{
#ifdef TRACE
 (*trace) << " entering WriteResultsGSI::WriteSolution " << std::endl;
#endif

 Dataset55(title, sol, step+1, time);
}

void WriteResultsGSI::WriteDataOnCell(const Vector<Double>&sol,const Integer step, const Double time, const std::string title)
{
#ifdef TRACE
  (*trace) << "entering WriteResultsGSI :: WriteDataOnCell" << std::endl;
#endif

  std::string aux;
  if (title == "elecfield") {
      aux="electric field";
      Dataset56_Scalar(aux, sol, step+1, time);
  }
  else 
    Dataset56_Scalar(title, sol, step+1, time);
}

void WriteResultsGSI::WriteVecDataOnCell(const Vector<Double>*vec,const Integer step, const Double time, const std::string title)
{
#ifdef TRACE
  (*trace) << "entering WriteResultsGSI :: WriteVecDataOnCell" << std::endl;
#endif

  std::string aux;
  if (title == "elecfield") {
      aux="electric field";
      Dataset56_Vec(aux, vec, step+1, time);
  }
  else
    // Soll ich das so lassen? Hier können doch auch andere Vektoren ausgegeben werden.
    Dataset56_Vec("electric field", vec, step+1, time);
}

void WriteResultsGSI::WriteMsg(const std::string msg) {
#ifdef TRACE
  (*trace) << "entering WriteResultsGSI :: WriteMsg" << std::endl;
#endif

  try {
    (*io_) << msg;
  }
  catch (CFS_IOException& e) {
      std::cerr << "Exception in WriteResultsGSI::WriteMsg" << std::endl \
                << e.description() << std::endl;
  }  
}

void WriteResultsGSI::Init(Grid * aptgrid)
{
#ifdef TRACE
  (*trace) << "entering WriteResultsGSI :: Init" << std::endl;
#endif

  ptgrid=aptgrid;
}


Integer WriteResultsGSI::DoCompute() {
#ifdef TRACE
  (*trace) << "entering WriteResultsGSI :: DoCompute" << std::endl;
#endif

  std::string in;
  
  try {
      (*io_) >> in;
  }
  catch (CFS_IOException& e) {
      std::cerr << "Exception in WriteResultsGSI::DoCompute" << std::endl \
                << e.description() << std::endl;
      return 0;
  }

  if(in == "---- COMPUTE ----")
    return 1;
  else
    return 0;
}

Integer WriteResultsGSI::TransferGrid() {
#ifdef TRACE
  (*trace) << "entering WriteResultsGSI :: TransferGrid" << std::endl;
#endif

  std::string in;

  try {
      (*io_) >> in;
  }
  catch (CFS_IOException& e) {
      std::cerr << "Exception in WriteResultsGSI::TransferGrid()" << std::endl \
                << e.description() << std::endl;
      return 0;
  }

  if(in == "---- TRANSFER GRID ----")
    return 1;
  else
    return 0;
}

void  WriteResultsGSI::Dataset666(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering WriteResultsGSI :: Dataset666" << std::endl;
#endif

  Integer dim;
  Integer maxnumnodes;
  
  if (!ptgrid)
    Error("ptgrid is not initialized", __FILE__,__LINE__);

  try {
      (*io_) << "---- DS666 ----";

      dim=ptgrid->GetDim();
      // maxnumnodes and visible nodes
      maxnumnodes=ptgrid-> GetMaxnumnodes(level);// + 4;
      maxnumnodes_=maxnumnodes;
      // Nachfragen!!!
      // if (flags->VisElements_) maxnumnodes+=4;
      Integer maxnumelem=ptgrid-> GetMaxnumElem(level);
      
      (*io_) << (int) dim << (int) maxnumnodes << (int) maxnumelem;
      
  }
  catch (CFS_IOException& e) {
      std::cerr << "Exception in WriteResultsGSI::Dataset666" << std::endl \
                << e.description() << std::endl;
  }
}

void  WriteResultsGSI::Dataset781(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering WriteResultsGSI :: Dataset781" << std::endl;
#endif

  Integer dim;
  Integer maxnumnodes;
  
  if (!ptgrid)
    Error("ptgrid is not initialized", __FILE__,__LINE__);

  try {
      (*io_) << "---- DS781 ----";
      
      dim=ptgrid->GetDim();
      maxnumnodes=ptgrid->GetMaxnumnodes(level);
      
      float *vec = new float[maxnumnodes*3];

      if (dim==2) 
      {
          for (Integer i=0; i<maxnumnodes; i++)
          {
              // (*io_) << (int) i+1 << (int) 0 << (int) 0 << (int) 11;
             
              Point<2> Point;
              ptgrid->GetCoordinateNode(i,level,Point);
              
              vec[i*3+2] = (float) 0.0;
              vec[i*3+1] = (float) Point[1];
              vec[i*3+0] = (float) Point[0];
          }
      }
      
      Integer igl=0;
      if (dim==3)
      {
          for (igl=0; igl<maxnumnodes; igl++)
          {
              //	 (*io_) << (int) igl+1 << (int) 0 << (int) 0 << (int) 11;
              
              Point<3> Point;
              ptgrid->GetCoordinateNode(igl,level,Point);
              
              vec[igl*3+2] = (float) Point[2];
              vec[igl*3+1] = (float) Point[1];
              vec[igl*3+0] = (float) Point[0];
          }
      }
      
      /*      
      Integer ilm=0;
      for (ilm=0; ilm<4; ilm++, igl++)
      {
          
          vec[igl*3+2] = (float) (*(visNodes_[ilm]))[2];
          vec[igl*3+1] = (float) (*(visNodes_[ilm]))[1];
          vec[igl*3+0] = (float) (*(visNodes_[ilm]))[0];        
      }
      */
      
      
      io_->writeFloatArray(vec, maxnumnodes*3);
  } 
  catch (CFS_IOException& e) {
      std::cerr << "Exception in WriteResultsGSI::Dataset781" << std::endl \
                << e.description() << std::endl;
  }


}

void  WriteResultsGSI::Dataset780(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering WriteResultsGSI :: Dataset780" << std::endl;
#endif

  Integer dim;
  Integer maxnumelem;
  Vector<Integer> connect;
  std::vector<Elem*> elemssd;
  Integer elmsgrp;
  std::vector<std::string>* subdoms;
  Integer i, j, k;
  
  if (!ptgrid)
    Error("ptgrid is not initialized", __FILE__,__LINE__);
  
  try {
      (*io_) << "---- DS780 ----";
      dim=ptgrid->GetDim();

      maxnumelem=ptgrid->GetMaxnumElem(level);

      elmsgrp=1;

      subdoms=ptgrid->GetAllSDs();
      k = 0;

      (*io_) << (int) subdoms->size();
  
      for (i=0; i<subdoms->size(); i++)
      {
          ptgrid->GetElemSD(elemssd,(*subdoms)[i],level);
          
          (*io_) << (int) elemssd.size();
          for (j=0; j < elemssd.size(); j++)
          {  
              k++; 
              connect=elemssd[j]->connect;
              
              (*io_) << (int) connect.size();
              
              Integer ii=0;
              
              /*              
              //Wegen diesem Zeug hier nachfragen
              if (elmsgrp == 3) {
                  
                  
                  Integer jkl;
                  for (jkl=1; jkl<5; jkl++) {
                      (*io_) << (int) maxnumnodes_+jkl-2;
                  }
	    
                  ii=4;
              }
              */

              for (; ii < connect.size(); ii++) 
              { 
                  (*io_) << (int) connect[ii]-1;
              }

          } // over elements of group
          elmsgrp++;
      } // over groups
  }
  catch (CFS_IOException& e) {
      std::cerr << "Exception in WriteResultsGSI::Dataset780" << std::endl \
                << e.description() << std::endl;
  }
}


void  WriteResultsGSI::Dataset55(const std::string & title, const Vector<Double> & x, const Integer step, const Double time)
{
#ifdef TRACE
  (*trace) << "entering WriteResultsGSI :: Dataset55" << std::endl;
#endif

  Integer i,n;
  std::vector<float> vec;

  if (!ptgrid)
     Error("ptgrid is not initialized", __FILE__,__LINE__);

  try {
      (*io_) << "---- DS55 ----";

      (*io_) << title;
      
      (*io_) << (int) step << (float) time;
 
      n=x.size();
      vec.resize(n);
      for (i=0; i<n; i++)
        vec[i] = (float) x[i];
 
      io_->writeFloatVector(vec);

      //  std::cout << "+++ DS55 +++ length: " << n << " " << title <<
      //  std::endl;
  }
  catch (CFS_IOException& e) {
      std::cerr << "Exception in WriteResultsGSI::Dataset55" << std::endl \
                << e.description() << std::endl;
  } 

}  

void  WriteResultsGSI::Dataset56_Scalar(const std::string & title, const Vector<Double> & x, const Integer step, const Double time)
{
#ifdef TRACE
  (*trace) << "entering WriteResultsGSI :: Dataset56_Scalar" << std::endl;
#endif

  Integer i,n;
  std::vector<float> vec;

  if (!ptgrid)
    Error("ptgrid is not initialized", __FILE__,__LINE__);

  try {
      (*io_) << "---- DS56S ----";
      
      (*io_) << title << (int) step << (float) time;
      
      n=x.size();  
      vec.resize(n);
      for (i=0; i<n; i++)
        vec[i] = (float) x[i];
      
      io_->writeFloatVector(vec);
      
      // std::cout << "+++ DS56S +++ length: " << n << std::endl;
  }
  catch (CFS_IOException& e) {
      std::cerr << "Exception in WriteResultsGSI::Dataset56_Scalar" << std::endl \
                << e.description() << std::endl;
  } 
}  

void  WriteResultsGSI::Dataset56_Vec(const std::string & title, const Vector<Double> * x, const Integer step, const Double time)
{
#ifdef TRACE
  (*trace) << "entering WriteResultsGSI :: Dataset56_Vec" << std::endl;
#endif

  Integer i,n;
  std::vector<float> vec;

  if (!ptgrid)
    Error("ptgrid is not initialized", __FILE__,__LINE__);

  try {
      (*io_) << "---- DS56V ----";

      (*io_) << title << (int) step << (float) time;  
      
      n=x[0].size();
      vec.resize(n*3);
      for (i=0; i<n; i++) 
      {    
          vec[i*3+0] = (float) x[0][i];
          vec[i*3+1] = (float) x[1][i];
          vec[i*3+2] = (float) x[2][i];
      }

      io_->writeFloatVector(vec);
  
      //  std::cout << "+++ DS56V +++ length: " << n*3 << std::endl;
  }
  catch (CFS_IOException& e) {
      std::cerr << "Exception in WriteResultsGSI::Dataset56_Vec" << std::endl \
                << e.description() << std::endl;
  } 

}


} // end of namespace
