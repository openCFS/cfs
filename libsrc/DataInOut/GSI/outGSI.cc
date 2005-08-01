#include <iostream>
#include <fstream>
#include <string>

#ifdef __sgi
#include <stdio.h>
#else
#include <cstdio>
#endif

#include "outGSI.hh"

#include "GSIBaseIO.hh"
#include "GSIRawIO.hh"
#include "GSIXDRIO.hh"
#include "GSIIOException.hh"

namespace CoupledField {

  // **************
  //  Constructor
  // **************
  WriteResultsGSI::WriteResultsGSI( const Char * const filename )
    : WriteResults(filename) {

    ENTER_FCN( "WriteResultsGSI::WriteResultsGSI" );

    fp_ = NULL;
    io_ = NULL;

    std::string myfilename(filename);
    myfilename += ".gsi";
  
    fp_ = fopen(myfilename.c_str(), "wb");
  
    io_ = new GridlibSocketInterface::RawIO(NULL, fp_, true);

  }


  // *************
  //  Destructor
  // *************
  WriteResultsGSI::~WriteResultsGSI() {

    ENTER_FCN( "WriteResultsGSI::~ WriteResultsGSI" );

    delete io_;
    fclose(fp_);
  }


  // ************
  //  WriteGrid
  // ************
  void WriteResultsGSI::WriteGrid() {

    ENTER_FCN( "WriteResultsGSI::WriteGrid" );


    if ( !io_ ) {
      Error( "File for output results is not initialized", __FILE__,
             __LINE__ );
    }

    try {
  
      // This is the number of timesteps/frequency steps 
      // at the moment this is one for test purposes
      (*io_) << (int) 1;
      (*io_) << "---- BEGIN GRID ----";
      //          (*io_) << (int) 1;
    }
    catch (GridlibSocketInterface::IOException& e) {
      std::cerr << "Exception in WriteResultsGSI::WriteGrid" << std::endl
                << e.GetErrorString() << std::endl;
    }
      
    Dataset666();
    Dataset781();
    Dataset780(); 

  }


  // ***********
  //  WriteMsg
  // ***********
  void WriteResultsGSI::WriteMsg(const std::string msg) {

    ENTER_FCN( "WriteResultsGSI::WriteMsg" );

    try {
      (*io_) << msg;
    }
    catch ( GridlibSocketInterface::IOException& e ) {
      std::cerr << "Exception in WriteResultsGSI::WriteMsg" << std::endl
                << e.GetErrorString() << std::endl;
    }
  }


  // ********
  //   Init
  // ********
  void WriteResultsGSI::Init( Grid * aptgrid ) {
    ENTER_FCN( "WriteResultsGSI::Init" );
    ptGrid_ = aptgrid;

    // Initialize history files
    InitHistoryFiles();
  }


  UInt WriteResultsGSI::DoCompute() {
    ENTER_FCN("WriteResultsGSI::DoCompute");

    std::string in;

    try {
      (*io_) >> in;
    }
    catch (GridlibSocketInterface::IOException& e) {
      std::cerr << "Exception in WriteResultsGSI::DoCompute" << std::endl
                << e.GetErrorString() << std::endl;
      return 0;
    }

    if(in == "---- COMPUTE ----")
      return 1;
    else
      return 0;
  }

  UInt WriteResultsGSI::TransferGrid() {

    ENTER_FCN("WriteResultsGSI::TransferGrid");


    std::string in;

    try {
      (*io_) >> in;
    }
    catch (GridlibSocketInterface::IOException& e) {
      std::cerr << "Exception in WriteResultsGSI::TransferGrid()"
                << std::endl
                << e.GetErrorString() << std::endl;
      return 0;
    }

    if(in == "---- TRANSFER GRID ----")
      return 1;
    else
      return 0;
  }

  void  WriteResultsGSI::Dataset666()
  {
    ENTER_FCN("WriteResultsGSI::Dataset666");

    if (!ptGrid_)
      Error("ptGrid_ is not initialized", __FILE__,__LINE__);

    UInt dim=ptGrid_->GetDim();
    UInt maxnumnodes=ptGrid_-> GetNumNodes();
    UInt maxnumelem=ptGrid_-> GetNumVolElems();

    try 
      {
        (*io_) << "---- DS666 ----";
        (*io_) << (int) dim << (int) maxnumnodes << (int) maxnumelem;

        //        std::cerr << "---- DS666 ---- writing dim " << dim << " maxnumnodes "
        //                  <<  maxnumnodes << " maxnumelem " << maxnumelem << std::endl;
      }
    catch (GridlibSocketInterface::IOException& e)
      {
        std::cerr << "Exception in WriteResultsGSI::Dataset666" << std::endl
                  << e.GetErrorString() << std::endl;
      }
  }

  void  WriteResultsGSI::Dataset781()
  {
    ENTER_FCN("WriteResultsGSI::Dataset781");

    if (!ptGrid_)
      Error("ptGrid_ is not initialized", __FILE__,__LINE__);

    UInt i;
    UInt dim=ptGrid_->GetDim();
    UInt maxnumnodes=ptGrid_->GetNumNodes();

    try 
      {
        (*io_) << "---- DS781 ----";
      
        float *vec = new float[maxnumnodes*3];
      
        if (dim==2) 
          {
            for (i=0; i<maxnumnodes; i++)
              {
                Point<2> Point;
                ptGrid_->GetNodeCoordinate(Point,i+1);
              
                vec[i*3+0] = (float) Point[0];
                vec[i*3+1] = (float) Point[1];
                vec[i*3+2] = (float) 0.0;
              }
          }
      
        if (dim==3)
          {
            for (i=0; i<maxnumnodes; i++)
              {
                Point<3> Point;
                ptGrid_->GetNodeCoordinate(Point,i+1);
              
                vec[i*3+0] = (float) Point[0];
                vec[i*3+1] = (float) Point[1];
                vec[i*3+2] = (float) Point[2];
              }
          }
      
        io_->WriteFloatArray(vec, maxnumnodes*3);
        //        std::cerr << "---- DS781 ---- writing maxnumnodes*3 floats "
        //                  << (maxnumnodes*3) << std::endl;
      }
    catch (GridlibSocketInterface::IOException& e)
      {
        std::cerr << "Exception in WriteResultsGSI::Dataset781" << std::endl
                  << e.GetErrorString() << std::endl;
      }
  }

  void  WriteResultsGSI::Dataset780()
  {
    ENTER_FCN("WriteResultsGSI::Dataset780");

    if (!ptGrid_)
      Error("ptGrid_ is not initialized", __FILE__,__LINE__);
  
    UInt maxnumelem;
    UInt dim=ptGrid_->GetDim();
  
    UInt elmsgrp;
  
    StdVector<RegionIdType> subdoms;
    StdVector<Elem*> elemssd;
    StdVector<UInt> connect;

    std::vector<int> connectsizes;
    std::vector<int> subdomains;
    std::vector<int> indexe;

    UInt i, j, k, l;

    elmsgrp=1;
    ptGrid_->GetVolRegionIds(subdoms);
    k = 0;
    l = 0;
    maxnumelem=ptGrid_->GetNumVolElems();

    try 
      {
        (*io_) << "---- DS780 ----";

        connectsizes.resize(maxnumelem);
        subdomains.resize(maxnumelem);
  
        for (i=0; i<subdoms.GetSize(); i++)
          {
            ptGrid_->GetVolElems(elemssd,subdoms[i]);
          
            //          (*io_) << (int) elemssd.size();
            for (j=0; j < elemssd.GetSize(); j++)
              {  
                connect=elemssd[j]->connect;
                connectsizes[k] = connect.GetSize();
                subdomains[k] = i;
                k++;

                //              k++; 
                //              connect=elemssd[j]->connect;
              
                //              (*io_) << (int) connect.size();
              
                UInt ii=0;
              
                /*              
                //Wegen diesem Zeug hier nachfragen
                if (elmsgrp == 3) {
                  
                  
                UInt jkl;
                for (jkl=1; jkl<5; jkl++) {
                (*io_) << (int) maxnumnodes_+jkl-2;
                }
            
                ii=4;
                }
                */

                for (; ii < connect.GetSize(); ii++) 
                  { 
                    indexe.push_back((int)connect[ii]-1);
                    l++;
                    //                  (*io_) << (int) connect[ii]-1;
                  }
              
                // write the element's subdomain
                //              (*io_) << (int) i;

              } // over elements of group
            elmsgrp++;
          } // over groups

        (*io_) << connectsizes;
        (*io_) << subdomains;
        (*io_) << indexe;


        //        std::cerr << "---- DS780 ---- subdoms " << (int) subdoms->GetSize()
        //                  << std::endl;
      }
    catch (GridlibSocketInterface::IOException& e) 
      {
        std::cerr << "Exception in WriteResultsGSI::Dataset780" << std::endl
                  << e.GetErrorString() << std::endl;
      }
  }


  //void  WriteResultsGSI::Dataset55(const std::string & title,
  // const Vector<Double> & x, const UInt step, const Double time)void  
  // nodal results
  void WriteResultsGSI::Dataset55_Transient(const std::string & title, 
                                            const Vector<Double> & x, 
                                            const UInt step, 
                                            const Double time, 
                                            const UInt nrNodes,
                                            const UInt nrDofs)
  {
    ENTER_FCN("WriteResultsGSI::Dataset55_Transient");

    if (!ptGrid_)
      Error("ptGrid_ is not initialized", __FILE__,__LINE__);

    UInt i,j,n;
    std::vector<float> vec;

    n = x.GetSize();
    vec.resize(n);

    try 
      {
        (*io_) << "---- DS55 TRANSIENT ----";
        (*io_) << title;      
        (*io_) << (int) step;
        (*io_) << (float) time;
        (*io_) << (int) nrNodes;
        (*io_) << (int) nrDofs;
 
        for (i=0; i<n; i++)
          {
            vec[i] = (float) x[i];
          }
 
        io_->WriteFloatVector(vec);
      }
    catch (GridlibSocketInterface::IOException& e) 
      {
        std::cerr << "Exception in WriteResultsGSI::Dataset55_Transient"
                  << std::endl
                  << e.GetErrorString() << std::endl;
      }
  }  

  void WriteResultsGSI::Dataset55_Harmonic(const std::string & title, 
                                           const Vector<Complex> & x, 
                                           const UInt step,
                                           const Double frequency, 
                                           const ComplexFormat format,
                                           const UInt nrNodes,
                                           const UInt nrDofs)
  {
    if (!ptGrid_)
      Error("ptGrid_ is not initialized", __FILE__,__LINE__);
  
    UInt i,j,n;
    n=nrNodes;
    std::string id;
    std::vector<float> vec1, vec2;

    vec1.resize(n*nrDofs);
    vec2.resize(n*nrDofs);

    if (format == REAL_IMAG)
      {
        id = "---- DS55 HARMONIC RI ----";

        for (i=0; i<n; i++)
          {
            for (j=0; j<nrDofs; j++)
              {
                vec1[i*nrDofs+j] = (float) x[i*nrDofs +j].real();
                vec2[i*nrDofs+j] = (float) x[i*nrDofs +j].imag();
              }
          }
      }
    else if (format == AMPLITUDE_PHASE) 
      {
        id = "---- DS55 HARMONIC AP ----";
        
        for (i=0; i<n; i++)
          {
            for (j=0; j<nrDofs; j++)
              {
                vec1[i*nrDofs+j] = (float) std::abs(x[i*nrDofs +j]);

                if (abs(x[i*nrDofs +j].imag()) > 1e-16)
                  vec2[i*nrDofs+j] = (float) std::arg(x[i*nrDofs +j])*180/PI;
                else
                  vec2[i*nrDofs+j] = (float) 0.0;
              }
          }
      }

    try 
      {
        (*io_) << id;
        (*io_) << title;      
        (*io_) << (int) step;
        (*io_) << (float) frequency;
        (*io_) << (int) nrNodes;
        (*io_) << (int) nrDofs;
 
        io_->WriteFloatVector(vec1);
        io_->WriteFloatVector(vec2);
      }
    catch (GridlibSocketInterface::IOException& e) 
      {
        std::cerr << "Exception in WriteResultsGSI::Dataset55_Harmonic"
                  << std::endl << e.GetErrorString() << std::endl;
      }  
  }


  // void  WriteResultsGSI::Dataset56_Scalar(const std::string & title,
  // const Vector<Double> & x, const UInt step, const Double time)
  void  WriteResultsGSI::Dataset56_Transient(const std::string & title, 
                                             const Vector<Double> & x, 
                                             const UInt step, 
                                             const Double time, 
                                             const UInt numElems,
                                             const UInt nrDofs)
  {
    ENTER_FCN("WriteResultsGSI::Dataset56_Transient");

    UInt i,n;
    std::vector<float> vec;

    n = x.GetSize();
    vec.resize(n);

    if (!ptGrid_)
      Error("ptGrid_ is not initialized", __FILE__,__LINE__);

    try {
      (*io_) << "---- DS56 TRANSIENT ----";
      (*io_) << title;      
      (*io_) << (int) step;
      (*io_) << (float) time;
      (*io_) << (int) numElems;
      (*io_) << (int) nrDofs;
 
      for (i=0; i<n; i++)
        {
          vec[i] = (float) x[i];
        }
 
      io_->WriteFloatVector(vec);
    }
    catch (GridlibSocketInterface::IOException& e) 
      {
        std::cerr << "Exception in WriteResultsGSI::Dataset56_Transient"
                  << std::endl << e.GetErrorString() << std::endl;
      }
  }  


  void WriteResultsGSI::Dataset56_Harmonic(const std::string & title, 
                                           const Vector<Complex> & x, 
                                           const UInt step,
                                           const Double frequency, 
                                           const ComplexFormat format, 
                                           const UInt numElems,
                                           const UInt nrDofs)
  {
    if (!ptGrid_)
      Error("ptGrid_ is not initialized", __FILE__,__LINE__);
  
    UInt i,j,n;
    n=numElems;
    std::string id;
    std::vector<float> vec1, vec2;

    vec1.resize(n*nrDofs);
    vec2.resize(n*nrDofs);

    if (format == REAL_IMAG)
      {
        id = "---- DS56 HARMONIC RI ----";

        for (i=0; i<n; i++)
          {
            for (j=0; j<nrDofs; j++)
              {
                vec1[i*nrDofs+j] = (float) x[i*nrDofs +j].real();
                vec2[i*nrDofs+j] = (float) x[i*nrDofs +j].imag();
              }
          }
      }
    else if (format == AMPLITUDE_PHASE) 
      {
        id = "---- DS56 HARMONIC AP ----";
       
        for (i=0; i<n; i++)
          {
            for (j=0; j<nrDofs; j++)
              {
                vec1[i*nrDofs+j] = (float) std::abs(x[i*nrDofs +j]);

                if (abs(x[i*nrDofs +j].imag()) > 1e-16)
                  vec2[i*nrDofs+j] = (float) std::arg(x[i*nrDofs +j])*180/PI;
                else
                  vec2[i*nrDofs+j] = (float) 0.0;
              }
          }
      }

    try 
      {
        (*io_) << id;
        (*io_) << title;      
        (*io_) << (int) step;
        (*io_) << (float) frequency;
        (*io_) << (int) numElems;
        (*io_) << (int) nrDofs;
 
        io_->WriteFloatVector(vec1);
        io_->WriteFloatVector(vec2);
      }
    catch (GridlibSocketInterface::IOException& e) 
      {
        std::cerr << "Exception in WriteResultsGSI::Dataset55_Harmonic"
                  << std::endl << e.GetErrorString() << std::endl;
      }
  }

  void WriteResultsGSI ::
  WriteNodeSolutionTransient( const NodeStoreSol<Double> & sol, 
                              const UInt step, 
                              const Double time )
  { 
    ENTER_FCN( "WriteResultsGSI::WriteNodeSolutionTransient" );
 
  
    Vector<Double> globalSolution;
    StdVector<SolutionType> solTypes;
    UInt numNodes =  ptGrid_->GetNumNodes();
    std::string title;

    sol.GetSolutionTypes(solTypes);
  
    for (UInt iSol=0; iSol<solTypes.GetSize(); iSol++)
      {
        sol.GetGlobalSolVector(solTypes[iSol],globalSolution);
        title = SolutionTypeToString(solTypes[iSol]);
        Dataset55_Transient(title, globalSolution, step, 
                            time, numNodes ,sol.GetDof(solTypes[iSol]));
      }
  
    if(solTypes.GetSize() > 0)
    {
      try 
        {
          (*io_) << "---- END OF RESULTS ----";
        }
      catch (GridlibSocketInterface::IOException& e) 
        {
          std::cerr << "Exception in WriteResultsGSI::WriteNodeSolutionTransient"
                    << std::endl << e.GetErrorString() << std::endl;
        }
    }
    
  }

  void
  WriteResultsGSI::WriteElemSolutionTransient( const ElemStoreSol<Double>& sol,
                                               const UInt step,
                                               const Double time )
  {
    ENTER_FCN( "WriteResultsGSI::WriteElemSolutionTransient" );

    Vector<Double> globalSolution;
    StdVector<SolutionType> solTypes;
    std::string title;

    UInt numElems =  ptGrid_->GetNumVolElems();  
  
    sol.GetSolutionTypes(solTypes);
    sol.TransformElemSolution(globalSolution,ptGrid_);
    title = SolutionTypeToString(solTypes[0]);
    Dataset56_Transient(title, globalSolution, step, 
                        time, numElems, sol.GetDof());

    if(solTypes.GetSize() > 0)
    {
      try 
        {
          (*io_) << "---- END OF RESULTS ----";
        }
      catch (GridlibSocketInterface::IOException& e) 
        {
          std::cerr << "Exception in WriteResultsGSI::WriteElemSolutionTransient"
                    << std::endl << e.GetErrorString() << std::endl;
        }
    }
  }

  void
  WriteResultsGSI::WriteNodeSolutionHarmonic( const NodeStoreSol<Complex>& sol,
                                              const UInt step,
                                              const Double frequency, 
                                              const ComplexFormat format )
  {

    ENTER_FCN( "WriteResultsGSI::WriteNodeSolutionHarmonic" );
  
    Vector<Complex> globalSolution;
    StdVector<SolutionType> solTypes;
    UInt numNodes =  ptGrid_->GetNumNodes();
    std::string title;  
    sol.GetSolutionTypes(solTypes);
  
    for (UInt iSol=0; iSol<solTypes.GetSize(); iSol++)
      {
        sol.GetGlobalSolVector(solTypes[iSol],globalSolution);
        title = SolutionTypeToString(solTypes[iSol]);
        Dataset55_Harmonic(title, globalSolution, step, frequency, 
                           format, numNodes ,sol.GetDof(solTypes[iSol]));
      }
  
    if(solTypes.GetSize() > 0)
    {
      try 
        {
          (*io_) << "---- END OF RESULTS ----";
        }
      catch (GridlibSocketInterface::IOException& e) 
        {
          std::cerr << "Exception in WriteResultsGSI::WriteNodeSolutionHarmonic"
                    << std::endl << e.GetErrorString() << std::endl;
        }
    }

  }

  void
  WriteResultsGSI::WriteElemSolutionHarmonic( const ElemStoreSol<Complex>& sol,
                                              const UInt step,
                                              const Double frequency, 
                                              const ComplexFormat format )
  {

    ENTER_FCN( "WriteResultsGSI::WriteElemSolutionHarmonic" );
  
    std::string title;  
    Vector<Complex> globalSolution;
    StdVector<SolutionType> solTypes;
    sol.GetSolutionTypes(solTypes);
    UInt numElems =  ptGrid_->GetNumVolElems(); 
  
    for (UInt iSol=0; iSol<solTypes.GetSize(); iSol++)
      {
        sol.GetGlobalSolVector(solTypes[iSol],globalSolution);
        title = SolutionTypeToString(solTypes[iSol]);
        Dataset56_Harmonic(title, globalSolution, step, frequency, 
                           format, numElems ,sol.GetDof(solTypes[iSol]));
      }
  
    if(solTypes.GetSize() > 0)
    {
      try 
        {
          (*io_) << "---- END OF RESULTS ----";
        }
      catch (GridlibSocketInterface::IOException& e) 
        {
          std::cerr << "Exception in WriteResultsGSI::WriteElemSolutionHarmonic"
                    << std::endl << e.GetErrorString() << std::endl;
        }
    }
  }

  std::string WriteResultsGSI::SolutionTypeToString(const SolutionType type)
    const
  {
    ENTER_FCN( "WriteResultsGSI::SolutionTypeToString" );

    switch (type)
      {
      case MECH_DISPLACEMENT:
        return "displacement";
        break;
      case MECH_ACCELERATION:
        return "acceleration";
        break;
      case MECH_VELOCITY:
        return "velocity";
        break;
      case MECH_FORCE:
        return "mechanic force";
        break;
      case MECH_STRESS:
        return "stress";
        break;
      case MECH_STRAIN:
        return "strain";
        break; 
      case MECH_ENERGY:
        return "energy";
        break;
      case ELEC_POTENTIAL:
        return "electric potential";
        break;
      case ELEC_FIELD_INTENSITY:
        return "electric field";
        break;
      case ELEC_FORCE_VWP: 
        return "electric force";
        break;
      case ELEC_CHARGE:
        return "electric charge";
        break;
      case ELEC_FLUX_DENSITY:
        return "electric flux density";
        break;
      case ELEC_ENERGY:
        return "electric energy";
        break;
      case SMOOTH_DISPLACEMENT:
        return "displacement";
        break;
      case ACOU_POTENTIAL:
        return "fluid potential";
        break;
      case ACOU_FORCE:
        return "acou force";
        break;
      case ACOU_POTENTIAL_DERIV_1:
        return "fluid potential, 1st deriv.";
        break;
      case ACOU_POTENTIAL_DERIV_2:
        return "fluid potential, 2nd deriv.";
        break;
      case MAG_POTENTIAL:
        return "mag. vector potential";
        break;
      case MAG_FLUX_DENSITY:
        return "mag. flux density";
        break;
      case MAG_EDDY_CURRENT:
        return "eddy current";
        break;
      case MAG_FORCE_VWP:
        return "mag. forceVWP";
        break;
      case MAG_FORCE_LORENTZ:
        return "mag. force Lorentz";
        break; 
      case MAG_ENERGY:
        return "mag. energy";
        break;
      default:
        (*error) << "Wrong type of solution or 'SolutionType2String' not "
                 << "implemented for this type of solution";
        Error( __FILE__, __LINE__ );
      }
  }

} // end of namespace
