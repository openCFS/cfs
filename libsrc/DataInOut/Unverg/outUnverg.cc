#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include "outUnverg.hh"

namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  WriteResultsUnverg::WriteResultsUnverg( const Char* const filename )
    : WriteResults( filename ) {

    ENTER_FCN( "WriteResultsUnverg::WriteResultsUnverg" );

    std::string name = namefile_ + ".unv";
    output = NULL;
    output = new std::ofstream(name.c_str());
  }


  // **************
  //   Destructor
  // **************
  WriteResultsUnverg::~WriteResultsUnverg() {
    ENTER_FCN( "WriteResultsUnverg::~WriteResultsUnverg" );
    delete output;
  }


  // *************
  //   WriteGrid
  // *************
  void WriteResultsUnverg::WriteGrid() {

    ENTER_FCN( "WriteResultsUnverg::WriteGrid" );

    if ( !output ) {
      Error( "File for output results is not initialized", __FILE__,
             __LINE__ );
    }

    Dataset666();
    Dataset781();
    Dataset780();
  }


  void  WriteResultsUnverg::Dataset666() {

    if (!ptGrid_) {
      Error("ptGrid_ is not initialized",  __FILE__, __LINE__ );
    }

    (*output) << std::setw(6) << -1 << std::endl << std::setw(6)
              << -666 << std::endl ;

    UInt dim=ptGrid_->GetDim();
    UInt maxnumnodes=ptGrid_->GetNumNodes();
    UInt maxnumelem=ptGrid_->GetNumVolElems();

    (*output) << std::setw(10) << 1 << std::setw(10) << 1 << std::setw(10)
              << dim << std::endl << std::setw(10) << maxnumnodes
              << std::setw(10) << maxnumelem << std::endl;

    (*output) << std::setw(6) << -1 << std::endl;
  }

  void  WriteResultsUnverg::Dataset781() {
    //
    if (!ptGrid_)
      Error("ptGrid_ is not initialized", __FILE__,__LINE__);

    (*output) << std::setw(6) << -1 << std::endl << std::setw(6) << 781
              << std::endl;

    UInt dim=ptGrid_->GetDim();
    UInt maxnumnodes=ptGrid_->GetNumNodes();

    (*output).setf(std::ios::scientific);
    (*output).precision(16);

    for ( UInt i = 0; i < maxnumnodes; i++ ) {
      (*output) << std::setw(10) << i+1 << std::setw(10) << 0
                << std::setw(10) << 0 << std::setw(10) << 11 << std::endl;

      (*output).setf(std::ios::uppercase);

      if ( dim==2 ) {
        Point<2> Point;
        ptGrid_->GetNodeCoordinate(Point,i+1);

        (*output) << "   " << 0.0 ;
        PrintPoint(Point,output);
        (*output) << std::endl;
      }
      else {
        Point<3> Point;
        ptGrid_->GetNodeCoordinate(Point,i+1);

        PrintPoint(Point,output);
        (*output) << std::endl;
      }
    }
 
    (*output) << std::setw(6) << -1 << std::endl;
  }

  void  WriteResultsUnverg::Dataset780()
  {
    //
    if (!ptGrid_)
      Error("ptGrid_ is not initialized", __FILE__,__LINE__);

    (*output) << std::setw(6) << -1 << std::endl << std::setw(6)
              << 780 << std::endl;

    UInt dim = ptGrid_->GetDim();

    StdVector<UInt> connect;
    StdVector<Elem*> elemssd;
    UInt elmsgrp=1;
    std::string errMsg;

    StdVector<RegionIdType> subdoms;
    ptGrid_->GetVolRegionIds(subdoms);
    UInt i, j, k;
    k = 0;
    for (i=0; i<subdoms.GetSize(); i++) {

      ptGrid_->GetVolElems(elemssd,subdoms[i]);

      for ( j = 0; j < elemssd.GetSize(); j++ ) {  
        k++;
        connect=elemssd[j]->connect;
        (*output) << std::setw(10) << elemssd[j]->elemNum << std::setw(10);

        if ( dim == 2 ) {
          switch(connect.GetSize()) {
          case 3: (*output) << 91; break;
          case 4: (*output) << 94; break;
          case 6: (*output) << 92; break;
          case 8: (*output) << 95; break;
          default:
            (*error) << "Please, put element type according to "
                     << "unverg-format for this number of nodes "
                     << "per element";
            Error( __FILE__, __LINE__ );
          }

          (*output) << std::setw(10) << 2 << std::setw(10) << 2
                    << std::setw(10) << 1 << std::setw(10) << 1
                    << std::setw(10) << elmsgrp << std::setw(10)
                    << connect.GetSize() << std::endl;
        }
        else {
          switch(connect.GetSize()) {
          case  4: (*output) << 111; break;  // tetraeder 1.ord
          case  6: (*output) << 112; break;  // prism     1.ord
          case  8: (*output) << 115; break;  // hexaeder  1.ord
          case 15: (*output) << 113; break;  // prism     2.ord
          case 20: (*output) << 116; break;  // hexaeder  2.ord
          default:
            (*error) << "Please, put element type according to "
                     << "unverg-format for " << connect.GetSize()
                     << " nodes per Element in 3D!";
            Error( __FILE__,__LINE__);
          }

          (*output) << std::setw(10) << 11 << std::setw(10) << 1
                    << std::setw(10) << 1  << std::setw(10) << 1
                    << std::setw(10) << elmsgrp << std::setw(10) 
                    << connect.GetSize() << std::endl;
        }

        if (dim == 2 && (connect.GetSize() == 6 || connect.GetSize() == 8)) {
          //quadratic elements
          UInt offset = UInt(connect.GetSize()/2);
          for (UInt ii=0; ii < offset; ii++) {
            (*output).width(10);
            (*output) << connect[ii];
            (*output).width(10);
            (*output) << connect[ii+offset]; 
          }
        }
        else {
          for (UInt ii=0; ii < connect.GetSize(); ii++) { 
            (*output).width(10);
            (*output) << connect[ii];
          }
        }

        (*output) << std::endl;
      } // over elements of group
      elmsgrp++;
    } // over groups
    (*output) << std::setw(6) << -1 << std::endl;
  }

  void  WriteResultsUnverg::NodeElemDataTransient(const UInt dataSetNr,
                                                  const std::string & title, 
                                                  const Vector<Double> & x, 
                                                  const UInt step, 
                                                  const Double time, 
                                                  const UInt nrNodes,
                                                  const UInt nrDofs)
  {
    //
    if (!ptGrid_)
      Error("ptGrid_ is not initialized", __FILE__,__LINE__);

    (*output) << std::setw(6) << -1 << std::endl 
              << std::setw(6) << dataSetNr << std::endl;

    (*output).setf(std::ios::scientific);
    (*output).precision(6);
    (*output).setf(std::ios::uppercase);

    // Standard for scalar values
    UInt dataCharac = 1;
    UInt valsPerNode = 1;
 
    // needed for undocumented value of
    // Dataset 55/56 in record8 , field4
    UInt specDataCharac = 0;

    // Vector type
    if (nrDofs > 1 && nrDofs <= 3) 
      {
        valsPerNode = 3;
        dataCharac = 2;
      } 
    // Tensor
    else if(nrDofs == 6) 
      { 
        valsPerNode = 6;
        specDataCharac = 2;
        dataCharac = 4; // symmetric tensor
        //      dataCharac = 3; // vector with 6 components
        //      dataCharac = 5; // unsymmetric tensor
      }
     
 
    (*output) << " " << title << " step" << std::setw(6) << step <<
      " time   " << time << std::endl;  
    (*output) << std::endl << std::endl << std::endl << std::endl;
    (*output) << std::setw(10) << 1 << std::setw(10) << 4 << std::setw(10) 
              << dataCharac  << std::setw(10) << specDataCharac
              << std::setw(10) << 2 << std::setw(10) << valsPerNode << std::endl;
    (*output) << std::setw(10) << 2 << std::setw(10) << 1 << std::setw(10) << 1 
              << std::setw(10) <<
      step << std::endl;
    (*output) << " " << time << std::endl;       

    UInt i,j,n;
    n=nrNodes;;  
    for (i=0; i<n; i++)
      {
     
        (*output) << std::setw(10) << i+1;
        if (dataSetNr == 56)
          (*output) << std::setw(10) << valsPerNode;

        (*output) << std::endl;
     
        // in the universal file either one or three results datas must exist
        if (nrDofs == 2)
          (*output) << 0.0;

        for (j=0; j<nrDofs; j++)
          {
            //std::cerr << "trying to write " << i << ", " << j << std::endl;
            (*output) << std::setw(14) << x[i*nrDofs +j];
          }
     
        (*output) << std::endl;
      }    
    (*output) << std::setw(6) << -1 << std::endl;
  }  

  void WriteResultsUnverg::NodeElemDataHarmonic(const UInt dataSetNr,
                                                const std::string & title, 
                                                const Vector<Complex> & x, 
                                                const UInt step,
                                                const Double frequency, 
                                                const ComplexFormat format,
                                                const UInt nrNodes,
                                                const UInt nrDofs)
  {
  
    UInt dataCharact = 1;
    if (!ptGrid_)
      Error("ptGrid_ is not initialized", __FILE__,__LINE__);
  
    (*output) << std::setw(6) << -1 << std::endl 
              << std::setw(6) << dataSetNr << std::endl;
  
    (*output).setf(std::ios::scientific);
    (*output).precision(6);
    (*output).setf(std::ios::uppercase);
  
    UInt valsPerNode = 1;
    if (nrDofs > 1)
      {
        dataCharact = 2;
        valsPerNode = 3;
      }
 

    if (format == REAL_IMAG)
      {
        // write out realpart
        (*output) << " " << title << " cw realpart" << std::setw(6) <<" frequency   " << frequency << std::endl;  
        (*output) << std::endl << std::endl << std::endl << std::endl;
        (*output) << std::setw(10) << 1 << std::setw(10) << 5 << std::setw(10) << dataCharact << std::setw(10) << 0
                  << std::setw(10) << 2 << std::setw(10) << valsPerNode << std::endl;
        (*output) << std::setw(10) << 2 << std::setw(10) << 1 << std::setw(10) << -2 << std::setw(10) <<
          step << std::endl;
        (*output) << " " << frequency << std::endl;       
      
        UInt i,j,n;
        n=nrNodes;
        for (i=0; i<n; i++)
          {
            (*output) << std::setw(10) << i+1 << std::endl;
          
            // in the universal file either one or three results datas must exist
            if (nrDofs == 2)
              (*output) << 0.0;
          
            for (j=0; j<nrDofs; j++)
              {
                //std::cerr << "trying to write " << i << ", " << j << std::endl;
                (*output) << std::setw(14) << x[i*nrDofs +j].real();
              }
          
            (*output) << std::endl;
          }    
        (*output) << std::setw(6) << -1 << std::endl;

        // write out imag part
        (*output) << std::setw(6) << -1 << std::endl << std::setw(6) << 55 << std::endl;
        (*output) << " " << title << " cw imagpart" << std::setw(6) <<" frequency   " << frequency << std::endl;  
        (*output) << std::endl << std::endl << std::endl << std::endl;
        (*output) << std::setw(10) << 1 << std::setw(10) << 5 << std::setw(10) << dataCharact << std::setw(10) << 0
                  << std::setw(10) << 2 << std::setw(10) << valsPerNode << std::endl;
        (*output) << std::setw(10) << 2 << std::setw(10) << 1 << std::setw(10) << -12 << std::setw(10) <<
          step << std::endl;
        (*output) << " " << frequency << std::endl;       
      
        for (i=0; i<n; i++)
          {
            (*output) << std::setw(10) << i+1 << std::endl;
          
            // in the universal file either one or three results datas must exist
            if (nrDofs == 2)
              (*output) << 0.0;
          
            for (j=0; j<nrDofs; j++)
              {
                //std::cerr << "trying to write " << i << ", " << j << std::endl;
                (*output) << std::setw(14) << x[i*nrDofs +j].imag();
              }
          
            (*output) << std::endl;
          }    
        (*output) << std::setw(6) << -1 << std::endl;
      }
  
    else if (format == AMPLITUDE_PHASE) {
      // write out amplitude
      (*output) << " " << title << " cw amplitude" << std::setw(6) <<" frequency   " << frequency << std::endl;  
      (*output) << std::endl << std::endl << std::endl << std::endl;
      (*output) << std::setw(10) << 1 << std::setw(10) << 5 << std::setw(10) << dataCharact << std::setw(10) << 0
                << std::setw(10) << 2 << std::setw(10) << valsPerNode << std::endl;
      (*output) << std::setw(10) << 2 << std::setw(10) << 1 << std::setw(10) << -1 << std::setw(10) <<
        step << std::endl;
      (*output) << " " << frequency << std::endl;       
    
      UInt i,j,n;
      n=nrNodes;
      for (i=0; i<n; i++)
        {
          (*output) << std::setw(10) << i+1 << std::endl;
        
          // in the universal file either one or three results datas must exist
          if (nrDofs == 2)
            (*output) << 0.0;
        
          for (j=0; j<nrDofs; j++)
            {
              //std::cerr << "trying to write " << i << ", " << j << std::endl;
              (*output) << std::setw(14) << std::abs(x[i*nrDofs +j]);
            }
        
          (*output) << std::endl;
        }    
      (*output) << std::setw(6) << -1 << std::endl;
    
      // write out phase
      (*output) << std::setw(6) << -1 << std::endl << std::setw(6) << 55 << std::endl;
      (*output) << " " << title << " cw phase" << std::setw(6) <<" frequency   " << frequency << std::endl;  
      (*output) << std::endl << std::endl << std::endl << std::endl;
      (*output) << std::setw(10) << 1 << std::setw(10) << 5 << std::setw(10) << dataCharact << std::setw(10) << 0
                << std::setw(10) << 2 << std::setw(10) << valsPerNode << std::endl;
      (*output) << std::setw(10) << 2 << std::setw(10) << 1 << std::setw(10) << -11 << std::setw(10) <<
        step << std::endl;
      (*output) << " " << frequency << std::endl;       
    
      for (i=0; i<n; i++)
        {
          (*output) << std::setw(10) << i+1 << std::endl;
        
          // in the universal file either one or three results datas must exist
          if (nrDofs == 2)
            (*output) << 0.0;
        
          for (j=0; j<nrDofs; j++)
            {
              if (abs(x[i*nrDofs +j].imag()) > 1e-16)
                (*output) << std::setw(14) << std::arg(x[i*nrDofs +j])*180/PI;
              else 
                (*output) << std::setw(14) << 0.0;
            }
        
          (*output) << std::endl;
        }    
      (*output) << std::setw(6) << -1 << std::endl;
    }
  
  }

  void  WriteResultsUnverg::Init(Grid * aptgrid)
  {
    ptGrid_=aptgrid;

    // Initialize history files
    InitHistoryFiles();
  }

  void  WriteResultsUnverg::WriteNodeSolutionTransient(const NodeStoreSol<Double> & sol, 
                                                       const UInt step, 
                                                       const Double time)
  {
  
    ENTER_FCN( "WriteResultsUnverg::WriteNodeSolutionTransient" );
    
    Vector<Double> globalSolution;
    StdVector<SolutionType> solTypes;
    UInt numNodes =  ptGrid_->GetNumNodes();
    std::string title;
    sol.GetSolutionTypes(solTypes);

    for (UInt iSol=0; iSol<solTypes.GetSize(); iSol++)
      {
        sol.GetGlobalSolVector(solTypes[iSol],globalSolution);

        title = SolutionTypeToString(solTypes[iSol]);

        NodeElemDataTransient(55,title, globalSolution, step, 
                              time, numNodes ,sol.GetDof(solTypes[iSol]));
      }
  }

  void  WriteResultsUnverg::WriteElemSolutionTransient(const ElemStoreSol<Double>& sol, 
                                                       const UInt step, 
                                                       const Double time)
  {
    ENTER_FCN( "WriteResultsUnverg::WriteElemSolutionTransient" );

    Vector<Double> globalSolution;
    StdVector<SolutionType> solTypes;
    std::string title;
    UInt numElems =  ptGrid_->GetNumVolElems();  
  
    sol.GetSolutionTypes(solTypes);
    sol.TransformElemSolution(globalSolution,ptGrid_);
    title = SolutionTypeToString(solTypes[0]);
    NodeElemDataTransient(56,title, globalSolution, step, 
                          time, numElems, sol.GetDof());
  }

  void  WriteResultsUnverg::WriteNodeSolutionHarmonic(const NodeStoreSol<Complex> & sol, 
                                                      const UInt step,
                                                      const Double frequency, 
                                                      const ComplexFormat format)
  {
  
    ENTER_FCN( "WriteResultsUnverg::WriteNodeSolutionHarmonic" );
  
    Vector<Complex> globalSolution;
  
    StdVector<SolutionType> solTypes;
    sol.GetSolutionTypes(solTypes);
  
    UInt numNodes =  ptGrid_->GetNumNodes();
    std::string title;

    for (UInt iSol=0; iSol<solTypes.GetSize(); iSol++)
      {
        sol.GetGlobalSolVector(solTypes[iSol],globalSolution);
        title = SolutionTypeToString(solTypes[iSol]);
        NodeElemDataHarmonic(55, title, globalSolution, step, frequency, 
                             format, numNodes ,sol.GetDof(solTypes[iSol]));
      }
  
  
  }


  void  WriteResultsUnverg::WriteElemSolutionHarmonic(const ElemStoreSol<Complex>& sol, 
                                                      const UInt step,
                                                      const Double frequency, 
                                                      const ComplexFormat format)
  {
    ENTER_FCN( "WriteResultsUnverg::WriteElemSolutionHarmonic" );
    Vector<Complex> globalSolution;
    StdVector<SolutionType> solTypes;
    UInt numElems =  ptGrid_->GetNumVolElems();  

    std::string title;

    sol.GetSolutionTypes(solTypes);
    sol.TransformElemSolution(globalSolution,ptGrid_);
    title = SolutionTypeToString(solTypes[0]);  

    NodeElemDataHarmonic(55, title, globalSolution, step, frequency, 
                         format, numElems ,sol.GetDof(solTypes[0]));
  }

  std::string WriteResultsUnverg::SolutionTypeToString(const SolutionType type) const
  {
    ENTER_FCN( "WriteResultsUnverg::SolutionTypeToString" );

    //   std::string warnMsg;

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
        break;
      case MECH_STRESS:
        return "stress";
        break;
      case MECH_STRAIN:
        Error("Not implemented", __FILE__, __LINE__);
        break;
      case ELEC_POTENTIAL:
        return "electric potential";
        break;
      case ELEC_FIELD_INTENSITY:
        return "electric field";
        break;
      case ELEC_FORCE_VWP: 
        Error("Not implemented", __FILE__, __LINE__);
        break;
      case ELEC_INTERFACE_FORCE:
        Error("Not implemented", __FILE__, __LINE__);
        break; 
      case ELEC_CHARGE:
        return "electric charge";
        break;
      case ELEC_FLUX_DENSITY:
        Error("Not implemented", __FILE__, __LINE__);
        break; 
      case ELEC_ENERGY:
        Error("Not implemented", __FILE__, __LINE__);
      case SMOOTH_DISPLACEMENT:
        return "displacement";
        break;
      case ACOU_POTENTIAL:
        return "fluid potential";
        break;
      case ACOU_PRESSUREXYZ:
        return "displacement";
        break;

      case ACOU_RHSVAL:
        return "fluid potential";
        break;
      case ACOU_PRESSURE:
        //       warnMsg = "Due to the restrictions in the .unv file format, the ";
        //       warnMsg += "acoustic pressure is written as acoustic (fluid) potential!";
        //       Warning(warnMsg.c_str(), __FILE__, __LINE__);
        return "fluid potential";
        break;
      case ACOU_FORCE:
        Error("Not implemented", __FILE__, __LINE__);
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
        Error("Not implemented", __FILE__, __LINE__);
        break;
      case MAG_FORCE_LORENTZ:
        Error("Not implemented", __FILE__, __LINE__);
        break;
      case MAG_ENERGY:
        Error("Not implemented", __FILE__, __LINE__);
        break;
      case HEAT_TEMPERATURE:
        return "temperature";
        break;

      default:
        Error( "Wrong type of solution or 'SolutionType2String' not implemented for \
this type of solution", __FILE__, __LINE__);
      }
    return std::string();
  }




  // void  WriteResultsUnverg::WriteDataOnCell(const Vector<Double> & sol, const UInt step, const Double time, const std::string title, const Integer nrDofs)
  // {
  //   std::string aux;
  //   if (title == "elecfield") {
  //     aux="electric field";
  //     Dataset56(aux, sol, step+1, time);
  //   }
  //   else 
  //     Dataset56(title, sol, step+1, time, nrDofs);
  //     //Warning("This cell-data is not printed, since this type is not supported by Capapost",__FILE__,__LINE__);
  // }

} // end of namespace
