#include <string>
#include <iomanip>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/exception.hpp>
namespace fs = boost::filesystem;

#include <fstream>
#include "SimOutputUnv.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  SimOutputUnv::SimOutputUnv(  const std::string& filename,
                               PtrParamNode outputNode ) 
    : SimOutput ( filename, outputNode ) {
    std::string sysPathSep;
    
    formatName_ = "unv";
    dirName_ = "results_" + formatName_;
    outputNode->GetValue("directory", dirName_, ParamNode::PASS );
    try 
    {
      sysPathSep = fs::path("/").native_directory_string();
    } catch (std::exception &ex)
    {
      EXCEPTION(ex.what());
    }

    fileName_ = dirName_ + sysPathSep + filename;
    stepNumOffset_ = 0;
    stepValOffset_ = 0.0;
    
    capabilities_.insert( MESH );
    capabilities_.insert( MESH_RESULTS );
   
  }


  // **************
  //   Destructor
  // **************
  SimOutputUnv::~SimOutputUnv() {
    delete output;
  }


  // *************
  //   WriteGrid
  // *************
  void SimOutputUnv::WriteGrid() {


    if ( !output ) {
      EXCEPTION( "File for output results is not initialized" );
    }

    Dataset666();
    Dataset781();
    Dataset780();
  }


  void  SimOutputUnv::Dataset666() {

    if (!ptGrid_) {
      EXCEPTION("ptGrid_ is not initialized" );
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

  void  SimOutputUnv::Dataset781() {
    //
    if (!ptGrid_)
      EXCEPTION("ptGrid_ is not initialized" );

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

      Point Point;
      ptGrid_->GetNodeCoordinate(Point,i+1);
      if( dim == 3 ) {
        (*output) << "   " << Point[0];
        (*output) << "   " << Point[1];
        (*output) << "   " << Point[2];
      } else {
        (*output) << "   " << 0.0;
        (*output) << "   " << Point[0];
        (*output) << "   " << Point[1];
      }
      (*output) << std::endl;
    }
 
    (*output) << std::setw(6) << -1 << std::endl;
  }

  void  SimOutputUnv::Dataset780()
  {
    //
    if (!ptGrid_)
      EXCEPTION("ptGrid_ is not initialized");

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
            EXCEPTION( "Please, put element type according to "
                       << "unverg-format for this number of nodes "
                       << "per element" );
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
            EXCEPTION( "Please, put element type according to "
                       << "unverg-format for " << connect.GetSize()
                       << " nodes per Element in 3D!" );
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

  void  SimOutputUnv::NodeElemDataTransient(const UInt dataSetNr,
                                                  const std::string & title, 
                                                  const Vector<Double> & x, 
                                                  const UInt step, 
                                                  const Double time, 
                                                  const UInt nrNodes,
                                                  const UInt nrDofs)
  {
    //
    if (!ptGrid_)
      EXCEPTION("ptGrid_ is not initialized" );

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

  void SimOutputUnv::NodeElemDataHarmonic(const UInt dataSetNr,
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
      EXCEPTION("ptGrid_ is not initialized" );
  
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
                (*output) << std::setw(14) << CPhase(x[i*nrDofs +j]);
              else 
                (*output) << std::setw(14) << 0.0;
            }
        
          (*output) << std::endl;
        }    
      (*output) << std::setw(6) << -1 << std::endl;
    }
  
  }

  void SimOutputUnv::Init(Grid * aptgrid, bool printGridOnly )
  {
    try 
    {
      fs::create_directory( dirName_ );
    } catch (std::exception &ex)
    {
      EXCEPTION(ex.what());
    }

    std::string name = fileName_ + ".unv";
    output = NULL;
    output = new std::ofstream(name.c_str());
    if(!output)
      EXCEPTION("Unv file ' " << name << "' could not be openend!" );

    ptGrid_=aptgrid;
    WriteGrid();

  }

  void SimOutputUnv::RegisterResult( shared_ptr<BaseResult> sol,
                                     UInt saveBegin, UInt saveInc,
                                     UInt saveEnd,
                                     bool isHistory ) {
    
  }

  void SimOutputUnv::BeginStep( UInt stepNum, Double stepVal ) {
    
    actStep_ = stepNum + stepNumOffset_;
    actStepVal_ = stepVal + stepValOffset_;
    resultMap_.clear();
  }

  void SimOutputUnv::AddResult( shared_ptr<BaseResult>  sol ) {
    resultMap_[sol->GetResultInfo()->resultName].Push_back( sol );
  }
  
  void SimOutputUnv::FinishStep( ) {

    // iterate over all result types
    ResultMapType::iterator it = resultMap_.begin();
    for( ; it != resultMap_.end(); it++ ) {
      
      // check if result is defined on nodes or elements
      ResultInfo & actInfo = *(it->second[0]->GetResultInfo());
      
      const StdVector<shared_ptr<BaseResult> > actResults =
        it->second;

      if(!ValidateNodesAndElements(actInfo)) continue;      
      
      ResultInfo::EntityUnknownType entityType = actInfo.definedOn;
      std::string title =  SolutionTypeToString( actInfo.resultType );
      
      // if result could not be mapped, omit it
      if( title == "") {
        WARN("Result '" << actInfo.resultName
             << "' coul not be mappted to unv result type and is "
             << "omitted!");
        continue;
      }
      
      StdVector<std::string> & dofNames = actInfo.dofNames;
      UInt numDofs = dofNames.GetSize();
      
      // Determine type of entity the result is defined on
      UInt mapTo = 0;
      UInt numEntities = 0;
      if( actInfo.definedOn == ResultInfo::NODE ) {
        mapTo = 55;
        numEntities = ptGrid_->GetNumNodes();;
      } else {
        mapTo = 56;
        numEntities = ptGrid_->GetNumVolElems();
      }
      
      if( actResults[0]->GetEntryType() == BaseMatrix::DOUBLE ) {
        Vector<Double> gSol;
        FillGlobalVec<Double>(gSol, actResults, entityType );
        
        // Special case: If result is mechStress, restort entries 
        // according to capa notation
        if( actInfo.resultType == MECH_STRESS ) { 
          SortStresses( gSol, dofNames);
          numDofs = 6;
        }
        
        NodeElemDataTransient( mapTo, title, gSol, actStep_, 
                               actStepVal_ , numEntities ,
                               numDofs );
      } else {
        Vector<Complex> gSol;
        FillGlobalVec<Complex>(gSol, actResults, entityType );

        // Special case: If result is mechStress, restort entries 
        // according to capa notation
        if( actInfo.resultType == MECH_STRESS ) { 
          SortStresses( gSol, dofNames);
          numDofs = 6;
        }
        
        NodeElemDataHarmonic( mapTo, title, gSol, actStep_, 
                              actStepVal_, actInfo.complexFormat,
                              numEntities, numDofs );
      }
    } // over all result types
  }
  
  
  void SimOutputUnv::FinishMultiSequenceStep() {

    // set offset for step value and number to last values
    stepNumOffset_ = actStep_;
    stepValOffset_ = actStepVal_;
  }
  
  
  template<class TYPE>
  void SimOutputUnv::SortStresses( Vector<TYPE> & vec,

                                   const StdVector<std::string>& dofNames ){


    //soring according to capa (unv) notation
    //our notation is Voigt: Txx Tyy Tzz Tyz Txz Txy


    // check size of dofNames:
    // 6: 3D
    // 4: axi
    // 3: plane 

    // ensure, that number of entries is multiples of number of dofs
    assert( vec.GetSize() % dofNames.GetSize() == 0 );
    UInt numDofs = dofNames.GetSize();
    UInt numEntities = (UInt) vec.GetSize() / dofNames.GetSize();

    // CAPA always expects a stress vector with 6 entries
    Vector<TYPE> sorted;
    sorted.Resize( numEntities * 6 );
    
    // a) 3D-case (6 entries for stress vector)
    if( numDofs == 6 ) {
      
      //Txx Txy Tyy Txz Tyz Tzz
      UInt j = 0;
      for( UInt i = 0; i < vec.GetSize(); i+=6 ) {
        sorted[j++] = vec[i];
        sorted[j++] = vec[i+5];
        sorted[j++] = vec[i+1];
        sorted[j++] = vec[i+4];
        sorted[j++] = vec[i+3];
        sorted[j++] = vec[i+2];
      }
    }
  
    // b) axisymmetric-case (4 entries for stress vector)  
    else if( numDofs == 4 ) {
      //Tphiphi 0 Trr 0 Trz Tzz
      
      
      UInt j = 0;
      for( UInt i = 0; i < vec.GetSize(); i+=4 ) {
        sorted[j++] = vec[i+3];
        sorted[j++] = 0.0;
        sorted[j++] = vec[i+0];
        sorted[j++] = 0.0;
        sorted[j++] = vec[i+2];
        sorted[j++] = vec[i+1];
      }
    }

    // c) planestress/strain-case (3 entries for stress vector)  
    else if( numDofs == 3 ) {
      //0 0 Tyy 0 Tyz Tzz
    
      UInt j = 0;
      for( UInt i = 0; i < vec.GetSize(); i+=3 ) {
        sorted[j++] = 0.0;
        sorted[j++] = 0.0;
        sorted[j++] = vec[i+0];
        sorted[j++] = 0.0;
        sorted[j++] = vec[i+2];
        sorted[j++] = vec[i+1];
      } 

    } else {
      EXCEPTION( "Can not resort a stress vector with " << numDofs << " entries" );
    }
   
    // store sorted vector back to original one
    vec = sorted;
  }
  
  std::string SimOutputUnv::SolutionTypeToString(const SolutionType type) const
  {

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
        EXCEPTION("Not implemented" );
        break;
      case MECH_PSEUDO_DENSITY:
        return "pseudo density";
        break;
      case MECH_RHS_LOAD:
        return "mechRhsLoad";
         break;
      case ELEC_POTENTIAL:
        return "electric potential";
        break;
      case ELEC_FIELD_INTENSITY:
        return "electric field";
        break;
      case ELEC_FORCE_VWP: 
        EXCEPTION("Not implemented" );
        break;
      case ELEC_INTERFACE_FORCE:
        EXCEPTION("Not implemented" );
        break; 
      case ELEC_CHARGE:
        return "electric charge";
        break;
      case ELEC_FLUX_DENSITY:
        EXCEPTION("Not implemented");
        break; 
      case ELEC_ENERGY:
        EXCEPTION("Not implemented");
        break;
      case ELEC_RHS_LOAD:
        return "elecRhsLoad";
        break;
      case SMOOTH_DISPLACEMENT:
        return "displacement";
        break;
      case ACOU_POTENTIAL:
        return "fluid potential";
        break;
      case ACOU_PRESSUREXYZ:
        return "displacement";
        break;

      case ACOU_RHS_LOAD:
        return "fluid potential";
        break;
      case ACOU_PRESSURE:
        //       warnMsg = "Due to the restrictions in the .unv file format, the ";
        //       warnMsg += "acoustic pressure is written as acoustic (fluid) potential!";
        //       WARN(warnMsg.c_str());
        return "fluid potential";
        break;
      case ACOU_FORCE:
        EXCEPTION("Not implemented" );
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
        EXCEPTION("Not implemented");
        break;
      case MAG_FORCE_LORENTZ:
        EXCEPTION("Not implemented");
        break;
      case MAG_ENERGY:
        EXCEPTION("Not implemented");
        break;
      case MAG_RHS_LOAD:
        return "magRhsLoad";
        break;

      case HEAT_TEMPERATURE:
        return "temperature";
        break;
      case HEAT_RHS_LOAD:
        return "temperatureRhsLoad";
        break;

      default:
        return "";
      }
    return std::string();
  }




  // void  SimOutputUnv::WriteDataOnCell(const Vector<Double> & sol, const UInt step, const Double time, const std::string title, const Integer nrDofs)
  // {
  //   std::string aux;
  //   if (title == "elecfield") {
  //     aux="electric field";
  //     Dataset56(aux, sol, step+1, time);
  //   }
  //   else 
  //     Dataset56(title, sol, step+1, time, nrDofs);
  //     //WARN("This cell-data is not printed, since this type is not supported by Capapost");
  // }

}
