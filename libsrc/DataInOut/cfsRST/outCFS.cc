#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>

#include "outCFS.hh"

namespace CoupledField
{

 WriteResultsCFS ::WriteResultsCFS(const Char * const filename)
    :WriteResults(filename)
  {
    ENTER_FCN( "WriteResultsUnverg::WriteResultsUnverg" );

    //create the result directory
    std::string nameDir = namefile_ + "_Results";
    std::string S="mkdir -p ";

    S += nameDir;
    system(S.c_str());

    //set pointer of output file to NULL
    output = NULL;

    currentSequence_ = 1;
    currentStep_     = 1;

    //open the first result file
    OpenNewResultFile();

  }

 WriteResultsCFS :: ~WriteResultsCFS()
  {
    ENTER_FCN( "WriteResultsCFS::~WriteResultsCFS" );
    if (output) {
      (*output) << "</multiSequenceStep>" << std::endl << std::endl;
      (*output) << "</cfsStepFile>" << std::endl << std::endl;

      output->close();
    }
  }


  void WriteResultsCFS :: OpenNewResultFile()
  {
    ENTER_FCN( "WriteResultsCFS::OpenNewResultFile" );

    //close current output file
   if (output) {
     (*output) << "</multiSequenceStep>" << std::endl << std::endl;
     (*output) << "</cfsStepFile>" << std::endl << std::endl;
     output->close();
     delete output;
   }
 
   //create new output file
   std::ostringstream mystream;
   mystream << namefile_ << "_Results/" << namefile_ << "_S" << currentSequence_ 
	    << "_" << currentStep_ << ".rst";

   std::string name = mystream.str();
   std::cout << "name = " << name << std::endl;
   output = new std::ofstream(name.c_str());
   WriteHeader();

  }


  void WriteResultsCFS :: WriteHeader()
  {
    ENTER_FCN( "WriteResultsCFS::WriteHeader" );


    (*output) << "<?xml version=\"1.0\"?>" << std::endl << std::endl;
    (*output) << "<cfsStepFile version=\"1.0\"" << std::endl << std::endl;

    (*output) << "<muliSequenceStep step=\"" << currentSequence_ << "\" step=\"" 
              <<  currentStep_ << "\">" << std::endl << std::endl;


  }


  void WriteResultsCFS::NodeElemDataTransient(const UInt dataSetNr,
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

  void WriteResultsCFS::NodeElemDataHarmonic(const UInt dataSetNr,
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

  void WriteResultsCFS::Init(Grid * aptgrid)
  {
    ptGrid_=aptgrid;

    // Initialize history files
    InitHistoryFiles();
  }

  void WriteResultsCFS::WriteNodeSolutionTransient(const NodeStoreSol<Double> & sol, 
                                                       const UInt step, 
                                                       const Double time)
  {
  
    ENTER_FCN( "WriteResultsCFS::WriteNodeSolutionTransient" );
    
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

  void WriteResultsCFS::WriteElemSolutionTransient(const ElemStoreSol<Double>& sol, 
                                                       const UInt step, 
                                                       const Double time)
  {
    ENTER_FCN( "WriteResultsCFS::WriteElemSolutionTransient" );

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

  void WriteResultsCFS::WriteNodeSolutionHarmonic(const NodeStoreSol<Complex> & sol, 
                                                      const UInt step,
                                                      const Double frequency, 
                                                      const ComplexFormat format)
  {
  
    ENTER_FCN( "WriteResultsCFS::WriteNodeSolutionHarmonic" );
  
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


  void WriteResultsCFS::WriteElemSolutionHarmonic(const ElemStoreSol<Complex>& sol, 
                                                      const UInt step,
                                                      const Double frequency, 
                                                      const ComplexFormat format)
  {
    ENTER_FCN( "WriteResultsCFS::WriteElemSolutionHarmonic" );
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

  std::string WriteResultsCFS::SolutionTypeToString(const SolutionType type) const
  {
    ENTER_FCN( "WriteResultsCFS::SolutionTypeToString" );

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
      default:
        Error( "Wrong type of solution or 'SolutionType2String' not implemented for \
this type of solution", __FILE__, __LINE__);
      }
    return std::string();
  }


} // end of namespace
