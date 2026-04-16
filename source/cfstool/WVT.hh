// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef WVT_HEADER_2013
#define WVT_HEADER_2013

#include <string>

#include <fstream>

#include "General/Environment.hh"
#include "DataInOut/SimInput.hh"

namespace CoupledField 
{
  class BaseBOperator;
}

namespace CFSTool
{

  /** This is the WVT class for weight vector theory-based post-processing.
   *
   *  This is the WVT class for the evaluation of post-processing results
   *  specific to coriolis mass flow meters (CMFs) according to the
   *  weight vector theory developed in \cite Hemp1994 and \cite Hemp1995. */
  class WVT 
  {
  public:

    WVT( const PtrParamNode& param,
         const PtrParamNode& info );

    ~WVT();
    
    typedef enum {PRIMARY_MODE, SECONDARY_MODE, MEAN_FLOW, EVAL_GRID, OUTPUT} FileType;
    static Enum<FileType> fileType;

    typedef enum {MF_GRID_DATA, MF_ANALYTIC_EXP, MF_SCATTERED_DATA} MeanFlowDataType;
    static Enum<MeanFlowDataType> meanFlowDataType;

    void PostProcess();
    
  private:

    void OpenIntPointsFile();
    void CloseIntPointsFile();

    void AcouPot2AcouVelFeFct(const shared_ptr<BaseFeFunction>& acouPotFeFct,
                              const shared_ptr<BaseBOperator>& gradOp,
                              shared_ptr<BaseFeFunction>& acouVelFeFct);

    static void WriteResultsToCSV(Double freq,
                                  Complex u_p_prime, 
                                  Double deltaPhiVol,
                                  Double deltaPhiSurf,
                                  Double vol,
                                  Double meanVel,
                                  Double meanVelCorrectionFactor);

    static void WriteResultsToMatlab(Double freq,
                                     Complex u_p_prime, 
                                     Double deltaPhiVol,
                                     Double deltaPhiSurf,
                                     Double vol,
                                     Double meanVel,
                                     Double meanVelCorrectionFactor);

    void Initialize();
    
  private:
    const PtrParamNode& param_;
    const PtrParamNode& info_;

    //! Dimension of grid
    UInt dim_;
    
    //! Number of DOFs for fluid velocity vectors
    UInt numDofs_;

    //! Names of DOFs for fluid velocity vectors
    StdVector<std::string> dofNames_;

    //! Global math parser instance.
    shared_ptr<MathParser> mp_;

    //! File   name  of   input  file   for  the   evaluation  on   which  the
    //! post-processing integrals  will be  evaluated.  If not  specified, the
    //! primary mode file will be used.
    std::string evalGridFile_;

    //! File name of input file for primary mode. Usually, this is the lateral
    //! mode.
    std::string priModeFile_;

    //! File name of input file for secondary mode. Usually, this is the coriolis
    //! mode.
    std::string secModeFile_;

    //! File name of input file for mean flow.
    std::string meanFlowFile_;

    //! File name of output file, if output for post-processing weight vectors
    //! should be written.
    std::string outFile_;
        
    bool writeOutputFile_;

    //! Name of named sensor node inside the primary mode HDF5 file. The harmonic
    //! mechanic displacement value is read for this named node either from the
    //! history branch or, if history for this node is not existant, the value
    //! is read from the mesh result branch.
    std::string sensorNodeName_;

    //! Integration order
    UInt integOrder_;

    //! Do primary and secondary mode input files come from a direct coupled FSI simulation?
    bool dirCoupled_;

    Complex u_p_;
    
    typedef std::map< FileType, shared_ptr<SimInput> > IOType;
    IOType inputs_;

    bool nodalResults_;

    //! What kind of mean flow data do we use? Data given on grid points either from
    //! HDF5 file or given as analytical expression. Or do we use data on a point cloud
    //! read from a CSV file and provided through the CoefFunctionScatteredData.
    MeanFlowDataType meanFlowType_;

    //! CoefFunction which holds mean flow data as point cloud and interpolates
    //! the data to the integration points using neareast neighbor mapping.
    shared_ptr<CoefFunction> meanFlowCoefScattered_;

    //! File in legacy VTK format for writing out a point cloud of
    //! volume integration points. These points can then be read into a CFD code
    //! for evaluation of the mean velocity. The file can also be read into
    //! ParaView for visualization.
    std::ofstream intPointsFile_;

    //! Shall we write the integration point file?
    bool writeIntPointsFile_;

    std::string intPointsFileName_;

    long intPointsFilePos_;
    long numIntPoints_;
  };

}

#endif
