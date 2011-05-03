// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <string>
#include <math.h>
#include <iomanip>
#include <complex>

#include "OLAS/algsys/basesystem.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Driver/stdSolveStep.hh"
#include "Driver/singleDriver.hh"

#include "acouCombustion.hh"


extern "C" {
  void read_in(const char *filename, double** data, int nNodes, int nVars, int* nodeNr);
}

namespace CoupledField
{

  AcouCombustionNoise::AcouCombustionNoise( Grid *aptgrid, PtrParamNode paramNode )
    :AcousticPDE( aptgrid, paramNode )
  {

    // pdename_ is also acoustic for this case
    pdename_          = "acoustic";
    pdematerialclass_ = FLUID;
    
    isHarmonic_=false;

    // here are the fixed parameters according to the input file
    varSpeedOfSoundIdx_   = 0;
    srcReynoldStressIdx_  = 1;
    srcMomentumIdx_       = 2;
    srcChemicalIdx_       = 3;
    srcDensityDeriv2Idx_  = 4;
    srcHeatReleaseIdx_    = 5;
    srcGasConstDeriv2Idx_ = 6;
    srcShearTermIdx_      = 7;

    //read din parameters from xml-file
    StdVector<Elem*> elemssd;
    StdVector<std::string> regionNames;
    StdVector<std::string> coupledRegionNames(1);


    // fetch combustion data node
    PtrParamNode combNode = myParam_->Get("combustionData");
    combNode->GetValue("coupledRegion",coupledRegionNames[0]);
    ptgrid_->GetRegion().Parse(regionNames, subdoms_);
    ptgrid_->GetRegion().Parse(coupledRegionNames, couplSubDomId_);

    std::cout << "CplRegion: " << coupledRegionNames[0] << std::endl;
    std::cout << "CplRegionID: " << couplSubDomId_ << std::endl;

    UInt varSpeedOfSound, srcReynoldStress, srcMomentum,
      srcChemical, srcDensityDeriv2, srcHeatRelease, srcShearTerm,
      srcGasConstDeriv2;

    combNode->GetValue("fileName",filenameCFD_);
    combNode->GetValue("varSpeedOfSound",varSpeedOfSound);
    combNode->GetValue("ReynoldStress",srcReynoldStress);
    combNode->GetValue("momentum",srcMomentum);
    combNode->GetValue("chemical",srcChemical);
    combNode->GetValue("densityDeriv2",srcDensityDeriv2);
    combNode->GetValue("heatRelease",srcHeatRelease);
    combNode->GetValue("gasConstDeriv2",srcGasConstDeriv2);
    combNode->GetValue("shearTerm",srcShearTerm);

    Info->WriteCombustionNoiseInfo(filenameCFD_, coupledRegionNames[0], 
				   varSpeedOfSound, srcReynoldStress, 
				   srcMomentum, srcChemical,  srcDensityDeriv2,
				   srcHeatRelease, srcGasConstDeriv2, srcShearTerm);

    // check for analysis
    AnalysisType analysisType = domain->GetSingleDriver()->GetAnalysisType();
    if ( analysisType==HARMONIC ) {
      isHarmonic_=true;
      ComputeRHSforHarm_=true;    
      Info->PrintF(pdename_, "Using FlowData from dataset as RHS nodal source\n" );
      Info->PrintF(pdename_, "Computing using nodal frequency files (No MpCCI used)\n" );
    }

    //get number of src-data in CFD-file
    combNode->GetValue("numDataInCFD", numDataInCFD_);

    UInt numSrc = srcReynoldStress + srcMomentum + srcChemical 
      + srcDensityDeriv2 + srcHeatRelease + srcGasConstDeriv2 + srcShearTerm;

    //check for source terms
    if ( numSrc == 1 ) {
      if ( srcDensityDeriv2 > 0 ) {
	srcType_ = DENSITY2DERIV;
      }
      else if ( srcReynoldStress > 0 ) {
	srcType_ = REYNOLDSTRESS;
      }
      else if ( srcMomentum  > 0 ) {
	srcType_ = MOMENTUM;
      }
      else if ( srcChemical > 0 ) {
	srcType_ = CHEMICAL;
      }
      else if ( srcHeatRelease > 0 ) {
	srcType_ = HEATRELEASE;
      }
      else if ( srcGasConstDeriv2 > 0 ) {
	srcType_ = GASCONSTDERIV2;
      }
      else if ( srcShearTerm > 0 ) {
	srcType_ = SHEARTERM;
      }
    }
    else if ( numSrc == 2 ) {
      if ( srcHeatRelease > 0 && srcGasConstDeriv2 > 0 ) {
         srcType_ = HEATRELEASE_GASCONSTDERIV2;
      }
    }
    else if ( numSrc == 3 ) {
      if ( srcReynoldStress > 0 && srcMomentum  > 0 && srcChemical > 0 ) { 
         srcType_ = REYNOLDSTRESS_MOMENTUM_CHEMICAL;
      }
      else if ( srcHeatRelease > 0 && srcGasConstDeriv2 > 0 && srcShearTerm > 0 ) {
         srcType_ = HEATRELEASE_GASCONSTDERIV2_SHEARTERM;
      }
    }
    else if ( numSrc == 0 && varSpeedOfSound == 1)  {
      srcType_ = NOSRC;
    }
    else {
      EXCEPTION("Check your xml-file according to the choice of your source terms");
    }

    // get number of nodes for coupled region!
    numNodesInCoupledRegion_ = ptgrid_->GetNumNodes(couplSubDomId_);

    // get memory for the CFD data
    dataCFD_ = new Double *[numDataInCFD_];

    for ( UInt i=0; i<numDataInCFD_; i++) {
      dataCFD_[i] = new Double[numNodesInCoupledRegion_];
    }
    //set aray to zero
    for ( UInt i=0; i<numDataInCFD_; i++) {
      for ( UInt j=0; j<numNodesInCoupledRegion_; j++) {
	dataCFD_[i][j] = 0.0;
      }
    }

    //get memory for node number
    nodesCFD_ = new Integer[numNodesInCoupledRegion_];


    //initialize the time step counter
    timeStep_ = 1;
  }

  AcouCombustionNoise::~AcouCombustionNoise()
  {
    if ( dataCFD_ != NULL ) {
      delete[]  dataCFD_[0];
      delete[]  dataCFD_;
    }
    if ( nodesCFD_ != NULL ) {
      delete[]  nodesCFD_;
    }
  }


  void AcouCombustionNoise::DefineSolveStep()
  {

    solveStep_ = new StdSolveStep(*this);
  }

  void AcouCombustionNoise::ComputeRHS(const Double atime)
  {

    Vector<Double> elemvec, nodalval;
    Matrix<Double> ptCoordNodes;
    StdVector<UInt> connecth;
    // For changing connecth to PDE
    StdVector<Integer> connect_PDE; 
    // unuased: BaseFE * ptEl;
    // unuased: UInt j;
    // unuased: UInt elsize = 0;
    StdVector<Elem*> elemssd;     
    // unuased: Double valmult;
    Vector<Double> valVec(1);

    // Get Data      
    ReadFlowData( filenameCFD_.c_str() );

    // Variables for ramping
    Double xfmin, yfmin, zfmin, xfmax, yfmax, zfmax, facRampXmin, facRampYmin, 
      facRampZmin, facRampXmax, facRampYmax, facRampZmax;
    Double bndoffsetXmin, bndoffsetYmin, bndoffsetZmin = 0.0, bndoffsetXmax, bndoffsetYmax, bndoffsetZmax = 0.0;

    // fetch combustion data node
    PtrParamNode combNode = myParam_->Get("combustionData");

    combNode->GetValue("xfmin",xfmin);
    combNode->GetValue("yfmin",yfmin);
    combNode->GetValue("xfmax",xfmax);
    combNode->GetValue("yfmax",yfmax);
    combNode->GetValue("facrampXmin",facRampXmin);
    combNode->GetValue("facrampYmin",facRampYmin);
    combNode->GetValue("facrampXmax",facRampXmax);
    combNode->GetValue("facrampYmax",facRampYmax);

    bndoffsetXmin=facRampXmin*xfmin;
    bndoffsetYmin=facRampYmin*yfmin; 
    bndoffsetXmax=facRampXmax*xfmax;
    bndoffsetYmax=facRampYmax*yfmax;

    if (dim_==3) {
      combNode->GetValue("zfmin",zfmin);
      combNode->GetValue("zfmax",zfmax);
      combNode->GetValue("facrampZmin",facRampZmin);
      combNode->GetValue("facrampZmax",facRampZmax);
      bndoffsetZmin=facRampZmin*zfmin;
      bndoffsetZmax=facRampZmax*zfmax;
    }
    
    UInt dof;
    Integer eqnNr;
    StdVector<UInt> connect(1);
    Double srcVal = 0.0;
    Double pointX, pointY, pointZ = 0.0;
    bool inBox;

    Vector<Double> sos(1);
    if ( variableSpeedOfSoundCN_ ) {
      // we will store here the variable speed of sound
      // so that it is available in the bilinearform
      speedOfSound_.Init();
    }

    dof=1;    
    Integer nodeNr;

    if (!isHarmonic_) {
      for (UInt idx=0; idx<numNodesInCoupledRegion_ ; idx++) {
	//get node number
        nodeNr = nodesCFD_[idx];
        //std::cout << "nodeNr: " << nodeNr << std::endl;

	// get the source term
	GetSrcTerm(srcVal, idx);
        //std::cout << "val = " << srcVal << std::endl;

        // TODO: Check if this is still needed
	// Double origVal = srcVal;

	//variable speed of sound
	if ( variableSpeedOfSoundCN_ ) {
	  sos[0] = dataCFD_[varSpeedOfSoundIdx_][idx];  
          //          std::cout << "sos = " << sos[0] << std::endl;
	}

	//	  node = idx + 1;
	  
	  // Ramping before adding to RHS vector
	  Matrix<Double> ptCoordNodes;
	  connecth.Resize(1);
	  connecth[0] = nodeNr;
	  ptgrid_->GetElemNodesCoord(ptCoordNodes, connecth);         
	  
	  pointX = ptCoordNodes[0][0];
	  pointY = ptCoordNodes[1][0];
	  if(dim_==3) {
	    pointZ = ptCoordNodes[2][0];
	  }

	  //check, if point is in specified (xml-file) box 
	  inBox = false;
	  if ( pointX > xfmin &&  pointX < xfmax ) {
	    if ( pointY > yfmin &&  pointY < yfmax ) {
	      inBox = true;
	      if ( dim_ == 3 ) {
		inBox = false;
		if ( pointZ > zfmin &&  pointZ < zfmax ) {
		  inBox = true;
		}
	      }
	    }
	  }
	  if ( !inBox ) {
	    srcVal = 0.0;
	  }
	  else {
	    // check for smoothing down
	    if ( pointX < bndoffsetXmin )
	      srcVal -= srcVal*(pointX-bndoffsetXmin)/(xfmin-bndoffsetXmin);
	    else
	      if ( pointX > bndoffsetXmax )
		srcVal -= srcVal*(pointX-bndoffsetXmax)/(xfmax-bndoffsetXmax);

	    if ( pointY < bndoffsetYmin )
	      srcVal  -= srcVal*(pointY-bndoffsetYmin)/(yfmin-bndoffsetYmin);
	    else      
	      if ( pointY > bndoffsetYmax )
		srcVal -= srcVal*(pointY-bndoffsetYmax)/(yfmax-bndoffsetYmax);
	    
	    if(dim_==3) {
	      if ( pointZ < bndoffsetZmin )
		srcVal -= srcVal*(pointZ-bndoffsetZmin)/(zfmin-bndoffsetZmin);
	      else      
		if ( pointZ > bndoffsetZmax )
		  srcVal -= srcVal*( pointZ-bndoffsetZmax)/(zfmax-bndoffsetZmax);
	    }
	  }

// 	  if ( abs(origVal-srcVal) > 1e-8 ) {
// 	    std::cout << "Orig. srcVal= " << origVal 
// 		      << "  Reduced srcVal= " << srcVal << std::endl;
// 	    std::cout << "Coords:\n" << ptCoordNodes << std::endl;
// 	  }

	  //add to RHS
	  eqnNr = eqnMap_->GetNodeEqn(nodeNr,dof );
	  //	  std::cout << "node: " << srcVal << std::endl;
	  algsys_->SetNodeRHS(srcVal, pdeId_, eqnNr);   

	  if ( variableSpeedOfSoundCN_ ) {
            //	    speedOfSound_.SetNodalResult(eqnNr, sos);
	    speedOfSound_.SetNodalResult(nodeNr, sos);
	  }
	  
	  if ( saveNodalSourcesRHS_ ) {
	    //		    Info->PrintSrcRhs(nodeNr, eqnNr , val);
	    valVec[0] = srcVal;
	    rhsNodalSrc_.SetNodalResult(nodeNr, valVec);
            //	    rhsNodalSrc_.SetNodalResult(eqnNr, valVec);
	  }   
	}
      
    } 
    else  {
      std::cout << "Using frequency source files..."<< std::endl;
      EXCEPTION("Currently harmonic case not implemented!");
    }    
    
    timeStep_ += 1;
  } 



  void AcouCombustionNoise::ReadFlowData(const char * aname)
  {
    
    std::string filename;
    char buf[128];
        
    filename = aname;
    filename+= "_01_";
    sprintf(buf, "%04i", timeStep_);
    filename+= buf;
    filename+= ".dat";

    std::cout<<"name of current CFD-File: "<< filename << std::endl;
    read_in(filename.c_str(), dataCFD_, numNodesInCoupledRegion_,  numDataInCFD_, nodesCFD_);

//     std::cout << "numData: " << numDataInCFD_ << std::endl;
//     for ( UInt i=0; i<numNodesInCoupledRegion_; i++ ) {
//       for ( UInt j=0; j< numDataInCFD_; j++ ) {
//         std::cout << dataCFD_[j][i] << "   ";
//       }
//       std::cout << std::endl;
//     }    

}

  void AcouCombustionNoise::GetSrcTerm( Double& val, UInt idx) {
 
  if ( srcType_ ==  REYNOLDSTRESS_MOMENTUM_CHEMICAL ) {
    val = dataCFD_[srcReynoldStressIdx_][idx] 
        + dataCFD_[srcMomentumIdx_][idx]
        + dataCFD_[srcChemicalIdx_][idx];
//     std::cout << "values: " << dataCFD_[srcReynoldStressIdx_][idx] 
// 	      << "  " << dataCFD_[srcMomentumIdx_][idx] 
// 	      << "  " <<  dataCFD_[srcChemicalIdx_][idx] << std::endl;
  }
  else if ( srcType_ ==  HEATRELEASE_GASCONSTDERIV2_SHEARTERM ) {
    val = dataCFD_[srcHeatReleaseIdx_][idx]
        + dataCFD_[srcGasConstDeriv2Idx_][idx]
        + dataCFD_[srcShearTermIdx_][idx];
  }
  else if ( srcType_ ==  HEATRELEASE_GASCONSTDERIV2 ) {
    val = dataCFD_[srcHeatReleaseIdx_][idx]
        + dataCFD_[srcGasConstDeriv2Idx_][idx];
  }
  else if ( srcType_ == DENSITY2DERIV ) {
    val = dataCFD_[srcDensityDeriv2Idx_][idx];
  }
  else if ( srcType_ == REYNOLDSTRESS ) {
    val = dataCFD_[srcReynoldStressIdx_][idx]; 
  }
  else if ( srcType_ == MOMENTUM ) {
    val = dataCFD_[srcMomentumIdx_][idx]; 
  }
  else if ( srcType_ == CHEMICAL ) {
    val = dataCFD_[srcChemicalIdx_][idx]; 
  }
  else if ( srcType_ ==  HEATRELEASE ) {
    val = dataCFD_[srcHeatReleaseIdx_][idx]; 
  }
  else if ( srcType_ == GASCONSTDERIV2 ) {
    val = dataCFD_[srcGasConstDeriv2Idx_][idx]; 
  }
  else if ( srcType_ ==  SHEARTERM ) {
    val = dataCFD_[srcShearTermIdx_][idx]; 
  }

  else if ( srcType_ == NOSRC ) {
     val = 0.0;
  }

 }


} // end of namespace
