#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "Domain/elem.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/WriteInfo.hh"
#include "PDE/StdPDE.hh"

#include "olas.hh"

#include "assemble.hh"

namespace CoupledField {
  
  Assemble::Assemble(BaseSystem * algsys, Grid * aptgrid)
    :algsys_(algsys),
     ptgrid_(aptgrid),
     stiffnessMatrix_(FALSE),
     dampingMatrix_(FALSE),
     massMatrix_(FALSE),
     convectionMatrix_(FALSE),
     actlevel_(0),
     integrators_(0),
     rhsIntegrators_(0),
     rhsSrcIntegrators_(0)
  {
    ENTER_FCN( "Assemble::Assemble" );
    firstTime_ = TRUE;
    oneIntIsNonlin_ = FALSE;
    nrMatrices_ = 5+1;  // 5 matrices, but index starts with 1!
    reassembleMat_.Resize(nrMatrices_);
    // reassembleMat_.resize(nrMatrices_);
    nonLinGeo = FALSE;
    deltaCoords_ = NULL;
    alternateMaterialData_=FALSE;

    olasParams_ = algsys_->GetOLASParams();
    olasReport_ = algsys_->GetOLASReport(); 

    startFreq_ = 0.0;

  }


  Assemble::~Assemble()
  {
    ENTER_FCN( "Assemble::~Assemble" );
    
    for (int i=0; i<integrators_.GetSize();i++)
      for (int j=0; j<integrators_[i]->GetSize(); j++)
        delete (*integrators_[i])[j];

    for (int i=0; i<surfintegrators_.GetSize();i++)
      for (int j=0; j<surfintegrators_[i]->GetSize(); j++)
        delete (*surfintegrators_[i])[j];

    for (int i=0; i<rhsIntegrators_.GetSize();i++)
      for (int j=0; j<rhsIntegrators_[i]->GetSize(); j++)
        delete (*rhsIntegrators_[i])[j];
  }

  void Assemble::SetPtr2EQNData(NodeEQN * aPtNodeEQN)
  {
    ENTER_FCN( "Assemble::SetPtr2EQNData" );
    if (aPtNodeEQN == NULL)
      Error( "The NodeEQN object is not created yet" , __FILE__, __LINE__);
    
    if (! aPtNodeEQN->IsInitialized())
      Error( "The NodeEQN object has to be initialized before assigning it to\
              an assemble object by calling 'CalcMapping()'", __FILE__, __LINE__);
    ptEQN_ = aPtNodeEQN;
  }


  Integer Assemble::SubDomIndex(const std::string & subDomName)
  {
    ENTER_FCN( "Assemble::SubDomIndex" );

    for (int i=0; i < subdoms_.GetSize();i++)
      {
        if (subDomName == subdoms_[i]) return i;
      }
  
    std::string errOut;
    errOut = "SubDomain " + subDomName + " not defined!";
    Info->Error(errOut, __FILE__, __LINE__);

    return 0;
  }
  
  Integer Assemble::SurfDomIndex(const std::string & surfDomName)
  {
    ENTER_FCN( "Assemble::SurfDomIndex" );

    for (int i=0; i < surfdoms_.GetSize();i++)
      if (surfDomName == surfdoms_[i])
        return i;
    
  
    std::string errOut;
    errOut = "Surface-Domain " + surfDomName + " not defined!";
    Info->Error(errOut, __FILE__, __LINE__);
    return -1;
  }


  void Assemble::GetElemCoords(const StdVector<Integer> connect, 
                               Matrix<Double> &coordMat, const Integer level)
  {
    ENTER_FCN( "Assemble:GetElemCoords" );
    
    ptgrid_->GetCoordNodesElemMat(connect, coordMat, level);
  
    if (nonLinGeo)
      {
        if (deltaCoords_ == NULL)
          Error("ElecPDE: set input_coupling_terms = smoothdisplacement or nonlin = no");

        StdVector<Integer> connect_PDE;
        ptEQN_->Mesh2PDENode(connect_PDE, connect);
        Double val;
        for (Integer i=0; i<coordMat.GetSizeRow(); i++)
          for (Integer j=0; j<coordMat.GetSizeCol(); j++) 
            {
              val = (*deltaCoords_)[i][connect_PDE[j] - 1];
              coordMat(i,j) += val;
            }
      }
  }

  
  // do the basic assembling stuff
  void Assemble::AssembleMatrices(const Integer level)
  {

    ENTER_FCN( "Assemble:AssembleMatrices" );
    
    Vector<Double> harmonicVec;
    Double dampTransform = 1.0;

    // initialize reassembling "indicator" vector
    if (firstTime_)
      for (Integer actMat=0; actMat < nrMatrices_; actMat++)
        reassembleMat_[actMat] = FALSE;

    
    for (Integer actDom=0; actDom < subdoms_.GetSize(); actDom++) {     
      StdVector<Elem*> elemssd;

      ptgrid_->GetElemSD(elemssd, subdoms_[actDom], level);

      //        std::cout<<"\n we have Number ofIntegrators= " << integrators_[actDom]->GetSize() <<std::endl;

      for(Integer actInteg=0; actInteg < integrators_[actDom]->GetSize(); actInteg++) {
        IntegratorDescriptor * actDescriptor = (*integrators_[actDom])[actInteg];

        if (alternateMaterialData_ == TRUE)
          actDescriptor->GetIntegrator()->SetMaterial(ptMaterial_);
                    
        // assemble only if nonlinear or first time
        if (reassembleMat_[actDescriptor->DestMat()] || firstTime_) {

          if ( actDescriptor->IsReducedInt()) {
            //set all elements to reduced integration!!
            SetFE2ReducedInt();
          }

          if ( actDescriptor->GetIntegrator()->IsFracDamping() ) {
            // get multiplicative pre factor depending on time step size
            Double damp;
            damp = ptPDE_->GetFracDampMatrixCoeff(actDom);
            actDescriptor->GetIntegrator()->SetFactor(damp);
          }

	  dampTransform = 1.0;
	  if ( ( actDescriptor->GetIntegrator()->IsRaylDamping() 
		 || actDescriptor->GetSecondaryMat() != NOTYPE ) 
	       && startFreq_ > 0 ) {
            // get multiplicative pre factor depending on time step size
	    if (startFreq_ > 0 && actFreq_ > 0 ) {
	      FEMatrixType destMat = actDescriptor->GetIntegrator()->GetBaseType();
	      if ( destMat == STIFFNESS) {
		dampTransform = startFreq_ / actFreq_;
	      }
	      else if ( destMat == MASS ) {
		dampTransform = actFreq_ / startFreq_;
	      }
	    }
	    
	    if ( actDescriptor->GetIntegrator()->IsRaylDamping() ) {
	      actDescriptor->GetIntegrator()->SetFactor(dampTransform);
	    }
	    std::cout << "startF=" << startFreq_ << "  actF=" << actFreq_ << std::endl;
	    std::cout << "dampTransform=" << dampTransform << std::endl;
          }
	  

          
          for (Integer actEl=0; actEl< elemssd.GetSize(); actEl++) {
            BaseFE * ptEl = elemssd[actEl]->ptElem;
            StdVector<Integer> connecth = elemssd[actEl]->connect;
                             
            Matrix<Double> ptCoord;
            GetElemCoords(connecth, ptCoord, level);
                    
            // map connect to PDE node numbers
            StdVector<Integer> connect_PDE;
                    
            ptEQN_->Node2EQN(connecth, connect_PDE);
                    
            Matrix<Double> elSol;
                    

            // this matrix is nonlinear and, therefore, has to be reassembled next time
            if ((oneIntIsNonlin_ || firstTime_ ) && analysisType_ != HARMONIC)
              // fetch solution at element nodes
              sol_->GetElemSolutionAsMatrix(elSol, connecth);
            
            actDescriptor->GetIntegrator()->SetElemPtr(ptEl);
            FEMatrixType destMat = actDescriptor->DestMat();
                
            // this matrix is nonlinear and, therefore, has to be reassembled next time
            if (actDescriptor->IsNonLin()) {
              oneIntIsNonlin_ = TRUE;
              reassembleMat_[actDescriptor->DestMat()] = TRUE;
              
              actDescriptor->GetIntegrator()->SetActElemSol(elSol);
            }       

            // ================================================================
            //                             assemble matrices
            // ================================================================

            actDescriptor->GetIntegrator()->CalcElementMatrix(ptCoord, elemmat);
                  
            piezoMaterialType matType = actDescriptor->GetPiezoMaterialType();
            //std::cout<<"\n assemble:: piezoMaterialType = " << matType <<std::endl;
            actDescriptor->SetPiezoMaterialType(matType);


            if (analysisType_ == HARMONIC) {
              TransformMatrix2Harmonic(harmonicVec,elemmat, actDescriptor->GetOrigMatrixType(),
                                       actDescriptor->GetPiezoMaterialType());
         
              algsys_->SetElementMatrix(&harmonicVec[0], connect_PDE.GetPointer(), connect_PDE.GetSize(), destMat);
            }
            else {
	      algsys_->SetElementMatrix(elemmat.GetDataPointer(), connect_PDE.GetPointer(), 
                                        connect_PDE.GetSize(), destMat);
	    }
#ifdef DEBUG
            // output matrices
            if (destMat == STIFFNESS)
              (*debug) << "Stiffness matrix of Element " << actEl << std::endl;
            if (destMat == MASS)
              (*debug) << "Mass      matrix of Element " << actEl << std::endl;
            if (destMat == DAMPING)
              (*debug) << "Damping   matrix of Element " << actEl << std::endl;
            (*debug) << elemmat << std::endl;
            if ( !elemmat.IsSymmetric() ) {
              (*debug) << " --> Matrix is not symmetric " << std::endl << std::endl;
            }
            else {
              (*debug) << " --> Matrix is symmetric " << std::endl << std::endl;
            }
#endif
            if (actDescriptor->GetSecondaryMat() != NOTYPE) {
	      Double damp = dampTransform * actDescriptor->GetSecMatFac();
              elemmat *= damp;
              if (analysisType_ == HARMONIC) {
                TransformMatrix2Harmonic(harmonicVec,elemmat,actDescriptor->GetOrigSecMatrixType(),
                                         actDescriptor->GetPiezoMaterialType());
                algsys_->SetElementMatrix(&harmonicVec[0], connect_PDE.GetPointer(),connect_PDE.GetSize(), destMat);
              }
              else
                algsys_->SetElementMatrix(elemmat.GetDataPointer(), connect_PDE.GetPointer(), 
                                          connect_PDE.GetSize(), actDescriptor->GetSecondaryMat()); 
            }
                  
          } //over all elements of subdomain            
                
        } //check, if we have to assemble
          
        if ( actDescriptor->IsReducedInt()) {
          //set all elements back to standard integration!!
          SetFE2StandardInt();
        }
          

      } //integrators
        
    } //subdomains
  
    //assemble matrices for surface integrators
    for (Integer actDom=0; actDom < surfdoms_.GetSize(); actDom++) {    

      //check is necessary, because the surface-integrator could also be a RHS-Src integrator
      if (surfintegrators_[actDom]->GetSize()) {        
        StdVector<Elem*> elemssd;
        elemssd=ptBCs_->getFacesBC(surfdoms_[actDom],level);
          
        for (Integer actEl=0; actEl< elemssd.GetSize(); actEl++) {
          BaseFE * ptEl = elemssd[actEl]->ptElem;
          StdVector<Integer> connecth = elemssd[actEl]->connect;

          Matrix<Double> ptCoord;
          GetElemCoords(connecth, ptCoord, level);

          // map connect to PDE node numbers
          StdVector<Integer> connect_PDE;
          ptEQN_->Node2EQN(connecth, connect_PDE);
              
          Matrix<Double> elSol;
              
          // this matrix is nonlinear and, therefore, has to be reassembled next time

          if (oneIntIsNonlin_ || firstTime_)
            // fetch solution at element nodes
            //GetSolOfElement(elSol, connect_PDE);
            sol_->GetElemSolutionAsMatrix(elSol, connecth);
                
          // ================================================================
          //                             assemble matrices
          // ================================================================
                
          for(Integer actInteg=0; actInteg < surfintegrators_[actDom]->GetSize(); actInteg++) {
            IntegratorDescriptor * actDescriptor = (*surfintegrators_[actDom])[actInteg];

            // assemble only if nonlinear or first time
            if (reassembleMat_[actDescriptor->DestMat()] || firstTime_) {
              actDescriptor->GetIntegrator()->SetElemPtr(ptEl);
              FEMatrixType destMat = actDescriptor->DestMat();
              // this matrix is nonlinear and, therefore, has to be reassembled next time
              if (actDescriptor->IsNonLin()) {
                oneIntIsNonlin_ = TRUE;
                reassembleMat_[actDescriptor->DestMat()] = TRUE;
                actDescriptor->GetIntegrator()->SetActElemSol(elSol);
              }
                      
              actDescriptor->GetIntegrator()->CalcElementMatrix(ptCoord, elemmat);

              if (analysisType_ ==HARMONIC) {
                TransformMatrix2Harmonic(harmonicVec,elemmat, actDescriptor->GetOrigMatrixType(),
                                         actDescriptor->GetPiezoMaterialType());
                algsys_->SetElementMatrix(&harmonicVec[0], connect_PDE.GetPointer(),connect_PDE.GetSize(), destMat);
              }

              else
                algsys_->SetElementMatrix(elemmat.GetDataPointer(), connect_PDE.GetPointer(), connect_PDE.GetSize(), destMat);
                        
              if (actDescriptor->GetSecondaryMat()  != NOTYPE ) {
                elemmat *= actDescriptor->GetSecMatFac();
                if (analysisType_ == HARMONIC) {
                  TransformMatrix2Harmonic(harmonicVec,elemmat,actDescriptor->GetOrigSecMatrixType(),
                                           actDescriptor->GetPiezoMaterialType());
                  algsys_->SetElementMatrix(&harmonicVec[0], connect_PDE.GetPointer(), connect_PDE.GetSize(), destMat);
                }
                else
                  algsys_->SetElementMatrix(elemmat.GetDataPointer(), connect_PDE.GetPointer(), 
                                            connect_PDE.GetSize(), actDescriptor->GetSecondaryMat()); 
              }
            }           
          }
        } // elements
      } // check for surface integrator
    } // subdomains
     
  firstTime_ = FALSE;

}


  // do the basic assembling stuff
  void Assemble::AssembleSrcRHS(const Integer level, const Double time)
  {
    ENTER_FCN( "Assemble:AssembleSrcRHS" );
    AssembleRHSNodalSources(level, time);
    AssembleRHSIntegralSources(level, time);
  }


  // do the basic assembling stuff
  void Assemble::AssembleRHSIntegralSources(const Integer level,const Double time)
  {
    ENTER_FCN( "Assemble:AssembleRHSIntegralSources" );

    Vector<Double> harmVec;

    //add the volume sources
    for (Integer actDom=0; actDom <  subdoms_.GetSize(); actDom++)
      { 
        if (rhsSrcIntegrators_[actDom]->GetSize())
          {
            StdVector<Elem*> elemssd;
            ptgrid_->GetElemSD(elemssd, subdoms_[actDom], level);
            
            Double val_tfunc = 1.0;
            Double valPhase = 0.0;
            if (analysisType_ == HARMONIC) 
              valPhase = rhsSrcPhase_[actDom];
            else {
              if (ptTimeFunc_->GetmaxTimeFnc() > 0 )
                val_tfunc=ptTimeFunc_->TimeFuncAtTime(time,fncname_rhs_[actDom]);
            }

            for (Integer actEl=0; actEl< elemssd.GetSize(); actEl++)
              {        
                BaseFE * ptEl = elemssd[actEl]->ptElem;
                StdVector<Integer> connecth = elemssd[actEl]->connect;
                
                Matrix<Double> ptCoord;
                GetElemCoords(connecth, ptCoord, level);
            
                // map connect to PDE node numbers
                StdVector<Integer> connect_PDE;
                ptEQN_->Node2EQN(connecth, connect_PDE);
                for(Integer actRhsInt=0; actRhsInt < rhsSrcIntegrators_[actDom]->GetSize(); actRhsInt++)
                  {
                    BaseIntDescriptor * actRhsID = (*rhsSrcIntegrators_[actDom])[actRhsInt];
                    
                    actRhsID->GetIntegrator()->SetElemPtr(ptEl);
                    
                    Vector<Double> elemVec;
                    actRhsID->GetIntegrator()->CalcElemVector(ptCoord, elemVec);
                    
                    if (analysisType_ == HARMONIC) {
                      TransformVector2Harmonic(harmVec,elemVec,valPhase);
                      algsys_->SetElementRHS(&harmVec[0], connect_PDE.GetPointer(), connect_PDE.GetSize());
                    }
                    else {
                      if (val_tfunc != 1.0)
                        elemVec *= val_tfunc;
                      algsys_->SetElementRHS(&elemVec[0], connect_PDE.GetPointer(), connect_PDE.GetSize());
                    }
                  }
              }
          }
      }

    //add the surface sources
    for (Integer actSurf=0; actSurf <  surfdoms_.GetSize(); actSurf++)
      { 
        if (rhsSrcSurfIntegrators_[actSurf]->GetSize())
          {
            StdVector<Elem*> elemssd;
            elemssd=ptBCs_->getFacesBC(surfdoms_[actSurf],level);
            
            Double val_tfunc = 1.0;
            Double valPhase = 0.0;
            if (analysisType_ == HARMONIC) 
              valPhase = rhsSrcSurfPhase_[actSurf];
            else {
              if (ptTimeFunc_->GetmaxTimeFnc() > 0 )
                val_tfunc=ptTimeFunc_->TimeFuncAtTime(time,fncname_rhsSurf_[actSurf]);
            }

            for (Integer actEl=0; actEl< elemssd.GetSize(); actEl++)
              {        
                BaseFE * ptEl = elemssd[actEl]->ptElem;
                StdVector<Integer> connecth = elemssd[actEl]->connect;
                
                Matrix<Double> ptCoord;
                GetElemCoords(connecth, ptCoord, level);
            
                // map connect to PDE node numbers
                StdVector<Integer> connect_PDE;
                ptEQN_->Node2EQN(connecth, connect_PDE);
                for(Integer actRhsInt=0; actRhsInt < rhsSrcSurfIntegrators_[actSurf]->GetSize(); actRhsInt++)
                  {
                    BaseIntDescriptor * actRhsID = (*rhsSrcSurfIntegrators_[actSurf])[actRhsInt];
                    
                    actRhsID->GetIntegrator()->SetElemPtr(ptEl);
                    
                    Vector<Double> elemVec;
                    actRhsID->GetIntegrator()->CalcElemVector(ptCoord, elemVec);
                    if (analysisType_ == HARMONIC) {
                      TransformVector2Harmonic(harmVec,elemVec,valPhase);
                      algsys_->SetElementRHS(&harmVec[0], connect_PDE.GetPointer(), 
                                             connect_PDE.GetSize());
                    }
                    else {
                      if (val_tfunc != 1.0)
                        elemVec *= val_tfunc;
                      algsys_->SetElementRHS(&elemVec[0], connect_PDE.GetPointer(), connect_PDE.GetSize());
                    }
                  }
              }
          }
      }

  }



  // do the basic assembling stuff
  void Assemble::AssembleRHSNodalSources(const Integer level, const Double time)
  {
    ENTER_FCN( "Assemble:AssembleRHSNodalSources" );
    
    Integer eqnNr, eqnDof;
    
    for (int actDom=0; actDom < loadDom_.GetSize(); actDom++)
      {
        std::string doftype = loadDom_[actDom];

        Integer dof = 1;
        if (dofsPerNode_ != 1)
          dof = GetBCDof( loadDof_[actDom] );   

        std::list<Integer> nodes;
        nodes = ptBCs_->GetNodesLevel(loadDom_[actDom], level);
        
        Double val = loadVals_[actDom];

        Double val_tfunc = 1.0;
        if (ptTimeFunc_->GetmaxTimeFnc() > 0 )
          val_tfunc=ptTimeFunc_->TimeFuncAtTime(time,fncname_loads_[actDom]);


        for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++)
          {
            Integer node = *p;
            
            val = loadVals_[actDom] * val_tfunc;

            ptEQN_->Node2EQN(node,dof,eqnNr,eqnDof);
            algsys_->SetNodeRHS(val, eqnNr, eqnDof);    
          }
      }
  }
  






  // do the basic assembling stuff
  void Assemble::AssembleNLRHS(const Integer level, const Double time)
  {
    ENTER_FCN( "Assemble:AssembleNLRHS" );

    Matrix<Double> elemmat;

    StdVector<Elem*> elemssd;
   
    for (int actDom=0; actDom < subdoms_.GetSize(); actDom++)
      { 
        if (rhsIntegrators_[actDom]->GetSize())
          {
            ptgrid_->GetElemSD(elemssd, subdoms_[actDom], level);
        
            for (int actEl=0; actEl< elemssd.GetSize(); actEl++)
              {
                BaseFE * ptEl = elemssd[actEl]->ptElem;
                StdVector<Integer> connecth = elemssd[actEl]->connect;


                Matrix<Double> ptCoord;
                GetElemCoords(connecth, ptCoord, level);


                // map connect to PDE node numbers
                StdVector<Integer> connect_PDE;
                ptEQN_->Node2EQN(connecth, connect_PDE);
                
                Matrix<Double> elSol;
                
                sol_->GetElemSolutionAsMatrix(elSol, connecth);
                
                // ================================================================
                //                             assemble RHS
                // ================================================================
                
                for(Integer actRhsInt=0; actRhsInt < rhsIntegrators_[actDom]->GetSize(); actRhsInt++)
                  {
                    BaseIntDescriptor * actRhsID = (*rhsIntegrators_[actDom])[actRhsInt];
                    
                    actRhsID->GetIntegrator()->SetElemPtr(ptEl);
                    
                    if (actRhsID->IsNonLin())
                      actRhsID->GetIntegrator()->SetActElemSol(elSol);
                    
                    
                    Vector<Double> elemVec;
                    actRhsID->GetIntegrator()->CalcElemVector(ptCoord, elemVec);
                    
                    algsys_->SetElementRHS(&elemVec[0], connect_PDE.GetPointer(), connect_PDE.GetSize());
                  }
              }
          }
      }
  }
  

  // do the basic assembling stuff
  void Assemble::AssembleSprings(const Integer level, const Double time)
  {
    ENTER_FCN( "Assemble::AssembleSprings" );
    
    Integer eqnNr, eqnDof;
    
    for (int actDom=0; actDom < springDom_.GetSize(); actDom++)
      {
        std::string doftype = springDom_[actDom];

        Integer dof = 1;
        if (dofsPerNode_ != 1)
          dof = GetBCDof( springDof_[actDom] );   

        std::list<Integer> nodes;
        nodes = ptBCs_->GetNodesLevel(springDom_[actDom], level);
        
        Double val_m = springMassVals_[actDom];
        Double val_c = springDampVals_[actDom];
        Double val_k = springStiffVals_[actDom];

        Double val_tfunc = 1.0;
        if (ptTimeFunc_->GetmaxTimeFnc() > 0 )
          val_tfunc=ptTimeFunc_->TimeFuncAtTime(time,fncname_springs_[actDom]);


        for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++)
          {
            Integer node = *p;
            
            val_m = springMassVals_[actDom]  * val_tfunc;
            val_c = springDampVals_[actDom]  * val_tfunc;
            val_k = springStiffVals_[actDom] * val_tfunc;
	    std::cout << "val_m=" << val_m << std::endl;
	    std::cout << "val_c=" << val_c << std::endl;
	    std::cout << "val_k=" << val_k << std::endl << std::endl;
            ptEQN_->Node2EQN(node,dof,eqnNr,eqnDof);
            //algsys_->AddToSystemMatrix(MASS, val_m, eqnNr, eqnDof);    
            //algsys_->AddToSystemMatrix(DAMPING, val_c, eqnNr, eqnDof);    
            //algsys_->AddToSystemMatrix(STIFFNESS,val_k, eqnNr, eqnDof);    
          }
      }
  }



  Integer Assemble::
  GetBCDof(const std::string dofString)
  {
    ENTER_FCN( "MechPDE::GetBCDof" );

    if (dofString == "ux")
      return 1;
    else
      if (dofString == "uy")
        return 2;
      else
        if (dofString == "uz")
          return 3;
        else {
          Error("The direction mentioned in the config-file is not implemented! ",__FILE__,__LINE__);
          return -1;
        }
  }




  void Assemble::InitMatrices()
  {
    ENTER_FCN( "Assemble::InitMatrices" );

    //set firstTime_ to TRUE, so that assembling of element matrices will be preformed
    firstTime_ = TRUE;

    // Initialize matrices in order to get BCs correct
    algsys_->InitMatrix(SYSTEM);

    if (stiffnessMatrix_)
      algsys_->InitMatrix(STIFFNESS);
    
    if (dampingMatrix_)
      algsys_->InitMatrix(DAMPING);

    if (convectionMatrix_)
      algsys_->InitMatrix(CONVECTION);
    
    if (massMatrix_)
      algsys_->InitMatrix(MASS); 
  }
 


  void Assemble::InitNonLinMatrices()
  {
    ENTER_FCN( "Assemble::InitNonLinMatrices" );

    // return, if matrices are not yet assembled
    if (!reassembleMat_.GetSize()) {
      // if ( reassembleMat_.GetSize() == 0 ) {
      InitMatrices();
      return;
    }
    
    // Initialize matrices in order to get BCs correct
    algsys_->InitMatrix(SYSTEM);

    if (stiffnessMatrix_ && reassembleMat_[STIFFNESS]) 
      algsys_->InitMatrix(STIFFNESS);
    
    
    if (dampingMatrix_ && reassembleMat_[DAMPING])
      algsys_->InitMatrix(DAMPING);

    if (convectionMatrix_ && reassembleMat_[CONVECTION])
      algsys_->InitMatrix(CONVECTION);
    
    if (massMatrix_ && reassembleMat_[MASS]) 
      algsys_->InitMatrix(MASS); 
    
  }
 


 

  void Assemble::CreateMatrices()
  {
    ENTER_FCN( "Assemble::CreateMatrices" );
    const Integer numconstraints = 0;  // currently not handled
    
    const Integer dofsPerEQN = ptEQN_->GetNumDofsPerEQN();

    FEMatrixType matrixsystype[5];
    matrixsystype[0] = OLAS::NOTYPE;
    matrixsystype[1] = OLAS::NOTYPE;
    matrixsystype[2] = OLAS::NOTYPE;
    matrixsystype[3] = OLAS::NOTYPE;
    matrixsystype[4] = OLAS::NOTYPE;
    
    matrixsystype[0] = SYSTEM;      // memory for the system matrix
    if (stiffnessMatrix_  == 1) matrixsystype[1] = STIFFNESS;   // memory for the stiffness matrix
    if (dampingMatrix_    == 1) matrixsystype[2] = DAMPING;     // memory for the damping matrix
    if (convectionMatrix_ == 1) matrixsystype[3] = CONVECTION;  // memory for the convection matrix
    if (massMatrix_       == 1) matrixsystype[4] = MASS;        // memory for the mass matrix
    
    Integer numBuildInDirichletEQNs_ = ptEQN_->GetNumBuildInDirichletEQNs();
    Integer numDir = numDirichletBCs_ - numBuildInDirichletEQNs_;
    //std::cerr << "number of Dirichlet to set by AlgSYS: " << numDir << std::endl;
    //put to algebraic system
    
    olasParams_->SetValue( "FEMatrixType1", matrixsystype[0] );
    olasParams_->SetValue( "FEMatrixType2", matrixsystype[1] );
    olasParams_->SetValue( "FEMatrixType3", matrixsystype[2] );
    olasParams_->SetValue( "FEMatrixType4", matrixsystype[3] );
    olasParams_->SetValue( "FEMatrixType5", matrixsystype[4] );
    olasParams_->SetValue( "NumDof", dofsPerEQN);
    olasParams_->SetValue( "NumDirichletBCs", numDir);
    olasParams_->SetValue( "NumConstraints", numconstraints );
    olasParams_->SetValue( "AuxiliaryMatrix", FALSE);

    // Create linear system. This involves finalising the assembly of the
    // matrix graph. The latter may trigger a re-ordering of the unknows /
    // equation numbers. We provide an array for giving back the new
    // ordering, so that we can incorporate it into our node -> eqn map
    Integer *order = new Integer[ptEQN_->GetNumEQNs()];
    algsys_->CreateLinSys(order);
    ptEQN_->ReorderMapping(order);
    delete [] order;
    
  }
  

  void Assemble::SetGeneralParams(const std::string & pdename, 
                                  const Integer dofsPerNode,
                                  const Integer numPDENodes, 
                                  const StdVector<std::string> subdoms,
                                  const StdVector<std::string> surfdoms,
                                  const std::string bcSequenceTag)
  {
    ENTER_FCN( "Assemble::SetGeneralParams" );

    pdename_       = pdename;
    dofsPerNode_   = dofsPerNode;
    numPDENodes_   = numPDENodes;
    subdoms_       = subdoms;
    surfdoms_      = surfdoms;
    bcSequenceTag_ = bcSequenceTag;
 


    // parameters for reading load values =========================================
    StdVector<std::string> keyVec, attrVec, valVec;

    // parameters for reading spring values =======================================
    StdVector<std::string> keyVecSpring;

    attrVec = "", "tag", "";
    valVec = "", bcSequenceTag_, "";

    loadDom_.Clear();
    loadDof_.Clear();
    loadVals_.Clear();
    fncname_loads_.Clear();

    keyVec = pdename_, "bcsAndLoads", "load", "name";
    params->GetList(keyVec, attrVec, valVec, loadDom_);

    keyVecSpring = pdename_, "bcsAndLoads", "spring", "name";
    params->GetList(keyVecSpring, attrVec, valVec, springDom_);

    if (dofsPerNode_ == 1) {
      loadDof_.Resize(loadDom_.GetSize());
      loadDof_.Init("ux");

      springDof_.Resize(loadDom_.GetSize());
      springDof_.Init("ux");
    }else {
      keyVec = pdename_, "bcsAndLoads", "load", "dof";
      params->GetList(keyVec, attrVec, valVec, loadDof_);

      keyVecSpring = pdename_, "bcsAndLoads", "spring", "dof";
      params->GetList(keyVecSpring, attrVec, valVec, springDof_);
    }
    
    keyVec = pdename_, "bcsAndLoads", "load", "value";
    params->GetList(keyVec, attrVec, valVec, loadVals_);

    keyVecSpring = pdename_, "bcsAndLoads", "spring", "val_m";
    params->GetList(keyVecSpring, attrVec, valVec, springMassVals_);

    keyVecSpring = pdename_, "bcsAndLoads", "spring", "val_c";
    params->GetList(keyVecSpring, attrVec, valVec, springDampVals_);

    keyVecSpring = pdename_, "bcsAndLoads", "spring", "val_k";
    params->GetList(keyVecSpring, attrVec, valVec, springStiffVals_);

    keyVec = pdename_, "bcsAndLoads", "load", "dynamics";
    params->GetList(keyVec, attrVec, valVec, fncname_loads_);

    keyVecSpring = pdename_, "bcsAndLoads", "spring", "dynamics";
    params->GetList(keyVecSpring, attrVec, valVec, fncname_springs_);


    // check if a fncname for load was given, although
    // in section transient there is non mentioned
    if (fncname_loads_.GetSize() > 0)
      {
        if (ptTimeFunc_->GetmaxTimeFnc() == 0  &&
            fncname_loads_[0] != "none") {
          std::string errmsg = "Loads: There was no time data file ";
          errmsg += "specified in the section 'transient', so 'dynamics' ";
          errmsg += "for PDE '";
          errmsg += pdename;
          errmsg += "' makes no sense!";
          Error(errmsg.c_str(), __FILE__, __LINE__);
        }
      }

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // this section must be changed so that it is possible to
    // specify two or more timeDataFile in section transient
    // and to check them here
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // check if a fncname for spring was given, although
    // in section transient there is non mentioned
    if (fncname_springs_.GetSize() > 0)
      {
        if (ptTimeFunc_->GetmaxTimeFnc() == 0  &&
            fncname_springs_[0] != "none") {
          std::string errmsg = "Springs: There was no time data file ";
          errmsg += "specified in the section 'transient', so 'dynamics' ";
          errmsg += "for PDE '";
          errmsg += pdename;
          errmsg += "' makes no sense!";
          Error(errmsg.c_str(), __FILE__, __LINE__);
        }
      }
    

    //     std::cerr << "load names = " <<  loadDom_ << std::endl;
    //     std::cerr << "load dofs = " <<  loadDof_ << std::endl;
    //     std::cerr << "load values = " <<  loadVals_ << std::endl;
    //     std::cerr << "load dynamics = " <<  fncname_loads_ << std::endl;

    // Check load consistency
    if ( loadDom_.GetSize() != loadDof_.GetSize() ||
         loadDom_.GetSize() != loadVals_.GetSize() )
      {
        std::string errmsg = "Loads: ";
        errmsg += "#name = " + Info->GenStr(loadDom_.GetSize());
        errmsg += ", #dof = " + Info->GenStr(loadDof_.GetSize());
        errmsg += ", #value = " + Info->GenStr(loadVals_.GetSize());
        errmsg += ", #dynamics = " + fncname_loads_.GetSize() + '\n';
        Info->Error( errmsg, __FILE__, __LINE__ );
      }

    // Check spring consistency
    if ( springDom_.GetSize() != springDof_.GetSize() ||
         springDom_.GetSize() != springMassVals_.GetSize() ||
         springDom_.GetSize() != springDampVals_.GetSize() ||
         springDom_.GetSize() != springStiffVals_.GetSize() )
      {
        std::string errmsg = "Springs: ";
        errmsg += "#name = " + Info->GenStr(springDom_.GetSize());
        errmsg += ", #dof = " + Info->GenStr(springDof_.GetSize());
        errmsg += ", #mass value = " + Info->GenStr(springMassVals_.GetSize());
        errmsg += ", #damping value = " + Info->GenStr(springDampVals_.GetSize());
        errmsg += ", #stiffness value = " + Info->GenStr(springStiffVals_.GetSize());
        errmsg += ", #dynamics = " + fncname_springs_.GetSize() + '\n';
        Info->Error( errmsg, __FILE__, __LINE__ );
      }

    // We need not have as many function/filenames as loads!
    for ( Integer k = fncname_loads_.GetSize(); k < loadDom_.GetSize(); k++ )
      {
        fncname_loads_.Push_back( "none" );
      }

    // ????????????????????????????????????????????????????????
    // We need not have as many function/filenames as springs!
    for ( Integer k = fncname_springs_.GetSize(); k < springDom_.GetSize(); k++ )
      {
        fncname_springs_.Push_back( "none" );
      }

#ifdef DEBUG
    // loads
    (*debug) << "Assemble::SetGeneralParams: We got " << loadDom_.GetSize()
             << " interfaces with loads" << std::endl;
    (*debug) << "Loads: #interfaces = " << loadDom_.GetSize()
             << ", #dof = " << loadDof_.GetSize()
             << ", #value = " << loadVals_.GetSize()
             << ", #dynamics = " << fncname_loads_.GetSize() << std::endl;
    for ( UInt k = 0; k < loadDom_.GetSize(); k++ )
      {
        (*debug) << "Loads: interface = " << loadDom_[k]
                 << ", dof = " << loadDof_[k]
                 << ", value = " << loadVals_[k];
        if ( k < fncname_loads_.GetSize() )
          {
            (*debug) << ", dynamics = " << fncname_loads_[k] << std::endl;
          }
        else
          {
            (*debug) << ", dynamics = " << std::endl;
          }
      }
    // springs
    (*debug) << "Assemble::SetGeneralParams: We got " << springDom_.GetSize()
             << " interfaces with springs" << std::endl;
    (*debug) << "Springs: #interfaces = " << springDom_.GetSize()
             << ", #dof = " << springDof_.GetSize()
             << ", #mass value = " << springMassVals_.GetSize()
             << ", #damping value = " << springDampVals_.GetSize()
             << ", #stiffness value = " << springStiffVals_.GetSize()
             << ", #dynamics = " << fncname_springs_.GetSize() << std::endl;
    for ( UInt k = 0; k < springDom_.GetSize(); k++ )
      {
        (*debug) << "Springs: interface = " << springDom_[k]
                 << ", dof = " << springDof_[k]
                 << ", mass value = " << springMassVals_[k]
                 << ", damping value = " << springDampVals_[k]
                 << ", stiffness value = " << springStiffVals_[k];
        if ( k < fncname_springs_.GetSize() )
          {
            (*debug) << ", dynamics = " << fncname_springs_[k] << std::endl;
          }
        else
          {
            (*debug) << ", dynamics = " << std::endl;
          }
      }
#endif

    // for every domain, we need an own integrator list ==========
    integrators_.Resize(subdoms_.GetSize());
    rhsSrcIntegrators_.Resize(subdoms_.GetSize());
    fncname_rhs_.Resize(subdoms_.GetSize());
    rhsIntegrators_.Resize(subdoms_.GetSize());
    for (int i=0; i<subdoms_.GetSize();i++) {
      integrators_[i]      = new StdVector<IntegratorDescriptor *>;
      rhsSrcIntegrators_[i]     = new StdVector<BaseIntDescriptor *>;
      rhsIntegrators_[i]        = new StdVector<BaseIntDescriptor *>;
    }

    surfintegrators_.Resize(surfdoms_.GetSize());
    rhsSrcSurfIntegrators_.Resize(surfdoms_.GetSize());
    fncname_rhsSurf_.Resize(surfdoms_.GetSize());
    for (int i=0; i<surfdoms_.GetSize();i++) {
      surfintegrators_[i]       = new StdVector<IntegratorDescriptor *>;
      rhsSrcSurfIntegrators_[i] = new StdVector<BaseIntDescriptor *>;
    }
    
  }
  

  //  void Assemble::SetupMatrixGraph(Integer numeq, Integer graphType)
  void Assemble::SetupMatrixGraph(Integer numeq)
  {
    ENTER_FCN( "Assemble::SetupMatrixGraph" );
    

    //initialize matrix graph
 
    algsys_->CreateGraph(numeq, graphType_);

    // set the graph - connectivity matrix
    BaseFE * ptElem; 
    Integer nsub, iel;
    FEType fe_type;
  
    StdVector<Integer> connecth;
    StdVector<Integer> connect_PDE;
    for (nsub=0; nsub<subdoms_.GetSize(); nsub++)
      {
        StdVector<Elem*> elemssd;
        ptgrid_->GetElemSD(elemssd,subdoms_[nsub],actlevel_);

        for (iel=0; iel < elemssd.GetSize(); iel++)
          {  
            ptElem=elemssd[iel]->ptElem;
            connecth = elemssd[iel]->connect;
            ptEQN_->Node2EQN(connecth, connect_PDE);

            fe_type=elemssd[iel]->ptElem->feType();
            //std::cerr << "SetElelemtPos" << std::endl << "----------" << std::endl;
            //std::cerr << "connect: " << std::endl << connecth << std::endl;
            //std::cerr << "connect_PDE " << std::endl;
            //std::cerr << connect_PDE << std::endl << std::endl;
            algsys_->SetElementPos(connect_PDE.GetPointer(),connect_PDE.GetSize(),fe_type);
          }
      }
  }



  Integer Assemble::GetNrBCDof(const std::string & dofStartString)
  {
    ENTER_FCN( "Analysis::GetNrBCDof" );
    
    Integer nrActDof = 0;
    
    if (dofStartString == "ux")
      nrActDof = 1;
    else
      if (dofStartString == "uy")
        nrActDof = 2;
      else
        if (dofsPerNode_ == 3)
          if (dofStartString == "uz")
            nrActDof = 3;
          else
            Error("Unknown dof-type in homog. BC; substring must start with ux, uy or uz!!",__FILE__,__LINE__);
        else
          Error("Unknown dof-type in homog. BC; substring must start with ux or uy!!",__FILE__,__LINE__);
    
    return nrActDof;
  }

  //set all FE-Elements to reduced integration  
  void Assemble::SetFE2ReducedInt()
  {
    ptQ1    -> SetReducedIntegration();
    ptQ2    -> SetReducedIntegration();
    ptTet1  -> SetReducedIntegration();
    ptL1    -> SetReducedIntegration();
    ptL2    -> SetReducedIntegration();
    ptTr1   -> SetReducedIntegration();
    ptTr2   -> SetReducedIntegration();
    ptHexa1 -> SetReducedIntegration();
    ptHexa2 -> SetReducedIntegration();
    //    ptPyra1 -> SetReducedIntegration();
    //    ptWedge1-> SetReducedIntegration();
  }

  //set all FE-Elements to reduced integration  
  void Assemble::SetFE2StandardInt()
  {
    ptQ1    -> SetStandardIntegration();
    ptQ2    -> SetStandardIntegration();
    ptTet1  -> SetStandardIntegration();
    ptL1    -> SetStandardIntegration();
    ptL2    -> SetStandardIntegration();
    ptTr1   -> SetStandardIntegration();
    ptTr2   -> SetStandardIntegration();
    ptHexa1 -> SetStandardIntegration();
    ptHexa2 -> SetStandardIntegration();
    //    ptPyra1 -> SetStandardIntegration();
    //    ptWedge1-> SetStandardIntegration();
  }


  // define integrators
  void Assemble::AddRhsIntegrator(BaseForm * integrator,
                                  const std::string & subDomName, 
                                  const Integer nonLin)
  {
    ENTER_FCN( "Assemble::AddRhsIntegrator" );
    BaseIntDescriptor * actRhsID = new  BaseIntDescriptor(integrator, nonLin);
    rhsIntegrators_[SubDomIndex(subDomName)]->Push_back(actRhsID);
  }


  //needed for static and transient analysis
  void Assemble::AddRhsSrcIntegrator(BaseForm * integrator,
                                     const std::string & subDomName, 
                                     const std::string fncname,
                                     const Integer nonLin)
  {
    ENTER_FCN( "Assemble::AddRhsSrcIntegrator" );
    BaseIntDescriptor * actRhsID = new  BaseIntDescriptor(integrator, nonLin);
    rhsSrcIntegrators_[SubDomIndex(subDomName)]->Push_back(actRhsID);
    fncname_rhs_[SubDomIndex(subDomName)] = fncname;
  }

  //needed for harmonic analysis
  void Assemble::AddRhsSrcIntegrator(BaseForm * integrator,
                                     const std::string & subDomName, 
                                     const Double phaseval,
                                     const Integer nonLin)
  {
    ENTER_FCN( "Assemble::AddRhsSrcIntegrator" );
    BaseIntDescriptor * actRhsID = new  BaseIntDescriptor(integrator, nonLin);
    rhsSrcIntegrators_[SubDomIndex(subDomName)]->Push_back(actRhsID);
    rhsSrcPhase_[SubDomIndex(subDomName)] = phaseval;
  }


  //needed for static and transient analysis
  void Assemble::AddRhsSrcSurfIntegrator(BaseForm * integrator,
                                         const std::string & subDomName, 
                                         const std::string fncname,
                                         const Integer nonLin)
  {
    ENTER_FCN( "Assemble::AddRhsSrcSurfIntegrator" );
    
    BaseIntDescriptor * actRhsID = new  BaseIntDescriptor(integrator, nonLin);
    rhsSrcSurfIntegrators_[SurfDomIndex(subDomName)]->Push_back(actRhsID);
    fncname_rhsSurf_[SurfDomIndex(subDomName)] = fncname;
  }

  //needed for harmonic analysis
  void Assemble::AddRhsSrcSurfIntegrator(BaseForm * integrator,
                                         const std::string & subDomName, 
                                         const Double phaseval,
                                         const Integer nonLin)
  {
    ENTER_FCN( "Assemble::AddRhsSrcSurfIntegrator" );
    BaseIntDescriptor * actRhsID = new  BaseIntDescriptor(integrator, nonLin);
    rhsSrcSurfIntegrators_[SurfDomIndex(subDomName)]->Push_back(actRhsID);
    rhsSrcSurfPhase_[SurfDomIndex(subDomName)] = phaseval;
  }


  // ==========================================================
  // STATIC ANALYSIS
  // ==========================================================


  StaticAssemble::StaticAssemble(BaseSystem * algsys, Grid * agrid)
    :Assemble(algsys, agrid)
  {
    ENTER_FCN( "StaticAssemble::StaticAssemble" );
    graphType_    = NODEGRAPH; 
    SetAnalysisType(STATIC);
  }
  

  // define integrators
  void StaticAssemble::AddIntegrator( BaseForm *integrator,
                                      const std::string & subDomName,
                                      const FEMatrixType destinationMatrix,
                                      const Integer nonLin ) {

    ENTER_FCN( "StaticAssemble::AddIntegrator" );

    FEMatrixType actMatType = destinationMatrix;
    
    if (actMatType == STIFFNESS)
      actMatType = SYSTEM;

    if (actMatType !=  SYSTEM)
      return;

    IntegratorDescriptor *actID = new IntegratorDescriptor(integrator,
                                                           actMatType, nonLin);
    integrators_[SubDomIndex(subDomName)]->Push_back(actID);
  }


  // define integrators
  void StaticAssemble::AddSurfIntegrator( BaseForm *integrator,
                                          const std::string &subDomName,
                                          const FEMatrixType destinationMatrix,
                                          const Integer nonLin ) {

    ENTER_FCN( "StaticAssemble::AddSurfIntegrator" );
   
    FEMatrixType actMatType = destinationMatrix;
    
    if (actMatType == STIFFNESS)
      actMatType = SYSTEM;

    if (actMatType !=  SYSTEM)
      return;

    IntegratorDescriptor * actID = new IntegratorDescriptor(integrator,
                                                            actMatType,
                                                            nonLin);
    surfintegrators_[SurfDomIndex(subDomName)]->Push_back(actID);
  }

  // define integrators
  void StaticAssemble::AddIntegrator(IntegratorDescriptor * actID,
                                     const std::string & subDomName) {

    ENTER_FCN( "StaticAssemble::AddIntegrator" );

    if (actID->DestMat() == STIFFNESS)
      actID->SetDestMat(SYSTEM);

    if (actID->DestMat() !=  SYSTEM)
      return;

    integrators_[SubDomIndex(subDomName)]->Push_back(actID);
  }
    

  // ==========================================================
  // TRANSIENT ANALYSIS
  // ==========================================================


  TransientAssemble::TransientAssemble(BaseSystem * algsys, Grid * agrid)
    :Assemble(algsys, agrid)
  {
    ENTER_FCN( "TransientAssemble::TransientAssemble" );
    graphType_    = NODEGRAPH; 
    SetAnalysisType(TRANSIENT);

    stiffnessMatrix_  = TRUE;
    massMatrix_       = TRUE;
  }

  // define integrators
  void TransientAssemble::AddIntegrator( BaseForm * integrator,
                                         const std::string & subDomName,
                                         const FEMatrixType destinationMatrix,
                                         const Integer nonLin ) {

    ENTER_FCN( "TransientAssemble::AddIntegrator" );
    
    if (destinationMatrix == SYSTEM)
      Info->Error( "In transient assembling, no SYSTEM matrix may be defined "
                   "directly", __FILE__, __LINE__ );
    
    IntegratorDescriptor * actID = new IntegratorDescriptor(integrator, destinationMatrix, nonLin);
    integrators_[SubDomIndex(subDomName)]->Push_back(actID);
  }


  // define integrators
  void TransientAssemble::AddSurfIntegrator(BaseForm * integrator,
                                            const std::string & subDomName,
                                            const FEMatrixType destinationMatrix, const Integer nonLin)
  {
    ENTER_FCN( "TransientAssemble::AddSurfIntegrator" );
    if (destinationMatrix == SYSTEM)
      Info->Error("In transient assembling, no SYSTEM matrix may be defined directly", __FILE__, __LINE__);

    IntegratorDescriptor * actID = new IntegratorDescriptor(integrator, destinationMatrix, nonLin);
    surfintegrators_[SurfDomIndex(subDomName)]->Push_back(actID);
  }


  // define integrators
  void TransientAssemble::AddIntegrator( IntegratorDescriptor * actID,
                                         const std::string & subDomName ) {
    ENTER_FCN( "TransientAssemble::AddIntegrator" );

    if (actID->DestMat() == SYSTEM)
      Info->Error("In transient assembling, no SYSTEM matrix may be defined directly", __FILE__, __LINE__);

    integrators_[SubDomIndex(subDomName)]->Push_back(actID);
  }


  // ==========================================================
  // RHS INTEGRATOR DESCRIPTOR
  // ==========================================================

  BaseIntDescriptor::BaseIntDescriptor() : 
    integrator(NULL), nonLin(FALSE), reducedIntegration_(FALSE)
  {
    ENTER_FCN( "BaseIntDescriptor::BaseIntDescriptor" );
  }

  BaseIntDescriptor::~BaseIntDescriptor()
  {
    ENTER_FCN( "BaseIntDescriptor::~BaseIntDescriptor" );
    if (integrator != NULL)
      delete integrator;
  }
  

  // define integrators
  BaseIntDescriptor::BaseIntDescriptor(BaseForm * aIntegrator,
                                       const Boolean aNonLin)
    : integrator(aIntegrator), nonLin(aNonLin), reducedIntegration_(FALSE)
  {
    ENTER_FCN( "BaseIntDescriptor::BaseIntDescriptor " );
  }
  

  // ==========================================================
  // INTEGRATOR DESCRIPTOR
  // ==========================================================

  IntegratorDescriptor::IntegratorDescriptor()
    :BaseIntDescriptor(),
     destinationMatrix(SYSTEM),
     secondaryMatrix(NOTYPE),
     secMatFac(0.0),
     piezoMaterialType_(realMaterialParameter)
  
  {
    ENTER_FCN( "IntegratorDescriptor::IntegratorDescriptor" );
    //    piezoMaterialType_=realMaterialParameter;
  }


  // define integrators
  IntegratorDescriptor::IntegratorDescriptor(BaseForm * aIntegrator,
                                             const enum FEMatrixType aDestMat,
                                             const Boolean aNonLin)
    :BaseIntDescriptor(aIntegrator, aNonLin),
     destinationMatrix(aDestMat),
     secondaryMatrix(NOTYPE),
     secMatFac(0.0),
     piezoMaterialType_(realMaterialParameter) {
    ENTER_FCN( "IntegratorDescriptor::IntegratorDescriptor" );
    //   piezoMaterialType_=realMaterialParameter;
  }

  // defines a secondary destination for the calculated element marices of an
  // integrator
  void IntegratorDescriptor::SetSecondaryMat( FEMatrixType aSecMat,
                                              Double aSecMatFac,
                                              AnalysisType analysisType ) {

    ENTER_FCN( "IntegratorDescriptor::SetSecondaryMat" );
    FEMatrixType MatType = aSecMat;

    if (analysisType == HARMONIC) {
      if (aSecMat == STIFFNESS || aSecMat == MASS || aSecMat == DAMPING )
        MatType = SYSTEM;
      else
        {
          std::string error_msg = "Matrix type ";
          error_msg += aSecMat + " not supported in harmonic analysis";
          Error(error_msg.c_str());
        }

      SetOrigSecMatrixType(aSecMat);
    }

    secondaryMatrix = MatType;
    secMatFac = aSecMatFac;
  }


  IntegratorDescriptor::~IntegratorDescriptor() {
    ENTER_FCN( "IntegratorDescriptor::~IntegratorDescriptor" );
  }


  // ==========================================================
  // HARMONIC ANALYSIS
  // ==========================================================


  HarmonicAssemble::HarmonicAssemble(BaseSystem * algsys, Grid * agrid)
    :Assemble(algsys, agrid)
  {
    ENTER_FCN( "HarmonicAssemble::HarmonicAssemble" );
    graphType_    = NODEGRAPH; 
    SetAnalysisType(HARMONIC);
  }

  /// set actual frequency (already multiplied by 2*pi)
  void HarmonicAssemble::SetFrequency(Double frequency)
  {
    ENTER_FCN( "HarmonicAssemble::SetFrequency" );

    actFreq_ = 2*PI*frequency;

  } 


  // define integrators
  void HarmonicAssemble::AddIntegrator(BaseForm * integrator,
                                       const std::string & subDomName,
                                       const FEMatrixType destinationMatrix,
                                       const Integer nonLin) {

    ENTER_FCN( "HarmonicAssemble::AddIntegrator" );

    FEMatrixType actMatType = destinationMatrix;
    FEMatrixType matType;  
    
    if (actMatType == STIFFNESS || actMatType == MASS ||
        actMatType == DAMPING )
      matType = SYSTEM;
    else {
      std::string error_msg = "Matrix type ";
      error_msg += actMatType + " not supported in harmonic analysis";
      Error(error_msg.c_str());
    }

    IntegratorDescriptor * actID = new IntegratorDescriptor(integrator,
                                                            matType, nonLin);
    actID->SetOrigMatrixType(actMatType);

    integrators_[SubDomIndex(subDomName)]->Push_back(actID);

  }

  // define integrators
  void HarmonicAssemble::AddSurfIntegrator(BaseForm * integrator,
                                           const std::string & subDomName,
                                           const FEMatrixType destinationMatrix,
                                           const Integer nonLin)
  {
    ENTER_FCN( "HarmonicAssemble::AddSurfIntegrator" );
   
    FEMatrixType actMatType = destinationMatrix;
    FEMatrixType matType;
    
    if (actMatType == STIFFNESS || actMatType == MASS || actMatType == DAMPING )
      matType = SYSTEM;
    else
      {
        std::string error_msg = "Matrix type ";
        error_msg += actMatType + " not supported in harmonic analysis";
        Error(error_msg.c_str());
      }

    IntegratorDescriptor * actID = new IntegratorDescriptor(integrator, matType, nonLin);
    actID->SetOrigMatrixType(actMatType);

    surfintegrators_[SurfDomIndex(subDomName)]->Push_back(actID);
  }

  // define integrators
  void HarmonicAssemble::AddIntegrator(IntegratorDescriptor * actID,
                                       const std::string & subDomName)
  {
    ENTER_FCN( "HarmonicAssemble::AddIntegrator" );

    actID->SetOrigMatrixType(actID->DestMat());
    if (actID->DestMat() == STIFFNESS || actID->DestMat() == MASS || actID->DestMat() == DAMPING )
      actID->SetDestMat(SYSTEM);
    else
      {
        std::string error_msg = "Matrix type ";
        error_msg += actID->DestMat() + " not supported in harmonic analysis";
        Error(error_msg.c_str());
      }

    integrators_[SubDomIndex(subDomName)]->Push_back(actID);
  }


  void HarmonicAssemble::
  TransformMatrix2Harmonic( Vector<Double>& harmMat,
                            Matrix<Double> origMat,
                            const FEMatrixType matrixType,
                            const piezoMaterialType piezoMatType ) {

    ENTER_FCN( "HarmonicAssemble::TransformMatrix2Harmonic" );

    Integer numRow = origMat.GetSizeRow();
    Integer numCol = origMat.GetSizeCol();
    harmMat.Resize(2*numRow*numCol);

    Integer k=0;

    if (piezoMatType == realMaterialParameter) {

      if (matrixType == STIFFNESS) {

	// std::cout<<"real_stiff - actfreq: "<< actFreq_<<std::endl;
	for (Integer row=0; row<numRow; row++)
	  for (Integer col=0; col<numCol; col++) {
	    harmMat[k] = origMat[row][col];
	    k++;
	  }
      }
    
      else if (matrixType == MASS)
        {
          //std::cout<<"real_mass"<<std::endl;

          Double factor = -actFreq_*actFreq_;
          for (Integer row=0; row<numRow; row++)
            for (Integer col=0; col<numCol; col++) {
              harmMat[k] = factor*origMat[row][col];
              k++;
            }
        }
    
      else if (matrixType == DAMPING)
        {       
          //    std::cout<<"real_damping"<<std::endl;
          Double factor = actFreq_;
        
          k=numRow*numCol;
          for (Integer row=0; row<numRow; row++)
            for (Integer col=0; col<numCol; col++) {
              harmMat[k] = factor*origMat[row][col];
              k++;
            }
        }
    } // end, if piezoMatType == real...

    else if(piezoMatType == imagMaterialParameter){  // the "imaginary parts"
   
      if (matrixType == STIFFNESS)
        {
          //std::cout<<"comlex_stiff"<<std::endl;
          k=numRow*numCol;
          for (Integer row=0; row<numRow; row++)
            for (Integer col=0; col<numCol; col++) {
              harmMat[k] = origMat[row][col];
              k++;
            }
        }
       
      else if (matrixType == DAMPING)
        {
          Double factor = actFreq_;
          //    std::cout<<"comlex_damping"<<std::endl;
        
          k=0;  
          for (Integer row=0; row<numRow; row++)
            for (Integer col=0; col<numCol; col++) {
              harmMat[k] = -factor*origMat[row][col];
              k++;
            }
        }

    } // end if piezoMatType == imag

    else {
      std::cerr<<"\n piezoMaterialType" << piezoMatType << "not specified "<<std::endl;
    }
  }



  void  HarmonicAssemble::TransformVector2Harmonic(Vector<Double>& harmVec, Vector<Double> origVec, const Double valPhase)
  {
    ENTER_FCN( "HarmonicAssemble::TransformVector2Harmonic" );

    Integer size = origVec.GetSize();
    harmVec.Resize(2*size);

    Double valReal = cos(valPhase);
    Double valImag = sin(valPhase);

    Integer k=0;
    std::cout<<" HarmonicAssemble::TransformVector2Harmonic valPhase = " <<valPhase <<std::endl;
    //real part
    for (Integer i=0; i<size; i++) {
      harmVec[k] = origVec[i]*valReal;
      k++;
    }

    //imaginary part
    for (Integer i=0; i<size; i++) {
      harmVec[k] = origVec[i]*valImag;
      k++;
    }

  }


} // end namespace CoupledField
