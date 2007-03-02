#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>

#include "magEdgePDE.hh"
#include <Forms/forms_header.hh>
#include <Estimator/spaceerror.hh>
#include <General/environment.hh>
#include <DataInOut/WriteInfo.hh>
 
namespace CoupledField
{

  MagEdgePDE::MagEdgePDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc, FileType *aptFileType, 
                         WriteResults *aptOut)
    :BasePDE(aptgrid, aptbcs, aptFileType, aptOut, aptTimeFunc)
  {
#ifdef TRACE
    (*trace) << "entering MagEdgePDE::MagEdgePDE " << std::endl;
#endif


    Error( "Not working at the moment", __FILE__, __LINE__ );

    dofsperedge_ =1;

    SetMatrixFactors();

    pdename_    = "magnetic";
    pdematerialclass_ = "magnetic";
  
    conf->getsubdompde(subdoms_,pdename_);
    ReadBCs(pdename_);

    // Map global numeration of element and nodes to local one
    AssignPDENodeNumbers(mesh2PDENode_, pde2MeshNode_, subdoms_);  
    AssignPDEElemNumbers(mesh2PDEElem_, pde2MeshElem_, subdoms_);
    numPDENodes_ = pde2MeshNode_.size();
    numElems_ = pde2MeshElem_.size();

    bFieldRe_ = NULL;
    magVecPotRe_ = NULL;
  
    ReadCoils();

    // Currently freq_ is unused
    // if (analysistype_==HARMONIC) {
    //  conf->get("frequency", freq_, pdename_);
    // }
  }


  void MagEdgePDE::SetMatrixFactors()
  {
#ifdef TRACE
    (*trace) << "entering MagEdgePDE::SetMatrixFactors" << std::endl;
#endif
 
    matrix_factor_[0] = 0.0;  // factor for stiffness matrix
    matrix_factor_[1] = 0.0;  // factor for damping matrix
    matrix_factor_[2] = 0.0;  // factor for convection matrix
    matrix_factor_[3] = 0.0;  // factor for mass matrix
  }


  void MagEdgePDE::DiscreteParamsPDE()
  {
#ifdef TRACE
    (*trace) << "entering MagEdgePDE::DiscreteParamsPDE" << std::endl;
#endif
    if (analysistype_ == STATIC)
      MatrixType_   = RSCALAREDGE; 
    else if (analysistype_ == HARMONIC)
      MatrixType_   = CSCALAREDGE;
    else 
      Error("Analysis Type not supported",__FILE__,__LINE__);
  
    SystemMatrix_ = true;
    GraphType_    = EDGEGRAPH;
  }

  void MagEdgePDE::SetAlgSys(const Integer as_sysid)
  {
#ifdef TRACE
    (*trace) << "entering  BasePDE::SetAlgSys" << std::endl;
#endif

    as_sysid_ = as_sysid;

    //allocate according algebraic system
    algsys_ = new StandardSystem();

    //set solver parameters  
    SetSolverParameters();
    algsys_->SetCycle(VCYCLE);
    algsys_->SetSmoothType(1);
    algsys_->SetSmoothStepFor(2);
    algsys_->SetSmoothStepBack(2);

    //ask the PDE discrete form
    DiscreteParamsPDE();

    //set the node graph
    SetupMatrixGraph(NumPDENodes_, GraphType_);

    //set the edge graph
    SetupEdgeGraph();

    //allocate the necessary matrices as well as solver and preconditioner
    CreateMatrices_Solver();
  }



  void MagEdgePDE::CreateMatrices_Solver()
  {
#ifdef TRACE
    (*trace) << "entering  MagEdgePDE::CreateMatrices_Solver" << std::endl;
#endif

    Integer matrixclass;
    Integer matrixsystype[5];    
    Integer graphtype; 
    Integer numconstraints;

    //eval number of dirichlet-edges
    EvalNumDirichlet();

    numconstraints = 0;  // currently not handled
    matrixclass    = MatrixType_;
    graphtype      = GraphType_;

    if (SystemMatrix_    != 1)
      Error("One needs at least a system matrix!",__FILE__,__LINE__);

    if (SystemMatrix_     == 1) matrixsystype[0] = SYSTEM;      // memory for the system matrix
    if (StiffnessMatrix_  == 1) matrixsystype[1] = STIFFNESS;   // memory for the stiffness matrix
    if (DampingMatrix_    == 1) matrixsystype[2] = DAMPING;     // memory for the damping matrix
    if (ConvectionMatrix_ == 1) matrixsystype[3] = CONVECTION;  // memory for the convection matrix
    if (MassMatrix_       == 1) matrixsystype[4] = MASS;        // memory for the mass matrix


    //account for real and imaginray part
    if (analysistype_ == HARMONIC)
      algsys_->CreateMatrix(matrixsystype, matrixclass, graphtype, dofsperedge_,  numEdgedir_*2,
                            numconstraints);
    else
      algsys_->CreateMatrix(matrixsystype, matrixclass, graphtype, dofsperedge_,  numEdgedir_,
                            numconstraints);


    //create solver and preconditioner
    algsys_->CreateSolver();
    algsys_->CreatePrecond(matrixclass);

    //now reset AlgebraicSystem 
    algsys_->InitRHS();
    algsys_->InitSol();

    for (Integer i=0;i<5;i++)
      {
        if (matrixsystype[i] !=0)
          algsys_->InitMatrix(i+1);
      }

  }


  void MagEdgePDE::SolveStepStatic(const Integer level, const Double aTime)
  {
#ifdef TRACE
    (*trace) << "entering MagEdgePDE::SolveStepStatic" << std::endl;
#endif

    Integer update = 0;
    Integer job = 4;

    Double * ptsol;

    //compute and assemble element matrices
    SetupMatrices(level);

    // calculate source term
    SetupRHS(level);
  
    //account for bcs
    SetBCs(level,update,0);

    CPUClock cpuClock;
    Double startTime = cpuClock.GetTime ();
    
    algsys_->CalcPrecond();
    Double preCondTime = cpuClock.GetTime ();
    (*cla) << std::endl << "TIME for PRECONDITIONER SETUP: " << preCondTime - startTime << std::endl;
    //    std::cout << "TIME for PRECONDITIONER SETUP: " << preCondTime - startTime << std::endl;
    
    algsys_->Solve();
    Double solveTime = cpuClock.GetTime ();
    (*cla) << std::endl << "TIME for SOLUTION: " << solveTime - preCondTime << std::endl;
    //    std::cout << "TIME for SOLUTION: " << solveTime - preCondTime << std::endl;

    algsys_->CalcComplexity();

    ptsol = algsys_->GetSolutionVal();

    // save solution
    for(Integer i=0; i < size_; i++)
      solRe_[i] = ptsol[i];


#ifdef DEBUG
    (*debug) << "SolveStepStatic: \n solution \n";
    for (Integer i=0; i < size_; i++)
      (*debug) << solRe_[i] << " ";
  
    (*debug) << std::endl;  
#endif
  }


  void MagEdgePDE::SolveStepHarmonic(const Integer level)
  {
#ifdef TRACE
    (*trace) << "entering MagEdgePDE::SolveStepHarmonic" << std::endl;
#endif

    Integer update = 0;
    Integer job = 4;

    Double * ptsol;

    //compute and assemble element matrices
    SetupMatrices(level);

    // calculate source term
    SetupRHS(level);
  
    //account for bcs
    SetBCs(level,update,0);


    CPUClock cpuClock;
    Double startTime = cpuClock.GetTime ();
    
    algsys_->CalcPrecond();
    Double preCondTime = cpuClock.GetTime ();
    (*cla) << std::endl << "TIME for PRECONDITIONER SETUP: " << preCondTime - startTime << std::endl;
    //    std::cout << "TIME for PRECONDITIONER SETUP: " << preCondTime - startTime << std::endl;
    
    algsys_->Solve();
    Double solveTime = cpuClock.GetTime ();
    (*cla) << std::endl << "TIME for SOLUTION: " << solveTime - preCondTime << std::endl;
    //    std::cout << "TIME for SOLUTION: " << solveTime - preCondTime << std::endl;

    algsys_->CalcComplexity();

    ptsol = algsys_->GetSolutionVal();

    // save solution
    for(Integer i=0; i < size_; i++)
      {
        solRe_[i] = ptsol[2*i];
        solIm_[i] = ptsol[2*i+1];
      }
  

#ifdef DEBUG
    (*debug) << "SolveStepHarmonic: \n solution real part \n";
    for (Integer i=0; i < size_; i++)
      (*debug) << solRe_[i] << " ";
    (*debug) << std::endl;  
  

    (*debug) << "SolveStepHarmonic: \n solution imaginary part \n";
    for (Integer i=0; i < size_; i++)
      (*debug) << solIm_[i] << " ";
    (*debug) << std::endl;  
#endif
  }




  void MagEdgePDE::SetupMatrices(const Integer level)
  {
#ifdef TRACE
    (*trace) << "entering MagEdgePDE::SetupMatrices" << std::endl;
#endif
  
    Matrix<Double> elemmat;  
    Matrix<Double> elemmat_aux;  

    Matrix<Double> ptCoord;
    BaseFE         * ptElem;

    Vector<Integer> connecth, connect_PDE;  
    Integer i, j;
    Double reluctivity, conductivity, harmfactor;

    //curently hard coded for tets
    Integer elemsize_edge = 6;
    std::vector<Integer> epos(elemsize_edge);
    std::vector<Integer> esign(elemsize_edge);

    //regularization parameter
    Double relaxFac = 1.0e-5;
    conf->ifget("relaxFac", relaxFac, pdename_);

    std::vector<Double> harmVec(2*elemsize_edge*elemsize_edge);

    for (i=0; i<subdoms_.size(); i++)
      {
        Double pemeability; 
        materialData_[i].GetPermeability(2,2,permeability);
        reluctivity  = 1.0/permeability;
        //      reluctivity  = 1.0/materialData_[i].GetPermiability();
        Double conductivity; 
        materialData_[i].GetConductivity(2,2,conductivity);

        //if STATIC, set conductivity to zero
        if (analysistype_==STATIC) conductivity = 0.0;

        // small conductivity is needed for regularization
        harmfactor = 1.0;
        if (conductivity <= 0 || analysistype_==STATIC)
          {
            conductivity = reluctivity * relaxFac;
            //if conductivity = 0.0, then we have in the harmonic case no mass matrix;
            //computed mass matrix is just used for for regulraization of stiffness matrix
            harmfactor = 0.0;
          }
        

        std::vector<Elem*> elemssd;
   
        ptgrid_->GetElemSD(elemssd,subdoms_[i],level);

        for (j=0; j < elemssd.size(); j++)
          {  
            ptElem=elemssd[j]->ptElem;

            BaseForm * bilinear_stiff = new CurlCurlEdgeInt(ptElem, reluctivity);
            BaseForm * bilinear_mass  = new MassEdgeInt(ptElem, conductivity);
            BaseForm * lapl           = new LaplaceInt(ptElem, reluctivity, isaxi_);

            connecth=elemssd[j]->connect;
          
            GetElemCoords(connecth, ptCoord, level);

            // CHANGE connecth
            Mesh2PDENode(connect_PDE, connecth, mesh2PDENode_);

            //get the edge numbers and their signs
            GetEdgeNumber(connect_PDE.get(), epos, esign, ptElem);

            // stiffness part
            bilinear_stiff->CalcElementMatrix(ptCoord, elemmat);
          
            CorrectEdgeDir(elemmat, esign);
          

            if(analysistype_==TRANSIENT || analysistype_==STATIC)
              algsys_->SetElementMatrix(elemmat.getinarray(), &epos[0], epos.size(),SYSTEM);

            // store stiffness matrix to special vector needed by algebraic system
            // for harmonic analysis
            Integer k=0;
            
            if(analysistype_==HARMONIC)
              {
                for(Integer iii=0; iii<elemmat.size_row(); iii++)
                  for(Integer jjj=0; jjj < elemmat.size_row(); jjj++)
                    {    
                      harmVec[k] = elemmat[iii][jjj];
                      k++;
                    }
              }
          

            // mass part
            bilinear_mass->CalcElementMatrix(ptCoord, elemmat);
          
            CorrectEdgeDir(elemmat, esign);


            if(analysistype_==TRANSIENT || analysistype_==STATIC)
              algsys_->SetElementMatrix(elemmat.getinarray(), &epos[0], epos.size(),SYSTEM);



            //auxilary matrix
            lapl->CalcElementMatrix(ptCoord, elemmat_aux);


            algsys_->SetAuxElementMatrix(elemmat_aux.getinarray(), 
                                         connect_PDE.get(), connecth.size());
          
            // store mass matrix to special vector needed by algebraic system
            // for harmonic analysis
            
            
            if(analysistype_==HARMONIC)
              {
                Integer kkk=0;
                Double regularizationFactor = reluctivity / conductivity *  relaxFac;
                
                if (k!=elemmat.size_row()*elemmat.size_col())
                  Error("k is wrong!!!!!!!!!!!!!!!!!!!!!",__FILE__,__LINE__);

                for(Integer iii=0; iii<elemmat.size_row(); iii++)
                  for(Integer jjj=0; jjj < elemmat.size_row(); jjj++)
                    {    
                      harmVec[k] = elemmat[iii][jjj] * 2 * PI * freq_ * harmfactor;
                      // regularization of element stiffness matrix
                      harmVec[kkk] += elemmat[iii][jjj] * regularizationFactor;

                      // k initially  set in for loop above!!
                      k++;
                      kkk++;
                    }
                algsys_->SetElementMatrix(&harmVec[0], &epos[0], epos.size(),SYSTEM);
              }

            delete bilinear_stiff;
            delete bilinear_mass;
            delete lapl;
          } 
      }
  }





  void MagEdgePDE::WriteResultsInFile()
  {
#ifdef TRACE
    (*trace) << "entering MagEdgePDE::WriteResultsInFile" << std::endl;
#endif

    Integer step=0;
    Double time=0;

    //   Vector<Double> Potentialh, *Eh, *Fh;
    //   TransformNodeSolution(Potentialh,solRe_);

    if (!bFieldRe_)
      Error("B-Field needs to be calculated before output!",__FILE__,__LINE__);

    if (!magVecPotRe_)
      Error("Magn. vector potential needs to be calculated before output!",__FILE__,__LINE__);


    if (OutFile_->IsGMV())
      {
        Error("GMV Currently not supported",__FILE__,__LINE__);
        //      OutFile_->WriteSolution(solRe_, step, time, "magnetic vector potential");
      
        // Write Out Vector Data
        //      for (ShortInt i=0; i<ptgrid_->GetDim(); i++) 
        //        {
        //          std::ostringstream b_fieldname;
        //          b_fieldname << "Bfield " << i;
        //          OutFile_->WriteDataOnCell(bFieldRe_[i], step, time, b_fieldname.str());
        //        }
      }
    else
      {
        const Integer dim = 3;

        // write magnetic flux      
        std::string fieldname = "mag. flux density";
        Array<Double> outMat(dim, NumElems_);
      
        if (analysistype_==HARMONIC)
          {
            for (Integer i=0; i<ptgrid_->GetDim(); i++) 
              for (Integer j=0; j < NumElems_; j++)
                // calc abs value of complex magnetic field
                // outMat[i][j] = sqrt(bFieldRe_[i][j]*bFieldRe_[i][j] + bFieldIm[i][j]*bFieldIm[i][j]);          
                // write just real part of solution 
                outMat[i][j] = bFieldRe_[i][j];

            OutFile_->WriteElemSolution(outMat, step, time, fieldname);
          
          
            // write magnetic vector potential 
            // fieldname = "eddy current";
            for (Integer i=0; i<ptgrid_->GetDim(); i++) 
              for (Integer j=0; j < NumElems_; j++)
                // write just real part of "eddy current" 
                outMat[i][j] = - 2 * PI * freq_ * magVecPotIm_[i][j];     
          
            //    OutFile_->WriteDataOnCell(outMat, step+1, time, fieldname);
          }
        else
          {
            for (Integer i=0; i<ptgrid_->GetDim(); i++) 
              for (Integer j=0; j < NumElems_; j++)
                outMat[i][j] = bFieldRe_[i][j];   

            OutFile_->WriteElemSolution(outMat, step, time, fieldname);

      
            // write magnetic vector potential 
            // fieldname = "eddy current";
            for (Integer i=0; i<ptgrid_->GetDim(); i++) 
              for (Integer j=0; j < NumElems_; j++)
                outMat[i][j] = magVecPotRe_[i][j];        

            //    OutFile_->WriteElemSolution(outMat, step+1, time, fieldname);
          }
  
      }
  }


  void MagEdgePDE::SetupEdgeGraph()
  {
#ifdef TRACE
    (*trace) << "entering MagEdgePDE::SetupEdgeGraph" << std::endl;
#endif

    algsys_->CreateEdgeGraph();

    //number of unknowns (edges)
    size_ = algsys_->GetEdgeSize();


    //now we can allocate the memory for the solution vector
    solRe_.Resize(size_);
    solRe_.Init(0);

    if (analysistype_ == HARMONIC)
      {
        solIm_.Resize(size_);
        solIm_.Init(0);
      }
  
    //curently hard coded for tets
    Integer elemsize_edge = 6;

    std::vector<Integer> epos(elemsize_edge);

    BaseFE * ptElem; 
    Integer nsub, iel, fe_type;
    Vector<Integer> connecth;
    Integer * pos;

    for (nsub=0; nsub<subdoms_.size(); nsub++)
      {
        std::vector<Elem*> elemssd;
        ptgrid_->GetElemSD(elemssd,subdoms_[nsub],actlevel_);

        for (iel=0; iel < elemssd.size(); iel++)
          {  
            ptElem=elemssd[iel]->ptElem;
          
            //Map Mesh Node numbers to PDE node numbers
            Mesh2PDENode(connecth,elemssd[iel]->connect,mesh2PDENode_);
            fe_type = ptElem->feType();
            if (fe_type != TET) 
              Error("Currently just TETs supported for MagEdgePDE",__FILE__,__LINE__);
          
            pos = connecth.get();

            ptElem->GetGlobalEdgeIndices(epos, pos, algsys_);
          

            //take the absolute value
            for (Integer j=0; j<elemsize_edge; j++)
              epos[j] = abs(epos[j]);

            algsys_->SetElementPosEdge(&epos[0],elemsize_edge,fe_type);
          }
      }

  }


  void MagEdgePDE::EvalNumDirichlet()
  {
#ifdef TRACE
    (*trace) << "entering MagEdgePDE::EvalNumDirichlet" << std::endl;
#endif

    Integer i,j,k;

    //curently hard coded for tets
    Integer elemsize_edge = 6;
    Integer * epos  = new Integer[elemsize_edge];
    Integer * dnode = new Integer[size_];

    //set array to zero
    for (i=0; i<size_; i++)
      dnode[i] = 0;
    
    BaseFE  * ptElem;
    std::vector<Elem*>  SurfD;
    std::vector<std::string> surfDirichlet;

    if (conf->ifgetliststr("SurfeDirichlet", surfDirichlet, "magnetic"))
    
      if (surfDirichlet.size())
        {
          SurfD = ptBCs_->getFacesBC(surfDirichlet[0],actlevel_);
          //    else
          //    Error("No Surfaces specified as Dirichlet Boundaries",__FILE__,__LINE__);
          //  else
          //    Error("No Surfaces specified as Dirichlet Boundaries",__FILE__,__LINE__);

          //hard coded: surface elements of tets are triangles
          Integer surfelemdim = 3;
          Vector<Integer> connecth;
          Integer * pos;
        
          for (Integer iel=0; iel< SurfD.size(); iel++)
            {
              ptElem=SurfD[iel]->ptElem;
            
              //Map Mesh Node numbers to PDE node numbers
              Mesh2PDENode(connecth,SurfD[iel]->connect,mesh2PDENode_);
              pos = connecth.get();
            
              epos[0] = algsys_->GetNode2Edge(pos[1], pos[0]);
              epos[1] = algsys_->GetNode2Edge(pos[2], pos[0]);
              epos[2] = algsys_->GetNode2Edge(pos[2], pos[1]);

              for (j=0; j<surfelemdim; j++)
                dnode[abs(epos[j])-1] = abs(epos[j]);
            
            }
        
          numEdgedir_=0;
          for (i=0; i<size_; i++)
            {
              if (dnode[i] != 0)
                numEdgedir_++;
            }
        
          EdgeDir_ = new Integer[numEdgedir_];
        
          k=0;
          for (i=0; i<size_; i++)
            {
              if (dnode[i])
                {
                  EdgeDir_[k] = dnode[i];
                  k++;
                }
            }
        }
      else
        numEdgedir_=0;


    delete [] dnode;
    delete [] epos;
  }




  void   MagEdgePDE::SetBCs(const Integer level, const Integer update, const Double time)
  {
#ifdef TRACE
    (*trace) << "entering   MagEdgePDE::SetBCs" << std::endl;
#endif

    //currently just homogeneous BC supported

    Double  val;
    Integer edge;

    for (Integer i=0; i<numEdgedir_; i++)
      {  
        edge = EdgeDir_[i];
        val=0; 
        if (analysistype_==STATIC)
          algsys_->SetDirichlet(i+1, edge, val, dofsperedge_, SYSTEM);
        else if (analysistype_==HARMONIC)
          {
            // set real part 
            algsys_->SetDirichlet(i*2+1, edge, val, dofsperedge_,SYSTEM);
            // set imag part 
            algsys_->SetDirichlet(i*2+2, edge, val, dofsperedge_+1, SYSTEM);        
          }
        else if (analysistype_==TRANSIENT)
          Error("BCs for transient not yet implemented!!",__FILE__,__LINE__);
      }

#ifdef TRACE
    (*trace) << " leaving  MagEdgePDE::SetBCs " << std::endl;
#endif 
  }

  void MagEdgePDE::GetEdgeNumber(Integer *pos, std::vector<Integer>& epos, 
                                 std::vector<Integer>& esign, BaseFE * ptElem)
  {
#ifdef TRACE
    (*trace) << "entering MagEdgePDE::GetEdgeNumber" << std::endl;
#endif

    ptElem->GetGlobalEdgeIndices(epos, pos, algsys_);

    Integer elemsize_edge = 6;
    for (Integer j=0; j<elemsize_edge; j++)
      {
        esign[j] = epos[j]/abs(epos[j]);
        epos[j]  = abs(epos[j]);
      }
  }




  void MagEdgePDE::SetupRHS(const Integer level)
  {
#ifdef TRACE
    (*trace) << "entering MagEdgePDE::SetupRHS" << std::endl;
#endif

    //   std::vector<std::string> coil;
    //   if (conf->ifgetliststr("Coil", coil, "magnetic"))
    
  
    std::vector<Double> elemVec;  
    std::vector<Double> elemVecHarm;  
    Matrix<Double> ptCoord;
    BaseFE         * ptElem;

    Vector<Integer> connecth, connect_PDE;  
    Integer i, j;
    Double reluctivity, conductivity;

    //curently hard coded for tets
    Integer elemsize_edge = 6;
    elemVecHarm.resize(elemsize_edge*2);
  
    std::vector<Integer> epos(elemsize_edge);
    std::vector<Integer> esign(elemsize_edge);

    for (Integer actCoil=0; actCoil < coilDomain_.size(); actCoil++)
      for (i=0; i<subdoms_.size(); i++)
        if (coilDomain_[actCoil] == subdoms_[i])
          {       
            std::vector<Elem*> elemssd;
   
            ptgrid_->GetElemSD(elemssd,subdoms_[i],level);

            const Double  currentDens   = coilDef_[actCoil].current / coilDef_[actCoil].coilArea;
            const Integer currDirection = coilDef_[actCoil].iDir;
            const Double  currPhase     = coilDef_[actCoil].currentPhase;

            std::vector<Double> * coilMidPt = &coilDef_[actCoil].coilMidPt;


            for (j=0; j < elemssd.size(); j++)
              {  
                ptElem=elemssd[j]->ptElem;
              
                BaseForm * rhsSource = new LinearEdgeInt(ptElem, currentDens, 
                                                         currDirection, coilMidPt);
              
                connecth=elemssd[j]->connect;
              
                GetElemCoords(connecth, ptCoord, level);
              
                // CHANGE connecth
                Mesh2PDENode(connect_PDE,connecth,mesh2PDENode_);
              
                //get the edge numbers and their signs
                GetEdgeNumber(connect_PDE.get(), epos, esign, ptElem);
              
              
                // stiffness part
                rhsSource->CalcElemVector(ptCoord, elemVec);
              
                // correct sign of entries in elemVec due to orientation of edge
                for(Integer ii=0; ii<elemVec.size(); ii++)
                  elemVec[ii] *= esign[ii];
              
                if (analysistype_==STATIC)
                  algsys_->SetElementRHS(&elemVec[0], &epos[0], epos.size());

                else if (analysistype_==HARMONIC)
                  {
                    for(Integer ii=0; ii<elemVec.size(); ii++)
                      {
                        elemVecHarm[ii] = elemVec[ii] * cos(currPhase);
                        elemVecHarm[ii+elemsize_edge] = elemVec[ii] * sin(currPhase);
                      }
                    algsys_->SetElementRHS(&elemVecHarm[0], &epos[0], epos.size());
                  }
                delete rhsSource;
              }
          }
  }


  // Calculation of the magnetic field
  void MagEdgePDE::PostProcess(const Integer level)
  {
#ifdef TRACE
    (*trace) << "entering MagEdgePDE::PostProcess" << std::endl;
#endif

    ShortInt dim = ptgrid_->GetDim();
    CurlEdgeOp * magFieldRe = new CurlEdgeOp(ptgrid_, this, &mesh2PDENode_, &solRe_, level, algsys_); 
    CurlEdgeOp * magFieldIm;

    std::vector<Double> lCoord(dim, 1./4);
    std::vector<Elem*> elemssd;
    Vector<Double> elemMagField;
    Vector<Double> elemMagVecPot;  
  

    bFieldRe_ = new Vector<Double>[dim];
    // Resize solution arrays
    for( Integer i=0; i<dim; i++)
      bFieldRe_[i].Resize(NumElems_);  


    magVecPotRe_ = new Vector<Double>[dim];
    // Resize solution arrays
    for( Integer i=0; i<dim; i++)
      magVecPotRe_[i].Resize(NumElems_);  
  

    if (analysistype_==HARMONIC)
      {  
        magFieldIm = new CurlEdgeOp(ptgrid_, this, &mesh2PDENode_, &solIm_, level, algsys_); 

        bFieldIm_ = new Vector<Double>[dim];
        // Resize solution arrays
        for( Integer i=0; i<dim; i++)
          bFieldIm_[i].Resize(NumElems_);  

        magVecPotIm_ = new Vector<Double>[dim];
        // Resize solution arrays
        for( Integer i=0; i<dim; i++)
          magVecPotIm_[i].Resize(NumElems_);
      }


    Integer actEl = 0;
  
    // loop over all subdomains
    for (Integer isd=0; isd<subdoms_.size(); isd++)
      {
        // get vector of Elem of subdomain with color: subdoms[isd]
        ptgrid_->GetElemSD(elemssd, subdoms_[isd], level);
      
        // loop over elements of subdomain
        for (Integer iel=0; iel< elemssd.size(); iel++)
          {
            magFieldRe->CalcElemCurlEdge(elemMagField, elemssd[iel], lCoord);
            magFieldRe->CalcElemMagVec(elemMagVecPot, elemssd[iel], lCoord);

            for (Integer k=0; k<dim; k++)
              {
                bFieldRe_[k][actEl]    = elemMagField[k];
                magVecPotRe_[k][actEl] = elemMagVecPot[k];
              }

            if (analysistype_==HARMONIC)
              {
                magFieldIm->CalcElemCurlEdge(elemMagField, elemssd[iel], lCoord);
                magFieldIm->CalcElemMagVec(elemMagVecPot, elemssd[iel], lCoord);
              
                for (Integer k=0; k<dim; k++)
                  {
                    bFieldIm_[k][actEl]    = elemMagField[k];
                    magVecPotIm_[k][actEl] = elemMagVecPot[k];
                  }
              }
  
            actEl++;
          }
      }

    delete magFieldRe;
    if (analysistype_==HARMONIC)
      delete magFieldIm;
  }



  // reads all information in the config file concerning coils 
  void MagEdgePDE::ReadCoils()
  {
#ifdef TRACE
    (*trace) << "entering MagEdgePDE::ReadCoils" << std::endl;
#endif  

    conf->ifgetliststr("list_coils", coilDomain_, "magnetic");

    Integer nrCoils = coilDomain_.size();

    if (nrCoils)
      {
        coilDef_.resize(nrCoils);
        
  
        for (Integer i=0; i < nrCoils; i++)
          {
            conf->get("iDir_"+coilDomain_[i], coilDef_[i].iDir,pdename_);
            conf->get("current_"+coilDomain_[i] ,coilDef_[i].current, pdename_);

            conf->get("coilArea_"+coilDomain_[i], coilDef_[i].coilArea, pdename_);
            // if a vector is given for the flow direction of the current
            if (coilDef_[i].iDir > 3)
              conf->getlist("coilMidPoint_"+coilDomain_[i], coilDef_[i].coilMidPt, pdename_);
            if (analysistype_==HARMONIC)
              {
                conf->get("currentPhase_"+coilDomain_[i], coilDef_[i].currentPhase, pdename_);
                coilDef_[i].currentPhase *= PI / 180;
              }


            //          conf->get("iDir", coilDef_[i].iDir, coilDomain_[i], pdename_);
            //          conf->get("current", coilDef_[i].current, coilDomain_[i], pdename_);
            //          conf->get("coilArea", coilDef_[i].coilArea, coilDomain_[i], pdename_);
            //          // if a vector is given for the flow direction of the current
            //          if (coilDef_[i].iDir > 3)
            //            conf->getlist("coilMidPoint", coilDef_[i].coilMidPt, coilDomain_[i], pdename_);
            //          if (analysistype_==HARMONIC)
            //            {
            //              conf->get("currentPhase", coilDef_[i].currentPhase, coilDomain_[i], pdename_);
            //              coilDef_[i].currentPhase *= PI / 180;
            //            }
          

            Info->PrintCoil(coilDomain_[i], coilDef_[i], analysistype_);
          }
      }  
  }



  void MagEdgePDE::CorrectEdgeDir(Matrix<Double> & elemmat, std::vector<Integer>& esign)
  {
#ifdef TRACE
    (*trace) << "entering MagEdgePDE::CorrectEdgeDir" << std::endl;
#endif  

    // correct sign of entries in elemmat due to orientation of edge
    for(Integer ii=0; ii<elemmat.getSize(); ii++)
      for(Integer jj=0; jj<elemmat.getSize(); jj++)
        elemmat[ii][jj] *= esign[ii] * esign[jj];
  }




  MagEdgePDE::~MagEdgePDE()
  {
    ;
  }

} // end of namespace


