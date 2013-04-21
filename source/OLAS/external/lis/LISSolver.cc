#include "LISSolver.hh"

#include <string>


#include "MatVec/SparseOLASMatrix.hh"
#include "MatVec/SCRS_Matrix.hh"
#include "MatVec/CRS_Matrix.hh"
#include <sstream>
#include <stdio.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

namespace CoupledField{

  static EnumTuple lisSolverTypeTuples[] =
  {
    EnumTuple( LISSolver::NOSOLVER, "noSolver" ),
    EnumTuple( LISSolver::CG, "CG" ),
    EnumTuple( LISSolver::BICG, "BiCG"),
    EnumTuple( LISSolver::CGS, "CGS"),
    EnumTuple( LISSolver::BICGSTAB, "BiCGSTAB" ),
    EnumTuple( LISSolver::BICGSTABL, "BiCGSTABL" ),
    EnumTuple( LISSolver::GPBICG, "CPBiCG"),
    EnumTuple( LISSolver::TFQMR, "TFQMR"),
    EnumTuple( LISSolver::ORTHOMIN, "Orthomin" ),
    EnumTuple( LISSolver::GMRES, "GMRES" ),
    EnumTuple( LISSolver::JACOBI, "Jacobi"),
    EnumTuple( LISSolver::GS, "GaussSeidel"),
    EnumTuple( LISSolver::SOR, "SOR" ),
    EnumTuple( LISSolver::BICGSAFE, "BiCGSave" ),
    EnumTuple( LISSolver::CR, "CR"),
    EnumTuple( LISSolver::BICR, "BiCR"),
    EnumTuple( LISSolver::CRS, "CRS"),
    EnumTuple( LISSolver::BICRSTAB, "BiCRSTAB"),
    EnumTuple( LISSolver::GPBICR, "GPBiCR"),
    EnumTuple( LISSolver::BICRSAFE, "BiCRSAFE"),
    EnumTuple( LISSolver::FGMRES, "FGMRES"),
    EnumTuple( LISSolver::IDRS, "IDRs"),
    EnumTuple( LISSolver::MINRES, "minres")
  };

  Enum<LISSolver::LISSolverType> LISSolver::lisSolverType = \
  Enum<LISSolver::LISSolverType>("LIS Solver Types",
      sizeof(lisSolverTypeTuples) / sizeof(EnumTuple),
      lisSolverTypeTuples);

  static EnumTuple lisPrecondTypeTuples[] =
  {
    EnumTuple( LISSolver::NONE, "none" ),
    EnumTuple( LISSolver::JACOBI_PRE, "jacobi" ),
    EnumTuple( LISSolver::ILU, "iluk"),
    EnumTuple( LISSolver::SSOR, "ssor"),
    EnumTuple( LISSolver::HYBRID, "hybrid" ),
    EnumTuple( LISSolver::IS, "is" ),
    EnumTuple( LISSolver::SAINV, "sainv"),
    EnumTuple( LISSolver::SAAMG, "saamg"),
    EnumTuple( LISSolver::ILUC, "iluc" ),
    EnumTuple( LISSolver::ILUT, "ilut" )
  };

  Enum<LISSolver::LISPrecondType> LISSolver::lisPrecondType = \
  Enum<LISSolver::LISPrecondType>("LIS Preconditioner Types",
      sizeof(lisPrecondTypeTuples) / sizeof(EnumTuple),
      lisPrecondTypeTuples);


LISSolver::LISSolver(PtrParamNode param, PtrParamNode olasInfo, BaseMatrix::EntryType type){

  int argc = -1;                            /* dummy arg count */
  char **argv = NULL;              /* dummy arg */
  Integer err = 0;
  err = lis_initialize(&argc,&argv); CHKERR(err);



  infoNode_ = olasInfo;
  xml_ = param;

}

LISSolver::~LISSolver(){
  Integer err = 0;
  err = lis_vector_destroy(x_); CHKERR(err);
  err = lis_vector_destroy(b_); CHKERR(err);

  err = lis_matrix_destroy(A_); CHKERR(err);

  lis_finalize();
}

void LISSolver::Setup(BaseMatrix &sysmat, PtrParamNode analysis_id){

  const StdMatrix& stdmat = dynamic_cast<const StdMatrix&>(sysmat);
  
  BaseMatrix::EntryType etype = stdmat.GetEntryType();
  BaseMatrix::StorageType stype = stdmat.GetStorageType();

  Integer err=0;

  //we create the matrix...
  err = lis_matrix_create(0,&A_); CHKERR(err);

  if(stype == BaseMatrix::SPARSE_SYM ){
    //const SCRS_Matrix<Double>& scrs = dynamic_cast<const SCRS_Matrix<Double>&>(som);
    //ok we need to think about the matrix conversion a smart way would be nice...
    EXCEPTION("LIS solver cannot yet handle SCRS matrices.");
  }else{
    if ( etype == BaseMatrix::DOUBLE ) {
      // non-symmteric real case
      const CRS_Matrix<Double>& crs = dynamic_cast<const CRS_Matrix<Double>&>(stdmat);

      if(crs.GetNumCols() != crs.GetNumRows()){
        EXCEPTION("IS solver only tested for quadratic matrices");
      }
      //gather info
      UInt nnz = crs.GetNnz();
      UInt dim = crs.GetNumRows();

      Integer * rowPtr = (Integer *)crs.GetRowPointer();
      Integer * colPtr = (Integer *)crs.GetColPointer();
      Double * dataPtr = const_cast<Double*>(crs.GetDataPointer());

      err = lis_matrix_set_size(A_,dim,0); CHKERR(err);
      err = lis_matrix_set_csr(nnz,rowPtr,colPtr,dataPtr,A_); CHKERR(err);
      err = lis_matrix_assemble(A_); CHKERR(err);

      // Create RHS vector;
      err = lis_vector_duplicate(A_,&b_); CHKERR(err);
    }
    else {
      // non-symmteric complex case
      const CRS_Matrix<Complex>& crs = dynamic_cast<const CRS_Matrix<Complex>&>(stdmat);

      if(crs.GetNumCols() != crs.GetNumRows()){
        EXCEPTION("IS solver only tested for quadratic matrices");
      }

      //gather info
      UInt dim = crs.GetNumRows();

      Integer * rowPtr = (Integer *)crs.GetRowPointer();
      Integer * colPtr = (Integer *)crs.GetColPointer();
      Complex * dataPtr = const_cast<Complex*>(crs.GetDataPointer());

      err = lis_matrix_set_size(A_,dim*2,0); CHKERR(err);
   
      err = lis_vector_create(0,&b_); CHKERR(err);
      err = lis_vector_set_size(b_,0,dim*2); CHKERR(err);
      
      for(UInt row=0; row<dim; row++) {
        for(Integer col=rowPtr[row]; col<rowPtr[row+1]; col++) {
#if 0	  
          LIS_INT i=row*2;
          LIS_INT j=colPtr[col]*2;
#endif	  
          LIS_INT i=row;
          LIS_INT j=colPtr[col];
          Complex val=dataPtr[col];
#if 0	  
          if(val.real()) lis_matrix_set_value(LIS_INS_VALUE,  i,  j, val.real(),A_);
          if(val.imag()) lis_matrix_set_value(LIS_INS_VALUE,i+1,  j,-val.imag(),A_);
          if(val.imag()) lis_matrix_set_value(LIS_INS_VALUE,  i,j+1, val.imag(),A_);
          if(val.real()) lis_matrix_set_value(LIS_INS_VALUE,i+1,j+1, val.real(),A_);
#endif	  
          if(val.real()) lis_matrix_set_value(LIS_INS_VALUE,  i,  j, val.real(),A_);
          if(val.imag())
          {
            lis_matrix_set_value(LIS_INS_VALUE,i+dim,  j,-val.imag(),A_);
            lis_matrix_set_value(LIS_INS_VALUE,  i,j+dim, val.imag(),A_);
          }
          else
          {
            lis_matrix_set_value(LIS_INS_VALUE,i+dim,  j,-1e0,A_);
            lis_matrix_set_value(LIS_INS_VALUE,  i,j+dim, 1e0,A_);
          }
          if(val.real()) lis_matrix_set_value(LIS_INS_VALUE,i+dim,j+dim, val.real(),A_);
        }
      }
      err = lis_matrix_set_type(A_,LIS_MATRIX_CSR); CHKERR(err);
      err = lis_matrix_assemble(A_); CHKERR(err);
    }
  }
  err = lis_vector_duplicate(b_,&x_); CHKERR(err);

  lis_vector_set_all(1.0,x_);

  //create the solver

  std::string config;
  createConfigString(xml_,config);
  err = lis_solver_create(&solver_); CHKERR(err);
  err = lis_solver_set_option(const_cast<char*>(config.c_str()),solver_);CHKERR(err);

  solver_->A = A_;
  err = lis_precon_create(solver_, &precond_);
  CHKERR(err);
  if( err ){
    std::cerr << "There was an error creating the preconditioner. Code: " << err << " ...Going to abort" << std::endl;
    EXCEPTION("ERROR");
  }


}

void LISSolver::Solve( const BaseMatrix &sysmat,
                       const BaseVector &rhs, BaseVector &sol, PtrParamNode analysis_id){
  ParamNode::ActionType at = ParamNode::APPEND;
  PtrParamNode out = infoNode_->Get(ParamNode::PN_PROCESS)->Get("solver", at);

  if(sysmat.GetEntryType() == BaseMatrix::DOUBLE) {
    for(Integer i=0, n=(Integer)rhs.GetSize(); i<n; i++){
      Double myEnt =0;
      rhs.GetEntry((UInt)i,myEnt);
    
      b_->value[i] = myEnt;
      //lis_vector_set_value(LIS_INS_VALUE,i,myEnt,b_);
    }
  } else {
    for(Integer i=0, n=(Integer)rhs.GetSize(); i<n; i++){
      Complex myEnt = 0;
      rhs.GetEntry((UInt)i,myEnt);

      // Set RHS value
#if 0
      lis_vector_set_value(LIS_INS_VALUE,i*2,myEnt.real(),b_);
      lis_vector_set_value(LIS_INS_VALUE,i*2+1,myEnt.imag(),b_);
#endif
      lis_vector_set_value(LIS_INS_VALUE,i,myEnt.real(),b_);
      lis_vector_set_value(LIS_INS_VALUE,i+n,myEnt.imag(),b_);
    }
  }
  
  Integer err = 0;
  std::cerr << "SOLVING" << std::endl;
  err = lis_solve_kernel(A_, b_, x_, solver_, precond_);
  if(err){
    EXCEPTION("Solver returned error code: " << err << " ...aborting")
  }
  std::cerr << "DONE" << std::endl;
  //lis_solve(A_, b_, x_, solver_);
  //copy solution
  if(sysmat.GetEntryType() == BaseMatrix::DOUBLE) {
    for(UInt i=0; i<sol.GetSize();i++){
      sol.SetEntry(i,x_->value[i]);
    }
  }
  else
  {
    for(UInt i=0, n=sol.GetSize(); i<n; i++){
      Complex myEnt = 0;
#if 0
      myEnt.real() = x_->value[i*2];
      myEnt.imag() = x_->value[i*2+1];
#endif
      myEnt = Complex(x_->value[i], x_->value[i+n]);
      
      sol.SetEntry(i,myEnt);
    }
  }
  Integer iterations;
  Double lastTime;
  Double norm;
  Integer solverCode;
  lis_solver_get_iters(solver_,&iterations);
  lis_solver_get_time(solver_,&lastTime);
  lis_solver_get_residualnorm(solver_,&norm);
  lis_solver_get_solver(solver_,&solverCode);


  out->Get("iterations")->SetValue(iterations);
  PtrParamNode timing = out->Get("timing");
  timing->Get("lastRun")->SetValue(lastTime);
  PtrParamNode norms = out->Get("norms");
  norms->Get("residualNorm")->SetValue(norm);
  out->Get("solverType")->SetValue(lisSolverType.ToString((LISSolverType)solverCode));
}


void LISSolver::createConfigString(PtrParamNode configNode, std::string& output){
  std::string solStr;
  std::string precondStr;
  createSolverString(configNode,solStr);
  createPrecondString(configNode,precondStr);

  std::stringstream globStream;

  //lets read teh basics
  PtrParamNode cNode = configNode->Get("maxIter",ParamNode::PASS);
  if(cNode){
    UInt maxIter = 1000;
    configNode->GetValue("maxIter",maxIter);
    globStream << " -maxiter " << maxIter;
  }

  cNode = configNode->Get("tolerance",ParamNode::PASS);
  if(cNode){
    Double tol = 1e-12;
    configNode->GetValue("tolerance",tol);
    globStream << " -tol " << tol;
  }

  cNode = configNode->Get("logging",ParamNode::PASS);
  if(cNode){
    bool log = false;
    configNode->GetValue("logging",log);
    if(log){
      globStream << " -print out";
    }else{
      globStream << " -print none";
    }
  }

  output = solStr + " " + precondStr + " " + globStream.str() + " -initx_ones false -initx_zeros false";
  std::cout << " the config string for LIS was: " << output << std::endl;
#ifdef _OPENMP
  std::cout << "max number of threads = " << omp_get_num_procs() << std::endl;
  std::cout << "number of threads = " << omp_get_max_threads() << std::endl;
#endif
  return;
}


void LISSolver::createSolverString(PtrParamNode solverNode, std::string& output){
  std::string solverString;
  std::string nodeName = solverNode->GetName();
  PtrParamNode sNode = solverNode->Get("solver",ParamNode::PASS);

  std::stringstream solstream;

  if(sNode){
    ParamNodeList sol = sNode->GetChildren();

    UInt numChilds = sol.GetSize();
    //new know that first child is id the other is the actual solver tag...
    //well, not very filsafe but good for now.
    if(numChilds==2){
      solverString = sol[1]->GetName();
      LISSolverType type;
      type = lisSolverType.Parse(solverString);
      UInt typeID = type;

      solstream << " -i "<< typeID;

      PtrParamNode curNode;
      UInt deg = 0;
      UInt restart = 0;
      Double omega = 0;
      switch(type){
      case BICGSTABL:
        deg = 0;
        sol[1]->GetValue("degree",deg,ParamNode::PASS);
        solstream << " -ell "<< deg;
        break;
      case ORTHOMIN:
        restart = 0;
        sol[1]->GetValue("restart",restart,ParamNode::PASS);
        solstream << " -restart "<< restart;
        break;
      case GMRES:
        restart = 0;
        sol[1]->GetValue("restart",restart,ParamNode::PASS);
        solstream << " -restart "<< restart;
        break;
      case SOR:
        omega = 0;
        sol[1]->GetValue("omega",omega,ParamNode::PASS);
        solstream << " -omega "<< omega;
        break;
      case FGMRES:
        restart = 0;
        sol[1]->GetValue("restart",restart,ParamNode::PASS);
        solstream << " -restart "<< restart;
        break;
      case IDRS:
        restart = 0;
        sol[1]->GetValue("irestart",restart,ParamNode::PASS);
        solstream << " -irestart "<< restart;
        break;
      case CG:
        //nothing to do
        break;
      case BICG:
        //nothing to do
        break;
      case CGS:
        //nothing to do
        break;
      case BICGSTAB:
        //nothing to do
        break;
      case GPBICG:
        //nothing to do
        break;
      case TFQMR:
        //nothing to do
        break;
      case JACOBI:
        //nothing to do
        break;
      case GS:
        //nothing to do
        break;
      case BICGSAFE:
        //nothing to do
        break;
      case CR:
        //nothing to do
        break;
      case BICR:
        //nothing to do
        break;
      case CRS:
        //nothing to do
        break;
      case BICRSTAB:
        //nothing to do
        break;
      case GPBICR:
        //nothing to do
        break;
      case BICRSAFE:
        //nothing to do
        break;
      case MINRES:
        //nothing to do
        break;
      default:
        EXCEPTION("Could not obtain a valid solver");
        break;
      }

    }else if(numChilds > 2){
      EXCEPTION("Expected 2 Childs of element solver got " << numChilds << "!");
    }
  }

  output = solstream.str();
}


void LISSolver::createPrecondString(PtrParamNode precondNode, std::string& output){
  std::string precondString;

  PtrParamNode sNode = precondNode->Get("precond",ParamNode::PASS);

  std::stringstream solstream;

  if(sNode){
    ParamNodeList sol = sNode->GetChildren();

    UInt numChilds = sol.GetSize();
    //new know that first child is id the other is the actual solver tag...
    //well, not very filsafe but good for now.
    if(numChilds==2){
      precondString = sol[1]->GetName();
      LISPrecondType type;
      type = lisPrecondType.Parse(precondString);
      UInt typeID = type;

      solstream << " -p "<< typeID;

      PtrParamNode curNode;
      Double fill = 0;
      UInt restart = 1;
      Double alpha = 1.0;
      UInt m = 3;
      Double drop = 0.05;
      bool unsym = 0;
      Double theta = 0;
      Double fillInRate = 5.0;
      switch(type){
      case ILU:
        fill = 0;
        sol[1]->GetValue("fill",fill,ParamNode::PASS);
        solstream << " -ilu_fill "<< fill;
        break;
      case SSOR:;
        restart = 1;
        sol[1]->GetValue("restart",restart,ParamNode::PASS);
        solstream << " -ssor_w "<< restart;
        break;
      case HYBRID:
        EXCEPTION("Hybrid preconditioner not supported right now...")
        break;
      case IS:
        alpha = 1.0;
        sol[1]->GetValue("alpha",alpha,ParamNode::PASS);
        solstream << " -is_alpha "<< alpha;
        m = 3;
        sol[1]->GetValue("m",m,ParamNode::PASS);
        solstream << " -is_m "<< alpha;
        break;
      case SAINV:
        drop = 0;
        sol[1]->GetValue("drop",drop,ParamNode::PASS);
        solstream << " -sainv_drop "<< drop;
        break;
      case SAAMG:
        unsym = 0;
        theta = 0;
        sol[1]->GetValue("unsym",unsym,ParamNode::PASS);
        sol[1]->GetValue("theta",theta,ParamNode::PASS);
        solstream << " -saamg_unsym "<< unsym;
        solstream << " -saamg_theta "<< theta;
        break;
      case ILUC:
        drop = 0.05;
        fillInRate = 5.0;
        sol[1]->GetValue("drop",drop,ParamNode::PASS);
        sol[1]->GetValue("fillInRate",fillInRate,ParamNode::PASS);
        solstream << " -iluc_drop "<< drop;
        solstream << " -iluc_rate "<< fillInRate;
        break;
      case ILUT:
        drop = 0.05;
        fillInRate = 5.0;
        sol[1]->GetValue("drop",drop,ParamNode::PASS);
        sol[1]->GetValue("fillInRate",fillInRate,ParamNode::PASS);
        solstream << " -ilut_drop "<< drop;
        solstream << " -ilut_rate "<< fillInRate;
        break;
      case NONE:
      case JACOBI_PRE:
        //nothing to do
        break;
      default:
        EXCEPTION("Could not obtain a valid preconditioner");
        break;
      }

    }else if(numChilds > 2){
      EXCEPTION("Expected 2 Childs of element solver got " << numChilds << "!");
    }
  }
  output = solstream.str();

}

}
