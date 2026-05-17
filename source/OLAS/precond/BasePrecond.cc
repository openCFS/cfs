#include "MatVec/StdMatrix.hh"
#include "MatVec/SBM_Matrix.hh"

#include "MatVec/BaseVector.hh"
#include "MatVec/SBM_Vector.hh"

#include "Utils/Timer.hh"
#include "BasePrecond.hh"

namespace CoupledField {

  static EnumTuple precondTypeTuples[] = 
  {
    EnumTuple( BasePrecond::NOPRECOND, "noPrecond" ),
    EnumTuple( BasePrecond::ID, "Id" ),
    
    // Initialization of Classical Preconditioners
    EnumTuple( BasePrecond::MG, "MG"),
    EnumTuple( BasePrecond::JACOBI, "Jacobi"),
    EnumTuple( BasePrecond::BLOCK_JACOBI, "BlockJacobi"),
    EnumTuple( BasePrecond::SSOR, "SSOR" ),
    EnumTuple( BasePrecond::ILU0, "ILU0" ),
    EnumTuple( BasePrecond::ILUTP, "ILUTP"),
    EnumTuple( BasePrecond::ILUK, "ILUK"),
    EnumTuple( BasePrecond::ILDL0, "ILDL0" ),
    EnumTuple( BasePrecond::ILDLK, "ILDLK" ),
    EnumTuple( BasePrecond::ILDLTP, "ILDLTP"),
    EnumTuple( BasePrecond::ILDLCN, "ILDLCN"),
    EnumTuple( BasePrecond::IC0, "IC0" ),
    EnumTuple( BasePrecond::SBM_DIAG, "SBMDiag" ),
    EnumTuple( BasePrecond::SBM_JACOBI, "SBMJacobi" ),
    // Initialization of Solver type preconditioners
    EnumTuple( BasePrecond::RICHARDSON, "richardson" ),
    EnumTuple( BasePrecond::DIAGSOLVER, "diagsolver"),
    EnumTuple( BasePrecond::CG, "cg"),
    EnumTuple( BasePrecond::GMRES, "gmres" ),
    EnumTuple( BasePrecond::MINRES, "minres" ),
    EnumTuple( BasePrecond::SYMMLQ, "symmlq"),
    EnumTuple( BasePrecond::LAPACK_LU, "lapackLU"),
    EnumTuple( BasePrecond::LAPACK_LL, "lapackLL" ),
    EnumTuple( BasePrecond::LU_SOLVER, "directLU" ),
    EnumTuple( BasePrecond::LDL_SOLVER, "directLDL"),
    EnumTuple( BasePrecond::LDL_SOLVER2, "directLDL2"),
    EnumTuple( BasePrecond::PARDISO_PRECOND, "pardiso" ),
    EnumTuple( BasePrecond::PARDISO_FACTREUSE, "pardisoFactReuse" ),
    EnumTuple( BasePrecond::CHOLMOD, "cholmod")
  };

  Enum<BasePrecond::PrecondType> BasePrecond::precondType = \
  Enum<BasePrecond::PrecondType>("Preconditioner Types",
      sizeof(precondTypeTuples) / sizeof(EnumTuple),
      precondTypeTuples); 



  
  void BasePrecond::PostInit() {
    // assert, that info node is present
    assert(infoNode_);

    std::string name = precondType.ToString(GetPrecondType());

    PtrParamNode b = infoNode_->Get(ParamNode::SUMMARY);

    setupTimer_   = b->Get("precond_" + name + "_setup/timer")->AsTimer();
    precondTimer_ = b->Get("precond_" + name + "_run/timer")->AsTimer();
  }
  
  void BasePrecond::GetPrecondSysMat( BaseMatrix& sysMat ) {
    WARN("Using fall-back default export for preconditioner matrix");
    sysMat.Init();
  }
    
  
  // ------------------------------------------------------------------------
  //  S B M   -  P R E C O N D I T I O N E R
  // ------------------------------------------------------------------------
  
  BaseSBMPrecond::BaseSBMPrecond( UInt numBlocks, PtrParamNode infoNode ) {
    numBlocks_ = numBlocks;
    infoNode_ = infoNode->Get("SBMPrecond");
  
  }
  
  BaseSBMPrecond::~BaseSBMPrecond() {
    
   
  }
  
  
  void BaseSBMPrecond::Apply( const BaseMatrix& sysmat, const BaseVector& r, 
                              BaseVector& z ) {
      const SBM_Matrix& sbmsysmat = dynamic_cast<const SBM_Matrix&>(sysmat);
      const SBM_Vector& sbmr      = dynamic_cast<const SBM_Vector&>(r);
      SBM_Vector& sbmz            = dynamic_cast<SBM_Vector&>(z);
      Apply(sbmsysmat,sbmr,sbmz);
    }
  

  void BaseSBMPrecond::Setup( BaseMatrix &A ) {
    SBM_Matrix& sbmA = dynamic_cast<SBM_Matrix&>(A);
    Setup(sbmA);
  }
 
}
