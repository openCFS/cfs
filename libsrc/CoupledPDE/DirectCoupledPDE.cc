#include "DirectCoupledPDE.hh"


namespace CoupledField{


DirectCoupledPDE::DirectCoupledPDE(Grid *aptgrid, BCs *aptBCs, 
				   FileType *aInFile, 
				   WriteResults *aOutFile, 
				   TimeFunc *aTimeFunc)
  : StdPDE(aptgrid, aptBCs, aInFile, aOutFile, aTimeFunc)
{
  ENTER_FCN( "DirectCoupledPDE::DirectCoupledPDE" );

}


DirectCoupledPDE::~DirectCoupledPDE()
{
  ENTER_FCN( "DirectCoupledPDE::~DirectCoupledPDE" );
  
}

void DirectCoupledPDE::WriteGeneralPDEdefines()
{

  ENTER_FCN( "DirectCoupledPDE::WriteGeneralPDEdefines" );
}

void DirectCoupledPDE::SetBCs(const Integer level, 
			      const Double atimestep)
{
  ENTER_FCN( "DirectCoupledPDE::SetBCs" );
}


// ======================================================
// POSTPROC SECTION
// ======================================================

void DirectCoupledPDE::PostProcess(const Integer level)
{
  ENTER_FCN( "DirectCoupledPDE::PostProcess" );
}


void DirectCoupledPDE::WriteResultsInFile(const Integer kstep,
					  const Double asteptime,
					  Integer stepOffset,
					  Double timeOffset)
{
  ENTER_FCN( "DirectCoupledPDE::WriteResultsInFile" );

}

// ======================================================
// COUPLING SECTION
// ======================================================

void DirectCoupledPDE::InitCoupling(PDECoupling * Coupling)
{
  ENTER_FCN( "DirectCoupledPDE::InitCoupling" );
}
void DirectCoupledPDE::ResetCoupling()
{
  ENTER_FCN( "DirectCoupledPDE::ResetCoupling ");
	     }
void DirectCoupledPDE::CalcInputCoupling()
{
  ENTER_FCN( "DirectCoupledPDE::CalcInputCoupling" );
}

void DirectCoupledPDE::CalcOutputCoupling()
{
  ENTER_FCN( "DirectCoupledPDE::CalcOutputCoupling" );
}


} // end of namesapce


