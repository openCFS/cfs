//#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>

#include <general_head.hh>
#include <utils_head.hh>
#include <datainout_head.hh>
#include <elements_head.hh>
#include <forms_head.hh>
#include <domain_head.hh>
#include "systemmatrix.hh"
#include "assemble.hh"

namespace CoupledField
{

template<class Dim, class T_Matrix>
Assemble<Dim, T_Matrix>::Assemble(Grid<Dim> * aptgrid)
{
#ifdef TRACE
   (*trace) << "entering Assemble::Assemble" << std::endl; 
#endif
   
   ptgrid=aptgrid;
   size=ptgrid->GetMaxnumnodes(0);
   if (InfoPrint) 
      if (A.IsSymmetric())
        (*infofile) << "we are working with symmetric matrix" << std::endl;
   A.Resize(size,size);
   M.Resize(size,size);
   S.Resize(size,size);
   b.Resize(size);

   IsCalcM=FALSE;
   IsCalcS=FALSE;
}

template<class Dim, class T_Matrix>
Assemble<Dim, T_Matrix>::~Assemble()
{
#ifdef TRACE
   (*trace) << "entering Assemble::~Assemble" << std::endl;
#endif
  ;
}

template<class Dim, class T_Matrix>
void Assemble<Dim, T_Matrix>::SetAb()
{
#ifdef TRACE
   (*trace) << "entering Assemble::SetAb" << std::endl;
#endif

  Integer n=5;
  A.Resize(n,n);
  b.Resize(n);

  Integer i,j;
  for (i=0; i<n; i++)
   for (j=0; j<n; j++)
     A[i][j]=1;
  
  for (i=0; i<n; i++) b[i]=n;
}

template<class Dim, class T_Matrix>
void Assemble<Dim,T_Matrix>::AssembleSysMatrix(const Double CoefL, const Double CoefM) 
{
#ifdef TRACE
  (*trace) << "entering Assemble::AssembleSysMatrix" << std::endl;
#endif
  Clock oClock;    
  oClock.ClockCount(Clock::beg);

  if (!IsCalcS) 
    { 
      AssembleGlobal< LaplaceInt<Dim> >(S);
      IsCalcS=TRUE;
#ifdef DEBUG
      (*debug) << "------- Stiffness Matrix ---------" << std::endl << S;
#endif
    }

  if (!IsCalcM) 
    {
      IsCalcM=TRUE;
      AssembleGlobal< MassInt<Dim> >(M);
#ifdef DEBUG
      (*debug) << "------- Mass Matrix ---------" << std::endl << M;
#endif
    }

  if (CoefM!=1.0) 
    A+=M*CoefM;
  else A+=M;

  if (CoefL!=1.0) 
    A+=S*CoefL;
  else A+=S;

  oClock.ClockCount(Clock::end, "Assemble matrix");

#ifdef DEBUG
  (*debug) << "----- System Matrix --------" << std::endl;
  (*debug) << A;
#endif

#ifdef TRACE
 (*trace) << " have left Assemble::AssembleSysMat " << std::endl;
#endif

}

template<class Dim, class T_Matrix>
void Assemble<Dim, T_Matrix>::SetRHS(const Vector<Double> & CoefM, const Vector<Double> & CoefS, const Vector<Double> & f)
{
#ifdef TRACE
  (*trace) << "entering Assemble::SetRHS" << std::endl;
#endif
  if (!CoefM.size() && !CoefS.size() && !f.size()) b.Init();
  if (CoefM.size()) b=M*CoefM;
  if (CoefS.size()) 
    if (CoefM.size()) b+=S*CoefS;
    else b=S*CoefS;
  if (f.size()) if (!CoefM.size() || !CoefS.size()) b+=f; 
  else b=f;
#ifdef DEBUG
  (*debug) << " ----- Right hand side -----" << std::endl;
  (*debug) << b;
#endif
}

template<class Dim, class T_Matrix>
void Assemble<Dim, T_Matrix>::SetDirichletBoundaryCondSysMat_PenaltyMethod()
{
#ifdef TRACE
  (*trace)<<"entering Assemble::SetDirichletBoundaryCondSysMat_PenaltyMethod"<< std::endl;
#endif
  Integer nn=A.getSize();
  Integer i;
  Double max=A[0][0];
  for (i=1; i < nn; i++)
    if (A(i,i)>max) max=A(i,i);
   
  BigConst=10e+10*max;

  Integer n;
  n=nodesDirBC.size();


  Integer aux; 
  for (i=0; i<n; i++) 
    {
      aux=nodesDirBC[i]-1;
      std::cout << "number of nodes " << aux << std::endl;
      A(aux,aux)+=BigConst;
    }
}

template<class Dim, class T_Matrix>
void Assemble<Dim, T_Matrix>::SetDirichletBoundaryCondZero_Cut()
{
#ifdef TRACE
  (*trace)<<"entering Assemble::SetDirichletBoundaryCondZero_Cut"<< std::endl;
#endif
 
  Integer n;
  n=nodesDirBC.size();
 
  Integer i,aux;
  for (i=0; i<n; i++)
    {     
      aux=nodesDirBC[i]-1;
      std::cout << "nodes for BC" << aux << std::endl;
      mark
      std::cout << "aux-i" << aux-i << std::endl;
      std::cout << A << std::endl;
      A.cut(aux-i,aux-i);
mark
      b.cut(aux-i);
mark
    }
}

template<class Dim, class T_Matrix>
void Assemble<Dim,T_Matrix>::SetNodesBoundaryCondition(FileType * aptFileType)
{
#ifdef TRACE
   (*trace) << "entering Assemble::SetNodesBoundaryCondition " << std::endl;
#endif
  aptFileType->ReadDirichletBC(nodesDirBC); 
}

template<class Dim, class T_Matrix>
void Assemble<Dim, T_Matrix>::SetDirichletBoundaryCondRHS_PenaltyMethod(const Double val_tf)
{
#ifdef TRACE
  (*trace) << "entering Assemble::SetDirichletBoundaryCondRHS_PenaltyMethod" << std::endl;
#endif
 
  if (!nodesDirBC.size()) Error("Array of nodes with BC is not defined. Use somewhere function SetNodesBoundaryCondition() before this function.");

//  if (val_tf==0) Error("Don't use penalty method to define zero boundary condition");

  Integer i;
  Integer aux;
  Integer n=nodesDirBC.size();
  for (i=0; i<n; i++)
    {
      aux=nodesDirBC[i]-1;
      b[aux]+=BigConst*val_tf;
    }
}

template<class Dim, class T_Matrix>
template <class typeBaseForm>
void Assemble<Dim,T_Matrix>::AssembleGlobal(T_Matrix & Mat) const
{
#ifdef TRACE
   (*trace) << "entering Assemble::AssembleGlobal" << std::endl;
#endif

  Integer i,ii,iii; 
  Integer irow,icln; 

  Integer numnodeelem;
  numnodeelem=ptgrid->GetNumNodesPerElem(0,0);

  Integer * help=new Integer[numnodeelem];
  Matrix<Double> elemmat;

  Dim * ptCoord=new Dim[numnodeelem];
  
  BaseElem * ptElem;
  switch(numnodeelem)
  {
    case 3:
       ptElem=new  Triangle1(GaussOrder5);
       break;

    case 4:
       ptElem=new Quad1(GaussOrder5);  
       break;

    default:
       Error("Number of nodes per element is strange",__FILE__,__LINE__);
  }
//  BaseElem * ptElem=new Quad1(GaussOrder5);  /////////////////////
//  BaseElem * ptElem=new Triangle1(GaussOrder5);
  
  typeBaseForm oElemMatrix(ptElem,1);

  Mat.Init();

  // This part we should do for every group of elements
  //   ptgrid->GetCoordOfNodesElem(0,0,ptCoord);
  //   oElemMatrix.CalcElemMatrix(ptCoord, elemmat);

  Integer numelem=ptgrid->GetMaxnumElem(0); 
  for (i=0; i<numelem; i++) 
    { 
      ptgrid->GetConnection(help,0,i,numnodeelem);
      ptgrid->GetCoordOfNodesElem(i,0,numnodeelem,ptCoord);
      oElemMatrix.CalcElemMatrix(ptCoord, elemmat);

      for (ii=0; ii<numnodeelem; ii++)
	for (iii=0; iii<numnodeelem; iii++)
	  {
            irow=help[ii]-1;
            icln=help[iii]-1;
            if (irow >= icln)
	      //       Mat(irow,icln)+=elemmat[ii][iii];
	      Mat.Add(irow,icln,elemmat[ii][iii]);
	  }
    }
  delete [] ptCoord;
  delete [] help;
  /// Convertion from SymMatrix to Matrix  
  if (!Mat.IsSymmetric()) 
    {
      Integer helpn=Mat.getSize();
      for (i=0; i<helpn; i++)
	for (ii=i+1; ii < helpn; ii++)
	  Mat(i,ii)=Mat(ii,i);
    }
}

template<class Dim, class T_Matrix>
void Assemble<Dim, T_Matrix>::Restore() 
{
#ifdef TRACE
  (*trace) << "entering Assemble::Restore" << std::endl;
#endif

  Integer i,aux;
  for (i=0; i<nodesDirBC.size(); i++)
    {
      aux=nodesDirBC[i]-1;
      x.add(0,aux);
    }
}


template<class Dim, class T_Matrix>
void Assemble<Dim, T_Matrix>::Print()
{
#ifdef TRACE
   (*trace) << "entering Assemble::Print" << std::endl;
#endif
  std::cout << A;
}

} // end of namespace
