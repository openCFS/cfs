#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "elec2dPDE.hh"
 
namespace CoupledField
{

Elec2dPDE::Elec2dPDE(Grid<Point2D> * grid, Material * aptMaterial, FileType * aptFileType)
:BasePDE(aptFileType,aptMaterial)
{
#ifdef TRACE
  (*trace) << "entering ElecPDE::ElecPDE " << std::endl;
#endif

  ptgrid  = grid;
  //  bilinear = new LaplaceInt<Point2D>(1);
}


void Elec2dPDE::SpecifySolver(Integer &solvertype, Integer &precondtype, Double &eps, Double &dampiter, 
                     Integer &maxnumit)
{

  eps      = 1e-6;                  // relative accuracy in the precond. energy
  dampiter = 0.7;                   // daming parameter for Jacobi, SSOR, ...
  maxnumit = 20;                    // max number of iterations

  solvertype = 2;                   // RICHARDSON = 1
                                    // CG         = 2

  precondtype = 3;                  // ID     = 1
                                    // JACOBI = 2
                                    // SSOR   = 3
                                    // MG     = 4
}

void Elec2dPDE::SpecifyMatrices(Integer &matrixtype, Integer * matrixsystype, Integer &graphtype, 
                                Integer &numdofpernode, Integer &numdirichlets, Integer &numconstraints)
{
  matrixtype = 1;     // NOCLASS = 0
                      // RSPARSE = 1
                      // CSPARSE = 2
                      // RBLOCK  = 3
                      // CBLOCK  = 4
                      // RFULL   = 5
                      // CFULL   = 6
                      // MIXED   = 7

  /* matrixsystype: NOTYPE     = 0
                    SYSTEM     = 1
                    STIFFNESS  = 2
                    DAMPING    = 3
                    CONVECTION = 4
                    MASS       = 5
  */

  matrixsystype[0] = 1;   // memory for the system matrix
  matrixsystype[1] = 0;   // memory for the stiffness matrix
  matrixsystype[2] = 0;   // memory for the damping matrix
  matrixsystype[3] = 0;   // memory for the convection matrix
  matrixsystype[4] = 0;   // memory for the mass matrix
  matrixsystype[5] = 0;   // external auxiliary matrix

  graphtype = 1; // NOGRAPH = 0
                 // NODE   = 1
                 // EDGE   = 2
                 // FACE   = 3 
                 // VOLUME = 4 

  numdofpernode  = 1;
  numdirichlets  = 1;
  numconstraints = 0;

}

void Elec2dPDE::CalcElemStiff()
{
#ifdef TRACE
  (*trace) << "entering Elec2dPDE::CalcElemeStiff" << std::endl;
#endif
 

}

void Elec2dPDE::AssembleGlobal(AbstractAlgebraicSys * algsys)
{
#ifdef TRACE
  (*trace) << "entering Elec2dPDE::AssembleGlobal" << std::endl;
#endif

  Integer i,ii,iii; 
  Integer irow,icln; 

  Integer numnodeelem;
mark
  numnodeelem=ptgrid->GetNumNodesPerElem(0,0);
mark
  Integer * help=new Integer[numnodeelem];
  Matrix<Double> elemmat;

  Point2D * ptCoord=new Point2D[numnodeelem];
  
  BaseElem * ptElem;
  switch(numnodeelem)
  {
    case 3:
       ptElem=new  Triangle1(GaussOrder3);
       break;

    case 4:
       ptElem=new Quad1(GaussOrder5);  
       break;

    default:
       Error("Number of nodes per element is strange",__FILE__,__LINE__);
  }

//   typeBaseForm oElemMatrix(ptElem,1);

//   Mat.Init();

//   // This part we should do for every group of elements
//   //   ptgrid->GetCoordOfNodesElem(0,0,ptCoord);
//   //   oElemMatrix.CalcElemMatrix(ptCoord, elemmat);

   Integer numelem=ptgrid->GetMaxnumElem(0); 
   std::cout << "numelem:" << numelem << std::endl;
   bilinear = new LaplaceInt<Point2D>(ptElem,1);

   std::cout << "nodesperelem: " << numnodeelem << std::endl;
   for (i=0; i<numelem; i++) 
     { 
       ptgrid->GetConnection(help,0,i,numnodeelem);
       ptgrid->GetCoordOfNodesElem(i,0,numnodeelem,ptCoord);
       for (int j=0;j<numnodeelem;j++)
       {
       std::cout << "x=" << ptCoord[j].x << " y=" << ptCoord[j].y << std::endl; 
       }
       bilinear->CalcElemMatrix(ptCoord, elemmat);
       std:: cout << "ElemStiff: "<< elemmat << std::endl;

     }
   delete [] ptCoord;
   delete [] help;

}


void Elec2dPDE::SolveStatic()
{
#ifdef TRACE
  (*trace) << "entering Elec2dPDE::SolveStatic" << std::endl;
#endif
 

}

} // end of namespace
