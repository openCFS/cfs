#include <string>
#include <fstream>

#include "Elements/elements_header.hh"
#include "grid.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"

namespace CoupledField
{

Grid::Grid(FileType * aptFileType)
{
  ENTER_FCN( "Grid::Grid" );

  ptFileType = aptFileType;

  ptQ1     = new Quad1FE();
  ptQ2     = new Quad2FE();
  ptTet1   = new Tetra1FE();
  ptL1     = new Line1FE();
  ptL2     = new Line2FE();
  ptTr1    = new Triangle1FE();
  ptTr2    = new Triangle2FE();
  ptHexa1  = new Hexa1FE();
  ptHexa2  = new Hexa2FE();
  ptPyra1  = new Pyra1FE();
  ptWedge1 = new Wedge1FE();
  ptWedge2 = new Wedge2FE();

  lastlevel_=0;

  params->GetList( "name", listSD_, "domain", "region" );

}



  
void Grid::SetIntTypeAllElems(IntegrationType aIntType)
{
  ENTER_FCN( "Grid::SetIntTypeAllElems" );

  ptQ1    -> SetIntegrationType(aIntType);
  ptQ2    -> SetIntegrationType(aIntType);
  ptTet1  -> SetIntegrationType(aIntType);
  ptL1    -> SetIntegrationType(aIntType);
  ptL2    -> SetIntegrationType(aIntType);
  ptTr1   -> SetIntegrationType(aIntType);
  ptTr2   -> SetIntegrationType(aIntType);
  ptHexa1 -> SetIntegrationType(aIntType);
  ptHexa2 -> SetIntegrationType(aIntType);
  ptPyra1 -> SetIntegrationType(aIntType);
  ptWedge1-> SetIntegrationType(aIntType);
  ptWedge2-> SetIntegrationType(aIntType);

  ptQ1    ->Init();
  ptQ2    ->Init();
  ptTet1  ->Init();
  ptL1    ->Init();
  ptL2    ->Init();
  ptTr1   -> Init();
  ptTr2   -> Init();
  ptHexa1 -> Init();
  ptHexa2 -> Init();
  ptPyra1 -> Init();
  ptWedge1-> Init();
  ptWedge2-> Init();

}




Grid::~Grid()
{
  ENTER_FCN( "Grid::~Grid" );
  if (ptQ1)     delete ptQ1;
  if (ptQ2)     delete ptQ2;
  if (ptTet1)   delete ptTet1;
  if (ptL1)     delete ptL1;
  if (ptL2)     delete ptL2;
  if (ptTr1)    delete ptTr1;
  if (ptTr2)    delete ptTr2;
  if (ptHexa1)  delete ptHexa1;
  if (ptHexa2)  delete ptHexa2;
  if (ptPyra1)  delete ptPyra1;
  if (ptWedge1) delete ptWedge1;
  if (ptWedge2) delete ptWedge2;
  
}


void Grid::CalcSurfNormalOutOfVol(Vector<Double> & n,
				  const Elem & surfElem,
				  const Elem & volElem)
{
  ENTER_IFCN( "Grid::CalcSurfNormalOutOfVol" );

  //compute normal vector
  Matrix<Double>  ptVolCoord, ptSurfCoord;

  GetCoordNodesElemMat(volElem.connect, ptVolCoord, lastlevel_);
  GetCoordNodesElemMat(surfElem.connect, ptSurfCoord, lastlevel_);
  

  Integer surfCorners = surfElem.ptElem->GetNumCorners();
  Integer volCorners = volElem.ptElem->GetNumCorners();
 
  // Check for dimension:
  // A 2D volume element has only one face
  // -> we are in 2D
  if (volElem.ptElem->GetNumFaces() == 1) {
  
    // 1. step: compute vector perpendicular to line element
    // but without defined sign
    Double dx  = ptSurfCoord[0][1] - ptSurfCoord[0][0];
    Double dy  = ptSurfCoord[1][1] - ptSurfCoord[1][0];
    Double len = sqrt(dx*dx + dy*dy);
    if (len <= 0.0)
      Error("length of normal vector is zero!");
    n.Resize(2);
    n[0] = dy/len;
    n[1] = -dx/len;
    
    // 2. step: compute direction
    
    Integer indexNode1=-1;
    Integer indexNode2=-1;
    
    for(Integer actNode=0; actNode < volCorners; actNode++)
      {
	if (volElem.connect[actNode] == surfElem.connect[0])
	  indexNode1 = actNode;
	if (volElem.connect[actNode] == surfElem.connect[1])
	  indexNode2 = actNode;
      }
    // if not clockwise orientation of nodes (difference of node indizes is -1)    
    if (indexNode1==-1 || indexNode2==-1)
      Error("Nodes of neighbouring element not found!", __FILE__, __LINE__);
    
    
    // counterclockwise orientation of nodes (difference of node indizes is +1)
    if ( (indexNode2-indexNode1 == -1 || 
	  (indexNode2-indexNode1)- volCorners == -1 ) ) {
      n *= -1;
    }
    
    else
      // counterclockwise orientation of nodes (difference of node indizes is +1)

      if (! (indexNode2-indexNode1 == 1 || 
	     (indexNode2-indexNode1)+volCorners == 1) )
	Error("Nodes of interface don't lie beneath each other in neighbouring element!", __FILE__, __LINE__);   
  }
  
  else {
    // 1. step: compute vector perpendicular to surface element
    // but without defined sign

    //compute the two vectors in the plane
    Vector<Double> vec1(3), vec2(3);
    for (Integer i=0; i<3; i++) {
      vec1[i] = ptSurfCoord[i][1]             - ptSurfCoord[i][0];
      vec2[i] = ptSurfCoord[i][surfCorners-1] - ptSurfCoord[i][0];
    }
    //compute cross product
    n.Resize(3);
    n[0] = vec1[1] * vec2[2] - vec1[2]*vec2[1];
    n[1] = vec1[2] * vec2[0] - vec1[0]*vec2[2];
    n[2] = vec1[0] * vec2[1] - vec1[1]*vec2[0];
    //normalize the length to 1
    Double length = n.NormL2();
    n /= length;
   
    // 2. step: compute direction
    
    // find first common vertex index
    Integer firstCommonIndex = -1;
    for (Integer i=0; i<volCorners; i++)
      if (volElem.connect[i] == surfElem.connect[0]){
	firstCommonIndex = i;
	break;
      }
    
    // calculate barycenter of volume element
    Vector<Double> barycenter(3);
    for (Integer i=0; i<volCorners; i++){
      barycenter[0] += ptVolCoord[0][i];
      barycenter[1] += ptVolCoord[1][i];
      barycenter[2] += ptVolCoord[2][i];
    }

    barycenter /= volCorners;

    // check, if scalar product with vector (going from barycenter to
    // common edge) and perpendicular vector  are pointing in same direction
    Vector<Double> innerVec(3);
    Double product = 0;
    innerVec[0] = ptVolCoord[0][firstCommonIndex] - barycenter[0];
    innerVec[1] = ptVolCoord[1][firstCommonIndex] - barycenter[1];
    innerVec[2] = ptVolCoord[2][firstCommonIndex] - barycenter[2];
    
    product = innerVec * n;
    if (product < 0) {
      n *= -1;
    }

  }
  
}

}
