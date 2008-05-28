#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <iomanip>
#include <stdarg.h>
#include <list>
#include <math.h>
#include <vector>

#include "elemIntegr.hh"


using namespace CoupledField;

  void  calcacousrcterms3d_quick_( const std::vector<int>& Topology, const std::vector<double>& NodalCoords, 
                                   const std::vector<double>& velocity, std::vector<double>& acouSrc, int nElems) {

#ifdef TRACE
    (*trace) << " entering CALCACOUSRCTERMS3D_QUICK_" << std::endl;
#endif


    const int ELEMNUMNOD=8;
    const int DIM=3;

    Matrix<Double> coordMat;
    Matrix<Double> NodaldTijdxj;
    Matrix<Double> NodalVel;

    coordMat.Resize(DIM, ELEMNUMNOD);
    NodaldTijdxj.Resize(DIM, ELEMNUMNOD);
    NodalVel.Resize(DIM, ELEMNUMNOD);    
    Vector<Double> elemvec4FASTEST;

//     Vector<Double> maxbnd;
//     Vector<Double> minbnd;
//     maxbnd.Resize(DIM);
//     minbnd.Resize(DIM);

//     minbnd[0]=0.002;
//     minbnd[1]=0.002; 
//     minbnd[2]=0.000;
//     maxbnd[0]=0.798;
//     maxbnd[1]=0.218; 
//     maxbnd[2]=0.218;

    ElemIntegr * ptElemIntegr=new ElemIntegr(ET_HEXA8);
  
        int k=0;
        for( int i=0; i<nElems; i++)
          {
            for( int n=0; n<ELEMNUMNOD; n++)
              {
                int idx = (Topology[k+n]-1)*3;

                for (int dimIdx=0;dimIdx<DIM;dimIdx++)
                  {
                    coordMat[dimIdx][n] = NodalCoords[idx+dimIdx];
                    NodalVel[dimIdx][n] = velocity[idx+dimIdx];
                    NodaldTijdxj[dimIdx][n] = 0;
                  }
              }
            
            ptElemIntegr->PerformIntegration(coordMat, NodaldTijdxj, NodalVel, elemvec4FASTEST);
            
            // Here we assemble the values contained in elemvec4FASTEST to the global acouSrc vector                 
            for( int n=0; n<ELEMNUMNOD; n++)
              {
                acouSrc[Topology[k+n]-1] += elemvec4FASTEST[n];
              }

            //nextElem:
            k += ELEMNUMNOD;
          }
        
        //delete objects
        if (ptElemIntegr) delete ptElemIntegr;
  
  }

//  void calcacousrcterms3d_withPress_quick_( int Topology[8], double NodalCoords[24], 
//                                            double pressure[8], double acouSrc[8], int nElems) {
  void  calcacousrcterms3d_withPress_quick_( const std::vector<int>& Topology, const std::vector<double>& NodalCoords, 
                                    const std::vector<double>& pressure, std::vector<double>& acouSrc, int nElems) {

#ifdef TRACE
    (*trace) << " entering CALCACOUSRCTERMS3D_WITHPRESS_QUICK_" << std::endl;
#endif

    const int ELEMNUMNOD=8;
    const int DIM=3;

    Matrix<Double> coordMat;
    Matrix<Double> NodaldTijdxj;
    Matrix<Double> NodalPress;

    coordMat.Resize(DIM, ELEMNUMNOD);
    NodaldTijdxj.Resize(DIM, ELEMNUMNOD);
    NodalPress.Resize(1, ELEMNUMNOD);    
    Vector<Double> elemvec4FASTEST;
//     Vector<Double> maxbnd;
//     Vector<Double> minbnd;
//     maxbnd.Resize(DIM);
//     minbnd.Resize(DIM);

//     minbnd[0]=0.001;
//     minbnd[1]=0.001; 
//     minbnd[2]=0.000;
//     maxbnd[0]=0.799;
//     maxbnd[1]=0.219; 
//     maxbnd[2]=0.219;

    ElemIntegr * ptElemIntegr=new ElemIntegr(ET_HEXA8);
  
        int k=0;
        for( int i=0; i<nElems; i++)
          {
            for( int n=0; n<ELEMNUMNOD; n++)
              {
                int idx = (Topology[k+n]-1)*3;

                for (int dimIdx=0;dimIdx<DIM;dimIdx++)
                  {
                    //      if ((NodalCoords[idx+dimIdx]>maxbnd[dimIdx])||
                    //    (NodalCoords[idx+dimIdx]<minbnd[dimIdx]))
                    //  {
                    //    goto nextElem;
                    //  }
                    coordMat[dimIdx][n] = NodalCoords[idx+dimIdx];
                    NodaldTijdxj[dimIdx][n] = 0;
                  }
         
                    NodalPress[0][n] = pressure[Topology[k+n]-1];
              }
            
            ptElemIntegr->PerformIntegration(coordMat, NodaldTijdxj, NodalPress, elemvec4FASTEST);
            
            // Here we assemble the values contained in elemvec4FASTEST to the global acouSrc vector                 
            for( int n=0; n<ELEMNUMNOD; n++)
              {
                acouSrc[Topology[k+n]-1] += elemvec4FASTEST[n];
              }




            // nextElem:
            k += ELEMNUMNOD;
          }

        //delete objects
        if (ptElemIntegr) delete ptElemIntegr;
  
  }


extern "C" 
{


  void  calcacousrcterms3d_( double coordMat_fortran[3][8], double Nodal_dTijdxj_fortran[3][8],
                           double Nodal_Vel_fortran[3][8], double elemVec[8]) {
#ifdef TRACE
    (*trace) << " entering CALC_ACOU_SOURCE_TERMS_3D_" << std::endl;
#endif


    const int ELEMNUMNOD=8;
    const int DIM=3;

    Matrix<Double> coordMat;
    Matrix<Double> NodaldTijdxj;
    Matrix<Double> NodalVel;

    coordMat.Resize(DIM, ELEMNUMNOD);
    NodaldTijdxj.Resize(DIM, ELEMNUMNOD);
    NodalVel.Resize(DIM, ELEMNUMNOD);    

    
    for (int n=0; n<ELEMNUMNOD; n++)
    {
	for (int k=0; k<DIM; k++)
        {
	    coordMat[k][n]=coordMat_fortran[k][n];
	    NodaldTijdxj[k][n]=Nodal_dTijdxj_fortran[k][n];
	    NodalVel[k][n]=Nodal_Vel_fortran[k][n];
        }
    }

    // Hex1-element
    ElemIntegr * ptElemIntegr=new ElemIntegr(ET_HEXA8);

    Vector<Double> elemvec4FASTEST;

    ptElemIntegr->PerformIntegration(coordMat, NodaldTijdxj, NodalVel, elemvec4FASTEST);
  
    // Here we send back the values contained in elemvec4FASTEST to FASTEST in the vector elemVec
    for (int n=0; n<ELEMNUMNOD; n++)
    {
        elemVec[n]=elemvec4FASTEST[n];
    }

    //delete objects
    if (ptElemIntegr) delete ptElemIntegr;

}



void  calcacousrcterms2d_( double coordMat_fortran[2][4], double Nodal_dTijdxj_fortran[2][4],
                           double Nodal_Vel_fortran[2][4], double elemVec[4]) {
#ifdef TRACE
    (*trace) << " CALC_ACOU_SOURCE_TERMS_2D_" << std::endl;
#endif


    const int ELEMNUMNOD=4;
    const int DIM=2;

    Matrix<Double> coordMat;
    Matrix<Double> NodaldTijdxj;
    Matrix<Double> NodalVel;

    coordMat.Resize(DIM, ELEMNUMNOD);
    NodaldTijdxj.Resize(DIM, ELEMNUMNOD);
    NodalVel.Resize(DIM, ELEMNUMNOD);    

    
    for (int n=0; n<ELEMNUMNOD; n++)
    {
	for (int k=0; k<DIM; k++)
        {
	    coordMat[k][n]=coordMat_fortran[k][n];
	    NodaldTijdxj[k][n]=Nodal_dTijdxj_fortran[k][n];
	    NodalVel[k][n]=Nodal_Vel_fortran[k][n];
        }
    }

    // quad1FE
    ElemIntegr * ptElemIntegr=new ElemIntegr(ET_QUAD4);

    Vector<Double> elemvec4FASTEST;
  
    ptElemIntegr->PerformIntegration(coordMat, NodaldTijdxj, NodalVel, elemvec4FASTEST);
  
    // Here we send back the values contained in elemvec4FASTEST to FASTEST in the vector elemVec
    for (int n=0; n<ELEMNUMNOD; n++)
    {
        elemVec[n]=elemvec4FASTEST[n];
    }
   
    //delete objects
    if (ptElemIntegr) delete ptElemIntegr;

}


void  calcacousrcterms2d_withPress_( double coordMat_fortran[2][4], double Nodal_dTijdxj_fortran[2][4],
                           double Nodal_Press_fortran[1][4], double elemVec[4]) {
#ifdef TRACE
    (*trace) << " CALC_ACOU_SOURCE_TERMS_2D_" << std::endl;
#endif


    const int ELEMNUMNOD=4;
    const int DIM=2;

    Matrix<Double> coordMat;
    Matrix<Double> NodaldTijdxj;
    Matrix<Double> NodalPress;

    coordMat.Resize(DIM, ELEMNUMNOD);
    NodaldTijdxj.Resize(DIM, ELEMNUMNOD);
    NodalPress.Resize(1, ELEMNUMNOD);    

    
    for (int n=0; n<ELEMNUMNOD; n++)
    {
	for (int k=0; k<DIM; k++)
        {
	    coordMat[k][n]=coordMat_fortran[k][n];
	    NodaldTijdxj[k][n]=Nodal_dTijdxj_fortran[k][n];
        }
        NodalPress[0][n]=Nodal_Press_fortran[0][n];
        //    std::cout<<"NodalPress[0]["<<n<<"]: "<<NodalPress[0][n]<<std::endl;
    }

    // quad1FE
    ElemIntegr * ptElemIntegr=new ElemIntegr(ET_QUAD4);

    Vector<Double> elemvec4FASTEST;
  
    ptElemIntegr->PerformIntegration(coordMat, NodaldTijdxj, NodalPress, elemvec4FASTEST);
  
    // Here we send back the values contained in elemvec4FASTEST to FASTEST in the vector elemVec
    for (int n=0; n<ELEMNUMNOD; n++)
    {
        elemVec[n]=elemvec4FASTEST[n];
    }
    //    std::cout<<"elemvec4FASTEST: "<<elemvec4FASTEST[0]<<" "<<elemvec4FASTEST[1]<<" "<<elemvec4FASTEST[2]<<" "<<elemvec4FASTEST[3]<<" "<<std::endl;
    //delete objects
    if (ptElemIntegr) delete ptElemIntegr;

}

//======================= Combustion Noise ===================================================

void getNumNodesAndDim(int type, int& numNodes, int& dim) {

  switch(type) 
    {
    case ET_LINE2:
      numNodes = 2;
      dim = 1;
      break;

    case ET_TRIA3:
      numNodes = 3;
      dim = 2;
      break;

    case ET_QUAD4:
      numNodes = 4;
      dim = 2;
      break;

    case ET_TET4:
      numNodes = 4;
      dim = 3;
      break;

    case ET_HEXA8:
      numNodes = 8;
      dim = 3;
      break;

    case ET_PYRA5:
      numNodes = 5;
      dim = 3;
      break;

    case ET_WEDGE6:
      numNodes = 6;
      dim = 3;
      break;

    default:
      Error("Element-Type not defined!",__FILE__,__LINE__);
    }
  
  
}


void  calccombustionsrctij_( double coordMat_fortran[24], double Nodal_Vel_fortran[24],
			     double Nodal_Rho_fortran[8], double elemVec[8], 
			     int elemType_fortran[1]) {
#ifdef TRACE
  (*trace) << " entering CALC_COMBUSTION_SOURCE_Tij_" << std::endl;
#endif

    int elemType = elemType_fortran[0];
    int ELEMNUMNOD = 0;
    int DIM = 0;

    getNumNodesAndDim(elemType, ELEMNUMNOD, DIM);

    Matrix<Double> coordMat;
    Matrix<Double> NodalVel;
    Vector<Double> NodalRho;
    coordMat.Resize(DIM, ELEMNUMNOD);
    NodalVel.Resize(DIM, ELEMNUMNOD);    
    NodalRho.Resize(ELEMNUMNOD);    

    int idx = 0;    
    for (int n=0; n<ELEMNUMNOD; n++) {
      for (int k=0; k<DIM; k++) {
	coordMat[k][n]=coordMat_fortran[idx];
	NodalVel[k][n]=Nodal_Vel_fortran[idx];
	idx++;
      }
      NodalRho[n]=Nodal_Rho_fortran[n];
    }

    ElemIntegr * ptElemIntegr=new ElemIntegr(elemType);

    Vector<Double> elemvec4CFD;
  
    ptElemIntegr->ComputeFromCombustionTij(coordMat, NodalVel, NodalRho, elemvec4CFD);
  
    // Here we send back the values contained in elemvec4CFD to CFD-Code in the vector elemVec
    for (int n=0; n<ELEMNUMNOD; n++) {
      elemVec[n]=elemvec4CFD[n];
    }
   
    //delete objects
    if (ptElemIntegr) delete ptElemIntegr;

}


void  calccombustionsrcvector_( double coordMat_fortran[24], double Nodal_Vec_fortran[24],
				double elemVec[8], int elemType_fortran[1]) {
#ifdef TRACE
    (*trace) << " entering CALC_COMBUSTION_SOURCE_Vector_" << std::endl;
#endif

    //
    // Nodal_Vec_fortran: \partial_tau ( (1-density) velocity )
    //

    int elemType = elemType_fortran[0];
    int ELEMNUMNOD = 0;
    int DIM = 0;

    getNumNodesAndDim(elemType, ELEMNUMNOD, DIM);

    Matrix<Double> coordMat;
    Matrix<Double> NodalVec;
    coordMat.Resize(DIM, ELEMNUMNOD);
    NodalVec.Resize(DIM, ELEMNUMNOD);    
    
    int idx = 0;
    for (int n=0; n<ELEMNUMNOD; n++) {
      for (int k=0; k<DIM; k++) {
	coordMat[k][n]=coordMat_fortran[idx];
	NodalVec[k][n]=Nodal_Vec_fortran[idx];
	idx++;
      }
    }

    ElemIntegr * ptElemIntegr=new ElemIntegr(elemType);

    Vector<Double> elemvec4CFD;
  
    ptElemIntegr->ComputeFromCombustionVector(coordMat, NodalVec, elemvec4CFD);
  
    // Here we send back the values contained in elemvec4CFD to CFD-Code in the vector elemVec
    for (int n=0; n<ELEMNUMNOD; n++) {
      elemVec[n]=elemvec4CFD[n];
    }
   
    //delete objects
    if (ptElemIntegr) delete ptElemIntegr;
}


void  calccombustionsrcscalar_( double coordMat_fortran[24], double Nodal_Scalar_fortran[8],
				double elemVec[8], int elemType_fortran[1]) {
#ifdef TRACE
    (*trace) << " entering CALC_COMBUSTION_SOURCE_Scalar_" << std::endl;
#endif

    //
    // Nodal_Scalar = Da * \partial_tau ( 1/density * \partial_C ( density omega_point ) )
    //

    int elemType = elemType_fortran[0];
    int ELEMNUMNOD = 0;
    int DIM = 0;

    getNumNodesAndDim(elemType, ELEMNUMNOD, DIM);

    Matrix<Double> coordMat;
    Vector<Double> NodalScalar;
    coordMat.Resize(DIM, ELEMNUMNOD);
    NodalScalar.Resize(ELEMNUMNOD);    

    int idx = 0;
    for (int n=0; n<ELEMNUMNOD; n++) {
      for (int k=0; k<DIM; k++) {
	coordMat[k][n]=coordMat_fortran[idx];
	idx++;
      }
      NodalScalar[n]=Nodal_Scalar_fortran[n];
    }

    ElemIntegr * ptElemIntegr=new ElemIntegr(elemType);

    Vector<Double> elemvec4CFD;
  
    ptElemIntegr->ComputeFromCombustionScalar(coordMat, NodalScalar, elemvec4CFD);
  
    // Here we send back the values contained in elemvec4CFD to CFD-Code in the vector elemVec
    for (int n=0; n<ELEMNUMNOD; n++) {
      elemVec[n]=elemvec4CFD[n];
    }
   
    //delete objects
    if (ptElemIntegr) delete ptElemIntegr;
}



void  calccombustionsrcvectorsurf_( double coordMat_fortran[12], double Nodal_Vec_fortran[12],
				    double elemVec[8], int elemType_fortran[1]) {
#ifdef TRACE
    (*trace) << " entering CALC_COMBUSTION_SOURCE_Tij_OnSurface_" << std::endl;
#endif

    //
    // Nodal_Vec = \partial_tau ( (1-density) velocity )
    //

    int elemType = elemType_fortran[0];
    int ELEMNUMNOD = 0;
    int DIM = 0;

    getNumNodesAndDim(elemType, ELEMNUMNOD, DIM);

    // due to fact, that element is a surface element
    if ( elemType == ET_QUAD4 || elemType == ET_TRIA3 ) {
      DIM = 3;
    } 
    else if ( elemType == ET_LINE2 ) {
      DIM = 2;
    }

    Matrix<Double> coordMat;
    Matrix<Double> NodalVec;
    coordMat.Resize(DIM, ELEMNUMNOD);
    NodalVec.Resize(DIM, ELEMNUMNOD);    

    int idx = 0;
    for (int n=0; n<ELEMNUMNOD; n++) {
      for (int k=0; k<DIM; k++) {
	coordMat[k][n]=coordMat_fortran[idx];
	NodalVec[k][n]=Nodal_Vec_fortran[idx];
	idx++;
      }
    }

    ElemIntegr * ptElemIntegr=new ElemIntegr(elemType);

    Vector<Double> elemvec4CFD;
  
    ptElemIntegr->ComputeFromCombustionVectorOnSurface(coordMat, NodalVec, elemvec4CFD);
  
    // Here we send back the values contained in elemvec4CFD to CFD-Code in the vector elemVec
    for (int n=0; n<ELEMNUMNOD; n++) {
      elemVec[n]=elemvec4CFD[n];
    }
   
    //delete objects
    if (ptElemIntegr) delete ptElemIntegr;

}


void  calccombustionsrctijsurf_( double coordMat_fortran[12], double Nodal_Vel_fortran[12],
			     double Nodal_Rho_fortran[8], double elemVec[8], 
			     int elemType_fortran[1]) {
#ifdef TRACE
    (*trace) << " entering CALC_COMBUSTION_SOURCE_Tij_OnSurface_" << std::endl;
#endif

    int elemType = elemType_fortran[0];
    int ELEMNUMNOD = 0;
    int DIM = 0;

    getNumNodesAndDim(elemType, ELEMNUMNOD, DIM);

    // due to fact, that element is a surface element
    if ( elemType == ET_QUAD4 || elemType == ET_TRIA3 ) {
      DIM = 3;
    } 
    else if ( elemType == ET_LINE2 ) {
      DIM = 2;
    }

    Matrix<Double> coordMat;
    Matrix<Double> NodalVel;
    Vector<Double> NodalRho;
    coordMat.Resize(DIM, ELEMNUMNOD);
    NodalVel.Resize(DIM, ELEMNUMNOD);    
    NodalRho.Resize(ELEMNUMNOD);    

    int idx = 0;
    for (int n=0; n<ELEMNUMNOD; n++) {
      for (int k=0; k<DIM; k++) {
	coordMat[k][n]=coordMat_fortran[idx];
	NodalVel[k][n]=Nodal_Vel_fortran[idx];
	idx++;
      }
      NodalRho[n]=Nodal_Rho_fortran[n];
    }

    ElemIntegr * ptElemIntegr=new ElemIntegr(elemType);

    Vector<Double> elemvec4CFD;
  
    ptElemIntegr->ComputeFromCombustionTijOnSurface(coordMat, NodalVel, NodalRho, elemvec4CFD);
  
    // Here we send back the values contained in elemvec4CFD to CFD-Code in the vector elemVec
    for (int n=0; n<ELEMNUMNOD; n++) {
      elemVec[n]=elemvec4CFD[n];
    }
   
    //delete objects
    if (ptElemIntegr) delete ptElemIntegr;

}



}//end externC
