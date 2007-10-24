#ifndef INTEGLIB_H_2006
#define INTEGLIB_H_2006

#include <vector>

#if defined (__cplusplus)
extern "C" {
#endif

    void  calcacousrcterms3d_( double coordMat_fortran[3][8], double Nodal_dTijdxj_fortran[3][8],
                               double Nodal_Vel_fortran[3][8], double elemVec[8]);
    

    void  calcacousrcterms2d_( double coordMat_fortran[2][4], double Nodal_dTijdxj_fortran[2][4],
                               double Nodal_Vel_fortran[2][4], double elemVec[4]);

    void  calcacousrcterms2d_withPress_( double coordMat_fortran[2][4], double Nodal_dTijdxj_fortran[2][4],
                                         double Nodal_Press_fortran[1][4], double elemVec[4]);

    // void calcacousrcterms3d_withPress_quick_( int Topology[8], double NodalCoords[24], 
    //                                     double pressure[8], double acouSrc[8], int nElems);
    
    
    void getNumNodesAndDim(int type, int& numNodes, int& dim);
    
    void  calccombustionsrctij_( double coordMat_fortran[24], double Nodal_Vel_fortran[24],
                                 double Nodal_Rho_fortran[8], double elemVec[8], 
                                 int elemType_fortran[1]);
    
    void  calccombustionsrcvector_( double coordMat_fortran[24], double Nodal_Vec_fortran[24],
                                    double elemVec[8], int elemType_fortran[1]);
    
    void  calccombustionsrcscalar_( double coordMat_fortran[24], double Nodal_Scalar_fortran[8],
                                    double elemVec[8], int elemType_fortran[1]);
    
    void  calccombustionsrcvectorsurf_( double coordMat_fortran[12], double Nodal_Vec_fortran[12],
                                        double elemVec[8], int elemType_fortran[1]);
    
    void  calccombustionsrctijsurf_( double coordMat_fortran[12], double Nodal_Vel_fortran[12],
                                     double Nodal_Rho_fortran[8], double elemVec[8], 
                                     int elemType_fortran[1]);

 #if defined (__cplusplus)
}
#endif

#ifdef CPLREADER
void  calcacousrcterms3d_quick_(const std::vector<int>& Topology, const std::vector<double>& NodalCoords, 
                                const std::vector<double>& velocity, std::vector<double>& acouSrc, int nElems);
void  calcacousrcterms3d_withPress_quick_(const std::vector<int>& Topology, const std::vector<double>& NodalCoords, 
                                const std::vector<double>& pressure, std::vector<double>& acouSrc, int nElems);
#endif

#endif
