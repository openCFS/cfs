
#include "coeff.h"
#include <cstdlib>
#include <stdio.h>      /* printf */
#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>
#include <boost/tokenizer.hpp>
#include <fstream>
#include <string>
#include <iomanip>

using namespace std;
using namespace boost;
int ijk2n(int i, int j, int k, int nx, int ny, int nz, std::string order="C") {
  assert(nx > 0 && ny > 0 && nz >0);
  assert(order == "C" || order == "F");
  // maps from three indices to one index
  int idx = i+nx*j+nx*ny*k;
  if (order =="F")
    idx = i * ny * nz + j * nz + k;
  assert(idx >= 0 && idx < nx*ny*nz);
  return idx;
}

// helper function for finite differences;
// does a range check and return index 'res' and spacing 'dx' within range of [0,nx)
// e.g. idx = end-1 and offset = +1 -> return  res=nx-1 and dx=1 as array element end=(end-1) + 1 exceeds the array bounds
void get_index(int idx, int offset, int end, double spacing, int& res, double& dx){
  assert(0 <= idx && idx < end);
  int new_idx = idx + offset;
  if (new_idx >= end || new_idx < 0){
    res = idx;
    dx = 0;
    return;
  }
  assert(0 <= new_idx && new_idx < end);
  res = new_idx;
  dx = spacing;
  // dx contains same value as before
}


void tricubic_solve(vector<double> & a, vector<double> x) {
  for (int i=0;i<64;i++) {
    a[i]=0.;
    for (int j=0;j<64;j++) {
      a[i]+=A[i][j]*x[j];
    }
  }
}

void tricubic_coeff(vector<double> & a, vector<double> f, vector<double> dfdx, vector<double> dfdy, vector<double> dfdz, vector<double> d2fdxdy,
    vector<double> d2fdxdz, vector<double> d2fdydz, vector<double> d3fdxdydz) {
  vector<double> x(64,0.);
  // Solve a linear system in order to get the interpolation coefficients
  for (int i=0;i<8;i++) {
    x[i]=f[i];
    x[8+i]=dfdx[i];
    x[16+i]=dfdy[i];
    x[24+i]=dfdz[i];
    x[32+i]=d2fdxdy[i];
    x[40+i]=d2fdxdz[i];
    x[48+i]=d2fdydz[i];
    x[56+i]=d3fdxdydz[i];
  }
  tricubic_solve(a,x);
}

void tricubic_partialderiv(vector<vector<vector<double> > > & dEda,vector<vector<vector<double> > > & dEdb,vector<vector<vector<double> > > & dEdc,
    vector<vector<vector<double> > > & dEdadb,vector<vector<vector<double> > > & dEdadc,vector<vector<vector<double> > > & dEdbdc,
    vector<vector<vector<double> > > & dEdadbdc,vector<double>a,vector<double>b,vector<double>c,vector<vector<vector<double> > >E){
  //Approximation of the derivatives with finite differences at the grid points

  // number of samples for each design parameter
  int m = a.size();
  int n = b.size();
  int o = c.size();

  double da = 1./(m-1);
  double db = 1./(n-1);
  double dc = 1./(o-1);

  //  std::cout << "tricubic_partialderiv: m,n,o: " << m << "," << n << "," << o << std::endl;
  // ip: i+1, im: i-1, jp: j+1, ...
  int ip, im, jp, jm, kp, km;
  // d_ip: da (i+1), d_im: da i(-1), ...
  double d_ip, d_im, d_jp, d_jm, d_kp, d_km;

  for (int i = 0; i < m;i++) {
    for (int j = 0; j < n;j++) {
      for (int k = 0; k < o;k++) {

        get_index(i,+1,m,da,ip,d_ip);
        get_index(i,-1,m,da,im,d_im);

        get_index(j,+1,n,db,jp,d_jp);
        get_index(j,-1,n,db,jm,d_jm);

        get_index(k,+1,o,dc,kp,d_kp);
        get_index(k,-1,o,dc,km,d_km);

        //f(i+1,j,k)-f(i-1,j,k) / 2
        dEda[i][j][k] = (E[ip][j][k] - E[im][j][k]) /(d_ip + d_im);
        dEdb[i][j][k] = (E[i][jp][k] - E[i][jm][k]) /(d_jp + d_jm);
        dEdc[i][j][k] = (E[i][j][kp] - E[i][j][km]) /(d_kp + d_km);
        //f(i+1,j+1,k)-f(i+1,j-1,k)-f(i-1,j+1,k)+f(i-1,j-1,k))/4
        dEdadb[i][j][k] = (E[ip][jp][k]-E[ip][jm][k] - E[im][jp][k]+E[im][jm][k])/((d_ip+d_im) * (d_jp+d_jm));
        dEdadc[i][j][k] = (E[ip][j][kp]-E[ip][j][km] - E[im][j][kp]+E[im][j][km])/((d_ip+d_im) * (d_kp+d_km));
        dEdbdc[i][j][k] = (E[i][jp][kp]-E[i][jp][km] - E[i][jm][kp]+E[i][jm][km])/((d_jp+d_jm) * (d_kp+d_km));
        //(f(i+1,j+1,z+1) - f(i+1,j,k-1) - f(i+1,j-1,k+1) + f(i+1,j-1,k-1) - f(i-1,j+1,k+1) + f(i-1,j+1,k-1)
        //+ f(i-1,j-1,k+1) - f(i-1,j-1,k-1))/8
        dEdadbdc[i][j][k] = (E[ip][jp][kp]-E[ip][jp][km]-E[ip][jm][kp]+E[ip][jm][km]-E[im][jp][kp]+
            E[im][jp][km]-E[im][jm][kp]+E[im][jm][km])/((d_ip+d_im) * (d_jp+d_jm) * (d_kp+d_km));
      }
    }
  }

//  /*********************cube faces*********************/
//
//  // face where s1 = 0 or s1 = end
//  int begin = 0;
//  int end = m-1;
//  for (int j=0;j<n;j++) {
//    for (int k=0;k<o;k++) {
//      //boundary treatment for partial derivatives dEda
//      dEda[begin][j][k] = (E[begin+1][j][k] - E[begin][j][k]) /(da);
//      dEda[end][j][k] = (E[end][j][k] - E[end-1][j][k]) /(da);
//      //      std::cout << "i,j,k:" << end << "," << j << "," << k << "  (" << E[end][j][k] << "-" << E[end-1][j][k] << ")/" << (da) << std::endl;
//      if (1 <= j && j < n-1) {
//        dEdb[begin][j][k] = (E[begin][j+1][k] - E[begin][j-1][k]) / (2*db);
//        dEdb[end][j][k] = (E[end][j+1][k]- E[end][j-1][k]) / (2*db);
////        std::cout << "dEdb:  i,j,k:" << begin << "," << j << "," << k << "  (" << E[begin][j+1][k] << "-" << E[begin][j-1][k] << ")/" << (2*db) << "=" << dEdb[begin][j][k] << std::endl;
//        dEdadb[begin][j][k] = (E[begin+1][j+1][k]-E[begin+1][j-1][k] - E[begin][j+1][k]+E[begin][j-1][k])/(2.*db * da);
//        dEdadb[end][j][k] = (E[end][j+1][k]-E[end][j-1][k] - E[end-1][j+1][k]+E[end-1][j-1][k])/(2.*db * da);
//      }
//      else if (j == 0){
//        dEdb[begin][j][k] = (E[begin][j+1][k] - E[begin][j][k]) / (db);
//        dEdb[end][j][k] = (E[end][j+1][k]- E[end][j][k]) / (db);
//        dEdadb[begin][j][k] = (E[begin+1][j+1][k]-E[begin+1][j][k] - E[begin][j+1][k]+E[begin][j][k])/(db * da);
//        dEdadb[end][j][k] = (E[end][j+1][k]-E[end][j][k] - E[end-1][j+1][k]+E[end-1][j][k])/(db * da);
//      }
//      else if (j == n-1){
//        dEdb[begin][j][k] = (E[begin][j][k] - E[begin][j-1][k]) / (db);
//        dEdb[end][j][k] = (E[end][j][k]- E[end][j-1][k]) / (db);
//        dEdadb[begin][j][k] = (E[begin+1][j][k]-E[begin+1][j-1][k] - E[begin][j][k]+E[begin][j-1][k])/(db * da);
//        dEdadb[end][j][k] = (E[end][j][k]-E[end][j-1][k] - E[end-1][j][k]+E[end-1][j-1][k])/(db * da);
//      }
//      if (1 <= k && k < o-1) {
//        dEdc[begin][j][k] = (E[begin][j][k+1] - E[begin][j][k-1]) / (2*dc);
//        dEdc[end][j][k] = (E[end][j][k+1] - E[end][j][k-1]) / (2*dc);
//        std::cout << "dEdc:  i,j,k:" << begin << "," << j << "," << k << "  (" << E[begin][j][k+1] << "-" << E[begin][j][k-1] << ")/" << (2*dc) << "=" << dEdc[begin][j][k] << std::endl;
//
//        dEdadc[begin][j][k] = (E[begin+1][j][k+1]-E[begin+1][j][k-1] - E[begin][j][k+1]+E[begin][j][k-1])/(2.*dc * da);
//        dEdadc[end][j][k] = (E[end][j][k+1]-E[end][j][k-1] - E[end-1][j][k+1]+E[end-1][j][k-1])/(2.*dc * da);
//        if (1 <= j && j < n-1){
//          dEdbdc[begin][j][k] = (E[begin][j+1][k+1]-E[begin][j+1][k-1] - E[begin][j-1][k+1]+E[begin][j-1][k-1])/(4. * db * dc);
//          dEdbdc[end][j][k] = (E[end][j+1][k+1]-E[begin][j+1][k-1] - E[end][j-1][k+1]+E[end][j-1][k-1])/(4. * db * dc);
//          dEdadbdc[begin][j][k] = (E[begin+1][j+1][k+1]-E[begin+1][j+1][k-1]-E[begin+1][j-1][k+1]+E[begin+1][j-1][k-1]-E[begin][j+1][k+1]+
//              E[begin][j+1][k-1]-E[begin][j-1][k+1]+E[begin][j-1][k-1])/(4.*db*dc * da);
//          dEdadbdc[end][j][k] = (E[end][j+1][k+1]-E[end][j+1][k-1]-E[end][j-1][k+1]+E[end][j-1][k-1]-E[end-1][j+1][k+1]+
//              E[end-1][j+1][k-1]-E[end-1][j-1][k+1]+E[end-1][j-1][k-1])/(4.*db*dc * da);
//        }
//        else if (j == 0){
//          dEdbdc[begin][j][k] = (E[begin][j+1][k+1]-E[begin][j+1][k-1] - E[begin][j][k+1]+E[begin][j][k-1])/(2. * db * dc);
//          dEdbdc[end][j][k] = (E[end][j+1][k+1]-E[begin][j+1][k-1] - E[end][j][k+1]+E[end][j][k-1])/(2. * db * dc);
//          dEdadbdc[begin][j][k] = (E[begin+1][j+1][k+1]-E[begin+1][j+1][k-1]-E[begin+1][j][k+1]+E[begin+1][j][k-1]-E[begin][j+1][k+1]+
//              E[begin][j+1][k-1]-E[begin][j][k+1]+E[begin][j][k-1])/(2.*dc * db * da);
//          dEdadbdc[end][j][k] = (E[end][j+1][k+1]-E[end][j+1][k-1]-E[end][j][k+1]+E[end][j][k-1]-E[end-1][j+1][k+1]+
//              E[end-1][j+1][k-1]-E[end-1][j][k+1]+E[end-1][j][k-1])/(2.*dc *db * da);
//        }
//        else if (j == n-1){
//          dEdbdc[begin][j][k] = (E[begin][j][k+1]-E[begin][j][k-1] - E[begin][j-1][k+1]+E[begin][j-1][k-1])/(2.*dc * db);
//          dEdbdc[end][j][k] = (E[end][j][k+1]-E[begin][j][k-1] - E[end][j-1][k+1]+E[end][j-1][k-1])/(2*dc * db);
//          dEdadbdc[begin][j][k] = (E[begin+1][j][k+1]-E[begin+1][j][k-1]-E[begin+1][j-1][k+1]+E[begin+1][j-1][k-1]-E[begin][j][k+1]+
//              E[begin][j][k-1]-E[begin][j-1][k+1]+E[begin][j-1][k-1])/(2.*dc * da * db);
//          dEdadbdc[end][j][k] = (E[end][j][k+1]-E[end][j][k-1]-E[end][j-1][k+1]+E[end][j-1][k-1]-E[end-1][j][k+1]+
//              E[end-1][j][k-1]-E[end-1][j-1][k+1]+E[end-1][j-1][k-1])/(2.*dc * da * db);
//        }
//      }
//      if (k == 0) {
//        dEdc[begin][j][k] = (E[begin][j][k+1] - E[begin][j][k]) / (dc);
//        dEdc[end][j][k] = (E[end][j][k+1] - E[end][j][k]) / (dc);
//
//      }
//    }
//  }
//
//  // face where s2 = 0 or s2 = end
//  begin = 0;
//  end = n-1;
//  for (int i=0;i<m;i++) {
//    for (int k=0;k<o;k++) {
//      if (1 <= i && i < m-1) {
//        dEda[i][begin][k] = (E[i+1][begin][k] - E[i-1][begin][k]) / (2*da);
//        dEda[i][end][k] = (E[i+1][end][k] - E[i-1][end][k]) / (2*da);
//        dEdadb[i][begin][k] = (E[i+1][begin+1][k]-E[i+1][begin][k] - E[i-1][begin+1][k]+E[i-1][begin][k])/(2.*da * db);
//        dEdadb[i][end][k] = (E[i+1][end][k]-E[i+1][end-1][k] - E[i-1][end][k]+E[i-1][end-1][k])/(2.*da * db);
//      }
//      if (1 <= k && k < o-1) {
//        dEdc[i][begin][k] = (E[i][begin][k+1] - E[i][begin][k-1]) / (2*dc);
//        dEdc[i][end][k] = (E[i][end][k+1] - E[i][end][k-1]) / (2*dc);
//        dEdbdc[i][begin][k] = (E[i][begin+1][k+1]-E[i][begin+1][k-1] - E[i][begin][k+1]+E[i][begin][k-1])/(2.*dc * db);
//        dEdbdc[i][end][k] = (E[i][end][k+1]-E[i][end][k-1] - E[i][end-1][k+1]+E[i][end-1][k-1])/(2.*dc * db);
//        if (1 <= i && i < m-1) {
//          dEdadc[i][begin][k] = (E[i+1][begin][k+1]-E[i+1][begin][k-1] - E[i-1][begin][k+1]+E[i-1][begin][k-1])/(4. * da * dc);
//          dEdadc[i][end][k] = (E[i+1][end][k+1]-E[i+1][end][k-1] - E[i-1][end][k+1]+E[i-1][end][k-1])/(4. * da * dc);
//          dEdadbdc[i][begin][k] = (E[i+1][begin+1][k+1]-E[i+1][begin+1][k-1]-E[i+1][begin][k+1]+E[i+1][begin][k-1]-E[i-1][begin+1][k+1]+
//              E[i-1][begin+1][k-1]-E[i-1][begin][k+1]+E[i-1][begin][k-1])/(4.*da*dc * db);
//          dEdadbdc[i][end][k] = (E[i+1][end][k+1]-E[i+1][end][k-1]-E[i+1][end-1][k+1]+E[i+1][end-1][k-1]-E[i-1][end][k+1]+
//              E[i-1][end][k-1]-E[i-1][end-1][k+1]+E[i-1][end-1][k-1])/(4.*da*dc * db);
//        }
//      }
//      //boundary treatment for partial derivatives dEdb
//      dEdb[i][begin][k] = (E[i][begin+1][k] - E[i][begin][k]) /(db);
//      dEdb[i][end][k] = (E[i][end][k] - E[i][end-1][k]) /(db);
//
//    }
//  }
//
//  // face where s3 = 0 or s3 = end
//  begin = 0;
//  end = o-1;
//  for (int i=0;i<m;i++) {
//    for (int j=0;j<n;j++) {
//      if (1 <= i && i < m-1) {
//        dEda[i][j][begin] = (E[i+1][j][begin] - E[i-1][j][begin]) / (2*da);
////        if (i == 1 && j == 0) {
////          std::cout <<"i=1 j=0 k=0:" << dEda[i][j][begin] << "=(" << E[i+1][j][begin] << "-" << E[i-1][j][begin] <<") / (2*" << da << ")" <<std::endl;
////        }
//        dEda[i][j][end] = (E[i+1][j][end] - E[i-1][j][end]) / (2*da);
//        dEdadc[i][j][begin] = (E[i+1][j][begin+1]-E[i+1][j][begin] - E[i-1][j][begin+1]+E[i-1][j][begin])/(2.*da * dc);
//        dEdadc[i][j][end] = (E[i+1][j][end]-E[i+1][j][end-1] - E[i-1][j][end]+E[i-1][j][end-1])/(2.*da * dc);
//      }
//      if (1 <= j && j < n-1) {
//        dEdb[i][j][begin] = (E[i][j+1][begin] - E[i][j-1][begin]) / (2*db);
//        dEdb[i][j][end] = (E[i][j+1][end] - E[i][j-1][end]) / (2*db);
//        dEdbdc[i][j][begin] = (E[i][j+1][begin+1]-E[i][j+1][begin] - E[i][j-1][begin+1]+E[i][j-1][begin])/(2.*db * dc);
//        dEdbdc[i][j][end] = (E[i][j+1][end]-E[i][j+1][end-1] - E[i][j-1][end]+E[i][j-1][end-1])/(2*db * dc);
//        if (1 <= i && i < m-1) {
//          dEdadb[i][j][begin] = (E[i+1][j+1][begin]-E[i+1][j-1][begin] - E[i-1][j+1][begin]+E[i-1][j-1][begin])/(4. * da * db);
//          dEdadb[i][j][end] = (E[i+1][j+1][end]-E[i+1][j-1][end] - E[i-1][j+1][end]+E[i-1][j-1][end])/(4. * da * db);
//          dEdadbdc[i][j][begin] = (E[i+1][j+1][begin+1]-E[i+1][j+1][begin]-E[i+1][j-1][begin+1]+E[i+1][j-1][begin]-E[i-1][j+1][begin+1]+
//              E[i-1][j+1][begin]-E[i-1][j-1][begin+1]+E[i-1][j-1][begin])/(4.*da*db * dc);
//          dEdadbdc[i][j][end] = (E[i+1][j+1][end]-E[i+1][j+1][end-1]-E[i+1][j-1][end]+E[i+1][j-1][end-1]-E[i-1][j+1][end]+
//              E[i-1][j+1][end-1]-E[i-1][j-1][end]+E[i-1][j-1][end-1])/(4.*da*db * dc);
//        }
//      }
//      //boundary treatment for partial derivatives dEdc
//      dEdc[i][j][begin] = (E[i][j][begin+1] - E[i][j][begin]) /(dc);
//      //      dEdc[i][j][end-1] = (E[i][j][end-1] - E[i][j][end-2]) /(dc);
//      dEdc[i][j][end] = (E[i][j][end] - E[i][j][end-1]) /(dc);
//    }
//  }
//
//  std::cout << "dEdadb[0][0][1]:" << dEdadb[0][0][1]<< std::endl;

  /*********************cube edges*********************/
  //  int i,j,k;
  //  i = 1;
  //  j = 0;
  //  k = 0;
  //  for (;i < m-1; i++) {
  //    dEda[i][j][k] =(E[i+1][j][k] - E[i-1][j][k]) / (2*da);
  //    dEdb[i][j][k] =(E[i][j+1][k] - E[i][j][k]) / (db);
  //    dEdc[i][j][k] =(E[i][j][k+1] - E[i][j][k]) / (dc);
  //    // dEdadb[i][j][k] = (E[i+1][j+1][k]-E[i+1][j-1][k] - E[i-1][j+1][k]+E[i-1][j-1][k])/(4. * da * db);
  //    dEdadb[i][j][k] = (E[i+1][j+1][k]-E[i+1][j][k] - E[i-1][j+1][k]+E[i-1][j][k])/(2.*db * da);
  //    // dEdadc[i][j][k] = (E[i+1][j][k+1]-E[i+1][j][k-1] - E[i-1][j][k+1]+E[i-1][j][k-1])/(4. * da * dc);
  //    dEdadc[i][j][k] = (E[i+1][j][k+1]-E[i+1][j][k-1] - E[i][j][k+1]+E[i][j][k-1])/(2.*dc * da);
  ////    std::cout << "(" << E[i+1][j][k] << "-" << E[i-1][j][k] << ")/" << (2*da) << std::endl;
  //  }
  //
  //  i = 0;
  //  j = 1;
  //  k = 0;
  //  for (;j < n-1; j++){
  //    dEda[i][j][k] =(E[i+1][j][k] - E[i][j][k]) / (da);
  //    dEdb[i][j][k] = (E[i][j+1][k] - E[i][j-1][k]) / (2*db);
  //    dEdc[i][j][k] =(E[i][j][k+1] - E[i][j][k]) / (dc);
  ////    std::cout << "i,j,k:0," << j << ",0  (" << E[i][j+1][k] << "-" << E[i][j-1][k] << ")/" << (2*db) << std::endl;
  //  }
  //
  //  i = 0;
  //  j = 0;
  //  k = 1;
  //  for (;k < o-1; k++){
  //    dEda[i][j][k] =(E[i+1][j][k] - E[i][j][k]) / (da);
  //    dEdb[i][j][k] =(E[i][j+1][k] - E[i][j][k]) / (db);
  //    dEdc[i][j][k] =(E[i][j][k+1] - E[i][j][k-1]) / (2*dc);
  //  }

  //  for (int i = 0; i < m; i++)
  //    std::cout << "i,j,k:" << 0 << "," << i << ",0" << "  dEdb:" << dEdb[0][i][0] << std::endl;
  //
  //  std::cout << std::endl;
  //  for (int i = 0; i < m; i++)
  //    std::cout << "i,j,k:1," << i << ",0" << "  dEdb:" << dEdb[1][i][0] << std::endl;
  //
  //  std::cout << std::endl;
  //  for (int i = 0; i < m; i++)
  //    std::cout << "i,j,k:2," << i << ",0" << "  dEdb:" << dEdb[2][i][0] << std::endl;
  //
  //  std::cout << std::endl;
  //  for (int i = 0; i < m; i++)
  //    std::cout << "i,j,k:3," << i << ",0" << "  dEdb:" << dEdb[3][i][0] << std::endl;
  // missing: boundary treatment for mixed partial derivatives, they are set to zero
}

void tricubic_offline(vector<vector<double> > & Coeff,vector<double> a, vector<double> b, vector<double> c, vector<vector<vector<double> > > & E, int m, int n, int o, double dx,double dy,double dz, vector<vector<double> > & deriv) {
  // Calculation of the coefficients of the interpolation polynomial for all possible intervals
  vector<double> res(64);
  vector<vector<vector <double > > > dEda(m,vector<vector<double> >(n,vector<double>(o,0.))),dEdb(m,vector<vector<double> >(n,vector<double>(o,0.))),
      dEdc(m,vector<vector<double> >(n,vector<double>(o,0.))),dEdadb(m,vector<vector<double> >(n,vector<double>(o,0.))),
      dEdadc(m,vector<vector<double> >(n,vector<double>(o,0.))),dEdbdc(m,vector<vector<double> >(n,vector<double>(o,0.))),
      dEdadbdc(m,vector<vector<double> >(n,vector<double>(o,0.)));
  tricubic_partialderiv(dEda,dEdb,dEdc,dEdadb,dEdadc,dEdbdc,dEdadbdc,a,b,c,E);
  vector<double> f(8,0);
  vector<double> Eda(8,0);
  vector<double> Edb(8,0);
  vector<double> Edc(8,0);
  vector<double> Edadb(8,0);
  vector<double> Edadc(8,0);
  vector<double> Edbdc(8,0);
  vector<double> Edadbdc(8,0);
  // for each interval, we store the data (function values and partial derivatives) needed to compute the interpolation polynomial
  assert(deriv.size() == (m-1)*(n-1)*(o-1));
//  deriv.resize((m-1)*(n-1)*(o-1));
//  for (int i = 0; i < (m-1)*(n-1)*(o-1); i++)
//    deriv[i].resize(64);

  for (int i=0;i<m-1;i++) {
    for (int j=0;j<n-1;j++) {
      for (int k =0;k<o-1;k++) {
        f[0]= E[i][j][k];
        f[1]= E[i+1][j][k];
        f[2]= E[i][j+1][k];
        f[3]= E[i+1][j+1][k];

//        if ((i == 0 && j == 0 && k== 1)){// || (i == 1 && j == 0 && k== 0)){
//          std::cout << "i,j,k:" << i << "," << j << "," << k << std::endl;
//          std::cout << "f[0]: " << f[0]<< std::endl;
//          std::cout << "f[1]: " << f[1] << std::endl;
//          std::cout << "f[2]: " << f[2] << std::endl;
//          std::cout << "f[3]: " << f[3] << std::endl;
//          std::cout << "f[4]: " << f[4] << std::endl;
//          std::cout << "f[5]: " << f[5] << std::endl;
//          std::cout << "f[6]: " << f[6] << std::endl;
//          std::cout << "f[7]: " << f[7] << std::endl;
//        }

        f[4]= E[i][j][k+1];
        f[5]= E[i+1][j][k+1];
        f[6]= E[i][j+1][k+1];
        f[7]= E[i+1][j+1][k+1];

        /* compute all the partial derivatives */
        //dfdx(i,j,k)= f(i+1,j,k)-f(i-1,j,k) / 2
        // 8 NN pixels
        Eda[0]= dEda[i][j][k];
        Eda[1]= dEda[i+1][j][k];
        Eda[2]= dEda[i][j+1][k];
        Eda[3]= dEda[i+1][j+1][k];

        Eda[4]= dEda[i][j][k+1];
        Eda[5]= dEda[i+1][j][k+1];
        Eda[6]= dEda[i][j+1][k+1];
        Eda[7]= dEda[i+1][j+1][k+1];

//        if ((i == 0 && j == 0 && k== 1)){// || (i == 1 && j == 0 && k== 0)){
//          std::cout << "dEda[i][j][k]: " << Eda[0] << std::endl;
//          std::cout << "dEda[i+1][j][k]: " << Eda[1] << std::endl;
//          std::cout << "dEda[i][j+1][k]: " << Eda[2] << std::endl;
//          std::cout << "dEda[i+1][j+1][k]: " << Eda[3] << std::endl;
//
//          std::cout << "dEda[i][j][k+1]: " << Eda[4] << std::endl;
//          std::cout << "dEda[i+1][j][k+1]: " << Eda[5] << std::endl;
//          std::cout << "dEda[i][j+1][k+1]: " << Eda[6] << std::endl;
//          std::cout << "dEda[i+1][j+1][k+1]: " << Eda[7] << std::endl;
//        }

        Edb[0]= dEdb[i][j][k];
        Edb[1]= dEdb[i+1][j][k];
        Edb[2]= dEdb[i][j+1][k];
        Edb[3]= dEdb[i+1][j+1][k];

        Edb[4]= dEdb[i][j][k+1];
        Edb[5]= dEdb[i+1][j][k+1];
        Edb[6]= dEdb[i][j+1][k+1];
        Edb[7]= dEdb[i+1][j+1][k+1];


        Edc[0]= dEdc[i][j][k];
        Edc[1]= dEdc[i+1][j][k];
        Edc[2]= dEdc[i][j+1][k];
        Edc[3]= dEdc[i+1][j+1][k];

        Edc[4]= dEdc[i][j][k+1];
        Edc[5]= dEdc[i+1][j][k+1];
        Edc[6]= dEdc[i][j+1][k+1];
        Edc[7]= dEdc[i+1][j+1][k+1];

//        std::cout << "i,j,k:" << i << "," << j << "," << k << " dEdc[i][j][k]:" << dEdc[i][j][k] << std::endl;

        //d2fdxdy(i,j,k)=(f(i+1,j+1,k)-f(i+1,j,k)-f(i,j+1,k)+f(i-1,j-1,k))/4
        Edadb[0]= dEdadb[i][j][k];
        Edadb[1]= dEdadb[i+1][j][k];
        Edadb[2]= dEdadb[i][j+1][k];
        Edadb[3]= dEdadb[i+1][j+1][k];

        Edadb[4]= dEdadb[i][j][k+1];
        Edadb[5]= dEdadb[i+1][j][k+1];
        Edadb[6]= dEdadb[i][j+1][k+1];
        Edadb[7]= dEdadb[i+1][j+1][k+1];


        //d2fdxdz(i,j,k)=(f(i+1,j,k+1)-f(i+1,j,k)-f(i,j,k+1)+f(i-1,j,k-1))/4
        Edadc[0]= dEdadc[i][j][k];
        Edadc[1]= dEdadc[i+1][j][k];
        Edadc[2]= dEdadc[i][j+1][k];
        Edadc[3]= dEdadc[i+1][j+1][k];

        Edadc[4]= dEdadc[i][j][k+1];
        Edadc[5]= dEdadc[i+1][j][k+1];
        Edadc[6]= dEdadc[i][j+1][k+1];
        Edadc[7]= dEdadc[i+1][j+1][k+1];


        //d2fdxdy(i,j,k)=(f(i,j+1,k+1)-f(i,j+1,k)-f(i,j,k+1)+f(i,j-1,k-1))/4
        Edbdc[0]= dEdbdc[i][j][k];
        Edbdc[1]= dEdbdc[i+1][j][k];
        Edbdc[2]= dEdbdc[i][j+1][k];
        Edbdc[3]= dEdbdc[i+1][j+1][k];

        Edbdc[4]= dEdbdc[i][j][k+1];
        Edbdc[5]= dEdbdc[i+1][j][k+1];
        Edbdc[6]= dEdbdc[i][j+1][k+1];
        Edbdc[7]= dEdbdc[i+1][j+1][k+1];


        //d3fdxdydz(i,j,k)= f(i+1,j+1,z+1) - f(i+1,j,k+1) - f(i,j+1,k+1) + f(i-1,j-1,k+1) - f(i+1,j+1,k-1) + f(i+1,j,k-1) + f(i,j+1,k-1) - f(i-1,j-1,k-1)
        Edadbdc[0]= dEdadbdc[i][j][k];
        Edadbdc[1]= dEdadbdc[i+1][j][k];
        Edadbdc[2]= dEdadbdc[i][j+1][k];
        Edadbdc[3]= dEdadbdc[i+1][j+1][k];

        Edadbdc[4]= dEdadbdc[i][j][k+1];
        Edadbdc[5]= dEdadbdc[i+1][j][k+1];
        Edadbdc[6]= dEdadbdc[i][j+1][k+1];
        Edadbdc[7]= dEdadbdc[i+1][j+1][k+1];

        int idx = ijk2n(i,j,k,m-1,n-1,o-1,"F");
        vector<double>& d = deriv[idx];
//        if ((i == 0 && j == 0 && k== 1))
//          std::cout << "deriv idx:" << idx << std::endl;
        // store info on interpolation data
        // clear d as it might be already initialized with 0 entries
        d.clear();
        d.reserve(64);
        d.insert(d.begin(),f.begin(),f.end());
        d.insert(d.end(),Eda.begin(),Eda.end());
        d.insert(d.end(),Edb.begin(),Edb.end());
        d.insert(d.end(),Edc.begin(),Edc.end());
        d.insert(d.end(),Edadb.begin(),Edadb.end());
        d.insert(d.end(),Edadc.begin(),Edadc.end());
        d.insert(d.end(),Edbdc.begin(),Edbdc.end());
        d.insert(d.end(),Edadbdc.begin(),Edadbdc.end());
        assert(d.size() == 64);
        //scale values
        double eps = 1e-6;
        for (int q=0;q<8;q++) {
//          if (i > 0)
//            assert(Eda[q] > eps);
//          if (j > 0){
//            std::cout << "i,j,k:" << i << "," << j << "," << k << " q:" << q << " Edb[q]:" << Edb[q] << std::endl;
//            assert(Edb[q] > eps);
//          }
////          if (k > 0) {
////            std::cout << "i,j,k:" << i << "," << j << "," << k << " q:" << q << " Edc[q]:" << Edc[q] << std::endl;
////            assert(Edc[q] > eps);
////          }
//          if (i > 0 && j > 0)
//            assert(Edadb[q] > eps);
//          if (i > 0 && k > 0)
//            assert(Edadc[q] > eps);
//          if (j > 0 && k > 0)
//            assert(Edbdc[q] > eps);
//          if (i > 0 && j > 0 && k > 0)
//            assert(Edadbdc[q] > eps);
          f[q]*=1.0;
          Eda[q]*=dx;
          Edb[q]*=dy;
          Edc[q]*=dz;
          Edadb[q]*=dx*dy;
          Edadc[q]*=dx*dz;
          Edbdc[q]*=dy*dz;
          Edadbdc[q]*=dx*dy*dz;
        }
        vector<double> res(64,0.);
        tricubic_coeff(res,f,Eda,Edb,Edc,Edadb,Edadc,Edbdc,Edadbdc);
        for (int ii = 0;ii<64;ii++) {
          Coeff[idx][ii] = res[ii];
          //          Coeff[(n-1)*(o-1)*i+(o-1)*j+k][ii] = res[ii];
        }
      }
    }
  }


}

double tricubic_eval(vector<double> a, double x, double y, double z) {
  //Evaluation of the interpolation polynomial at point x,y,z
  double ret= 0.;
  for (int i=0;i<4;i++) {
    for (int j=0;j<4;j++) {
      for (int k=0;k<4;k++) {
//        if (j == 1 && k == 0) {
//          std::cout << "i,j,k:" << i << "," << j << "," << k << "  a_ijk:" << a[ijk2n(i,j,k,4,4,4)] << std::endl;
//          std::cout << "res=" << ret << "+" << a[ijk2n(i,j,k,4,4,4)] << "*" << pow(x,i) << "*" << pow(y,j) << "*" << pow(z,k)  << std::endl;
//        }
        ret+=a[ijk2n(i,j,k,4,4,4)]*pow(x,i)*pow(y,j)*pow(z,k);
      }
    }
  }
  return ret;
}

void write_to_xml(fstream & f, vector<vector<vector<double> > > Coeff, vector<double> aa, vector<double> bb, vector<double> cc, int mode, vector<vector<vector<double> > > Coeff_deriv) {
  f<<"<homRectC1 notation=\"voigt\">"<<endl;
  f<<"<a>"<<endl;
  f<<"<matrix dim1=\""<<aa.size()<<"\" dim2=\"1\">"<<endl;
  f<<"<real>"<<endl;
  for (int i = 0;i<aa.size();i++) {
    f<<aa[i]<<" ";
  }
  f<<endl;
  f<<"</real>"<<endl;
  f<<"</matrix>"<<endl;
  f<<"</a>"<<endl;

  f<<"<b>\n<matrix dim1=\""<<bb.size()<<"\" dim2=\"1\">\n<real>\n";
  for (int i = 0;i<bb.size();i++) {
    f<<bb[i]<<" ";
  }
  f<<"\n</real>\n</matrix>\n</b>\n";
  f<<"<c> \n<matrix dim1=\""<<cc.size()<<"\" dim2=\"1\">\n<real>\n";
  for (int i = 0;i<cc.size();i++) {
    f<<cc[i]<<" ";
  }
  f<<"\n</real>\n</matrix>\n</c>\n";
  // number of intervals
  int n_int = (aa.size()-1)*(bb.size()-1)*(cc.size()-1);

  std::string name_coeff;
  for (int cidx = 0; cidx < Coeff.size(); cidx++) {
    switch(cidx){
      case 0:
        name_coeff = "coeff11";
        break;
      case 1:
        name_coeff = "coeff12";
        break;
      case 2:
        name_coeff = "coeff13";
        break;
      case 3:
        name_coeff = "coeff22";
        break;
      case 4:
        name_coeff = "coeff23";
        break;
      case 5:
        name_coeff = "coeff33";
        break;
      case 6:
        name_coeff = "coeff44";
        break;
      case 7:
        name_coeff = "coeff55";
        break;
      case 8:
        name_coeff = "coeff66";
        break;
      default:
        assert(false);
    }
    std::string name_deriv = name_coeff + "_deriv";
    // interpolation coefficients
    f << "<" << name_coeff << ">\n<matrix dim1=\""<<n_int<<"\" dim2=\"64\">\n<real>\n";
    vector<vector<double>>& data = Coeff[cidx];
    for (int i=0;i<n_int;i++) {
      for (int j=0;j<64;j++) {
        f<<setprecision(12)<<data[i][j]<<" ";
      }
      f<<endl;
    }
    f<<"</real>\n</matrix>\n</" << name_coeff << ">\n";

    // interpolation data (derivatives)
    f << "<" << name_deriv << ">\n<matrix dim1=\""<<n_int<<"\" dim2=\"64\">\n<real>\n";
    vector<vector<double>>& data_deriv = Coeff_deriv[cidx];
    for (int i=0;i<n_int;i++) {
      for (int j=0;j<64;j++) {
        f<<setprecision(12)<<data_deriv[i][j]<<" ";
      }
      f<<endl;
    }
    f<<"</real>\n</matrix>\n</" << name_deriv << ">\n";

  }
}

void write_to_xml_vol(fstream & f, vector<vector<double>  > Coeff, vector<vector<double>  > vol_deriv, vector<double> aa, vector<double> bb, vector<double> cc, bool combined = true) {
  if (!combined) {
    f<<"<vol3D>"<<endl;
    f<<"<a>"<<endl;
    f<<"<matrix dim1=\""<<aa.size()<<"\" dim2=\"1\">"<<endl;
    f<<"<real>"<<endl;
    for (int i = 0;i<aa.size();i++) {
      f<<aa[i]<<" ";
    }
    f<<endl;
    f<<"</real>"<<endl;
    f<<"</matrix>"<<endl;
    f<<"</a>"<<endl;

    f<<"<b>"<<endl;
    f<<"<matrix dim1=\""<<bb.size()<<"\" dim2=\"1\">"<<endl;
    f<<"<real>"<<endl;
    for (int i = 0;i<bb.size();i++) {
      f<<bb[i]<<" ";
    }
    f<<"\n</real>"<<endl;
    f<<"</matrix>"<<endl;
    f<<"</b>"<<endl;
    f<<"<c> \n<matrix dim1=\""<<cc.size()<<"\" dim2=\"1\">\n<real>\n";
    for (int i = 0;i<cc.size();i++) {
      f<<cc[i]<<" ";
    }
    f<<endl;
    f<<"</real>"<<endl;
    f<<"</matrix>"<<endl;
    f<<"</c>"<<endl;
  }
  int n_int = (aa.size()-1)*(bb.size()-1)*(cc.size()-1);
  f<<"<volcoeff>"<<endl;
  f<<"<matrix dim1=\""<<n_int<<"\" dim2=\"64\">"<<endl;
  f<<"<real>"<<endl;
  for (int i=0;i<n_int;i++) {
    for (int j=0;j<64;j++) {
      f<<setprecision(12)<<Coeff[i][j]<<" ";
    }
    f<<endl;
  }
  f<<"</real>"<<endl;
  f<<"</matrix>"<<endl;
  f<<"</volcoeff>"<<endl;

  // write out the derivatives
  f<<"<volcoeff_deriv>"<<endl;
  f<<"<matrix dim1=\"" << n_int << "\" dim2=\"64\">" << endl;
  f<<"<real>"<<endl;
  for (int i=0;i<n_int;i++) {
    for (int j=0;j<64;j++) {
      f<<setprecision(12) << vol_deriv[i][j] << " ";
    }
    f<<endl;
  }
  f<<"</real>"<<endl;
  f<<"</matrix>"<<endl;
  f<<"</volcoeff_deriv>"<<endl;

  if (!combined)
    f<<"</vol3D>"<<endl;
  //f.close();
}


double Calc3DCrossVolume(double stiff1, double stiff2, double stiff3, bool derivative, int der) {
  double vol;
  if (!derivative) {
    if (stiff1 >= stiff2 && stiff1 >= stiff3) {
      vol = stiff1*stiff1 + stiff2*stiff2 + stiff3*stiff3 - stiff1*stiff3*stiff3 - stiff1*stiff2*stiff2;
    } else if (stiff2 >= stiff1 && stiff2 >= stiff3) {
      vol = stiff1*stiff1 + stiff2*stiff2 + stiff3*stiff3 - stiff2*stiff3*stiff3 - stiff2*stiff1*stiff1;
    } else { 
      // case: (stiff3 >= stiff1 && stiff3 >= stiff2)
      vol = stiff1*stiff1 + stiff2*stiff2 + stiff3*stiff3 - stiff3*stiff2*stiff2 - stiff3*stiff1*stiff1;
    }
    return vol;
  } else {
    switch(der)
    {
      case 1:
        if (stiff1 >= stiff2 && stiff1 >= stiff3) {
          vol = 2*stiff1 - stiff3*stiff3 - stiff2*stiff2;
        } else if (stiff2 >= stiff1 && stiff2 >= stiff3) {
          vol = 2*stiff1 - 2* stiff2*stiff1;
        } else {
          //case: (stiff3 >= stiff1 && stiff3 >= stiff2) {
          vol = 2*stiff1 - 2* stiff3*stiff1;
        }
        return vol;
      case 2:
        if (stiff1 >= stiff2 && stiff1 >= stiff3) {
          vol = 2*stiff2 - 2*stiff1*stiff2;
        } else if (stiff2 >= stiff1 && stiff2 >= stiff3) {
          vol = 2*stiff2 - stiff3*stiff3 - stiff1*stiff1;
        } else {
          //case: (stiff3 >= stiff1 && stiff3 >= stiff2) {
          vol = 2*stiff2 - 2*stiff3*stiff2;
        }
        return vol;
      case 3:
        if (stiff1 >= stiff2 && stiff1 >= stiff3) {
          vol = 2*stiff3 - 2*stiff1*stiff3;
        } else if (stiff2 >= stiff1 && stiff2 >= stiff3) {
          vol = 2*stiff3 - 2*stiff2*stiff3;
        } else {
          // (stiff3 >= stiff1 && stiff3 >= stiff2) {
          vol = 2*stiff3 - stiff2*stiff2 - stiff1*stiff1;
        }
        return vol;
      default:
        return 0.0;
    }
  }
}

std::string toString(std::vector<double> vec){
  std::string str = "";
  for(std::vector<double>::iterator it = vec.begin(); it < vec.end(); it++)
    str +=  " " + std::to_string(*it);

  return str;
}

int main(int argc, char * argv[]) {
  int mode;
  cout<<"Select number of interpolation goal: 3D cross shaped volume 1 or homogenized tensor 2 or debug mode 3 or homogenized tensor with volume 4 or mode 5 create plot data for volume."<<endl;
  cin>>mode;
  /*if (argc !=2 ) {
    cout <<"Error number of input parameter wrong"<<endl;
    return -1;
  }*/
  if (mode==1) {
    cout<<"Enter number of discretized intervals of interval [0,1] for volume grid"<<endl;
    double disc;
    cin>>disc;
    int m = disc;
    int n = disc;
    int o = disc;
    double da = 1./(m);
    double db = 1./(n);
    double dc = 1./(o);

    //    std::cout << "main() mode 1: m,n,o: " << m << "," << n << "," << o << std::endl;

    //Create vector with the different sizes
    vector<double> aa(m,0.);
    vector<double> bb(n,0.);
    vector<double> cc(o,0.);
    for (int i=0;i<m;i++) {
      aa[i]=(static_cast<double>(i)/static_cast<double>(m-1));
      bb[i]=(static_cast<double>(i)/static_cast<double>(n-1));
      cc[i]=(static_cast<double>(i)/static_cast<double>(o-1));
    }
    cout<<"data written"<<endl;
    //Calculate coefficient matrix for the interpolation of the different entries in E
    vector<vector<double> > Coeff((m-1)*(n-1)*(o-1),vector<double>(64,0.));
    vector<vector<vector <double > > > E(m,vector<vector<double> >(n,vector<double>(o,0.)));
    for (int i=0;i<m;i++) {
      for (int j=0;j<n;j++){
        for (int k=0;k<o;k++) {
          E[i][j][k] = Calc3DCrossVolume(aa[i], bb[j], cc[k], false,0.);
        }
      }
    }
    /*cout<<"a= [ ";
    for (int i=0;i<m;i++) {
      cout<<" "<<aa[i];
    }
    cout<<endl;
    for (int i=0;i<m;i++) {
        cout<<"E("<<i<<") = "<<endl;
        for (int j=0;j<n;j++) {
          for (int k=0;k<o;k++) {
            cout<<" "<<E[i][j][k];
          }
          cout<<endl;
        }
    }*/

    cout<<"Volume table calculated."<<endl;
    vector<vector<double>> vol_deriv((m-1)*(n-1)*(o-1),vector<double>(64,0.));
    tricubic_offline(Coeff,aa, bb, cc, E,m, n, o,da,db,dc,vol_deriv);
    cout<<"Insert output xml file name:"<<endl;
    string name;
    cin>>name;
    string file(name);
    fstream f;
    f.open(file.c_str(), ios::out);
    //Write interpolation coefficients in the xml file above
    write_to_xml_vol(f, Coeff, vol_deriv, aa, bb, cc, false);
    f.close();

  } else if (mode == 2 || mode == 3 || mode == 4 || mode == 5) {
    cout<<"Enter name of the homogenized data file:"<<endl;
    string input;
    string input_vol = "none";
    cin>>input;
    if (mode == 4 || mode == 5) {
      cout<<"Enter name of the homogenized volume data file:"<<endl;
      cin>>input_vol;
      cout<<"Read "<<input_vol<<endl;
    }
    double result;
    double x1 = -1;
    double x2 = -1;
    double x3 = -1;
    if (mode == 3) {
      cout<<"Evaluate interpolation polynomial at x1 x2 x3. Enter x1 x2 x3!"<<endl;
      cin >> x1;
      cin >> x2;
      cin >> x3;
    }
    cout<<"Read "<<input<<endl;
    string zeile;
    const char *inp = input.c_str();
    const char *inp_vol = input_vol.c_str();
    ifstream fin(inp);
    ifstream fin_vol;
    if (inp_vol != "none") {
      fin_vol.open(inp_vol);
      //Read header of the input file from homogenization volume
      getline(fin_vol, zeile, '\n');
    }
    //Read header of the input file from homogenization
    getline(fin, zeile, '\n');
    typedef tokenizer<char_separator<char> > tokenizer;
    tokenizer tok(zeile);
    tokenizer::iterator head = tok.begin();
    int m = atoi((*head).c_str())+1;
    int n = atoi((*(++head)).c_str())+1;
    int o = atoi((*(++head)).c_str())+1;

    //    std::cout << "m,n,o:" << m << "," << n << "," << o << std::endl;

    //Read data of the material catalogue into a data structure
    int ncol;
    if (mode == 2 || mode == 3 || mode == 4) {
      ncol = 3 + 9;
    } else {
      ncol = 2 + 6 + 1;
    }
    vector<vector<double> > data(m*n*o,vector<double>(ncol,0.));
    int count1 = 0;
    int count2 = 0;
    while(!(fin.eof())) {
      getline(fin, zeile, '\n');
      istringstream ss(zeile);
      count2=0;
      for (double d; ss >> d; ) {
        data[count1][count2] = d;
        count2++;
      }
      count1++;
    }
    vector<vector<double> > data_vol(m*n*o,vector<double>(4,0.));
    if (fin_vol.is_open()) {
      int count1 = 0;
      int count2 = 0;
      while(!(fin_vol.eof())) {
        getline(fin_vol, zeile, '\n');
        istringstream ss(zeile);
        count2=0;
        for (double d; ss >> d; ) {
          data_vol[count1][count2] = d;
          count2++;
        }
        count1++;
      }
    }
    cout<<"Files read."<<endl;
    //Create vector with the different sizes
    vector<double> aa(m,0.);
    vector<double> bb(n,0.);
    vector<double> cc(o,0.);
    for (int i=0;i<m;i++) {
      aa[i]=(static_cast<double>(i)/static_cast<double>(m-1));
      bb[i]=(static_cast<double>(i)/static_cast<double>(n-1));
      cc[i]=(static_cast<double>(i)/static_cast<double>(o-1));
    }

    double da = 1./(m-1);
    double db = 1./(n-1);
    double dc = 1./(o-1);

    //    std::cout << "main() mode 2,3,4,5: m,n,o: " << m << "," << n << "," << o << std::endl;

    /*cout<<"vektor: ";
    for (int i=0;i<m;i++) {
      cout<<aa[i]<<" ";
    }
    cout<<endl;*/
    /*aa.assign(taa,taa+m);
    bb.assign(tbb,tbb+n);
    cc.assign(tcc,tcc+o);*/
    //Calculate coefficient matrix for the interpolation of the different entries in E
    int nrow;
    if (mode == 2 || mode == 3 || mode == 4 || mode == 5) {
      nrow = 9;
    } else {
      nrow = 6 + 1;
    }
    vector<vector<vector<double > > > Coeff(nrow, vector<vector<double> >((m-1)*(n-1)*(o-1),vector<double>(64,0.)));
    vector<vector<vector<double > > > Coeff_deriv(nrow, vector<vector<double> >((m-1)*(n-1)*(o-1),vector<double>(64,0.)));
    vector<vector<vector <double > > > E(m,vector<vector<double> >(n,vector<double>(o,0.)));
    vector<vector<double > > Coeff_vol((m-1)*(n-1)*(o-1),vector<double>(64,0.));
    vector<vector<double>> vol_deriv((m-1)*(n-1)*(o-1),vector<double>(64,0.));
    if (fin_vol.is_open()) {
      vector<vector<vector <double > > > E_vol(m,vector<vector<double> >(n,vector<double>(o,0.)));
      for (int i=0;i<m*n*o;i++) {
        E_vol[data_vol[i][0]][data_vol[i][1]][data_vol[i][2]] = data_vol[i][3];
      }
      tricubic_offline(Coeff_vol, aa, bb, cc, E_vol, m, n, o, da, db, dc,vol_deriv);
    }
    for (int ll = 0;ll<nrow;ll++) {
//    for (int ll = 5;ll<6;ll++) {
      //    for (int ll = 0;ll<4;ll++) {
      // k22 -> Coeff[3]
//    for (int ll = 3;ll<4;ll++) {
      for (int i=0;i<m*n*o;i++) {
        E[data[i][0]][data[i][1]][data[i][2]] = data[i][3+ll];
      }

      tricubic_offline(Coeff[ll], aa, bb, cc, E, m, n, o, da, db, dc, Coeff_deriv[ll]);
    }
    /*for (int i=0;i<o;i++) {
        cout<<"E("<<i<<") = "<<endl;
        for (int j=0;j<n;j++) {
          for (int k=0;k<m;k++) {
            cout<<" "<<E[i][j][k];
          }
          cout<<endl;
        }
    }*/

    /*for (int i=0;i<m;i++) {
      cout<<" "<<aa[i];
    }
    cout<<endl;*/


    cout<<"Insert output xml file name:"<<endl;
    string name;
    cin>>name;

    // open output xml
    string file(name);
    fstream f;
    f.open(file.c_str(), ios::out);
    //Write interpolation coefficients in the xml file above
    if (mode != 5) {
      write_to_xml(f,Coeff,aa,bb,cc, mode,Coeff_deriv);
      write_to_xml_vol(f, Coeff_vol, vol_deriv, aa,bb, cc);
      f<<"</homRectC1>";
    }
    f.close();
    /*for (int i=0;i<m-1;i++) {
      for (int j=0;j<n-1;j++) {
        for (int k=0;k<o-1;k++) {
          for (int l=0;l<64;l++) {
            cout<<" "<<Coeff[(n-1)*(o-1)*i+(o-1)*j+k][i];
          }
          cout<<endl;
        }
      }
    }*/

    if (mode == 3) {
      //Select the correct interval for x1, x2 and x3
      int a1=-1;
      int b1=-1;
      int c1 =-1;
      for (int i=0;i<m;i++) {
        if (aa[i] <= x1 && x1 < aa[i+1]) {
          a1=i;
        } else if( x1 == aa[m-1]) {
          a1=m-2;
          break;
        } else if (x1 > aa[m-1]) {
          cout<<"x1 out of bounds"<<endl;
          break;
        }
      }

      for (int i=0;i<n;i++) {
        if (bb[i] <= x2 && x2 < bb[i+1]) {
          b1=i;
        } else if( x2 == bb[n-1]) {
          b1=n-2;
          break;
        } else if (x2 > bb[n-1]) {
          cout<<"x2 out of bounds"<<endl;
          break;
        }
      }
      for (int i=0;i<o;i++) {
        if (cc[i] <= x3 && x3 < cc[i+1]) {
          c1=i;
        } else if( x3 == cc[o-1]) {
          c1=o-2;
          break;
        } else if (x3 > cc[o-1]) {
          cout<<"x3 out of bounds"<<endl;
          break;
        }
      }

      double coeff1=(x1-aa[a1])/da;
      double coeff2=(x2-bb[b1])/db;
      double coeff3=(x3-cc[c1])/dc;

      cout << "a1,b1,c1:" << a1 << "," << b1 << "," << c1 << endl;
      cout << "aa,bb,cc:" <<  aa[a1] << "," << bb[b1] << "," << cc[c1] << endl;
      cout << "coeff1,coeff2,coeff3: " << coeff1 << " " << coeff2 << " " << coeff3 << endl;

      vector<double> a(64,0.);
      vector<double> b(64,0.);
      vector<double> c(64,0.);
      vector<double> d(64,0.);
      vector<double> e(64,0.);
      vector<double> f(64,0.);
      vector<double> g(64,0.);
      vector<double> h(64,0.);
      vector<double> x(64,0.);
      vector<double> v(64,0.);
      for (int i=0;i<64;i++) {
        a[i] = Coeff[0][(n-1)*(o-1)*a1+(o-1)*b1+c1][i];
        b[i] = Coeff[1][(n-1)*(o-1)*a1+(o-1)*b1+c1][i];
        c[i] = Coeff[2][(n-1)*(o-1)*a1+(o-1)*b1+c1][i];
        d[i] = Coeff[3][(n-1)*(o-1)*a1+(o-1)*b1+c1][i];
        e[i] = Coeff[4][(n-1)*(o-1)*a1+(o-1)*b1+c1][i];
        f[i] = Coeff[5][(n-1)*(o-1)*a1+(o-1)*b1+c1][i];
        g[i] = Coeff[6][(n-1)*(o-1)*a1+(o-1)*b1+c1][i];
        h[i] = Coeff[7][(n-1)*(o-1)*a1+(o-1)*b1+c1][i];
        x[i] = Coeff[8][(n-1)*(o-1)*a1+(o-1)*b1+c1][i];
        v[i] = Coeff_vol[(n-1)*(o-1)*a1+(o-1)*b1+c1][i];
      }

      //      cout << "a: " << toString(a) << endl;
      //Evaluate the interpolation interval at point x1,x2,x3
//            double k11 = tricubic_eval(a, coeff1, coeff2, coeff3);
//            double k22 = tricubic_eval(d, coeff1, coeff2, coeff3);
//            double k33 = tricubic_eval(f, coeff1, coeff2, coeff3);
      //      cout<<"k11: " << k11 << endl;

      //      cout<<"e12: "<<tricubic_eval(b, coeff1, coeff2, coeff3)<<endl;
      //      cout<<"e13: "<<tricubic_eval(c, coeff1, coeff2, coeff3)<<endl;
      //        cout<<"e22: " << k22 << endl;
      //      cout<<"e23: "<<tricubic_eval(e, coeff1, coeff2, coeff3)<<endl;
      //        cout<<"e33: " << k33 << endl;
      //      cout<<"e44: "<<tricubic_eval(g, coeff1, coeff2, coeff3)<<endl;
      //      cout<<"e55: "<<tricubic_eval(h, coeff1, coeff2, coeff3)<<endl;
      //      cout<<"e66: "<<tricubic_eval(x, coeff1, coeff2, coeff3)<<endl;
      //      cout<<"vol: "<<tricubic_eval(v, coeff1, coeff2, coeff3)<<endl;
      return 0;
    }
    if (mode == 5) {
      //vector<vector<double > > vol_values(101*101,vector<double>(3,0.));
      fstream f,f2;
      f.open("output_vol.dat", ios::out);
      f2.open("output_kreuze_vol.dat", ios::out);
      for (int ii = 0;ii <= 100;ii++) {
        for (int jj = 0;jj <= 100;jj++) {
          x1 = ii/100.;
          x2 = jj/100.;
          x3 = 0.;
          //Select the correct interval for x1, x2 and x3
          int a1=-1;
          int b1=-1;
          int c1 =-1;
          for (int i=0;i<m;i++) {
            if (aa[i] <= x1 && x1 < aa[i+1]) {
              a1=i;
            } else if( x1 == aa[m-1]) {
              a1=m-2;
              break;
            } else if (x1 > aa[m-1]) {
              cout<<"x1 out of bounds"<<endl;
              break;
            }
          }

          for (int i=0;i<n;i++) {
            if (bb[i] <= x2 && x2 < bb[i+1]) {
              b1=i;
            } else if( x2 == bb[n-1]) {
              b1=n-2;
              break;
            } else if (x2 > bb[n-1]) {
              cout<<"x2 out of bounds"<<endl;
              break;
            }
          }
          for (int i=0;i<o;i++) {
            if (cc[i] <= x3 && x3 < cc[i+1]) {
              c1=i;
            } else if( x3 == cc[o-1]) {
              c1=o-2;
              break;
            } else if (x3 > cc[o-1]) {
              cout<<"x3 out of bounds"<<endl;
              break;
            }
          }
          // Map x1,x2,x3 into the chosen intervall
          double coeff1=(x1-aa[a1])/da;
          double coeff2=(x2-bb[b1])/db;
          double coeff3=(x3-cc[c1])/dc;
          f<<x1<<" "<<x2<< " " <<x3<<" "<<tricubic_eval(Coeff_vol[(n-1)*(o-1)*a1+(o-1)*b1+c1], coeff1, coeff2, coeff3)<<endl;
          f2<<x1<<" "<<x2<< " " <<x3<<" "<<Calc3DCrossVolume(x1,x2,x3,false,1)<<endl;
        }
      }
      f.close();
      f2.close();
      return 0;
    }
    return 0;
  } else {
    cout<<"wrong interpolation goal number chosen."<<endl;
    return -1;
  }
}
