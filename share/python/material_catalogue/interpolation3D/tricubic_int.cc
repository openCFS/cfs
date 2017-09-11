
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
int ijk2n(int i, int j, int k) {
  return(i+4*j+16*k);
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


  int m = a.size();
  int n = b.size();
  int o = c.size();
  for (int i=1;i<m-1;i++) {
    for (int j=1;j<n-1;j++) {
      for (int k=1;k<o-1;k++) {
        //f(i+1,j,k)-f(i-1,j,k) / 2
        dEda[i][j][k] = (E[i+1][j][k] - E[i-1][j][k]) /2.;
        dEdb[i][j][k] = (E[i][j+1][k] - E[i][j-1][k]) /2.;
        dEdc[i][j][k] = (E[i][j][k+1] - E[i][j][k-1]) /2.;
        //f(i+1,j+1,k)-f(i+1,j-1,k)-f(i-1,j+1,k)+f(i-1,j-1,k))/4
        dEdadb[i][j][k] = (E[i+1][j+1][k]-E[i+1][j-1][k] - E[i-1][j+1][k]+E[i-1][j-1][k])/4.;
        dEdadc[i][j][k] = (E[i+1][j][k+1]-E[i+1][j][k-1] - E[i-1][j][k+1]+E[i-1][j][k-1])/4.;
        dEdbdc[i][j][k] = (E[i][j+1][k+1]-E[i][j+1][k-1] - E[i][j-1][k+1]+E[i][j-1][k-1])/4.;
        //(f(i+1,j+1,z+1) - f(i+1,j,k-1) - f(i+1,j-1,k+1) + f(i+1,j-1,k-1) - f(i-1,j+1,k+1) + f(i-1,j+1,k-1)
        //+ f(i-1,j-1,k+1) - f(i-1,j-1,k-1))/8
        dEdadbdc[i][j][k] = (E[i+1][j+1][k+1]-E[i+1][j+1][k-1]-E[i+1][j-1][k+1]+E[i+1][j-1][k-1]-E[i-1][j+1][k+1]+
            E[i-1][j+1][k-1]+E[i-1][j-1][k+1]-E[i-1][j-1][k-1])/8.;
      }
    }
  }

}

void tricubic_offline(vector<vector<double> > & Coeff,vector<double> a, vector<double> b, vector<double> c, vector<vector<vector<double> > > & E, int m, int n, int o, double dx,double dy,double dz) {
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
  for (int i=0;i<m-1;i++) {
    for (int j=0;j<n-1;j++) {
      for (int k =0;k<o-1;k++) {
                    f[0]= E[i][j][k];
                    f[1]= E[i+1][j][k];
                    f[2]= E[i][j+1][k];
                    f[3]= E[i+1][j+1][k];

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
                    //scale values
                    for (int q=0;q<8;q++) {
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
                      Coeff[(n-1)*(o-1)*i+(o-1)*j+k][ii] = res[ii];
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
        ret+=a[ijk2n(i,j,k)]*pow(x,i)*pow(y,j)*pow(z,k);
      }
    }
  }
  return ret;
}

void write_to_xml(fstream & f, vector<vector<vector<double> > > Coeff, vector<double> aa, vector<double> bb, vector<double> cc, int mode) {
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
  int ende = (aa.size()-1)*(bb.size()-1)*(cc.size()-1);
  f<<"<coeff11>\n<matrix dim1=\""<<ende<<"\" dim2=\"64\">\n<real>\n";
  for (int i=0;i<ende;i++) {
    for (int j=0;j<64;j++) {
      f<<setprecision(12)<<Coeff[0][i][j]<<" ";
    }
    f<<endl;
  }
  f<<"</real>\n</matrix>\n</coeff11>\n";
  f<<"<coeff12>\n<matrix dim1=\""<<ende<<"\" dim2=\"64\">\n<real>\n";
  for (int i=0;i<ende;i++) {
    for (int j=0;j<64;j++) {
      f<<setprecision(12)<<Coeff[1][i][j]<<" ";
    }
    f<<endl;
  }
  f<<"</real>\n</matrix>\n</coeff12>\n";
  f<<"<coeff13>\n<matrix dim1=\""<<ende<<"\" dim2=\"64\">\n<real>\n";
  for (int i=0;i<ende;i++) {
    for (int j=0;j<64;j++) {
      f<<setprecision(12)<<Coeff[2][i][j]<<" ";
    }
    f<<endl;
  }
  f<<"</real>\n</matrix>\n</coeff13>\n";
  f<<"<coeff22>\n<matrix dim1=\""<<ende<<"\" dim2=\"64\">\n<real>\n";
  for (int i=0;i<ende;i++) {
    for (int j=0;j<64;j++) {
      f<<setprecision(12)<<Coeff[3][i][j]<<" ";
    }
    f<<endl;
  }
  f<<"</real>\n</matrix>\n</coeff22>\n";
  f<<"<coeff23>\n<matrix dim1=\""<<ende<<"\" dim2=\"64\">\n<real>\n";
  for (int i=0;i<ende;i++) {
    for (int j=0;j<64;j++) {
      f<<setprecision(12)<<Coeff[4][i][j]<<" ";
    }
    f<<endl;
  }
  f<<"</real>\n</matrix>\n</coeff23>\n";
  f<<"<coeff33>\n<matrix dim1=\""<<ende<<"\" dim2=\"64\">\n<real>\n";
  for (int i=0;i<ende;i++) {
    for (int j=0;j<64;j++) {
      f<<setprecision(12)<<Coeff[5][i][j]<<" ";
    }
    f<<endl;
  }
  f<<"</real>\n</matrix>\n</coeff33>\n";
  f<<"<coeff44>\n<matrix dim1=\""<<ende<<"\" dim2=\"64\">\n<real>\n";
  for (int i=0;i<ende;i++) {
    for (int j=0;j<64;j++) {
      f<<setprecision(12)<<Coeff[6][i][j]<<" ";
    }
    f<<endl;
  }
  f<<"</real>\n</matrix>\n</coeff44>\n";
    f<<"<coeff55>\n<matrix dim1=\""<<ende<<"\" dim2=\"64\">\n<real>\n";
    for (int i=0;i<ende;i++) {
      for (int j=0;j<64;j++) {
        f<<setprecision(12)<<Coeff[7][i][j]<<" ";
      }
      f<<endl;
    }
    f<<"</real>\n</matrix>\n</coeff55>\n";
    f<<"<coeff66>\n<matrix dim1=\""<<ende<<"\" dim2=\"64\">\n<real>\n";
    for (int i=0;i<ende;i++) {
      for (int j=0;j<64;j++) {
        f<<setprecision(12)<<Coeff[8][i][j]<<" ";
      }
      f<<endl;
    }
    f<<"</real>\n</matrix>\n</coeff66>\n";
}

void write_to_xml_vol(fstream & f, vector<vector<double>  > Coeff, vector<double> aa, vector<double> bb, vector<double> cc, bool combined = true) {
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
  int ende = (aa.size()-1)*(bb.size()-1)*(cc.size()-1);
  f<<"<volcoeff>"<<endl;
  f<<"<matrix dim1=\""<<ende<<"\" dim2=\"64\">"<<endl;
  f<<"<real>"<<endl;
  for (int i=0;i<ende;i++) {
    for (int j=0;j<64;j++) {
      f<<setprecision(12)<<Coeff[i][j]<<" ";
    }
    f<<endl;
  }
  f<<"</real>"<<endl;
  f<<"</matrix>"<<endl;
  f<<"</volcoeff>"<<endl;
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
    tricubic_offline(Coeff,aa, bb, cc, E,m, n, o,da,db,dc);
    cout<<"Insert output xml file name:"<<endl;
    string name;
    cin>>name;
    string file(name);
    fstream f;
    f.open(file.c_str(), ios::out);
    //Write interpolation coefficients in the xml file above
    write_to_xml_vol(f,Coeff,aa,bb,cc,false);
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

    double da = 1./(m);
    double db = 1./(n);
    double dc = 1./(o);
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
    vector<vector<vector <double > > > E(m,vector<vector<double> >(n,vector<double>(o,0.)));
    vector<vector<double > > Coeff_vol((m-1)*(n-1)*(o-1),vector<double>(64,0.));
    if (fin_vol.is_open()) {
      vector<vector<vector <double > > > E_vol(m,vector<vector<double> >(n,vector<double>(o,0.)));
      for (int i=0;i<m*n*o;i++) {
            E_vol[data_vol[i][0]][data_vol[i][1]][data_vol[i][2]] = data_vol[i][3];
      }
      tricubic_offline(Coeff_vol, aa, bb, cc, E_vol, m, n, o, da, db, dc);
    }
    for (int ll = 0;ll<nrow;ll++) {
      for (int i=0;i<m*n*o;i++) {
            E[data[i][0]][data[i][1]][data[i][2]] = data[i][3+ll];
      }
      tricubic_offline(Coeff[ll], aa, bb, cc, E, m, n, o, da, db, dc);
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
      write_to_xml(f,Coeff,aa,bb,cc, mode);
      write_to_xml_vol(f, Coeff_vol, aa,bb, cc);
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
      // Map x1,x2,x3 into the chosen interval
      double coeff1=(x1-aa[a1])/da;
      double coeff2=(x2-bb[b1])/db;
      double coeff3=(x3-cc[c1])/dc;

      vector<double> a(64,0.);
      vector<double> b(64,0.);
      vector<double> c(64,0.);
      vector<double> d(64,0.);
      vector<double> e(64,0.);
      vector<double> f(64,0.);
      vector<double> g(64,0.);
      vector<double> h(64,0.);
      vector<double> x(64,0.);
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
      }

      //Evaluate the interpolation interval at point x1,x2,x3
      result= tricubic_eval(a, coeff1, coeff2, coeff3);
      cout<<"e11: "<<result<<endl;
      cout<<"e12: "<<tricubic_eval(b, coeff1, coeff2, coeff3)<<endl;
      cout<<"e13: "<<tricubic_eval(c, coeff1, coeff2, coeff3)<<endl;
      cout<<"e22: "<<tricubic_eval(d, coeff1, coeff2, coeff3)<<endl;
      cout<<"e23: "<<tricubic_eval(e, coeff1, coeff2, coeff3)<<endl;
      cout<<"e33: "<<tricubic_eval(f, coeff1, coeff2, coeff3)<<endl;
      cout<<"e44: "<<tricubic_eval(g, coeff1, coeff2, coeff3)<<endl;
      cout<<"e55: "<<tricubic_eval(h, coeff1, coeff2, coeff3)<<endl;
      cout<<"e66: "<<tricubic_eval(x, coeff1, coeff2, coeff3)<<endl;
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

