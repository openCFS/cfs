#include <iostream>
#include "math.h"

#include "scpip30.hh"



/* ********************************************************************************* */
/* Subroutine */ int scpfct_(int* n, int* mie, int* meq, double* x, double* f_org__, double* g_org__, double* h_org__)
{
/*     input:  n,mie,meq,x */
/*     output: f_org,g_org,h_org */
    /* Parameter adjustments */
    --h_org__;
    --g_org__;
    --x;

    /* Function Body */
    *f_org__ = -x[1];
    h_org__[1] = exp(x[1]) - x[2];
    h_org__[2] = exp(x[2]) - x[3];
    return 0;
} /* scpfct_ */

/* ********************************************************************************* */
/* Subroutine */ int scpgrd_(int* n, int* mie, int* meq, double* x, double* f_org__, double* g_org__, double* h_org__, 
  int* active, double* df, int* iern, int* iecn, double* iederv, int* ieleng, int* eqrn, int* eqcn, double* eqcoef, int* eqleng)
{

    /* Local variables */
    static int numbercon;

/*     input:  n,mie,meq,x,f_org,g_org,h_org,active */
/*     output: df,iern,iecn,iederv,ieleng,eqrn,eqcn,eqcoef,eqleng */
    /* Parameter adjustments */
    --eqcoef;
    --eqcn;
    --eqrn;
    --iederv;
    --iecn;
    --iern;
    --df;
    --active;
    --h_org__;
    --g_org__;
    --x;

    /* Function Body */
    df[1] = -1.;
    df[2] = 0.;
/* gradient of objective fully */
    df[3] = 0.;
    *ieleng = 0;
    numbercon = 0;
    if (active[1] == 1) {
  ++numbercon;
/* constraint 1, component 1: */
/* one more active constraint */
  ++(*ieleng);
/* derivative always > 0 */
  iecn[*ieleng] = numbercon;
/* function number */
  iern[*ieleng] = 1;
/* component number */
  iederv[*ieleng] = exp(x[1]);
/* constraint 1, component 2: */
/* value */
  ++(*ieleng);
/* derivative always > 0 */
  iecn[*ieleng] = numbercon;
/* function number */
  iern[*ieleng] = 2;
/* component number */
  iederv[*ieleng] = -1.;
/* value */
    }
    if (active[2] == 1) {
  ++numbercon;
/* constraint 2, component 2: */
/* one more active constraint */
  ++(*ieleng);
/* derivative always > 0 */
  iecn[*ieleng] = numbercon;
/* function number */
  iern[*ieleng] = 2;
/* component number */
  iederv[*ieleng] = exp(x[2]);
/* constraint 2, component 3: */
/* value */
  ++(*ieleng);
/* derivative always > 0 */
  iecn[*ieleng] = numbercon;
/* function number */
  iern[*ieleng] = 3;
/* component number */
  iederv[*ieleng] = -1.;
/* value */
    }
    return 0;
} /* scpgrd_ */




/* Table of constant values */

static int c__100 = 100;

int f2c_test()
{
    /* Builtin functions */
    int s_wsle(), do_lio(), e_wsle();
    /* Subroutine */ int s_stop();

    /* Local variables */
    static int iecn[100], idim, mode, eqcn[100];
    static Double y_ie__[10];
    static int info[23], rdim, iern[100];
    static Double y_eq__[10];
    static int ierr, eqrn[100];
    static Double spdw[1000];
    static int spiw[1000], nout, n;
    static Double f_org__;
    static int i_scp__[1000];
    static Double h_org__[10], g_org__[10];
    static int i_sub__[1000], iemax;
    static Double r_scp__[1000];
    static int icntl[13];
    static Double r_sub__[1000];
    static int eqmax;
    static Double rinfo[5], rcntl[6], x0[100], df[100], eqcoef[100];
    static int ieleng, active[10], eqleng;
    static Double iederv[100];
    static int mactiv;
    static int linsys, mie, meq;
    static Double x_l__[100], y_l__[100], x_u__[100], y_u__[100];
    static int isubdim, rsubdim, spdwdim, spiwdim, spstrat;


/*     Sample driving program for problem #34 in Hock/Schittkowski */
/* mi inequality constraints */
/* me equality constraints */
/* nn maximum number of variab */
/* Jacobian of inequalities. e */
/* for equalities */
/* ielpar maximum number of en */
/* maximum for working arrays */
    n = 3;
    mie = 2;
    meq = 0;
    iemax = mie > 1 ? mie : 1;
    eqmax = meq > 1 ? meq : 1;
    x0[0] = 0.;
    x0[1] = 1.05;
    x0[2] = 2.9;
    x_l__[0] = 0.;
    x_l__[1] = 0.;
    x_l__[2] = 0.;
    x_u__[0] = 100.;
    x_u__[1] = 100.;
    x_u__[2] = 10.;
    icntl[0] = 1;
    icntl[1] = 1;
    icntl[2] = 100;
    icntl[3] = 1;
    icntl[4] = 10;
    icntl[5] = 0;
    icntl[10] = 0;
    icntl[11] = 0;
    icntl[12] = 0;
    rcntl[0] = 1e-7;
    rcntl[1] = 1e30;
    rcntl[2] = 1e30;
    rcntl[3] = 0.;
    rcntl[4] = 0.;
    rcntl[5] = 0.;
    ierr = 0;
    nout = 7;
    rdim = 1000;
    rsubdim = 1000;
    idim = 1000;
    isubdim = 1000;
    spiwdim = 1000;
    spdwdim = 1000;
    mode = 2;
    spstrat = 1;
    linsys = 1;

    for(;;) // we break out in success and error
    {    
      scpip30_(&n, &mie, &meq, &iemax, &eqmax, x0, x_l__, x_u__, &f_org__, 
        h_org__, g_org__, df, y_ie__, y_eq__, y_l__, y_u__, icntl, rcntl, 
        info, rinfo, &nout, r_scp__, &rdim, r_sub__, &rsubdim, i_scp__, &
        idim, i_sub__, &isubdim, active, &mode, &ierr, iern, iecn, iederv,
         &c__100, &ieleng, eqrn, eqcn, eqcoef, &c__100, &eqleng, &mactiv, 
        spiw, &spiwdim, spdw, &spdwdim, &spstrat, &linsys);

      std::cout << "scpip30 returned with " << ierr << std::endl;
        
      if (ierr == 0) {
        std::cout << "solution obtained" << std::endl;
        break;
      } 
      if (ierr == -1) {
        scpfct_(&n, &mie, &meq, x0, &f_org__, g_org__, h_org__);
      } 
      if (ierr == -2) {
        scpgrd_(&n, &mie, &meq, x0, &f_org__, g_org__, h_org__, active, df, 
                iern, iecn, iederv, &ieleng, eqrn, eqcn, eqcoef, &eqleng);
      } 
      if (ierr == -3) {
        /* now you have the chance to adopt spiwdim and spdwdim with the information */
        /* of info(6) and info(7). we just jump back here. */
      }
      if(ierr < -3 || ierr > 0) {
        std::cout << "scpip finished with an error! ierr = " << ierr << std::endl;
        break;
      }
    }
    
    return ierr;      
} /* MAIN__ */

















void scpfct (int n, int mie, int meq, const double* x,
             Double f_org, double* g_org, double* h_org)
{
  f_org = -1.0 * x[0];
  h_org[0] = exp(x[0]) - x[1];
  h_org[1] = exp(x[1]) - x[2];
}


/*
 * void scpgrd (int n, int mie, int meq,const double* x, Double f_org, const double* g_org, const double* h_org,
             const int* active, double* df, int* iern, int* iecn, double* iederv, int& ieleng,
             int* eqrn, int* eqcn, double* eqcoef, int& eqleng)
{             

  df[0] = -1.0; // gradient of objective fully stored!
  df[1] = 0.0;  
  df[2] = 0.0;

  ieleng = 0;
  int numbercon = 0;

  if(active[0] == 1)
  {
    numbercon = numbercon+1;    // one more active constraint
    // constraint 0, component 0
    ieleng = ieleng + 1;        // derivative always > 0
    iecn[ieleng] = numbercon;   // function number
    iern[ieleng] = 0;           // component number
    iederv[ieleng] = exp(x[0]); // value

    // constraint 0, component 1
    ieleng = ieleng + 1;
    iecn[ieleng] = numbercon;
    iern[ieleng] = 1;
    iederv[ieleng] = -1.0;
  }

  if(active[1] == 1) 
  {
    numbercon=numbercon+1;
    // constraint 1, component 1
    ieleng = ieleng + 1;
    iecn[ieleng] = numbercon;
    iern[ieleng] = 1;
    iederv[ieleng] = exp(x[1]);

    // constraint 1, component 2
    ieleng = ieleng + 1;
    iecn[ieleng] = numbercon;
    iern[ieleng] = 2;     
    iederv[ieleng] = -1.0;
  }
} */



   void scpgrd (int n, int mie, int meq,const double* x, Double f_org, const double* g_org, const double* h_org,
             const int* active, double* df, int* iern, int* iecn, double* iederv, int& ieleng,
             int* eqrn, int* eqcn, double* eqcoef, int& eqleng)
{             
    /* Function Body */
    df[0] = -1.;
    df[1] = 0.;
/* gradient of objective fully */
    df[2] = 0.;
    ieleng = 0;
     static int numbercon = 0;
    if (active[0] == 1) {
  ++numbercon;
/* constraint 1, component 1: */
/* one more active constraint */
  ++ieleng;
/* derivative always > 0 */
  iecn[ieleng] = numbercon;
/* function number */
  iern[ieleng] = 1;
/* component number */
  iederv[ieleng] = exp(x[0]);
/* constraint 1, component 2: */
/* value */
  ++(ieleng);
/* derivative always > 0 */
  iecn[ieleng] = numbercon;
/* function number */
  iern[ieleng] = 2;
/* component number */
  iederv[ieleng] = -1.;
/* value */
    }
    if (active[1] == 1) {
  ++numbercon;
/* constraint 2, component 2: */
/* one more active constraint */
  ++(ieleng);
/* derivative always > 0 */
  iecn[ieleng] = numbercon;
/* function number */
  iern[ieleng] = 2;
/* component number */
  iederv[ieleng] = exp(x[1]);
/* constraint 2, component 3: */
/* value */
  ++(ieleng);
/* derivative always > 0 */
  iecn[ieleng] = numbercon;
/* function number */
  iern[ieleng] = 3;
/* component number */
  iederv[ieleng] = -1.;
/* value */
    }
}




void test()
{
  // this is a plain fortan -> c translation of v30main.f written by Ch. Zillober
  // Sample driving program for problem #34 in Hock/Schittkowski

  int nn,mi,me,ielpar,eqlpar,rd,rd2,rd3,id,id2,id3;
  
  nn = 100; // maximum number of variables,
  mi = 10;  // inequality constraints
  me = 10;  // equality constraints

  ielpar=mi*10; // maximum number of entries in Jacobian of inequalities  
  eqlpar=me*10; // same for equalities

  rd = rd2 = rd3 = id = id2 = id3 = 1000; // maximum for working arrays

  int n,mie,meq,iemax,eqmax,nout,rdim,rsubdim,idim,isubdim,
      mode,ierr,ieleng,eqleng,mactiv,spiwdim,spdwdim,spstrat,linsys;
      
  int icntl[13];
  int info[23];    
  int* i_scp = new int[id];
  int* i_sub = new int[id2];
  int* active = new int[mi];
  int* iern = new int[ielpar];
  int* iecn = new int[ielpar];
  int* eqrn = new int[eqlpar];
  int* eqcn = new int[eqlpar];
  int* spiw = new int[id3];
      
  Double f_org;

  Double rcntl[6];
  Double rinfo[5];       
         
  double* x0 = new double[nn];
  double* x_l = new double[nn];
  double* x_u = new double[nn];
  double* h_org = new double[mi];
  double* g_org = new double[me];
  double* df = new double[nn];  
  double* y_ie = new double[mi];  
  double* y_eq = new double[me];
  double* y_l = new double[nn];
  double* y_u = new double[nn];      
  double* r_scp = new double[rd];
  double* r_sub = new double[rd2];
  double* iederv = new double[ielpar];
  double* eqcoef = new double[eqlpar];
  double* spdw = new double[rd3];      
      
  n = 3;
  mie = 2;
  meq = 0;

  iemax = mie > 1 ? mie : 1;
  eqmax = meq > 1 ? meq : 1;

  x0[0] = 0.0;
  x0[1] = 1.05;
  x0[2] = 2.9;

  x_l[0] = 0.0;
  x_l[1] = 0.0;
  x_l[2] = 0.0;

  x_u[0] = 100;
  x_u[1] = 100;
  x_u[2] = 10;

  icntl[0] = 1;
  icntl[1] = 1;
  icntl[2] = 100;
  icntl[3] = 1;
  icntl[4] = 10;
  icntl[5] = 0;
  icntl[10] = 0;
  icntl[11] = 0;
  icntl[12] = 0;

  rcntl[0] = 1.e-7;
  rcntl[1] = 1.e30;
  rcntl[2] = 1.e30;
  rcntl[3] = 0.0;
  rcntl[4] = 0.0;
  rcntl[5] = 0.0;

  ierr = 0;
  nout = 7;
  rdim = rd;
  rsubdim = rd2;
  idim = id;
  isubdim = id2;
  spiwdim = id3;
  spdwdim = rd3;

  mode = 2;

  spstrat = 1;
  linsys = 1;
/*
  int n,mie,meq,iemax,eqmax,nout,rdim,rsubdim,idim,isubdim,
      mode,ierr,ieleng,eqleng,mactiv,spiwdim,spdwdim,spstrat,linsys,i;
 int nn,mi,me,ielpar,eqlpar,rd,rd2,rd3,id,id2,id3;*/

  bool loop = true;
  while(loop)
  {
    scpip30_(&n,&mie,&meq,&iemax,&eqmax,x0,x_l,x_u,&f_org,h_org,
            g_org,df,y_ie,y_eq,y_l,y_u,icntl,rcntl,info,
            rinfo,&nout,r_scp,&rdim,r_sub,&rsubdim,i_scp,&idim,
            i_sub,&isubdim,active,&mode,&ierr,iern,iecn,iederv,
            &ielpar,&ieleng,eqrn,eqcn,eqcoef,&eqlpar,&eqleng,&mactiv,
            spiw,&spiwdim,spdw,&spdwdim,&spstrat,&linsys);
    std::cout << "scpip returns " << ierr << std::endl;
    if(ierr == 0) 
    {
      std::cout << "solution obtained" << std::endl;
      loop = false;
    }
    if(ierr == -1) 
    {
      scpfct(n,mie,meq,x0,f_org,g_org,h_org);
    }
    if(ierr == -2) 
    {
       scpgrd(n,mie,meq,x0,f_org,g_org,h_org,active,df,iern,
               iecn,iederv,ieleng,eqrn,eqcn,eqcoef,eqleng);
    }
    if(ierr == -3)
    {
       // now you have the chance to adopt spiwdim and spdwdim with the information
       // of info(6) and info(7). we just jump back here.
    }
    if(ierr < -3 || ierr > 0)
    {
      std::cout << "scpip finished with an error! ierr = " << ierr << std::endl;
      loop = false;
    }
  }  
}

int main(int argc, char* argv[])
{
  //test();
  f2c_test();
}
