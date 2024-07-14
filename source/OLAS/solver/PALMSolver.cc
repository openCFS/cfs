/*

The code in this section is taken from the work of Ding Lu: http://www.unige.ch/~dlu/palm.html

Publications using this work should cite:
A Padé approximate linearization algorithm for solving the quadratic eigenvalue problem with low-rank damping
by Ding Lu, Xin Huang, Zhaojun Bai, and Yangfeng Su
Int. J. Numer. Methods Eng., 2015. 103(11): 840–858. (paper)

The original code was adapted to CFS by Florian Schlägl.

*/

#include "PALMSolver.hh"
#include "Utils/Timer.hh"
#include "MatVec/BLASLAPACKInterface.hh"
DEFINE_LOG(palm, "palm")

namespace CoupledField{

  PALMEigenSolver::PALMEigenSolver(shared_ptr<SolStrategy> strat, PtrParamNode xml, PtrParamNode solverList,  PtrParamNode precondList, PtrParamNode eigenInfo)
  : BaseEigenSolver(strat, xml, solverList, precondList, eigenInfo)
  {
    this->generalized_ = false;
    this->bloch_ = false;
    this->numFreq_ = 0;
    this->freqShift_ = 0;
    xml_ = xml;
    eigenSolverName_ = EigenSolverType::PALM;

    // Data from XML
    xml_->GetValue("PadeOrder", p, ParamNode::INSERT);

    isscaling = iseigv = isresid = true;
    isfactor = resdone = isscaled = false;

    normm = normd = normk = 0;

    /* C) Parameter Initialization () */

    threshold = 0.1;
    artol = 1.0E-16;
    nconv = -1;
    delta = 1.0;
    order = 1;

    // set a lot of things to NULL to make Eclipse and compiler happy
    EigVal = EigVec = NULL;
    res = NULL;
    rk = 0;

    // set in lufactor
    permr = permc = NULL;

    // (re)set in CalcEigenValues
    nzmax = n = nl = 0;
    c_ = k_ = m_ = NULL;
    c_CRS = k_CRS = m_CRS = NULL;

    // (re)set in CalcEigenValues
    a = b = NULL;
    cidx = ridx = NULL;
    nev = 0;

    BaseEigenSolver::PostInit();
  }

  PALMEigenSolver::~PALMEigenSolver()
  {
    LOG_DBG(palm) << "destructor";
    delete[] EigVal;
    delete[] EigVec;
    delete[] res;
    delete[] permr;
    delete[] permc;
    delete[] a;
    delete[] b;
    delete[] cidx;
    delete[] ridx;
    c_ = k_ = m_ = NULL;
    c_CRS = k_CRS = m_CRS = NULL;
    xml_ = NULL;
  }

  void PALMEigenSolver::ToInfo()
  {

  }

  void PALMEigenSolver::Setup(const BaseMatrix & K, const BaseMatrix & C, const BaseMatrix & M)
  {
    // Array containing the nonzero elements of either the full matrix A or the upper or lower
    // triangular part of the matrix A, as specified by uplo.
    m_ = &dynamic_cast<const StdMatrix&>(M);
    const SCRS_Matrix<double>* m_const = dynamic_cast<const SCRS_Matrix<double>*>(m_);
    m_CRS = const_cast<SCRS_Matrix<double>*>(m_const);

    c_ = &dynamic_cast<const StdMatrix&>(C);
    const SCRS_Matrix<double>* c_const = dynamic_cast<const SCRS_Matrix<double>*>(c_);
    c_CRS = const_cast<SCRS_Matrix<double>*>(c_const);


    k_ = &dynamic_cast<const StdMatrix&>(K);
    const SCRS_Matrix<double>* k_const = dynamic_cast<const SCRS_Matrix<double>*>(k_);
    k_CRS = const_cast<SCRS_Matrix<double>*>(k_const);

    this->SetProblemType(M,false);

    n = K.GetNumRows();

    /* C) Parameter Initialization */
    nzmax = m_->GetNnz()+k_->GetNnz()+c_->GetNnz();
  }

  void PALMEigenSolver::LowRankDcomp()
  {

    /*  A) Get the row indices I and column indices J,
     *      such that X(I,J) forms the nonzero part
     *      the matrix X.
     */

    int nnz = c_CRS->GetNnz();
    ridx = new UInt[ nnz ]; // Row indices.
    cidx = new UInt[ nnz ]; // Column indices.
    //nr    // Number of nonzero rows.
    //nc    // Number of nonzero columns.
    int t;    // Internal variables.

    double temp_new = 0;
    int nr = 0;
    int nc = 0;
    UInt col_count = c_CRS->GetNumCols();
    //bool col_notzero[col_count] = {false};
    StdVector<bool> col_notzero; col_notzero.Resize(col_count,false);
    UInt row_temp=0;
    for(UInt i = 0; i<(c_CRS->GetNnz()+c_CRS->GetNumCols())/2; i++) // loop thought all entries
    {
      if(i == c_CRS->GetRowPointer()[row_temp+1])
        row_temp++;
      UInt j = c_CRS->GetColPointer()[i];
      c_CRS->GetMatrixEntry(row_temp, j, temp_new);
      LOG_DBG3(palm) << i << ":" << row_temp << "," << j << " = " << temp_new;
      if(temp_new != 0)
      {
        LOG_DBG3(palm) << "  ^ nonzero found ^ ";
        if(row_temp != c_CRS->GetColPointer()[i])
        {
          LOG_DBG3(palm) << "  off diagonal";
          col_notzero[c_CRS->GetColPointer()[i]] = true;
          col_notzero[row_temp] = true;
          ridx[nr] = row_temp;
          ridx[nr+1] = c_CRS->GetColPointer()[i];
          cidx[nr] = c_CRS->GetColPointer()[i];
          cidx[nr+1] = row_temp;
          nr += 2;
        }
        else
        {
          LOG_DBG3(palm) << "  diagonal";
          col_notzero[c_CRS->GetColPointer()[i]] = true;
          col_notzero[row_temp] = true;
          ridx[nr] = row_temp;
          cidx[nr] = c_CRS->GetColPointer()[i];
          nr++;
        }
      }
      LOG_DBG3(palm) << " nr = " << nr ;
    }

    for(UInt i=0; i<col_notzero.GetSize(); i++)
    {
      if(col_notzero[i] == true)
        nc++;
    }

    /*  B) Form a dense matrix B = A(I,J). */

    double *mat = new double[ nr * nc ];


    for( int i=0; i<nr; i++ ) {
      for( int j=0; j<nc; j++ ) {
        double temp = 0;
        if(ridx[i]<cidx[i])
          c_CRS->GetMatrixEntry(ridx[i], cidx[j], temp);
        else
          c_CRS->GetMatrixEntry(cidx[i], ridx[j], temp);
        mat[ j * nc + i] = temp;
      }
    }
    /*  C) Apply SVD to B.   */
    t = (nr<nc)? nc : nr;  // t = max( nc , nr )

    int lwork = 5*t;  // SVD WORKING SPACE
    int info;
    double *MS = new double(t);
    double *MU = new double[nr*nc];
    double *MVT = new double[nr*nc];
    double *work = new double[lwork];

    char s = 'S';
    dgesvd( &s, &s, &nr, &nc, mat, &nr, MS, MU, &nr, MVT, &nc, work, &lwork, &info);
    /*  D) generate E and F.   */

    //  Check numerical rank.
    rk = (nr<nc)? nr : nc;
    for( int i=0; i<t; i++ ) {
      if( MS[ i ] < MS[ 0 ]*10E-16) {
        rk = i + 1;
        break;
      }
    }

    double tnm;
    int inc = 1;

    t = (nr<nc)? nr : nc;  // t = min( nc , nr )

    for(int i=0; i<rk;i++){
      tnm = sqrt( MS[ i ] );
      dscal(&nr, &tnm, &MU[ i*nr ], &inc); // Scale column of u
      dscal(&nc, &tnm, &MVT[i], &t);  // Scale row of vt
    }

    // Construct the low-rank matrix in CSR format.
    UInt* temp_colptr_E = new UInt[rk];
    for(int i = 0; i<rk; i++)
      temp_colptr_E[i] = i*rk;

    UInt* temp_rowptr_F = new UInt[rk];
    for(int i = 0; i<rk; i++)
      temp_rowptr_F[i] = i*rk;

    E = Matrix<double>(n, rk);
    F = Matrix<double>(rk, n);

    for(int i=0; i<rk; i++)
    {
      E.SetEntry(ridx[i],i , MU[i]);
      F.SetEntry(i, cidx[i], MVT[i]);
    }
    LOG_DBG(palm) << "LowRankDcomp: END";
    return;
  }

  double PALMEigenSolver::
  Norm(const SCRS_Matrix<double>* Matr)
  {

    double* nm_col = new double[Matr->GetNumCols()] ;
    double nm = 0;
    int row_temp=0;
    for(UInt i = 0; i<(Matr->GetNnz()+Matr->GetNumCols())/2; i++)
    {
      if(i == Matr->GetRowPointer()[row_temp+1])
        row_temp++;

      double temp = 0;
      Matr->GetMatrixEntry(row_temp, Matr->GetColPointer()[i], temp);
      nm_col[Matr->GetColPointer()[i]] = nm_col[Matr->GetColPointer()[i]] + abs(temp);

      if(i>0)
      {
        if(nm_col[Matr->GetColPointer()[i]] > nm)
          nm = nm_col[Matr->GetColPointer()[i]];
      }
    }
    return nm;
  }

  double PALMEigenSolver::
  scaling() // compute norm estimate and scaling factors
  {
    if(isscaled || !isscaling) return delta; // precheck

    // compute matrix norms
    normm = Norm(m_CRS);
    normd = Norm(c_CRS);
    normk = Norm(k_CRS);

    double t; // internal variable


    // compute scaling factors delta= 1/ max(sigma^2*nm, 2*p*sigma*nd, nk)
    delta = ABS(sigma)*ABS(sigma) * normm;
    t = 2*p*normd*ABS(sigma);
    delta = (delta<t)? t : delta;
    delta = (delta<normk)? normk : delta;
    delta = 1.0/delta;

    isscaled = true;

    return delta;
  }

  void PALMEigenSolver::
  Qsigma(Complex* val, int* rowind, int* colptr, int &nnz, int &info)
  {
    nnz = 0;

    for(int row = 0; row<n; row++)
    {
      colptr[row] = nnz;
      for(int col = 0; col < n; col++)
      {
        double temp_m = 0, temp_k = 0, temp_c = 0;
        int temp_pos_m =0, temp_pos_k = 0, temp_pos_c=0;

        if(row < col)
        {
          temp_pos_m = Calculate_DataArrayPosition(row,col,m_CRS->GetRowPointer(), m_CRS->GetColPointer());
          temp_pos_k = Calculate_DataArrayPosition(row,col,k_CRS->GetRowPointer(), k_CRS->GetColPointer());
          temp_pos_c = Calculate_DataArrayPosition(row,col,c_CRS->GetRowPointer(), c_CRS->GetColPointer());
        }
        else
        {
          temp_pos_m = Calculate_DataArrayPosition(col,row,m_CRS->GetRowPointer(), m_CRS->GetColPointer());
          temp_pos_k = Calculate_DataArrayPosition(col,row,k_CRS->GetRowPointer(), k_CRS->GetColPointer());
          temp_pos_c = Calculate_DataArrayPosition(col,row,c_CRS->GetRowPointer(), c_CRS->GetColPointer());
        }
        if(temp_pos_m != -1)
          temp_m = m_CRS->GetDataPointer()[temp_pos_m];
        if(temp_pos_k != -1)
          temp_k = k_CRS->GetDataPointer()[temp_pos_k];
        if(temp_pos_c != -1)
          temp_c = c_CRS->GetDataPointer()[temp_pos_c];

        val[nnz] = sigma * sigma * temp_m + sigma * temp_c + temp_k;
        if(val[nnz] != Complex(0.0,0.0))
        {
          rowind[nnz] = col;
          nnz++;

        }
      }
    }
    colptr[ n ] = nnz;


    return;
  }// compute Q(sigma)

  int PALMEigenSolver::Calculate_DataArrayPosition(UInt row, UInt col,const UInt* rowPtr,const UInt* colInd)
  {
    UInt row_begin = rowPtr[row];
    UInt row_end = rowPtr[row+1];
    for ( UInt i = row_begin; i < row_end; i++ ) {
      if ( colInd[i] == col ) {
        return i;
      }
    }
    return -1;

  }//Calculates the Position in the CRS_Matrix. If Element is not found returns -1;

  void PALMEigenSolver::
  lufactor() // compute LU factorization of Q(s)
  {

    /* reserve storage for Q(s), in CSC format  */

    Complex *matms = new Complex[ nzmax ];
    int     *rowindms = new int[ nzmax ];
    int    *colptrms = new int[ n+1 ];
    int nnzms, info;


    /* compute Q(sigma)    */

    Qsigma(matms, rowindms, colptrms, nnzms, info);

    /* prepare data for SUPERLU.   */

    SuperMatrix  A;
    zCreate_CompCol_Matrix(&A, n, n, nnzms, (doublecomplex*)matms,
        rowindms, colptrms, SLU_NC, SLU_Z, SLU_GE);

    /*
     * LU factorization to A
     */

    /* Defining local variables. */

    int          *etree;
    SuperMatrix  AC;

    /* Setting default values for gstrf parameters. */

    int     panel_size = sp_ienv( 1 );
    int     relax = sp_ienv( 2 );
    superlu_options_t  options;

    set_default_options( &options );
    options.DiagPivotThresh = threshold;

    /* Reserving memory for etree (used in matrix decomposition). */

    etree = new int[ n ];
    permc = new int[ n ];
    permr = new int[ n ];

    StatInit( &stat);  //StatInit(panel_size, relax);

    get_perm_c(order, &A, permc); // Defining the column permutation of matrix A
    // (using minimum degree ordering on A'*A).
    sp_preorder(&options, &A, permc, etree, &AC);

    zgstrf(&options, &AC, relax, panel_size, etree, NULL, 0, permc, permr, &L, &U, &Glu, &stat, &info);

    isfactor = true;
    /* Deleting AC and etree, and working space */

    Destroy_CompCol_Permuted( &AC );
    delete[]  matms;
    delete[]  rowindms;
    delete[]  colptrms;
    delete[]  etree;

    // POSSIBLE IMPROVEMENTS: check info...
  }

  void PALMEigenSolver ::
  Mult(const SCRS_Matrix<double>* M ,Complex *v, Complex *w, Complex alpha)
  {
    UInt row_temp=0;
    for(UInt i = 0; i<(M->GetNnz()+M->GetNumCols())/2; i++)
    {
      if(i == M->GetRowPointer()[row_temp+1])
        row_temp++;

      double temp = 0;
      M->GetMatrixEntry(row_temp, M->GetColPointer()[i], temp);
      w[row_temp] = w[row_temp]+ alpha*temp*v[M->GetColPointer()[i]];
      if(row_temp != M->GetColPointer()[i])
        w[M->GetColPointer()[i]] = w[M->GetColPointer()[i]] + alpha*temp*v[row_temp];
    }

    return;
  }// Matrix Vector product w=M*v

  void PALMEigenSolver::MultM(Matrix<double>* M ,Complex *v, Complex *w, Complex alpha, char trans)
  {

    if(trans == 'N')
    {
      for(int i=0; i<rk; i++)
      {
        double temp = 0;
        M->GetEntry(ridx[i], i, temp);
        w[cidx[i]] = w[cidx[i]]+ alpha*temp*v[i];
      }
    }
    else
    {
      for(UInt col = 0; col < M->GetNumCols(); col++)
      {
        for(UInt row = 0; row < M->GetNumRows(); row++)
        {
          double temp = 0;
          M->GetEntry(row, col, temp);
          w[row] = w[row]+ alpha*temp*v[col];
        }
      }
    }

    return;
  }// Matrix-Vector Product temp=F^T*w(1:n)

  void PALMEigenSolver::  // compute w = op(v)
  MultMv(Complex* v, Complex* w)
  {

    Complex zzero(0.0,0.0), alpha;
    int             idx;
    Complex* temp = new Complex[rk]; // internal working space

    for(int i=0;i<n+rk*p;i++) w[i] = zzero; // initialize w to all zero;

    /* A) compute w(1:n) = Q^-1 (-sigma^2 M*v(1:n)-sqrt(sigma)F(Ir\cron a^T*E)v(n+1:end) ) */

    /* A.1) compute w(1:n) = -sigma^2*M*v(1:n) */

    alpha = -sigma*sigma;
    Mult(m_CRS,v,w, alpha);

    /* A.2) compute temp = (I_r\cron a^T*E)*v(n+1:end) */

    for(int i=0;i<rk; i++){
      idx = n+i*p;
      temp[i]=-a[0]*b[0]*v[idx];
      for(int j=1;j<p; j++){
        temp[i] = temp[i] -a[j]*b[j]*v[idx+j];
      }
    }

    /* A.3) compute w(1:n) = w(1:n) -sqrt(sigma) F * temp. */

    alpha = -sqrt(sigma/delta);
    MultM(&E,temp,w,alpha, 'N');

    /* B) compute Q^{-1} w(1:n) */

    if(!isfactor) lufactor();  // call LU factorization if not factored.

    int         info, one=1;
    SuperMatrix B;
    zCreate_Dense_Matrix(&B, n, one, (doublecomplex *)w, n, SLU_DN, SLU_Z, SLU_GE);
    StatInit(&stat);
    trans_t trans = NOTRANS;
    zgstrs(trans, &L, &U, permc, permr, &B, &stat, &info);
    Destroy_SuperMatrix_Store(&B);

    /* C) compute w(n+1:end) = (Ir\cron E)*v(n+1:end) - sqrt(sigma)(Ir\cron a)F^T w(1:n) */

    for(int i=0;i<rk;i++) temp[i]=zzero;
    alpha = 1.0;

    /* C.1) compute temp=F^T*w(1:n) */

    MultM(&F, w,temp,alpha, 'T');

    /* C.2) compute w2 = sqrt(sigma)(I_r\cron a)*temp. */

    for(int i=0;i<rk;i++){
      for(int j = 0;j<p;j++){
        w[n+i*p+j] = -b[j]*v[n+i*p+j]-sqrt(sigma*delta)*a[j]*temp[i];
      }
    }

    return;
  }

  void PALMEigenSolver::
  CalcEigenValues( BaseVector &sol, BaseVector &err, UInt N, Double shiftPoint )// EigenValue Solver
  {
    LOG_DBG(palm) << "CalcEigenValues: START";
    shared_ptr<Timer> timer = info_->Get("palm/timer")->AsTimer();
    timer->Start();

    nev = N;
    sigma = Complex(0.0,shiftPoint);

    LowRankDcomp();

    nl = n + p*rk;

    /* B) Compute Pade parameter a and b. */

    b = new double[p];
    a = new double[p];

    for(int i=0; i<p; i++){
      b[ i ] = cos( (i + 1) * M_PI / (2 * p + 1) );
      a[ i ] = sin( (i + 1) * M_PI / (2 * p + 1) );
      a[ i ] = a[ i ] / b[ i ] * sqrt(2) / sqrt(2 * p + 1);
      b[ i ] = b[ i ] * b[ i ];
    }
    /* PARAMETER AND WORKSPACE SETTING UP */
    LOG_DBG(palm) << "CalcEigenValues: Mark 1";
    int ido = 0;  // First call to ARPACK
    const char bmat = 'I'; // Bmat
    std::string which("LM"); // Part of spectral

    int ncv = 2*nev + 1; // (default=2*nev+1) Number of Arnoldi vectors at each iteration.
    int lworkl = ncv * ( 3 * ncv + 6 );
    // Dimension of internal vectors
    int lworkv = 2 * ncv; // Dimension of internal vectors
    int lrwork = ncv;  // Dimension of internal vectors
    Complex *V = new Complex[ nl * ncv + 1 ];
    // Arnoldi vectors
    Complex *resid = new Complex[ nl ];
    // Working space
    Complex *workd = new Complex[ 3 * nl + 1 ];
    // Working space
    Complex *workl = new Complex[ lworkl + 1 ];
    // Working space
    Complex *workv = new Complex [lworkv + 1 ];
    // Working space
    double  *rwork = new double[ lrwork + 1];
    // Working space
    LOG_DBG(palm) << "CalcEigenValues: Mark 2";
    int  iparam[ 12 ] = {0}; // Vector that handles ARPACK parameters.
    int  ipntr[ 15 ]; // Vector than handles ARPACK pointers.

    nconv = 0;  // No eigenvalues found yet
    iparam[ 1 ] = 1;  // Using AutoShift: CAN USER-DEFINED
    iparam[ 3 ] = 300; // Maximum number of IRAM
    iparam[ 4 ] = 1;  // Block size must be 1.
    iparam[ 7 ] = 1;  // mode bmat=1.

    int info = 0; // Use random starting vectors if = 0,
    // otherwise user-defined.

    scaling();   // Call scaling function,
    // it will do nothing if isscaling = false.

    if ( !isfactor ) lufactor(); // LU factorization, if LU is not done.

    bool BasisOK = false;

    /* MAIN ARPACK LOOP: REVERSE COMMUNICATION */

    while ( !BasisOK ) {

      /* Call ARPACK routine for updating Arnoldi basis */

      znaupd( &ido, &bmat, &nl, which.c_str(), &nev, &artol, resid, &ncv,
          &V[1], &nl, &iparam[1], &ipntr[1], &workd[1],
          &workl[1], &lworkl, &rwork[1], &info, sizeof(bmat), sizeof(which.c_str()));

      if ( ido == 99 ) {
        nconv = iparam[ 5 ];
        if ( info >= 0 ) BasisOK = true;
      }

      if (( ido == -1 ) || ( ido == 1 )) {
        // Performing Matrix vector multiplication: y<- op*x.
        MultMv( &workd[ ipntr[1] ], &workd[ ipntr[2] ]);
      }

    }
    LOG_DBG(palm) << "CalcEigenValues: Mark 3";
    /* Computing Ritz values */

    EigVal = new Complex[ nev + 1 ];  // Eigenvalues
    EigVec = new Complex[ nl*nev + 1 ]; // Eigenvectors

    bool rvec = (iseigv) ? 1 : 0; // 0 not to compute eigenvectors, 1 compute eigenvectors.

    char HowMny = 'A';   // 'A' for Ritz vectors, 'P' for Shur vectors.

    bool *iselect = new bool[ ncv ]; // Internal working space.
    LOG_DBG(palm) << "CalcEigenValues: Mark 4";
    /* Call ARPACK routine for the Ritz values and/or Ritz vectors. */
    zneupd( &rvec, &HowMny, iselect, EigVal, EigVec, &nl, &sigma, &workv[1],
        &bmat, &nl, which.c_str(), &nev, &artol, resid, &ncv,
        &V[1], &nl, &iparam[1], &ipntr[1], &workd[1], &workl[1],
        &lworkl, &rwork[1], &info);

    LOG_DBG(palm) << "CalcEigenValues: Mark 5";
    /* Convert LEP eigenvalues back to QEP.  */
    for( int i = 0; i < nev; i++) {
      LOG_DBG(palm) << i << ":" << EigVal[ i ];
      EigVal[ i ] = sigma * sqrt( 1.0 / EigVal[ i ] + 1.0 );
    }
    //LOG_DBG(palm) << "CalcEigenValues: Mark 6";
    /* Release storage.    */
    delete[] resid;
    delete[] workd;
    delete[] workl;
    delete[] workv;
    delete[] rwork;

    /* Return if only eigenvalues are required.  */
    if(!iseigv) return;

    /* Return if relative residual norms are not required */
    if(!isresid) return;

    ResidComp();  // Relative residual norms computation

    sol.Resize(nev);
    err.Resize(nev);

    for(int i=0; i<nev; i++)
    {
      sol.SetEntry(i, EigVal[i]);
      err.SetEntry(i, (Complex)res[i]);
    }
    timer->Stop();
    LOG_DBG(palm) << "CalcEigenValues: END";
  }

  void PALMEigenSolver::
  ResidComp()
  {
    if(resdone) return;    // Quick return.


    res = new double[ nev ];

    /* Internal variables.    */
    Complex *temp = new Complex[ n ];
    Complex ev;
    Complex alpha = 1.0;
    int te = 1;
    double   relerr;

    /* Compute matrix 1-norm if necessary.  */
    normm = ( normm==0 )? Norm(m_CRS) : normm;
    normd = ( normd==0 )? Norm(c_CRS) : normd;
    normk = ( normk==0 )? Norm(k_CRS) : normk;

    /* Calculate relative residual norms.  */

    for ( int i=0; i<nev; i++ ) {
      ev = EigVal[i];            // Eigenvalue i.

      Mult(m_CRS, &EigVec[i*nl], temp, ev*ev); // Compute residual Q(e_i)x_i.
      Mult(c_CRS, &EigVec[i*nl], temp, ev);
      Mult(k_CRS, &EigVec[i*nl], temp, alpha);

      relerr = dznrm2_(&n, temp, &te); // Residual 2-norm r = |Q(e_i)x_i|_2

      res[i] = relerr/(ABS(ev)*ABS(ev)*normm+ABS(ev)*normd+normk);
      // Relative residual norm
    }

    resdone = true;  // Relative residual computed.

  }

  void PALMEigenSolver::GetEigenMode(unsigned int modeNr, Vector<Complex>& mode, bool right)
  {
    mode.Resize(n);
    for(int i=0; i<n; i++)
    {
      mode.SetEntry(i, EigVec[modeNr*nl+i]);
    }
  }
}


