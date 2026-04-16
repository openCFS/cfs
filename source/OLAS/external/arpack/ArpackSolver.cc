#include <limits>
#include <boost/type_traits/is_complex.hpp>

#include "ArpackSolver.hh"
#include "General/Environment.hh"
#include "General/Exception.hh"
#include "MatVec/Vector.hh"
#include "Utils/ToolsFull.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/LogConfigurator.hh"


DEFINE_LOG(as, "arpackSolver")


namespace CoupledField {

  ArpackSolver::ArpackSolver(PtrParamNode xml) :
    eigenValues_(NULL),
    eigenVectors_(NULL)
  {
    computeMode_ = ArpackMatInterface::ComputeMode::SHIFT_INVERT;
    tolerance_ = xml->Has("tolerance") ? xml->Get("tolerance")->As<double>() : 1.e-8;
    maxIterations_ = xml->Has("maxIt") ? xml->Get("maxIt")->As<unsigned int>() : 5000;
    arnoldiFactor_ = xml->Has("arnoldiFactor") ? xml->Get("arnoldiFactor")->As<double>() : 2.0;
    numArnoldiVec_ = 0; // set later, based on arnoldiFactor_
    // cast to char* or receive compiler warning
    xml_which_ = xml->Has("which") ? xml->Get("which")->As<std::string>() : "";

    // overwritten via ArpackEigenSolver ?? why the hell largest magnitude??
    which_ = xml_which_.empty() ? (char*) "LM" : const_cast<char*>(xml_which_.c_str());

    type_ = (char*) "G";
    size_ = 0;
    interface_ = NULL;

    logging_ = xml->Has("logging") ? xml->Get("logging")->As<bool>() : false;
    if(logging_)
      DebugOn();
    else
      DebugOff();

    counter_calll_aupd = 0;
    counter_solve_OP_x = 0;
    counter_solve_OP_B_x = 0;
    counter_B_x = 0;
  }
  
  ArpackSolver::~ArpackSolver() {

    delete eigenValues_;
    delete eigenVectors_;
    interface_ = NULL;

  }

  void ArpackSolver::Setup( ArpackMatInterface *matInterface, UInt size, char* which, char* type,
                            ArpackMatInterface::ComputeMode computeMode, bool complex ) {

    size_ = size;
    computeMode_ = computeMode;

    // new reference to matrix interface
    interface_ = matInterface;

    // Save frequency parameters
    which_ = which;
    type_ = type;

    // set default value for Arnoldi vectors
    numArnoldiVec_ = std::max(int(size_ + 1), int(size_ * arnoldiFactor_));

    // check, if number of Arnoldi vectors is larger than
    // size of system
    // if( numArnoldiVec_ > size_) {
    //   numFreq_ = size / 2;
    //   numArnoldiVec_ = numFreq_*2;
    //   WARN( "Number of Arnoldi vectors will be re-set to " << numFreq_ << " as the number of eigenfrequencies may only be 1/2 the number of unknowns of the system");
    // }
    
    // Setup() is called for each wave_vector
    if(eigenValues_ == NULL || eigenVectors_ == NULL)
    {
      // eigenvalues, tolerances and vectors
      if(complex) {
        eigenValues_  = new Vector<Complex>();
        eigenVectors_ = new Vector<Complex>();
      } else {
        eigenValues_  = new Vector<double>();
        eigenVectors_ = new Vector<double>();
      }
    }
  }


  void ArpackSolver::QuadSetup( ArpackMatInterface *matInterface, UInt size, char* which,
                                ArpackMatInterface::ComputeMode computeMode ) {

    size_ = size;
    computeMode_ = computeMode;

    // new reference to matrix interface
    interface_ = matInterface;

    // Save frequency parameters
    which_ = which;

    // adjust defaults for number of iterations and tolerance
    tolerance_     = 1e-10;
    maxIterations_ = 10000;

    // set default value for Arnoldi vectors
    numArnoldiVec_ = size_*2;

    // eigenvalues, tolerances and vectors
    eigenValues_  = new Vector<Complex>();
    eigenVectors_ = new Vector<Complex>();

    DebugOff();

  }


  void ArpackSolver::ToInfo(PtrParamNode info)
  {
    info->Get("frequencies")->SetValue(numEV_);
    info->Get("type")->SetValue(which_);
    info->Get("shift")->SetValue(valueShift_);
    if (computeMode_ == ArpackMatInterface::ComputeMode::SHIFT_INVERT) {
      info->Get("computeMode")->SetValue("shift_and_invert");
    }
    else if (computeMode_ == ArpackMatInterface::BUCKLING) {
      info->Get("computeMode")->SetValue("buckling");
    }
    else {
      info->Get("computeMode")->SetValue("regular");
    }
    info->Get("tol")->SetValue(tolerance_);
    info->Get("maxIter")->SetValue(maxIterations_);
    info->Get("arnoldiVectors")->SetValue(numArnoldiVec_);
    info->Get("scaling")->SetValue(interface_->GetDiagScaling());
    info->Get("logging")->SetValue(logging_);
  }

  template <class TYPE>
  UInt ArpackSolver::FindEigenvalues(UInt numEV, TYPE valueShift)
  {
    // based on zndrv4.f and dsdrv4.f
    numEV_ = numEV;
    valueShift_ = (Complex) valueShift;

    // set default value for Arnoldi vectors
    numArnoldiVec_ = std::max(int(numEV_ + 2), int(numEV_ * arnoldiFactor_));

    // check, if number of Arnoldi vectors is larger than size of system
    if( numArnoldiVec_ > size_) {
      UInt newNumEV = static_cast<UInt>( size_ / 2.0 );
      numEV_ = newNumEV;
      numArnoldiVec_ = size_;
      std::stringstream out;
      WARN( "Number of eigenvalues (" << numEV << ") has been reset to "
            << newNumEV << " as the number of eigenvalues may only be"
            << " 1/2 the number of unknowns of the system" );
    }

    bool complex = boost::is_complex<TYPE>::value;

    // eigenvalues, tolerances and vectors
    eigenValues_->Resize(numArnoldiVec_);
    eigenVectors_->Resize(numEV_*size_);
    eigenTolerances_.Resize(numArnoldiVec_);

    bool converged = false;
    int ido = 0; // necessary initial value for first aupd() call
    int info = 0;

    // temp vector to store B*x
    StdVector<TYPE> tempV(size_); // TYPE* tempV = new Double[size_];
    // pointers to the position of x and y in workD
    TYPE* vecX, *vecY;

    // temp working space required in *aupd
    StdVector<TYPE> workD(3*size_);
    StdVector<TYPE> residual(size_);
    // int lenWorkL = (complex ? numArnoldiVec_ * (3 * numArnoldiVec_ + 5): numArnoldiVec_ * (numArnoldiVec_ + 8));
    int lenWorkL = (complex ? numArnoldiVec_ * (10 * numArnoldiVec_ + 20): numArnoldiVec_ * (numArnoldiVec_ + 8));
    StdVector<TYPE> workL(lenWorkL);
    StdVector<TYPE> matrixV(size_*numArnoldiVec_);

    // this double array is only for the complex part
    StdVector<double> rwork(size_);


    StdVector<int> iparams(21); // initialized with 0's
    StdVector<int> ipntr(21);

    iparams[0] = 1; // ISHIFT = 1:
    iparams[2] = maxIterations_; // NXITER
    iparams[3] = 1; // NB blocksize to be used in the recurrence
    // ARPACK manual, chap 3.5
    if (computeMode_ == ArpackMatInterface::BUCKLING) {
      assert(!complex);
      iparams[6] = 4;
    } else {
      iparams[6] = 3;
    }

    counter_calll_aupd = 0;
    counter_solve_OP_x = 0;
    counter_solve_OP_B_x = 0;
    counter_B_x = 0;

    // Reverse Communication Interface
    // https://www.caam.rice.edu/software/ARPACK/UG/node80.html#SECTION001013000000000000000
    UInt itNum;
    for (itNum=0; itNum<maxIterations_; itNum++) {
      // additionally rwork is ignored for the real case
      CallAUPD(&ido, type_, &size_, which_, &numEV_, &tolerance_, residual.GetPointer(),
            &numArnoldiVec_, matrixV.GetPointer(), &size_, iparams.GetPointer(), ipntr.GetPointer(),
            workD.GetPointer(), workL.GetPointer(), &lenWorkL, rwork.GetPointer(), &info);

      counter_calll_aupd++;

      switch (ido)
      {
      case -1:
        // NOTE: addresses in workd are given as FORTRAN addreses
        // solve (A-shift*B)*y = w, where w = A*x (buckling) or w = B*x (else).
        // x is returned by dsaupd in workd(ipntr[1]:ipntr[1]+size-1)
        // y is expected by dsaupd in workd(ipntr[2]:ipntr[2]+size-1)

        vecX = workD.GetPointer() + (ipntr[0] -1);
        vecY = workD.GetPointer() + (ipntr[1] -1);

        if (computeMode_ == ArpackMatInterface::BUCKLING) {
          // w = A*x
          interface_->MultAV(vecX, tempV.GetPointer());
        } else {
          // w = B*x
          interface_->MultBV(vecX, tempV.GetPointer());
        }
        // solve (A-shift*B)*y = w
        interface_->MultShiftOpV(tempV.GetPointer(), vecY);
        counter_solve_OP_x++;
        break;

      case 1:
        // solve (A-shift*B)*y = x, where x = A*z (buckling) or x = B*z (else).
        // x is returned by dsaupd in workd(ipntr[3]:ipntr[3]+size-1),
        //   i.e. action by A or B has already been made.
        // y is expected by dsaupd in workd(ipntr[2]:ipntr(2]+size-1)
        vecX = workD.GetPointer() + (ipntr[2]-1);
        vecY = workD.GetPointer() + (ipntr[1]-1);
        interface_->MultShiftOpV(vecX,vecY);
        counter_solve_OP_B_x++;
        break;

      case 2:
        // x is returned by dsaupd in workd(ipntr[1]:ipntr[1]+size-1)
        // y is expected by dsaupd in workd(ipntr[2]:ipntr[2]+size-1)
        vecX = workD.GetPointer() + (ipntr[0]-1);
        vecY = workD.GetPointer() + (ipntr[1]-1);
        if (computeMode_ == ArpackMatInterface::BUCKLING) {
          // calculate y = A*x
          interface_->MultAV(vecX,vecY);
        } else {
          // calculate y = B*x
          interface_->MultBV(vecX,vecY);
        }
        counter_B_x++;
        break;

      case 3:
        // algorithm requires shifts - > not implemented
        EXCEPTION( "User required shifts not yet implemented!\nARPACK requested ido=3 in DSAUPD!" );
        break;

      case 99:
        converged = true;
        break;

      default:
        // something went wrong -> arpack requires undefind ido
        EXCEPTION( "ARPACK requested not implemented ido=" << ido << " in DSAUPD!" );
        break;
      }

      LOG_DBG3(as) << "FE: complex=" << complex << " it=" << itNum << " ido=" << ido << " info=" << info;
      LOG_DBG3(as) << " x=" << ToString<TYPE>(vecX, size_);
      LOG_DBG3(as) << " y=" << ToString<TYPE>(vecY, size_);

      // if convergence has been achieved, exit loop
      if (converged)
        break;

    }

    // check whether everything went well
    if (info<0)
      EXCEPTION("Error reported in Ritz value calculation: info=" << info << " iparam=" << iparams.ToString() << "\n"  << ArpackError(info) );


    if ( itNum>maxIterations_ && info != 99 ) {
      EXCEPTION("Maximum number of iterations has been reached, but no solution"
          << "could be obtained which satisfies requested tolerance setting.\n"
          << "Increase maxIt (=" << maxIterations_ << ") in xml file!" );
    }

    // proceed with calculation of eigenvectors, values, and tolerances

    // NOTE: in the FORTRAN routine dseupd d and eigenVector_ array are
    //       treated as 2-dimensional arrays of dimension
    //       (numFreq_,2) and (size_,numArnoldiVec_), respectively.
    //       However, taking account of FORTRAN storage scheme for arrays
    //       (column first) in data transfer, we use a simple 1d vector
    //       memory allocation here!

    bool rvec = true;
    StdVector<double> select(numArnoldiVec_); // is logical in Fortran! // Double *select = new Double [numArnoldiVec_];
    StdVector<TYPE> d(numArnoldiVec_*2);  // in dsdrv4.f (maxncv,2) and in zndrv4.f (maxncv)
    TYPE vShift = valueShift;

    //eigenValues_->Resize(numArnoldiVec_);
    //eigenVectors_->Resize(numEV_*size_);
    //eigenTolerances_.Resize(numArnoldiVec_);
    Vector<TYPE>& eval = dynamic_cast<Vector<TYPE>&>(*eigenValues_);
    Vector<TYPE>& evec = dynamic_cast<Vector<TYPE>&>(*eigenVectors_);

    // only complex part
    StdVector<TYPE> workev(2 * numArnoldiVec_);

    // cast to char* or receive compiler warning
    CallEUPD(&rvec, (char*) "A", select.GetPointer(), d.GetPointer(), evec.GetPointer(), &size_,
        &vShift, workev.GetPointer() , (char*) "G", &size_, which_, &numEV_,
        &tolerance_, residual.GetPointer(), &numArnoldiVec_, matrixV.GetPointer(),
        &size_, iparams.GetPointer(), ipntr.GetPointer(), workD.GetPointer(),
        workL.GetPointer(), &lenWorkL, rwork.GetPointer(), &info);

    LOG_DBG3(as) << "FE post d=" << d.ToString();

    if(info != 0)
      EXCEPTION("error postprocessing arpack: info=" << info);

    Integer numEVConverged = iparams[4];
    for (Integer i=0; i< numEVConverged; i++)
    {
      eval[i] = d[i];                 // cf. note above!
      // calculate error norm of obtained solution
      vecX = matrixV.GetPointer() + i*size_;
      vecY = workD.GetPointer() + size_;
      interface_->MultBV(vecX, vecY);
      vecY = workD.GetPointer();
      interface_->MultAV(vecX, vecY);
      double vNrm2 = 0.0;
      for (Integer j=0; j<size_; j++) {
        vecX[j] = vecY[j] - d[i]*vecY[j+size_];
        vNrm2 += abs(vecX[j])*abs(vecX[j]);
      }
      eigenTolerances_[i] = sqrt(vNrm2)/abs(d[i]);    // cf. note above!
    }

    return numEVConverged;
  }


  /**
  * Call of ARPACK subroutine.
  * The first letter of a subroutine denotes the data type:
  * s  single precision real arithmetic,
  * d   double precision real arithmetic,
  * c   single precision complex arithmetic,
  * z   double precision complex arithmetic.
  * The second one denotes the type of eigensystem:
  * n   non-symmetric,
  * s   symmetric.
  */
  template <>
  void ArpackSolver::CallAUPD<double>(Integer* ido, char* bmat, Integer* n, char* which, Integer* nev, Double* tol,
      Double *resid, Integer *ncv, Double *V, Integer *ldv, Integer *iparam, Integer *ipntr,
      Double *workd, Double *workl, Integer *lworkl, Double *workDbleD, Integer *info)
  {
    /**
    * This method has been designed to compute approximations to a
    * few eigenpairs of a linear operator OP that is real and symmetric
    * with respect to a real positive semi-definite symmetric matrix B.
    * https://www.caam.rice.edu/software/ARPACK/UG/node136.html
    */
    dsaupd(ido, bmat, n, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, info);
  }

  // @see CallAUPD<double>()
  template <>
  void ArpackSolver::CallAUPD<Complex>(Integer* ido, char* bmat, Integer* n, char* which, Integer* nev, Double* tol,
      Complex *resid, Integer *ncv, Complex *V, Integer *ldv, Integer *iparam, Integer *ipntr,
      Complex *workd, Complex *workl, Integer *lworkl, Double *workDbleD, Integer *info)
  {
    /**
    * This is intended to be used to find a few eigenpairs of a
    * complex linear operator OP with respect to a semi-inner product defined
    * by a hermitian positive semi-definite real matrix B.
    * https://www.caam.rice.edu/software/ARPACK/UG/node138.html
    */
    znaupd(ido, bmat, n, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, workDbleD, info);
  }

  // @see CallAUPD<double>()
  template <>
  void ArpackSolver::CallEUPD<double>(bool *rvec, char *howMny, Double *select, double *d, double *z, Integer *ldz,
      double *shift, double *zwork, char *bmat, Integer* size, char *which, Integer *nev, Double *tol,
      double *resid, Integer *ncv, double *V, Integer *ldv, Integer *iparam, Integer *ipntr,
      double *workd, double *workl, Integer *lworkl, Double *workDbleD, Integer *info)
  {
    /**
    * On the final return from dsaupd (indicated by ido = 99), the error flag info must be checked.
    * If info = 0 then no fatal errors have occurred and it is time to post-process using dseupd
    * to get eigenvalues of the original problem and the corresponding eigenvectors if desired.
    * https://www.caam.rice.edu/software/ARPACK/UG/node40.html
    */
    dseupd(rvec, howMny, select, d, z, ldz, shift, bmat, size, which, nev, tol,
        resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, info);
  }

  // @see CallAUPD<double>()
  template <>
  void ArpackSolver::CallEUPD<Complex>(bool *rvec, char *howMny, Double *select, Complex *d, Complex *z, Integer *ldz,
      Complex *shift, Complex *zwork, char *bmat, Integer* size, char *which, Integer *nev, Double *tol,
      Complex *resid, Integer *ncv, Complex *V, Integer *ldv, Integer *iparam, Integer *ipntr,
      Complex *workd, Complex *workl, Integer *lworkl, Double *workDbleD, Integer *info)
  {
    /**
    * On the final return from znaupd (indicated by ido = 99), the error flag info must be checked.
    * If info = 0 then no fatal errors have occurred and it is time to post-process using zneupd
    * to get eigenvalues of the original problem and the corresponding eigenvectors if desired.
    * https://www.caam.rice.edu/software/ARPACK/UG/node44.html
    */
    zneupd(rvec, howMny, select, d, z, ldz, shift, zwork, bmat, size, which, nev, tol,
        resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, workDbleD, info);
  }

  template <class TYPE>
  UInt ArpackSolver::FindQuadEigenvalues(UInt numEV, TYPE valueShift)
  {
    numEV_ = numEV;
    valueShift_ =(Complex) valueShift;

    // adjust to higher number of wanted eigenvalues
    Double logNF = std::log10(numEV_*1.0);
    if (logNF > 1.0) {
      Integer iterInc = static_cast<Integer>( 500 * std::floor(5*pow(10,logNF-1)) );
      maxIterations_ += iterInc;
    }

    // set default value for Arnoldi vectors
    numArnoldiVec_ = std::max(int(numEV_ + 2), int(numEV_ * arnoldiFactor_));

    // check, if number of arnoldi vectors is larger than size of system
    if( numArnoldiVec_ > size_ ) {
      UInt newNumFreq = static_cast<UInt>( size_ / 2.0 );
      numEV_ = newNumFreq;
      numArnoldiVec_ = size_;
      std::stringstream out;
      WARN( "Number of eigenvalues (" << numEV << ") has been reset to "
            << newNumFreq << " as the number of eigenvalues may only be"
            << " 1/2 the number of unknowns of the system" );
    }

    eigenValues_->Resize(numArnoldiVec_);
    eigenVectors_->Resize(numEV_*size_);
    eigenTolerances_.Resize(numArnoldiVec_);
    Vector<Complex>& eval = dynamic_cast<Vector<Complex>& >(*eigenValues_);
    Vector<Complex>& evec = dynamic_cast<Vector<Complex>& >(*eigenVectors_);

    bool converged = false;
    Integer  ido = 0, info = 0;

    // temp vector to store B*x
    Complex *tempV = new Complex [size_];
    // pointers to the position of x and y in workD
    Complex *vecX, *vecY;

    // temp working space required in znaupd
    Complex *workD = new Complex [3*size_];
    Complex *residual = new Complex [size_];
    Integer lenWorkL = numArnoldiVec_*(3*numArnoldiVec_+5);
    Complex *workL = new Complex [lenWorkL];
    Complex *matrixV = new Complex [size_*numArnoldiVec_];
    Double  *workDbleD = new Double [numArnoldiVec_];

    InitQuadTempSpace(tempV, residual, workD, workL, matrixV, workDbleD);

    Integer iparams[21], ipntr[21];
    for (Integer i=0; i<21; i++) {
         iparams[i] = 0;
         ipntr[i] = 0;
    }

    iparams[0] = 1;
    iparams[2] = maxIterations_;
    iparams[6] = 3;

    UInt itNum;
    for (itNum=1; itNum<=maxIterations_; itNum++) {

      znaupd(&ido, (char*) "G", &size_, which_, &numEV_, &tolerance_, residual,
                      &numArnoldiVec_, matrixV, &size_, iparams, ipntr, workD, workL,
                      &lenWorkL, workDbleD, &info);

        switch (ido)  {

        case -1:
            // NOTE: addresses in workd are given as FORTRAN addreses
            // solve (A+shift*D+shift**2*B)*y = B*x
            // x is returned by znaupd in workd(ipntr[1]:ipntr[1]+size-1)
            // y is expected by znaupd in workd(ipntr[2]:ipntr[2]+size-1)

           vecX = workD + --ipntr[0];
           vecY = workD + --ipntr[1];
           interface_->MultBVQuad(vecX, tempV);
           interface_->MultShiftOpVQuad(tempV, vecY);
           break;

        case 1:
            // solve (A+shift*D+shift**2*B)*y = B*x
            // B*x is returned by znaupd in workd(ipntr[3]:ipntr[3]+size-1)
            // y is expected by znaupd in workd(ipntr[2]:ipntr(2]+size-1)
           vecX = workD + --ipntr[2];
           vecY = workD + --ipntr[1];
           interface_->MultShiftOpVQuad(vecX,vecY);
           break;

        case 2:
            // calculate y = B*x
            // x is returned by znaupd in workd(ipntr[1]:ipntr[1]+size-1)
            // y is expected by znaupd in workd(ipntr[2]:ipntr[2]+size-1)
            vecX = workD + --ipntr[0];
            vecY = workD + --ipntr[1];
            interface_->MultBVQuad(vecX,vecY);
            break;

        case 3:
            // algorithm requires shifts - > not implemented
            EXCEPTION("User required shifts not yet implemented!\n"
                     << "ARPACK requested ido=3 in ZNAUPD!");
            break;

        case 99:
            converged = true;
            break;

        default:
            // something went wrong -> arpack requires undefind ido
            EXCEPTION("ARPACK requested not implemented ido=" << ido
                     << " in ZNAUPD!");
            break;

        }

        LOG_DBG3(as) << "FQE: it=" << itNum << " ido=" << ido << " info=" << info;
        LOG_DBG3(as) << " x=" << StdVector<Complex>::ToString(size_, vecX);
        LOG_DBG3(as) << " y=" << StdVector<Complex>::ToString(size_, vecY);

        // if convergence has been achieved, exit loop
        if (converged)
            break;


    }

    // check whether everything went well
    if (info != 0) {
      EXCEPTION("IPARAM(5) = " << iparams[4] <<
      ". Error reported in ritz value calculation:\n" << ArpackError(info) );
    }

    if ( itNum>maxIterations_ && info != 99 ) {
      EXCEPTION("Maximum number of iterations has been reached, but no solution"
          << "could be obtained which satisfies requested tolerance setting.\n"
          << "Increase maxIt (=" << maxIterations_ << ") in xml file!");
    }

    // proceed with calculation of eigenvectors, values, and tolerances

    // NOTE: in the FORTRAN routine zneupd d and eigenVector_ array are
    //       treated as 2-dimensional arrays of dimension
    //       (numFreq_,2) and (size_,numArnoldiVec_), respectively.
    //       However, taking account of FORTRAN storage scheme for arrays
    //       (column first) in data transfer, we use a simple 1d vector
    //       memory allocation here!

    bool rvec = true;
    Double *select = new Double [numArnoldiVec_];
    Complex *d = new Complex [numArnoldiVec_*2];
    Complex *zwv = new Complex [numArnoldiVec_*2];

    // we define a second tmp-array to store all eigenvectors,
    // also for the negative frequencies.
    // zEigenVectors_ will only store data for positive freqs.
    Complex *zEV = new Complex [numEV_*size_];

    for (int i=0; i< numArnoldiVec_*2; i++) {
      d[i] = (Complex) 0.0;
      zwv[i] = (Complex) 0.0;
    }
    for (int i=0; i< numEV_*size_; i++) {
      zEV[i]  = (Complex) 0.0;
      evec[i] = (Complex) 0.0;
    }

    Complex cShift = (Complex) valueShift_;

    zneupd(&rvec, (char*) "A", select, d, zEV, &size_, &cShift, zwv, (char*) "G",
                    &size_, which_, &numEV_, &tolerance_, residual, &numArnoldiVec_,
                    matrixV, &size_, iparams, ipntr, workD,
                    workL, &lenWorkL, workDbleD, &info);

    LOG_DBG3(as) << "FQE post d=" << ToString(d, numArnoldiVec_*2);

    Integer numEVConverged = iparams[4];

    Complex *z1 = new Complex [size_/2];
    Complex *z2 = new Complex [size_/2];
    Complex *z3 = new Complex [size_/2];
    Complex *zA = new Complex [size_];
    Complex *zB = new Complex [size_];
    Complex *zC = new Complex [size_];
    Complex *ztmp;

    // we only count and treat/store eigenvalues with positive frequency value
    Integer *sortedEV = new Integer [numEVConverged];
    Double  *sortedFreq = new Double [numEVConverged];
    for (Integer ncv=0; ncv < numEVConverged; ncv++) {
      sortedEV[ncv] = -1;
      sortedFreq[ncv] = 0.0;
    }

    Integer found = 0, pos;
    for (Integer ncv=0; ncv < numEVConverged; ncv++) {
      Double freq = d[ncv].imag();
      if (freq >= 0.0) {
        if (found==0) {
           sortedFreq[found] = freq;
           sortedEV[found] = ncv;
        } else {
          pos = found;
          for (Integer k=found-1; k >=0; k--) {
            if (sortedFreq[k] > freq)
               pos--;
          }
          for (Integer k=found-1; k >=pos; k--) {
            sortedFreq[k+1] = sortedFreq[k];
            sortedEV[k+1] = sortedEV[k];
          }
          sortedFreq[pos] = freq;
          sortedEV[pos] = ncv;
        }
        found++;
      }
    }

    Integer *evPos = new Integer [numEVConverged];
    for (Integer pos=0; pos<numEVConverged; pos++) {
      evPos[pos] = -1;
    }
    for (Integer pos=0; pos<found; pos++) {
      if (sortedEV[pos] >= 0) {
        evPos[sortedEV[pos]] = pos;
      }
    }

    for (Integer ncv=0; ncv < numEVConverged; ncv++) {

      if (d[ncv].imag() >= 0.0) {

        pos = evPos[ncv];
        if (pos >= 0) {
          eval[pos] = d[ncv];

          LOG_DBG3(as) << "CQV post ncv=" << ncv << " pos=" << pos << " eval=" << eval[pos];

          for (int j=0; j<size_; j++) {
            evec[pos*size_+j] = zEV[ncv*size_+j];
          }

          ztmp = evec.GetPointer() + pos*size_;
          for (int j=0; j<size_/2; j++) {
            z3[j] = ztmp[j];
          }
          interface_->MultZDampV(ztmp, z1);
          ztmp += size_/2;
          interface_->MultZStiffV(ztmp, z2);
          for (int j=0; j<size_/2; j++) {
            z1[j] = -z1[j] - z2[j];
          }

          if ( interface_->UseStiffInNMat() ) {
            for (int j=0; j<size_/2; j++) {
              z2[j] = z3[j];
            }
            interface_->MultZStiffV(z2, z3);
          }
          interface_->ScaleDiag(z3);
          for (int j=0; j<size_/2; j++) {
              zA[j]         = z1[j];
              zA[j+size_/2] = z3[j];
          }

          ztmp = evec.GetPointer() + pos*size_;
          for (int j=0; j<size_/2; j++) {
             z3[j] = ztmp[j+size_/2];
          }
          interface_->MultZMassV(ztmp, z1);

          if ( interface_->UseStiffInNMat() ) {
            for (int j=0; j<size_/2; j++) {
                z2[j] = z3[j];
            }
            interface_->MultZStiffV(z2, z3);
          }
          interface_->ScaleDiag(z3);
          for (int j=0; j<size_/2; j++) {
            zB[j]         = z1[j];
            zB[j+size_/2] = z3[j];
          }

          Double vNrm2 = 0.0,v2Nrm2 = 0.0;
          for (int j=0; j<size_/2; j++) {
            zC[j] = zA[j]-d[ncv]*zB[j];
            vNrm2 += pow(abs(zC[j]),2);
          }
          for (int j=size_/2; j<size_; j++) {
            zC[j] = zA[j]-d[ncv]*zB[j];
            v2Nrm2 += pow(abs(zC[j]),2);
          }
          vNrm2 += v2Nrm2;
          eigenTolerances_[pos] = sqrt(vNrm2)/abs(d[ncv]);

          // "normalize" eigenvectors:
          ztmp = evec.GetPointer() + pos*size_;
          Double  epsTest = 1.0e-6;
          // 1) for a unit vector
          Complex vLen = (Complex) 0.0;
          for (int j=size_/2; j<size_; j++) {
            vLen += ztmp[j]*std::conj(ztmp[j]);
          }
          vLen = sqrt(vLen);
          for (int j=0; j<size_; j++) {
            ztmp[j] /= vLen;
          }

          // 2) start at size_/2, i.e. use "essential part" of the eigenvector
          //    and look for first non-noisy entry
          int jNorm = size_/2, jFound=0;
          while (jFound==0 && jNorm<size_) {
            if (std::abs(ztmp[jNorm]) > epsTest)
              jFound = 1;
            else
              jNorm++;
          }

          // 3) rescale complete vector so that found entry will be scaled
          //    with equal real and imag parts. this is done to get rid of
          //    numerical noise in comparing eigenvectors
          //    NOTE: this should be modified, as soon as complex valued
          //    eigenvectors are available for the result files

          Complex zbase = Complex (1.,1.);
          Complex zNorm = ztmp[jNorm]*zbase;
          for (int j=0; j<size_; j++) {
            ztmp[j] /= zNorm;
          }

    // 4) calculated conj(z)^T * M * z
          //          ztmp = zEigenVectors_ + pos*size_ + size_/2;
          //          interface_->MultZMassV(ztmp, z1);
          //          Complex normalize = (Complex) 0.0;
          //          for (Integer j=0; j<size_/2; j++) {
          //             normalize += std::conj(ztmp[j])*z1[j];
          //          }
          //	  normalize = sqrt(normalize);
          //          ztmp = zEigenVectors_ + pos*size_;
          //          for (Integer j=0; j<size_; j++) {
          //             ztmp[j] /= normalize;
          //          }
    //-------------------------------------
        }
      }
    }

    return found;
  }

  void ArpackSolver::InitQuadTempSpace(Complex *tempV,Complex *residual,
         Complex* workD, Complex* workL, Complex* matrixV, Double* workDbleD) {

      int i;
      for (i=0; i < size_;i++) {
           tempV[i] = (Complex) 0.0;
           residual[i] = (Complex) 0.0;
      }

      for (i=0; i < 3*size_;i++)
           workD[i] = (Complex) 0.0;

      int len = numArnoldiVec_*(3*numArnoldiVec_+5);
      for (i=0; i < len; i++)
          workL[i] = (Complex) 0.0;

      len = size_*numArnoldiVec_;
      for (i=0; i < len; i++)
           matrixV[i] = (Complex) 0.0;

      for (i=0; i < numArnoldiVec_; i++)
           workDbleD[i] = 0.0;
  }

  Double ArpackSolver::Eigenvalue(UInt i) {
      Double v;
      eigenValues_->GetEntry(i, v);
      return v;
  }

  Complex ArpackSolver::CmplxEigenvalue(UInt i) {
      Complex v;
      eigenValues_->GetEntry(i, v);
      return v;
  }

  Double ArpackSolver::Tolerance(UInt i) {
      return eigenTolerances_[i];
  }

  Double* ArpackSolver::GetEigenvector(UInt modeNr) {
      return dynamic_cast<Vector<Double>&>(*eigenVectors_).GetPointer() + modeNr*size_;
  }
  
  Complex* ArpackSolver::GetComplexEigenvector(UInt modeNr) {
    return dynamic_cast<Vector<Complex>&>(*eigenVectors_).GetPointer() + modeNr*size_;
  }
  
  std::string ArpackSolver::ArpackError(Integer errNo) {
    std::string msg;

    switch (errNo ) {
      case 0:
        msg = "No error";
        break;

      case 1: 
        msg = "Maximum number of iterations taken. "
          " All possible eigenvalues of OP has been found. IPARAM(5) "
          " returns the number of wanted converged Ritz values";
        break;

      case 2: 
        msg = "No longer an informational error. Deprecated starting"
          "with release 2 of ARPACK.";
        break;

      case 3: 
        msg = "No shifts could be applied during a cycle of the"
          "Implicitly restarted Arnoldi iteration. One possibility"
          "is to increase the size of NCV relative to NEV";
        break;

      case  -1: 
        msg = "N must be positive";
        break;

      case -2: 
        msg = "NEV must be positive";
        break;

      case -3: 
        msg = "NCV-NEV >= 2 and less than or equal to N";
        break;
      
      case -4: 
        msg = "The maximum number of Arnoldi update iterations allowed"
        "must be greater than zero.";
        break;
      
      case -5: 
        msg = "WHICH must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'";
        break;
      
      case -6: 
        msg = "BMAT must be one of 'I' or 'G'";
        break;
      
      case -7: 
        msg = "Length of private work array WORKL is not sufficient.";
        break;
      
      case -8: 
        msg = "Error return from trid. eigenvalue calculation;" 
          " Informatinal error from LAPACK routine dsteqr.";
        break;
      
      case -9: 
        msg = "Starting vector is zero.";
        break;
      
      case -10: 
        msg = "IPARAM(7) must be 1,2,3,4,5";
        break;
      
      case -11: 
        msg = "IPARAM(7) = 1 and BMAT = 'G' are incompatable.";
        break;
      
      case -12: 
        msg = "IPARAM(1) must be equal to 0 or 1.";
        break;
      
      case -13: 
        msg = "NEV and WHICH = 'BE' are incompatable.";
        break;
      
      case -9999: 
      msg = "Could not build an Arnoldi factorization."
        " IPARAM(5) returns the size of the current Arnoldi"
        " factorization. The user is advised to check that"
        " enough workspace and array storage has been allocated.";
        
        break;

      default:
        msg = "An unknown error occured";
    }
    return msg;

  }

  void ArpackSolver::DebugOff( ) {

      ::debug.logfil = 6;
      ::debug.ndigit = 0;
      ::debug.mgetv0 = 0;
      ::debug.msaupd = 0;
      ::debug.msaup2 = 0;
      ::debug.msaitr = 0;
      ::debug.mseigt = 0;
      ::debug.msapps = 0;
      ::debug.msgets = 0;
      ::debug.mseupd = 0;
      ::debug.mnaupd = 0;
      ::debug.mnaup2 = 0;
      ::debug.mnaitr = 0;
      ::debug.mneigt = 0;
      ::debug.mnapps = 0;
      ::debug.mngets = 0;
      ::debug.mneupd = 0;
      ::debug.mcaupd = 0;
      ::debug.mcaup2 = 0;
      ::debug.mcaitr = 0;
      ::debug.mceigt = 0;
      ::debug.mcapps = 0;
      ::debug.mcgets = 0;
      ::debug.mceupd = 0;

  }

  void ArpackSolver::DebugOn( ) {

      ::debug.logfil = 6;
      ::debug.ndigit = -3;
      ::debug.mgetv0 = 1;
      ::debug.msaupd = 2;
      ::debug.msaup2 = 2;
      ::debug.msaitr = 1;
      ::debug.mseigt = 1;
      ::debug.msapps = 1;
      ::debug.msgets = 1;
      ::debug.mseupd = 2;
      ::debug.mnaupd = 2;
      ::debug.mnaup2 = 2;
      ::debug.mnaitr = 1;
      ::debug.mneigt = 1;
      ::debug.mnapps = 1;
      ::debug.mngets = 1;
      ::debug.mneupd = 1;
      ::debug.mcaupd = 2;
      ::debug.mcaup2 = 2;
      ::debug.mcaitr = 1;
      ::debug.mceigt = 1;
      ::debug.mcapps = 1;
      ::debug.mcgets = 1;
      ::debug.mceupd = 2;

  }

  // Explicit template instantiation
  template UInt ArpackSolver::FindEigenvalues<Double>(UInt, Double);
  template UInt ArpackSolver::FindEigenvalues<Complex>(UInt, Complex);

  template UInt ArpackSolver::FindQuadEigenvalues<Double>(UInt, Double);
  template UInt ArpackSolver::FindQuadEigenvalues<Complex>(UInt, Complex);

} // end of namespace
