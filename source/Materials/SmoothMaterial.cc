// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "SmoothMaterial.hh"
#include "Domain/Domain.hh"
#include "Domain/CoefFunction/CoefFunctionCompound.hh"
#include "Domain/CoefFunction/CoefFunctionConst.hh"
#include "Domain/CoefFunction/CoefXpr.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Utils/ToolsFull.hh"

DEFINE_LOG(mat_smooth, "mat_smooth")

namespace CoupledField
{

  // ***********************
  //   Default Constructor
  // ***********************
  SmoothMaterial::SmoothMaterial(MathParser* mp,
                                     CoordSystem * defaultCoosy) 
  : BaseMaterial(SMOOTH, mp, defaultCoosy)
  {
    // set the allowed material parameters
    isAllowed_.insert( SMOOTH_STIFFNESS_TENSOR );
    isAllowed_.insert( SMOOTH_KMODULUS );
    isAllowed_.insert( SMOOTH_LAME_MU );
    isAllowed_.insert( SMOOTH_LAME_LAMBDA );
    isAllowed_.insert( SMOOTH_EMODULUS );
    isAllowed_.insert( SMOOTH_EMODULUS_1 );
    isAllowed_.insert( SMOOTH_EMODULUS_2 );
    isAllowed_.insert( SMOOTH_EMODULUS_3 );
    isAllowed_.insert( SMOOTH_POISSON );
    isAllowed_.insert( SMOOTH_POISSON_3 );
    isAllowed_.insert( SMOOTH_POISSON_12 );
    isAllowed_.insert( SMOOTH_POISSON_23 );
    isAllowed_.insert( SMOOTH_POISSON_13 );
    isAllowed_.insert( SMOOTH_GMODULUS );
    isAllowed_.insert( SMOOTH_GMODULUS_3 );
    isAllowed_.insert( SMOOTH_GMODULUS_23 );
    isAllowed_.insert( SMOOTH_GMODULUS_13 );
    isAllowed_.insert( SMOOTH_GMODULUS_12 );
  }

  void SmoothMaterial::Finalize() {

    // Calculation of stiffness tensor
    ComputeFullStiffTensor();
  }

  PtrCoefFct SmoothMaterial::GetScalCoefFnc(MaterialType matType,
                                              Global::ComplexPart matDataType) const
  {
    PtrCoefFct mFunct;
    CoefMap::const_iterator it = scalarCoef_.find(matType);
    if( it !=  scalarCoef_.end() ) {
      // --------------------------------------
      //  Coefficient Function already defined
      // --------------------------------------
      mFunct = it->second->GetComplexPart( matDataType );

    }
    else {
      // Conversion is available for isotropic materials only
      std::map<MaterialType, SymmetryType>::const_iterator symType =
          symmetryType_.find(SMOOTH_STIFFNESS_TENSOR);
      assert(symType != symmetryType_.end());
      if (symType->second == ISOTROPIC) {
        // First try to find all available parameters
        CoefMap::const_iterator eModIt = scalarCoef_.find(SMOOTH_EMODULUS),
                                bulkIt = scalarCoef_.find(SMOOTH_KMODULUS),
                                poissonIt = scalarCoef_.find(SMOOTH_POISSON),
                                shearIt = scalarCoef_.find(SMOOTH_GMODULUS);
        std::map<std::string, PtrCoefFct> vars;

        switch (matType) {

          case SMOOTH_KMODULUS:
            if (eModIt != scalarCoef_.end() && shearIt != scalarCoef_.end()) {
              vars["E"] = eModIt->second;
              vars["G"] = shearIt->second;
              shared_ptr< CoefFunctionCompound<Double> >
                  kFunc(new CoefFunctionCompound<Double>(mp_));
              kFunc->SetScalar("E_R*G_R/(3*(3*G_R-E_R))", vars);
              mFunct = kFunc;
            }
            else if (eModIt != scalarCoef_.end() && poissonIt != scalarCoef_.end()) {
              vars["E"] = eModIt->second;
              vars["nu"] = poissonIt->second;
              shared_ptr< CoefFunctionCompound<Double> >
                  kFunc(new CoefFunctionCompound<Double>(mp_));
              kFunc->SetScalar("E_R/(3*(1-2*nu_R))", vars);
              mFunct = kFunc;
            }
            else if (shearIt != scalarCoef_.end() && poissonIt != scalarCoef_.end()) {
              vars["G"] = shearIt->second;
              vars["nu"] = poissonIt->second;
              shared_ptr< CoefFunctionCompound<Double> >
                  kFunc(new CoefFunctionCompound<Double>(mp_));
              kFunc->SetScalar("2*G_R*(1+nu_R)/(3*(1-2*nu_R))", vars);
              mFunct = kFunc;
            }
            break;

          case SMOOTH_EMODULUS:
            if (bulkIt != scalarCoef_.end() && shearIt != scalarCoef_.end()) {
              vars["K"] = bulkIt->second;
              vars["G"] = shearIt->second;
              shared_ptr< CoefFunctionCompound<Double> >
                  kFunc(new CoefFunctionCompound<Double>(mp_));
              kFunc->SetScalar("9*K_R*G_R/(3*K_R+G_R)", vars);
              mFunct = kFunc;
            }
            else if (bulkIt != scalarCoef_.end() && poissonIt != scalarCoef_.end()) {
              vars["K"] = bulkIt->second;
              vars["nu"] = poissonIt->second;
              shared_ptr< CoefFunctionCompound<Double> >
                  kFunc(new CoefFunctionCompound<Double>(mp_));
              kFunc->SetScalar("3*K_R*(1-2*nu_R)", vars);
              mFunct = kFunc;
            }
            else if (shearIt != scalarCoef_.end() && poissonIt != scalarCoef_.end()) {
              vars["G"] = shearIt->second;
              vars["nu"] = poissonIt->second;
              shared_ptr< CoefFunctionCompound<Double> >
                  kFunc(new CoefFunctionCompound<Double>(mp_));
              kFunc->SetScalar("2*G_R*(1+nu_R)", vars);
              mFunct = kFunc;
            }
            break;

          case SMOOTH_GMODULUS:
            if (bulkIt != scalarCoef_.end() && eModIt != scalarCoef_.end()) {
              vars["K"] = bulkIt->second;
              vars["E"] = eModIt->second;
              shared_ptr< CoefFunctionCompound<Double> >
                  kFunc(new CoefFunctionCompound<Double>(mp_));
              kFunc->SetScalar("3*K_R*E_R/(9*K_R-E_R)", vars);
              mFunct = kFunc;
            }
            else if (bulkIt != scalarCoef_.end() && poissonIt != scalarCoef_.end()) {
              vars["K"] = bulkIt->second;
              vars["nu"] = poissonIt->second;
              shared_ptr< CoefFunctionCompound<Double> >
                  kFunc(new CoefFunctionCompound<Double>(mp_));
              kFunc->SetScalar("3*K_R*(1-2*nu_R)/(2*(1+nu_R))", vars);
              mFunct = kFunc;
            }
            else if (eModIt != scalarCoef_.end() && poissonIt != scalarCoef_.end()) {
              vars["E"] = eModIt->second;
              vars["nu"] = poissonIt->second;
              shared_ptr< CoefFunctionCompound<Double> >
                  kFunc(new CoefFunctionCompound<Double>(mp_));
              kFunc->SetScalar("E_R/(2*(1+nu_R))", vars);
              mFunct = kFunc;
            }
            break;

          case SMOOTH_POISSON:
            if (bulkIt != scalarCoef_.end() && eModIt != scalarCoef_.end()) {
              vars["K"] = bulkIt->second;
              vars["E"] = eModIt->second;
              shared_ptr< CoefFunctionCompound<Double> >
                  kFunc(new CoefFunctionCompound<Double>(mp_));
              kFunc->SetScalar("(3*K_R-E_R)/(6*K_R)", vars);
              mFunct = kFunc;
            }
            else if (bulkIt != scalarCoef_.end() && shearIt != scalarCoef_.end()) {
              vars["K"] = bulkIt->second;
              vars["G"] = shearIt->second;
              shared_ptr< CoefFunctionCompound<Double> >
                  kFunc(new CoefFunctionCompound<Double>(mp_));
              kFunc->SetScalar("(3*K_R-2*G_R)/(2*(3*K_R+G_R))", vars);
              mFunct = kFunc;
            }
            else if (eModIt != scalarCoef_.end() && shearIt != scalarCoef_.end()) {
              vars["E"] = eModIt->second;
              vars["G"] = shearIt->second;
              shared_ptr< CoefFunctionCompound<Double> >
                  kFunc(new CoefFunctionCompound<Double>(mp_));
              kFunc->SetScalar("E_R/(2*G_R)-1", vars);
              mFunct = kFunc;
            }
            break;

          default:
            break;
        }
      }

      if (!mFunct) {
        EXCEPTION("Material Data Type '" <<
                  MaterialTypeEnum.ToString(matType) << "' not available for "
                  << "material '" << name_ << "'");
      }
    }
    mFunct->SetCoordinateSystem(this->coosy_);
    return mFunct;
  }

  PtrCoefFct SmoothMaterial::GetSubTensorCoefFnc( MaterialType matType, 
                                                    SubTensorType tensorType,
                                                    Global::ComplexPart matDataType,
                                                    bool transposed ) const
  {
    PtrCoefFct mFunct;

    CoefMap::const_iterator it = tensorCoef_.find(matType);
    if ( it != tensorCoef_.end() ) {
      CoefXprMechSubTensor subTensorXpr(mp_,  it->second );
      subTensorXpr.SetSubTensorType( tensorType, transposed );
      mFunct = CoefFunction::Generate( mp_, matDataType, subTensorXpr );
    }
    else {
      matTypeNotInDataBase(matType, "tensor");
    }

    return mFunct;
  }

  PtrCoefFct SmoothMaterial::GetSubVectorCoefFnc( MaterialType matType,
                                                    SubTensorType tensorType,
                                                    Global::ComplexPart matDataType ) const
  {
    PtrCoefFct mFunct;

    CoefMap::const_iterator it = vectorCoef_.find(matType);
    if ( it !=  vectorCoef_.end() ) {
      CoefXprMechSubVector subTensorXpr(mp_,  it->second );
      subTensorXpr.SetSubTensorType( tensorType );
      if ( subTensorXpr.IsComplex() ) {
        mFunct = CoefFunction::Generate( mp_, matDataType, subTensorXpr );
      }
      else {
        mFunct = CoefFunction::Generate( mp_, Global::REAL, subTensorXpr );
      }
    }
    else {
      matTypeNotInDataBase(matType, "vector");
    }

    return mFunct;
  }

  void SmoothMaterial::CalcIsotropicStiffnessTensorFromEAndPoisson(Matrix<Double>& out, Double emod, Double poi)
  {
    Complex EModul(emod);
    Complex poisson(poi); 
    Complex LameLambda = (poisson*EModul);
    LameLambda /=  ((Complex(1.0,0) + poisson)*(Complex(1.0,0)  - Complex(2.0,0)*poisson));
    Complex LameMu = (EModul)/(Complex(2.0,0)*(Complex(1.0)+poisson));
    
    Matrix<Complex> elasticityTensor;
    CalcComplexIsotropicStiffnessTensor(elasticityTensor, LameLambda, LameMu);
    out = elasticityTensor.GetPart(Global::REAL);
  }

  void SmoothMaterial::CalcIsotropicStiffnessTensorFromLame(Matrix<Double>& out, Double lambda, Double mu)
  { 
    Matrix<Complex> elasticityTensor;
    CalcComplexIsotropicStiffnessTensor(elasticityTensor, static_cast<Complex>(lambda), static_cast<Complex>(mu));
    
    out = elasticityTensor.GetPart(Global::REAL);
  }

  double SmoothMaterial::CalcIsotropyError(const Matrix<double>& tensor, SubTensorType stt)
  {
    double v = CalcIsotropicPoissonsRatio(tensor, stt);
    double E = CalcIsotropicYoungsModulus(tensor, stt);

    // this is FULL
    Matrix<double> full_hom;
    CalcIsotropicStiffnessTensorFromEAndPoisson(full_hom, E, v);
    // eventually reduce
    Matrix<double> hom;
    ComputeSubTensor(hom, stt, full_hom);

    LOG_DBG(mat_smooth) << "MM::CIE E=" << E << " v=" << v << " err=" << hom.DiffNormL1(tensor);
    LOG_DBG2(mat_smooth) << "MM::CIE tensor=" << tensor.ToString();
    LOG_DBG2(mat_smooth) << "MM::CIE full_hom=" << full_hom.ToString();
    LOG_DBG2(mat_smooth) << "MM::CIE hom=" << hom.ToString();

    return hom.DiffNormL1(tensor);
  }

  double SmoothMaterial::CalcIsotropicPoissonsRatio(const Matrix<double>& tensor, SubTensorType subTensor)
  {
    assert(tensor.GetNumCols() == 3 || tensor.GetNumCols() == 6);

    double E11 = tensor[0][0];
    double E12 = tensor[0][1];

    switch(subTensor)
    {
    case FULL:
    case PLANE_STRAIN:
      return E12 / (E11 + E12);

    case PLANE_STRESS:
      return E12 / E11;

    default:
      EXCEPTION("fail");
      return 0.0;
    }
  }

  double SmoothMaterial::CalcIsotropicYoungsModulus(const Matrix<double>& tensor, SubTensorType subTensor)
  {
    double E11 = tensor[0][0];
    double v = CalcIsotropicPoissonsRatio(tensor, subTensor);

    switch(subTensor)
    {
    case FULL:
    case PLANE_STRAIN:
      E11 *= (1.0 + v) * (1.0 - 2.0 * v) / (1.0 - v);
      break;

    case PLANE_STRESS:
      E11 *= (1.0 - v*v);
      break;

    default:
      assert(false);
      break;
    }

    return E11;
  }

  void SmoothMaterial::CalcComplexIsotropicStiffnessTensor(Matrix<Complex>& out,
                                                             Complex LameLambda,
                                                             Complex LameMu)
  {
    out.Resize(6);
    out.Init();
        
    out[0][0] = LameLambda + Complex(2.0,0) * LameMu;
    out[1][1] = LameLambda + Complex(2.0,0) * LameMu;
    out[2][2] = LameLambda + Complex(2.0,0) * LameMu;
    
    out[0][1] = LameLambda;
    out[1][0] = LameLambda;
    out[0][2] = LameLambda;
    out[1][2] = LameLambda;
    out[2][0] = LameLambda;
    out[2][1] = LameLambda;

    out[3][3] = LameMu;
    out[4][4] = LameMu;
    out[5][5] = LameMu;
  }

  template<class T>
  void SmoothMaterial::ComputeSubTensor(Matrix<T>& matMatrix, SubTensorType subTensor, const Matrix<T>& mat_smooth)
  {

    switch(subTensor)
    {
    case AXI:
    {
      UInt nrElemsAxi = 4;
      matMatrix.Resize( nrElemsAxi, nrElemsAxi );
      matMatrix.Init();

      // indices of rows and lines for xy-plane (rr, zz, rz, phiphi)
      UInt rowPtr[] = {1,2,6,3};  
      for ( UInt i=0; i<nrElemsAxi; i++ )
        for ( UInt j=0; j<nrElemsAxi; j++ )
          matMatrix[i][j] = mat_smooth[rowPtr[i]-1][rowPtr[j]-1];
      break;
    }
    case PLANE_STRAIN:
    {
      UInt nrElems = 3;
      matMatrix.Resize(nrElems, nrElems);
      matMatrix.Init();

      // indices of rows and lines for xy-plane (xx, yy, xy)
      UInt rowPtr[] = {1,2,6}; 
      for ( UInt i=0; i<nrElems; i++ )
        for ( UInt j=0; j<nrElems; j++ )
          matMatrix[i][j] = mat_smooth[rowPtr[i]-1][rowPtr[j]-1];
      break;
    }
    case PLANE_STRESS:
    {
      UInt nrElems = 3;
      matMatrix.Resize( nrElems, nrElems );
      matMatrix.Init();

      // This is a bad test for singularity!
      assert( std::abs(mat_smooth[0][0]) > 1.09E-15 && "Singular material tensor when computing plane stress case" );

      // calculate plane stress matrix for xy-plane
      matMatrix[0][0] = mat_smooth[0][0] - mat_smooth[2][0]*mat_smooth[0][2]/mat_smooth[2][2];
      matMatrix[0][1] = mat_smooth[0][1] - mat_smooth[2][1]*mat_smooth[0][2]/mat_smooth[2][2];
      matMatrix[0][2] = mat_smooth[0][5];
      matMatrix[1][0] = mat_smooth[1][0] - mat_smooth[2][0]*mat_smooth[1][2]/mat_smooth[2][2];
      matMatrix[1][1] = mat_smooth[1][1] - mat_smooth[2][1]*mat_smooth[1][2]/mat_smooth[2][2];
      matMatrix[1][2] = mat_smooth[1][5];
      matMatrix[2][0] = mat_smooth[5][0];
      matMatrix[2][1] = mat_smooth[5][1];
      matMatrix[2][2] = mat_smooth[5][5];
      break;
    }
    case FULL:
      matMatrix = mat_smooth; // copy as nothing changes
      break;
    default:
      // PLAIN is unspecific
      subTensorNotAvailable(NO_MATERIAL, subTensor); // shall be clear
      break;
    }
  }
 
  void SmoothMaterial::ComputeFullStiffTensor() {
    if (tensorCoef_.find(SMOOTH_STIFFNESS_TENSOR) != tensorCoef_.end()) {
      return;
    }

    SymmetryType symType = GetSymmetryType(SMOOTH_STIFFNESS_TENSOR);

    StdVector<std::string> tensorReal(36), tensorImag(36);
    tensorReal.Init("0");
    tensorImag.Init("0");

    switch (symType) {
      case ISOTROPIC: {
          std::string lambdaR, lambdaI, muR, muI;

          if (scalarCoef_.find(SMOOTH_EMODULUS) != scalarCoef_.end() &&
              scalarCoef_.find(SMOOTH_POISSON) != scalarCoef_.end())
          {
            if (!scalarCoef_[SMOOTH_EMODULUS]->IsAnalytic() ||
                !scalarCoef_[SMOOTH_POISSON]->IsAnalytic()) {
              EXCEPTION("Cannot calculate stiffness tensor from non-analytical data");
            }

            shared_ptr<CoefFunctionAnalytic> eModul =
                dynamic_pointer_cast<CoefFunctionAnalytic>(scalarCoef_[SMOOTH_EMODULUS]);
            shared_ptr<CoefFunctionAnalytic> poisson =
                dynamic_pointer_cast<CoefFunctionAnalytic>(scalarCoef_[SMOOTH_POISSON]);
            assert(eModul && poisson);

            std::string eR, eI, nuR, nuI;
            eModul->GetStrScalar(eR, eI);
            poisson->GetStrScalar(nuR, nuI);

            lambdaR = Bracket(nuR) + "*" + Bracket(eR) + "/((1+" +
                      Bracket(nuR) + ")*(1-2*" + Bracket(nuR) + "))";
            lambdaI = Bracket(nuI) + "*" + Bracket(eI) + "/((1+" +
                      Bracket(nuI) + ")*(1-2*" + Bracket(nuI) + "))";
            muR = Bracket(eR) + "/(2*(1+" + Bracket(nuR) + "))";
            muI = Bracket(eI) + "/(2*(1+" + Bracket(nuI) + "))";
          }
          else if (scalarCoef_.find(SMOOTH_KMODULUS) != scalarCoef_.end() &&
                   scalarCoef_.find(SMOOTH_GMODULUS) != scalarCoef_.end())
          {
            if (!scalarCoef_[SMOOTH_KMODULUS]->IsAnalytic() ||
                !scalarCoef_[SMOOTH_GMODULUS]->IsAnalytic()) {
              EXCEPTION("Cannot calculate stiffness tensor from non-analytical data");
            }

            shared_ptr<CoefFunctionAnalytic> bulkMod =
                dynamic_pointer_cast<CoefFunctionAnalytic>(scalarCoef_[SMOOTH_KMODULUS]);
            shared_ptr<CoefFunctionAnalytic> shearMod =
                dynamic_pointer_cast<CoefFunctionAnalytic>(scalarCoef_[SMOOTH_GMODULUS]);
            assert(bulkMod && shearMod);

            std::string kR, kI;
            bulkMod->GetStrScalar(kR, kI);
            shearMod->GetStrScalar(muR, muI);

            lambdaR = Bracket(kR) + "-2/3*" + Bracket(muR);
            lambdaI = Bracket(kI) + "-2/3*" + Bracket(muI);
          }
          else if (scalarCoef_.find(SMOOTH_LAME_LAMBDA) != scalarCoef_.end() &&
                   scalarCoef_.find(SMOOTH_LAME_MU) != scalarCoef_.end())
          {
            if (!scalarCoef_[SMOOTH_LAME_LAMBDA]->IsAnalytic() ||
                !scalarCoef_[SMOOTH_LAME_MU]->IsAnalytic()) {
              EXCEPTION("Cannot calculate stiffness tensor from non-analytical data");
            }

            shared_ptr<CoefFunctionAnalytic> l =
                dynamic_pointer_cast<CoefFunctionAnalytic>(scalarCoef_[SMOOTH_LAME_LAMBDA]);
            shared_ptr<CoefFunctionAnalytic> m =
                dynamic_pointer_cast<CoefFunctionAnalytic>(scalarCoef_[SMOOTH_LAME_MU]);
            assert(l && m);

            l->GetStrScalar(lambdaR, lambdaI);
            m->GetStrScalar(muR, muI);
          }
          else {
            EXCEPTION("Unsupported definition of property '"
                << MaterialTypeEnum.ToString(SMOOTH_STIFFNESS_TENSOR) << "'.");
          }

          std::string diagR = Bracket(lambdaR) + "+2*" + Bracket(muR);
          std::string diagI = Bracket(lambdaI) + "+2*" + Bracket(muI);

          tensorReal[0] = diagR;
          tensorImag[0] = diagI;
          tensorReal[7] = diagR;
          tensorImag[7] = diagI;
          tensorReal[14] = diagR;
          tensorImag[14] = diagI;

          tensorReal[1] = lambdaR;
          tensorImag[1] = lambdaI;
          tensorReal[2] = lambdaR;
          tensorImag[2] = lambdaI;
          tensorReal[6] = lambdaR;
          tensorImag[6] = lambdaI;
          tensorReal[8] = lambdaR;
          tensorImag[8] = lambdaI;
          tensorReal[12] = lambdaR;
          tensorImag[12] = lambdaI;
          tensorReal[13] = lambdaR;
          tensorImag[13] = lambdaI;

          tensorReal[21] = muR;
          tensorImag[21] = muI;
          tensorReal[28] = muR;
          tensorImag[28] = muI;
          tensorReal[35] = muR;
          tensorImag[35] = muI;
        }
        break;

      case TRANS_ISOTROPIC: {
          if (scalarCoef_.find(SMOOTH_EMODULUS) == scalarCoef_.end() ||
              scalarCoef_.find(SMOOTH_EMODULUS_3) == scalarCoef_.end() ||
              scalarCoef_.find(SMOOTH_GMODULUS) == scalarCoef_.end() ||
              scalarCoef_.find(SMOOTH_GMODULUS_3) == scalarCoef_.end() ||
              scalarCoef_.find(SMOOTH_POISSON) == scalarCoef_.end() ||
              scalarCoef_.find(SMOOTH_POISSON_3) == scalarCoef_.end()) {
            EXCEPTION("Transversly isotropic definition of stiffness tensor is incomplete");
          }
          if (!scalarCoef_[SMOOTH_EMODULUS]->IsAnalytic() ||
              !scalarCoef_[SMOOTH_EMODULUS_3]->IsAnalytic() ||
              !scalarCoef_[SMOOTH_GMODULUS]->IsAnalytic() ||
              !scalarCoef_[SMOOTH_GMODULUS_3]->IsAnalytic() ||
              !scalarCoef_[SMOOTH_POISSON]->IsAnalytic() ||
              !scalarCoef_[SMOOTH_POISSON_3]->IsAnalytic()) {
            EXCEPTION("Cannot calculate stiffness tensor from non-analytical data");
          }

          shared_ptr<CoefFunctionAnalytic> E =
              dynamic_pointer_cast<CoefFunctionAnalytic>(scalarCoef_[SMOOTH_EMODULUS]);
          shared_ptr<CoefFunctionAnalytic> E3 =
              dynamic_pointer_cast<CoefFunctionAnalytic>(scalarCoef_[SMOOTH_EMODULUS_3]);
          shared_ptr<CoefFunctionAnalytic> G =
              dynamic_pointer_cast<CoefFunctionAnalytic>(scalarCoef_[SMOOTH_GMODULUS]);
          shared_ptr<CoefFunctionAnalytic> G3 =
              dynamic_pointer_cast<CoefFunctionAnalytic>(scalarCoef_[SMOOTH_EMODULUS_3]);
          shared_ptr<CoefFunctionAnalytic> nu =
              dynamic_pointer_cast<CoefFunctionAnalytic>(scalarCoef_[SMOOTH_POISSON]);
          shared_ptr<CoefFunctionAnalytic> nu13 =
              dynamic_pointer_cast<CoefFunctionAnalytic>(scalarCoef_[SMOOTH_POISSON_3]);
          assert(E && E3 && G && G3 && nu && nu13);

          std::string eR, eI, e3R, e3I, gR, gI, g3R, g3I, nuR, nuI, nu13R, nu13I;
          E->GetStrScalar(eR, eI);
          E3->GetStrScalar(e3R, e3I);
          G->GetStrScalar(gR, gI);
          G3->GetStrScalar(g3R, g3I);
          nu->GetStrScalar(nuR, nuI);
          nu13->GetStrScalar(nu13R, nu13I);

          std::string nu31R = Bracket(Bracket(e3R) + "/" + Bracket(eR)) + "*" + Bracket(nu13R);
          std::string nu31I = Bracket(Bracket(e3I) + "/" + Bracket(eI)) + "*" + Bracket(nu13I);

          std::string auxR = "(1+" + Bracket(nuR) + ")*(1-" + Bracket(nuR) +
                             "-2*" + Bracket(nu13R) + "*" + Bracket(nu31R) + ")";
          std::string auxI = "(1+" + Bracket(nuI) + ")*(1-" + Bracket(nuI) +
                             "-2*" + Bracket(nu13I) + "*" + Bracket(nu31I) + ")";

          tensorReal[0] = Bracket(eR) + "*(1-" + Bracket(nu13R) + "*" +
                          Bracket(nu31R) + ")/" + Bracket(auxR);
          tensorImag[0] = Bracket(eI) + "*(1-" + Bracket(nu13I) + "*" +
                          Bracket(nu31I) + ")/" + Bracket(auxI);
          tensorReal[7] = tensorReal[0];
          tensorImag[7] = tensorImag[0];
          tensorReal[14] = Bracket(e3R) + "*(1-" + Bracket(nuR) + "*" + Bracket(nuR) +
                           "/" + Bracket(auxR);
          tensorImag[14] = Bracket(e3I) + "*(1-" + Bracket(nuI) + "*" + Bracket(nuI) +
                           "/" + Bracket(auxI);

          tensorReal[1] = Bracket(eR) + "*(" + Bracket(nuR) + "+" + Bracket(nu13R) +
                          "*" + Bracket(nu31R) + "/" + Bracket(auxR);
          tensorImag[1] = Bracket(eI) + "*(" + Bracket(nuI) + "+" + Bracket(nu13I) +
                          "*" + Bracket(nu31I) + "/" + Bracket(auxI);
          tensorReal[2] = Bracket(eR) + "*((1+" + Bracket(nuR) + ")*" +
                          Bracket(nu31R) + ")/" + Bracket(auxR);
          tensorImag[2] = Bracket(eI) + "*((1+" + Bracket(nuI) + ")*" +
                          Bracket(nu31I) + ")/" + Bracket(auxI);
          tensorReal[6] = tensorReal[1];
          tensorImag[6] = tensorImag[1];
          tensorReal[8] = tensorReal[2];
          tensorImag[8] = tensorImag[2];
          tensorReal[12] = tensorReal[2];
          tensorImag[12] = tensorImag[2];
          tensorReal[13] = tensorReal[2];
          tensorImag[13] = tensorImag[2];

          tensorReal[21] = g3R;
          tensorImag[21] = g3I;
          tensorReal[28] = g3R;
          tensorImag[28] = g3I;
          tensorReal[35] = gR;
          tensorImag[35] = gI;
        }
        break;

      case ORTHOTROPIC: {
          if (scalarCoef_.find(SMOOTH_EMODULUS_1) == scalarCoef_.end() ||
              scalarCoef_.find(SMOOTH_EMODULUS_2) == scalarCoef_.end() ||
              scalarCoef_.find(SMOOTH_EMODULUS_3) == scalarCoef_.end() ||
              scalarCoef_.find(SMOOTH_GMODULUS_23) == scalarCoef_.end() ||
              scalarCoef_.find(SMOOTH_GMODULUS_13) == scalarCoef_.end() ||
              scalarCoef_.find(SMOOTH_GMODULUS_12) == scalarCoef_.end() ||
              scalarCoef_.find(SMOOTH_POISSON_12) == scalarCoef_.end() ||
              scalarCoef_.find(SMOOTH_POISSON_23) == scalarCoef_.end() ||
              scalarCoef_.find(SMOOTH_POISSON_13) == scalarCoef_.end()) {
            EXCEPTION("Orthotropic definition of stiffness tensor is incomplete");
          }
          if (!scalarCoef_[SMOOTH_EMODULUS_1]->IsAnalytic() ||
              !scalarCoef_[SMOOTH_EMODULUS_2]->IsAnalytic() ||
              !scalarCoef_[SMOOTH_EMODULUS_3]->IsAnalytic() ||
              !scalarCoef_[SMOOTH_GMODULUS_23]->IsAnalytic() ||
              !scalarCoef_[SMOOTH_GMODULUS_13]->IsAnalytic() ||
              !scalarCoef_[SMOOTH_GMODULUS_12]->IsAnalytic() ||
              !scalarCoef_[SMOOTH_POISSON_12]->IsAnalytic() ||
              !scalarCoef_[SMOOTH_POISSON_23]->IsAnalytic() ||
              !scalarCoef_[SMOOTH_POISSON_13]->IsAnalytic()) {
            EXCEPTION("Cannot calculate stiffness tensor from non-analytical data");
          }

          shared_ptr<CoefFunctionAnalytic> E1 =
                        dynamic_pointer_cast<CoefFunctionAnalytic>(scalarCoef_[SMOOTH_EMODULUS_1]);
          shared_ptr<CoefFunctionAnalytic> E2 =
                        dynamic_pointer_cast<CoefFunctionAnalytic>(scalarCoef_[SMOOTH_EMODULUS_2]);
          shared_ptr<CoefFunctionAnalytic> E3 =
                        dynamic_pointer_cast<CoefFunctionAnalytic>(scalarCoef_[SMOOTH_EMODULUS_3]);
          shared_ptr<CoefFunctionAnalytic> nu12 =
                        dynamic_pointer_cast<CoefFunctionAnalytic>(scalarCoef_[SMOOTH_POISSON_12]);
          shared_ptr<CoefFunctionAnalytic> nu23 =
                        dynamic_pointer_cast<CoefFunctionAnalytic>(scalarCoef_[SMOOTH_POISSON_23]);
          shared_ptr<CoefFunctionAnalytic> nu13 =
                        dynamic_pointer_cast<CoefFunctionAnalytic>(scalarCoef_[SMOOTH_POISSON_13]);
          shared_ptr<CoefFunctionAnalytic> G12 =
                        dynamic_pointer_cast<CoefFunctionAnalytic>(scalarCoef_[SMOOTH_GMODULUS_12]);
          shared_ptr<CoefFunctionAnalytic> G23 =
                        dynamic_pointer_cast<CoefFunctionAnalytic>(scalarCoef_[SMOOTH_GMODULUS_23]);
          shared_ptr<CoefFunctionAnalytic> G13 =
                        dynamic_pointer_cast<CoefFunctionAnalytic>(scalarCoef_[SMOOTH_GMODULUS_13]);
          assert(E1 && E2 && E3 && nu12 && nu23 && nu13 && G12 && G23 && G13);

          std::string e1R, e1I, e2R, e2I, e3R, e3I, nu12R, nu12I, nu23R, nu23I,
                      nu13R, nu13I, g12R, g12I, g23R, g23I, g13R, g13I;
          E1->GetStrScalar(e1R, e1I);
          E2->GetStrScalar(e2R, e2I);
          E3->GetStrScalar(e3R, e3I);
          nu12->GetStrScalar(nu12R, nu12I);
          nu23->GetStrScalar(nu23R, nu23I);
          nu13->GetStrScalar(nu13R, nu13I);
          G12->GetStrScalar(g12R, g12I);
          G23->GetStrScalar(g23R, g23I);
          G13->GetStrScalar(g13R, g13I);

          // http://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_orthotropic.cfm
          std::string nu21R = Bracket( Bracket(e2R) + "/" + Bracket(e1R)) + "*" +
                              Bracket(nu12R);
          std::string nu21I = Bracket( Bracket(e2I) + "/" + Bracket(e1I)) + "*" +
                              Bracket(nu12I);
          std::string nu32R = Bracket( Bracket(e3R) + "/" + Bracket(e2R)) + "*" +
                              Bracket(nu23R);
          std::string nu32I = Bracket( Bracket(e3I) + "/" + Bracket(e2I)) + "*" +
                              Bracket(nu23I);
          std::string nu31R = Bracket( Bracket(e3R) + "/" + Bracket(e1R)) + "*" +
                              Bracket(nu13R);
          std::string nu31I = Bracket( Bracket(e3I) + "/" + Bracket(e1I)) + "*" +
                              Bracket(nu13I);
          std::string auxR = "(1-" + Bracket(nu12R) + "*" + Bracket(nu21R) + "-" +
                             Bracket(nu23R) + "*" + Bracket(nu32R) + "-" +
                             Bracket(nu13R) + "*" + Bracket(nu31R) + "-2*" +
                             Bracket(nu12R) + "*" + Bracket(nu23R) + "*" +
                             Bracket(nu31R) + ")/(" + Bracket(e1R) + "*" +
                             Bracket(e2R) + "*" + Bracket(e3R) + ")";
          std::string auxI = "(1-" + Bracket(nu12I) + "*" + Bracket(nu21I) + "-" +
                             Bracket(nu23I) + "*" + Bracket(nu32I) + "-" +
                             Bracket(nu13I) + "*" + Bracket(nu31I) + "-2*" +
                             Bracket(nu12I) + "*" + Bracket(nu23I) + "*" +
                             Bracket(nu31I) + ")/(" + Bracket(e1I) + "*" +
                             Bracket(e2I) + "*" + Bracket(e3I) + ")";

          tensorReal[0] = "(1-" + Bracket(nu23R) + "*" + Bracket(nu32R) + ")/(" +
                          Bracket(e2R) + "*" + Bracket(e3R) + "*" + Bracket(auxR) + ")";
          tensorImag[0] = "(1-" + Bracket(nu23I) + "*" + Bracket(nu32I) + ")/(" +
                          Bracket(e2I) + "*" + Bracket(e3I) + "*" + Bracket(auxI) + ")";
          tensorReal[7] = "(1-" + Bracket(nu13R) + "*" + Bracket(nu31R) + ")/(" +
                          Bracket(e1R) + "*" + Bracket(e3R) + "*" + Bracket(auxR) + ")";
          tensorImag[7] = "(1-" + Bracket(nu13I) + "*" + Bracket(nu31I) + ")/(" +
                          Bracket(e1I) + "*" + Bracket(e3I) + "*" + Bracket(auxI) + ")";
          tensorReal[14] = "(1-" + Bracket(nu12R) + "*" + Bracket(nu21R) + ")/(" +
                          Bracket(e1R) + "*" + Bracket(e2R) + "*" + Bracket(auxR) + ")";
          tensorImag[14] = "(1-" + Bracket(nu12I) + "*" + Bracket(nu21I) + ")/(" +
                          Bracket(e1I) + "*" + Bracket(e2I) + "*" + Bracket(auxI) + ")";

          tensorReal[1] = "(" + Bracket(nu21R) + "+" + Bracket(nu31R) + "*" +
                          Bracket(nu23R) + ")/(" + Bracket(e2R) + "*" +
                          Bracket(e3R) + "*" + Bracket(auxR) + ")";
          tensorImag[1] = "(" + Bracket(nu21I) + "+" + Bracket(nu31I) + "*" +
                          Bracket(nu23I) + ")/(" + Bracket(e2I) + "*" +
                          Bracket(e3I) + "*" + Bracket(auxI) + ")";
          tensorReal[2] = "(" + Bracket(nu31R) + "+" + Bracket(nu21R) + "*" +
                          Bracket(nu32R) + ")/(" + Bracket(e2R) + "*" +
                          Bracket(e3R) + "*" + Bracket(auxR) + ")";
          tensorImag[2] = "(" + Bracket(nu31I) + "+" + Bracket(nu21I) + "*" +
                          Bracket(nu32I) + ")/(" + Bracket(e2I) + "*" +
                          Bracket(e3I) + "*" + Bracket(auxI) + ")";
          tensorReal[6] = "(" + Bracket(nu12R) + "+" + Bracket(nu13R) + "*" +
                          Bracket(nu32R) + ")/(" + Bracket(e1R) + "*" +
                          Bracket(e3R) + "*" + Bracket(auxR) + ")";
          tensorImag[6] = "(" + Bracket(nu12I) + "+" + Bracket(nu13I) + "*" +
                          Bracket(nu32I) + ")/(" + Bracket(e1I) + "*" +
                          Bracket(e3I) + "*" + Bracket(auxI) + ")";
          tensorReal[8] = "(" + Bracket(nu32R) + "+" + Bracket(nu31R) + "*" +
                          Bracket(nu12R) + ")/(" + Bracket(e1R) + "*" +
                          Bracket(e3R) + "*" + Bracket(auxR) + ")";
          tensorImag[8] = "(" + Bracket(nu32I) + "+" + Bracket(nu31I) + "*" +
                          Bracket(nu12I) + ")/(" + Bracket(e1I) + "*" +
                          Bracket(e3I) + "*" + Bracket(auxI) + ")";
          tensorReal[12] = "(" + Bracket(nu13R) + "+" + Bracket(nu12R) + "*" +
                           Bracket(nu13R) + ")/(" + Bracket(e1R) + "*" +
                           Bracket(e2R) + "*" + Bracket(auxR) + ")";
          tensorImag[12] = "(" + Bracket(nu13I) + "+" + Bracket(nu12I) + "*" +
                           Bracket(nu13I) + ")/(" + Bracket(e1I) + "*" +
                           Bracket(e2I) + "*" + Bracket(auxI) + ")";
          tensorReal[13] = "(" + Bracket(nu23R) + "+" + Bracket(nu13R) + "*" +
                           Bracket(nu21R) + ")/(" + Bracket(e1R) + "*" +
                           Bracket(e2R) + "*" + Bracket(auxR) + ")";
          tensorImag[13] = "(" + Bracket(nu23I) + "+" + Bracket(nu13I) + "*" +
                           Bracket(nu21I) + ")/(" + Bracket(e1I) + "*" +
                           Bracket(e2I) + "*" + Bracket(auxI) + ")";

          tensorReal[21] = g23R;
          tensorImag[21] = g23I;
          tensorReal[28] = g13R;
          tensorImag[28] = g13I;
          tensorReal[35] = g12R;
          tensorImag[35] = g12I;
        }
        break;

      case GENERAL:
        EXCEPTION("Property '" << MaterialTypeEnum.ToString(SMOOTH_STIFFNESS_TENSOR)
            << "' is undefined.");
        break;

      default:
        EXCEPTION("Unknown symmetry type '"
            << SymmetryTypeEnum.ToString(symType) << "' for property '"
            << MaterialTypeEnum.ToString(SMOOTH_STIFFNESS_TENSOR) << "'.");
    }

    PtrCoefFct c = CoefFunction::Generate(mp_, Global::COMPLEX, 6, 6,
                                           tensorReal, tensorImag);
    SetCoefFct(SMOOTH_STIFFNESS_TENSOR, c);
  }

  // Compute the constant full stiffness tensor from isotropic,
  // transversely isotropic or orthotropic data
  Matrix<Complex> SmoothMaterial::GetFullStiffTensor(
      BaseMaterial::SymmetryType symType, BaseMaterial::CoefMap &coefMap)
  {
    Matrix<Complex> tensor(6, 6);
    LocPointMapped lpm;

    switch (symType) {
      case ISOTROPIC: {
          Complex LameLambda, LameMu;

          if (coefMap.find(SMOOTH_EMODULUS) != coefMap.end() &&
              coefMap.find(SMOOTH_POISSON) != coefMap.end())
          {
            Complex eModul, poisson;
            coefMap[SMOOTH_EMODULUS]->GetScalar(eModul, lpm);
            coefMap[SMOOTH_POISSON]->GetScalar(poisson, lpm);

            LameLambda = (poisson*eModul) / ((1.0+poisson) * (1.0-2.0*poisson));
            LameMu = eModul / (2.0 * (1.0 + poisson));
          }
          else if (coefMap.find(SMOOTH_KMODULUS) != coefMap.end() &&
                   coefMap.find(SMOOTH_GMODULUS) != coefMap.end())
          {
            Complex k, g;
            coefMap[SMOOTH_KMODULUS]->GetScalar(k, lpm);
            coefMap[SMOOTH_GMODULUS]->GetScalar(g, lpm);

            LameLambda = k - 2.0/3.0 * g;
            LameMu = g;
          }
          else {
            EXCEPTION("Unsupported definition of property '"
                << MaterialTypeEnum.ToString(SMOOTH_STIFFNESS_TENSOR) << "'.");
          }

          tensor[0][0] = LameLambda + 2.0 * LameMu;
          tensor[1][1] = LameLambda + 2.0 * LameMu;
          tensor[2][2] = LameLambda + 2.0 * LameMu;

          tensor[0][1] = LameLambda;
          tensor[0][2] = LameLambda;
          tensor[1][0] = LameLambda;
          tensor[1][2] = LameLambda;
          tensor[2][0] = LameLambda;
          tensor[2][1] = LameLambda;

          tensor[3][3] = LameMu;
          tensor[4][4] = LameMu;
          tensor[5][5] = LameMu;
        }
        break;

      case TRANS_ISOTROPIC: {
          Complex E, E3, G, G3, nu, nu13;
          coefMap[SMOOTH_EMODULUS]->GetScalar(E, lpm);
          coefMap[SMOOTH_EMODULUS_3]->GetScalar(E3, lpm);
          coefMap[SMOOTH_GMODULUS]->GetScalar(G, lpm);
          coefMap[SMOOTH_GMODULUS_3]->GetScalar(G3, lpm);
          coefMap[SMOOTH_POISSON]->GetScalar(nu, lpm);
          coefMap[SMOOTH_POISSON_3]->GetScalar(nu13, lpm);

          Complex nu31 = (E3/E)*nu13;

          Complex aux = (1.0 + nu) * (1.0 - nu - 2.0*nu13*nu31);

          tensor[0][0] = E * (1.0-nu13*nu31) / aux;
          tensor[1][1] = E * (1.0-nu13*nu31) / aux;
          tensor[2][2] = E3 * (1.0-nu*nu) / aux;

          tensor[0][1] = E * (nu+nu13*nu31) / aux;
          tensor[0][2] = E * ((1.0+nu)*nu31) / aux;
          tensor[1][0] = tensor[0][1];
          tensor[1][2] = tensor[0][2];
          tensor[2][0] = tensor[0][2];
          tensor[2][1] = tensor[0][2];

          tensor[3][3] = G3;
          tensor[4][4] = G3;
          tensor[5][5] = G;
        }
        break;

      case ORTHOTROPIC: {
          // http://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_orthotropic.cfm
          Complex E1, E2, E3, nu12, nu23, nu13, G23, G31, G12;
          coefMap[SMOOTH_EMODULUS_1]->GetScalar(E1, lpm);
          coefMap[SMOOTH_EMODULUS_2]->GetScalar(E2, lpm);
          coefMap[SMOOTH_EMODULUS_3]->GetScalar(E3, lpm);
          coefMap[SMOOTH_POISSON_12]->GetScalar(nu12, lpm);
          coefMap[SMOOTH_POISSON_23]->GetScalar(nu23, lpm);
          coefMap[SMOOTH_POISSON_13]->GetScalar(nu13, lpm);
          coefMap[SMOOTH_GMODULUS_23]->GetScalar(G23, lpm);
          coefMap[SMOOTH_GMODULUS_13]->GetScalar(G31, lpm);
          coefMap[SMOOTH_GMODULUS_12]->GetScalar(G12, lpm);

          Complex nu21 = (E2/E1)*nu12;
          Complex nu32 = (E3/E2)*nu23;
          Complex nu31 = (E3/E1)*nu13;

          Complex aux = (1.0 - nu12*nu21 - nu23*nu32 - nu13*nu31 -
                         2.0*nu12*nu23*nu31) / (E1*E2*E3);

          tensor[0][0] = (1.0-nu23*nu32)/(E2*E3*aux);
          tensor[1][1] = (1.0-nu13*nu31)/(E1*E3*aux);
          tensor[2][2] = (1.0-nu12*nu21)/(E1*E2*aux);

          tensor[0][1] = (nu21+nu31*nu23)/(E2*E3*aux);
          tensor[0][2] = (nu31+nu21*nu32)/(E2*E3*aux);
          tensor[1][0] = (nu12+nu13*nu32)/(E1*E3*aux);
          tensor[1][2] = (nu32+nu31*nu12)/(E1*E3*aux);
          tensor[2][0] = (nu13+nu12*nu23)/(E1*E2*aux);
          tensor[2][1] = (nu23+nu13*nu21)/(E1*E2*aux);

          tensor[3][3] = G23;
          tensor[4][4] = G31;
          tensor[5][5] = G12;
        }
        break;

      case GENERAL:
        if (coefMap.find(SMOOTH_STIFFNESS_TENSOR) == coefMap.end()) {
          EXCEPTION("Property '" << MaterialTypeEnum.ToString(SMOOTH_STIFFNESS_TENSOR)
              << "' is undefined.");
        }

        coefMap[SMOOTH_STIFFNESS_TENSOR]->GetTensor(tensor, lpm);
        break;

      default:
        EXCEPTION("Unknown symmetry type '"
            << SymmetryTypeEnum.ToString(symType) << "' for property '"
            << MaterialTypeEnum.ToString(SMOOTH_STIFFNESS_TENSOR) << "'.");
    }

    return tensor;
  }

  StdVector<double> SmoothMaterial::CalcOrthotropeYoungsModulus(const Matrix<double>& tensor, BaseMaterial* mat_smooth, SubTensorType stt, double vol)
  {
    Matrix<double> D;
    tensor.Invert(D);

    assert(tensor.GetNumRows() == (stt == FULL ? 6 : 3));
    StdVector<double> res(stt == PLANE_STRESS ? 2 : 3);

    switch(stt)
    {
    case FULL:
      for(UInt i = 0; i < 3; i++)
        res[i] = 1.0/D[i][i];
      break;

    case PLANE_STRAIN:
    {
      assert(mat_smooth != NULL && vol > 0 && vol <= 1.000001);

      // core properties
      double E_core, v_core;
      mat_smooth->GetScalar(E_core, SMOOTH_EMODULUS, Global::REAL);
      mat_smooth->GetScalar(v_core, SMOOTH_POISSON, Global::REAL);

      double E3 = E_core * vol;

      res[0] = 1.0 / (D[0][0] + v_core*v_core / E3);
      res[1] = 1.0 / (D[1][1] + v_core*v_core / E3);
      res[2] = E3;
      break;
    }

    case PLANE_STRESS:
    {
      StdVector<double> v = CalcOrthotropePoissonsRatio(tensor, mat_smooth, stt, vol);
      assert(v.GetSize() == 2);

      // E1 = C(1,1) * (1-v12*v21)
      res[0] = tensor[1-1][1-1] * (1.0 - v[0] * v[1]);
      // E2 = C(2,2) * (1-v12*v21)
      res[1] = tensor[2-1][2-1] * (1.0 - v[0] * v[1]);
      break;
    }
    default:
      assert(false);
      break;
    }
    return res;
  }

  StdVector<double> SmoothMaterial::CalcOrthotropePoissonsRatio(const Matrix<double>& tensor, BaseMaterial* mat_smooth, SubTensorType stt, double vol)
  {
    Matrix<double> D;
    tensor.Invert(D);

    StdVector<double> res;
    StdVector<double> E; // v_21, v_12, v_31, v_13, v_32, v_23



    switch(stt)
    {
    case FULL:
      E = CalcOrthotropeYoungsModulus(tensor, mat_smooth, stt, vol); // does the asserts!
      res.Push_back(-1.0 * D[1-1][2-1] * E[2-1]); // v_21
      res.Push_back(-1.0 * D[2-1][1-1] * E[1-1]); // v_12
      res.Push_back(-1.0 * D[1-1][3-1] * E[3-1]); // v_31
      res.Push_back(-1.0 * D[3-1][1-1] * E[1-1]); // v_13
      res.Push_back(-1.0 * D[2-1][3-1] * E[3-1]); // v_32
      res.Push_back(-1.0 * D[3-1][2-1] * E[2-1]); // v_23
      break;

    case PLANE_STRAIN:
    {
      E = CalcOrthotropeYoungsModulus(tensor, mat_smooth, stt, vol); // does the asserts!
      // core properties
      double E_core, v_core;
      mat_smooth->GetScalar(E_core, SMOOTH_EMODULUS, Global::REAL);
      mat_smooth->GetScalar(v_core, SMOOTH_POISSON, Global::REAL);

      // % the two 2D poisson ratios and Es have been derived by Bastian and Fabian S.
      // v12 = (-B(1,2)-nusteg^2/rho/Esteg)*E1
      // v21 = E2/E1*v12
      double v_12 = (-D[1-1][2-1] - v_core * v_core / E[3-1]) * E[1-1];
      double v_21 = (E[2-1] / E[1-1]) * v_12;

      // % the remaining poisson ratios
      // v13 = E1*nusteg/rho/Esteg
      // v31 = E3/E1*v13
      double v_13 = E[1-1] * v_core / E[3-1];
      double v_31 = (E[3-1] / E[1-1]) * v_13;
      // v23 = E2*nusteg/rho/Esteg
      // v32 = E3/E2*v23
      double v_23 = E[2-1] * v_core / E[3-1];
      double v_32 = (E[3-1] / E[2-1]) * v_23;

      res.Push_back(v_21);
      res.Push_back(v_12);
      res.Push_back(v_31);
      res.Push_back(v_13);
      res.Push_back(v_32);
      res.Push_back(v_23);
      break;
    }

    case PLANE_STRESS:
      // v21 = c(1,2)/c(1,1)
      res.Push_back(tensor[1-1][2-1] / tensor[1-1][1-1]);

      // v12 = C(1,2)/C(2,2)
      res.Push_back(tensor[1-1][2-1] / tensor[2-1][2-1]);
      break;

    default:
      assert(false);
      break;
    }

    return res;
  }

  double SmoothMaterial::CalcOrthotropeError(const Matrix<double>& tensor)
  {
    // measure zeros in 1-norm
    double err = 0.0;
    if(tensor.GetNumRows() == 3)
    {
       err += abs(tensor[1-1][3-1]);
       err += abs(tensor[2-1][3-1]);
    }
    else
    {
      for(UInt r = 1; r <= 5; r++)
      {
        for(UInt c = 4; c <= 6; c++)
        {
          // only above the diagonal
          if(r >= c) continue;
          err += abs(tensor[r-1][c-1]);
        }
      }
    }
    return err;
  }

  StdVector<std::pair<std::string, double> > SmoothMaterial::CalcIsotropicProperties(const Matrix<double>& tensor, SubTensorType stt)
  {
    double E = CalcIsotropicYoungsModulus(tensor, stt);
    double v = CalcIsotropicPoissonsRatio(tensor, stt);
    double G = E / (2.0 + 2.0*v); // comparing to the tensor gives an idea about the error
    double err = CalcIsotropyError(tensor, stt);

    StdVector<std::pair<std::string, double> > res;
    res.Push_back(std::make_pair("E", E));
    res.Push_back(std::make_pair("v", v));
    res.Push_back(std::make_pair("G", G));
    res.Push_back(std::make_pair("err", err));

    return res;
  }

  StdVector<std::pair<std::string, double> > SmoothMaterial::CalcOrthotropeProperties(const Matrix<double>& tensor, BaseMaterial* bm, SubTensorType stt, double vol)
  {
    LOG_DBG2(mat_smooth) << "GOP tensor=" << tensor.ToString();
    Matrix<double> D;
    tensor.Invert(D);

    StdVector<std::pair<std::string, double> > res;

    StdVector<double> E = CalcOrthotropeYoungsModulus(tensor, bm, stt, vol);

    res.Push_back(std::make_pair("E_1", E[0]));
    res.Push_back(std::make_pair("E_2", E[1]));
    if(stt == FULL) {
      res.Push_back(std::make_pair("E_3", E[2]));
    }

    StdVector<double> v = CalcOrthotropePoissonsRatio(tensor, bm, stt, vol);
    // v_21=0, v_12=1, v_31=2, v_13=3, v_32=4, v_23=5
    res.Push_back(std::make_pair("v_21", v[0]));
    res.Push_back(std::make_pair("v_12", v[1]));
    if(stt == FULL || stt == PLANE_STRAIN) {
      res.Push_back(std::make_pair("v_31", v[2]));
      res.Push_back(std::make_pair("v_13", v[3]));
      res.Push_back(std::make_pair("v_32", v[4]));
      res.Push_back(std::make_pair("v_23", v[5]));
    }

    if(stt == FULL) {
      res.Push_back(std::make_pair("G_23", 1.0 / D[4-1][4-1]));
      res.Push_back(std::make_pair("G_13", 1.0 / D[5-1][5-1]));
      res.Push_back(std::make_pair("G_12", 1.0 / D[6-1][6-1]));
    } else {
      res.Push_back(std::make_pair("G_12", 1.0 / D[3-1][3-1]));
    }

    double err = CalcOrthotropeError(tensor);
    res.Push_back(std::make_pair("err", err));

    return res;
  }

  // required in ErsatzMaterial. The complex version is explictely called here
  template void SmoothMaterial::ComputeSubTensor<Double>(Matrix<Double>& matMatrix, SubTensorType subTensor, const Matrix<Double>& mat_smooth);

} // end of namespace
