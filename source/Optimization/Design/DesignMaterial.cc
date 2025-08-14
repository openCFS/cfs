#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <boost/iterator/counting_iterator.hpp>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/ParamHandling/XmlReader.hh"
#include "Domain/CoefFunction/CoefFunctionOpt.hh"
#include "Domain/Domain.hh"
#include "Domain/Mesh/Grid.hh"
#include "FeBasis/H1/H1Elems.hh"
#include "FeBasis/H1/H1ElemsLagExpl.hh"
#include "General/defs.hh"
#include "General/Exception.hh"
#include "Optimization/Design/DesignMaterial.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/MaterialTensor.hh"
#include "Optimization/Design/FeatureMappingDesign.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/TransferFunction.hh"
#include "PDE/SinglePDE.hh"
#include "Utils/CubicInterpolate.hh"
#include "Utils/BiCubicInterpolate.hh"
#include "Utils/TriCubicInterpolate.hh"
#include "Utils/StdVector.hh"

#ifdef USE_SGPP
#include "sgpp_base.hpp"
#include "sgpp_optimization.hpp"
#endif

DEFINE_LOG(dm, "designMaterial")

using namespace CoupledField;
using std::string;

Enum<DesignMaterial::Type> DesignMaterial::type;
Enum<DesignMaterial::TransIsoType> DesignMaterial::transIsoType;
Enum<DesignMaterial::RotationType> DesignMaterial::rotationType;

DesignMaterial::DesignMaterial(PtrParamNode pn, OptimizationMaterial::System material, StdVector<DesignID>& design, DesignSpace* space)
#ifdef USE_SGPP
  :
  alpha1_(sgpp::base::DataVector(0)),
  alpha2_(sgpp::base::DataVector(0)),
  alpha3_(sgpp::base::DataVector(0)),
  alpha4_(sgpp::base::DataVector(0)),
  alpha5_(sgpp::base::DataVector(0)),
  alpha6_(sgpp::base::DataVector(0)),
  alpha7_(sgpp::base::DataVector(0))
#endif
{
  type_ = type.Parse(pn->Get("type")->As<string>());

  dim = domain->GetDim();

  transIsoType_ = transIsoType.Parse(pn->Get("isoplane")->As<string>());

  massIsDesign_ = pn->Get("optimizeMass")->As<bool>();

  massFactor_ = pn->Get("massFactor")->As<double>();

  trace_ = pn->Get("trace")->As<double>();

  dampingIsDesign_ = pn->Get("optimizeDamping")->As<bool>();

  space_ = space;

  if(pn->Has("rotationtype"))
    rotationType_ = rotationType.Parse(pn->Get("rotationtype")->As<std::string>());

  bias_ = pn->Get("bias")->As<bool>();
  // TODO: set biasMaterials_

  // initialize, maybe overwritten later
  interpolation_ = DesignMaterial::NOTYPE;

  // collect all designs here, to check whether all are given
  unsigned int r = RequiredParameters(material);
  StdVector<DesignElement::Type> d;
  d.Reserve(r);
  // copy the ones from DesignSpace
  for (unsigned int i = 0; i < design.GetSize(); ++i) {
    d.Push_back(design[i].design);
  }

  shearIsDesign_ = d.Find(DesignElement::SHEAR1) == -1 ? 0 : 1;

  // read non-design parameters
  ParamNodeList params = pn->GetList("param");
  for (unsigned int i = 0; i < params.GetSize(); i++) {
    DesignElement::Type dt = DesignElement::type.Parse(
        params[i]->Get("name")->As<string>());
    SetParameter(dt, params[i]->Get("value")->As<Double>(), true);
    if (d.Find(dt) < 0) {
      d.Push_back(dt);
    }
    if (d.Find(DesignElement::SHEAR1) < 0 && type_ != HOM_RECT_C1) {
      SetParameter(DesignElement::SHEAR1, 0.5, true);
    }
  }
  if (!CheckRequiredDesigns(d)) {
    throw Exception("Not all Parameters for chosen DesignMaterial given. See DesignMaterial::CheckRequiredDesigns().");
  }else if(d.GetSize() > r){ // design.GetSize() < r is impossible as CheckRequiredDesigns passed
    domain->GetInfoRoot()->Get("optimization/designSpace/header")->SetWarning("There are designs specified that are not used!");
  }

  if(type_ == HOM_RECT || type_ == D_HOM_RECT)
  {
    PtrParamNode hr = pn->Get("homRect");
    hom_rect_samples_.Resize(9, 6);
    FillHomRectSamples(hr, 0, "0.0", "0.0");
    FillHomRectSamples(hr, 1, "0.5", "0.0");
    FillHomRectSamples(hr, 2, "0.5", "0.5");
    FillHomRectSamples(hr, 3, "0.0", "0.5");
    FillHomRectSamples(hr, 4, "0.25", "0.0");
    FillHomRectSamples(hr, 5, "0.5", "0.25");
    FillHomRectSamples(hr, 6, "0.25", "0.5");
    FillHomRectSamples(hr, 7, "0.0", "0.25");
    FillHomRectSamples(hr, 8, "0.25", "0.25");
  }

  std::string interpolation_str;
  if (type_ == HOM_RECT_C1 || type_ == HOM_ISO_C1 || type_ == HEAT) {
    if (type_ == HEAT && dim !=3)
      EXCEPTION("DesignMaterial Heat only implemented and tested for 3d!");
    string p_node = "";
    if (type_ == HOM_RECT_C1)
      p_node = "homRectC1";
    if (type_ == HOM_ISO_C1)
      p_node = "homIsoC1";
    if (type_ == HEAT)
      p_node = "heat";

    PtrParamNode hr = pn->Get(p_node);
    std::string file = hr->Get("file")->As<std::string>();
    // read interpolation method
    interpolation_str = hr->Get("interpolation")->As<std::string>();

    LOG_DBG3(dm) << "DM: sequence=" << Optimization::context->sequence << " reading coeff file " << file;

    // full C1 interpolation with text coefficients for 2D (without shearing angle)
    if (interpolation_str == "c1_text_2d") {
      interpolation_ = C1;
      // read coefficients from text file
      std::string filename = hr->Get("file")->As<std::string>();
      std::ifstream file(filename.c_str());
      std::string word;
      // read notation
      file >> word;
      MaterialTensorNotation notation = ((word == "voigt") ? VOIGT : HILL_MANDEL);
      // read a
      unsigned int a_count;
      file >> a_count;
      hom_rect_a_ = Vector<double>(a_count);
      for (unsigned int i = 0; i < a_count; i++) {
        file >> hom_rect_a_[i];
      }
      // read b
      unsigned int b_count;
      file >> b_count;
      hom_rect_b_ = Vector<double>(b_count);
      for (unsigned int i = 0; i < b_count; i++) {
        file >> hom_rect_b_[i];
      }
      // calculate number of coefficient rows
      const unsigned int count = (a_count-1) * (b_count-1);
      // create and fill matrices
      Matrix<double>* matrices[] = {&hom_rect_coeff11_, &hom_rect_coeff12_,
                                    &hom_rect_coeff22_, &hom_rect_coeff33_};
      for (unsigned int k = 0; k < 4; k++) {
        *matrices[k] = Matrix<double>(count, 16);
        for (unsigned int i = 0; i < count; i++) {
          for (unsigned int j = 0; j < 16; j++) {
            file >> (*matrices[k])[i][j];
          }
        }
      }
      // apply notation
      hom_rect_coeff33_ *= (notation == HILL_MANDEL) ? 0.5 : 1.0;
      file.close();
      if (type_ == HOM_ISO_C1) {
        throw Exception("Not implemented for non orthotropic tensors!");
      }
    // full C1 interpolation with text coefficients for 3D (with shearing angle)
    } else if (interpolation_str == "c1_text_3d") {
      interpolation_ = C1;
      // read coefficients from text file
      std::string filename = hr->Get("file")->As<std::string>();
      std::ifstream file(filename.c_str());
      std::string word;
      // read notation
      file >> word;
      MaterialTensorNotation notation = ((word == "voigt") ? VOIGT : HILL_MANDEL);
      // read a
      unsigned int a_count;
      file >> a_count;
      hom_rect_a_ = Vector<double>(a_count);
      for (unsigned int i = 0; i < a_count; i++) {
        file >> hom_rect_a_[i];
      }
      // read b
      unsigned int b_count;
      file >> b_count;
      hom_rect_b_ = Vector<double>(b_count);
      for (unsigned int i = 0; i < b_count; i++) {
        file >> hom_rect_b_[i];
      }
      // read c
      unsigned int c_count;
      file >> c_count;
      hom_rect_c_ = Vector<double>(c_count);
      for (unsigned int i = 0; i < c_count; i++) {
        file >> hom_rect_c_[i];
      }
      // calculate number of coefficient rows
      const unsigned int count = (a_count-1) * (b_count-1) * (c_count-1);
      // create and fill matrices
      Matrix<double>* matrices[] = {&hom_rect_coeff11_, &hom_rect_coeff12_, &hom_rect_coeff13_,
                                    &hom_rect_coeff22_, &hom_rect_coeff23_, &hom_rect_coeff33_};
      for (unsigned int k = 0; k < 6; k++) {
        *matrices[k] = Matrix<double>(count, 64);
        for (unsigned int i = 0; i < count; i++) {
          for (unsigned int j = 0; j < 64; j++) {
            file >> (*matrices[k])[i][j];
          }
        }
      }
      // apply notation
      hom_rect_coeff33_ *= (notation == HILL_MANDEL) ? 0.5 : 1.0;
      file.close();
      if (type_ == HOM_ISO_C1) {
        throw Exception("Not implemented for non orthotropic tensors!");
      }
    // full C1 interpolation with XML coefficients
    } else if (interpolation_str == "c1") {

      bool is_non_ortho = (type_ == HOM_ISO_C1) ? true : false;

      interpolation_ = C1;
      PtrParamNode root = XmlReader::ParseFile(file);
      MaterialTensorNotation notation = tensorNotation.Parse(root->Get("notation")->As<string>());

      // number of sample intervals
      int num_intervals = root->Get("coeff11/matrix/dim1")->As<int>();
      // number of interpolation coefficients per interval
      int num_interp_coeffs = root->Get("coeff11/matrix/dim2")->As<int>();

      // samples for parameters
      assert(root->Has("param1"));
      ParamTools::AsVector<double>(root->Get("param1/matrix/real"), hom_rect_a_);
      if (root->Has("param2"))
        ParamTools::AsVector<double>(root->Get("param2/matrix/real"), hom_rect_b_);
      else {
        hom_rect_b_.Resize(0);
        hom_rect_b_.Init();
      }
      if (root->Has("param3"))
        ParamTools::AsVector<double>(root->Get("param3/matrix/real"), hom_rect_c_);
      else {
        hom_rect_c_.Resize(0);
        hom_rect_c_.Init();
      }

      ReadCoeff(root, "coeff11", num_intervals, num_interp_coeffs, hom_rect_coeff11_);
      ReadCoeff(root, "coeff12", num_intervals, num_interp_coeffs, hom_rect_coeff12_);
      ReadCoeff(root, "coeff22", num_intervals, num_interp_coeffs, hom_rect_coeff22_);

      // heat tensor has only size of 2 x 2 (2d) and is symmetric
      if (type_ != HEAT || dim == 3) {
        ReadCoeff(root, "coeff13", num_intervals, num_interp_coeffs, hom_rect_coeff13_);
        ReadCoeff(root, "coeff23", num_intervals, num_interp_coeffs, hom_rect_coeff23_);
        ReadCoeff(root, "coeff33", num_intervals, num_interp_coeffs, hom_rect_coeff33_);
      }

      if (dim == 2) {
        hom_rect_coeff13_ *= (notation == HILL_MANDEL) ? 1.0/sqrt(2.0) : 1.0;
        hom_rect_coeff23_ *= (notation == HILL_MANDEL) ? 1.0/sqrt(2.0) : 1.0;
        if (type_ != HEAT)
          hom_rect_coeff33_ *= (notation == HILL_MANDEL) ? 0.5 : 1.0;
      }

      if (dim == 3) {
        // heat tensor has only size of 3 x 3 (3d) and is symmetric
        if (type_ != HEAT){
          ReadCoeff(root, "coeff44", num_intervals, num_interp_coeffs, hom_rect_coeff44_);
          ReadCoeff(root, "coeff55", num_intervals, num_interp_coeffs, hom_rect_coeff55_);
          ReadCoeff(root, "coeff66", num_intervals, num_interp_coeffs, hom_rect_coeff66_);
          hom_rect_coeff44_ *= (notation == HILL_MANDEL) ? 0.5 : 1.0;
          hom_rect_coeff55_ *= (notation == HILL_MANDEL) ? 0.5 : 1.0;
          hom_rect_coeff66_ *= (notation == HILL_MANDEL) ? 0.5 : 1.0;
        }

        // returns true if coefficients for entry 14 are given in catalogue
        is_non_ortho = ReadCoeff(root, "coeff14", num_intervals, num_interp_coeffs, hom_rect_coeff14_);

        if (is_non_ortho) {
          // only necessary for the non-ortho case
          ReadCoeff(root, "coeff15", num_intervals, num_interp_coeffs, hom_rect_coeff15_);
          ReadCoeff(root, "coeff16", num_intervals, num_interp_coeffs, hom_rect_coeff16_);
          ReadCoeff(root, "coeff24", num_intervals, num_interp_coeffs, hom_rect_coeff24_);
          ReadCoeff(root, "coeff25", num_intervals, num_interp_coeffs, hom_rect_coeff25_);
          ReadCoeff(root, "coeff26", num_intervals, num_interp_coeffs, hom_rect_coeff26_);
          ReadCoeff(root, "coeff34", num_intervals, num_interp_coeffs, hom_rect_coeff34_);
          ReadCoeff(root, "coeff35", num_intervals, num_interp_coeffs, hom_rect_coeff35_);
          ReadCoeff(root, "coeff36", num_intervals, num_interp_coeffs, hom_rect_coeff36_);
          ReadCoeff(root, "coeff45", num_intervals, num_interp_coeffs, hom_rect_coeff45_);
          ReadCoeff(root, "coeff46", num_intervals, num_interp_coeffs, hom_rect_coeff46_);
          ReadCoeff(root, "coeff56", num_intervals, num_interp_coeffs, hom_rect_coeff56_);
          hom_rect_coeff14_ *= (notation == HILL_MANDEL) ? 1.0/sqrt(2.0) : 1.0;
          hom_rect_coeff15_ *= (notation == HILL_MANDEL) ? 1.0/sqrt(2.0) : 1.0;
          hom_rect_coeff16_ *= (notation == HILL_MANDEL) ? 1.0/sqrt(2.0) : 1.0;
          hom_rect_coeff24_ *= (notation == HILL_MANDEL) ? 1.0/sqrt(2.0) : 1.0;
          hom_rect_coeff25_ *= (notation == HILL_MANDEL) ? 1.0/sqrt(2.0) : 1.0;
          hom_rect_coeff26_ *= (notation == HILL_MANDEL) ? 1.0/sqrt(2.0) : 1.0;
          hom_rect_coeff34_ *= (notation == HILL_MANDEL) ? 1.0/sqrt(2.0) : 1.0;
          hom_rect_coeff35_ *= (notation == HILL_MANDEL) ? 1.0/sqrt(2.0) : 1.0;
          hom_rect_coeff36_ *= (notation == HILL_MANDEL) ? 1.0/sqrt(2.0) : 1.0;
          hom_rect_coeff45_ *= (notation == HILL_MANDEL) ? 0.5 : 1.0;
          hom_rect_coeff46_ *= (notation == HILL_MANDEL) ? 0.5 : 1.0;
          hom_rect_coeff56_ *= (notation == HILL_MANDEL) ? 0.5 : 1.0;
        }
      }

      // read micro load factor coefficients and threshold above which micro load factors are extrapolated
      bool hasMLF = ReadCoeff(root, "microloadfactor", num_intervals, num_interp_coeffs, hom_mlf_coeff_);

      if (hasMLF)
        extrapolationThreshold_ = hr->Get("mlfExtrapolationThreshold")->As<double>();

      // copy Vector to StdVector
      StdVector<double> a(hom_rect_a_.GetSize()), b(hom_rect_b_.GetSize()), c(hom_rect_c_.GetSize());
      std::copy(hom_rect_a_.GetPointer(), hom_rect_a_.GetPointer() + hom_rect_a_.GetSize(), a.Begin());
      std::copy(hom_rect_b_.GetPointer(), hom_rect_b_.GetPointer() + hom_rect_b_.GetSize(), b.Begin());
      std::copy(hom_rect_c_.GetPointer(), hom_rect_c_.GetPointer() + hom_rect_c_.GetSize(), c.Begin());

      // create interpolators for tensor entries
      interpolator11_ = CreateInterpolator(a, b, c, hom_rect_coeff11_);
      interpolator12_ = CreateInterpolator(a, b, c, hom_rect_coeff12_);
      interpolator13_ = CreateInterpolator(a, b, c, hom_rect_coeff13_);
      interpolator23_ = CreateInterpolator(a, b, c, hom_rect_coeff23_);
      interpolator22_ = CreateInterpolator(a, b, c, hom_rect_coeff22_);
      interpolator33_ = CreateInterpolator(a, b, c, hom_rect_coeff33_);

      if (dim == 3) {
        // heat tensor has only size of 3 x 3 (3d) and is symmetric
        if (type_ != HEAT){
          interpolator44_ = CreateInterpolator(a, b, c, hom_rect_coeff44_);
          interpolator55_ = CreateInterpolator(a, b, c, hom_rect_coeff55_);
          interpolator66_ = CreateInterpolator(a, b, c, hom_rect_coeff66_);
        }
        if (is_non_ortho) {
          // only necessary for the non-ortho case
          interpolator14_ = CreateInterpolator(a, b, c, hom_rect_coeff14_);
          interpolator15_ = CreateInterpolator(a, b, c, hom_rect_coeff15_);
          interpolator16_ = CreateInterpolator(a, b, c, hom_rect_coeff16_);
          interpolator24_ = CreateInterpolator(a, b, c, hom_rect_coeff24_);
          interpolator25_ = CreateInterpolator(a, b, c, hom_rect_coeff25_);
          interpolator26_ = CreateInterpolator(a, b, c, hom_rect_coeff26_);
          interpolator34_ = CreateInterpolator(a, b, c, hom_rect_coeff34_);
          interpolator35_ = CreateInterpolator(a, b, c, hom_rect_coeff35_);
          interpolator36_ = CreateInterpolator(a, b, c, hom_rect_coeff36_);
          interpolator45_ = CreateInterpolator(a, b, c, hom_rect_coeff45_);
          interpolator46_ = CreateInterpolator(a, b, c, hom_rect_coeff46_);
          interpolator56_ = CreateInterpolator(a, b, c, hom_rect_coeff56_);
        }
      }

      // create interpolator for micro load factor
      if (hasMLF)
        interpolatorMLF_ = CreateInterpolator(a, b, c, hom_mlf_coeff_);

      LOG_DBG3(dm) << "a = " << hom_rect_a_;
      LOG_DBG3(dm) << "b = " << hom_rect_b_;
      LOG_DBG3(dm) << "c = " << hom_rect_c_;
      LOG_DBG3(dm) << "Size of coeff11 = " << hom_rect_coeff11_.GetNumRows() << " x "<< hom_rect_coeff11_.GetNumCols();
      LOG_DBG3(dm) << "Size of coeff12 = " << hom_rect_coeff12_.GetNumRows() << " x "<< hom_rect_coeff12_.GetNumCols();
      LOG_DBG3(dm) << "Size of coeff13 = " << hom_rect_coeff13_.GetNumRows() << " x "<< hom_rect_coeff13_.GetNumCols();
      LOG_DBG3(dm) << "Size of coeff22 = " << hom_rect_coeff22_.GetNumRows() << " x "<< hom_rect_coeff22_.GetNumCols();
      LOG_DBG3(dm) << "Size of coeff23 = " << hom_rect_coeff23_.GetNumRows() << " x "<< hom_rect_coeff23_.GetNumCols();
      LOG_DBG3(dm) << "Size of coeff33 = " << hom_rect_coeff33_.GetNumRows() << " x "<< hom_rect_coeff33_.GetNumCols();
    } else {

#ifdef USE_SGPP
      // sparse grid interpolation
      if (interpolation_str == "sgpp") {
        interpolation_ = SG;
        unsigned int dimension = (shearIsDesign_ ? 3 : 2);
        std::string basis_str = hr->Get("sgppBasis")->As<std::string>();
        // sparse grid basis to be used
        if ((basis_str == "bspline") || (basis_str == "modbspline")) {
          // B-spline basis
          bspline_degree_ = hr->Get("bsplineDegree")->As<unsigned int>();
          if (basis_str == "bspline") {
            // B-splines
            sgpp_basis_ = BSPLINE;
            grid_ = sgpp::base::Grid::createBsplineGrid(dimension, bspline_degree_);
          } else {
            // modified B-splines
            sgpp_basis_ = MODBSPLINE;
            grid_ = sgpp::base::Grid::createModBsplineGrid(dimension, bspline_degree_);
          }
          // optional verbosity for initial hierarchization
        } else if (basis_str == "modlinear") {
          // modified linear basis
          sgpp_basis_ = MODLINEAR;
          grid_ = sgpp::base::Grid::createModLinearGrid(dimension);
        } else {
          // (un)modified linear basis
          sgpp_basis_ = LINEAR;
          grid_ = sgpp::base::Grid::createLinearGrid(dimension);
        }
        // read coefficients
        InitializeSparseGrid(file.c_str());
        LOG_DBG(dm)<<"Sparse grid initialized and filled.";

      // full B-spline interpolation with XML coefficients
      } else if (interpolation_str == "full_bspline") {
        interpolation_ = FULL_BSPLINE;
        bspline_degree_ = hr->Get("bsplineDegree")->As<unsigned int>();
        PtrParamNode root = XmlReader::ParseFile(file);
        MaterialTensorNotation notation = tensorNotation.Parse(root->Get("notation")->As<string>());
        // read coefficients from XML
        ParamTools::AsMatrix<double>(root->Get("coeff11/matrix"), full_bspline_coeff11_);
        ParamTools::AsMatrix<double>(root->Get("coeff12/matrix"), full_bspline_coeff12_);
        if (grid_->getStorage().getDimension() == 3) {
          ParamTools::AsMatrix<double>(root->Get("coeff13/matrix"), full_bspline_coeff13_);
        }
        ParamTools::AsMatrix<double>(root->Get("coeff22/matrix"), full_bspline_coeff22_);
        if (grid_->getStorage().getDimension() == 3) {
          ParamTools::AsMatrix<double>(root->Get("coeff23/matrix"), full_bspline_coeff23_);
        }
        ParamTools::AsMatrix<double>(root->Get("coeff33/matrix"), full_bspline_coeff33_);
        // apply notation
        full_bspline_coeff33_ *= (notation == HILL_MANDEL) ? 0.5 : 1.0;

      // full B-spline interpolation with text coefficients for 2D (without shearing angle)
      } else if (interpolation_str == "full_bspline_text_2d") {
        interpolation_ = FULL_BSPLINE;
        bspline_degree_ = hr->Get("bsplineDegree")->As<unsigned int>();
        // read coefficients from text file
        std::string filename = hr->Get("file")->As<std::string>();
        std::ifstream file(filename.c_str());
        std::string word;
        // read notation
        file >> word;
        MaterialTensorNotation notation = ((word == "voigt") ? VOIGT : HILL_MANDEL);
        // read number of B-splines
        unsigned int count;
        file >> count;
        // create matrices
        Matrix<double>* matrices[] = {&full_bspline_coeff11_, &full_bspline_coeff12_,
                                      &full_bspline_coeff22_, &full_bspline_coeff33_};
        for (unsigned int k = 0; k < 4; k++) {
          *matrices[k] = Matrix<double>(count, 1);
        }
        // create matrices
        for (unsigned int i = 0; i < count; i++) {
          for (unsigned int k = 0; k < 4; k++) {
            file >> (*matrices[k])[i][0];
          }
        }
        // apply notation
        full_bspline_coeff33_ *= (notation == HILL_MANDEL) ? 0.5 : 1.0;
        file.close();

      // full B-spline interpolation with text coefficients for 3D (with shearing angle)
      } else if (interpolation_str == "full_bspline_text_3d") {
        interpolation_ = FULL_BSPLINE;
        bspline_degree_ = hr->Get("bsplineDegree")->As<unsigned int>();
        // read coefficients from text file
        std::string filename = hr->Get("file")->As<std::string>();
        std::ifstream file(filename.c_str());
        std::string word;
        // read notation
        file >> word;
        MaterialTensorNotation notation = ((word == "voigt") ? VOIGT : HILL_MANDEL);
        // read number of B-splines
        unsigned int count;
        file >> count;
        // create matrices
        Matrix<double>* matrices[] = {&full_bspline_coeff11_, &full_bspline_coeff12_, &full_bspline_coeff13_,
                                      &full_bspline_coeff22_, &full_bspline_coeff23_, &full_bspline_coeff33_};
        for (unsigned int k = 0; k < 6; k++) {
          *matrices[k] = Matrix<double>(count, 1);
        }
        // fill matrices
        for (unsigned int i = 0; i < count; i++) {
          for (unsigned int k = 0; k < 6; k++) {
            file >> (*matrices[k])[i][0];
          }
        }
        // apply notation
        full_bspline_coeff33_ *= (notation == HILL_MANDEL) ? 0.5 : 1.0;
        file.close();
      }
#else //USE_SGPP
    EXCEPTION("CFS is compiled without SGpp toolbox! Please recompile with SGpp or choose another interpolation method.")
#endif //USE_SGPP
    }
  }
  else if (type_ == MSFEM_C1)
  {
    if(!space_->IsRegular())
      throw Exception("MSFEM requires regular design");
    if(dim > 2)
      throw Exception("MSFEM not implemented for 3D");

    // Read interpolation coefficients of MSFEM element stiffness matrices from material catalogue
    PtrParamNode hr = pn->Get("MSFEMC1");
    std::string file = hr->Get("file")->As<std::string>();
    PtrParamNode root = XmlReader::ParseFile(file);
    Matrix<double> msfem_a, msfem_b;
    ParamTools::AsMatrix<double>(root->Get("param1/matrix"), msfem_a);
    msfem_a.GetCol(msfem_a_, 0);
    ParamTools::AsMatrix<double>(root->Get("param2/matrix"), msfem_b);
    msfem_b.GetCol(msfem_b_, 0);
    StdVector<std::string> index;
    index = "11","12","13","14","15","16","17","18","22","23","24","25",
        "26","27","28","33","34","35","36","37","38","44","45","46","47","48","55","56","57","58","66","67","68","77","78","88";
    msfem_coeff_.Resize(36);
    for (int i = 0;i<36;i++) {
      std::stringstream ss;
      ss<<"coeff"<<index[i]<<"/matrix";
      std::string tmp = ss.str();
      ParamTools::AsMatrix<double>(root->Get(tmp), msfem_coeff_[i]);
    }
    LOG_DBG3(dm) << "a = " << msfem_a_;
    LOG_DBG3(dm) << "b = " << msfem_b_;
    LOG_DBG3(dm) << "Size of msfem_coeff = " << msfem_coeff_.GetSize();
  } else {
    LOG_DBG3(dm) << "a = " << hom_rect_a_;
    LOG_DBG3(dm) << "b = " << hom_rect_b_;
    LOG_DBG3(dm) << "c = " << hom_rect_c_;
    LOG_DBG3(dm) << "Size of coeff11 = " << hom_rect_coeff11_.GetNumRows() << " x "<< hom_rect_coeff11_.GetNumCols();
    LOG_DBG3(dm) << "Size of coeff12 = " << hom_rect_coeff12_.GetNumRows() << " x "<< hom_rect_coeff12_.GetNumCols();
    LOG_DBG3(dm) << "Size of coeff13 = " << hom_rect_coeff13_.GetNumRows() << " x "<< hom_rect_coeff13_.GetNumCols();
    LOG_DBG3(dm) << "Size of coeff22 = " << hom_rect_coeff22_.GetNumRows() << " x "<< hom_rect_coeff22_.GetNumCols();
    LOG_DBG3(dm) << "Size of coeff23 = " << hom_rect_coeff23_.GetNumRows() << " x "<< hom_rect_coeff23_.GetNumCols();
    LOG_DBG3(dm) << "Size of coeff33 = " << hom_rect_coeff33_.GetNumRows() << " x "<< hom_rect_coeff33_.GetNumCols();
  }
}


DesignMaterial::DesignMaterial(DesignMaterial::Type type, DesignSpace* space)
{
  space_ = space;
  this->type_ = type;
  dim = domain->GetDim();
}

DesignMaterial::~DesignMaterial() {
  if(interpolator11_) { delete interpolator11_; }
  if(interpolator12_) { delete interpolator12_; }
  if(interpolator13_) { delete interpolator13_; }
  if(interpolator14_) { delete interpolator14_; }
  if(interpolator15_) { delete interpolator15_; }
  if(interpolator16_) { delete interpolator16_; }
  if(interpolator22_) { delete interpolator22_; }
  if(interpolator23_) { delete interpolator23_; }
  if(interpolator24_) { delete interpolator24_; }
  if(interpolator25_) { delete interpolator25_; }
  if(interpolator26_) { delete interpolator26_; }
  if(interpolator33_) { delete interpolator33_; }
  if(interpolator34_) { delete interpolator34_; }
  if(interpolator35_) { delete interpolator35_; }
  if(interpolator36_) { delete interpolator36_; }
  if(interpolator44_) { delete interpolator44_; }
  if(interpolator45_) { delete interpolator45_; }
  if(interpolator46_) { delete interpolator46_; }
  if(interpolator55_) { delete interpolator55_; }
  if(interpolator56_) { delete interpolator56_; }
  if(interpolator66_) { delete interpolator66_; }
  if(interpolatorMLF_) { delete interpolatorMLF_; }
}

ApproxData* DesignMaterial::CreateInterpolator(StdVector<double>& a, StdVector<double>& b, StdVector<double>& c, Matrix<double>& coeff) {
  StdVector<double> dummy_data;
  ApproxData* interpolator = nullptr;

  if (coeff.NormL2() == 0)
    return NULL;

  if (c.GetSize() > 0) {
    assert(a.GetSize() > 0 && b.GetSize() > 0);
    dummy_data.Resize(a.GetSize()*b.GetSize()*c.GetSize());
    interpolator = new TriCubicInterpolate(dummy_data, a, b, c, coeff);
  } else if (b.GetSize() > 0) {
    assert(a.GetSize() > 0);
    dummy_data.Resize(a.GetSize()*b.GetSize());
    interpolator = new BiCubicInterpolate(dummy_data, a, b, coeff);
  } else {
    assert(a.GetSize() > 0);
    dummy_data.Resize(a.GetSize());
    interpolator = new CubicInterpolate(dummy_data, a, coeff);
  }
  return interpolator;
}

inline bool DesignMaterial::ReadCoeff(PtrParamNode pn, const string& name, int nRows, int nCols, Matrix<double>& coeff) const {
  if (pn->Has(name)) {
    ParamTools::AsTensor<double>(pn->Get(name + "/matrix/real"), nRows, nCols, coeff);
    return true;
  } else {
    coeff.Resize(nRows,nCols);
    coeff.Init();
    return false;
  }
}

void DesignMaterial::ToInfo(PtrParamNode in) const
{
  in->Get("type")->SetValue(type.ToString(type_));
  in->Get("bias")->SetValue(bias_); // TODO: Add material names in case

  // dont't easily add GetParameters() map, as there are constants and element temp values mixed
}

void DesignMaterial::FillHomRectSamples(PtrParamNode homRect, unsigned int idx, const string& a, const string& b)
{
  // the internal tensor representation in hom_rect_samples_ is HILL-MANDEL!
  MaterialTensorNotation notation = tensorNotation.Parse(homRect->Get("notation")->As<string>());

  PtrParamNode data = homRect->GetByVal("data", "a", a, "b", b);
  hom_rect_samples_[idx][0] = data->Get("e11")->As<double>();
  hom_rect_samples_[idx][5] = data->Get("e12")->As<double>();
  hom_rect_samples_[idx][1] = data->Get("e22")->As<double>();
  hom_rect_samples_[idx][2] = data->Get("e33")->As<double>() * (notation == HILL_MANDEL ? 0.5 : 1.0);
  hom_rect_samples_[idx][4] = 0.0;
  hom_rect_samples_[idx][3] = 0.0;
}

bool DesignMaterial::CollectMaterialParametersForElement(DesignSpace* space, const Elem* elem)
{
  int base = space->Find(elem, false);
  if(base < 0)
    return false;

  // we must not clear the parameters here as only designs are rewritten but not fixed parameters
  for(unsigned int index = base; index < space->data.GetSize(); index += space->elements)
  {
    DesignElement* de = &space->data[index];
    assert(de->elem->elemNum == elem->elemNum);
    double val = de->GetDesign(DesignElement::SMART);
    LOG_DBG2(dm) << "CMPFE e=" << elem->elemNum << " de=" << de->ToString() << " v=" << val; // << " thread:" << omp_get_thread_num();
    SetParameter(de->GetType(), val, false); // not global element data
  }
  return true;
}

void DesignMaterial::SetParameter(const DesignElement::Type key, const double value, bool global)
{
  if(global)
    for(unsigned int i = 0; i < params_.GetNumSlots(); i++)
      params_.Mine(i)[key] = value;
  else
    params_.Mine()[key] = value;
}

double DesignMaterial::GetParameter(const std::map<DesignElement::Type, double>& map, const DesignElement::Type p)
{
  assert(map.size() == params_.Mine().size());
  assert(HasParameter(p));

  std::map<DesignElement::Type, double>::const_iterator iter = map.find(p);
  assert(iter != map.end());
  return iter->second;
}

const std::map<DesignElement::Type, double>& DesignMaterial::GetParameters() const
{
  return params_.ConstMine();
}

unsigned int DesignMaterial::RequiredParameters( OptimizationMaterial::System material)
{
  unsigned int r = MassIsDesign() ? 1 : 0;
  if (DampingIsDesign()) {
    r += 2;
  }
  switch (type_)
  {
  case FMO:
    assert(material == OptimizationMaterial::MECH || material == OptimizationMaterial::PIEZOCOUPLING);
    if (dim == 2)
      return r + (material == OptimizationMaterial::MECH ? 6 : 15);
    else
      return r + (material == OptimizationMaterial::MECH ? 21 : 15);
  case SGP_GRADIENTCHECK:
  case SGP_MATLAB:
    assert(material == OptimizationMaterial::MECH || material == OptimizationMaterial::PIEZOCOUPLING);
    if (dim == 2)
      return r + (material == OptimizationMaterial::MECH ? 9 : 15);
    else
      return r + (material == OptimizationMaterial::MECH ? 25 : 15);
  case ISOTROPIC:
  case LAME_ISOTROPIC:
    return r + 2;
  case TRANSVERSAL_ISOTROPIC:
  case TRANSVERSAL_ISOTROPIC_BOXED:
  { //FIXME: not initialized anywhere yet, but this is ugly as hell..
    std::string subType;
    domain->GetParamRoot()->GetValue("sequenceStep/pdeList/mechanic/subType", subType);
    return (subType == "planeStress") ? r + 4 : r + 5;
  }
  case DENSITY_TIMES_TRANSVERSAL_ISOTROPIC:
  case DENSITY_TIMES_TRANSVERSAL_ISOTROPIC_BOXED:
  case LAMINATES:
    return r + 5;
  case D_LAMINATES:
    return r + 6;
  case HOM_RECT:
    return r + 3;
  case HOM_RECT_C1:
    if (dim == 2)
      return r + 4;
    else
      return r + 6;
  case HOM_ISO_C1:
  case HEAT:
    assert(dim == 2 || dim == 3);
    return r + dim;
  case MSFEM_C1:
    if (dim == 2)
      return r + 3;
    else
      return r + 6;
  case D_HOM_RECT:
    return r + 4;
  case DENSITY_TIMES_2D_TENSOR_CONSTANT_TRACE:
  case ORTHOTROPIC:
    return r + 6;
  case DENSITY_TIMES_2D_TENSOR:
  case DENSITY_TIMES_ROTATED_2D_TENSOR:
  case DENSITY_TIMES_ORTHOTROPIC:
    return r + 7;
  case DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC:
  case DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED:
    return r + 5 + (dim == 3 ? 3 : 1);
  case DENSITY_TIMES_ROT_PA12:
    return r + 5 + (dim == 3 ? 3 : 1);
  case D_INTERP_IN718_TENSOR:
    return r + 2;
  case D_INTERP_IN718_TENSOR_ROT:
    return r + 3 + (dim == 3 ? 3 : 1);
  case FEATURE_MAPPING_ANISO:
  case NO_TYPE:
    return 0;
  }

  assert(false);
  return 0;
}

bool DesignMaterial::CheckRequiredDesigns(
    StdVector<DesignElement::Type>& design) {
  if (MassIsDesign() && design.Find(DesignElement::MASS) < 0) {
    return (false);
  }
  if (DampingIsDesign()
      && (design.Find(DesignElement::DAMPINGALPHA) < 0
          || design.Find(DesignElement::DAMPINGBETA) < 0)) {
    return (false);
  }

  switch (type_) {
  case FMO:
  case SGP_MATLAB:
  case SGP_GRADIENTCHECK:
    if (dim == 2) {
      return (design.Find(DesignElement::MECH_11) >= 0
          && design.Find(DesignElement::MECH_22) >= 0
          && design.Find(DesignElement::MECH_33) >= 0
          && design.Find(DesignElement::MECH_23) >= 0
          && design.Find(DesignElement::MECH_13) >= 0
          && design.Find(DesignElement::MECH_12) >= 0);
    } else {
      return (design.Find(DesignElement::MECH_11) >= 0
          && design.Find(DesignElement::MECH_12) >= 0
          && design.Find(DesignElement::MECH_13) >= 0
          && design.Find(DesignElement::MECH_14) >= 0
          && design.Find(DesignElement::MECH_15) >= 0
          && design.Find(DesignElement::MECH_16) >= 0
          && design.Find(DesignElement::MECH_22) >= 0
          && design.Find(DesignElement::MECH_23) >= 0
          && design.Find(DesignElement::MECH_24) >= 0
          && design.Find(DesignElement::MECH_25) >= 0
          && design.Find(DesignElement::MECH_26) >= 0
          && design.Find(DesignElement::MECH_33) >= 0
          && design.Find(DesignElement::MECH_34) >= 0
          && design.Find(DesignElement::MECH_35) >= 0
          && design.Find(DesignElement::MECH_36) >= 0
          && design.Find(DesignElement::MECH_44) >= 0
          && design.Find(DesignElement::MECH_45) >= 0
          && design.Find(DesignElement::MECH_46) >= 0
          && design.Find(DesignElement::MECH_55) >= 0
          && design.Find(DesignElement::MECH_56) >= 0
          && design.Find(DesignElement::MECH_66) >= 0);
    }
  case ISOTROPIC:
    return (design.Find(DesignElement::EMODUL) >= 0
        && design.Find(DesignElement::POISSON) >= 0);
  case LAME_ISOTROPIC:
    return (design.Find(DesignElement::LAMELAMBDA) >= 0
        && design.Find(DesignElement::LAMEMU) >= 0);
  case TRANSVERSAL_ISOTROPIC:
  case TRANSVERSAL_ISOTROPIC_BOXED:
  { //FIXME: not initialized anywhere yet, but this is ugly as hell..
    std::string subType;
    domain->GetParamRoot()->GetValue("sequenceStep/pdeList/mechanic/subType", subType);
    return (
        design.Find(DesignElement::EMODULISO) >= 0
            && design.Find(DesignElement::EMODUL) >= 0
            && design.Find(DesignElement::POISSON) >= 0
            && design.Find(DesignElement::GMODUL) >= 0
            && (subType == "planeStress") ? true : design.Find(DesignElement::POISSONISO) >= 0);
  }
  case DENSITY_TIMES_TRANSVERSAL_ISOTROPIC:
  case DENSITY_TIMES_TRANSVERSAL_ISOTROPIC_BOXED:
    return (design.Find(DesignElement::DENSITY) >= 0
        && design.Find(DesignElement::EMODULISO) >= 0
        && design.Find(DesignElement::EMODUL) >= 0
        && design.Find(DesignElement::POISSON) >= 0
        && design.Find(DesignElement::GMODUL) >= 0);
  case DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC:
  case DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED:
    return (design.Find(DesignElement::DENSITY) >= 0
        && design.Find(DesignElement::EMODULISO) >= 0
        && design.Find(DesignElement::EMODUL) >= 0
        && design.Find(DesignElement::POISSON) >= 0
        && design.Find(DesignElement::GMODUL) >= 0
        && (dim != 2 || design.Find(DesignElement::ROTANGLE) >= 0)
        && (dim != 3 || (design.Find(DesignElement::POISSONISO) >= 0
                        && design.Find(DesignElement::ROTANGLEFIRST) >= 0
                        && design.Find(DesignElement::ROTANGLESECOND) >= 0
                        && design.Find(DesignElement::ROTANGLETHIRD) >= 0) ) );
  case DENSITY_TIMES_ROT_PA12:
    return (design.Find(DesignElement::DENSITY) >= 0
        && design.Find(DesignElement::EMODULISO) >= 0
        && design.Find(DesignElement::POISSON) >= 0
        && design.Find(DesignElement::GMODUL) >= 0
        && (dim != 2 || design.Find(DesignElement::ROTANGLE) >= 0)
        && (dim != 3 || (design.Find(DesignElement::ROTANGLEFIRST) >= 0
                        && design.Find(DesignElement::ROTANGLESECOND) >= 0
                        && design.Find(DesignElement::ROTANGLETHIRD) >= 0) ) );
  case D_INTERP_IN718_TENSOR:
    return (design.Find(DesignElement::DENSITY) >= 0
        && design.Find(DesignElement::INTERPOLATION) >=0);
  case D_INTERP_IN718_TENSOR_ROT:
    return (design.Find(DesignElement::DENSITY) >= 0
        && design.Find(DesignElement::EMODUL) >=0
        && design.Find(DesignElement::INTERPOLATION) >=0
        && (dim != 2 || design.Find(DesignElement::ROTANGLE) >= 0)
        && (dim != 3 || (design.Find(DesignElement::ROTANGLEFIRST) >= 0
                        && design.Find(DesignElement::ROTANGLESECOND) >= 0
                        && design.Find(DesignElement::ROTANGLETHIRD) >= 0) ) );
  case ORTHOTROPIC:
    return(design.Find(DesignElement::MECH_11) >= 0
        && design.Find(DesignElement::MECH_22) >= 0
        && design.Find(DesignElement::MECH_33) >= 0
        && design.Find(DesignElement::MECH_12) >= 0
        && design.Find(DesignElement::ROTANGLE) >= 0
        && design.Find(DesignElement::LOWER_EIG_BOUND) >= 0);
  case DENSITY_TIMES_ORTHOTROPIC:
    return(design.Find(DesignElement::DENSITY) >= 0
        && design.Find(DesignElement::MECH_11) >= 0
        && design.Find(DesignElement::MECH_22) >= 0
        && design.Find(DesignElement::MECH_33) >= 0
        && design.Find(DesignElement::MECH_12) >= 0
        && design.Find(DesignElement::ROTANGLE) >= 0
        && design.Find(DesignElement::LOWER_EIG_BOUND) >= 0);
  case DENSITY_TIMES_2D_TENSOR:
    return (design.Find(DesignElement::DENSITY) >= 0
        && design.Find(DesignElement::MECH_11) >= 0
        && design.Find(DesignElement::MECH_22) >= 0
        && design.Find(DesignElement::MECH_33) >= 0
        && design.Find(DesignElement::MECH_23) >= 0
        && design.Find(DesignElement::MECH_13) >= 0
        && design.Find(DesignElement::MECH_12) >= 0);
  case DENSITY_TIMES_2D_TENSOR_CONSTANT_TRACE:
    return (design.Find(DesignElement::DENSITY) >= 0
        && design.Find(DesignElement::MECH_11) >= 0
        && design.Find(DesignElement::MECH_22) >= 0
        && design.Find(DesignElement::MECH_23) >= 0
        && design.Find(DesignElement::MECH_13) >= 0
        && design.Find(DesignElement::MECH_12) >= 0);
  case DENSITY_TIMES_ROTATED_2D_TENSOR:
    return (design.Find(DesignElement::DENSITY) >= 0
        && design.Find(DesignElement::MECH_11) >= 0
        && design.Find(DesignElement::MECH_33) >= 0
        && design.Find(DesignElement::MECH_23) >= 0
        && design.Find(DesignElement::MECH_13) >= 0
        && design.Find(DesignElement::MECH_12) >= 0
        && design.Find(DesignElement::ROTANGLE) >= 0);
  case LAMINATES:
    return (design.Find(DesignElement::STIFF1) >= 0
        && design.Find(DesignElement::STIFF2) >= 0
        && design.Find(DesignElement::ROTANGLE) >= 0
        && design.Find(DesignElement::EMODUL) >= 0
        && design.Find(DesignElement::POISSON) >= 0);
  case D_LAMINATES:
    return(design.Find(DesignElement::STIFF1) >= 0
        && design.Find(DesignElement::STIFF2) >= 0
        && design.Find(DesignElement::DENSITY) >= 0
        && design.Find(DesignElement::ROTANGLE) >= 0
        && design.Find(DesignElement::EMODUL) >= 0
        && design.Find(DesignElement::POISSON) >= 0);
  case HOM_RECT:
    return (design.Find(DesignElement::STIFF1) >= 0
        && design.Find(DesignElement::STIFF2) >= 0
        && design.Find(DesignElement::ROTANGLE) >= 0);
  case D_HOM_RECT:
    return(design.Find(DesignElement::STIFF1) >= 0
            && design.Find(DesignElement::STIFF2) >= 0
            && design.Find(DesignElement::ROTANGLE) >= 0
            && design.Find(DesignElement::DENSITY) >= 0);
  case HOM_RECT_C1:
    if (dim == 3) {
      return (design.Find(DesignElement::STIFF1) >= 0
          && design.Find(DesignElement::STIFF2) >= 0
          && design.Find(DesignElement::STIFF3) >= 0
          && design.Find(DesignElement::ROTANGLETHIRD) >= 0
          && design.Find(DesignElement::ROTANGLESECOND) >= 0);
    } else {
      return (design.Find(DesignElement::STIFF1) >= 0
          && design.Find(DesignElement::STIFF2) >= 0
          && design.Find(DesignElement::SHEAR1) >= 0
          && design.Find(DesignElement::ROTANGLE) >= 0);
    }
  case HOM_ISO_C1:
    if (dim == 3) {
      return (design.Find(DesignElement::STIFF1) >= 0
          && design.Find(DesignElement::ROTANGLETHIRD) >= 0
          && design.Find(DesignElement::ROTANGLESECOND) >= 0);
    } else {
      return (design.Find(DesignElement::STIFF1) >= 0
          && design.Find(DesignElement::ROTANGLE) >= 0);
    }
  case MSFEM_C1:
    if (dim == 3) {
      return (design.Find(DesignElement::STIFF1) >= 0
          && design.Find(DesignElement::STIFF2) >= 0
          && design.Find(DesignElement::STIFF3) >= 0
          && design.Find(DesignElement::ROTANGLETHIRD) >= 0
          && design.Find(DesignElement::ROTANGLESECOND) >= 0);
    } else {
      return (design.Find(DesignElement::STIFF1) >= 0
          && design.Find(DesignElement::STIFF2) >= 0
        && design.Find(DesignElement::ROTANGLE) >= 0);
    }
  case HEAT:
    if (dim == 3) {
      return (design.Find(DesignElement::STIFF1) >= 0
          && design.Find(DesignElement::STIFF2) >= 0
          && design.Find(DesignElement::STIFF3) >= 0);
    } else {
      assert(dim == 2);
      return (design.Find(DesignElement::STIFF1) >= 0
          && design.Find(DesignElement::STIFF2) >= 0);
    }
  case FEATURE_MAPPING_ANISO:  
  case NO_TYPE:
    // we have no param mat variables
    return false;
  }
  assert(false);
  return false;
}

inline void DesignMaterial::GetIsoMaterialTensor(MaterialTensor<double>& mt, SubTensorType subTensor, DesignElement::Type direction) {
  double E  = GetParameter(DesignElement::EMODUL);
  double nu = GetParameter(DesignElement::POISSON);

  Matrix<double>& t = mt.GetMatrix(VOIGT);

  switch (direction) {
    case DesignElement::NO_DERIVATIVE: {
      double lambda = nu * E / ((1.0 + nu) * (1.0 - 2.0 * nu));
      double mu = E / (2.0 * (1.0 + nu));
      double diag = lambda + 2.0 * mu;
      SetIsoMatrix(t, subTensor, diag, lambda, mu);
      break;
    }
    case DesignElement::EMODUL: {
      double dlambda_dE = nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
      double dmu_dE = 1.0 / (2.0 * (1.0 + nu));
      double ddiag_dE = dlambda_dE + 2.0 * dmu_dE;
      SetIsoMatrix(t, subTensor, ddiag_dE, dlambda_dE, dmu_dE);
      break;
    }
    case DesignElement::POISSON: {
      double dlambda_dnu = (1.0 + 2.0 * nu * nu) * E
          / ((1.0 + nu) * (1.0 + nu) * (1.0 - 2.0 * nu) * (1.0 - 2.0 * nu));
      double dmu_dnu = E / (-2.0 * (1.0 + nu) * (1.0 + nu));
      double ddiag_dnu = dlambda_dnu + 2.0 * dmu_dnu;
      SetIsoMatrix(t, subTensor, ddiag_dnu, dlambda_dnu, dmu_dnu);
      break;
    }
    default:
      ZeroMatrix(t, subTensor); // any derivative in any direction other than EMODUL or POISSON is zero
      break;
  }
}

inline double DesignMaterial::GetIsoMaterialMass(DesignElement::Type direction) {
  double E =  GetParameter(DesignElement::EMODUL);
  double nu = GetParameter(DesignElement::POISSON);
  switch (direction) {
  case DesignElement::NO_DERIVATIVE: {
    double lambda = nu * E / ((1.0 + nu) * (1.0 - 2.0 * nu));
    double mu = E / (2.0 * (1.0 + nu));
    double diag = lambda + 2.0 * mu;
    return (GetIsoMass(diag, mu));
  }
  case DesignElement::EMODUL: {
    double dlambda_dE = nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    double dmu_dE = 1.0 / (2.0 * (1.0 + nu));
    double ddiag_dE = dlambda_dE + 2.0 * dmu_dE;
    return (GetIsoMass(ddiag_dE, dmu_dE));
  }
  case DesignElement::POISSON: {
    double dlambda_dnu = (1.0 + 2.0 * nu * nu) * E
        / ((1.0 + nu) * (1.0 + nu) * (1.0 - 2.0 * nu) * (1.0 - 2.0 * nu));
    double dmu_dnu = E / (-2.0 * (1.0 + nu) * (1.0 + nu));
    double ddiag_dnu = dlambda_dnu + 2.0 * dmu_dnu;
    return (GetIsoMass(ddiag_dnu, dmu_dnu));
  }
  default:
    return (0.0); // any derivative in any direction other than EMODUL or POISSON is zero
  }
}

inline void DesignMaterial::GetLameMaterialTensor(MaterialTensor<double>& mt,
    SubTensorType subTensor, DesignElement::Type direction) {

  Matrix<double>& t = mt.GetMatrix(VOIGT);

  switch (direction) {
    case DesignElement::NO_DERIVATIVE: {
      double lambda = GetParameter(DesignElement::LAMELAMBDA);
      double mu = GetParameter(DesignElement::LAMEMU);
      double diag = lambda + 2.0 * mu;
      SetIsoMatrix(t, subTensor, diag, lambda, mu);
      break;
    }
    case DesignElement::LAMELAMBDA:
      SetIsoMatrix(t, subTensor, 1.0, 1.0, 0.0);
      break;
    case DesignElement::LAMEMU:
      SetIsoMatrix(t, subTensor, 2.0, 0.0, 1.0);
      break;
    default:
      ZeroMatrix(t, subTensor); // any derivative in any direction other than EMODUL or POISSON is zero
      break;
    }
}

inline double DesignMaterial::GetLameMaterialMass(DesignElement::Type direction) {
  switch (direction) {
  case DesignElement::NO_DERIVATIVE: {
    double lambda = GetParameter(DesignElement::LAMELAMBDA);
    double mu = GetParameter(DesignElement::LAMEMU);
    double diag = lambda + 2.0 * mu;
    return (GetIsoMass(diag, mu));
  }
  case DesignElement::LAMELAMBDA:
    return (GetIsoMass(1.0, 0.0));
  case DesignElement::LAMEMU:
    return (GetIsoMass(2.0, 1.0));
  default:
    return (0.0); // any derivative in any direction other than EMODUL or POISSON is zero
  }
}

inline void DesignMaterial::GetTransIsoMaterialTensor(MaterialTensor<double>& mt, SubTensorType subTensor, DesignElement::Type direction, bool core)
{
  LOG_DBG2(dm) << "GetTransIsoMaterialTensor called with direction=" << (direction == DesignElement::NO_DERIVATIVE ? "no_derivative" : DesignElement::type.ToString(direction));
  assert(type_ != DENSITY_TIMES_ROT_PA12 || subTensor == FULL);

  Matrix<double>& t = mt.GetMatrix(VOIGT);

  double E3(0.0);
  double E = GetParameter(DesignElement::EMODULISO);
  if (type_ != DENSITY_TIMES_ROT_PA12)
    E3 = GetParameter(DesignElement::EMODUL);
  else
  {
    E3 = 137.4 + 2.4*E;
    E = 145.0 - 5.8*E;
  }
  double G3 = GetParameter(DesignElement::GMODUL);
  double nu13 = GetParameter(DesignElement::POISSON); //used as theta in the boxed formulations

  if (subTensor == PLANE_STRESS) {
    double dens(1.0), factor(1.0), ninv2(0.0), D(0.0), D3(0.0), nD3(0.0);
    if((type_ == DENSITY_TIMES_TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_TRANSVERSAL_ISOTROPIC_BOXED
        || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED))
    {
      dens = GetParameter(DesignElement::DENSITY);
      TransferFunction* tf = space_->GetTransferFunction(DesignElement::DENSITY, App::MECH);
      factor = (direction == DesignElement::DENSITY) ? tf->Derivative(dens) : tf->Transform(dens);
    } else {
      if(direction == DesignElement::DENSITY)
        factor = 0.0;
    }
    if (type_ == TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC) {
      ninv2 = E3 - nu13 * nu13 * E;
      //assert(ninv2>0.0); //positivity of the elasticity tensor is violated. Use constraint "parametrized-plane-stress-pos-def" > 0
      ninv2 = 1 / (ninv2 * ninv2);
    }
    switch (direction) {
    case DesignElement::NO_DERIVATIVE:
    case DesignElement::DENSITY:
    case DesignElement::ROTANGLE: {
      if (type_ == TRANSVERSAL_ISOTROPIC
          || type_ == DENSITY_TIMES_TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC) {
        D = E * E3 / (E3 - nu13 * nu13 * E);
        D3 = E3 * E3 / (E3 - nu13 * nu13 * E);
        nD3 = nu13 * E * E3 / (E3 - nu13 * nu13 * E);
      } else if(type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED || type_ == TRANSVERSAL_ISOTROPIC_BOXED || type_ == DENSITY_TIMES_TRANSVERSAL_ISOTROPIC_BOXED){
        D = E / (1 - nu13);
        D3 = E3 / (1 - nu13);
        nD3 = sqrt(E * E3 * nu13) / (1 - nu13);
      }else{
        throw Exception("Not yet implemented!");
      }
      SetTransIsoMatrix(t, subTensor, factor*D, 0, 0, factor*D3, factor*nD3, factor*G3);
      break;
    }
    case DesignElement::EMODULISO: {
      if (type_ == TRANSVERSAL_ISOTROPIC
          || type_ == DENSITY_TIMES_TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC) {
        D = E3 * E3 * ninv2;
        D3 = nu13 * nu13 * E3 * E3 * ninv2;
        nD3 = nu13 * E3 * E3 * ninv2;
      } else if(type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED || type_ == TRANSVERSAL_ISOTROPIC_BOXED || type_ == DENSITY_TIMES_TRANSVERSAL_ISOTROPIC_BOXED){
        D = 1 / (1 - nu13);
        nD3 = sqrt(E3 * nu13 / E) / (2 - 2 * nu13);
      }else{
        throw Exception("Not yet implemented!");
      }
      SetTransIsoMatrix(t, subTensor, factor * D, 0, 0, factor * D3, factor * nD3, 0);
      break;
    }
    case DesignElement::EMODUL: {
      if (type_ == TRANSVERSAL_ISOTROPIC
          || type_ == DENSITY_TIMES_TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC) {
        D = -E * E * nu13 * nu13 * ninv2;
        D3 = -E3 * (-E3 + 2 * nu13 * nu13 * E) * ninv2;
        nD3 = -nu13 * nu13 * nu13 * E * E * ninv2;
      } else if(type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED || type_ == TRANSVERSAL_ISOTROPIC_BOXED || type_ == DENSITY_TIMES_TRANSVERSAL_ISOTROPIC_BOXED){
        D3 = 1 / (1 - nu13);
        nD3 = sqrt(E * nu13 / E3) / (2 - 2 * nu13);
      }else{
        throw Exception("Not yet implemented!");
      }
      SetTransIsoMatrix(t, subTensor, factor * D, 0, 0, factor * D3, factor * nD3, 0);
      break;
    }
    case DesignElement::POISSON: {
      if (type_ == TRANSVERSAL_ISOTROPIC
          || type_ == DENSITY_TIMES_TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC) {
        D = 2 * nu13 * E * E * E3 * ninv2;
        D3 = 2 * nu13 * E * E3 * E3 * ninv2;
        nD3 = E * E3 * (nu13 * nu13 * E + E3) * ninv2;
      } else if(type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED || type_ == TRANSVERSAL_ISOTROPIC_BOXED || type_ == DENSITY_TIMES_TRANSVERSAL_ISOTROPIC_BOXED){
        D = 1 / ((1 - nu13) * (1 - nu13));
        D3 = E3 * D;
        nD3 = sqrt(E * E3 / nu13) * (nu13 + 1) * D * 0.5;
        D = E * D;
      }else{
        throw Exception("Not yet implemented!");
      }
      SetTransIsoMatrix(t, subTensor, factor * D, 0, 0, factor * D3, factor * nD3, 0);
      break;
    }
    case DesignElement::GMODUL: {
      SetTransIsoMatrix(t, subTensor, 0, 0, 0, 0, 0, factor);
      break;
    }
    default:
      ZeroMatrix(t, subTensor);
      return;
    } // switch direction
    if(type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED){
      assert(!core); // not verified for core = true
      double rotAngle = GetParameter(DesignElement::ROTANGLE);
      LOG_DBG2(dm)<< "GetTransIsoMaterialTensor: E before rotation = " << t.ToString();
      RotateTensor(mt, direction, true, rotAngle);
      LOG_DBG2(dm)<< "GetTransIsoMaterialTensor: E after rotation = " << t.ToString();

      //    static int count(0);
      //    if (count % 10 == 0 && count/100 % 10 == 0){
      ////      std::cout << "(" << (count/100 % 10)*(count % 10)+1 << ")" << t.ToString() << std::endl;
      //      std::cout << t(0,0) << " " << t(0,1) << " " << t(1,1) << " " << t(0,2) << " " << t(1,2) << " " << t(2,2) << std::endl;
      //    }
      //    count++;
    }
    return;
  } // PLANE_STRESS
  
  // 3D AND PLANE_STRAIN
  double nu = GetParameter(DesignElement::POISSONISO);
  double nu3;
  double n3;
  double c;
  double dens = 1.0, factor = 1.0;
  if (!core && (type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED || type_ == DENSITY_TIMES_ROT_PA12)) {
    dens = GetParameter(DesignElement::DENSITY);
    TransferFunction* tf = space_->GetTransferFunction(DesignElement::DENSITY, App::MECH);
    factor = (direction == DesignElement::DENSITY) ? tf->Derivative(dens) : tf->Transform(dens);
  } else {
    if(direction == DesignElement::DENSITY)
      factor = 0.0;
  }
  if (type_ == TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC) {
    nu3 = nu13 * E3/E;
    n3 = nu13*nu3;
    c = (1.0-nu-2.0*n3); // this is the interesting thing, this must not get 0, however this would imply a volume (trace of tensor) of infinity, so it is hopefully not occuring
    if(c < 1e-8) {
      c = 1e-8;
    }
  } else if (type_ == TRANSVERSAL_ISOTROPIC_BOXED || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED || type_ == DENSITY_TIMES_ROT_PA12){
    nu3 = sqrt(0.5*(1.0-nu)*E3/E)*nu13;
    n3 = nu3*nu3*E/E3;
    c = (1.0-nu-2.0*n3);
  } else {
    throw Exception("Not yet implemented!");
  }

  LOG_DBG3(dm) << "GTIMTensor: nu3=" << nu3;

  double f = E / ((1.0 + nu) * c);
  double dE = 0.0, dE3 = 0.0, dnu = 0.0, dnu3 = 0.0, dn3 = 0.0, dG3 = 0.0;
  
  bool tensorset = false;
  switch(direction){
  case DesignElement::NO_DERIVATIVE:
  case DesignElement::DENSITY: // almost the same as no derivative, we only changed the factor above
  case DesignElement::ROTANGLE:
  case DesignElement::ROTANGLETHIRD:
  case DesignElement::ROTANGLESECOND:
  {
    double D = (1.0-n3)*f;
    double D3 = (1.0-nu)*E3/c;
    double nD = (nu+n3)*f;
    double nD3 = (1.0+nu)*nu3*f;
    double G = 0.5*E/(1.0+nu);
    SetTransIsoMatrix(t, subTensor, factor * D, factor * nD, factor * G, factor * D3, factor * nD3, factor*G3);
    tensorset = true;
    break;
  }
  case DesignElement::EMODULISO:
    dE = 1.0;
    if(type_ == TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC){
      dnu3 = -E3*nu13/(E*E);
      dn3 = nu3/E3 * (2.0*E*dnu3 + nu3);
    } else if (type_ == TRANSVERSAL_ISOTROPIC_BOXED || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED) {
      dnu3 = -sqrt(0.125*(1.0-nu)*E3/E)*nu13/E;
    } else if (type_ == DENSITY_TIMES_ROT_PA12) {
      double D =  -(29*((nu - 1)*nu13*nu13 + 2))/(10*(nu*nu - 1)*(nu13*nu13 - 1));
      double D3 = -12/(5*(nu13*nu13 - 1));
      double nD =  (29*(nu*(nu13*nu13 - 2) - nu13*nu13))/(10*(nu*nu - 1)*(nu13*nu13 - 1));
      double nD3 = -(0.3*nu13*(4*E - 954.1))/(E*(nu13*nu13 - 1)*sqrt(E3*.5*(1-nu)/E));
      double G = -29/(10*(nu + 1));
      SetTransIsoMatrix(t, subTensor, factor * D, factor * nD, factor * G, factor * D3, factor * nD3, 0.0);
      tensorset = true;
    }
    break;
  case DesignElement::EMODUL:
    dE3 = 1.0;
    if (type_ == TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC) {
      dnu3 = nu13/E;
      dn3 = nu3*E/E3 * (2.0*dnu3 - nu3/E3);
    } else if (type_ == TRANSVERSAL_ISOTROPIC_BOXED || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED) {
      dnu3 = sqrt(0.125*(1.0-nu)/(E*E3))*nu13;
    } else if(type_ == DENSITY_TIMES_ROT_PA12) {
      dE3 = 0.0;
    }
    break;
  case DesignElement::POISSONISO:
    dnu = 1.0;
    if (type_ == TRANSVERSAL_ISOTROPIC_BOXED || type_ == DENSITY_TIMES_ROT_PA12) { // else = 0
      dnu3 = -sqrt(0.125 * E3 / (E * (1.0 - nu))) * nu13;
      dn3 = -0.5 * nu13 * nu13;
    }
    break;
  case DesignElement::POISSON:
    if (type_ == TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC) {
      dnu3 = 1.0;
      dn3 = 2.0 * nu3 * E / E3 * dnu3;
    } else if (type_ == TRANSVERSAL_ISOTROPIC_BOXED || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED || type_ == DENSITY_TIMES_ROT_PA12) {
      dnu3 = sqrt(0.5 * (1.0 - nu) * E3 / E);
      dn3 = (1.0 - nu) * nu13;
    }
    break;
  case DesignElement::GMODUL:
    dG3 = 1.0;
    break;
  default:
    ZeroMatrix(t, subTensor);
    return;
  } // switch(direction)

  if(!tensorset){ // several cases are handled already and set the tensor
    double dc = -dnu-2.0*dn3;
    double df = ( dE - E*dnu/(1.0+nu) - E*dc/c ) / ((1.0+nu)*c);
    double dD = (1.0-n3)*df - dn3*f;
    double dnD = (nu+n3)*df + (dnu+dn3)*f;
    double dD3 = ( (1.0-nu)*dE3 - dnu*E3 - (1.0-nu)*E3*dc/c ) / c;
    double dnD3 = (1.0+nu)*nu3*df + (1.0+nu)*dnu3*f + dnu*nu3*f;
    double dG = 0.5 * ( (1.0+nu)*dE - E*dnu ) / ( (1.0+nu)*(1.0+nu) );
    SetTransIsoMatrix(t, subTensor, factor * dD, factor * dnD, factor * dG, factor * dD3, factor * dnD3, factor * dG3);
  }
  
  if(!core && (type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED || type_ == DENSITY_TIMES_ROT_PA12)){
    // for all rotated types, rotate the material tensor
    LOG_DBG2(dm) << "GetTransIsoMaterialTensor: tensor before rotation=" << t.ToString();
    RotateTensor(mt, direction);
    LOG_DBG2(dm)<< "GetTransIsoMaterialTensor: tensor after rotation = " << t.ToString();
  }
  LOG_DBG2(dm) << "GetTransIsoMaterialTensor: tensor result is " << t.ToString();
}


inline double DesignMaterial::GetTransIsoMaterialMass(DesignElement::Type direction)
{
  double E = GetParameter(DesignElement::EMODULISO);
  double E3 = GetParameter(DesignElement::EMODUL);
  double nu = GetParameter(DesignElement::POISSONISO);
  double nu13 = GetParameter(DesignElement::POISSON);
  double nu3 = nu13 * E3 / E;
  double n3 = nu3 * nu3 * E / E3;
  double c = (1.0 - nu - 2.0 * n3); // this is the interesting thing, this must not get 0, however this would imply a volume (trace of tensor) of infinity, so it is hopefully not occuring
  if (type_ == TRANSVERSAL_ISOTROPIC) {
    if (c < 1e-8) {
      c = 1e-8;
    }
  } else {
    nu3 = sqrt(0.5 * (1.0 - nu) * E3 / E) * nu13;
    n3 = nu3 * nu3 * E / E3;
    c = (1.0 - nu - 2.0 * n3);
  }
  double f = E / ((1.0 + nu) * c);
  double dE = 0.0, dE3 = 0.0, dnu = 0.0, dnu3 = 0.0, dn3 = 0.0, dG3 = 0.0;

  switch (direction) {
  case DesignElement::NO_DERIVATIVE: {
    double D = (1.0 - n3) * f;
    double D3 = (1.0 - nu) * E3 / c;
    double G3 = GetParameter(DesignElement::GMODUL);
    double G = 0.5 * E / (1.0 + nu);
    return (GetTransIsoMass(D, G, D3, G3));
  }
  case DesignElement::EMODULISO:
    dE = 1.0;
    if (type_ == TRANSVERSAL_ISOTROPIC) {
      dnu3 = -E3 * nu13 / (E * E);
      dn3 = nu3 / E3 * (2.0 * E * dnu3 + nu3);
    } else {
      dnu3 = -sqrt(0.125 * (1.0 - nu) * E3 / E) * nu13 / E;
    }
    break;
  case DesignElement::EMODUL:
    dE3 = 1.0;
    if (type_ == TRANSVERSAL_ISOTROPIC) {
      dnu3 = nu13 / E;
      dn3 = nu3 * E / E3 * (2.0 * dnu3 - nu3 / E3);
    } else {
      dnu3 = sqrt(0.125 * (1.0 - nu) / (E * E3)) * nu13;
    }
    break;
  case DesignElement::POISSONISO:
    dnu = 1.0;
    if (type_ == TRANSVERSAL_ISOTROPIC_BOXED) { // else = 0
      dnu3 = -sqrt(0.125 * E3 / (E * (1.0 - nu))) * nu13;
      dn3 = -0.5 * nu13 * nu13;
    }
    break;
  case DesignElement::POISSON:
    if (type_ == TRANSVERSAL_ISOTROPIC) {
      dnu3 = 1.0;
      dn3 = 2.0 * nu3 * E / E3 * dnu3;
    } else {
      dnu3 = sqrt(0.5 * (1.0 - nu) * E3 / E);
      dn3 = (1.0 - nu) * nu13;
    }
    break;
  case DesignElement::GMODUL:
    dG3 = 1.0;
    break;
  default:
    return (0.0);
  } // switch(direction)
  double dc = -dnu - 2.0 * dn3;
  double df = (dE - E * dnu / (1.0 + nu) - E * dc / c) / ((1.0 + nu) * c);
  double dD = (1.0 - n3) * df - dn3 * f;
  double dD3 = ((1.0 - nu) * dE3 - dnu * E3 - (1.0 - nu) * E3 * dc / c) / c;
  double dG = 0.5 * ((1.0 + nu) * dE - E * dnu) / ((1.0 + nu) * (1.0 + nu));
  return (GetTransIsoMass(dD, dG, dD3, dG3));
}

inline double DesignMaterial::GetDensityTimesTensorMass(DesignElement::Type direction){
  double dens = GetParameter(DesignElement::DENSITY);
  TransferFunction* tf = space_->GetTransferFunction(DesignElement::DENSITY, App::MECH);
  switch (direction){
  case DesignElement::NO_DERIVATIVE:
  {
    return tf->Transform(dens);
  }
  case DesignElement::DENSITY:
  {
    return tf->Derivative(dens);
  }
  default:
    return 0.0;
  }
}

inline void DesignMaterial::GetOrthotropicMaterialTensor(MaterialTensor<double>& mt, SubTensorType subTensor, DesignElement::Type direction){
  double e11 = GetParameter(DesignElement::MECH_11);
  double e22 = GetParameter(DesignElement::MECH_22);
  double e33 = GetParameter(DesignElement::MECH_33); //This is already Hill-Mandel notation -> scaling for Voigt notation
  double e12 = GetParameter(DesignElement::MECH_12);
  double rotAngle = GetParameter(DesignElement::ROTANGLE);
  double lowerEigBound = GetParameter(DesignElement::LOWER_EIG_BOUND);
  double dens(1.0), factor(1.0);
  if(type_ == DENSITY_TIMES_ORTHOTROPIC){
    dens = GetParameter(DesignElement::DENSITY);
  }
  TransferFunction* tf = space_->GetTransferFunction(DesignElement::DENSITY, App::MECH); // Identity TransferFunction if not defined
  factor = tf->Transform(dens); // true in most cases, otherwise set again

  Matrix<double>& t = mt.GetMatrix(HILL_MANDEL);

  if(subTensor == PLANE_STRESS){ //This is the only implemented case for now
    switch(direction){
    case DesignElement::NO_DERIVATIVE:
    case DesignElement::ROTANGLE:
    {
      SetTransIsoMatrix(t, subTensor, factor*(e11*e11+e12*e12)+lowerEigBound, 0, 0,
          factor*(e12*e12+e22*e22)+lowerEigBound, factor*(e11*e12+e12*e22), factor*e33+lowerEigBound);
      break;
    }
    case DesignElement::DENSITY:
    {
      factor = tf->Derivative(dens);
      SetTransIsoMatrix(t, subTensor, factor*(e11*e11+e12*e12), 0, 0, factor*(e12*e12+e22*e22), factor*(e11*e12+e12*e22), factor*e33);
      break;
    }
    case DesignElement::MECH_11:
    {
      SetTransIsoMatrix(t, subTensor, 2.0*factor*e11, 0, 0, 0, factor*e12, 0);
      break;
    }
    case DesignElement::MECH_22:
    {
      SetTransIsoMatrix(t, subTensor, 0, 0, 0, 2.0*factor*e22, factor*e12, 0);
      break;
    }
    case DesignElement::MECH_12:
    {
      SetTransIsoMatrix(t, subTensor, 2.0*factor*e12, 0, 0, 2.0*factor*e12, factor*(e11+e22), 0);
      break;
    }
    case DesignElement::MECH_33:
    {
      SetTransIsoMatrix(t, subTensor, 0, 0, 0, 0, 0, factor);
      break;
    }
    default:
      ZeroMatrix(t, subTensor);
      return;
    } // switch direction

    LOG_DBG2(dm)<< "GetOrthotropicMaterialTensor: E before rotation = " << t.ToString();
    // transform to Voigt notation for rotation
    mt.ToVoigt();
    RotateTensor(mt, direction, true, rotAngle);
    LOG_DBG2(dm)<< "GetOrthotropicMaterialTensor: E after rotation = " << t.ToString();

    return;
  } // PLANE_STRESS
  else
    throw Exception("subTensor not implemented yet");
}

inline void DesignMaterial::GetDensityTimes2dTensorTensor(MaterialTensor<double>& mt, SubTensorType subTensor, DesignElement::Type direction, bool core)
{
  // DumpParams();
  double e11 = 0;
  double e22 = 0;
  double e33 = 0;
  double e23 = 0;
  double e13 = 0;
  double e12 = 0;
  if (direction == DesignElement::NO_DERIVATIVE
      || direction == DesignElement::DENSITY
      || direction == DesignElement::ROTANGLE) {
    e11 = GetParameter(DesignElement::MECH_11);
    if (type_ == DENSITY_TIMES_2D_TENSOR_CONSTANT_TRACE) {
      e22 = GetParameter(DesignElement::MECH_22);
      e33 = 0.5 * (trace_ - e11 - e22);
    }
    // TODO: this was hard coded for
    //bvu: commented out was it does not make sense for me!!!
    // else if (type_ == DENSITY_TIMES_ROTATED_2D_TENSOR) {
    //      e22 = 15 - e11;
    //      e11 += 1.0;
    //      e33 = GetParameter(DesignElement::MECH_33);}
    else{
      e22 = GetParameter(DesignElement::MECH_22);
      e33 = GetParameter(DesignElement::MECH_33);
    }
    e23 = GetParameter(DesignElement::MECH_23);
    e13 = GetParameter(DesignElement::MECH_13);
    e12 = GetParameter(DesignElement::MECH_12);
  }

  Matrix<double>& t = mt.GetMatrix(VOIGT);
  switch (direction) {
    case DesignElement::NO_DERIVATIVE:
    case DesignElement::ROTANGLE:
    case DesignElement::DENSITY: // Treated after switch
      Set2dMatrix(t, e11, e22, e33, e23, e13, e12);
      break;
    case DesignElement::MECH_11:
      Set2dMatrix(t, 1.0, type_ == DENSITY_TIMES_ROTATED_2D_TENSOR ? -1.0 : 0.0,
          type_ == DENSITY_TIMES_2D_TENSOR_CONSTANT_TRACE ? -0.5 : 0.0, 0.0,
          0.0, 0.0);
      break;
    case DesignElement::MECH_22:
      Set2dMatrix(t, 0.0, 1.0,
          type_ == DENSITY_TIMES_2D_TENSOR_CONSTANT_TRACE ? -0.5 : 0.0, 0.0,
          0.0, 0.0);
      break;
    case DesignElement::MECH_33:
      Set2dMatrix(t, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0);
      break;
    case DesignElement::MECH_23:
      Set2dMatrix(t, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
      break;
    case DesignElement::MECH_13:
      Set2dMatrix(t, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
      break;
    case DesignElement::MECH_12:
      Set2dMatrix(t, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
      break;
    default:
      ZeroMatrix(t, subTensor);
      break;
  }

  // we're only interested in the tensor of core material
  if (core)
    return;

  if (type_ == DENSITY_TIMES_ROTATED_2D_TENSOR) {
    double rotAngle = GetParameter(DesignElement::ROTANGLE);
    LOG_DBG2(dm)<< "GetDensityTimes2dTensorTensor: E before rotation = " << t.ToString();
    RotateTensor(mt, direction, true, rotAngle);
    LOG_DBG2(dm)<< "GetDensityTimes2dTensorTensor: E after rotation = " << t.ToString();
//    static int count(0);
//    if (count % 10 == 0 && count/100 % 10 == 0){
////      std::cout << "(" << (count/100 % 10)*(count % 10)+1 << ")" << t.ToString() << std::endl;
//      std::cout << t(0,0) << " " << t(0,1) << " " << t(1,1) << " " << t(0,2) << " " << t(1,2) << " " << t(2,2) << std::endl;
//    }
//    count++;
  }
  double dens = GetParameter(DesignElement::DENSITY);
  // standard (SIMP) case: we want product of core tensor and transformed pseudo-density
  TransferFunction* tf = space_->GetTransferFunction(DesignElement::DENSITY, App::MECH);
  t *= (direction == DesignElement::DENSITY) ? tf->Derivative(dens) : tf->Transform(dens);
}

inline void DesignMaterial::GetElasticFMOTensor(MaterialTensor<double>& mt, SubTensorType subTensor, DesignElement::Type direction)
{
  // We use the anisotropic tensor only for solving FMO problems. We assume the design to be in Hill-Mandel
  // notation and therefore we need to transform it for using it in CFS
  const std::map<DesignElement::Type, double>& map = GetParameters();

  bool set = direction == DesignElement::NO_DERIVATIVE; //|| direction == DesignElement::ROTANGLE;

  double e11 = set ? GetParameter(map, DesignElement::MECH_11) : 0;
  double e12 = set ? GetParameter(map, DesignElement::MECH_12) : 0;
  double e13 = set ? GetParameter(map, DesignElement::MECH_13) : 0;
  double e22 = set ? GetParameter(map, DesignElement::MECH_22) : 0;
  double e23 = set ? GetParameter(map, DesignElement::MECH_23) : 0;
  double e33 = set ? GetParameter(map, DesignElement::MECH_33) : 0;

  double e14 = 0, e15 = 0, e16 = 0, e24 = 0, e25 = 0, e26 = 0, e34 = 0, e35 = 0, e36 = 0, e44 = 0, e45 = 0, e46 = 0, e55 = 0, e56 = 0, e66 = 0;
  if (subTensor == FULL) {
    // 3D tensor
    e14 = set ? GetParameter(map, DesignElement::MECH_14) : 0;
    e15 = set ? GetParameter(map, DesignElement::MECH_15) : 0;
    e16 = set ? GetParameter(map, DesignElement::MECH_16) : 0;
    e24 = set ? GetParameter(map, DesignElement::MECH_24) : 0;
    e25 = set ? GetParameter(map, DesignElement::MECH_25) : 0;
    e26 = set ? GetParameter(map, DesignElement::MECH_26) : 0;
    e34 = set ? GetParameter(map, DesignElement::MECH_34) : 0;
    e35 = set ? GetParameter(map, DesignElement::MECH_35) : 0;
    e36 = set ? GetParameter(map, DesignElement::MECH_36) : 0;
    e44 = set ? GetParameter(map, DesignElement::MECH_44) : 0;
    e45 = set ? GetParameter(map, DesignElement::MECH_45) : 0;
    e46 = set ? GetParameter(map, DesignElement::MECH_46) : 0;
    e55 = set ? GetParameter(map, DesignElement::MECH_55) : 0;
    e56 = set ? GetParameter(map, DesignElement::MECH_56) : 0;
    e66 = set ? GetParameter(map, DesignElement::MECH_66) : 0;
  }

  // we need Hill-Mandel notation
  // exception: in SGP_GRADIENTCHECK case, take the tensor as it is specified in the xml file
  // if we are called by ApplyPhysicalDesign, mt is in Voigt notation.
  // if we are called by IntegrateDesignVariable, mt is in Hill-Mandel notation.
  if (mt.GetNotation() == VOIGT && type_ != SGP_GRADIENTCHECK && space_->GetMethod() != ErsatzMaterial::SPAGHETTI_PARAM_MAT)
  {
    // for ToHillMandel(), mt.matrix_ has to be initialized
    mt.GetMatrix(VOIGT).Resize(subTensor == FULL? 6 : 3, subTensor == FULL? 6 : 3);
    mt.GetMatrix(VOIGT).Init();
    mt.ToHillMandel();
  }


  MaterialTensorNotation fmo_notation = HILL_MANDEL;
  // SGP_GRADIENTCHECK: we want the tensor entries, as they are specified in the xml file (without any further transformation!)

  if (space_->GetMethod() == ErsatzMaterial::SPAGHETTI_PARAM_MAT || type_ == SGP_GRADIENTCHECK)
    fmo_notation = VOIGT;

  Matrix<double>& E = mt.GetMatrix(fmo_notation);

  switch (direction) {
  case DesignElement::NO_DERIVATIVE:
  //case DesignElement::ROTANGLE:
    if (subTensor != FULL) {
      Set2dMatrix(E, e11, e22, e33, e23, e13, e12);
    } else {
      // temporarily: orthotropic 3D tensor
      Set3dMatrix(E,subTensor, e11,e12,e13,e14,e15,e16,e22,e23,e24,e25,e26,e33,e34,e35,e36,e44,e45,e46,e55,e56,e66);
    }
    break;
  case DesignElement::MECH_11:
    if (subTensor != FULL) {
      Set2dMatrix(E, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    } else {
      Set3dMatrix(E,subTensor, 1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.);
    }
    break;
  case DesignElement::MECH_12:
    if (subTensor != FULL) {
      Set2dMatrix(E, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    } else {
      Set3dMatrix(E,subTensor, 0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.);
    }
    break;
  case DesignElement::MECH_13:
    if (subTensor != FULL) {
      Set2dMatrix(E, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    } else {
      Set3dMatrix(E,subTensor, 0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.);
    }
    break;
  case DesignElement::MECH_14:
    assert(subTensor == FULL);
    Set3dMatrix(E,subTensor, 0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.);
    break;
  case DesignElement::MECH_15:
    assert(subTensor == FULL);
    Set3dMatrix(E,subTensor, 0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.);
    break;
  case DesignElement::MECH_16:
    assert(subTensor == FULL);
    Set3dMatrix(E,subTensor, 0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.);
    break;
  case DesignElement::MECH_22:
    if (subTensor != FULL) {
      Set2dMatrix(E, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0);
    } else {
      Set3dMatrix(E,subTensor, 0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.);
    }
    break;
  case DesignElement::MECH_23:
    if (subTensor != FULL) {
      Set2dMatrix(E, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    } else {
      Set3dMatrix(E,subTensor, 0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.);
    }
    break;
  case DesignElement::MECH_24:
    assert(subTensor == FULL);
    Set3dMatrix(E,subTensor, 0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.);
    break;
  case DesignElement::MECH_25:
    assert(subTensor == FULL);
    Set3dMatrix(E,subTensor, 0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.);
    break;
  case DesignElement::MECH_26:
    assert(subTensor == FULL);
    Set3dMatrix(E,subTensor, 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.);
    break;
  case DesignElement::MECH_33:
    if (subTensor != FULL) {
      Set2dMatrix(E, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0);
    } else {
      Set3dMatrix(E,subTensor, 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.);
    }
    break;
  case DesignElement::MECH_34:
    assert(subTensor == FULL);
    Set3dMatrix(E,subTensor,   0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.);
    break;
  case DesignElement::MECH_35:
    assert(subTensor == FULL);
    Set3dMatrix(E,subTensor,   0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.);
    break;
  case DesignElement::MECH_36:
    assert(subTensor == FULL);
    Set3dMatrix(E,subTensor,   0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.);
    break;
  case DesignElement::MECH_44:
    assert(subTensor == FULL);
    Set3dMatrix(E,subTensor,   0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.);
    break;
  case DesignElement::MECH_45:
    assert(subTensor == FULL);
    Set3dMatrix(E,subTensor,   0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.);
    break;
  case DesignElement::MECH_46:
    assert(subTensor == FULL);
    Set3dMatrix(E,subTensor,   0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.);
    break;
  case DesignElement::MECH_55:
    assert(subTensor == FULL);
    Set3dMatrix(E,subTensor,   0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.);
    break;
  case DesignElement::MECH_56:
    assert(subTensor == FULL);
    Set3dMatrix(E,subTensor,   0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.);
    break;
  case DesignElement::MECH_66:
    assert(subTensor == FULL);
    Set3dMatrix(E,subTensor, 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.);
    break;
  default:
    // for piezo FMO the derivative w.r.t. dielec_11, ... is zero
    ZeroMatrix(E, subTensor);
    break;
  }

  if (type_ != SGP_GRADIENTCHECK || fmo_notation == HILL_MANDEL)
    mt.ToVoigt();

  LOG_DBG2(dm) << "GEFMOT: E  =  " << mt.GetMatrix(VOIGT).ToString();
}

inline void DesignMaterial::GetInterpolatedHomTensor(MaterialTensor<double>& mt, SubTensorType subTensor, const Elem* elem, DesignElement::Type direction)
{
  const std::map<DesignElement::Type, double>& map = GetParameters();

  // Get design variables
  double a,b,c;
  if (type_ == HOM_ISO_C1) {
    a = GetParameter(map, DesignElement::STIFF1);
    b = a;
    c = a;
  } else {
    a = GetParameter(map, DesignElement::STIFF1);
    b = GetParameter(map, DesignElement::STIFF2);
    // iff subTensor == FULL we have 3D
    c = GetParameter(map, subTensor == FULL ? DesignElement::STIFF3 : DesignElement::SHEAR1);
  }
  double rotAngle = subTensor != FULL ? GetParameter(map, DesignElement::ROTANGLE) : 0;

  LocPoint p;
  p.coord.Resize(3);
  if (type_ == HOM_RECT || type_ == D_HOM_RECT) {
    p.coord[0] = -1.0 + 4 * a; // assume max 0.5
    p.coord[1] = -1.0 + 4 * b; // assume max 0.5
  }
  if (type_ == HOM_RECT_C1 || type_ == HOM_ISO_C1 || type_ == HEAT) {
    p.coord[0] = a;
    p.coord[1] = b;
    p.coord[2] = c;
  }

/* #ifndef NDEBUG
  Vector<double> peps(p.coord);
  double eps = 1e-8;
  Matrix<double> Eeps(E);
  Matrix<double> Etmp(E);
  Vector<double> Diff(3);
 #endif
 */
  LOG_DBG2(dm) << "GetInterpolatedHomTensor: dir=" << (direction == DesignElement::NO_DERIVATIVE ? "no_derivative" : DesignElement::type.ToString(direction))
               << " not=" << mt.GetNotation() << " rotAngle=" << rotAngle << " a=" << a << " b=" << b <<" c="<<(subTensor == FULL ? c : 0.0)<< " -> " << p.coord.ToString();

  // only relevant for HOM_RECT and D_HOM_RECT
  FeH1LagrangeQuad9 fe;

  // need this for default switch case and logging
  Matrix<double>& E = mt.GetMatrix(VOIGT);
  switch (direction) {
    case DesignElement::NO_DERIVATIVE:
    case DesignElement::ROTANGLE:
    case DesignElement::ROTANGLETHIRD:
    case DesignElement::ROTANGLESECOND:
    case DesignElement::DENSITY:
      if (type_ == HOM_RECT || type_ == D_HOM_RECT) {
        Vector<double> shape;
        fe.GetShFnc(shape, p, elem);
        GetHomRectTensor(mt, shape);
        LOG_DBG2(dm)<< "GetInterpolatedHomTensor: shape=" << shape.ToString();
      }
      if (type_ == HOM_RECT_C1 || type_ == HEAT) {
        if (interpolation_ == C1) {
          GetHomC1Tensor(mt,p.coord,direction,subTensor);
        }
  #ifdef USE_SGPP
          else if (interpolation_ == SG) {
          GetHomRectSGPPTensor(mt,p.coord,direction,subTensor);
        } else if (interpolation_ == FULL_BSPLINE) {
          GetHomRectFullBsplineTensor(mt,p.coord,direction,subTensor);
        }
  #endif //USE_SGPP
      }
      if (type_ == HOM_ISO_C1) {
        if (interpolation_ == C1) {
          GetHomIsoC1Tensor(mt,p.coord,direction,subTensor);
        }
      }
      break;
    case DesignElement::STIFF1:
    case DesignElement::STIFF2:
    case DesignElement::STIFF3:
    case DesignElement::SHEAR1:
      if(type_ == HOM_RECT || type_ == D_HOM_RECT) {
        Matrix<double>& jac = fe.GetLocDerivShFnc(p, elem);
        LOG_DBG3(dm) << "GetInterpolatedHomTensor: jac=" << jac.ToString();
        Vector<double> d_shape;
        jac.GetCol(d_shape, direction == DesignElement::STIFF1 ? 0 : 1);// a or by
        GetHomRectTensor(mt, d_shape);
        // correct scaling to local FE coordinates - why??
        mt.GetMatrix(VOIGT) *= 4;
        LOG_DBG2(dm) << "GetInterpolatedHomTensor: d_shape=" << d_shape.ToString();
      } else if(type_ == HOM_RECT_C1 || type_ == HEAT) {
        if (interpolation_ == C1) {
          GetHomC1Tensor(mt,p.coord,direction,subTensor);
        }
  #ifdef USE_SGPP
          else if (interpolation_ == SG) {
          GetHomRectSGPPTensor(mt,p.coord,direction,subTensor);
        } else if (interpolation_ == FULL_BSPLINE) {
          GetHomRectFullBsplineTensor(mt,p.coord,direction,subTensor);
        }
  #endif //USE_SGPP
      } else if (type_ == HOM_ISO_C1) {
        if (interpolation_ == C1) {
          GetHomIsoC1Tensor(mt, p.coord, direction, subTensor);
        }
      }
      break;
    default:
      ZeroMatrix(E, subTensor);
  }

  if (type_ != HEAT) {
    LOG_DBG2(dm)<< "GetInterpolatedHomTensor: E before rotation = " << E.ToString();
    RotateTensor(mt, direction);
    LOG_DBG2(dm)<< "GetInterpolatedHomTensor: E after rotation =  " << E.ToString();
  }

  if (type_ == D_HOM_RECT)
  {
    double dens = GetParameter(map, DesignElement::DENSITY);
    TransferFunction* tf = space_->GetTransferFunction(DesignElement::DENSITY, App::MECH);
    E *= (direction == DesignElement::DENSITY) ? tf->Derivative(dens) : tf->Transform(dens);
  }
  /*   for(double y = 0; y <= 0.5; y += 0.25)
   {
   for(double x = 0; x <= 0.5; x += 0.25 )
   {
   p[0] = -1.0 + 4 * x; //
   p[1] = -1.0 + 4 * y;
   fe.GetShFnc(shape, p, NULL);
   hom_rect_samples_.GetCol(data, DesignElement::TENSOR11 - DesignElement::TENSOR11);
   assert(shape.GetSize() == data.GetSize());
   double val = shape * data;
   std::cout << "x=" << x << " y=" << y << " xi=" << p[0] << " eta=" << p[1] << " -> " << val << std::endl; //" s=" << shape.ToString() << " d=" << data.ToString() << std::endl;
   }
   }*/
}


double DesignMaterial::GetMicroLoadFactor(double vol, DesignElement::Type direction)
{
  // do a quadratic extrapolation for all x > x0
  double x0 = extrapolationThreshold_;
  double ev0 = interpolatorMLF_->EvaluateFunc(x0);
  double dev0 = interpolatorMLF_->EvaluateDeriv(x0);
  double ddev0 = interpolatorMLF_->EvaluateSecondDeriv(x0);

  double ev = -1;
  if(direction == DesignElement::NO_DERIVATIVE)
    if(vol < x0)
      ev = interpolatorMLF_->EvaluateFunc(vol);
    else
      ev = ddev0/2 * std::pow(vol-x0, 2) + dev0 * (vol-x0) + ev0;
  else
    if(vol < x0)
      ev = interpolatorMLF_->EvaluateDeriv(vol);
    else
      ev = ddev0 * (vol-x0) + dev0;

  return ev;
}


void DesignMaterial::GetHomRectTensor(MaterialTensor<double>& mt, const Vector<double>& shape) const
{
  Matrix<double>& E = mt.GetMatrix(VOIGT);

  E.Resize(3, 3);
  E.Init(); // for off-diagonal
  Vector<double> data;
  hom_rect_samples_.GetCol(data, 0);
  
  E[1 - 1][1 - 1] = shape * data;
  LOG_DBG(dm)<< "AHRT 11=" << E[1-1][1-1] << " data=" << data.ToString();
  hom_rect_samples_.GetCol(data, 5);
  E[1 - 1][2 - 1] = shape * data;
  E[2 - 1][1 - 1] = E[1 - 1][2 - 1];
  LOG_DBG(dm)<< "AHRT 12=" << E[1-1][2-1] << " data=" << data.ToString();
  hom_rect_samples_.GetCol(data, 1);
  E[2 - 1][2 - 1] = shape * data;
  LOG_DBG(dm)<< "AHRT 22=" << E[2-1][2-1] << " data=" << data.ToString();
  hom_rect_samples_.GetCol(data, 2);
  E[3 - 1][3 - 1] = shape * data;
  LOG_DBG(dm)<< "AHRT 33=" << E[3-1][3-1] << " data=" << data.ToString();
}

void DesignMaterial::GetHomC1Tensor(MaterialTensor<double>& mt, Vector<double>& p,
    DesignElement::Type direction, SubTensorType subTensor) {

  // assume mech material
  int nrows = dim == 2? 3: 6;
  if (type_ == HEAT)
    nrows = dim;

  Matrix<double>& E = mt.GetMatrix(VOIGT);

  E.Resize(nrows, nrows);
  E.Init(); // for off-diagonal
  if (subTensor == FULL) {
    // Calculation of the interpolated tensor values
    if (direction == DesignElement::NO_DERIVATIVE || direction == DesignElement::ROTANGLE
        || direction == DesignElement::ROTANGLETHIRD || direction == DesignElement::ROTANGLESECOND) {
      E[0][0] = interpolator11_->EvaluateFunc(p);
      E[0][1] = interpolator12_->EvaluateFunc(p);
      E[1][0] = E[0][1];
      if (interpolator13_ != NULL) {
        E[0][2] = interpolator13_->EvaluateFunc(p);
        E[2][0] = E[0][2];
      }
      E[1][1] = interpolator22_->EvaluateFunc(p);
      if (interpolator23_ != NULL) {
        E[1][2] = interpolator23_->EvaluateFunc(p);
        E[2][1] = E[1][2];
      }
      E[2][2] = interpolator33_->EvaluateFunc(p);
      LOG_DBG(dm) << "E11= " << E[0][0] << " E12= " << E[0][1] << " E13= " << E[0][2] << " E22= " << E[1][1] << " E23= " << E[1][2]<< " E33= " << E[2][2];
      if (type_ != HEAT) {
        E[3][3] = interpolator44_->EvaluateFunc(p);
        E[4][4] = interpolator55_->EvaluateFunc(p);
        E[5][5] = interpolator66_->EvaluateFunc(p);
        LOG_DBG(dm) << " E44= "<<E[3][3]<<" E55= "<<E[4][4]<<" E66= "<<E[5][5];
        E[0][3] = interpolator14_ ? interpolator14_->EvaluateFunc(p) : 0.0;
        E[0][4] = interpolator15_ ? interpolator15_->EvaluateFunc(p) : 0.0;
        E[0][5] = interpolator16_ ? interpolator16_->EvaluateFunc(p) : 0.0;
        E[1][3] = interpolator24_ ? interpolator24_->EvaluateFunc(p) : 0.0;
        E[1][4] = interpolator25_ ? interpolator25_->EvaluateFunc(p) : 0.0;
        E[1][5] = interpolator26_ ? interpolator26_->EvaluateFunc(p) : 0.0;
        E[2][3] = interpolator34_ ? interpolator34_->EvaluateFunc(p) : 0.0;
        E[2][4] = interpolator35_ ? interpolator35_->EvaluateFunc(p) : 0.0;
        E[2][5] = interpolator36_ ? interpolator36_->EvaluateFunc(p) : 0.0;
        E[3][4] = interpolator45_ ? interpolator45_->EvaluateFunc(p) : 0.0;
        E[3][5] = interpolator46_ ? interpolator46_->EvaluateFunc(p) : 0.0;
        E[4][5] = interpolator56_ ? interpolator56_->EvaluateFunc(p) : 0.0;
        E[3][0] = E[0][3];
        E[4][0] = E[0][4];
        E[5][0] = E[0][5];
        E[3][1] = E[1][3];
        E[4][1] = E[1][4];
        E[5][1] = E[1][5];
        E[3][2] = E[2][3];
        E[4][2] = E[2][4];
        E[5][2] = E[2][5];
        E[4][3] = E[3][4];
        E[5][3] = E[3][5];
        E[5][4] = E[4][5];
      }
    } else {
      // Calculation of the interpolated tensor derivatives
      E[0][0] = EvaluateC1Interpolation_Deriv(p, interpolator11_, direction);
      E[0][1] = EvaluateC1Interpolation_Deriv(p, interpolator12_, direction);
      E[1][0] = E[0][1];
      if (interpolator13_ != NULL) {
        E[0][2] = EvaluateC1Interpolation_Deriv(p, interpolator13_, direction);
        E[2][0] = E[0][2];
      }
      E[1][1] = EvaluateC1Interpolation_Deriv(p, interpolator22_, direction);
      if (interpolator23_ != NULL) {
        E[1][2] = EvaluateC1Interpolation_Deriv(p, interpolator23_, direction);
        E[2][1] = E[1][2];
      }
      E[2][2] = EvaluateC1Interpolation_Deriv(p, interpolator33_, direction);
      LOG_DBG(dm) << "Derivative " << ((direction == DesignElement::STIFF1)?"1":((direction == DesignElement::STIFF2) ? "2":"3")) << " E11= " << E[0][0] << " E12= " << E[0][1] << " E13= " << E[0][2] << " E22= " << E[1][1] << " E33= " << E[2][2] << " E23= " << E[1][2];
      if (type_ != HEAT) {
        E[3][3] = EvaluateC1Interpolation_Deriv(p, interpolator44_, direction);
        E[4][4] = EvaluateC1Interpolation_Deriv(p, interpolator55_, direction);
        E[5][5] = EvaluateC1Interpolation_Deriv(p, interpolator66_, direction);
        LOG_DBG(dm)<<" E44= "<<E[3][3]<<" E55= "<<E[4][4]<<" E66= "<<E[5][5];

        E[0][3] = EvaluateC1Interpolation_Deriv(p, interpolator14_, direction);
        E[0][4] = EvaluateC1Interpolation_Deriv(p, interpolator15_, direction);
        E[0][5] = EvaluateC1Interpolation_Deriv(p, interpolator16_, direction);
        E[1][3] = EvaluateC1Interpolation_Deriv(p, interpolator24_, direction);
        E[1][4] = EvaluateC1Interpolation_Deriv(p, interpolator25_, direction);
        E[1][5] = EvaluateC1Interpolation_Deriv(p, interpolator26_, direction);
        E[2][3] = EvaluateC1Interpolation_Deriv(p, interpolator34_, direction);
        E[2][4] = EvaluateC1Interpolation_Deriv(p, interpolator35_, direction);
        E[2][5] = EvaluateC1Interpolation_Deriv(p, interpolator36_, direction);
        E[3][4] = EvaluateC1Interpolation_Deriv(p, interpolator45_, direction);
        E[3][5] = EvaluateC1Interpolation_Deriv(p, interpolator46_, direction);
        E[4][5] = EvaluateC1Interpolation_Deriv(p, interpolator56_, direction);
        E[3][0] = E[0][3];
        E[4][0] = E[0][4];
        E[5][0] = E[0][5];
        E[3][1] = E[1][3];
        E[4][1] = E[1][4];
        E[5][1] = E[1][5];
        E[3][2] = E[2][3];
        E[4][2] = E[2][4];
        E[5][2] = E[2][5];
        E[4][3] = E[3][4];
        E[5][3] = E[3][5];
        E[5][4] = E[4][5];
      }
    }
  } else {
    if (direction == DesignElement::NO_DERIVATIVE || direction == DesignElement::ROTANGLE
        || direction == DesignElement::ROTANGLETHIRD || direction == DesignElement::ROTANGLESECOND) {
      E[0][0] = interpolator11_->EvaluateFunc(p);
      E[0][1] = interpolator12_->EvaluateFunc(p);
      E[1][1] = interpolator22_->EvaluateFunc(p);
      E[1][0] = E[0][1];
      if (type_ != HEAT) {
        if (interpolator13_ != NULL) {
          E[0][2] = interpolator13_->EvaluateFunc(p);
          E[2][0] = E[0][2];
        }
        if (interpolator23_ != NULL) {
          E[1][2] = interpolator23_->EvaluateFunc(p);
          E[2][1] = E[1][2];
        }
        E[2][2] = interpolator33_->EvaluateFunc(p);
      }
      LOG_DBG(dm) << p.ToString();
      LOG_DBG(dm) << "E11= "<<E[0][0]<<" E12= "<<E[0][1]<<" E13= "<<E[0][2]<<" E22= "<< E[1][1]<<" E23= "<<E[1][2]<<" E33= "<<E[2][2];
    } else {
      // derivatives
      assert(dim == 2);
      E[0][0] = EvaluateC1Interpolation_Deriv(p, interpolator11_, direction);
      E[0][1] = EvaluateC1Interpolation_Deriv(p, interpolator12_, direction);
      E[1][1] = EvaluateC1Interpolation_Deriv(p, interpolator22_, direction);
      E[1][0] = E[0][1];
      if (type_ != HEAT) {
        if (interpolator13_ != NULL) {
          E[0][2] = EvaluateC1Interpolation_Deriv(p, interpolator13_, direction);
          E[2][0] = E[0][2];
        }
        if (interpolator23_ != NULL) {
          E[1][2] = EvaluateC1Interpolation_Deriv(p, interpolator23_, direction);
          E[2][1] = E[1][2];
        }
        E[2][2] = EvaluateC1Interpolation_Deriv(p, interpolator33_, direction);
      }
      LOG_DBG(dm) << p.ToString();
      double E_33 = type_ == HEAT? 0 : E[2][2];
      LOG_DBG(dm) << "Derivative "<<((direction == DesignElement::STIFF1)?"1":(direction == DesignElement::STIFF2)?"2":"3")
          <<" E11= "<<E[0][0]<<" E12= "<<E[0][1]<<" E13= "<<E[0][2]<<" E22= "<< E[1][1]<<" E23= "<<E[1][2]<<" E33= "<<E_33;
    }
  }

}

void DesignMaterial::GetHomIsoC1Tensor(MaterialTensor<double>& mt, Vector<double>& p,
    DesignElement::Type direction, SubTensorType subTensor) const {

  assert(type_ == HOM_ISO_C1);
  assert(dynamic_cast<CubicInterpolate*>(interpolator11_) != NULL); // check type

  Matrix<double>& E = mt.GetMatrix(VOIGT);

  if (subTensor == FULL) {
    E.Resize(6, 6);
    E.Init(); // for off-diagonal
    // Calculation of the interpolated tensor values
    if (direction == DesignElement::NO_DERIVATIVE || direction == DesignElement::ROTANGLE
        || direction == DesignElement::ROTANGLETHIRD || direction == DesignElement::ROTANGLESECOND) {
      E[1 - 1][1 - 1] = interpolator11_->EvaluateFunc(p);
      E[1 - 1][2 - 1] = interpolator12_->EvaluateFunc(p);
      E[1 - 1][3 - 1] = interpolator13_->EvaluateFunc(p);
      E[1 - 1][4 - 1] = interpolator14_->EvaluateFunc(p);
      E[1 - 1][5 - 1] = interpolator15_->EvaluateFunc(p);
      E[1 - 1][6 - 1] = interpolator16_->EvaluateFunc(p);
      E[2 - 1][2 - 1] = interpolator22_->EvaluateFunc(p);
      E[2 - 1][3 - 1] = interpolator23_->EvaluateFunc(p);
      E[2 - 1][4 - 1] = interpolator24_->EvaluateFunc(p);
      E[2 - 1][5 - 1] = interpolator25_->EvaluateFunc(p);
      E[2 - 1][6 - 1] = interpolator26_->EvaluateFunc(p);
      E[3 - 1][3 - 1] = interpolator33_->EvaluateFunc(p);
      E[3 - 1][4 - 1] = interpolator34_->EvaluateFunc(p);
      E[3 - 1][5 - 1] = interpolator35_->EvaluateFunc(p);
      E[3 - 1][6 - 1] = interpolator36_->EvaluateFunc(p);
      E[4 - 1][4 - 1] = interpolator44_->EvaluateFunc(p);
      E[4 - 1][5 - 1] = interpolator45_->EvaluateFunc(p);
      E[4 - 1][6 - 1] = interpolator46_->EvaluateFunc(p);
      E[5 - 1][5 - 1] = interpolator55_->EvaluateFunc(p);
      E[5 - 1][6 - 1] = interpolator56_->EvaluateFunc(p);
      E[6 - 1][6 - 1] = interpolator66_->EvaluateFunc(p);
      E[2 - 1][1 - 1] = E[1 - 1][2 - 1];
      E[3 - 1][1 - 1] = E[1 - 1][3 - 1];
      E[4 - 1][1 - 1] = E[1 - 1][4 - 1];
      E[5 - 1][1 - 1] = E[1 - 1][5 - 1];
      E[6 - 1][1 - 1] = E[1 - 1][6 - 1];
      E[3 - 1][2 - 1] = E[2 - 1][3 - 1];
      E[4 - 1][2 - 1] = E[2 - 1][4 - 1];
      E[5 - 1][2 - 1] = E[2 - 1][5 - 1];
      E[6 - 1][2 - 1] = E[2 - 1][6 - 1];
      E[4 - 1][3 - 1] = E[3 - 1][4 - 1];
      E[5 - 1][3 - 1] = E[3 - 1][5 - 1];
      E[6 - 1][3 - 1] = E[3 - 1][6 - 1];
      E[5 - 1][4 - 1] = E[4 - 1][5 - 1];
      E[6 - 1][4 - 1] = E[4 - 1][6 - 1];
      E[6 - 1][5 - 1] = E[5 - 1][6 - 1];
      LOG_DBG(dm)<<"E11= "<<E[0][0]<<" E12= "<<E[0][1]<<" E22= "<< E[1][1]<<" E33= "<<E[2][2]<<" E23= "<<E[1][2]<<" E13= "<<E[0][2]<<" E44= "<<E[3][3]<<" E55= "<<E[4][4]<<" E66= "<<E[5][5];
    } else {
      // Calculation of the interpolated tensor derivatives
      E[1 - 1][1 - 1] = EvaluateC1Interpolation_Deriv(p, interpolator11_, direction);
      E[1 - 1][2 - 1] = EvaluateC1Interpolation_Deriv(p, interpolator12_, direction);
      E[1 - 1][3 - 1] = EvaluateC1Interpolation_Deriv(p, interpolator13_, direction);
      E[1 - 1][4 - 1] = EvaluateC1Interpolation_Deriv(p, interpolator14_, direction);
      E[1 - 1][5 - 1] = EvaluateC1Interpolation_Deriv(p, interpolator15_, direction);
      E[1 - 1][6 - 1] = EvaluateC1Interpolation_Deriv(p, interpolator16_, direction);
      E[2 - 1][2 - 1] = EvaluateC1Interpolation_Deriv(p, interpolator22_, direction);
      E[2 - 1][3 - 1] = EvaluateC1Interpolation_Deriv(p, interpolator23_, direction);
      E[2 - 1][4 - 1] = EvaluateC1Interpolation_Deriv(p, interpolator24_, direction);
      E[2 - 1][5 - 1] = EvaluateC1Interpolation_Deriv(p, interpolator25_, direction);
      E[2 - 1][6 - 1] = EvaluateC1Interpolation_Deriv(p, interpolator26_, direction);
      E[3 - 1][3 - 1] = EvaluateC1Interpolation_Deriv(p, interpolator33_, direction);
      E[3 - 1][4 - 1] = EvaluateC1Interpolation_Deriv(p, interpolator34_, direction);
      E[3 - 1][5 - 1] = EvaluateC1Interpolation_Deriv(p, interpolator35_, direction);
      E[3 - 1][6 - 1] = EvaluateC1Interpolation_Deriv(p, interpolator36_, direction);
      E[4 - 1][4 - 1] = EvaluateC1Interpolation_Deriv(p, interpolator44_, direction);
      E[4 - 1][5 - 1] = EvaluateC1Interpolation_Deriv(p, interpolator45_, direction);
      E[4 - 1][6 - 1] = EvaluateC1Interpolation_Deriv(p, interpolator46_, direction);
      E[5 - 1][5 - 1] = EvaluateC1Interpolation_Deriv(p, interpolator55_, direction);
      E[5 - 1][6 - 1] = EvaluateC1Interpolation_Deriv(p, interpolator56_, direction);
      E[6 - 1][6 - 1] = EvaluateC1Interpolation_Deriv(p, interpolator66_, direction);
      E[2 - 1][1 - 1] = E[1 - 1][2 - 1];
      E[3 - 1][1 - 1] = E[1 - 1][3 - 1];
      E[4 - 1][1 - 1] = E[1 - 1][4 - 1];
      E[5 - 1][1 - 1] = E[1 - 1][5 - 1];
      E[6 - 1][1 - 1] = E[1 - 1][6 - 1];
      E[3 - 1][2 - 1] = E[2 - 1][3 - 1];
      E[4 - 1][2 - 1] = E[2 - 1][4 - 1];
      E[5 - 1][2 - 1] = E[2 - 1][5 - 1];
      E[6 - 1][2 - 1] = E[2 - 1][6 - 1];
      E[4 - 1][3 - 1] = E[3 - 1][4 - 1];
      E[5 - 1][3 - 1] = E[3 - 1][5 - 1];
      E[6 - 1][3 - 1] = E[3 - 1][6 - 1];
      E[5 - 1][4 - 1] = E[4 - 1][5 - 1];
      E[6 - 1][4 - 1] = E[4 - 1][6 - 1];
      E[6 - 1][5 - 1] = E[5 - 1][6 - 1];
      LOG_DBG(dm)<<"Derivative "<<"1"<<" E11= "<<E[0][0]<<" E12= "<<E[0][1]<<" E22= "<< E[1][1]<<" E33= "<<E[2][2]<<" E23= "<<E[1][2]<<" E13= "<<E[0][2]<<" E44= "<<E[3][3]<<" E55= "<<E[4][4]<<" E66= "<<E[5][5];
    }
  } else {
    E.Resize(3,3);
    E.Init();
    if (direction == DesignElement::NO_DERIVATIVE || direction == DesignElement::ROTANGLE
        || direction == DesignElement::ROTANGLETHIRD || direction == DesignElement::ROTANGLESECOND) {
      E[1-1][1-1] = interpolator11_->EvaluateFunc(p);
      E[1-1][2-1] = interpolator12_->EvaluateFunc(p);
      E[1-1][3-1] = interpolator13_->EvaluateFunc(p);
      E[2-1][2-1] = interpolator22_->EvaluateFunc(p);
      E[2-1][3-1] = interpolator23_->EvaluateFunc(p);
      E[3-1][3-1] = interpolator33_->EvaluateFunc(p);
      E[2-1][1-1] = E[1-1][2-1];
      E[3-1][1-1] = E[1-1][3-1];
      E[3-1][2-1] = E[2-1][3-1];
      LOG_DBG(dm)<<p;
      LOG_DBG(dm)<<"E11= "<<E[0][0]<<" E12= "<<E[0][1]<<" E13= "<<E[0][2]<<" E22= "<< E[1][1]<<" E23= "<<E[1][2]<<" E33= "<<E[2][2];
    } else {
      E[1-1][1-1] = EvaluateC1Interpolation_Deriv(p, interpolator11_, direction);
      E[1-1][2-1] = EvaluateC1Interpolation_Deriv(p, interpolator12_, direction);
      E[1-1][3-1] = EvaluateC1Interpolation_Deriv(p, interpolator13_, direction);
      E[2-1][2-1] = EvaluateC1Interpolation_Deriv(p, interpolator22_, direction);
      E[2-1][3-1] = EvaluateC1Interpolation_Deriv(p, interpolator23_, direction);
      E[3-1][3-1] = EvaluateC1Interpolation_Deriv(p, interpolator33_, direction);
      E[2-1][1-1] = E[1-1][2-1];
      E[3-1][1-1] = E[1-1][3-1];
      E[3-1][2-1] = E[2-1][3-1];
      LOG_DBG(dm)<<p;
      LOG_DBG(dm)<<"Derivative "<<"1"
          <<" E11= "<<E[0][0]<<" E12= "<<E[0][1]<<" E13= "<<E[0][2]<<" E22= "<< E[1][1]<<" E23= "<<E[1][2]<<" E33= "<<E[2][2];
    }
  }

}// end function GetHomIsoC1Tensor

#ifdef USE_SGPP

void DesignMaterial::InitializeSparseGrid(const char * filename) {
  std::ifstream file(filename);
  if (!file) {
    EXCEPTION("File " << filename << " does not exist.\n");
  }
  std::string word;
  // read first word
  file >> word;
  if (word != "sparsegrid") {
    // ==> old format
    file.close();
    bool dataIsSparse;
    Matrix<double> data;
    dataIsSparse = ReadDetailedStats(filename, data);
    // create regular grid
    sgpp::base::GridGenerator& gridGen = grid_->getGenerator();
    gridGen.regular(level_);
    if (dataIsSparse) {
      FillSparseGridWithSparseGridData(data);
    } else {
      FillSparseGridWithFullGridData(data);
    }
    return;
  }

  // new format
  unsigned int dim = grid_->getStorage().getDimension();
  unsigned int N, d, m;
  file >> N >> d >> m >> word;
  // standard assumptions
  assert((d == 2) || (d == 3));
  assert((m == 5) || (m == 6) || (m == 7));
  // logically equivalent to "(shearIsDesign_) ==> ((d=3) and (m=7))"
  assert(!shearIsDesign_ || ((d == 3) && ((m == 6) || (m == 7))));
  bool hierarchized = (word == "hierarchized");
  file >> word;
  MaterialTensorNotation notation = ((word == "voigt") ? VOIGT : HILL_MANDEL);
  // initialize coefficient vectors
  alpha1_.resize(N);
  alpha2_.resize(N);
  alpha3_.resize(N);
  alpha4_.resize(N);
  if (dim == 3) {
    alpha5_.resize(N);
    alpha6_.resize(N);
  }
  alpha7_.resize(N);
  double leveld;
  double indexd;
  std::vector<unsigned int> level(d, 0);
  std::vector<unsigned int> index(d, 0);
  sgpp::base::GridIndex grid_point(dim);
  level_ = 0; // we set that to the maximal level of the grid points
  double duck; // dummy
  unsigned int j = 0;
  for (unsigned int i = 0; i < N; i++) {
    if (dim == 3) {
      // shearing angle should be optimized ==> Read ALL The Data!
      for (unsigned int idx = 0; idx < 3; idx++) {
        file >> leveld >> indexd;
        level[idx] = static_cast<unsigned int>(leveld);
        index[idx] = static_cast<unsigned int>(indexd);
        grid_point.set(idx, level[idx], index[idx]);
      }
      level_ = std::max(level_, std::max(level[0], std::max(level[1], level[2])));
    } else {
      // shearing angle should not be optimized ==> maybe we have to pick the right data
      if (d == 2) {
        for (unsigned int idx = 0; idx < 2; idx++) {
          file >> leveld >> indexd;
          level[idx] = static_cast<unsigned int>(leveld);
          index[idx] = static_cast<unsigned int>(indexd);
          grid_point.set(idx, level[idx], index[idx]);
        }
      } else {
        for (unsigned int idx = 0; idx < 3; idx++) {
          file >> leveld >> indexd;
          level[idx] = static_cast<unsigned int>(leveld);
          index[idx] = static_cast<unsigned int>(indexd);
        }
        if (level[2] != 1) {
          // data is 3D, but this data point is not in the relevant z=0.5 plane ==> skip
          for (unsigned int j = 0; j < m; j++) {
            file >> duck;
          }
          continue;
        }
      }
      grid_point.set(0, level[0], index[0]);
      grid_point.set(1, level[1], index[1]);
      level_ = std::max(level_, std::max(level[0], level[1]));
    }
    // add grid point
    grid_->getStorage().insert(grid_point);
    if (dim == 3) {
      // shearing angle should be optimized ==> Read ALL The Data!
      // (except for the final value, the volume)
      file >> alpha1_[j] >> alpha2_[j] >> alpha3_[j] >> alpha4_[j]
           >> alpha5_[j] >> alpha6_[j];
      if (m == 7) {
        file >> alpha7_[j];
      }
      if (notation == HILL_MANDEL) {
        alpha6_[j] /= 2.0;
      }
    } else {
      // shearing angle should not be optimized ==> maybe we have to pick the right data
      switch (m) {
      case 4:
        file >> alpha1_[j] >> alpha2_[j] >> alpha3_[j] >> alpha4_[j];
        break;
      case 5:
        file >> alpha1_[j] >> alpha2_[j] >> alpha3_[j] >> alpha4_[j] >> alpha7_[j];
        break;
      case 6:
        file >> alpha1_[j] >> alpha2_[j] >> duck >> alpha3_[j] >> duck >> alpha4_[j];
        break;
      case 7:
        file >> alpha1_[j] >> alpha2_[j] >> duck >> alpha3_[j] >> duck >> alpha4_[j] >> alpha7_[j];
        break;
      }
      if (notation == HILL_MANDEL) {
        alpha4_[j] /= 2.0;
      }
    }
    if (dim == 3) {
      LOG_DBG3(dm) << grid_point.getCoord(0) << " " << grid_point.getCoord(1) << " " << grid_point.getCoord(2) << " -> "
          << alpha1_[j] << " " << alpha2_[j] << " " << alpha3_[j] << " " << alpha4_[j] << " "
          << alpha5_[j] << " " << alpha6_[j] << " " << alpha7_[j];
    } else {
      LOG_DBG3(dm) << grid_point.getCoord(0) << " " << grid_point.getCoord(1) << " -> "
          << alpha1_[j] << " " << alpha2_[j] << " " << alpha3_[j] << " " << alpha4_[j] << " " << alpha7_[j];
    }
    j++;
  }
  LOG_DBG(dm) << "DM::ISG: level = " << level_ << "\n";
  file.close();
  if (dim == 2) {
    // coefficient vectors were too big, because we skipped grid points
    alpha1_.resize(grid_->getStorage().getSize());
    alpha2_.resize(grid_->getStorage().getSize());
    alpha3_.resize(grid_->getStorage().getSize());
    alpha4_.resize(grid_->getStorage().getSize());
    alpha7_.resize(grid_->getStorage().getSize());
  }

  // DEBUG
  /*std::cout << "m = " << m << ", d = " << d << ", shearIsDesign_ = " << shearIsDesign_ << "\n";
  std::cout << "alpha2_before_hierarchization = [";
  for (size_t i = 0; i < alpha2_.getSize(); i++) {
    if (i > 0) {
      std::cout << ", ";
    }
    std::cout << alpha2_[i];
  }
  std::cout << "]\n";*/

  // hierarchization if needed
  if (!hierarchized) {
    HierarchizeSparseGridCoefficients();
  }

  // DEBUG
  /*std::cout << "alpha2_after_hierarchization = [";
  for (size_t i = 0; i < alpha2_.getSize(); i++) {
    if (i > 0) {
      std::cout << ", ";
    }
    std::cout << alpha2_[i];
  }
  std::cout << "]\n";*/

//  EvaluateFullGrid();
}

void DesignMaterial::EvaluateFullGrid() {
  Matrix<double> E;
  double stepsize;
  LocPoint p;
  p.coord.Resize(3);

  p.coord[2] = .5;
  stepsize = 1./128;

  std::ofstream file;
  file.open("gridData.5.txt");
  file.precision(8);
  file.flags(std::ios::scientific);

  for (unsigned int ii = 1; ii < 1/stepsize; ii++) {
    p.coord[0] = ii * stepsize;
    for (unsigned int jj = 1; jj < 1/stepsize; jj++) {
      p.coord[1] = jj * stepsize;
      GetHomRectSGPPTensor(E,p.coord,DesignElement::NO_DERIVATIVE,FULL);
      file << p.coord[0] << "\t" << p.coord[1] << "\t" << p.coord[2] << "\t" << E(0,0) << "\t" << E(0,1) << "\t" << E(0,2) << "\t" << E(1,1) << "\t" << E(1,2) << "\t" << E(2,2) << std::endl;
    }
  }
  file.close();

  p.coord[2] = .125;

  file.open("gridData.125.txt");
  file.precision(8);
  file.flags(std::ios::scientific);

  for (unsigned int ii = 1; ii < 1/stepsize; ii++) {
    p.coord[0] = ii * stepsize;
    for (unsigned int jj = 1; jj < 1/stepsize; jj++) {
      p.coord[1] = jj * stepsize;
      GetHomRectSGPPTensor(E,p.coord,DesignElement::NO_DERIVATIVE,FULL);
      file << p.coord[0] << "\t" << p.coord[1] << "\t" << p.coord[2] << "\t" << E(0,0) << "\t" << E(0,1) << "\t" << E(0,2) << "\t" << E(1,1) << "\t" << E(1,2) << "\t" << E(2,2) << std::endl;
    }
  }
  file.close();
}

void DesignMaterial::FillSparseGridWithFullGridData(Matrix<double>& data) {
  sgpp::base::GridStorage& gridStorage = grid_->getStorage();

  // create coefficient vectors
  alpha1_.resize(gridStorage.getSize());
  alpha2_.resize(gridStorage.getSize());
  alpha3_.resize(gridStorage.getSize());
  alpha4_.resize(gridStorage.getSize());
  alpha5_.resize(gridStorage.getSize());
  alpha6_.resize(gridStorage.getSize());
  alpha7_.resize(gridStorage.getSize());

  // put data values in coefficient vectors
  unsigned int dim1, dim2, dim3, index1, index2, index3, row;
  dim1 = catalogueSize_[0];
  dim2 = catalogueSize_[1];
  dim3 = catalogueSize_[2];
  sgpp::base::GridIndex* gp;
  unsigned int sz = gridStorage.getSize();
  sz = sz + 1;
  for (unsigned int i=0; i < gridStorage.getSize(); i++) {
    gp = gridStorage.get(i);
    if (catalogueSize_.GetSize() == 2) {
      index1 = gp->getCoord(0)*(dim1);
      index2 = gp->getCoord(1)*(dim2);
      row = (index1-1)*dim2 + index2 - 1;
    } else {
      index1 = gp->getCoord(0)*(dim1);
      index2 = gp->getCoord(1)*(dim2);
      if(shearIsDesign_) {
        index3 = gp->getCoord(2)*(dim3+1);
      } else {
        index3 = .5*(dim3+1);
      }
      row = (index1-1)*dim2*dim3 + (index2-1)*dim3 + index3 - 1;
    }
    alpha1_[i] = data[row][0];
    alpha2_[i] = data[row][1];
    if (grid_->getStorage().getDimension() == 3) {
      alpha3_[i] = data[row][2];
      alpha4_[i] = data[row][3];
      alpha5_[i] = data[row][4];
      alpha6_[i] = data[row][5];
      if (data.GetNumCols() == 7) {
        alpha7_[i] = data[row][6];
      }
    } else {
      if (catalogueSize_.GetSize() == 2) {
        alpha3_[i] = data[row][2];
        alpha4_[i] = data[row][3];
        if (data.GetNumCols() == 5) {
          alpha7_[i] = data[row][4];
        }
      } else {
        alpha3_[i] = data[row][3];
        alpha4_[i] = data[row][5];
        if (data.GetNumCols() == 7) {
          alpha7_[i] = data[row][6];
        }
      }
    }
    LOG_DBG3(dm) << "DM:FSGF: " << gp->getCoord(0) << " " << gp->getCoord(1) << " " << gp->getCoord(2) << " -> "
        << alpha1_[i] << " " << alpha2_[i] << " " << alpha3_[i] << " " << alpha4_[i] << " " << alpha7_[i];
  }
  // hierarchize data vectors
  HierarchizeSparseGridCoefficients();
}

void DesignMaterial::FillSparseGridWithSparseGridData(Matrix<double>& data) {
  sgpp::base::GridStorage& gridStorage = grid_->getStorage();
  // Catalogue has to be in correct order! High error risk!
  // Better use new format instead

  // create coefficient vectors
  alpha1_.resize(gridStorage.getSize());
  alpha2_.resize(gridStorage.getSize());
  alpha3_.resize(gridStorage.getSize());
  alpha4_.resize(gridStorage.getSize());
  if (gridStorage.getDimension() == 3) {
    alpha5_.resize(gridStorage.getSize());
    alpha6_.resize(gridStorage.getSize());
  }
  alpha7_.resize(gridStorage.getSize());

  // put data values in coefficient vectors
  sgpp::base::GridIndex* gp;
  unsigned int sz = gridStorage.getSize();
  sz = sz + 1;
  for (unsigned int i=0; i < gridStorage.getSize(); i++) {
    if (gridStorage.getDimension() == 3) {
      alpha1_[i] = data[i][0];
      alpha2_[i] = data[i][1];
      alpha3_[i] = data[i][2];
      alpha4_[i] = data[i][3];
      alpha5_[i] = data[i][4];
      alpha6_[i] = data[i][5];
      if (data.GetNumCols() == 7) {
        alpha7_[i] = data[i][6];
      }
    } else {
      alpha1_[i] = data[i][0];
      alpha2_[i] = data[i][1];
      if (catalogueSize_.GetSize() == 2) {
        alpha3_[i] = data[i][2];
        alpha4_[i] = data[i][3];
        if (data.GetNumCols() == 5) {
          alpha7_[i] = data[i][4];
        }
      } else {
        alpha3_[i] = data[i][3];
        alpha4_[i] = data[i][5];
        if (data.GetNumCols() == 7) {
          alpha7_[i] = data[i][6];
        }
      }
    }
    gp = gridStorage.get(i);
    if (gridStorage.getDimension() == 3) {
      LOG_DBG3(dm) << "DM:FSGS: " << gp->getCoord(0) << " " << gp->getCoord(1) << " " << gp->getCoord(2) << " -> "
          << alpha1_[i] << " " << alpha2_[i] << " " << alpha3_[i] << " " << alpha4_[i] << " " << alpha7_[i];
    } else {
      LOG_DBG3(dm) << "DM:FSGS: " << gp->getCoord(0) << " " << gp->getCoord(1) << " -> "
          << alpha1_[i] << " " << alpha2_[i] << " " << alpha3_[i] << " " << alpha4_[i] << " " << alpha7_[i];
    }
  }
  // hierarchize data vectors
  HierarchizeSparseGridCoefficients();
}

void DesignMaterial::HierarchizeSparseGridCoefficients() {
  std::cout << "++ Hierarchizing data vectors of Sparse Grid\n" << std::flush;
  if ((sgpp_basis_ == LINEAR) || (sgpp_basis_ == MODLINEAR)) {
    std::unique_ptr<sgpp::base::OperationHierarchisation> hierOp =
        sgpp::op_factory::createOperationHierarchisation(*grid_);
    hierOp->doHierarchisation(alpha1_);
    hierOp->doHierarchisation(alpha2_);
    hierOp->doHierarchisation(alpha3_);
    hierOp->doHierarchisation(alpha4_);
    if (grid_->getStorage().getDimension() == 3) {
      hierOp->doHierarchisation(alpha5_);
      hierOp->doHierarchisation(alpha6_);
    }
    hierOp->doHierarchisation(alpha7_);
  } else {
    sgpp::base::DataMatrix alphas(alpha1_.getSize(), (grid_->getStorage().getDimension() == 3 ? 7 : 5));
    alphas.setColumn(0, alpha1_);
    alphas.setColumn(1, alpha2_);
    alphas.setColumn(2, alpha3_);
    alphas.setColumn(3, alpha4_);
    if (grid_->getStorage().getDimension() == 3) {
      alphas.setColumn(4, alpha5_);
      alphas.setColumn(5, alpha6_);
      alphas.setColumn(6, alpha7_);
    } else {
      alphas.setColumn(4, alpha7_);
    }
    
    std::unique_ptr<sgpp::optimization::OperationMultipleHierarchisation> hierOp =
      sgpp::op_factory::createOperationMultipleHierarchisation(*grid_);
    hierOp->doHierarchisation(alphas);
    
    alphas.getColumn(0, alpha1_);
    alphas.getColumn(1, alpha2_);
    alphas.getColumn(2, alpha3_);
    alphas.getColumn(3, alpha4_);
    if (grid_->getStorage().getDimension() == 3) {
      alphas.getColumn(4, alpha5_);
      alphas.getColumn(5, alpha6_);
      alphas.getColumn(6, alpha7_);
    } else {
      alphas.getColumn(4, alpha7_);
    }
  }
}

void DesignMaterial::GetHomRectSGPPTensor(MaterialTensor<double>& mt, Vector<double>& p,
     DesignElement::Type direction, SubTensorType subTensor) {
  // Method uses SGPP interpolation
  sgpp::base::DataVector point(p.GetPointer(), p.GetSize());
  LOG_DBG2(dm) << p;

  Matrix<double>& E = mt.GetMatrix(HILL_MANDEL);
  E.Resize(3,3);
  E.Init(); // for off-diagonal

  if ((sgpp_basis_ == LINEAR) || (sgpp_basis_ == MODLINEAR)) {
    if (direction == DesignElement::NO_DERIVATIVE || direction == DesignElement::ROTANGLE || direction == DesignElement::ROTANGLETHIRD || direction == DesignElement::ROTANGLESECOND) {
      if (!shearIsDesign_) { // no shearing
        E[1-1][1-1] = op_eval_->eval(alpha1_, point);
        E[1-1][2-1] = op_eval_->eval(alpha2_, point);
        E[2-1][1-1] = E[1-1][2-1];
        E[2-1][2-1] = op_eval_->eval(alpha3_, point);
        E[3-1][3-1] = op_eval_->eval(alpha4_, point);
      } else { // shearing
        E[1-1][1-1] = op_eval_->eval(alpha1_, point);
        E[1-1][2-1] = op_eval_->eval(alpha2_, point);
        E[1-1][3-1] = op_eval_->eval(alpha3_, point);
        E[2-1][1-1] = E[1-1][2-1];
        E[2-1][2-1] = op_eval_->eval(alpha4_, point);
        E[2-1][3-1] = op_eval_->eval(alpha5_, point);
        E[3-1][1-1] = E[1-1][3-1];
        E[3-1][2-1] = E[2-1][3-1];
        E[3-1][3-1] = op_eval_->eval(alpha6_, point);
      }
      LOG_DBG(dm)<<" E11= "<<E[0][0]<<" E12= "<<E[0][1]<<" E13= "<<E[0][2]<<" E22= "<< E[1][1]<<" E23= "<<E[1][2]<<" E33= "<<E[2][2];
    } else {
      if (!shearIsDesign_) { // no shearing
        E[1-1][1-1] = EvaluateSGPPInterpolation_Deriv(alpha1_, point, direction);
        E[1-1][2-1] = EvaluateSGPPInterpolation_Deriv(alpha2_, point, direction);
        E[2-1][1-1] = E[1-1][2-1];
        E[2-1][2-1] = EvaluateSGPPInterpolation_Deriv(alpha3_, point, direction);
        E[3-1][3-1] = EvaluateSGPPInterpolation_Deriv(alpha4_, point, direction);
      } else { // shearing
        E[1-1][1-1] = EvaluateSGPPInterpolation_Deriv(alpha1_, point, direction);
        E[1-1][2-1] = EvaluateSGPPInterpolation_Deriv(alpha2_, point, direction);
        E[1-1][3-1] = EvaluateSGPPInterpolation_Deriv(alpha3_, point, direction);
        E[2-1][1-1] = E[1-1][2-1];
        E[2-1][2-1] = EvaluateSGPPInterpolation_Deriv(alpha4_, point, direction);
        E[2-1][3-1] = EvaluateSGPPInterpolation_Deriv(alpha5_, point, direction);
        E[3-1][1-1] = E[1-1][3-1];
        E[3-1][2-1] = E[2-1][3-1];
        E[3-1][3-1] = EvaluateSGPPInterpolation_Deriv(alpha6_, point, direction);
      }
      LOG_DBG(dm)<<"Derivative "<<((direction == DesignElement::STIFF1)?"1":(direction == DesignElement::STIFF2)?"2":"3")
          <<" E11= "<<E[0][0]<<" E12= "<<E[0][1]<<" E13= "<<E[0][2]<<" E22= "<< E[1][1]<<" E23= "<<E[1][2]<<" E33= "<<E[2][2];
    }
  } else {
    if (direction == DesignElement::NO_DERIVATIVE || direction == DesignElement::ROTANGLE || direction == DesignElement::ROTANGLETHIRD || direction == DesignElement::ROTANGLESECOND) {
      if (!shearIsDesign_) { // no shearing
        E[1-1][1-1] = op_naive_eval_->eval(alpha1_, point);
        E[1-1][2-1] = op_naive_eval_->eval(alpha2_, point);
        E[2-1][1-1] = E[1-1][2-1];
        E[2-1][2-1] = op_naive_eval_->eval(alpha3_, point);
        E[3-1][3-1] = op_naive_eval_->eval(alpha4_, point);
      } else { // shearing
        E[1-1][1-1] = op_naive_eval_->eval(alpha1_, point);
        E[1-1][2-1] = op_naive_eval_->eval(alpha2_, point);
        E[1-1][3-1] = op_naive_eval_->eval(alpha3_, point);
        E[2-1][1-1] = E[1-1][2-1];
        E[2-1][2-1] = op_naive_eval_->eval(alpha4_, point);
        E[2-1][3-1] = op_naive_eval_->eval(alpha5_, point);
        E[3-1][1-1] = E[1-1][3-1];
        E[3-1][2-1] = E[2-1][3-1];
        E[3-1][3-1] = op_naive_eval_->eval(alpha6_, point);
      }
      LOG_DBG(dm)<<" E11= "<<E[0][0]<<" E12= "<<E[0][1]<<" E13= "<<E[0][2]<<" E22= "<< E[1][1]<<" E23= "<<E[1][2]<<" E33= "<<E[2][2];
    } else {
      if (!shearIsDesign_) { // no shearing
        E[1-1][1-1] = EvaluateSGPPInterpolation_Deriv(alpha1_, point, direction);
        E[1-1][2-1] = EvaluateSGPPInterpolation_Deriv(alpha2_, point, direction);
        E[2-1][1-1] = E[1-1][2-1];
        E[2-1][2-1] = EvaluateSGPPInterpolation_Deriv(alpha3_, point, direction);
        E[3-1][3-1] = EvaluateSGPPInterpolation_Deriv(alpha4_, point, direction);
      } else { // shearing
        E[1-1][1-1] = EvaluateSGPPInterpolation_Deriv(alpha1_, point, direction);
        E[1-1][2-1] = EvaluateSGPPInterpolation_Deriv(alpha2_, point, direction);
        E[1-1][3-1] = EvaluateSGPPInterpolation_Deriv(alpha3_, point, direction);
        E[2-1][1-1] = E[1-1][2-1];
        E[2-1][2-1] = EvaluateSGPPInterpolation_Deriv(alpha4_, point, direction);
        E[2-1][3-1] = EvaluateSGPPInterpolation_Deriv(alpha5_, point, direction);
        E[3-1][1-1] = E[1-1][3-1];
        E[3-1][2-1] = E[2-1][3-1];
        E[3-1][3-1] = EvaluateSGPPInterpolation_Deriv(alpha6_, point, direction);
      }
      LOG_DBG(dm)<<"Derivative "<<((direction == DesignElement::STIFF1)?"1":(direction == DesignElement::STIFF2)?"2":"3")
          <<" E11= "<<E[0][0]<<" E12= "<<E[0][1]<<" E13= "<<E[0][2]<<" E22= "<< E[1][1]<<" E23= "<<E[1][2]<<" E33= "<<E[2][2];
    }
  }
}

void DesignMaterial::GetHomRectFullBsplineTensor(MaterialTensor<double>& mt, Vector<double>& p,
     DesignElement::Type direction, SubTensorType subTensor) const {

  Matrix<double>& E = mt.GetMatrix(HILL_MANDEL);
  E.Resize(3,3);
  E.Init(); // for off-diagonal
  const int margin = (bspline_degree_+1)/2 + 1;
  const double a = static_cast<double>((bspline_degree_+1)/2);
  sgpp::base::BsplineBasis<int, int> bspline(bspline_degree_);

  if (!shearIsDesign_) {
    E[1-1][1-1] = 0.0;
    E[1-1][2-1] = 0.0;
    E[2-1][1-1] = 0.0;
    E[2-1][2-1] = 0.0;
    E[3-1][3-1] = 0.0;
    // loop through all relevant B-splines (the ones that don't vanish in p)
    for (int i1 = static_cast<int>(64.0 * p[0] - a + 1.0);
             i1 <= static_cast<int>(64.0 * p[0] + a); i1++) {
      // B-spline out of bounds, skip
      if ((i1 < 1-margin) || (i1 > 64+margin)) continue;
      for (int i2 = static_cast<int>(64.0 * p[1] - a + 1.0);
               i2 <= static_cast<int>(64.0 * p[1] + a); i2++) {
        // B-spline out of bounds, skip
        if ((i2 < 1-margin) || (i2 > 64+margin)) continue;
        double bspl_val;
        // evaluate tensor product B-spline
        if (direction == DesignElement::NO_DERIVATIVE || direction == DesignElement::ROTANGLE || direction == DesignElement::ROTANGLETHIRD || direction == DesignElement::ROTANGLESECOND) {
          // 6 because 64 is 2 to the 6th
          bspl_val = bspline.eval(6, i1, p[0]) * bspline.eval(6, i2, p[1]);
        } else if (direction == DesignElement::STIFF1) {
          bspl_val = bspline.evalDx(6, i1, p[0]) * bspline.eval(6, i2, p[1]);
        } else if (direction == DesignElement::STIFF2) {
          bspl_val = bspline.eval(6, i1, p[0]) * bspline.evalDx(6, i2, p[1]);
        } else {
          std::cout << "ERROR\n";
          bspl_val = 0.0;
        }
        // calculate coefficient of B-spline
        const int i = (64+2*margin) * (i1+margin-1) + (i2+margin-1);
        // update result
        E[1-1][1-1] += full_bspline_coeff11_(i, 0) * bspl_val;
        E[1-1][2-1] += full_bspline_coeff12_(i, 0) * bspl_val;
        E[2-1][2-1] += full_bspline_coeff22_(i, 0) * bspl_val;
        E[3-1][3-1] += full_bspline_coeff33_(i, 0) * bspl_val;
      }
    }
//    // take positive part (if we're not calculating derivatives)
//    if (direction == DesignElement::NO_DERIVATIVE || direction == DesignElement::ROTANGLE || direction == DesignElement::ROTANGLETHIRD || direction == DesignElement::ROTANGLESECOND) {
//      E[1-1][1-1] = std::max(0.0, E[1-1][1-1]);
//      E[1-1][2-1] = std::max(0.0, E[1-1][2-1]);
//      E[2-1][2-1] = std::max(0.0, E[2-1][2-1]);
//      E[3-1][3-1] = std::max(0.0, E[3-1][3-1]);
//    }
    // symmetric entry
    E[2-1][1-1] = E[1-1][2-1];
  } else {
    E[1-1][1-1] = 0.0;
    E[1-1][2-1] = 0.0;
    E[1-1][3-1] = 0.0;
    E[2-1][1-1] = 0.0;
    E[2-1][2-1] = 0.0;
    E[2-1][3-1] = 0.0;
    E[3-1][1-1] = 0.0;
    E[3-1][2-1] = 0.0;
    E[3-1][3-1] = 0.0;
    // loop through all relevant B-splines (the ones that don't vanish in p)
    for (int i1 = static_cast<int>(64.0 * p[0] - a + 1.0);
             i1 <= static_cast<int>(64.0 * p[0] + a); i1++) {
      // B-spline out of bounds, skip
      if ((i1 < 1-margin) || (i1 > 64+margin)) continue;
      for (int i2 = static_cast<int>(64.0 * p[1] - a + 1.0);
               i2 <= static_cast<int>(64.0 * p[1] + a); i2++) {
        // B-spline out of bounds, skip
        if ((i2 < 1-margin) || (i2 > 64+margin)) continue;
        for (int i3 = static_cast<int>(64.0 * p[2] - a + 1.0);
                 i3 <= static_cast<int>(64.0 * p[2] + a); i3++) {
          // B-spline out of bounds, skip
          if ((i3 < 1-margin) || (i3 > 63+margin)) continue;
          double bspl_val;
          // evaluate tensor product B-spline
          if (direction == DesignElement::NO_DERIVATIVE || direction == DesignElement::ROTANGLE || direction == DesignElement::ROTANGLETHIRD || direction == DesignElement::ROTANGLESECOND) {
            // 6 because 64 is 2 to the 6th
            bspl_val = bspline.eval(6, i1, p[0]) * bspline.eval(6, i2, p[1]) * bspline.eval(6, i3, p[2]);
          } else if (direction == DesignElement::STIFF1) {
            bspl_val = bspline.evalDx(6, i1, p[0]) * bspline.eval(6, i2, p[1]) * bspline.eval(6, i3, p[2]);
          } else if (direction == DesignElement::STIFF2) {
            bspl_val = bspline.eval(6, i1, p[0]) * bspline.evalDx(6, i2, p[1]) * bspline.eval(6, i3, p[2]);
          } else if ((direction == DesignElement::STIFF3) || (direction == DesignElement::SHEAR1)) {
            bspl_val = bspline.eval(6, i1, p[0]) * bspline.eval(6, i2, p[1]) * bspline.evalDx(6, i3, p[2]);
          } else {
            std::cout << "ERROR\n";
            bspl_val = 0.0;
          }
          // calculate coefficient of B-spline
          const int i = (63+2*margin)*(64+2*margin) * (i1+margin-1) + (63+2*margin) * (i2+margin-1) + (i3+margin-1);
          // update result
          E[1-1][1-1] += full_bspline_coeff11_(i, 0) * bspl_val;
          E[1-1][2-1] += full_bspline_coeff12_(i, 0) * bspl_val;
          E[1-1][3-1] += full_bspline_coeff13_(i, 0) * bspl_val;
          E[2-1][2-1] += full_bspline_coeff22_(i, 0) * bspl_val;
          E[2-1][3-1] += full_bspline_coeff23_(i, 0) * bspl_val;
          E[3-1][3-1] += full_bspline_coeff33_(i, 0) * bspl_val;
        }
      }
    }
//    // take positive part (if we're not calculating derivatives)
//    if (direction == DesignElement::NO_DERIVATIVE || direction == DesignElement::ROTANGLE || direction == DesignElement::ROTANGLETHIRD || direction == DesignElement::ROTANGLESECOND) {
//      E[1-1][1-1] = std::max(0.0, E[1-1][1-1]);
//      E[1-1][2-1] = std::max(0.0, E[1-1][2-1]);
//      E[2-1][2-1] = std::max(0.0, E[2-1][2-1]);
//      E[3-1][3-1] = std::max(0.0, E[3-1][3-1]);
//    }
    // symmetric entries
    E[2-1][1-1] = E[1-1][2-1];
    E[3-1][1-1] = E[1-1][3-1];
    E[3-1][2-1] = E[2-1][3-1];
  }
}

double DesignMaterial::EvaluateSGPPInterpolation_Deriv(sgpp::base::DataVector& alpha, sgpp::base::DataVector& point, DesignElement::Type direction) const {
  // Approximates the derivative with finite differences
  unsigned int dimension;
  double h = 1./pow(2,level_+1) * 1e-6;
  switch (direction) {
    case DesignElement::STIFF1:
      dimension = 1;
      break;
    case DesignElement::STIFF2:
      dimension = 2;
      break;
    case DesignElement::STIFF3:
    case DesignElement::SHEAR1:
      dimension = 3;
      break;
    default:
      return 0;
    }
  sgpp::base::DataVector pointL = point;
  sgpp::base::DataVector pointU = point;
  pointL[dimension-1] -= h;
  pointU[dimension-1] += h;

  double valL = op_eval_->eval(alpha,pointL);
  double valU = op_eval_->eval(alpha,pointU);

  return (valU - valL) / (2*h);
}

#endif //USE_SGPP


double DesignMaterial::CalcHomVolume(Vector<double>& p, DesignElement::Type direction, bool derivative) {

#ifdef USE_SGPP
  // Method uses SGPP interpolation
  sgpp::base::DataVector point(p.GetPointer(), p.GetSize());
  LOG_DBG2(dm) << p;

  double vol;

  if ((sgpp_basis_ == LINEAR) || (sgpp_basis_ == MODLINEAR)) {
    std::unique_ptr<sgpp::base::OperationEval> opEval = sgpp::op_factory::createOperationEval(*grid_);
    if (!derivative) {
      vol = opEval->eval(alpha7_, point);
      LOG_DBG(dm) << "DM::CHV: volume= " << vol;
    } else {
      vol = EvaluateSGPPInterpolation_Deriv(alpha7_, point, direction);
    }
  } else {
    if (!derivative) {
      std::unique_ptr<sgpp::base::OperationNaiveEval> op_naive_eval_ = sgpp::op_factory::createOperationNaiveEval(*grid_);
      vol = op_naive_eval_->eval(alpha7_, point);
      LOG_DBG(dm) << "DM::CHV: volume= " << vol;
    } else {
      vol = EvaluateSGPPInterpolation_Deriv(alpha7_, point, direction);
    }
  }
  return vol;

#else //USE_SGPP
  // should never be reached: CalcHomVolume is only called in Function.cc when interpolation method
  // is set to SGPP which can only occur when CFS is compiled with SGPP
  return -1;

#endif //USE_SGPP
}

bool DesignMaterial::GetErsatzElementMatrixMSFEM(Matrix<double>& A,
    const Elem* elem,DesignElement::Type direction) {
    assert((type_ == MSFEM_C1));

    // collect all parameters
    if(!CollectMaterialParametersForElement(space_, elem))
      throw Exception("no elem data for MSFEM defined");

    // read design variables
    double a = GetParameter(DesignElement::STIFF1);
    double b = GetParameter(DesignElement::STIFF2);
    //double rotAngle = 0.;
    //if (HasParameter(DesignElement::ROTANGLE)) {
    //  rotAngle = params_[DesignElement::ROTANGLE];
    //}
    Vector<double> p(2);
    //if (HasParameter(DesignElement::ROTANGLE)) {
    //  p.Resize(3);
    //}
    p[0] = a;
    p[1] = b;
    //if (HasParameter(DesignElement::ROTANGLE)) {
    //  p[2] = rotAngle;
    //}
    // length of the discretized design interval
    int m = msfem_a_.GetSize();
    int n = msfem_b_.GetSize();
    //int o = -1;
    //if (HasParameter(DesignElement::ROTANGLE)) {
    //  o = msfem_c_.GetNumRows();
    //}

    int j = GetInterpolationIndex(msfem_a_, p[0]);
    int k = GetInterpolationIndex(msfem_b_, p[1]);
    //int l = -1;
    //if (HasParameter(DesignElement::ROTANGLE)) {
    //  l = GetInterpolationIndex(msfem_c_, p[2]);
    //}

    int count = 0;
    A.Resize(8,8);
    for (int ii = 0;ii<8;ii++) {
      for (int jj = ii;jj<8;jj++) {
        //if (HasParameter(DesignElement::ROTANGLE)) {
        //  if (direction == DesignElement::NO_DERIVATIVE) {
        //    A[ii][jj] = EvaluateC1Interpolation3D(p, msfem_coeff_[count], j, k, l, m, n, o);
        //  } else {
        //    A[ii][jj] = EvaluateC1Interpolation3D(p, msfem_coeff_[count], j, k, l, m, n, o);
        //  }
        //} else {
          if (direction == BaseDesignElement::NO_DERIVATIVE) {
            A[ii][jj] = EvaluateC1Interpolation2D(p, msfem_coeff_[count], j, k, m, n);
          } else {
            A[ii][jj] = EvaluateC1Interpolation2D_Deriv(p, msfem_coeff_[count], j, k, m, n, direction);
          }
        //}
        if (ii!=jj) {
          A[jj][ii] = A[ii][jj];
        }
        count++;
      }
    }
    return true;
}

int DesignMaterial::GetInterpolationIndex(const Vector<double>& interval, double& point) const {
  int sz = interval.GetSize();
  assert(sz > 0);
  int idx = -1;
  // set index for values close to boundaries manually
  if (close(point, interval[0])) {
    idx = 0;
  } else if (point > interval[sz - 1]) {
    idx = sz - 2;
    point = interval[sz - 1]; //FIXME
  } else if (point < interval[0]) {
    idx = 0;
    point = interval[0]; //FIXME
  } else { // interval[0] < point && point <= interval[sz-1]
    idx = std::upper_bound(interval.GetPointer(), interval.GetPointer() + interval.GetSize() - 1, point) - interval.GetPointer() - 1;
  }
  assert(idx > -1);
  assert(idx < sz-1);
  return idx;
}

double DesignMaterial::EvaluateC1Interpolation_Deriv(Vector<double>& p, ApproxData* interpolator, DesignElement::Type direction) const {
  if (interpolator == NULL)
    return 0.0;

  switch (direction) {
  case DesignElement::STIFF1:
    return interpolator->EvaluateDeriv(p, 0);
  case DesignElement::STIFF2:
    return interpolator->EvaluateDeriv(p, 1);
  case DesignElement::STIFF3:
  case DesignElement::SHEAR1:
    return interpolator->EvaluateDeriv(p, 2);
  default:
    return 0.0;
  }
}

double DesignMaterial::EvaluateC1Interpolation1D_Deriv(double p, const Matrix<double> & coeff,
    int & j, DesignElement::Type direction) const {

  assert(direction != DesignElement::NO_TYPE);

  double t;
  if (type_ == MSFEM_C1) {
    t = (p - msfem_a_[j]) / (msfem_a_[j+1] - msfem_a_[j]);
  } else {
    t = (p - hom_rect_a_[j]) / (hom_rect_a_[j+1] - hom_rect_a_[j]);
  }
  LOG_DBG3(dm) << "EC1ID: t = " << t << "\n";

  double tmp_t = t;
  double pow_t[4];
  // precalculate powers
  pow_t[1] = 1;
  for (int i = 2; i < 4; i++) {
    pow_t[i] = tmp_t;
    tmp_t *= t;
  }

  double deriv = 0;
  for (int i = 1; i < 4; i++) {
    deriv += coeff[j][i] * i * pow_t[i];
  }
  deriv /= type_ == MSFEM_C1 ? msfem_a_[j+1] - msfem_a_[j] : hom_rect_a_[j+1] - hom_rect_a_[j];
  LOG_DBG3(dm) << "EC1ID: Result =" << deriv;
  return deriv;
}

double DesignMaterial::EvaluateC1Interpolation1D_Deriv2(double p, const Matrix<double> & coeff,
    int & j, DesignElement::Type direction) const {

  assert(direction != DesignElement::NO_TYPE);

  double t;
  if (type_ == MSFEM_C1) {
    t = (p - msfem_a_[j]) / (msfem_a_[j+1] - msfem_a_[j]);
  } else {
    t = (p - hom_rect_a_[j]) / (hom_rect_a_[j+1] - hom_rect_a_[j]);
  }
  LOG_DBG3(dm) << "EC1ID: t = " << t << "\n";

  double pow_t[4];
  pow_t[2] = 1;
  pow_t[3] = t;

  double deriv = 0;
  for (int i = 2; i < 4; i++) {
    deriv += coeff[j][i] * i * (i-1) * pow_t[i];
  }
  deriv /= type_ == MSFEM_C1 ? msfem_a_[j+1] - msfem_a_[j] : hom_rect_a_[j+1] - hom_rect_a_[j];
  LOG_DBG3(dm) << "EC1ID: Result =" << deriv;
  return deriv;
}

double DesignMaterial::EvaluateC1Interpolation2D(Vector<double>& p,
    const Matrix<double> & coeff, int & j, int & k, int & m, int & n) const {
  LOG_DBG3(dm) << "EC1I: p=[" << p[0] << "," << p[1] << "]";
  double u,t;
  if (type_ == MSFEM_C1) {
    t = (p[0] - msfem_a_[j]) / (msfem_a_[j+1] - msfem_a_[j]);
    u = (p[1] - msfem_b_[k]) / (msfem_b_[k+1] - msfem_b_[k]);
  } else {
    t = (p[0] - hom_rect_a_[j]) / (hom_rect_a_[j+1] - hom_rect_a_[j]);
    u = (p[1] - hom_rect_b_[k]) / (hom_rect_b_[k+1] - hom_rect_b_[k]);
  }
  LOG_DBG3(dm) << "EC1I: u = " << u << " t = " << t << "\n";
  LOG_DBG3(dm) << "EC1I: j = " << j << " k = " << k << "\n";

  double tmp_t = t, tmp_u = u;
  double pow_u[4], pow_t[4];
  pow_u[0] = 1;
  pow_t[0] = 1;
  // precalculate powers
  for (int i = 1; i < 4; i++) {
    pow_t[i] = tmp_t;
    pow_u[i] = tmp_u;
    tmp_t *= t;
    tmp_u *= u;
  }

  double res = 0;
  for (int i = 0;i<4;i++) {
    for (int l=0;l<4;l++) {
      res += coeff[(n-1)*j+k][i*4+l]*pow_t[i]*pow_u[l];
    }
  }
  LOG_DBG3(dm) << "EC1I: Result =" << res;
  return res;
}

double DesignMaterial::EvaluateC1Interpolation2D_Deriv(Vector<double>& p,
    const Matrix<double> & coeff, int & j, int & k,
    int & m, int & n, DesignElement::Type direction) const {
  double u,t;
  if (type_ == MSFEM_C1) {
    t = (p[0] - msfem_a_[j]) / (msfem_a_[j+1] - msfem_a_[j]);
    u = (p[1] - msfem_b_[k]) / (msfem_b_[k+1] - msfem_b_[k]);
  } else {
    t = (p[0] - hom_rect_a_[j]) / (hom_rect_a_[j+1] - hom_rect_a_[j]);
    u = (p[1] - hom_rect_b_[k]) / (hom_rect_b_[k+1] - hom_rect_b_[k]);
  }
  LOG_DBG3(dm) << "EC1ID: u = " << u << " t= " << t << "\n";

  double deriv = 0;
  if (direction == DesignElement::STIFF1) {
    for (int i = 1; i < 4; i++) {
      for (int l = 0; l < 4; l++) {
        deriv += coeff[(n - 1) * j + k][(i) * 4 + l] * i * pow(t, i - 1) * pow(u, l);
      }
    }
    deriv /= type_ == MSFEM_C1 ? msfem_a_[j+1] - msfem_a_[j] : hom_rect_a_[j+1] - hom_rect_a_[j];
  }
  if (direction == DesignElement::STIFF2) {
    for (int i = 0; i < 4; i++) {
      for (int l = 1; l < 4; l++) {
        deriv += coeff[(n - 1) * j + k][(i) * 4 + l] * l * pow(t, i) * pow(u, l - 1);
      }
    }
    deriv /= type_ == MSFEM_C1 ? msfem_b_[k+1] - msfem_b_[k] : hom_rect_b_[k+1] - hom_rect_b_[k];
  }
  LOG_DBG3(dm) << "EC1ID: Result =" << deriv;
  return deriv;
}

double DesignMaterial::EvaluateC1Interpolation3D(Vector<double>& p,
    const Matrix<double> & coeff, int & j, int & k, int & l,
    int & m, int & n, int &o) const {
  // FIXME
  // dirty fix: does nothing if program works correctly
  j = (j > m - 2) ? m - 2 : j;
  k = (k > n - 2) ? n - 2 : k;
  l = (l > o - 2) ? o - 2 : l;
  LOG_DBG(dm)<<"p=["<<p[0]<<","<<p[1]<<", "<<p[2]<<"]";
  double t,u,v;
  if (type_ == MSFEM_C1) {
    t = (p[0] - msfem_a_[j]) / (msfem_a_[j+1] - msfem_a_[j]);
    u = (p[1] - msfem_b_[k]) / (msfem_b_[k+1] - msfem_b_[k]);
    v = (p[2] - msfem_c_[l]) / (msfem_c_[l+1] - msfem_c_[l]);
  } else {
    t = (p[0] - hom_rect_a_[j]) / (hom_rect_a_[j+1] - hom_rect_a_[j]);
    u = (p[1] - hom_rect_b_[k]) / (hom_rect_b_[k+1] - hom_rect_b_[k]);
    v = (p[2] - hom_rect_c_[l]) / (hom_rect_c_[l+1] - hom_rect_c_[l]);
  }
  LOG_DBG(dm)<<"u = "<<u<<" t= "<<t<<" v= "<<v;
  double res = 0;
  for (int ii = 0; ii<4; ii++) {
    for (int jj=0; jj<4; jj++) {
      for (int kk=0; kk<4; kk++) {
        res += coeff[(n-1)*(o-1)*j+(o-1)*k+l][ii+4*jj+16*kk]*pow(t,ii)*pow(u,jj)*pow(v,kk);
      }
    }
  }
  LOG_DBG(dm) << "Result =" << res;
  return res;
}

double DesignMaterial::EvaluateC1Interpolation3D_Deriv(Vector<double>& p,
    const Matrix<double> & coeff, int & j, int & k, int & l,
    int & m, int & n, int & o, DesignElement::Type direction) const {
  double t,u,v;
  if (type_ == MSFEM_C1) {
    t = (p[0] - msfem_a_[j]) / (msfem_a_[j+1] - msfem_a_[j]);
    u = (p[1] - msfem_b_[k]) / (msfem_b_[k+1] - msfem_b_[k]);
    v = (p[2] - msfem_c_[l]) / (msfem_c_[l+1] - msfem_c_[l]);
  } else {
    t = (p[0] - hom_rect_a_[j]) / (hom_rect_a_[j+1] - hom_rect_a_[j]);
    u = (p[1] - hom_rect_b_[k]) / (hom_rect_b_[k+1] - hom_rect_b_[k]);
    v = (p[2] - hom_rect_c_[l]) / (hom_rect_c_[l+1] - hom_rect_c_[l]);
  }
  LOG_DBG(dm)<<"Deriv: u = "<<u<<" t= "<<t<<" v= "<<v<<" j= "<<j<<" k= "<<k<<" l= "<<l;
  LOG_DBG(dm)<<"p_deriv: ["<<p[0]<<", "<<", "<<p[1]<<", "<<p[2];
  double deriv = 0;
  if (direction == DesignElement::STIFF1) {
    for (int ii = 1; ii < 4; ii++) {
      for (int jj = 0; jj < 4; jj++) {
        for (int kk = 0; kk < 4; kk++) {
          deriv += coeff[(n - 1) * (o - 1) * j + (o - 1) * k + l][ii + 4 * jj
              + 16 * kk] * ii * pow(t, ii - 1) * pow(u, jj) * pow(v, kk);
        }
      }
    }
    deriv /= type_ == MSFEM_C1 ? msfem_a_[j+1] - msfem_a_[j] : hom_rect_a_[j+1] - hom_rect_a_[j];
  }
  if (direction == DesignElement::STIFF2) {
    for (int ii = 0; ii < 4; ii++) {
      for (int jj = 1; jj < 4; jj++) {
        for (int kk = 0; kk < 4; kk++) {
          deriv += coeff[(n - 1) * (o - 1) * j + (o - 1) * k + l][ii + 4 * jj
              + 16 * kk] * jj * pow(t, ii) * pow(u, jj - 1) * pow(v, kk);
        }
      }
    }
    deriv /= type_ == MSFEM_C1 ? msfem_b_[k+1] - msfem_b_[k] : hom_rect_b_[k+1] - hom_rect_b_[k];
  }
  if (direction == DesignElement::STIFF3 || direction == DesignElement::SHEAR1) {
    for (int ii = 0; ii < 4; ii++) {
      for (int jj = 0; jj < 4; jj++) {
        for (int kk = 1; kk < 4; kk++) {
          deriv += coeff[(n - 1) * (o - 1) * j + (o - 1) * k + l][ii + 4 * jj
              + 16 * kk] * kk * pow(t, ii) * pow(u, jj) * pow(v, kk - 1);
        }
      }
    }
    deriv /= type_ == MSFEM_C1 ? msfem_c_[l+1] - msfem_c_[l] : hom_rect_c_[l+1] - hom_rect_c_[l];
  }
  LOG_DBG(dm)<< "Deriv Result =" << deriv;
  return deriv;
}

double DesignMaterial::EvaluateC1Interpolation1D(double p, const Matrix<double>& coeff, int& j) const {
  LOG_DBG3(dm) << "EC1I: p=" << p << "";
  double t;
  if (type_ == MSFEM_C1)
    t = (p - msfem_a_[j]) / (msfem_a_[j+1] - msfem_a_[j]);
  else
    t = (p - hom_rect_a_[j]) / (hom_rect_a_[j+1] - hom_rect_a_[j]);
  LOG_DBG3(dm) << "EC1I: t = " << t << "\n";
  LOG_DBG3(dm) << "EC1I: j = " << j << "\n";

  double tmp_t = t;
  double pow_t[4];
  // precalculate powers
  pow_t[0] = 1;
  for (int i = 1; i < 4; i++) {
    pow_t[i] = tmp_t;
    tmp_t *= t;
  }

  double res = 0;
  for (int i = 0; i < 4; i++) {
    res += coeff[j][i] * pow_t[i];
  }
  LOG_DBG3(dm) << "EC1I: Result =" << res;
  return res;
}

bool DesignMaterial::ReadDetailedStats(const char * filename, Matrix<double>& ret) {
  bool isHierarchized;
  unsigned int dim1, dim2, dim3, nRows, nCols;

  // open file
  std::ifstream file;
  file.open(filename);
  if (file.fail()) {
    LOG_DBG(dm)<<"Cannot open file "<<filename;
    std::cout<<"Cannot open file "<<filename<<std::endl;
    std::terminate();
  }

  // read first line
  string line;
  StdVector<string> strVec;
  getline(file,line,'\n');
  SplitStringListWhitespace(line, strVec);

  // read notation
  std::transform(strVec[3].begin(),strVec[3].end(),strVec[3].begin(), ::tolower);
  MaterialTensorNotation notation = !strVec[3].compare(0,strVec[3].size(),"voigt") ? VOIGT : HILL_MANDEL;

  // initialize matrix
  // if second entry starts with "L", we have hierarchized data
  if (strVec[1].compare(0,1,"L") == 0) {
    isHierarchized = true;
    unsigned int catalogueDimension = strVec[0].compare("3D") == 0 ? 3 : 2;
    level_ = boost::lexical_cast<int>(strVec[1].substr(1));
    for (unsigned int i=0; i<catalogueDimension; i++) {
      catalogueSize_.Push_back(pow(2,level_));
    }
    nRows = boost::lexical_cast<int>(strVec[2]);
  } else {
    isHierarchized = false;
    dim1 = boost::lexical_cast<int>(strVec[0]);
    dim2 = boost::lexical_cast<int>(strVec[1]);
    dim3 = boost::lexical_cast<int>(strVec[2]);
    catalogueSize_.Push_back(dim1);
    catalogueSize_.Push_back(dim2);
    if (dim3 == 0) {
      nRows = dim1*dim2;
    } else {
      catalogueSize_.Push_back(dim3);
      nRows = dim1*dim2*dim3;
    }
    level_ = log2(dim1);
  }

#ifdef USE_SGPP
  bool assertiontest = grid_->getStorage().getDimension() == 3;
#else
  bool assertiontest = shearIsDesign_;
#endif
  if (assertiontest) {
    assert(catalogueSize_.GetSize() == 3);
  }

  if (strVec.GetSize() <= catalogueSize_.GetSize() + 4 + 1) {
    // only tensor entries E11 E12 E22 E33 (and volume) are given
    nCols = 4;
  } else {
    // all tensor entries are given
    nCols = 6;
  }
  Matrix<double> data;
  data.Resize(nRows,nCols);
  data.Init();
  LOG_DBG(dm)<<"ReadDetailedStats: nRows = "<<nRows<<", nCols = "<<nCols<<", level = "<<level_;

  for (unsigned int i=0; i<nRows; i++) {
    // read line and write data into matrix
    getline(file,line,'\n');
    strVec.Clear();
    SplitStringListWhitespace(line, strVec);
    for (unsigned int j=0; j<nCols; j++) {
      data[i][j] = boost::lexical_cast<double>(strVec[j+catalogueSize_.GetSize()]);
    }
    if (notation == HILL_MANDEL) {
      if (nCols == 4) {
        data[i][3] /= 2.0;
      } else {
        data[i][2] /= std::sqrt(2.0);
        data[i][4] /= std::sqrt(2.0);
        data[i][5] /= 2.0;
      }
    }
  }

  // If too many tensor entries are given we extract the needed ones
  if (catalogueSize_.GetSize() == 2 && nCols == 6) {
    // extract the columns for entries E11 E12 E22 E33
    std::vector<UInt> rows(boost::counting_iterator<int>( 0 ), boost::counting_iterator<int>( nRows ));
    std::vector<UInt> cols(4);
    cols[0] = 0;
    cols[1] = 1;
    cols[2] = 3;
    cols[3] = 5;
    data.GetSubMatrixByInd(ret,rows,cols);
  } else {
    ret = data;
  }

  LOG_DBG3(dm)<<"Data: \n"<<ret;

  return isHierarchized;
}

inline void DesignMaterial::GetIN718Tensor(MaterialTensor<double>& mt, SubTensorType subTensor, DesignElement::Type direction){
  double a = (direction == DesignElement::INTERPOLATION) ? 1.0 : GetParameter(DesignElement::INTERPOLATION);
  double ma = (direction == DesignElement::INTERPOLATION) ? -1.0 : 1.0-a;

  Matrix<double>& t = mt.GetMatrix(HILL_MANDEL);
  switch (subTensor) {
  case FULL:
    t.Resize(6, 6);
    t.Init();
//    SetOrthotropicTensor(t, subTensor, a*255.68181818+ma*294.03409091, a*99.43181818+ma*80.0, a*99.43181818+ma*114.34659091,
//        a*255.68181818+ma*166.19318182, a*99.43181818+ma*80.0, a*255.68181818+ma*294.03409091, a*78.125+ma*70.0, a*78.125+ma*70.0, a*78.125+ma*60.0);
//    SetOrthotropicTensor(t, subTensor, a*1.91761363636+ma*1.54548810664, a*0.745738636364+ma*0.529863106641,
//        a*0.745738636364+ma*0.622605363985, a*1.91761363636+ma*1.54548810664, a*0.745738636364+ma*0.622605363985,
//        a*1.91761363636+ma*2.87356321839, a*0.5859375+ma*0.57, a*0.5859375+ma*0.57, a*0.5859375+ma*0.5078125);
//    SetOrthotropicTensor(t, subTensor, a*1.91761363636+ma*(1.54548810664-E/5), a*0.745738636364+ma*0.529863106641,
//        a*0.745738636364+ma*0.622605363985, a*1.91761363636+ma*(1.54548810664-E/5), a*0.745738636364+ma*0.622605363985,
//        a*1.91761363636+ma*(2.87356321839+E), a*0.5859375+ma*(0.57-E/10), a*0.5859375+ma*(0.57-E/10), a*0.5859375+ma*(0.5078125-E/10));
    // a*VoronoiCS + (1-a)*VoronoiNoCS
    t[0][0] = a*260.910 +ma*257.89;
    t[0][1] = a*122.880 +ma*125.21;
    t[0][2] = a*147.040 +ma*147.49;
    t[0][3] = a*-0.0011639 +ma*0.064781;
    t[0][4] = a*0.0015323 +ma*-0.014303;
    t[0][5] = a*-0.14353 +ma*0.22157;
    t[1][0] = t[0][1];
    t[1][1] = a*260.910 +ma*257.91;
    t[1][2] = a*147.040 +ma*147.48;
    t[1][3] = a*0.040643 +ma*-0.0021714;
    t[1][4] = a*0.000039308 +ma*-0.00072162;
    t[1][5] = a*0.18457 +ma*-0.24112;
    t[2][0] = t[0][2];
    t[2][1] = t[1][2];
    t[2][2] = a*236.720 +ma*235.19;
    t[2][3] = a*0.00024251 +ma*-0.45645;
    t[2][4] = a*-0.00068456 +ma*0.00069223;
    t[2][5] = a*-0.010199 +ma*-0.00041267;
    t[3][0] = t[0][3];
    t[3][1] = t[1][3];
    t[3][2] = t[2][3];
    t[3][3] = a*105.770 +ma*106.68;
    t[3][4] = a*-0.052229 +ma*-0.0017872;
    t[3][5] = a*0.1147 +ma*-0.10683;
    t[4][0] = t[0][4];
    t[4][1] = t[1][4];
    t[4][2] = t[2][4];
    t[4][3] = t[3][4];
    t[4][4] = a*105.840 +ma*107.01;
    t[4][5] = a*-0.00092836 +ma*0.00061316;
    t[5][0] = t[0][5];
    t[5][1] = t[1][5];
    t[5][2] = t[2][5];
    t[5][3] = t[3][5];
    t[5][4] = t[4][5];
    t[5][5] = a*68.952 +ma*67.691;
    t *= 1000;
    break;
  case PLANE_STRAIN:
  case PLANE_STRESS:
    transIsoType_ = TRANSISO_XZ;
    SetTransIsoMatrix(t, subTensor, a*1.91761363636+ma*1.54548810664, 0.0, 0.0, a*1.91761363636+ma*2.87356321839, a*0.745738636364+ma*0.529863106641, a*0.5859375+ma*0.5078125);
    break;
  default:
    throw Exception("subTensor not implemented yet");
  }

  if (type_ == D_INTERP_IN718_TENSOR || type_ == D_INTERP_IN718_TENSOR_ROT)
  {
    double dens = GetParameter(DesignElement::DENSITY);
    TransferFunction* tf = space_->GetTransferFunction(DesignElement::DENSITY, App::MECH);
    t *= (direction == DesignElement::DENSITY) ? tf->Derivative(dens) : tf->Transform(dens);
  }
  if(type_ == D_INTERP_IN718_TENSOR_ROT){
    // for all rotated types, rotate the material tensor
    LOG_DBG3(dm) << "GetIN718Tensor: tensor before rotation=" << t.ToString();
    LOG_DBG2(dm)<< "GetIN718Tensor: E before rotation = " << t.ToString();
    // RotateTensor needs Hill Mandel matrix
    RotateTensor(mt, direction);
    LOG_DBG2(dm)<< "GetIN718Tensor: E after rotation = " << t.ToString();
  }
}

inline void DesignMaterial::GetLaminatesTensor(MaterialTensor<double>& mt, SubTensorType subTensor, DesignElement::Type direction)
{
  const std::map<DesignElement::Type, double>& map = GetParameters();

  Matrix<double>& t = mt.GetMatrix(VOIGT);

  switch(subTensor)
  {
  case PLANE_STRAIN:    //see Allaire: Shape optimization by the homogenization method, pp. 127 [(2.64),(2.65)]
  {
    t.Resize(3,3);
    t.Init();
    double eps = 5.0625e-4;
    double stiff1 = GetParameter(map, DesignElement::STIFF1);
    double stiff2 = GetParameter(map, DesignElement::STIFF2);
    double E = GetParameter(map, DesignElement::EMODUL);
    double nu = GetParameter(map, DesignElement::POISSON);
    double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
    double mu = E / (2 * (1 + nu));
    Matrix<double> D(3, 3);
    Matrix<double> Dinv(3, 3);
    D.SetEntry(0, 0, 1 / (4 * (mu + lambda)) + 1 / (4 * mu));
    D.SetEntry(0, 1, 1 / (4 * (mu + lambda)) - 1 / (4 * mu));
    D.SetEntry(1, 0, D(0, 1));
    D.SetEntry(1, 1, D(0, 0));
    D.SetEntry(2, 2, 1 / (2 * mu));
    D *= 1 / (eps - 1);
    D.AddToEntry(0, 0, stiff2 / (2 * mu + lambda));
    D.AddToEntry(2, 2, stiff2 / (2 * mu));
    D.AddToEntry(1, 1, stiff1 * (1 - stiff2) / (2 * mu + lambda));
    D.AddToEntry(2, 2, stiff1 * (1 - stiff2) / (2 * mu));
    D.Invert(Dinv);
    switch (direction) {
    case DesignElement::NO_DERIVATIVE:
    case DesignElement::ROTANGLE:
      D.Init();
//      mu /= eps;
//      lambda /= eps;
      D.SetEntry(0, 0, 2 * mu + lambda);
      D.SetEntry(1, 1, D(0, 0));
      D.SetEntry(2, 2, 2 * mu);
      D.SetEntry(0, 1, lambda);
      D.SetEntry(1, 0, lambda);
      t.Add((1 - stiff2) * (1 - stiff1), Dinv);
      t.Add(1.0, D);
      break;
    case DesignElement::STIFF1:
      t.SetEntry(1, 1, 1 / (2 * mu + lambda));
      t.SetEntry(2, 2, 1 / (2 * mu));
      t *= (stiff2 - 1) * (1 - stiff2) * (1 - stiff1);
      Dinv.Mult(t, D);
      D.Mult(Dinv, t);
      t.Add(stiff2-1, Dinv);
      break;
    case DesignElement::STIFF2:
      t.SetEntry(0, 0, 1 / (2 * mu + lambda));
      t.SetEntry(1, 1, -stiff1 / (2 * mu + lambda));
      t.SetEntry(2, 2, (1 - stiff1) / (2 * mu));
      t *= (stiff2 - 1) * (1 - stiff1);
      Dinv.Mult(t, D);
      D.Mult(Dinv, t);
      t.Add(stiff1-1, Dinv);
      break;
    default:
      ZeroMatrix(t, subTensor);
      return;
    }
    break;
  }
  case PLANE_STRESS:      //see Bendsoe, Sigmund: Topology Optimization S. 166
  {
    double E33 = 0.02;
    double stiff1 = GetParameter(map, DesignElement::STIFF1);
    double stiff2 = GetParameter(map, DesignElement::STIFF2);
    double E = GetParameter(map, DesignElement::EMODUL);
    double nu = GetParameter(map, DesignElement::POISSON);
    double n = (stiff2 + stiff1 * stiff2 * (nu * nu - 1) - 1);
    switch (direction) {
    case DesignElement::NO_DERIVATIVE:
    case DesignElement::ROTANGLE:
    case DesignElement::DENSITY: // Treated after switch
    {
      double E11 = -(E*stiff1)/n;
      double E22 = stiff2*E+stiff2*stiff2*nu*nu*E11;
      double E12 = stiff2*nu*E11;
      Set2dMatrix(t, E11, E22, E33, 0.0, 0.0, E12);
      break;
    }
    case DesignElement::STIFF1: {
      double E11 = -(E * (stiff2 - 1)) / (n * n);
      double E22 = stiff2 * stiff2 * nu * nu * E11;
      double E12 = stiff2 * nu * E11;
      Set2dMatrix(t, E11, E22, 0.0, 0.0, 0.0, E12);
      break;
    }
    case DesignElement::STIFF2: {
      double E11 = (E * stiff1 * (stiff1 * (nu * nu - 1) + 1)) / (n * n);
      double E22 = E - 2 * stiff2 * nu * nu * E * stiff1 / n
          + stiff2 * stiff2 * nu * nu * E11;
      double E12 = -nu * E * stiff1 / n + stiff2 * nu * E11;
      Set2dMatrix(t, E11, E22, 0.0, 0.0, 0.0, E12);
    break;
    }
    default:
      ZeroMatrix(t, subTensor);
      return;
    }
    break;
  }
  default:
    throw Exception("subTensor not implemented yet");
  }

  if (type_ == D_LAMINATES)
  {
    double dens = GetParameter(map, DesignElement::DENSITY);
    TransferFunction* tf = space_->GetTransferFunction(DesignElement::DENSITY, App::MECH);
    t *= (direction == DesignElement::DENSITY) ? tf->Derivative(dens) : tf->Transform(dens);
  }

  double rotAngle = GetParameter(map, DesignElement::ROTANGLE);
  LOG_DBG2(dm)<< "GetLaminatesTensor: E before rotation = " << t.ToString();
  RotateTensor(mt, direction, true, rotAngle);
  LOG_DBG2(dm)<< "GetLaminatesTensor: E after rotation = " << t.ToString();
  return;
}

inline void DesignMaterial::ZeroMatrix(Matrix<double>& t, SubTensorType subTensor) {
  assert(Optimization::context->ToApp() == App::MECH);
  assert(subTensor != NO_TENSOR);
  switch (subTensor) {
  case FULL:
    t.Resize(6, 6);
    LOG_DBG(dm)<<"Zero Tensor: "<<t.ToString();
    break;
  case PLANE_STRAIN:
  case PLANE_STRESS:
  case PLANE:
    t.Resize(3, 3);
    break;
  default:
    throw Exception("subTensor not implemented yet");
  }
  t.Init();
}

inline void DesignMaterial::Set2dMatrix(Matrix<double>& t, double t11, double t22,
    double t33, double t23, double t13, double t12) {
  t.Resize(3, 3);
  t.Init();
  t[0][0] = t11;
  t[0][1] = t12;
  t[0][2] = t13;
  t[1][0] = t12;
  t[1][1] = t22;
  t[1][2] = t23;
  t[2][0] = t13;
  t[2][1] = t23;
  t[2][2] = t33;
}

inline void DesignMaterial::Set2dMatrix(Matrix<double>& t, double t11, double t12, double t13, double t21, double t22, double t23, double t31, double t32, double t33) {
  t.Resize(3, 3);
  t.Init();
  t[0][0] = t11;
  t[0][1] = t12;
  t[0][2] = t13;
  t[1][0] = t21;
  t[1][1] = t22;
  t[1][2] = t23;
  t[2][0] = t31;
  t[2][1] = t32;
  t[2][2] = t33;
}

inline void DesignMaterial::Set3dMatrix(Matrix<double>& t,
    SubTensorType subTensor, double e11, double e12, double e13, double e14, double e15, double e16, double e22, double e23, double e24, double e25, double e26,
    double e33, double e34, double e35, double e36, double e44, double e45, double e46, double e55, double e56, double e66) {
  assert(subTensor == FULL);
  switch (subTensor) {
    case FULL:
      t.Resize(6, 6);
      t.Init();
      t[0][0] = e11;
      t[0][1] = e12;
      t[0][2] = e13;
      t[0][3] = e14;
      t[0][4] = e15;
      t[0][5] = e16;
      t[1][0] = e12;
      t[2][0] = e13;
      t[3][0] = e14;
      t[4][0] = e15;
      t[5][0] = e16;
      t[1][1] = e22;
      t[1][2] = e23;
      t[1][3] = e24;
      t[1][4] = e25;
      t[1][5] = e26;
      t[2][1] = e23;
      t[3][1] = e24;
      t[4][1] = e25;
      t[5][1] = e26;
      t[2][2] = e33;
      t[2][3] = e34;
      t[2][4] = e35;
      t[2][5] = e36;
      t[3][2] = e34;
      t[4][2] = e35;
      t[5][2] = e36;
      t[3][3] = e44;
      t[3][4] = e45;
      t[3][5] = e46;
      t[4][3] = e45;
      t[5][3] = e46;
      t[4][4] = e55;
      t[4][5] = e56;
      t[5][4] = e56;
      t[5][5] = e66;
      break;
    case PLANE_STRAIN:
    case PLANE_STRESS:
      SetTransIsoMatrix(t, subTensor, e11, 0.0, 0.0, e22, e12, e66);
      break;
    default:
      throw Exception("subTensor not implemented yet");
    }
}

inline void DesignMaterial::SetTransIsoMatrix(Matrix<double>& t, SubTensorType subTensor, double iD, double inD, double iG, double oD, double onD, double oG) {
  switch (subTensor) {
  case FULL:
    t.Resize(6, 6);
    t.Init();
    switch (transIsoType_) {
    case TRANSISO_XY:
      t[0][0] = iD;
      t[0][1] = inD;
      t[0][2] = onD;
      t[1][0] = inD;
      t[1][1] = iD;
      t[1][2] = onD;
      t[2][0] = onD;
      t[2][1] = onD;
      t[2][2] = oD;
      t[3][3] = oG;
      t[4][4] = oG;
      t[5][5] = iG;
      break;
    case TRANSISO_YZ:
      t[0][0] = oD;
      t[0][1] = onD;
      t[0][2] = onD;
      t[1][0] = onD;
      t[1][1] = iD;
      t[1][2] = inD;
      t[2][0] = onD;
      t[2][1] = inD;
      t[2][2] = iD;
      t[3][3] = iG;
      t[4][4] = oG;
      t[5][5] = oG;
      break;
    case TRANSISO_XZ:
      t[0][0] = iD;
      t[0][1] = onD;
      t[0][2] = inD;
      t[1][0] = onD;
      t[1][1] = oD;
      t[1][2] = onD;
      t[2][0] = inD;
      t[2][1] = onD;
      t[2][2] = iD;
      t[3][3] = oG;
      t[4][4] = iG;
      t[5][5] = oG;
      break;
    }
    break;
  case PLANE_STRAIN:
  case PLANE_STRESS:
    t.Resize(3, 3);
    t.Init();
    switch (transIsoType_) {
    case TRANSISO_XY:
      t[0][0] = iD;
      t[0][1] = inD;
      t[1][0] = inD;
      t[1][1] = iD;
      t[2][2] = iG;
      break;
    case TRANSISO_YZ:
      t[0][0] = oD;
      t[0][1] = onD;
      t[1][0] = onD;
      t[1][1] = iD;
      t[2][2] = oG;
      break;
    case TRANSISO_XZ:
      t[0][0] = iD;
      t[0][1] = onD;
      t[1][0] = onD;
      t[1][1] = oD;
      t[2][2] = oG;
      break;
    }
    break;
  default:
    throw Exception("subTensor not implemented yet");
  }
}

inline void DesignMaterial::SetIsoMatrix(Matrix<double>& t, SubTensorType subTensor,
    double D, double nd, double G) {
  SetTransIsoMatrix(t, subTensor, D, nd, G, D, nd, G);
}

void DesignMaterial::RotateTensor(MaterialTensor<double>& mt, DesignElement::Type direction, bool angles, double rx, double ry, double rz){
  // rotation matrix is found in Dissertation of B. Schmidt: Topology Preserving Multi-Layer Shape and Material Optimization p. 62
  // and also found in Wikipedia Drehmatrix (german)
  // rotates the material by ROTANGLEFIRST around the first (z-)axis, by ROTANGLESECOND around the second (y-)axis and by ROTANGLETHIRD around the third (x-)axis in this given order or rz,ry,rx
  // direction of rotation around an axis is positive (ccw), i.e. right hand rule applies
  // this is identical to BaseMaterial::RotateTensorByRotationAngles

  Matrix<double>& t = mt.GetMatrix(VOIGT);

  LOG_DBG3(dm) << "RT: input tensor=" << t.ToString() << " d=" << direction  << " angles=" << angles << " rx=" << rx << " ry=" << ry << " rz=" << rz;

  int dim = t.GetNumRows() > 3 ? 3 : 2;

  double theta3 = 0.0, theta2 = 0.0, theta1 = 0.0;
  if(dim == 3){
    if (angles) {
      theta3 = rx;
      theta2 = ry;
      theta1 = rz;
    } else {
      theta3 = GetParameter(DesignElement::ROTANGLETHIRD);
      theta2 = GetParameter(DesignElement::ROTANGLESECOND);
      theta1 = 0.0;
      if(HasParameter(DesignElement::ROTANGLEFIRST))
        theta1 = GetParameter(DesignElement::ROTANGLEFIRST);
    }
  }else{ // dim == 2
    if (angles) {
      // this is correct
      theta1 = rx;
    } else {
      theta1 = GetParameter(DesignElement::ROTANGLE);
    }
  }

  Matrix<Double> R(dim, dim);
  SetRotationMatrix(R, theta1, theta2, theta3);

  // see also baseMaterial.cc for this
  int dimQ = dim == 3 ? 6 : 3;
  Matrix<Double> Q(dimQ, dimQ);
  Q.Resize(dimQ, dimQ);
  int l = dimQ-1;
  Q[0][0] = R[0][0]*R[0][0];
  Q[0][1] = R[0][1]*R[0][1];
  Q[0][l] = 2.0*R[0][0]*R[0][1];
  Q[1][0] = R[1][0]*R[1][0];
  Q[1][1] = R[1][1]*R[1][1];
  Q[1][l] = 2.0*R[1][0]*R[1][1];
  Q[l][0] = R[0][0]*R[1][0];
  Q[l][1] = R[0][1]*R[1][1];
  Q[l][l] = R[0][0]*R[1][1] + R[0][1]*R[1][0];

  if(dim == 3){
    Q[0][2] = R[0][2]*R[0][2];
    Q[0][3] = 2.0*R[0][1]*R[0][2];
    Q[0][4] = 2.0*R[0][0]*R[0][2];
    Q[1][2] = R[1][2]*R[1][2];
    Q[1][3] = 2.0*R[1][1]*R[1][2];
    Q[1][4] = 2.0*R[1][0]*R[1][2];
    Q[2][0] = R[2][0]*R[2][0];
    Q[2][1] = R[2][1]*R[2][1];
    Q[2][2] = R[2][2]*R[2][2];
    Q[2][3] = 2.0*R[2][1]*R[2][2];
    Q[2][4] = 2.0*R[2][0]*R[2][2];
    Q[2][5] = 2.0*R[2][0]*R[2][1];
    Q[3][0] = R[1][0]*R[2][0];
    Q[3][1] = R[1][1]*R[2][1];
    Q[3][2] = R[1][2]*R[2][2];
    Q[3][3] = R[1][1]*R[2][2] + R[1][2]*R[2][1];
    Q[3][4] = R[1][0]*R[2][2] + R[1][2]*R[2][0];
    Q[3][5] = R[1][0]*R[2][1] + R[1][1]*R[2][0];
    Q[4][0] = R[0][0]*R[2][0];
    Q[4][1] = R[0][1]*R[2][1];
    Q[4][2] = R[0][2]*R[2][2];
    Q[4][3] = R[0][1]*R[2][2] + R[0][2]*R[2][1];
    Q[4][4] = R[0][0]*R[2][2] + R[0][2]*R[2][0];
    Q[4][5] = R[0][0]*R[2][1] + R[0][1]*R[2][0];
    Q[5][2] = R[0][2]*R[1][2];
    Q[5][3] = R[0][1]*R[1][2] + R[0][2]*R[1][1];
    Q[5][4] = R[0][0]*R[1][2] + R[0][2]*R[1][0];
  }
  LOG_DBG3(dm) << "RT: Corresponding Q is " << Q.ToString();
  if(direction != DesignElement::ROTANGLETHIRD && direction != DesignElement::ROTANGLESECOND && direction != DesignElement::ROTANGLEFIRST && direction != DesignElement::ROTANGLE) {
    // calculate Q*t*Q' and store back to t. unfortunately MultT is the wrong way
    Matrix<double> help(dimQ, dimQ);
    Q.Mult(t, help); // help = Q * t
    Matrix<double> QT(dimQ, dimQ);
    Q.Transpose(QT); // QT = Q^T
    help.Mult(QT, t); // t = help * QT == Q * t * Q^T
    LOG_DBG3(dm) << "RT: final tensor " << t.ToString();
  } else { // we need a derivative
    Matrix<double> dR(dim, dim);
    SetRotationMatrix(dR, theta1, theta2, theta3, direction); // this now produces the derivative

    Matrix<double> dQ(dimQ, dimQ);
    // this part can be produced from the definition of Q above by sed 's/Q/dQ;s/R\(\[\d\]\[\d\]\)\*R\(\[\d\]\[\d\]\)/(dR\1*R\2+R\1*dR\2)/g', effectively using the product rule
    dQ[0][0] = (dR[0][0]*R[0][0]+R[0][0]*dR[0][0]);
    dQ[0][1] = (dR[0][1]*R[0][1]+R[0][1]*dR[0][1]);
    dQ[0][l] = 2.0*(dR[0][0]*R[0][1]+R[0][0]*dR[0][1]);
    dQ[1][0] = (dR[1][0]*R[1][0]+R[1][0]*dR[1][0]);
    dQ[1][1] = (dR[1][1]*R[1][1]+R[1][1]*dR[1][1]);
    dQ[1][l] = 2.0*(dR[1][0]*R[1][1]+R[1][0]*dR[1][1]);
    dQ[l][0] = (dR[0][0]*R[1][0]+R[0][0]*dR[1][0]);
    dQ[l][1] = (dR[0][1]*R[1][1]+R[0][1]*dR[1][1]);
    dQ[l][l] = (dR[0][0]*R[1][1]+R[0][0]*dR[1][1]) + (dR[0][1]*R[1][0]+R[0][1]*dR[1][0]);

    if(dim == 3) {
      dQ[0][2] = (dR[0][2]*R[0][2]+R[0][2]*dR[0][2]);
      dQ[0][3] = 2.0*(dR[0][1]*R[0][2]+R[0][1]*dR[0][2]);
      dQ[0][4] = 2.0*(dR[0][0]*R[0][2]+R[0][0]*dR[0][2]);
      dQ[1][2] = (dR[1][2]*R[1][2]+R[1][2]*dR[1][2]);
      dQ[1][3] = 2.0*(dR[1][1]*R[1][2]+R[1][1]*dR[1][2]);
      dQ[1][4] = 2.0*(dR[1][0]*R[1][2]+R[1][0]*dR[1][2]);
      dQ[2][0] = (dR[2][0]*R[2][0]+R[2][0]*dR[2][0]);
      dQ[2][1] = (dR[2][1]*R[2][1]+R[2][1]*dR[2][1]);
      dQ[2][2] = (dR[2][2]*R[2][2]+R[2][2]*dR[2][2]);
      dQ[2][3] = 2.0*(dR[2][1]*R[2][2]+R[2][1]*dR[2][2]);
      dQ[2][4] = 2.0*(dR[2][0]*R[2][2]+R[2][0]*dR[2][2]);
      dQ[2][5] = 2.0*(dR[2][0]*R[2][1]+R[2][0]*dR[2][1]);
      dQ[3][0] = (dR[1][0]*R[2][0]+R[1][0]*dR[2][0]);
      dQ[3][1] = (dR[1][1]*R[2][1]+R[1][1]*dR[2][1]);
      dQ[3][2] = (dR[1][2]*R[2][2]+R[1][2]*dR[2][2]);
      dQ[3][3] = (dR[1][1]*R[2][2]+R[1][1]*dR[2][2]) + (dR[1][2]*R[2][1]+R[1][2]*dR[2][1]);
      dQ[3][4] = (dR[1][0]*R[2][2]+R[1][0]*dR[2][2]) + (dR[1][2]*R[2][0]+R[1][2]*dR[2][0]);
      dQ[3][5] = (dR[1][0]*R[2][1]+R[1][0]*dR[2][1]) + (dR[1][1]*R[2][0]+R[1][1]*dR[2][0]);
      dQ[4][0] = (dR[0][0]*R[2][0]+R[0][0]*dR[2][0]);
      dQ[4][1] = (dR[0][1]*R[2][1]+R[0][1]*dR[2][1]);
      dQ[4][2] = (dR[0][2]*R[2][2]+R[0][2]*dR[2][2]);
      dQ[4][3] = (dR[0][1]*R[2][2]+R[0][1]*dR[2][2]) + (dR[0][2]*R[2][1]+R[0][2]*dR[2][1]);
      dQ[4][4] = (dR[0][0]*R[2][2]+R[0][0]*dR[2][2]) + (dR[0][2]*R[2][0]+R[0][2]*dR[2][0]);
      dQ[4][5] = (dR[0][0]*R[2][1]+R[0][0]*dR[2][1]) + (dR[0][1]*R[2][0]+R[0][1]*dR[2][0]);
      dQ[5][2] = (dR[0][2]*R[1][2]+R[0][2]*dR[1][2]);
      dQ[5][3] = (dR[0][1]*R[1][2]+R[0][1]*dR[1][2]) + (dR[0][2]*R[1][1]+R[0][2]*dR[1][1]);
      dQ[5][4] = (dR[0][0]*R[1][2]+R[0][0]*dR[1][2]) + (dR[0][2]*R[1][0]+R[0][2]*dR[1][0]);
    }
    LOG_DBG3(dm) << "RT: Corresponding dQ is " << dQ.ToString();

    // we now, have to calculate dQ*t*Q' + Q*t*dQ'
    // we calculate dQ*t*Q' + (dQ*t*Q')' = dQ*t*Q' + Q''*t'*dQ' = dQ*t*Q' + Q*t*dQ'
    Matrix<Double> help(dimQ, dimQ);
    dQ.Mult(t, help); // help = dQ * t
    Q.Transpose(dQ); // dQ is no longer needed, we overwrite it
    help.Mult(dQ, t);
    t.Transpose(help);
    t.Add(1.0, help);    // and add the rest
    //FIXME: this section is ugly and should be fixed if expression templates work reliably
  }

}
void DesignMaterial::SetOneAxisRotationMatrix(Matrix<double>& R, double theta, int axis, bool derivative) {
  // rotation matrix in 2d around z axis or in 3d around chosen coordinate axis

  // see https://en.wikipedia.org/wiki/Rotation_matrix

  double stheta, ctheta, dx;
  if(!derivative) {
    stheta = sin(theta);
    ctheta = cos(theta);
    dx = 1;
  }  else {
    stheta = cos(theta);
    ctheta = -sin(theta);
    dx = 0;
  }

  Matrix<double> Q;
  Q.Resize(2,2);
  Q[0][0] = ctheta;
  Q[0][1] = -stheta;
  Q[1][0] = stheta;
  Q[1][1] = ctheta;

  R.Resize(dim, dim);
  if (dim == 2)
    R = Q;
  else {
    // For three dimensions shift entries modulo 3
    R.Init();
    R[axis][axis] = dx;
    for (int row = 0; row < 2; ++row) {
      int new_row = (row+axis+1) % 3;
      for (int col = 0; col < 2; ++col) {
        int new_col = (col+axis+1) % 3;
        R[new_row][new_col] = Q[row][col];
      }
    }
  }
  LOG_DBG2(dm) << "SOARM: rotation matrix axis=" << axis << " angle=" << theta << " d=" << derivative << " Q=" << Q.ToString() << " -> R=" << R.ToString();
}

void DesignMaterial::SetRotationMatrix(Matrix<double>& R, double theta1, double theta2, double theta3, DesignElement::Type direction) {
  // rotation axes ares given by rotationType (default XYZ, i.e. first z, then y, then x)
  // direction of rotation around an axis is positive (ccw), i.e. right hand rule applies
  R.Resize(dim, dim);

  if(dim == 2) {
    SetOneAxisRotationMatrix(R, theta1, NONE, direction == DesignElement::ROTANGLE);
  } else {
    Matrix<double> R1, R2, R3;
    R1.Resize(dim, dim);
    R2.Resize(dim, dim);
    R3.Resize(dim, dim);

    int axis1 = -1, axis2 = -1, axis3 = -1;
    switch(rotationType_) {
    case RotationType::ZXZ:
      axis3 = 2;
      axis2 = 0;
      axis1 = 2;
      break;
    case RotationType::ZYZ:
      axis3 = 2;
      axis2 = 1;
      axis1 = 2;
      break;
    case RotationType::YZY:
      axis3 = 1;
      axis2 = 2;
      axis1 = 1;
      break;
    case RotationType::YXY:
      axis3 = 1;
      axis2 = 0;
      axis1 = 1;
      break;
    case RotationType::XYX:
      axis3 = 0;
      axis2 = 1;
      axis1 = 0;
      break;
    case RotationType::XZX:
      axis3 = 0;
      axis2 = 2;
      axis1 = 0;
      break;
    case RotationType::XYZ:
      axis3 = 0;
      axis2 = 1;
      axis1 = 2;
      break;
    case RotationType::YXZ:
      axis3 = 0;
      axis2 = 1;
      axis1 = 2;
      break;
    case RotationType::XZY:
      axis3 = 0;
      axis2 = 2;
      axis1 = 1;
      break;
    case RotationType::ZXY:
      axis3 = 2;
      axis2 = 0;
      axis1 = 1;
      break;
    case RotationType::ZYX:
      axis3 = 2;
      axis2 = 1;
      axis1 = 0;
      break;
    case RotationType::YZX:
      axis3 = 1;
      axis2 = 2;
      axis1 = 0;
      break;
    default: // same as RotationType::XYZ
      axis3 = 0;
      axis2 = 1;
      axis1 = 2;
      break;
    }

    assert((axis1 != -1) && (axis2 != -1) && (axis3 != -1));

    SetOneAxisRotationMatrix(R1, theta1, axis1, direction == DesignElement::ROTANGLEFIRST);
    SetOneAxisRotationMatrix(R2, theta2, axis2, direction == DesignElement::ROTANGLESECOND);
    SetOneAxisRotationMatrix(R3, theta3, axis3, direction == DesignElement::ROTANGLETHIRD);

    // we use R as temporary cache: R=R2*R1, R1=R3*R, R<-R1
    R2.Mult(R1,R);
    R3.Mult(R,R1);
    R = R1;
  }
  LOG_DBG3(dm)  << "SRM:Rotation matrix for t1=" << theta1 << ", t2=" << theta2 << ", t3=" << theta3 << " with derivative w.r.t. to " << direction << " is \n" << R.ToString();
//  std::cout << "SRM:Rotation matrix for t1=" << theta1 << ", t2=" << theta2 << ", t3=" << theta3 << " with derivative w.r.t. to " << direction << " is \n" << R.ToString() << std::endl;
}

void DesignMaterial::RotatePiezoCouplingTensor(Matrix<double>& E, double phi, DesignElement::Type direction)
{
  // R(phi) * [e] * Q(phi)^T
  // derivative: dR(phi)/dphi * ([e] * Q(phi)^T) + R(phi) * ([e] * dQ(phi)/dphi)^T

  // rotation is clockwise

  Matrix<double> R(2,2);
  R[0][0] = cos(phi);
  R[0][1] = sin(phi);
  R[1][0] = -R[0][1];
  R[1][1] = R[0][0];

  Matrix<double> QT(3,3);

  QT[0][0] = R[0][0]*R[0][0];
  QT[0][1] = R[1][0]*R[1][0];
  QT[0][2] = R[0][0]*R[1][0];

  QT[1][0] = R[0][1]*R[0][1];
  QT[1][1] = R[1][1]*R[1][1];
  QT[1][2] = R[0][1]*R[1][1];

  QT[2][0] = 2.0*R[0][0]*R[0][1];
  QT[2][1] = 2.0*R[1][0]*R[1][1];
  QT[2][2] = R[0][0]*R[1][1] + R[0][1]*R[1][0];

  Matrix<double> help(2,3);
  E.Mult(QT, help); // help = E * Q^T
  if(direction != DesignElement::ROTANGLE)
  {
    R.Mult(help, E); // E = R * (E * Q^T)

    LOG_DBG3(dm) << "RPCT phi=" << phi << " R=" << R.ToString() << " QT=" << QT.ToString();
    return;
  }
  else
  {
    Matrix<double> dR(2,2);
    dR[0][0] = -sin(phi);
    dR[0][1] = cos(phi);
    dR[1][0] = -cos(phi);
    dR[1][1] = -sin(phi);

    Matrix<double> dQT(3,3);

    dQT[0][0] = -2.0*R[0][0]*R[0][1];   // -2 cos(a) sin(a)
    dQT[0][1] = -dQT[0][0];             // 2 cos(a) sin(a)
    dQT[0][2] = -1.0*cos(2.0*phi);      // sin(a)^2 - cos(a)^2

    dQT[1][0] = dQT[0][1];              // 2 cos(a) sin(a)
    dQT[1][1] = dQT[0][0];              // -2 cos(a) sin(a)
    dQT[1][2] = -dQT[0][2];             // cos(a)^2 - sin(a)^2

    dQT[2][0] = 2.0*dQT[1][2];          // 2 cos(a)^2 - 2 sin(a)^2
    dQT[2][1] = 2.0*dQT[0][2];          // 2 sin(a)^2 - 2 cos(a)^2
    dQT[2][2] = 2.0 * dQT[0][0];        // -4 cos(a) sin(a)

    Matrix<double> left(2,3);
    dR.Mult(help, left); // left = dR * (E * Q^T)

    E.Mult(dQT, help); // help = E * dQ^T

    Matrix<double> right(2,3);
    R.Mult(help, right); // right = R * (help) = R * (E * dQ^T)
    E = left + right;
    LOG_DBG3(dm) << "RPCT phi=" << phi << " R=" << R.ToString() << " dR=" << dR.ToString() << " QT=" << QT.ToString() << " dQT=" << dQT.ToString();
    return;
  }
}


void DesignMaterial::RotateElecTensor(MaterialTensor<double>& mt, double phi, DesignElement::Type direction)
{
  // R(phi) * [e] * R(phi)^T
  // derivative: dR(phi)/dphi * ([e] * R(phi)^T) + R(phi) * ([e] * dR(phi)/dphi)^T

  // rotation is counterclockwise

  Matrix<double>& E = mt.GetMatrix(NO_NOTATION);

  Matrix<double> RT(2,2);
  RT[0][0] = cos(phi);
  RT[0][1] = -sin(phi);
  RT[1][0] = -RT[0][1];
  RT[1][1] = RT[0][0];

  Matrix<double> help(2,2);
  E.Mult(RT, help); // help = E * R^T

  if(direction != DesignElement::ROTANGLE)
  {
    RT.MultT(help, E); // E = R * (E * R^T)
    return;
  }
  else
  {
    Matrix<double> dRT(2,2);
    dRT[0][0] = -sin(phi);
    dRT[0][1] = -cos(phi);
    dRT[1][0] = -dRT[0][1];
    dRT[1][1] = dRT[0][0];

    Matrix<double> left(2,2);
    dRT.MultT(help, left); // left = dR * (E * R^T)

    E.Mult(dRT, help); // help = E * dR^T

    RT.MultT(help, dRT); // overwrite dR to use temporary: dR = R * (help) = R * (E * dR^T)
    E = left + dRT;
    return;
  }
}


inline double DesignMaterial::GetTransIsoMass(double iD, double iG, double oD, double oG){
  switch(dim){
  case 2:
    switch (transIsoType_) {
    case TRANSISO_XY:
      return (2.0 * iD + iG);
    case TRANSISO_YZ:
    case TRANSISO_XZ:
      return (iD + oD + oG);
    default:
      throw Exception("transIsoType not implemented yet");
    }
    break;
  case 3:
    return (2.0 * iD + oD + iG + 2.0 * oG);
  default:
    throw Exception("strange dimension");
  }
}

inline double DesignMaterial::GetIsoMass(double D, double G) {
  return (GetTransIsoMass(D, G, D, G));
}

bool DesignMaterial::GetTensor(MaterialTensor<double>& mt, DesignElement::Type type, SubTensorType subTensor, const Elem* elem, DesignElement::Type direction, MaterialTensorNotation notation, bool core)
{
  switch(type)
  {
  case DesignElement::MECH_TRACE:
  case DesignElement::MECH_ALL:
  case DesignElement::ALL_DESIGNS:
  {
    bool ret = GetMechTensor(mt, subTensor, elem, direction, core);
    assert(mt.GetNotation() == VOIGT);
    if (notation == HILL_MANDEL)
      mt.ToHillMandel();
    return ret;
  }
  case DesignElement::DIELEC_TRACE:
  case DesignElement::DIELEC_ALL:
    if(direction == DesignElement::NO_DERIVATIVE && !HasParameter(DesignElement::DIELEC_11))
      return false;
    return GetElecTensor(mt, elem, direction);
  case DesignElement::PIEZO_ALL:
    if(direction == DesignElement::NO_DERIVATIVE && !HasParameter(DesignElement::PIEZO_11))
        return false;
    return GetPiezoCouplingTensor(mt, elem, direction);
  default:
    return false;
    break;
  }
  return false;
}

bool DesignMaterial::GetMechTensor(MaterialTensor<Complex>& mtc, SubTensorType subTensor, const Elem* elem, DesignElement::Type direction, bool core, const CoefFunctionOpt* coef)
{
  assert(mtc.GetNotation() == VOIGT);

  // we assume we have no complex material (special form of damping)
  MaterialTensor<double> mtd(VOIGT);
  if(!GetMechTensor(mtd, subTensor, elem, direction, core))
    return false;

  Matrix<double>&  md = mtd.GetMatrix(VOIGT);
  Matrix<Complex>& mc = mtc.GetMatrix(VOIGT);

  mc.Resize(md.GetNumRows(), md.GetNumCols());
  mc.SetPart(Global::REAL, md, true); // zero other part
  return true;
}

bool DesignMaterial::GetMechTensor(MaterialTensor<double>& mt, SubTensorType subTensor, const Elem* elem, DesignElement::Type direction, bool core, const CoefFunctionOpt* coef)
{
  // FIXME! Check whether assertion makes sense
  //assert(!(notation == HILL_MANDEL && type_ != FMO && type_ != LAMINATES && type_ != D_LAMINATES && type_ != HOM_RECT && type_ != D_HOM_RECT && type_ != HOM_RECT_C1 && type_ != HOM_ISO_C1  && type_ !=  DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC && type_ != DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED && type_ != ORTHOTROPIC && type_ != DENSITY_TIMES_ROT_PA12));
  // FIXME!! with parallel assembling GetMechTensor seems to be not thread save
  // make the code save and remove the lock in calling DesingSpace!
  if(!CollectMaterialParametersForElement(space_, elem))
    return false;

  // optimization only debugged for specific design material cases
  if (domain->GetOptimization()->GetOptimizerType() == Optimization::SGP_SOLVER && !(type_ == SGP_GRADIENTCHECK || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_2D_TENSOR || type_ == DENSITY_TIMES_ROTATED_2D_TENSOR))
    EXCEPTION("SGP solver only verified for specific design material cases!");

  switch (type_) {
  case FMO:
  case SGP_MATLAB:
  case SGP_GRADIENTCHECK:
    GetElasticFMOTensor(mt, subTensor, direction);
    break;

  case ORTHOTROPIC:
  case DENSITY_TIMES_ORTHOTROPIC:
    GetOrthotropicMaterialTensor(mt, subTensor, direction);
    break;
  case ISOTROPIC:
    GetIsoMaterialTensor(mt, subTensor, direction);
    break;
  case LAME_ISOTROPIC: // LAME_ISOTROPIC
    GetLameMaterialTensor(mt, subTensor, direction);
    break;
  case TRANSVERSAL_ISOTROPIC:
  case TRANSVERSAL_ISOTROPIC_BOXED:
  case DENSITY_TIMES_TRANSVERSAL_ISOTROPIC:
  case DENSITY_TIMES_TRANSVERSAL_ISOTROPIC_BOXED:
  case DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC:
  case DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED:
  case DENSITY_TIMES_ROT_PA12:
    GetTransIsoMaterialTensor(mt, subTensor, direction, core);
    break;
  case DENSITY_TIMES_2D_TENSOR:
  case DENSITY_TIMES_2D_TENSOR_CONSTANT_TRACE:
  case DENSITY_TIMES_ROTATED_2D_TENSOR:
    GetDensityTimes2dTensorTensor(mt, subTensor, direction, core);
    break;
  case LAMINATES:
  case D_LAMINATES:
    GetLaminatesTensor(mt, subTensor, direction);
    break;
  case HOM_RECT:
  case D_HOM_RECT:
  case HOM_RECT_C1:
  case HOM_ISO_C1:
  case HEAT:
    GetInterpolatedHomTensor(mt, subTensor, elem, direction);
    break;
  case D_INTERP_IN718_TENSOR:
  case D_INTERP_IN718_TENSOR_ROT:
    GetIN718Tensor(mt, subTensor, direction);
    break;
  case MSFEM_C1:
  {
    Matrix<double>& t = mt.GetMatrix(VOIGT);
    ZeroMatrix(t, subTensor);
    break;
  }
  case FEATURE_MAPPING_ANISO:
    assert(dynamic_cast<FeatureMappingDesign*>(space_) != nullptr);
    dynamic_cast<FeatureMappingDesign*>(space_)->GetAnisoMechTensor(elem, mt, direction, coef, this);
    break;
  case NO_TYPE:
    assert(false);
    break;
  }

  LOG_DBG2(dm) << "GMT: e=" << elem->elemNum << " t=" << type.ToString(type_) << " d=" << DesignElement::type.ToString(direction) << " p=" << core << " -> " << mt.GetMatrix(VOIGT).ToString();

  assert(mt.GetMatrix(VOIGT).GetNumRows() >= 3 && mt.GetMatrix(VOIGT).GetNumCols() >= 3);
  
  return true;
}


bool DesignMaterial::GetElecTensor(MaterialTensor<double>& mt, const Elem* elem, DesignElement::Type direction)
{
  if(!CollectMaterialParametersForElement(space_, elem))
    return false;

  // only 2D!
  bool set = direction == DesignElement::NO_DERIVATIVE || direction == DesignElement::ROTANGLE;

  double e11 = set ? GetParameter(DesignElement::DIELEC_11) : 0;
  double e22 = set ? GetParameter(DesignElement::DIELEC_22) : 0;
  double e12 = set ? GetParameter(DesignElement::DIELEC_12) : 0;
  double rotAngle = set ? GetParameter(DesignElement::ROTANGLE) : 0;

  Matrix<double>& E = mt.GetMatrix(NO_NOTATION);

  E.Resize(2,2);
  E.Init();

  switch(direction)
  {
  case DesignElement::NO_DERIVATIVE:
  case DesignElement::ROTANGLE:
    // negative for the piezo case
    E[0][0] = -e11; E[0][1] = -e12;
    E[1][0] = -e12; E[1][1] = -e22;
    break;
  case DesignElement::DIELEC_11:
    E[0][0] = -1.0;
    break;
  case DesignElement::DIELEC_22:
    E[1][1] = -1.0;
    break;
  case DesignElement::DIELEC_12:
    E[0][1] = -1.0;
    break;
  default:
    // sensitivity is zero!
    break;
  }

  LOG_DBG2(dm) << "GET: E before rotation = " << E.ToString();
  RotateElecTensor(mt, rotAngle, direction);
  LOG_DBG2(dm) << "GET: E after rotation =  " << E.ToString();

  return true;
}


bool DesignMaterial::GetPiezoCouplingTensor(MaterialTensor<double>& mt, const Elem* elem, DesignElement::Type direction)
{
  if(!CollectMaterialParametersForElement(space_, elem))
    return false;

  const std::map<DesignElement::Type, double>& map = GetParameters();

  // only 2D!
  bool set = direction == DesignElement::NO_DERIVATIVE || direction == DesignElement::ROTANGLE;
  double e11 = set ? GetParameter(map, DesignElement::PIEZO_11) : 0;
  double e12 = set ? GetParameter(map, DesignElement::PIEZO_12) : 0;
  double e13 = set ? GetParameter(map, DesignElement::PIEZO_13) : 0;
  double e21 = set ? GetParameter(map, DesignElement::PIEZO_21) : 0;
  double e22 = set ? GetParameter(map, DesignElement::PIEZO_22) : 0;
  double e23 = set ? GetParameter(map, DesignElement::PIEZO_23) : 0;
  double rotAngle = set ? GetParameter(map, DesignElement::ROTANGLE) : 0;

  Matrix<double>& E = mt.GetMatrix(NO_NOTATION);
  E.Resize(2,3);
  E.Init();

  switch(direction)
  {
  case DesignElement::NO_DERIVATIVE:
  case DesignElement::ROTANGLE:
    E[0][0] = e11; E[0][1] = e12; E[0][2] = e13;
    E[1][0] = e21; E[1][1] = e22; E[1][2] = e23;
    break;
  case DesignElement::PIEZO_11:
    E[0][0] = 1.0;
    break;
  case DesignElement::PIEZO_12:
    E[0][1] = 1.0;
    break;
  case DesignElement::PIEZO_13:
    E[0][2] = 1.0;
    break;
  case DesignElement::PIEZO_21:
    E[1][0] = 1.0;
    break;
  case DesignElement::PIEZO_22:
    E[1][1] = 1.0;
    break;
  case DesignElement::PIEZO_23:
    E[1][2] = 1.0;
    break;

  default:
    // empty, sensitivity is zero
    break;
  }

  LOG_DBG2(dm) << "GPCT: E before rotation = " << E.ToString() << " ra=" << rotAngle << " d=" << DesignElement::type.ToString(direction);
  RotatePiezoCouplingTensor(E, rotAngle, direction);
  LOG_DBG2(dm) << "GPCT: E after rotation =  " << E.ToString();

  return true;
}

double DesignMaterial::GetMechMass(const Elem* elem, DesignElement::Type direction)
{
  if(!CollectMaterialParametersForElement(space_, elem))
    throw Exception("no mass data found");

  if(massIsDesign_)
  {
    switch(direction)
    {
    case DesignElement::MASS:
      return massFactor_;
    case DesignElement::NO_DERIVATIVE:
      return GetParameter(DesignElement::MASS) * massFactor_;
    default:
      return 0.0;
    }
  }
  else
  {
    switch (type_)
    {
    case ISOTROPIC:
      return GetIsoMaterialMass(direction) * massFactor_;
    case LAME_ISOTROPIC: // LAME_ISOTROPIC
      return GetLameMaterialMass(direction) * massFactor_;
    case TRANSVERSAL_ISOTROPIC:
    case TRANSVERSAL_ISOTROPIC_BOXED:
      return GetTransIsoMaterialMass(direction) * massFactor_;
    case DENSITY_TIMES_ROT_PA12:
      return (GetDensityTimesTensorMass(direction) * massFactor_);
    default: // case default
      throw Exception("DesignMaterial Type not implemented yet");
    }
  }
}

bool DesignMaterial::GetMaterialDamping(double& alpha, double& beta, DesignElement::Type direction)
{
  if (DampingIsDesign()) {
    switch (direction) {
    case DesignElement::DAMPINGALPHA:
      alpha = 1.0;
      beta = 0.0;
      break;
    case DesignElement::DAMPINGBETA:
      alpha = 0.0;
      beta = 1.0;
      break;
    case DesignElement::NO_DERIVATIVE:
      alpha = GetParameter(DesignElement::DAMPINGALPHA);
      beta = GetParameter(DesignElement::DAMPINGBETA);
      break;
    default:
      alpha = 0.0;
      beta = 0.0;
      break;
    }
    return (true);
  } else {
    return (false);
  }
}

void DesignMaterial::DumpParams()
{
  for(const auto& iter : GetParameters())
    std::cout << "params[" << DesignElement::type.ToString(iter.first) << "] = " << iter.second << std::endl;
}

void DesignMaterial::SetEnums() {
  type.SetName("DesignMaterial::Type");
  type.Add(NO_TYPE, "no-type");
  type.Add(FMO, "fmo");
  type.Add(ORTHOTROPIC, "orthotropic");
  type.Add(DENSITY_TIMES_ORTHOTROPIC, "density-times-orthotropic");
  type.Add(ISOTROPIC, "isotropic");
  type.Add(LAME_ISOTROPIC, "lame-isotropic");
  type.Add(TRANSVERSAL_ISOTROPIC, "transversal-isotropic");
  type.Add(TRANSVERSAL_ISOTROPIC_BOXED, "transversal-isotropic-boxed");
  type.Add(DENSITY_TIMES_TRANSVERSAL_ISOTROPIC, "density-times-transversal-isotropic");
  type.Add(DENSITY_TIMES_TRANSVERSAL_ISOTROPIC_BOXED, "density-times-transversal-isotropic-boxed");
  type.Add(DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC, "density-times-rotated-transversal-isotropic");
  type.Add(DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED, "density-times-rotated-transversal-isotropic-boxed");
  type.Add(DENSITY_TIMES_ROT_PA12, "density-times-rotated-pa12");
  type.Add(DENSITY_TIMES_2D_TENSOR, "density-times-2dtensor");
  type.Add(DENSITY_TIMES_2D_TENSOR_CONSTANT_TRACE,
      "density-times-2dtensor-constant-trace");
  type.Add(DENSITY_TIMES_ROTATED_2D_TENSOR, "density-times-rotated-2dtensor");
  type.Add(LAMINATES, "laminates");
  type.Add(D_LAMINATES, "density-times-laminates");
  type.Add(HOM_RECT, "hom-rect");
  type.Add(D_HOM_RECT, "density-times-hom-rect");
  type.Add(HOM_RECT_C1, "hom-rect-C1");
  type.Add(HOM_ISO_C1, "hom-iso-C1");
  type.Add(MSFEM_C1, "msfem-C1");
  type.Add(D_INTERP_IN718_TENSOR, "density-times-interpolated-in718-tensor");
  type.Add(D_INTERP_IN718_TENSOR_ROT, "density-times-rotated-interpolated-in718-tensor");
  type.Add(SGP_MATLAB, "sgp-matlab");
  type.Add(SGP_GRADIENTCHECK, "sgp-gradient-check");
  type.Add(HEAT, "heat");
  type.Add(FEATURE_MAPPING_ANISO, "featureMappingAniso"); 

  transIsoType.SetName("DesignMaterial::TransIsoType");
  transIsoType.Add(TRANSISO_XY, "xy");
  transIsoType.Add(TRANSISO_YZ, "yz");
  transIsoType.Add(TRANSISO_XZ, "xz");

  rotationType.SetName("RotationType");
  rotationType.Add(ZXZ, "zxz");
  rotationType.Add(ZYZ, "zyz");
  rotationType.Add(YZY, "yzy");
  rotationType.Add(YXY, "yxy");
  rotationType.Add(XYX, "xyx");
  rotationType.Add(XZX, "xzx");
  rotationType.Add(XYZ, "xyz");
  rotationType.Add(YXZ, "yxz");
  rotationType.Add(XZY, "xzy");
  rotationType.Add(ZXY, "zxy");
  rotationType.Add(ZYX, "zyx");
  rotationType.Add(YZX, "yzx");
}

