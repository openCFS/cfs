#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>

#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "DesignMaterial.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "General/defs.hh"
#include "General/exception.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/TransferFunction.hh"
#include "PDE/SinglePDE.hh"
#include "Utils/StdVector.hh"
#include "Elements/2D/quad9fe.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/ParamHandling/Xerces.hh"
#include "MatVec/matrix.hh"

#ifdef USE_SGPP
#include "sgpp_base.hpp"
#include "sgpp_opt.hpp"
#endif

DECLARE_LOG(dm)
DEFINE_LOG(dm, "designMaterial")

using namespace CoupledField;
using std::string;

Enum<DesignMaterial::Type> DesignMaterial::type;
Enum<DesignMaterial::TransIsoType> DesignMaterial::transIsoType;
Enum<DesignMaterial::Notation> DesignMaterial::notation;

DesignMaterial::DesignMaterial(PtrParamNode pn, OptimizationMaterial::System material, StdVector<DesignID>& design, ErsatzMaterial* em)
{
  type_ = type.Parse(pn->Get("type")->As<string>());

  dim = domain->GetGrid()->GetDim();

  transIsoType_ = transIsoType.Parse(pn->Get("isoplane")->As<string>());

  massIsDesign_ = pn->Get("optimizeMass")->As<bool>();

  massFactor_ = pn->Get("massFactor")->As<double>();

  penalty_ = pn->Get("penalty")->As<double>();

  trace_ = pn->Get("trace")->As<double>();

  dampingIsDesign_ = pn->Get("optimizeDamping")->As<bool>();

  em_ = em;

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
    SetParameter(dt, params[i]->Get("value")->As<Double>());
    if (d.Find(dt) < 0) {
      d.Push_back(dt);
    }
    if (d.Find(DesignElement::SHEAR1) < 0 && type_ != HOM_RECT_C1) {
      SetParameter(DesignElement::SHEAR1, 0.5);
    }
  }
  if (!CheckRequiredDesigns(d)) {
    throw Exception("Not all Parameters for chosen DesignMaterial given");
  }else if(d.GetSize() > r){ // design.GetSize() < r is impossible as CheckRequiredDesigns passed
    info->Get("optimization/header/designSpace")->Get(ParamNode::WARNING)->SetValue("There are designs specified that are not used!");
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

  if((type_ == REDBAS_PARAM) || (type_ == REDBAS_FREE) || (type_ == REDBAS_MAPPING))
  {

    std::cout << "Reading files " << std::endl;

    std::string file = pn->Get("modRed/file")->As<std::string>();
    Xerces xerces(file);
    PtrParamNode root = xerces.CreateParamNodeInstance();
    PtrParamNode hr = root->Get("modRed");
    dimension_ = hr->Get("dimension")->As<UInt>();
    all_param_ = hr->Get("allParameters")->As<bool>();
    PtrParamNode mean = hr->Get("meantensor");
    PtrParamNode tensor = mean->Get("tensor");
    Matrix<double> mat(3,1);
    ParamTools::AsTensor<double>(tensor->Get("real"),3, 1, mat);
    mean_tensor_ = mat;

    StdVector<Matrix<double> > matrices(30);
    StdVector<Vector<double> > vectors(12);
    mod_red_matrices_ = matrices;
    mod_red_vectors_ = vectors;
//    PtrParamNode matnode = hr->Get("matrices");
//    PtrParamNode vecnode = hr->Get("vectors");

    StdVector<std::string> tensor_comp(3);// =["11", "12", "33"];
    tensor_comp[0]="11";
    tensor_comp[1]="12";
    tensor_comp[2]="33";

    for (int tensor_int =0; tensor_int <3; tensor_int ++)
    {

        FillModRedMatrices(hr, tensor_comp, tensor_int, dimension_);
        FillModRedVectors(hr, tensor_comp,tensor_int, dimension_);
    }

    std::cout << "Files read" << std::endl;
  }


  if((type_ == GREEDY_FREE) || (type_ == GREEDY_PARAM) || (type_ == GREEDY_MAPPING))
  {
    std::string file = pn->Get("modRed/file")->As<std::string>();
    Xerces xerces(file);
    PtrParamNode root = xerces.CreateParamNodeInstance();
    PtrParamNode hr = root->Get("modRed");

    //Here dimension is the number of terms in the tensor-product expansion
    dimension_ = hr->Get("dimension")->As<UInt>();
    dimension_tot_ = hr->Get("TotalDimension")->As<UInt>();
   // dimension_tot_ = hr->Get("dimensionTotale")->As<UInt>();
    Nl_ = hr->Get("Nl")->As<UInt>();
    Na_ = hr->Get("Na")->As<UInt>();
    lmin_ = hr->Get("lmin")->As<double>();
    lmax_ = hr->Get("lmax")->As<double>();

    all_param_ = hr->Get("allParameters")->As<bool>();

    //If type_ == GREEDY_FREE, MAPPING necessarily, we must have all_param_
    assert((type_==GREEDY_PARAM) || (all_param_));

    PtrParamNode mean = hr->Get("meantensor");
    PtrParamNode tensor = mean->Get("tensor");
    Matrix<double> mat(3,1);
    ParamTools::AsTensor<double>(tensor->Get("real"),3, 1, mat);
    mean_tensor_ = mat;

    int num_param=3;
    if (all_param_) num_param = 4;

    StdVector<Matrix<double> > mat_param(num_param);
    matrices_param_ = mat_param;

///////////////Initialisation of matl2////////////////////////

    /*mat.Resize(Nl_+1,3*dimension_);

    PtrParamNode matl2 = hr->Get("matl2");
    tensor = matl2->Get("tensor");
    Matrix<double> mattrans(Nl_+1,3*dimension_);
    ParamTools::AsTensor<double>(tensor->Get("real"),Nl_+1, 3*dimension_, mat);
    matrices_param_[num_param-1] = mat;*/


    Matrix<double> mattot(Nl_+1, 3*dimension_tot_);
    PtrParamNode matl2 = hr->Get("matl2");
    tensor = matl2->Get("tensor");
    ParamTools::AsTensor<double>(tensor->Get("real"),Nl_+1, 3*dimension_tot_, mattot);

    Matrix<double> mattrans(Nl_+1,dimension_);
    mat.Resize(Nl_+1, 3*dimension_);
    mat.Init();
    for (int i=0; i<3; i++)
    {
      mattot.GetSubMatrix(mattrans, 0, i*dimension_tot_);
      mat.SetSubMatrix(mattrans, 0, i*dimension_);
    }
    matrices_param_[num_param-1] = mat;

/////////Initialisation of matl1 /////////////////////////////

    /*mat.Resize(Nl_+1,3*dimension_);
    PtrParamNode matl1 = hr->Get("matl1");
    tensor = matl1->Get("tensor");
    ParamTools::AsTensor<double>(tensor->Get("real"),Nl_+1, 3*dimension_, mat);
    matrices_param_[num_param-2] = mat;*/

    mattot.Resize(Nl_+1, 3*dimension_tot_);
    mattot.Init();
    PtrParamNode matl1 = hr->Get("matl1");
    tensor = matl1->Get("tensor");
    ParamTools::AsTensor<double>(tensor->Get("real"),Nl_+1, 3*dimension_tot_, mattot);

    mattrans.Resize(Nl_+1,dimension_);
    mattrans.Init();
    mat.Resize(Nl_+1, 3*dimension_);
    mat.Init();
    for (int i=0; i<3; i++)
    {
      mattot.GetSubMatrix(mattrans, 0, i*dimension_tot_);
      mat.SetSubMatrix(mattrans, 0, i*dimension_);
    }
    matrices_param_[num_param-2] = mat;



///////////////Inintialisation of matphi//////////////////////////
 /* mat.Resize(2*Na_+1,3*dimension_);
  PtrParamNode matphi = hr->Get("matphi");
  tensor = matphi->Get("tensor");
  ParamTools::AsTensor<double>(tensor->Get("real"),2*Na_+1, 3*dimension_, mat);

  matrices_param_[num_param-3] = mat;*/

    mattot.Resize(2*Na_+1, 3*dimension_tot_);
    mattot.Init();
    PtrParamNode matphi = hr->Get("matphi");
    tensor = matphi->Get("tensor");
    ParamTools::AsTensor<double>(tensor->Get("real"),2*Na_+1, 3*dimension_tot_, mattot);

    mattrans.Resize(2*Na_+1,dimension_);
    mattrans.Init();
    mat.Resize(2*Na_+1, 3*dimension_);
    mat.Init();
    for (int i=0; i<3; i++)
    {
      mattot.GetSubMatrix(mattrans, 0, i*dimension_tot_);
      mat.SetSubMatrix(mattrans, 0, i*dimension_);
    }
    matrices_param_[num_param-3] = mat;



    if (all_param_)
    {

          /*PtrParamNode mattheta = hr->Get("mattheta");
          tensor = mattheta->Get("tensor");
          mat.Resize(2*Na_+1,3*dimension_);
          ParamTools::AsTensor<double>(tensor->Get("real"),2*Na_+1, 3*dimension_, mat);
          matrices_param_[num_param-4] = mat;*/

          mattot.Resize(2*Na_+1, 3*dimension_tot_);
          mattot.Init();
          PtrParamNode mattheta = hr->Get("mattheta");
          tensor = mattheta->Get("tensor");
          ParamTools::AsTensor<double>(tensor->Get("real"),2*Na_+1, 3*dimension_tot_, mattot);

          mattrans.Resize(2*Na_+1,dimension_);
          mattrans.Init();
          mat.Resize(2*Na_+1, 3*dimension_);
          mat.Init();
          for (int i=0; i<3; i++)
          {
            mattot.GetSubMatrix(mattrans, 0, i*dimension_tot_);
            mat.SetSubMatrix(mattrans, 0, i*dimension_);
          }
          matrices_param_[num_param-4] = mat;

    }


    StdVector<Matrix<double> > matrices(30);
    StdVector<Vector<double> > vectors(12);
    mod_red_matrices_ = matrices;
    mod_red_vectors_ = vectors;

        StdVector<std::string> tensor_comp(3);// =["11", "12", "33"];
           tensor_comp[0]="11";
           tensor_comp[1]="12";
           tensor_comp[2]="33";

           for (int tensor_int =0; tensor_int <3; tensor_int ++)
           {

               FillModRedMatrices(hr, tensor_comp, tensor_int, dimension_, dimension_tot_);
               FillModRedVectors(hr, tensor_comp,tensor_int, dimension_, dimension_tot_);
           }
    std::cout << "Matrices filled" << std::endl;
  }

  if (type_ == HOM_RECT_C1) {
    PtrParamNode hr = pn->Get("homRectC1");
    std::string file = hr->Get("file")->As<std::string>();
    // read interpolation method
    std::string interpolation_str = hr->Get("interpolation")->As<std::string>();

    // full C1 interpolation with text coefficients for 2D (without shearing angle)
    if (interpolation_str == "c1_text_2d") {
      interpolation_ = C1;
      // read coefficients from text file
      std::string filename = hr->Get("file")->As<std::string>();
      std::ifstream file(filename.c_str());
      std::string word;
      // read notation
      file >> word;
      Notation notation = ((word == "voigt") ? VOIGT : HILL_MANDEL);
      // read a
      unsigned int a_count;
      file >> a_count;
      hom_rect_a_ = Matrix<double>(a_count, 1);
      for (unsigned int i = 0; i < a_count; i++) {
        file >> hom_rect_a_[i][0];
      }
      // read b
      unsigned int b_count;
      file >> b_count;
      hom_rect_b_ = Matrix<double>(b_count, 1);
      for (unsigned int i = 0; i < b_count; i++) {
        file >> hom_rect_b_[i][0];
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
      hom_rect_coeff33_ = hom_rect_coeff33_ * (notation == VOIGT ? 2.0 : 1.0);
      file.close();

    // full C1 interpolation with text coefficients for 3D (with shearing angle)
    } else if (interpolation_str == "c1_text_3d") {
      interpolation_ = C1;
      // read coefficients from text file
      std::string filename = hr->Get("file")->As<std::string>();
      std::ifstream file(filename.c_str());
      std::string word;
      // read notation
      file >> word;
      Notation notation = ((word == "voigt") ? VOIGT : HILL_MANDEL);
      // read a
      unsigned int a_count;
      file >> a_count;
      hom_rect_a_ = Matrix<double>(a_count, 1);
      for (unsigned int i = 0; i < a_count; i++) {
        file >> hom_rect_a_[i][0];
      }
      // read b
      unsigned int b_count;
      file >> b_count;
      hom_rect_b_ = Matrix<double>(b_count, 1);
      for (unsigned int i = 0; i < b_count; i++) {
        file >> hom_rect_b_[i][0];
      }
      // read c
      unsigned int c_count;
      file >> c_count;
      hom_rect_c_ = Matrix<double>(c_count, 1);
      for (unsigned int i = 0; i < c_count; i++) {
        file >> hom_rect_c_[i][0];
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
      hom_rect_coeff33_ = hom_rect_coeff33_ * (notation == VOIGT ? 2.0 : 1.0);
      file.close();

    // full C1 interpolation with XML coefficients
    } else if (interpolation_str == "c1") {
      interpolation_ = C1;
      Xerces xerces(file);
      PtrParamNode root = xerces.CreateParamNodeInstance();
      Notation notation = root->Get("notation")->As<string>() == "voigt" ? VOIGT : HILL_MANDEL;
      if (dim == 2) {
        //Check for old format of coefficient file. Delete at some point in the future.
        if (!root->Has("a/matrix/real")) {
          ParamTools::AsMatrix<double>(root->Get("a/matrix"),hom_rect_a_);
          ParamTools::AsMatrix<double>(root->Get("b/matrix"),hom_rect_b_);
          ParamTools::AsMatrix<double>(root->Get("coeff11/matrix"), hom_rect_coeff11_);
          ParamTools::AsMatrix<double>(root->Get("coeff12/matrix"), hom_rect_coeff12_);
          ParamTools::AsMatrix<double>(root->Get("coeff22/matrix"), hom_rect_coeff22_);
          ParamTools::AsMatrix<double>(root->Get("coeff33/matrix"), hom_rect_coeff33_);
          hom_rect_coeff33_ = hom_rect_coeff33_ * (notation == VOIGT ? 2.0 : 1.0);
          if (root->Has("c")) {
            ParamTools::AsMatrix<double>(root->Get("c/matrix"), hom_rect_c_);
            ParamTools::AsMatrix<double>(root->Get("coeff13/matrix"), hom_rect_coeff13_);
            ParamTools::AsMatrix<double>(root->Get("coeff23/matrix"), hom_rect_coeff23_);
            hom_rect_coeff13_ = hom_rect_coeff13_ * (notation == VOIGT ? sqrt(2.0) : 1.0);
            hom_rect_coeff23_ = hom_rect_coeff23_ * (notation == VOIGT ? sqrt(2.0) : 1.0);
          } else {
            unsigned int dim1 = hom_rect_coeff11_.GetNumRows();
            unsigned int dim2 = hom_rect_coeff11_.GetNumCols();
            hom_rect_c_.Resize(0,0);
            hom_rect_c_.Init();
            hom_rect_coeff13_.Resize(dim1,dim2);
            hom_rect_coeff13_.Init();
            hom_rect_coeff23_.Resize(dim1,dim2);
            hom_rect_coeff23_.Init();
          }
        } else {
          int dim1 = root->Get("coeff11/matrix/dim1")->As<int>();
          int dim2 = root->Get("coeff11/matrix/dim2")->As<int>();
          int dim3 = root->Get("a/matrix/dim1")->As<int>();
          int dim4 = root->Get("b/matrix/dim1")->As<int>();
          ParamTools::AsTensor<double>(root->Get("coeff11/matrix/real"), dim1, dim2, hom_rect_coeff11_);
          ParamTools::AsTensor<double>(root->Get("coeff12/matrix/real"), dim1, dim2, hom_rect_coeff12_);
          ParamTools::AsTensor<double>(root->Get("coeff22/matrix/real"), dim1, dim2, hom_rect_coeff22_);
          ParamTools::AsTensor<double>(root->Get("coeff33/matrix/real"), dim1, dim2, hom_rect_coeff33_);
          ParamTools::AsTensor<double>(root->Get("a/matrix/real"), dim3, 1, hom_rect_a_);
          ParamTools::AsTensor<double>(root->Get("b/matrix/real"), dim4, 1, hom_rect_b_);
          hom_rect_coeff33_ = hom_rect_coeff33_ * (notation == VOIGT ? 2.0 : 1.0);
          if (root->Has("c")) {
            int dim5 = root->Get("c/matrix/dim1")->As<int>();
            ParamTools::AsTensor<double>(root->Get("c/matrix/real"), dim5, 1, hom_rect_c_);
            ParamTools::AsTensor<double>(root->Get("coeff13/matrix/real"), dim1, dim2, hom_rect_coeff13_);
            ParamTools::AsTensor<double>(root->Get("coeff23/matrix/real"), dim1, dim2, hom_rect_coeff23_);
            hom_rect_coeff13_ = hom_rect_coeff13_ * (notation == VOIGT ? sqrt(2.0) : 1.0);
            hom_rect_coeff23_ = hom_rect_coeff23_ * (notation == VOIGT ? sqrt(2.0) : 1.0);
          } else {
            hom_rect_c_.Resize(0,0);
            hom_rect_c_.Init();
            hom_rect_coeff13_.Resize(dim1,dim2);
            hom_rect_coeff13_.Init();
            hom_rect_coeff23_.Resize(dim1,dim2);
            hom_rect_coeff23_.Init();
          }
        }
      } else if (dim == 3) {
        int dim1 = root->Get("coeff11/matrix/dim1")->As<int>();
        int dim2 = root->Get("coeff11/matrix/dim2")->As<int>();
        int dim3 = root->Get("a/matrix/dim1")->As<int>();
        int dim4 = root->Get("b/matrix/dim1")->As<int>();
        int dim5 = root->Get("c/matrix/dim1")->As<int>();
        ParamTools::AsTensor<double>(root->Get("coeff11/matrix/real"), dim1, dim2, hom_rect_coeff11_);
        ParamTools::AsTensor<double>(root->Get("coeff12/matrix/real"), dim1, dim2, hom_rect_coeff12_);
        ParamTools::AsTensor<double>(root->Get("coeff22/matrix/real"), dim1, dim2, hom_rect_coeff22_);
        ParamTools::AsTensor<double>(root->Get("coeff33/matrix/real"), dim1, dim2, hom_rect_coeff33_);
        ParamTools::AsTensor<double>(root->Get("coeff13/matrix/real"), dim1, dim2, hom_rect_coeff13_);
        ParamTools::AsTensor<double>(root->Get("coeff23/matrix/real"), dim1, dim2, hom_rect_coeff23_);
        ParamTools::AsTensor<double>(root->Get("coeff44/matrix/real"), dim1, dim2, hom_rect_coeff44_);
        ParamTools::AsTensor<double>(root->Get("coeff55/matrix/real"), dim1, dim2, hom_rect_coeff55_);
        ParamTools::AsTensor<double>(root->Get("coeff66/matrix/real"), dim1, dim2, hom_rect_coeff66_);
        ParamTools::AsTensor<double>(root->Get("a/matrix/real"), dim3, 1, hom_rect_a_);
        ParamTools::AsTensor<double>(root->Get("b/matrix/real"), dim4, 1, hom_rect_b_);
        ParamTools::AsTensor<double>(root->Get("c/matrix/real"), dim5, 1, hom_rect_c_);
        // the internal tensor representation in 3D hom_rect_samples_ is VOIGT
        hom_rect_coeff44_ = hom_rect_coeff44_ * (notation == VOIGT ? 1.0 : 0.5);
        hom_rect_coeff55_ = hom_rect_coeff55_ * (notation == VOIGT ? 1.0 : 0.5);
        hom_rect_coeff66_ = hom_rect_coeff66_ * (notation == VOIGT ? 1.0 : 0.5);
        // the tensor is orthotropic
      }
    } else {

#ifdef USE_SGPP
      // sparse grid interpolation
      if (interpolation_str == "sgpp") {
        interpolation_ = SGPP;
        unsigned int dimension = (shearIsDesign_ ? 3 : 2);
        std::string basis_str = hr->Get("sgppBasis")->As<std::string>();
        // sparse grid basis to be used
        if ((basis_str == "bspline") || (basis_str == "modbspline")) {
          // B-spline basis
          bspline_degree_ = hr->Get("bsplineDegree")->As<unsigned int>();
          if (basis_str == "bspline") {
            // B-splines
            sgpp_basis_ = BSPLINE;
            grid_ = sg::base::Grid::createBsplineGrid(dimension, bspline_degree_);
          } else {
            // modified B-splines
            sgpp_basis_ = MODBSPLINE;
            grid_ = sg::base::Grid::createModBsplineGrid(dimension, bspline_degree_);
          }
          // optional verbosity for initial hierarchisation
          sg::opt::tools::printer.setVerbosity(2);
        } else if (basis_str == "modlinear") {
          // modified linear basis
          sgpp_basis_ = MODLINEAR;
          grid_ = sg::base::Grid::createModLinearGrid(dimension);
        } else {
          // (un)modified linear basis
          sgpp_basis_ = LINEAR;
          grid_ = sg::base::Grid::createLinearGrid(dimension);
        }
        // read coefficients
        InitializeSparseGrid(file.c_str());
        LOG_DBG(dm)<<"Sparse grid initialized and filled.";

      // full B-spline interpolation with XML coefficients
      } else if (interpolation_str == "full_bspline") {
        interpolation_ = FULL_BSPLINE;
        bspline_degree_ = hr->Get("bsplineDegree")->As<unsigned int>();
        Xerces xerces(file);
        PtrParamNode root = xerces.CreateParamNodeInstance();
        Notation notation = ((root->Get("notation")->As<string>() == "voigt") ? VOIGT : HILL_MANDEL);
        // read coefficients from XML
        ParamTools::AsMatrix<double>(root->Get("coeff11/matrix"), full_bspline_coeff11_);
        ParamTools::AsMatrix<double>(root->Get("coeff12/matrix"), full_bspline_coeff12_);
        if (shearIsDesign_) {
          ParamTools::AsMatrix<double>(root->Get("coeff13/matrix"), full_bspline_coeff13_);
        }
        ParamTools::AsMatrix<double>(root->Get("coeff22/matrix"), full_bspline_coeff22_);
        if (shearIsDesign_) {
          ParamTools::AsMatrix<double>(root->Get("coeff23/matrix"), full_bspline_coeff23_);
        }
        ParamTools::AsMatrix<double>(root->Get("coeff33/matrix"), full_bspline_coeff33_);
        // apply notation
        full_bspline_coeff33_ = full_bspline_coeff33_ * ((notation == VOIGT) ? 2.0 : 1.0);

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
        Notation notation = ((word == "voigt") ? VOIGT : HILL_MANDEL);
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
        full_bspline_coeff33_ = full_bspline_coeff33_ * (notation == VOIGT ? 2.0 : 1.0);
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
        Notation notation = ((word == "voigt") ? VOIGT : HILL_MANDEL);
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
        full_bspline_coeff33_ = full_bspline_coeff33_ * (notation == VOIGT ? 2.0 : 1.0);
        file.close();
      }
#else //USE_SGPP
    EXCEPTION("CFS is compiled without SGpp toolbox! Please recompile with SGpp or choose another interpolation method.")
#endif //USE_SGPP
    }
  }

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

void DesignMaterial::FillHomRectSamples(PtrParamNode homRect, unsigned int idx,
    const string& a, const string& b) {
  // the internal tensor representation in hom_rect_samples_ is HILL-MANDEL!
  Notation notation =
    homRect->Get("notation")->As<string>() == "voigt" ? VOIGT : HILL_MANDEL;

  PtrParamNode data = homRect->GetByVal("data", "a", a, "b", b);
  hom_rect_samples_[idx][DesignElement::TENSOR11 - DesignElement::TENSOR11] = data->Get("e11")->As<double>();
  hom_rect_samples_[idx][DesignElement::TENSOR12 - DesignElement::TENSOR11] = data->Get("e12")->As<double>();
  hom_rect_samples_[idx][DesignElement::TENSOR22 - DesignElement::TENSOR11] = data->Get("e22")->As<double>();
  hom_rect_samples_[idx][DesignElement::TENSOR33 - DesignElement::TENSOR11] = data->Get("e33")->As<double>() * (notation == VOIGT ? 2.0 : 1.0);
  hom_rect_samples_[idx][DesignElement::TENSOR13 - DesignElement::TENSOR11] = 0.0;
  hom_rect_samples_[idx][DesignElement::TENSOR23 - DesignElement::TENSOR11] = 0.0;
}

void DesignMaterial::FillModRedMatrices(PtrParamNode matnode, const StdVector<std::string>& tensor_comp, const int& tensor_int, const UInt& dimbas)
{

      PtrParamNode mat1 = matnode->Get("matuxux" + tensor_comp[tensor_int]);
      PtrParamNode tensor = mat1->Get("tensor");
      Matrix<double> mat(dimbas,dimbas);
      ParamTools::AsTensor<double>(tensor->Get("real"),dimbas, dimbas, mat);
      mod_red_matrices_[10*tensor_int +0] = mat;

      PtrParamNode mat2= matnode->Get("matuxvx" + tensor_comp[tensor_int] );
      tensor = mat2->Get("tensor");
      ParamTools::AsTensor<double>(tensor->Get("real"),dimbas, dimbas, mat);
       mod_red_matrices_[10*tensor_int + 1] = mat;

      PtrParamNode mat3= matnode->Get("matvxvx" + tensor_comp[tensor_int]);
      tensor = mat3->Get("tensor");
      ParamTools::AsTensor<double>(tensor->Get("real"),dimbas, dimbas, mat);
      mod_red_matrices_[10*tensor_int + 2] = mat;


      PtrParamNode mat4 = matnode->Get("matuxuy" + tensor_comp[tensor_int]);
      tensor = mat4->Get("tensor");
      ParamTools::AsTensor<double>(tensor->Get("real"),dimbas, dimbas, mat);
      mod_red_matrices_[ 10*tensor_int +3] = mat;

      PtrParamNode mat5 = matnode->Get("matuxvy" + tensor_comp[tensor_int]);
      tensor = mat5->Get("tensor");
      ParamTools::AsTensor<double>(tensor->Get("real"),dimbas, dimbas, mat);
       mod_red_matrices_[10*tensor_int + 4] = mat;

       PtrParamNode mat6 = matnode->Get("matuyvx" + tensor_comp[tensor_int]);
       tensor = mat6->Get("tensor");
       ParamTools::AsTensor<double>(tensor->Get("real"),dimbas, dimbas, mat);
       mod_red_matrices_[ 10*tensor_int + 5] = mat;

      PtrParamNode mat7 = matnode->Get("matvxvy" + tensor_comp[tensor_int]);
      tensor = mat7->Get("tensor");
      ParamTools::AsTensor<double>(tensor->Get("real"),dimbas, dimbas, mat);
      mod_red_matrices_[10*tensor_int + 6] = mat;


      PtrParamNode mat8 = matnode->Get("matuyuy" + tensor_comp[tensor_int]);
      tensor = mat8->Get("tensor");
      ParamTools::AsTensor<double>(tensor->Get("real"),dimbas, dimbas, mat);
       mod_red_matrices_[10*tensor_int + 7] = mat;

       PtrParamNode mat9 = matnode->Get("matuyvy" + tensor_comp[tensor_int]);
       tensor = mat9->Get("tensor");
       ParamTools::AsTensor<double>(tensor->Get("real"),dimbas, dimbas, mat);
       mod_red_matrices_[10*tensor_int + 8] = mat;

      PtrParamNode mat10 = matnode->Get("matvyvy" + tensor_comp[tensor_int]);
      tensor = mat10->Get("tensor");
      ParamTools::AsTensor<double>(tensor->Get("real"),dimbas, dimbas, mat);
      mod_red_matrices_[10*tensor_int + 9] = mat;

}

void DesignMaterial::FillModRedVectors(PtrParamNode vecnode, const StdVector<std::string>& tensor_comp, const int& tensor_int, const UInt& dimbas)
{

      PtrParamNode vecux = vecnode->Get("vecux" + tensor_comp[tensor_int] );
      PtrParamNode tensor = vecux->Get("tensor");
      Matrix<double> mat(dimbas,1);
      Vector<double> vec(dimbas);
      ParamTools::AsTensor<double>(tensor->Get("real"),dimbas, 1, mat);
      mat.GetCol(vec,0);
      mod_red_vectors_[4*tensor_int +0] = vec;

      PtrParamNode vecuy = vecnode->Get("vecuy" + tensor_comp[tensor_int]);
      tensor = vecuy->Get("tensor");
      ParamTools::AsTensor<double>(tensor->Get("real"),dimbas, 1, mat);
      mat.GetCol(vec,0);
      mod_red_vectors_[4*tensor_int +1] = vec;

      PtrParamNode vecvx = vecnode->Get("vecvx" + tensor_comp[tensor_int]);
      tensor = vecvx->Get("tensor");
      ParamTools::AsTensor<double>(tensor->Get("real"),dimbas, 1, mat);
      mat.GetCol(vec,0);
      mod_red_vectors_[ 4*tensor_int +2] = vec;

      PtrParamNode vecvy = vecnode->Get("vecvy" + tensor_comp[tensor_int]);
      tensor = vecvy->Get("tensor");
      ParamTools::AsTensor<double>(tensor->Get("real"),dimbas, 1, mat);
      mat.GetCol(vec,0);
      mod_red_vectors_[ 4*tensor_int +3] = vec;

}


void DesignMaterial::FillModRedMatrices(PtrParamNode matnode, const StdVector<std::string>& tensor_comp, const int& tensor_int, const UInt& dimbas, const UInt& dimbastot)
{

      PtrParamNode mat1 = matnode->Get("matuxux" + tensor_comp[tensor_int]);
      PtrParamNode tensor = mat1->Get("tensor");
      Matrix<double> mattot(3*dimbastot,3*dimbastot);
      ParamTools::AsTensor<double>(tensor->Get("real"),3*dimbastot, 3*dimbastot, mattot);
      Matrix<double> mat(3*dimbas, 3*dimbas);
      Matrix<double> mattrans(dimbas, dimbas);
      for (int i=0; i< 3; i++)
      {
        for (int j=0; j<3; j++)
        {
          mattot.GetSubMatrix(mattrans, i*dimbastot, j*dimbastot);
          mat.SetSubMatrix(mattrans, i*dimbas, j*dimbas);
        }
      }
      mod_red_matrices_[10*tensor_int +0] = mat;


      PtrParamNode mat2= matnode->Get("matuxvx" + tensor_comp[tensor_int] );
      tensor = mat2->Get("tensor");
      ParamTools::AsTensor<double>(tensor->Get("real"),3*dimbastot, 3*dimbastot, mattot);
      for (int i=0; i< 3; i++)
      {
        for (int j=0; j<3; j++)
        {
          mattot.GetSubMatrix(mattrans, i*dimbastot, j*dimbastot);
          mat.SetSubMatrix(mattrans, i*dimbas, j*dimbas);
        }
      }
       mod_red_matrices_[10*tensor_int + 1] = mat;


      PtrParamNode mat3= matnode->Get("matvxvx" + tensor_comp[tensor_int]);
      tensor = mat3->Get("tensor");
      ParamTools::AsTensor<double>(tensor->Get("real"),3*dimbastot, 3*dimbastot, mattot);
      for (int i=0; i< 3; i++)
      {
        for (int j=0; j<3; j++)
        {
          mattot.GetSubMatrix(mattrans, i*dimbastot, j*dimbastot);
          mat.SetSubMatrix(mattrans, i*dimbas, j*dimbas);
        }
      }
      mod_red_matrices_[10*tensor_int + 2] = mat;


      PtrParamNode mat4 = matnode->Get("matuxuy" + tensor_comp[tensor_int]);
      tensor = mat4->Get("tensor");
      ParamTools::AsTensor<double>(tensor->Get("real"),3*dimbastot, 3*dimbastot, mattot);
      for (int i=0; i< 3; i++)
      {
        for (int j=0; j<3; j++)
        {
          mattot.GetSubMatrix(mattrans, i*dimbastot, j*dimbastot);
          mat.SetSubMatrix(mattrans, i*dimbas, j*dimbas);
        }
      }
      mod_red_matrices_[ 10*tensor_int +3] = mat;

      PtrParamNode mat5 = matnode->Get("matuxvy" + tensor_comp[tensor_int]);
      tensor = mat5->Get("tensor");
      ParamTools::AsTensor<double>(tensor->Get("real"),3*dimbastot, 3*dimbastot, mattot);
      for (int i=0; i< 3; i++)
      {
        for (int j=0; j<3; j++)
        {
          mattot.GetSubMatrix(mattrans, i*dimbastot, j*dimbastot);
          mat.SetSubMatrix(mattrans, i*dimbas, j*dimbas);
        }
      }
       mod_red_matrices_[10*tensor_int + 4] = mat;

       PtrParamNode mat6 = matnode->Get("matuyvx" + tensor_comp[tensor_int]);
       tensor = mat6->Get("tensor");
       ParamTools::AsTensor<double>(tensor->Get("real"),3*dimbastot, 3*dimbastot, mattot);
       for (int i=0; i< 3; i++)
       {
         for (int j=0; j<3; j++)
         {
           mattot.GetSubMatrix(mattrans, i*dimbastot, j*dimbastot);
           mat.SetSubMatrix(mattrans, i*dimbas, j*dimbas);
         }
       }
       mod_red_matrices_[ 10*tensor_int + 5] = mat;

      PtrParamNode mat7 = matnode->Get("matvxvy" + tensor_comp[tensor_int]);
      tensor = mat7->Get("tensor");
      ParamTools::AsTensor<double>(tensor->Get("real"),3*dimbastot, 3*dimbastot, mattot);
      for (int i=0; i< 3; i++)
             {
               for (int j=0; j<3; j++)
               {
                 mattot.GetSubMatrix(mattrans, i*dimbastot, j*dimbastot);
                 mat.SetSubMatrix(mattrans, i*dimbas, j*dimbas);
               }
             }
      mod_red_matrices_[10*tensor_int + 6] = mat;

      PtrParamNode mat8 = matnode->Get("matuyuy" + tensor_comp[tensor_int]);
      tensor = mat8->Get("tensor");
      ParamTools::AsTensor<double>(tensor->Get("real"),3*dimbastot, 3*dimbastot, mattot);
      for (int i=0; i< 3; i++)
      {
        for (int j=0; j<3; j++)
        {
          mattot.GetSubMatrix(mattrans, i*dimbastot, j*dimbastot);
          mat.SetSubMatrix(mattrans, i*dimbas, j*dimbas);
        }
      }
       mod_red_matrices_[10*tensor_int + 7] = mat;

       PtrParamNode mat9 = matnode->Get("matuyvy" + tensor_comp[tensor_int]);
       tensor = mat9->Get("tensor");
       ParamTools::AsTensor<double>(tensor->Get("real"),3*dimbastot, 3*dimbastot, mattot);
       for (int i=0; i< 3; i++)
       {
         for (int j=0; j<3; j++)
         {
           mattot.GetSubMatrix(mattrans, i*dimbastot, j*dimbastot);
           mat.SetSubMatrix(mattrans, i*dimbas, j*dimbas);
         }
       }
       mod_red_matrices_[10*tensor_int + 8] = mat;

      PtrParamNode mat10 = matnode->Get("matvyvy" + tensor_comp[tensor_int]);
      tensor = mat10->Get("tensor");
      ParamTools::AsTensor<double>(tensor->Get("real"),3*dimbastot, 3*dimbastot, mattot);
      for (int i=0; i< 3; i++)
      {
        for (int j=0; j<3; j++)
        {
          mattot.GetSubMatrix(mattrans, i*dimbastot, j*dimbastot);
          mat.SetSubMatrix(mattrans, i*dimbas, j*dimbas);
        }
      }
      mod_red_matrices_[10*tensor_int + 9] = mat;

}

void DesignMaterial::FillModRedVectors(PtrParamNode vecnode, const StdVector<std::string>& tensor_comp, const int& tensor_int, const UInt& dimbas, const UInt& dimbastot)
{

      PtrParamNode vecux = vecnode->Get("vecux" + tensor_comp[tensor_int] );
      PtrParamNode tensor = vecux->Get("tensor");
      Matrix<double> mattot(3*dimbastot,1);
      ParamTools::AsTensor<double>(tensor->Get("real"),3*dimbastot, 1, mattot);
      Matrix<double> mat(3*dimbas,1);
      Matrix<double> mattrans(dimbas,1);
      for (int i=0; i< 3; i++)
      {
        mattot.GetSubMatrix(mattrans, i*dimbastot, 0);
        mat.SetSubMatrix(mattrans, i*dimbas, 0);
      }
      Vector<double> vec(dimbas);
      mat.GetCol(vec,0);
      mod_red_vectors_[4*tensor_int +0] = vec;

      PtrParamNode vecuy = vecnode->Get("vecuy" + tensor_comp[tensor_int]);
      tensor = vecuy->Get("tensor");
      ParamTools::AsTensor<double>(tensor->Get("real"),3*dimbastot, 1, mattot);
      for (int i=0; i< 3; i++)
            {
              mattot.GetSubMatrix(mattrans, i*dimbastot, 0);
              mat.SetSubMatrix(mattrans, i*dimbas, 0);
            }
      mat.GetCol(vec,0);
      mod_red_vectors_[4*tensor_int +1] = vec;

      PtrParamNode vecvx = vecnode->Get("vecvx" + tensor_comp[tensor_int]);
      tensor = vecvx->Get("tensor");
      ParamTools::AsTensor<double>(tensor->Get("real"),3*dimbastot, 1, mattot);
      for (int i=0; i< 3; i++)
            {
              mattot.GetSubMatrix(mattrans, i*dimbastot, 0);
              mat.SetSubMatrix(mattrans, i*dimbas, 0);
            }
      mat.GetCol(vec,0);
      mod_red_vectors_[ 4*tensor_int +2] = vec;

      PtrParamNode vecvy = vecnode->Get("vecvy" + tensor_comp[tensor_int]);
      tensor = vecvy->Get("tensor");
      ParamTools::AsTensor<double>(tensor->Get("real"),3*dimbastot, 1, mattot);
      for (int i=0; i< 3; i++)
            {
              mattot.GetSubMatrix(mattrans, i*dimbastot, 0);
              mat.SetSubMatrix(mattrans, i*dimbas, 0);
            }
      mat.GetCol(vec,0);
      mod_red_vectors_[ 4*tensor_int +3] = vec;

}

unsigned int DesignMaterial::RequiredParameters(
    OptimizationMaterial::System material) {
  unsigned int r = MassIsDesign() ? 1 : 0;
  if (DampingIsDesign()) {
    r += 2;
  }
  switch (type_) {
  case FMO:
    assert(
        material == OptimizationMaterial::MECH
            || material == OptimizationMaterial::PIEZOCOUPLING);
    return r + (material == OptimizationMaterial::MECH ? 6 : 15);
  case ISOTROPIC:
  case LAME_ISOTROPIC:
    return r + 2;
  case TRANSVERSAL_ISOTROPIC:
  case TRANSVERSAL_ISOTROPIC_BOXED:
    return domain->GetSinglePDE("mechanic")->GetSubTensorType() == PLANE_STRESS ? r + 4 : r + 5;
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
  case D_INTERP_TENSOR:
    return r + 2;
  case D_INTERP_TENSOR_ROT:
    return r + 3 + (dim == 3 ? 3 : 1);
  case REDBAS_PARAM:
  case REDBAS_FREE:
  case GREEDY_PARAM:
  case GREEDY_FREE:
    return r + 4;
  case GREEDY_MAPPING:
  case REDBAS_MAPPING:
    return r + 2;
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
    return (design.Find(DesignElement::TENSOR11) >= 0
        && design.Find(DesignElement::TENSOR22) >= 0
        && design.Find(DesignElement::TENSOR33) >= 0
        && design.Find(DesignElement::TENSOR23) >= 0
        && design.Find(DesignElement::TENSOR13) >= 0
        && design.Find(DesignElement::TENSOR12) >= 0);
  case ISOTROPIC:
    return (design.Find(DesignElement::EMODUL) >= 0
        && design.Find(DesignElement::POISSON) >= 0);
  case LAME_ISOTROPIC:
    return (design.Find(DesignElement::LAMELAMBDA) >= 0
        && design.Find(DesignElement::LAMEMU) >= 0);
  case TRANSVERSAL_ISOTROPIC:
  case TRANSVERSAL_ISOTROPIC_BOXED:
    return (
        design.Find(DesignElement::EMODULISO) >= 0
            && design.Find(DesignElement::EMODUL) >= 0
            && design.Find(DesignElement::POISSON) >= 0
            && design.Find(DesignElement::GMODUL) >= 0
            && domain->GetSinglePDE("mechanic")->GetSubTensorType() == PLANE_STRESS ? true : design.Find(DesignElement::POISSONISO) >= 0);
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
        && (dim != 2 || design.Find(DesignElement::ROTANGLE))
        && (dim != 3 || (design.Find(DesignElement::POISSONISO) >= 0 && design.Find(DesignElement::ROTANGLEX) >= 0 && design.Find(DesignElement::ROTANGLEY) >= 0 ) ) );
  case DENSITY_TIMES_ROT_PA12:
    return (design.Find(DesignElement::DENSITY) >= 0
        && design.Find(DesignElement::EMODULISO) >= 0
        && design.Find(DesignElement::POISSON) >= 0
        && design.Find(DesignElement::GMODUL) >= 0
        && (dim != 2 || design.Find(DesignElement::ROTANGLE) >= 0)
        && (dim != 3 || (design.Find(DesignElement::ROTANGLEX) >= 0 && design.Find(DesignElement::ROTANGLEY) >= 0 ) ) );
  case D_INTERP_TENSOR:
    return (design.Find(DesignElement::DENSITY) >= 0
        && design.Find(DesignElement::INTERPOLATION) >=0);
  case D_INTERP_TENSOR_ROT:
    return (design.Find(DesignElement::DENSITY) >= 0
        && design.Find(DesignElement::EMODUL) >=0
        && design.Find(DesignElement::INTERPOLATION) >=0
        && (dim != 2 || design.Find(DesignElement::ROTANGLE) >= 0)
        && (dim != 3 || (design.Find(DesignElement::ROTANGLEX) >= 0 && design.Find(DesignElement::ROTANGLEY) >= 0 && design.Find(DesignElement::ROTANGLEZ) >= 0 ) ) );
  case ORTHOTROPIC:
    return(design.Find(DesignElement::TENSOR11) >= 0
        && design.Find(DesignElement::TENSOR22) >= 0
        && design.Find(DesignElement::TENSOR33) >= 0
        && design.Find(DesignElement::TENSOR12) >= 0
        && design.Find(DesignElement::ROTANGLE) >= 0
        && design.Find(DesignElement::LOWER_EIG_BOUND) >= 0);
  case DENSITY_TIMES_ORTHOTROPIC:
    return(design.Find(DesignElement::DENSITY) >= 0
        && design.Find(DesignElement::TENSOR11) >= 0
        && design.Find(DesignElement::TENSOR22) >= 0
        && design.Find(DesignElement::TENSOR33) >= 0
        && design.Find(DesignElement::TENSOR12) >= 0
        && design.Find(DesignElement::ROTANGLE) >= 0
        && design.Find(DesignElement::LOWER_EIG_BOUND) >= 0);
  case DENSITY_TIMES_2D_TENSOR:
    return (design.Find(DesignElement::DENSITY) >= 0
        && design.Find(DesignElement::TENSOR11) >= 0
        && design.Find(DesignElement::TENSOR22) >= 0
        && design.Find(DesignElement::TENSOR33) >= 0
        && design.Find(DesignElement::TENSOR23) >= 0
        && design.Find(DesignElement::TENSOR13) >= 0
        && design.Find(DesignElement::TENSOR12) >= 0);
  case DENSITY_TIMES_2D_TENSOR_CONSTANT_TRACE:
    return (design.Find(DesignElement::DENSITY) >= 0
        && design.Find(DesignElement::TENSOR11) >= 0
        && design.Find(DesignElement::TENSOR22) >= 0
        && design.Find(DesignElement::TENSOR23) >= 0
        && design.Find(DesignElement::TENSOR13) >= 0
        && design.Find(DesignElement::TENSOR12) >= 0);
  case DENSITY_TIMES_ROTATED_2D_TENSOR:
    return (design.Find(DesignElement::DENSITY) >= 0
        && design.Find(DesignElement::TENSOR11) >= 0
        && design.Find(DesignElement::TENSOR33) >= 0
        && design.Find(DesignElement::TENSOR23) >= 0
        && design.Find(DesignElement::TENSOR13) >= 0
        && design.Find(DesignElement::TENSOR12) >= 0
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
  case REDBAS_PARAM:
  case GREEDY_PARAM:
    return(design.Find(DesignElement::ROTANGLE) >= 0
        && design.Find(DesignElement::ROTANGLE2) >= 0
        && design.Find(DesignElement::SCALING1) >= 0
        && design.Find(DesignElement::SCALING2) >= 0);
  case REDBAS_FREE:
  case GREEDY_FREE:
      return(design.Find(DesignElement::G11) >= 0
          && design.Find(DesignElement::G12) >= 0
          && design.Find(DesignElement::G21) >= 0
          && design.Find(DesignElement::G22) >= 0);
  case GREEDY_MAPPING:
  case REDBAS_MAPPING:
      return(design.Find(DesignElement::G_MAP_X) >= 0
          && design.Find(DesignElement::G_MAP_Y) >= 0
          );
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
          && design.Find(DesignElement::ROTANGLEX) >= 0
          && design.Find(DesignElement::ROTANGLEY) >= 0);
    } else {
      return (design.Find(DesignElement::STIFF1) >= 0
          && design.Find(DesignElement::STIFF2) >= 0
          && design.Find(DesignElement::SHEAR1) >= 0
          && design.Find(DesignElement::ROTANGLE) >= 0);
    }
  }
  assert(false);
  return false;
}

void DesignMaterial::SetParameter(const DesignElement::Type p,
    const double value) {
  params_[p] = value;
  // LOG_DBG3(dm) << "SP p=" << DesignElement::type.ToString(p) << " v=" << value;
}

void DesignMaterial::GetIsoMaterialTensor(Matrix<double>& t,
    SubTensorType subTensor, DesignElement::Type direction) {
  double E = params_[DesignElement::EMODUL];
  double nu = params_[DesignElement::POISSON];
  switch (direction) {
  case DesignElement::NO_DERIVATIVE: {
    double lambda = nu * E / ((1.0 + nu) * (1.0 - 2.0 * nu));
    double mu = E / (2.0 * (1.0 + nu));
    double diag = lambda + 2.0 * mu;
    SetIsoTensor(t, subTensor, diag, lambda, mu);
    break;
  }
  case DesignElement::EMODUL: {
    double dlambda_dE = nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    double dmu_dE = 1.0 / (2.0 * (1.0 + nu));
    double ddiag_dE = dlambda_dE + 2.0 * dmu_dE;
    SetIsoTensor(t, subTensor, ddiag_dE, dlambda_dE, dmu_dE);
    break;
  }
  case DesignElement::POISSON: {
    double dlambda_dnu = (1.0 + 2.0 * nu * nu) * E
        / ((1.0 + nu) * (1.0 + nu) * (1.0 - 2.0 * nu) * (1.0 - 2.0 * nu));
    double dmu_dnu = E / (-2.0 * (1.0 + nu) * (1.0 + nu));
    double ddiag_dnu = dlambda_dnu + 2.0 * dmu_dnu;
    SetIsoTensor(t, subTensor, ddiag_dnu, dlambda_dnu, dmu_dnu);
    break;
  }
  default:
    ZeroTensor(t, subTensor); // any derivative in any direction other than EMODUL or POISSON is zero
    break;
  }
}

double DesignMaterial::GetIsoMaterialMass(DesignElement::Type direction) {
  double E = params_[DesignElement::EMODUL];
  double nu = params_[DesignElement::POISSON];
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

void DesignMaterial::GetLameMaterialTensor(Matrix<double>& t,
    SubTensorType subTensor, DesignElement::Type direction) {
  switch (direction) {
  case DesignElement::NO_DERIVATIVE: {
    double lambda = params_[DesignElement::LAMELAMBDA];
    double mu = params_[DesignElement::LAMEMU];
    double diag = lambda + 2.0 * mu;
    SetIsoTensor(t, subTensor, diag, lambda, mu);
    break;
  }
  case DesignElement::LAMELAMBDA:
    SetIsoTensor(t, subTensor, 1.0, 1.0, 0.0);
    break;
  case DesignElement::LAMEMU:
    SetIsoTensor(t, subTensor, 2.0, 0.0, 1.0);
    break;
  default:
    ZeroTensor(t, subTensor); // any derivative in any direction other than EMODUL or POISSON is zero
    break;
  }
}

double DesignMaterial::GetLameMaterialMass(DesignElement::Type direction) {
  switch (direction) {
  case DesignElement::NO_DERIVATIVE: {
    double lambda = params_[DesignElement::LAMELAMBDA];
    double mu = params_[DesignElement::LAMEMU];
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

void DesignMaterial::GetTransIsoMaterialTensor(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction, Notation notation){
  LOG_DBG2(dm) << "GetTransIsoMaterialTensor called with direction=" << (direction == DesignElement::NO_DERIVATIVE ? "no_derivative" : DesignElement::type.ToString(direction)) << " and notation=" << notation;
  assert(type_ != DENSITY_TIMES_ROT_PA12 || subTensor == FULL);
  double E3(0.0);
  double E = params_[DesignElement::EMODULISO];
  if (type_ != DENSITY_TIMES_ROT_PA12)
    E3 = params_[DesignElement::EMODUL];
  else
  {
    E3 = 137.4 + 2.4*E;
    E = 145.0 - 5.8*E;
  }
  double G3 = params_[DesignElement::GMODUL];
  double nu13 = params_[DesignElement::POISSON]; //used as theta in the boxed formulations

  if (subTensor == PLANE_STRESS) {
    double dens(1.0), ninv2(0.0), D(0.0), D3(0.0), nD3(0.0);
    if((type_ == DENSITY_TIMES_TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_TRANSVERSAL_ISOTROPIC_BOXED
        || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED) && notation != HILL_MANDEL_NO_DENSITY)
    {
      dens = params_[DesignElement::DENSITY];
      dens = std::pow(dens, penalty_); // SIMP
 //     dens = dens/(1.0 + 10.0*(1.0-dens)); // RAMP
    }
    if (type_ == TRANSVERSAL_ISOTROPIC
        || type_ == DENSITY_TIMES_TRANSVERSAL_ISOTROPIC) {
      ninv2 = E3 - nu13 * nu13 * E;
      //assert(ninv2>0.0); //positivity of the elasticity tensor is violated. Use constraint "parametrized-plane-stress-pos-def" > 0
      ninv2 = 1 / (ninv2 * ninv2);
    }
    switch (direction) {
    case DesignElement::NO_DERIVATIVE:
    case DesignElement::ROTANGLE: {
      if (type_ == TRANSVERSAL_ISOTROPIC
          || type_ == DENSITY_TIMES_TRANSVERSAL_ISOTROPIC) {
        D = E * E3 / (E3 - nu13 * nu13 * E);
        D3 = E3 * E3 / (E3 - nu13 * nu13 * E);
        nD3 = nu13 * E * E3 / (E3 - nu13 * nu13 * E);
      } else if(type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED || type_ == TRANSVERSAL_ISOTROPIC_BOXED || DENSITY_TIMES_TRANSVERSAL_ISOTROPIC_BOXED){
        D = E / (1 - nu13);
        D3 = E3 / (1 - nu13);
        nD3 = sqrt(E * E3 * nu13) / (1 - nu13);
      }else{
        throw Exception("Not yet implemented!");
      }
      SetTransIsoTensor(t, subTensor, dens*D, 0, 0, dens*D3, dens*nD3, (type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED) ? dens*2.0*G3 : dens*G3);
      break;
    }
    case DesignElement::DENSITY:
    {
      if(type_ == DENSITY_TIMES_TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_TRANSVERSAL_ISOTROPIC_BOXED
          || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED)
      {
        dens = params_[DesignElement::DENSITY];
        if (penalty_ == 1.0) {
          dens = 1.0;
        } else
          dens = penalty_ * std::pow(dens, penalty_ - 1);
        if(notation == HILL_MANDEL_NO_DENSITY)
          dens = 0.0;
      }
      if (type_ == TRANSVERSAL_ISOTROPIC
          || type_ == DENSITY_TIMES_TRANSVERSAL_ISOTROPIC) {
        D = E * E3 / (E3 - nu13 * nu13 * E);
        D3 = E3 * E3 / (E3 - nu13 * nu13 * E);
        nD3 = nu13 * E * E3 / (E3 - nu13 * nu13 * E);
      } else if(type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED || type_ == TRANSVERSAL_ISOTROPIC_BOXED || DENSITY_TIMES_TRANSVERSAL_ISOTROPIC_BOXED){
        D = E / (1 - nu13);
        D3 = E3 / (1 - nu13);
        nD3 = sqrt(E * E3 * nu13) / (1 - nu13);
      }else{
        throw Exception("Not yet implemented!");
      }
      SetTransIsoTensor(t, subTensor, dens*D, 0, 0, dens*D3, dens*nD3, (type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED) ? dens*2.0*G3 : dens*G3);
      break;
    }
    case DesignElement::EMODULISO: {
      if (type_ == TRANSVERSAL_ISOTROPIC
          || type_ == DENSITY_TIMES_TRANSVERSAL_ISOTROPIC) {
        D = E3 * E3 * ninv2;
        D3 = nu13 * nu13 * E3 * E3 * ninv2;
        nD3 = nu13 * E3 * E3 * ninv2;
      } else if(type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED || type_ == TRANSVERSAL_ISOTROPIC_BOXED || DENSITY_TIMES_TRANSVERSAL_ISOTROPIC_BOXED){
        D = 1 / (1 - nu13);
        nD3 = sqrt(E3 * nu13 / E) / (2 - 2 * nu13);
      }else{
        throw Exception("Not yet implemented!");
      }
      SetTransIsoTensor(t, subTensor, dens * D, 0, 0, dens * D3, dens * nD3, 0);
      break;
    }
    case DesignElement::EMODUL: {
      if (type_ == TRANSVERSAL_ISOTROPIC
          || type_ == DENSITY_TIMES_TRANSVERSAL_ISOTROPIC) {
        D = -E * E * nu13 * nu13 * ninv2;
        D3 = -E3 * (-E3 + 2 * nu13 * nu13 * E) * ninv2;
        nD3 = -nu13 * nu13 * nu13 * E * E * ninv2;
      } else if(type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED || type_ == TRANSVERSAL_ISOTROPIC_BOXED || DENSITY_TIMES_TRANSVERSAL_ISOTROPIC_BOXED){
        D3 = 1 / (1 - nu13);
        nD3 = sqrt(E * nu13 / E3) / (2 - 2 * nu13);
      }else{
        throw Exception("Not yet implemented!");
      }
      SetTransIsoTensor(t, subTensor, dens * D, 0, 0, dens * D3, dens * nD3, 0);
      break;
    }
    case DesignElement::POISSON: {
      if (type_ == TRANSVERSAL_ISOTROPIC
          || type_ == DENSITY_TIMES_TRANSVERSAL_ISOTROPIC) {
        D = 2 * nu13 * E * E * E3 * ninv2;
        D3 = 2 * nu13 * E * E3 * E3 * ninv2;
        nD3 = E * E3 * (nu13 * nu13 * E + E3) * ninv2;
      } else if(type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED || type_ == TRANSVERSAL_ISOTROPIC_BOXED || DENSITY_TIMES_TRANSVERSAL_ISOTROPIC_BOXED){
        D = 1 / ((1 - nu13) * (1 - nu13));
        D3 = E3 * D;
        nD3 = sqrt(E * E3 / nu13) * (nu13 + 1) * D * 0.5;
        D = E * D;
      }else{
        throw Exception("Not yet implemented!");
      }
      SetTransIsoTensor(t, subTensor, dens * D, 0, 0, dens * D3, dens * nD3, 0);
      break;
    }
    case DesignElement::GMODUL: {
      SetTransIsoTensor(t, subTensor, 0, 0, 0, 0, 0, (type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED) ? 2*dens : dens);
      break;
    }
    default:
      ZeroTensor(t, subTensor);
      return;
    } // switch direction
    if(type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED){
      double rotAngle = params_[DesignElement::ROTANGLE];
      RotateHMStiffnessTensor(t, subTensor, direction, rotAngle, notation);
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
  double nu = params_[DesignElement::POISSONISO];
  double nu3;
  double n3;
  double c;
  double dens = 1.0, factor = 1.0;
  if((type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED || type_ == DENSITY_TIMES_ROT_PA12) && notation != HILL_MANDEL_NO_DENSITY){
    dens = params_[DesignElement::DENSITY];
    TransferFunction* tf = em_->GetDesign()->GetTransferFunction(DesignElement::DENSITY, Optimization::MECH);
    factor = (direction == DesignElement::DENSITY) ? tf->Derivative(dens) : tf->Transform(dens);
  } else {
    if(direction == DesignElement::DENSITY)
      factor = 0.0;
  }
  if(type_ == TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC){
    nu3 = nu13 * E3/E;
    n3 = nu3*nu3*E/E3;
    c = (1.0-nu-2.0*n3); // this is the interesting thing, this must not get 0, however this would imply a volume (trace of tensor) of infinity, so it is hopefully not occuring
    if(c < 1e-8) {
      c = 1e-8;
    }
  }else if(type_ == TRANSVERSAL_ISOTROPIC_BOXED || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED || type_ == DENSITY_TIMES_ROT_PA12){
    nu3 = sqrt(0.5*(1.0-nu)*E3/E)*nu13;
    n3 = nu3*nu3*E/E3;
    c = (1.0-nu-2.0*n3);
  }else{
    throw Exception("Not yet implemented!");
  }
  double f = E / ((1.0 + nu) * c);
  double dE = 0.0, dE3 = 0.0, dnu = 0.0, dnu3 = 0.0, dn3 = 0.0, dG3 = 0.0;
  
  bool tensorset = false;
  switch(direction){
  case DesignElement::NO_DERIVATIVE:
  case DesignElement::DENSITY: // almost the same as no derivative, we only changed the factor above
  case DesignElement::ROTANGLE:
  case DesignElement::ROTANGLEX:
  case DesignElement::ROTANGLEY:
  {
    double D = (1.0-n3)*f;
    double D3 = (1.0-nu)*E3/c;
    double nD = (nu+n3)*f;
    double nD3 = (1.0+nu)*nu3*f;
    double G = 0.5*E/(1.0+nu);
    SetTransIsoTensor(t, subTensor, factor * D, factor * nD, factor * G, factor * D3, factor * nD3, factor*G3);
    tensorset = true;
    break;
  }
  case DesignElement::EMODULISO:
    dE = 1.0;
    if(type_ == TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC){
      dnu3 = -E3*nu13/(E*E);
      dn3 = nu3/E3 * (2.0*E*dnu3 + nu3);
    }else if (type_ == TRANSVERSAL_ISOTROPIC_BOXED || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED){
      dnu3 = -sqrt(0.125*(1.0-nu)*E3/E)*nu13/E;
    }else if(type_ == DENSITY_TIMES_ROT_PA12){
      double D =  -(29*((nu - 1)*nu13*nu13 + 2))/(10*(nu*nu - 1)*(nu13*nu13 - 1));
      double D3 = -12/(5*(nu13*nu13 - 1));
      double nD =  (29*(nu*(nu13*nu13 - 2) - nu13*nu13))/(10*(nu*nu - 1)*(nu13*nu13 - 1));
      double nD3 = -(0.3*nu13*(4*E - 954.1))/(E*(nu13*nu13 - 1)*sqrt(E3*.5*(1-nu)/E));
      double G = -29/(10*(nu + 1));
      SetTransIsoTensor(t, subTensor, factor * D, factor * nD, factor * G, factor * D3, factor * nD3, 0.0);
      tensorset = true;
    }
    break;
  case DesignElement::EMODUL:
    dE3 = 1.0;
    if(type_ == TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC){
      dnu3 = nu13/E;
      dn3 = nu3*E/E3 * (2.0*dnu3 - nu3/E3);
    }else if (type_ == TRANSVERSAL_ISOTROPIC_BOXED || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED){
      dnu3 = sqrt(0.125*(1.0-nu)/(E*E3))*nu13;
    }else if(type_ == DENSITY_TIMES_ROT_PA12){
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
    if(type_ == TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC){
      dnu3 = 1.0;
      dn3 = 2.0 * nu3 * E / E3 * dnu3;
    }else if (type_ == TRANSVERSAL_ISOTROPIC_BOXED || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED || type_ == DENSITY_TIMES_ROT_PA12){
      dnu3 = sqrt(0.5 * (1.0 - nu) * E3 / E);
      dn3 = (1.0 - nu) * nu13;
    }
    break;
  case DesignElement::GMODUL:
    dG3 = 1.0;
    break;
  default:
    ZeroTensor(t, subTensor);
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
    SetTransIsoTensor(t, subTensor, factor * dD, factor * dnD, factor * dG, factor * dD3, factor * dnD3, factor * dG3);
  }
  
  if(type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC || type_ == DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED || type_ == DENSITY_TIMES_ROT_PA12){
    // for all rotated types, rotate the material tensor
    LOG_DBG3(dm) << "GetTransIsoMaterialTensor: tensor before rotation=" << t.ToString();
    RotateVoigtTensor(t, direction);
  }
  if(notation == HILL_MANDEL || notation == HILL_MANDEL_NO_DENSITY){
    double sq2 = sqrt(2);
    for(UInt i=0;i<3;++i){
      for(UInt j=3;j<6;++j){
        t(i,j)*=sq2;
        t(j,i)*=sq2;
        t(i+3,j)*=2;
      }
    }
  }
  LOG_DBG2(dm) << "GetTransIsoMaterialTensor: tensor result is " << t.ToString();
  
}
  
 
double DesignMaterial::GetTransIsoMaterialMass(DesignElement::Type direction){
  double E = params_[DesignElement::EMODULISO];
  double E3 = params_[DesignElement::EMODUL];
  double nu = params_[DesignElement::POISSONISO];
  double nu13 = params_[DesignElement::POISSON];
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
    double G3 = params_[DesignElement::GMODUL];
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

double DesignMaterial::GetDensityTimesTensorMass(DesignElement::Type direction){
  double dens = params_[DesignElement::DENSITY];
  switch (direction){
  case DesignElement::NO_DERIVATIVE:
  {
// for mass identity transfer function (which should normally be used) is assumed. This is messy because of lack of TransferFunctions!
    return dens;
  }
  case DesignElement::DENSITY:
  {
    return 1.0;
  }
  default:
    return 0.0;
  }
}

void DesignMaterial::GetOrthotropicMaterialTensor(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction, Notation notation){
  double e11 = params_[DesignElement::TENSOR11];
  double e22 = params_[DesignElement::TENSOR22];
  double e33 = params_[DesignElement::TENSOR33]; //This is already Hill-Mandel notation -> no scaling
  double e12 = params_[DesignElement::TENSOR12];
  double rotAngle = params_[DesignElement::ROTANGLE];
  double lowerEigBound = params_[DesignElement::LOWER_EIG_BOUND];
  double dens(1.0);
  if(type_ == DENSITY_TIMES_ORTHOTROPIC)
  {
    dens = params_[DesignElement::DENSITY];
    dens = std::pow(dens, penalty_);
  }

  if(subTensor == PLANE_STRESS){ //This is the only implemented case for now
    switch(direction){
    case DesignElement::NO_DERIVATIVE:
    case DesignElement::ROTANGLE:
    {
      SetTransIsoTensor(t, subTensor, dens*(e11*e11+e12*e12)+lowerEigBound, 0, 0,
          dens*(e12*e12+e22*e22)+lowerEigBound, dens*(e11*e12+e12*e22), dens*e33+lowerEigBound);
      break;
    }
    case DesignElement::DENSITY:
    {
      if(type_ == DENSITY_TIMES_ORTHOTROPIC)
      {
        dens = params_[DesignElement::DENSITY];
        if(penalty_ == 1.0){
          dens = 1.0;
        }else dens = penalty_*std::pow(dens, penalty_-1);
      }
      SetTransIsoTensor(t, subTensor, dens*(e11*e11+e12*e12), 0, 0, dens*(e12*e12+e22*e22), dens*(e11*e12+e12*e22), dens*e33);
      break;
    }
    case DesignElement::TENSOR11:
    {
      SetTransIsoTensor(t, subTensor, 2.0*dens*e11, 0, 0, 0, dens*e12, 0);
      break;
    }
    case DesignElement::TENSOR22:
    {
      SetTransIsoTensor(t, subTensor, 0, 0, 0, 2.0*dens*e22, dens*e12, 0);
      break;
    }
    case DesignElement::TENSOR12:
    {
      SetTransIsoTensor(t, subTensor, 2.0*dens*e12, 0, 0, 2.0*dens*e12, dens*(e11+e22), 0);
      break;
    }
    case DesignElement::TENSOR33:
    {
      SetTransIsoTensor(t, subTensor, 0, 0, 0, 0, 0, dens);
      break;
    }
    default:
      ZeroTensor(t, subTensor);
      return;
    } // switch direction
      RotateHMStiffnessTensor(t, subTensor, direction, rotAngle, notation);
      //    static int count(0);
      //    if (count % 10 == 0 && count/100 % 10 == 0){
      ////      std::cout << "(" << (count/100 % 10)*(count % 10)+1 << ")" << t.ToString() << std::endl;
      //      std::cout << t(0,0) << " " << t(0,1) << " " << t(1,1) << " " << t(0,2) << " " << t(1,2) << " " << t(2,2) << std::endl;
      //    }
      //    count++;
    return;
  } // PLANE_STRESS
  else
    throw Exception("subTensor not implemented yet");
}

void DesignMaterial::GetDensityTimes2dTensorTensor(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction)
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
    e11 = params_[DesignElement::TENSOR11];
    if (type_ == DENSITY_TIMES_2D_TENSOR_CONSTANT_TRACE) {
      e22 = params_[DesignElement::TENSOR22];
      e33 = 0.5 * (trace_ - e11 - e22);
    } else if (type_ == DENSITY_TIMES_ROTATED_2D_TENSOR) {
      e22 = 15 - e11;
      e11 += 1.0;
      e33 = params_[DesignElement::TENSOR33];
    } else {
      e22 = params_[DesignElement::TENSOR22];
      e33 = params_[DesignElement::TENSOR33];
    }
    e23 = params_[DesignElement::TENSOR23];
    e13 = params_[DesignElement::TENSOR13];
    e12 = params_[DesignElement::TENSOR12];
  }
  double d = params_[DesignElement::DENSITY];
  switch (direction) {
  case DesignElement::NO_DERIVATIVE:
  case DesignElement::ROTANGLE:
    d = pow(d, penalty_);
    Set2dVoigtTensor(t, d * e11, d * e22, d * e33, d * e23, d * e13, d * e12);
    break;
  case DesignElement::DENSITY:
    if (penalty_ == 1.0) {
      d = 1.0;
    } else {
      d = penalty_ * pow(d, penalty_ - 1);
    }
    Set2dVoigtTensor(t, d * e11, d * e22, d * e33, d * e23, d * e13, d * e12);
    break;
  case DesignElement::TENSOR11:
    d = pow(d, penalty_);
    Set2dVoigtTensor(t, d, type_ == DENSITY_TIMES_ROTATED_2D_TENSOR ? -d : 0.0,
        type_ == DENSITY_TIMES_2D_TENSOR_CONSTANT_TRACE ? -0.5 * d : 0.0, 0.0,
        0.0, 0.0);
    break;
  case DesignElement::TENSOR22:
    d = pow(d, penalty_);
    Set2dVoigtTensor(t, 0.0, d,
        type_ == DENSITY_TIMES_2D_TENSOR_CONSTANT_TRACE ? -0.5 * d : 0.0, 0.0,
        0.0, 0.0);
    break;
  case DesignElement::TENSOR33:
    d = pow(d, penalty_);
    Set2dVoigtTensor(t, 0.0, 0.0, d, 0.0, 0.0, 0.0);
    break;
  case DesignElement::TENSOR23:
    d = pow(d, penalty_);
    Set2dVoigtTensor(t, 0.0, 0.0, 0.0, d, 0.0, 0.0);
    break;
  case DesignElement::TENSOR13:
    d = pow(d, penalty_);
    Set2dVoigtTensor(t, 0.0, 0.0, 0.0, 0.0, d, 0.0);
    break;
  case DesignElement::TENSOR12:
    d = pow(d, penalty_);
    Set2dVoigtTensor(t, 0.0, 0.0, 0.0, 0.0, 0.0, d);
    break;
  default:
    ZeroTensor(t, subTensor);
    break;
  }
  if (type_ == DENSITY_TIMES_ROTATED_2D_TENSOR) {
    double rotAngle = params_[DesignElement::ROTANGLE];
    RotateHMStiffnessTensor(t, subTensor, direction, rotAngle);
//    static int count(0);
//    if (count % 10 == 0 && count/100 % 10 == 0){
////      std::cout << "(" << (count/100 % 10)*(count % 10)+1 << ")" << t.ToString() << std::endl;
//      std::cout << t(0,0) << " " << t(0,1) << " " << t(1,1) << " " << t(0,2) << " " << t(1,2) << " " << t(2,2) << std::endl;
//    }
//    count++;
    return;
  }
}


void DesignMaterial::GetElasticFMOTensor(Matrix<double>& E, DesignElement::Type direction, Notation notation)
{
  // We use the anisotropic tensor only for solving FMO problems. We assume the design to be in Hill-Mandel
  // notation and therefore we need to transform it for using it in CFS

  bool set = direction == DesignElement::NO_DERIVATIVE || direction == DesignElement::ROTANGLE;

  double e11 = set ? params_[DesignElement::TENSOR11] : 0;
  double e22 = set ? params_[DesignElement::TENSOR22] : 0;
  double e33 = set ? params_[DesignElement::TENSOR33] : 0;
  double e23 = set ? params_[DesignElement::TENSOR23] : 0;
  double e13 = set ? params_[DesignElement::TENSOR13] : 0;
  double e12 = set ? params_[DesignElement::TENSOR12] : 0;
  double rotAngle = set ? params_[DesignElement::ROTANGLE] : 0;

  switch (direction) {
  case DesignElement::NO_DERIVATIVE:
  case DesignElement::ROTANGLE:
    Set2dVoigtTensor(E, e11, e22, e33, e23, e13, e12);
    break;
  case DesignElement::TENSOR11:
    Set2dVoigtTensor(E, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    break;
  case DesignElement::TENSOR22:
    Set2dVoigtTensor(E, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0);
    break;
  case DesignElement::TENSOR33:
    Set2dVoigtTensor(E, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0);
    break;
  case DesignElement::TENSOR23:
    Set2dVoigtTensor(E, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    break;
  case DesignElement::TENSOR13:
    Set2dVoigtTensor(E, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    break;
  case DesignElement::TENSOR12:
    Set2dVoigtTensor(E, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    break;
  default:
    // for piezo FMO the derivative w.r.E. dielec_11, ... is zero
    ZeroTensor(E, PLANE_STRAIN);
  }

  LOG_DBG2(dm) << "GEFMOT: E before rotation = " << E.ToString(2) << " d=" << DesignElement::type.ToString(direction);
  RotateHMStiffnessTensor(E, PLANE, direction, rotAngle, notation);
  LOG_DBG2(dm) << "GEFMOT: E after rotation =  " << E.ToString(2) << " n=" << (notation == VOIGT ? "v" : "n") ;

}

void DesignMaterial::GetHomRectTensor(Matrix<double>& E,
    SubTensorType subTensor, DesignElement::Type direction, Notation notation) {
  // only relevant for hom_rect
  Quad9FE fe;
  double a = params_[DesignElement::STIFF1];
  double b = params_[DesignElement::STIFF2];
  double c = params_[DesignElement::SHEAR1];
  if (subTensor == FULL) {
    c = params_[DesignElement::STIFF3];
  }
  double rotAngle = 0.;
  if (subTensor != FULL) {
    rotAngle = params_[DesignElement::ROTANGLE];
  }
  Vector<double> p(3);
  if (type_ == HOM_RECT || type_ == D_HOM_RECT) {
    p[0] = -1.0 + 4 * a; // assume max 0.5
    p[1] = -1.0 + 4 * b; // assume max 0.5
  }
  if (type_ == HOM_RECT_C1) {
    p[0] = a;
    p[1] = b;
    p[2] = c;
  }
#ifndef NDEBUG
  Vector<double> peps(p);
  double eps = 1e-8;
  Matrix<double> Eeps(E);
  Matrix<double> Etmp(E);
  Vector<double> Diff(3);
#endif
  LOG_DBG2(dm)<< "GHRT: dir=" << (direction == DesignElement::NO_DERIVATIVE ? "no_derivative" : DesignElement::type.ToString(direction))
  << " not=" << notation << " rotAngle=" << rotAngle << " a=" << a << " b=" << b <<" c="<<(subTensor == FULL ? c : 0.0)<< " -> " << p.ToString();
  switch (direction) {
  case DesignElement::NO_DERIVATIVE:
  case DesignElement::ROTANGLE:
  case DesignElement::ROTANGLEX:
  case DesignElement::ROTANGLEY:
  case DesignElement::DENSITY:
  {
    if (type_ == HOM_RECT || type_ == D_HOM_RECT) {
      Vector<double> shape;
      fe.GetShFnc(shape, p, NULL);
      ApplyHomRectTensor(E, shape);
      LOG_DBG2(dm)<< "GHRT: shape=" << shape.ToString();
    }
    if (type_ == HOM_RECT_C1) {
      if (interpolation_ == C1) {
        ApplyHomRectC1Tensor(E,p,direction,subTensor);
      }
#ifdef USE_SGPP
      else if (interpolation_ == SGPP) {
        ApplyHomRectSGPPTensor(E,p,direction,subTensor);
      } else if (interpolation_ == FULL_BSPLINE) {
        ApplyHomRectFullBsplineTensor(E,p,direction,subTensor);
      }
#endif //USE_SGPP
    }
    break;
  }
  case DesignElement::STIFF1:
  case DesignElement::STIFF2:
  case DesignElement::STIFF3:
  case DesignElement::SHEAR1:
  {
    if(type_ == HOM_RECT || type_ == D_HOM_RECT) {
      Matrix<double> jac;
      Matrix<double> dummy; // not used -> strange function ?! :(
      fe.GetLocDerivShFnc(jac, p, dummy, NULL);
      LOG_DBG3(dm) << "GHRT: jac=" << jac.ToString(2);
      Vector<double> d_shape;
      jac.GetCol(d_shape, direction == DesignElement::STIFF1 ? 0 : 1);// a or by
      ApplyHomRectTensor(E, d_shape);
      // correct scaling to local FE coordinates
      E *= 4;
      LOG_DBG2(dm) << "GHRT: d_shape=" << d_shape.ToString();
    } else if(type_ == HOM_RECT_C1) {
      if (interpolation_ == C1) {
        ApplyHomRectC1Tensor(E,p,direction,subTensor);
      }
#ifdef USE_SGPP
      else if (interpolation_ == SGPP) {
        ApplyHomRectSGPPTensor(E,p,direction,subTensor);
      } else if (interpolation_ == FULL_BSPLINE) {
        ApplyHomRectFullBsplineTensor(E,p,direction,subTensor);
      }
#endif //USE_SGPP

      if (subTensor == FULL) {
#ifndef NDEBUG
        if (direction == DesignElement::STIFF1) {
          peps[0] += eps;
        } else if (direction == DesignElement::STIFF2) {
          peps[1] += eps;
        } else if (direction == DesignElement::STIFF3 || direction == DesignElement::SHEAR1) {
          peps[2] += eps;
        }
        ApplyHomRectC1Tensor(Eeps,peps,DesignElement::NO_DERIVATIVE,subTensor);
        ApplyHomRectC1Tensor(Etmp,p,DesignElement::NO_DERIVATIVE,subTensor);
        LOG_DBG(dm)<<"Eeps11: "<<std::setprecision(10)<<Eeps[0][0]<<", E11: "<<Etmp[0][0]<<" Diff: "<<Eeps[0][0]-Etmp[0][0];
        LOG_DBG(dm)<<"Eeps13: "<<std::setprecision(10)<<Eeps[0][2]<<", E13: "<<Etmp[0][2]<<" Diff: "<<Eeps[0][2]-Etmp[0][2];
        double e11 = (Eeps[0][0]-Etmp[0][0])/eps;
        double e12 = (Eeps[0][1]-Etmp[0][1])/eps;
        double e13 = (Eeps[0][2]-Etmp[0][2])/eps;
        double e22 = (Eeps[1][1]-Etmp[1][1])/eps;
        double e23 = (Eeps[1][2]-Etmp[1][2])/eps;
        double e33 = (Eeps[2][2]-Etmp[2][2])/eps;
        double e44 = (Eeps[3][3]-Etmp[3][3])/eps;
        double e55 = (Eeps[4][4]-Etmp[4][4])/eps;
        double e66 = (Eeps[5][5]-Etmp[5][5])/eps;
        LOG_DBG(dm)<<"FD Derivative "<<((direction == DesignElement::STIFF1)?"1":((direction == DesignElement::STIFF2) ? "2":"3"))<<" E11= "<<e11<<" E12= "<<e12<<" E22= "<< e22<<
        " E33= "<<e33<<" E23= "<<e23<<" E13= "<<e13<<" E44= "<<e44<<" E55= "<<e55<<" E66= "<<e66;
        LOG_DBG(dm)<<"deriv p= "<<p[0]<<", "<<p[1]<<", "<<p[2];
        LOG_DBG(dm)<<"FD Derivative - Derivative: "<<((direction == DesignElement::STIFF1)?"1":((direction == DesignElement::STIFF2) ? "2":"3"))<<" diff E11= "<<E[0][0]-e11<<" diff E12= "<<E[0][1]-e12<<" diff E22= "<< E[1][1]-e22<<
        " diff E33= "<<E[2][2]-e33<<" diff E23= "<<E[1][2]-e23<<" diff E13= "<<E[0][2]-e13<<" diff E44= "<<E[3][3]-e44<<" diff E55= "<<E[4][4]-e55<<" diff E66= "<<E[5][5]-e66;
#endif
      }
    }
    break;
  }
  default:
  ZeroTensor(E, subTensor == FULL ? FULL : PLANE);
}
  LOG_DBG2(dm)<< "GHRT: E before rotation = " << E.ToString(2);
  if (subTensor == FULL) {
    RotateVoigtTensor(E, direction);
  } else {
    RotateHMStiffnessTensor(E, PLANE, direction, rotAngle, notation);
  }
  LOG_DBG2(dm)<< "GHRT: E after rotation =  " << E.ToString(2);

  if (type_ == D_HOM_RECT)
  {
    double dens = params_[DesignElement::DENSITY];
    if (direction == DesignElement::DENSITY)
    {
      if(penalty_ == 1.0)
        dens = 1.0;
      else
        dens = penalty_*std::pow(dens, penalty_-1.0);
    }
    else
    {
      dens = std::pow(dens, penalty_);
    }
    E *= dens;
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

void DesignMaterial::GetRedBasCorrector(StdVector<Vector<double> >& corrector_, const Matrix<double>& G)
{
    //std::cout << "Computing corrector " << std::endl;

    assert ( (type_ == REDBAS_FREE) || (type_ == REDBAS_PARAM) || (type_ == REDBAS_MAPPING));

   corrector_ = StdVector<Vector<double> >(3 );
  //std::cout << "corrector"<< std::endl;
    double G11 = G[1-1][1-1];
    double G12 = G[1-1][2-1];
    double G21 = G[2-1][1-1];
    double G22 = G[2-1][2-1];

   //std::cout << "G11=" << G11 << " G12=" << G12 << " G21=" << G21 << " G22=" << G22;

  Matrix<double> mat(dimension_, dimension_);
  mat.Init();

  //mat =    G11*G11*matuxux11 + G11*G12*(matuxuy11 + matuxuy11t) + G12*G12*matuyuy11
  //    +    G21*G21*matvxvx11 + G21*G22*(matvxvy11 + matvxvy11t) + G22*G22*matvyvy11
  //    + G11*G21*(matuxvx12 + matuxvx12t) + G12*G21*(matuyvx12 + matuyvx12t) + G11*G22*(matuxvy12 + matuxvy12t) + G12*G22*(matuyvy12 + matuyvy12t)
  //    +    (G11*G11*matvxvx33 + G11*G12*(matvxvy33 + matvxvy33t) + G12*G12*matvyvy33)
  //       + (G21*G21*matuxux33 + G21*G22*(matuxuy33 + matuxuy33t) + G22*G22*matuyuy33)
  //    + ( G11*G21*(matuxvx33 + matuxvx33t) + G12*G21*(matuxvy33 + matuxvy33t) + G11*G22*(matuyvx33 + matuyvx33t) + G12*G22*(matuyvy33 + matuyvy33t))


  mat.Add(G11*G11, mod_red_matrices_[10*0 + 0]);
   mat.Add(G11*G12, mod_red_matrices_[10*0 + 3]);
   //mat.AddT(G11*G12, mod_red_matrices_[10*0 + 3]);
   mat.Add(G12*G12, mod_red_matrices_[10*0 + 7]);

   mat.Add(G21*G21, mod_red_matrices_[10*0 + 2]);
   mat.Add(G21*G22 ,mod_red_matrices_[10*0 + 6]);
   //mat.AddT(G21*G22,mod_red_matrices_[10*0 + 6]);
   mat.Add(G22*G22, mod_red_matrices_[10*0 + 9]);

   mat.Add(G11*G21, mod_red_matrices_[10*1 + 1] );
   //mat.AddT(G11*G21, mod_red_matrices_[10*1 + 1] );
   mat.Add(G12*G21, mod_red_matrices_[10*1 + 5]);
   //mat.AddT(G12*G21, mod_red_matrices_[10*1 + 5]);
   mat.Add(G11*G22, mod_red_matrices_[10*1 + 4]);
   //mat.AddT(G11*G22, mod_red_matrices_[10*1 + 4]);
   mat.Add(G12*G22, mod_red_matrices_[10*1 + 8]);
   //mat.AddT(G12*G22, mod_red_matrices_[10*1 + 8]);

   mat.Add(G21*G21, mod_red_matrices_[10*2 + 0]);
   mat.Add(G21*G22, mod_red_matrices_[10*2 + 3] );
   //mat.AddT(G21*G22, mod_red_matrices_[10*2 + 3] );
   mat.Add(G22*G22,mod_red_matrices_[10*2 + 7] );

   mat.Add(G11*G11,mod_red_matrices_[10*2 + 2] );
   mat.Add(G11*G12,mod_red_matrices_[10*2 + 6] );
   //mat.AddT(G11*G12,mod_red_matrices_[10*2 + 6] );
   mat.Add(G12*G12, mod_red_matrices_[10*2 + 9]);

   mat.Add(G11*G21, mod_red_matrices_[10*2 + 1]);
   //mat.AddT(G11*G21, mod_red_matrices_[10*2 + 1]);
   mat.Add(G11*G22, mod_red_matrices_[10*2 + 5]);
   //mat.AddT(G11*G22, mod_red_matrices_[10*2 + 5]);
   mat.Add(G12*G21, mod_red_matrices_[10*2 + 4]);
   //mat.AddT(G12*G21, mod_red_matrices_[10*2 + 4]);
   mat.Add(G12*G22, mod_red_matrices_[10*2 + 8]);
   //mat.AddT(G12*G22, mod_red_matrices_[10*2 + 8]);


  for (int corrector_int =0; corrector_int<3; corrector_int++)
  {

    Vector<double> vec(dimension_);

        if (corrector_int == 0)
        {
          //vec = G11*vecux11 + G12vecuy11 + G21*vecvx12 + G22*vecvy12
          vec = mod_red_vectors_[4*0 +0]*G11 + mod_red_vectors_[4*0 +1]*G12 +mod_red_vectors_[4*1 +2]*G21 + mod_red_vectors_[4*1 +3]*G22;
        }
        else if (corrector_int == 1)
        {
          //vec = G21*vecvx11 + G22*vecvy11 + G11*vecux12 + G12*vecuy12
          vec = mod_red_vectors_[4*0 +2]*G21 + mod_red_vectors_[4*0 +3]*G22 +mod_red_vectors_[4*1 +0]*G11 + mod_red_vectors_[4*1 +1]*G12;
        }
        else if (corrector_int == 2)
        {
          //vec=0.5*(G21*vecux33 + G22*vecuy33 + G11*vecvx33 + G12*vecvy33)
          vec = mod_red_vectors_[ 4*2 +0]*G21 + mod_red_vectors_[4*2 +1]*G22 +mod_red_vectors_[4*2 +2]*G11 + mod_red_vectors_[4*2 +3]*G12;
        }

        Vector<double> sol(dimension_);
        mat.DirectSolve(sol,vec);
        corrector_[corrector_int] = sol;

  }


}

void DesignMaterial::GetModRedHomTensor(Matrix<double>& E, const Matrix<double>& G, const StdVector<Vector<double> >& corrector_, Notation notation)
{

    //std::cout << "Computing tensor" << std::endl;
  E.Resize(3,3);
  E.Init();

 // std::cout << "GetModRedTensorHom" << std::endl;
  double G11 = G(0,0);
  double G12 = G(0,1);
  double G21 = G(1,0);
  double G22 = G(1,1);

  UInt dimbas =0;
  if ((type_ == REDBAS_PARAM) || (type_ == REDBAS_FREE) || (type_ == REDBAS_MAPPING)) dimbas = dimension_;
  if ((type_ == GREEDY_PARAM) || (type_ == GREEDY_FREE) || (type_ == GREEDY_MAPPING)) dimbas = 3*dimension_;

  Matrix<double> mat(dimbas, dimbas);
  mat.Init();

  mat.Add(G11*G11, mod_red_matrices_[10*0 + 0]);
   mat.Add(G11*G12, mod_red_matrices_[10*0 + 3]);
   //mat.AddT(G11*G12, mod_red_matrices_[10*0 + 3]);
   mat.Add(G12*G12, mod_red_matrices_[10*0 + 7]);

   mat.Add(G21*G21, mod_red_matrices_[10*0 + 2]);
   mat.Add(G21*G22 ,mod_red_matrices_[10*0 + 6]);
   //mat.AddT(G21*G22,mod_red_matrices_[10*0 + 6]);
   mat.Add(G22*G22, mod_red_matrices_[10*0 + 9]);

   mat.Add(G11*G21, mod_red_matrices_[10*1 + 1] );
   //mat.AddT(G11*G21, mod_red_matrices_[10*1 + 1] );
   mat.Add(G12*G21, mod_red_matrices_[10*1 + 5]);
   //mat.AddT(G12*G21, mod_red_matrices_[10*1 + 5]);
   mat.Add(G11*G22, mod_red_matrices_[10*1 + 4]);
   //mat.AddT(G11*G22, mod_red_matrices_[10*1 + 4]);
   mat.Add(G12*G22, mod_red_matrices_[10*1 + 8]);
   //mat.AddT(G12*G22, mod_red_matrices_[10*1 + 8]);

   mat.Add(G21*G21, mod_red_matrices_[10*2 + 0]);
   mat.Add(G21*G22, mod_red_matrices_[10*2 + 3] );
   //mat.AddT(G21*G22, mod_red_matrices_[10*2 + 3] );
   mat.Add(G22*G22,mod_red_matrices_[10*2 + 7] );

   mat.Add(G11*G11,mod_red_matrices_[10*2 + 2] );
   mat.Add(G11*G12,mod_red_matrices_[10*2 + 6] );
   //mat.AddT(G11*G12,mod_red_matrices_[10*2 + 6] );
   mat.Add(G12*G12, mod_red_matrices_[10*2 + 9]);

   mat.Add(G11*G21, mod_red_matrices_[10*2 + 1]);
   //mat.AddT(G11*G21, mod_red_matrices_[10*2 + 1]);
   mat.Add(G11*G22, mod_red_matrices_[10*2 + 5]);
   //mat.AddT(G11*G22, mod_red_matrices_[10*2 + 5]);
   mat.Add(G12*G21, mod_red_matrices_[10*2 + 4]);
   //mat.AddT(G12*G21, mod_red_matrices_[10*2 + 4]);
   mat.Add(G12*G22, mod_red_matrices_[10*2 + 8]);
   //mat.AddT(G12*G22, mod_red_matrices_[10*2 + 8]);



  Vector<double> mat0(dimbas);
  Vector<double> mat1(dimbas);
  Vector<double> mat2(dimbas);

  mat.Mult(corrector_[0], mat0);
  mat.Mult(corrector_[1], mat1);
  mat.Mult(corrector_[2], mat2);

  Vector<double> vec0(dimbas);
  Vector<double> vec1(dimbas);
  Vector<double> vec2(dimbas);

  vec0 = mod_red_vectors_[4*0+0]*G11 + mod_red_vectors_[4*0+1]*G12 + mod_red_vectors_[4*1+2]*G21 + mod_red_vectors_[4*1+3]*G22;
  vec1 = mod_red_vectors_[4*1+0]*G11 + mod_red_vectors_[4*1+1]*G12 + mod_red_vectors_[4*0+2]*G21 + mod_red_vectors_[4*0+3]*G22;
  vec2 = mod_red_vectors_[4*2+0]*G21 + mod_red_vectors_[4*2+1]*G22 + mod_red_vectors_[4*2+2]*G11 + mod_red_vectors_[4*2+3]*G12;


 // double p11 = corrector_[0]*vec1;
 // double p12 = corrector_[0]*vec2;
 // double p13 = corrector_[0]*vec3;

 // double p21 = corrector_[1]*vec1;
 // double p22 = corrector_[1]*vec2;
 // double p23 = corrector_[1]*vec3;

 // double p31 = corrector_[2]*vec1;
 // double p32 = corrector_[2]*vec2;
 // double p33 = corrector_[2]*vec3;

   //E11 = E11mean + E11*dxg(u1) + E12*dyg(v1)
   E[1-1][1-1] = mean_tensor_[1-1][0]

               + corrector_[0].Inner(mat0) - 2.0*(corrector_[0].Inner(vec0));



  // E12 = E12mean - E11*dxg(u2) - E12*dyg(v2)
   E[1-1][2-1] = mean_tensor_[2-1][0]

               + corrector_[1].Inner(mat0)  - corrector_[1].Inner(vec0) - corrector_[0].Inner(vec1);


   E[2-1][1-1] =  E[1-1][2-1] ;


  //E13 = - E11*dxg(u3) - E12*dyg(v3)
   E[1-1][3-1] = 0

             + corrector_[2].Inner(mat0) - corrector_[2].Inner(vec0) - corrector_[0].Inner(vec2);


   E[3-1][1-1] = E[1-1][3-1];

  //E22 = E11mean - E12*dxg(u2) - E11*dyg(v2)
   E[2-1][2-1] = mean_tensor_[1-1][0]

              + corrector_[1].Inner(mat1) - 2.0*(corrector_[1].Inner(vec1));


  //E23 = - E12*dxg(u3) - E11*dyg(v3)
   E[2-1][3-1] = 0

             + corrector_[2].Inner(mat1) - corrector_[2].Inner(vec1) - corrector_[1].Inner(vec2);

   E[3-1][2-1] = E[2-1][3-1];

  //E33 = E33mean - E33*0.5*(dyg(u3) + dxg(v3))
   E[3-1][3-1] = mean_tensor_[3-1][0]

               + corrector_[2].Inner(mat2) - 2.0*(corrector_[2].Inner(vec2));

   E[1-1][3-1] = sqrt(2)*E[1-1][3-1];
   E[3-1][1-1] = E[1-1][3-1];
   E[2-1][3-1] = sqrt(2)*E[2-1][3-1];
   E[3-1][2-1] = E[2-1][3-1];
   E[3-1][3-1] = 2*E[3-1][3-1];

   if(notation != HILL_MANDEL)
   {
     E[1-1][3-1] = 1.0/sqrt(2)*E[1-1][3-1];
     E[3-1][1-1] = E[1-1][3-1];
     E[2-1][3-1] = 1.0/sqrt(2)*E[2-1][3-1];
     E[3-1][2-1] = E[2-1][3-1];
     E[3-1][3-1] = 0.5*E[3-1][3-1];
   }
}


void DesignMaterial::GetModRedHomTensor(Matrix<double>& E, const Matrix<double>& G, const Matrix<double>& Gderiv, const StdVector<Vector<double> >& corrector_, Notation notation)
{

     //std::cout << "Computing derivative of tensor" << std::endl;
    //Should be done with Hill-Mandel notation
  E.Resize(3,3);
  E.Init();

  double G11 = G[1-1][1-1];
  double G12 = G[1-1][2-1];
  double G21 = G[2-1][1-1];
  double G22 = G[2-1][2-1];

  double G11d = Gderiv[1-1][1-1];
  double G12d = Gderiv[1-1][2-1];
  double G21d = Gderiv[2-1][1-1];
  double G22d = Gderiv[2-1][2-1];

  UInt dimbas =0;
  if ((type_ == REDBAS_PARAM) || (type_ == REDBAS_FREE) || (type_ == REDBAS_MAPPING)) dimbas = dimension_;
  if ((type_ == GREEDY_PARAM) || (type_ == GREEDY_FREE) || (type_ == GREEDY_MAPPING)) dimbas = 3*dimension_;

  Matrix<double> matd(dimbas, dimbas);
  matd.Init();

  matd.Add(G11d*G11 + G11*G11d, mod_red_matrices_[10*0 + 0]);
   matd.Add(G11d*G12 + G11*G12d, mod_red_matrices_[10*0 + 3]);
   //matd.AddT(G11d*G12 + G11*G12d, mod_red_matrices_[10*0 + 3]);
   matd.Add(G12d*G12 + G12*G12d, mod_red_matrices_[10*0 + 7]);

   matd.Add(G21d*G21 + G21*G21d, mod_red_matrices_[10*0 + 2]);
   matd.Add(G21d*G22 + G21*G22d ,mod_red_matrices_[10*0 + 6]);
   //matd.AddT(G21d*G22 + G21*G22d,mod_red_matrices_[10*0 + 6]);
   matd.Add(G22d*G22 + G22*G22d, mod_red_matrices_[10*0 + 9]);

   matd.Add(G11d*G21 + G11*G21d, mod_red_matrices_[10*1 + 1] );
   //matd.AddT(G11d*G21 + G11*G21d, mod_red_matrices_[10*1 + 1] );
   matd.Add(G12d*G21 + G12*G21d, mod_red_matrices_[10*1 + 5]);
   //matd.AddT(G12d*G21 + G12*G21d, mod_red_matrices_[10*1 + 5]);
   matd.Add(G11d*G22 + G11*G22d, mod_red_matrices_[10*1 + 4]);
   //matd.AddT(G11d*G22 + G11*G22d, mod_red_matrices_[10*1 + 4]);
   matd.Add(G12d*G22 + G12*G22d, mod_red_matrices_[10*1 + 8]);
   //matd.AddT(G12d*G22 + G12*G22d, mod_red_matrices_[10*1 + 8]);

   matd.Add(G21d*G21 + G21*G21d, mod_red_matrices_[10*2 + 0]);
   matd.Add(G21d*G22 + G21*G22d, mod_red_matrices_[10*2 + 3] );
   //matd.AddT(G21d*G22 + G21*G22d, mod_red_matrices_[10*2 + 3] );
   matd.Add(G22d*G22 + G22*G22d,mod_red_matrices_[10*2 + 7] );

   matd.Add(G11d*G11 + G11*G11d,mod_red_matrices_[10*2 + 2] );
   matd.Add(G11d*G12 + G11*G12d,mod_red_matrices_[10*2 + 6] );
   //matd.AddT(G11d*G12 + G11*G12d,mod_red_matrices_[10*2 + 6] );
   matd.Add(G12d*G12 + G12*G12d, mod_red_matrices_[10*2 + 9]);

   matd.Add(G11d*G21 + G11*G21d, mod_red_matrices_[10*2 + 1]);
   //matd.AddT(G11d*G21 + G11*G21d, mod_red_matrices_[10*2 + 1]);
   matd.Add(G11d*G22 + G11*G22d, mod_red_matrices_[10*2 + 5]);
   //matd.AddT(G11d*G22 + G11*G22d, mod_red_matrices_[10*2 + 5]);
   matd.Add(G12d*G21 + G12*G21d, mod_red_matrices_[10*2 + 4]);
   //matd.AddT(G12d*G21 + G12*G21d, mod_red_matrices_[10*2 + 4]);
   matd.Add(G12d*G22 + G12*G22d, mod_red_matrices_[10*2 + 8]);
   //matd.AddT(G12d*G22 + G12*G22d, mod_red_matrices_[10*2 + 8]);


   Vector<double> mat0(dimbas);
   Vector<double> mat1(dimbas);
   Vector<double> mat2(dimbas);

   matd.Mult(corrector_[0], mat0);
   matd.Mult(corrector_[1], mat1);
   matd.Mult(corrector_[2], mat2);


   Vector<double> vec0(dimbas);
    Vector<double> vec1(dimbas);
    Vector<double> vec2(dimbas);

    vec0 = mod_red_vectors_[4*0+0]*G11d + mod_red_vectors_[4*0+1]*G12d + mod_red_vectors_[4*1+2]*G21d + mod_red_vectors_[4*1+3]*G22d;
    vec1 = mod_red_vectors_[4*1+0]*G11d + mod_red_vectors_[4*1+1]*G12d + mod_red_vectors_[4*0+2]*G21d + mod_red_vectors_[4*0+3]*G22d;
    vec2 = mod_red_vectors_[4*2+0]*G21d + mod_red_vectors_[4*2+1]*G22d + mod_red_vectors_[4*2+2]*G11d + mod_red_vectors_[4*2+3]*G12d;



  //E11 =- E11*dxgd(u1) - E12*dygd(v1) - E11*dxgd(u1) - E12*dygd(v1) + corr0'*mat*corr0;
  E[1-1][1-1] = 0

      +corrector_[0].Inner(mat0) - 2.0*(corrector_[0].Inner(vec0));

 // E12 = - E11*dxgd(u2) - E12*dygd(v2) - E12*dxgd(u1) - E22*dygd(v1) + corr1'*mat*corr0;
  E[1-1][2-1] =  0

               +corrector_[1].Inner(mat0) - corrector_[1].Inner(vec0) - corrector_[0].Inner(vec1);

  E[2-1][1-1] =  E[1-1][2-1] ;


 //E13 = - E11*dxgd(u3) - E12*dygd(v3) - E33*0.5*(dygd(u1)+dxgd(v1)) + corr2'*mat*corr0;
  E[1-1][3-1] =  0

                + corrector_[2].Inner(mat0) - corrector_[2].Inner(vec0) - corrector_[0].Inner(vec2);


  E[3-1][1-1] = E[1-1][3-1];

 //E22 = - E12*dxgd(u2) - E11*dygd(v2)- E12*dxgd(u2) - E11*dygd(v2) + corr1'*mat*corr1;
  E[2-1][2-1] =  0

                + corrector_[1].Inner(mat1) - 2.0*(corrector_[1].Inner(vec1));

 //E23 = - E12*dxgd(u3) - E11*dygd(v3) - E33*05*(dygd(u2) + dxgd(v2)) + corr1'*mat*corr2
  E[2-1][3-1] =  0

              + corrector_[2].Inner(mat1) - corrector_[2].Inner(vec1) - corrector_[1].Inner(vec2);


  E[3-1][2-1] = E[2-1][3-1];

 //E33 = - E33*0.5*(dygd(u3) + dxgd(v3))- E33*0.5*(dygd(u3) + dxgd(v3)) + corr2'*mat*corr2;
  E[3-1][3-1] = 0

             +corrector_[2].Inner(mat2) - 2.0*(corrector_[2].Inner(vec2));

                E[1-1][3-1] = sqrt(2)*E[1-1][3-1];
                E[3-1][1-1] = E[1-1][3-1];
                E[2-1][3-1] = sqrt(2)*E[2-1][3-1];
                E[3-1][2-1] = E[2-1][3-1];
                E[3-1][3-1] = 2*E[3-1][3-1];

          if (notation != HILL_MANDEL)
          {
              E[1-1][3-1] = 1.0/sqrt(2)*E[1-1][3-1];
              E[3-1][1-1] = E[1-1][3-1];
              E[2-1][3-1] = 1.0/sqrt(2)*E[2-1][3-1];
              E[3-1][2-1] = E[2-1][3-1];
              E[3-1][3-1] = 0.5*E[3-1][3-1];
          }



}

void DesignMaterial::GetModRedGTensor(Matrix<double>& G, DesignElement::Type direction, const bool& all_param)
{

  G.Resize(2,2);
  G.Init();

  assert((type_== REDBAS_PARAM) || (type_== REDBAS_FREE) || (type_== GREEDY_PARAM) || (type_== GREEDY_FREE));

  if((type_ == REDBAS_PARAM) | (type_ == GREEDY_PARAM))
  {
  //Quad9FE fe;
  // std::cout << "GetModRedTensor" << std::endl;

  double theta = params_[DesignElement::ROTANGLE];
  double phi = params_[DesignElement::ROTANGLE2];
  double l1 = params_[DesignElement::SCALING1];
  double l2 = params_[DesignElement::SCALING2];

  switch(direction)
  {
     case DesignElement::NO_DERIVATIVE:
     {
       if (!all_param)
       {
          G(0,0) = l1*cos(phi); //matrix element G11
          G(0,1) = l1*sin(phi); //matrix element G12
          G(1,0) = -l2*sin(phi); //matrix element G21
          G(1,1) = l2*cos(phi); //matrix element G22
       }
       else
       {
           G(0,0) = l1*cos(theta)*cos(phi) -l2*sin(theta)*sin(phi); //matrix element G11
           G(0,1) = l1*cos(theta)*sin(phi) +l2*sin(theta)*cos(phi); //matrix element G12
           G(1,0) = -l1*sin(theta)*cos(phi) -l2*cos(theta)*sin(phi); //matrix element G21
           G(1,1) = -l1*sin(theta)*sin(phi) + l2*cos(theta)*cos(phi); //matrix element G22
       }
       break;
     }

     case DesignElement::ROTANGLE:
     {
         if (!all_param)
         {
           throw Exception("GetModRedGTensor should not be called when all_param_ = false and Design::Element direction = ROTANGLE");
         }
         else
         {
             G(0,0) = -l1*sin(theta)*cos(phi) -l2*cos(theta)*sin(phi); //matrix element G11
             G(0,1) = -l1*sin(theta)*sin(phi) +l2*cos(theta)*cos(phi); //matrix element G12
             G(1,0) = -l1*cos(theta)*cos(phi) + l2*sin(theta)*sin(phi); //matrix element G21
             G(1,1) = -l1*cos(theta)*sin(phi) - l2*sin(theta)*cos(phi); //matrix element G22
         }
       break;
     }

     case DesignElement::ROTANGLE2:
     {
       if (!all_param)
       {
         G(0,0) = -l1*sin(phi); //matrix element G11
         G(0,1) = l1*cos(phi); //matrix element G12
         G(1,0) = - l2*cos(phi); //matrix element G21
         G(1,1) = - l2*sin(phi); //matrix element G22
       }
       else
       {

           G(0,0) = -l1*cos(theta)*sin(phi) -l2*sin(theta)*cos(phi); //matrix element G11
           G(0,1) = l1*cos(theta)*cos(phi) - l2*sin(theta)*sin(phi); //matrix element G12
           G(1,0) = l1*sin(theta)*sin(phi) -l2*cos(theta)*cos(phi); //matrix element G21
           G(1,1) = -l1*sin(theta)*cos(phi) - l2*cos(theta)*sin(phi); //matrix element G22

       }
       break;
     }

     case DesignElement::SCALING1:
     {
       if (!all_param)
       {
         G(0,0) = cos(phi); //matrix element G11
         G(0,1) = sin(phi); //matrix element G12
         G(1,0) = 0.0; //matrix element G21
         G(1,1) = 0.0; //matrix element G22
       }
       else
       {



         G(0,0) = cos(theta)*cos(phi); //matrix element G11
         G(0,1) = cos(theta)*sin(phi); //matrix element G12
         G(1,0) = -sin(theta)*cos(phi); //matrix element G21
         G(1,1) = -sin(theta)*sin(phi); //matrix element G22
       }
       break;
     }

     case DesignElement::SCALING2:
     {
       if (!all_param)
       {
         G(0,0) =0.0;
         G(0,1) =0.0;
         G(1,0) =-sin(phi);
         G(1,1) = cos(phi);
       }
       else
       {

         G(0,0) =-sin(theta)*sin(phi);
         G(0,1) =sin(theta)*cos(phi);
         G(1,0) =-cos(theta)*sin(phi);
         G(1,1) = cos(theta)*cos(phi);
       }
       break;
     }

     default: //Zero matrix
       G(0,0) =0.0;
       G(0,1) =0.0;
       G(1,0) =0.0;
       G(1,1) =0.0;
     break;
     }

  }
  else if ((type_ == REDBAS_FREE) | (type_ == GREEDY_FREE))
  {
     double G11 = params_[DesignElement::G11];
     double G12 = params_[DesignElement::G12];
     double G21 = params_[DesignElement::G21];
     double G22 = params_[DesignElement::G22];

     switch(direction)
       {
          case DesignElement::NO_DERIVATIVE:
          {
             G(0,0) = G11; //matrix element G11
             G(0,1) = G12; //matrix element G12
             G(1,0) = G21; //matrix element G21
             G(1,1) = G22; //matrix element G22

            break;
          }
          case DesignElement::G11:
          {
              G(0,0) = 1.0; //matrix element G11
              G(0,1) = 0.0; //matrix element G12
              G(1,0) = 0.0; //matrix element G21
              G(1,1) = 0.0; //matrix element G22
            break;
          }
          case DesignElement::G12:
          {
             G(0,0) = 0.0; //matrix element G11
             G(0,1) = 1.0; //matrix element G12
             G(1,0) = 0.0; //matrix element G21
             G(1,1) = 0.0; //matrix element G22
            break;
          }
          case DesignElement::G21:
          {

            G(0,0) = 0.0; //matrix element G11
            G(0,1) = 0.0; //matrix element G12
            G(1,0) = 1.0; //matrix element G21
            G(1,1) = 0.0; //matrix element G22

            break;
          }
          case DesignElement::G22:
          {
            G(0,0) =0.0;
            G(0,1) =0.0;
            G(1,0) =0.0;
            G(1,1) =1.0;
            break;
          }
          default: //Identity matrix
            G(0,0) =0.0;
            G(0,1) =0.0;
            G(1,0) =0.0;
            G(1,1) =0.0;
            break;

          }
  }


}

void DesignMaterial::GetMappingTensor(Matrix<double>& E, DesignElement::Type direction, Notation notation)
{

    //std::string str = direction->ToString();
    //std::cout << "direction = " << str << std::endl;

    assert((type_== GREEDY_MAPPING) || (type_ == REDBAS_MAPPING));

    //assert( (direction==DesignElement::GX_0) || (direction==DesignElement::GY_0) || (direction==DesignElement::GX_PX)  || (direction==DesignElement::GY_PX)   || (direction==DesignElement::GX_PY)  || (direction==DesignElement::GY_PY)   || (direction==DesignElement::GX_PXY)  || (direction==DesignElement::GY_PXY)  || (direction==DesignElement::NO_DERIVATIVE)   );

    Matrix<double> G(2,2);
    //We begin with calculating the gradient of the mapping in the element
    GetMappingGradient(G);



     StdVector<Vector<double> > corrector_(3);
     if (type_ == GREEDY_MAPPING)
     {
         Vector<double> paramvec(4);
         GetModRedParamVector(paramvec);
         GetGreedyCorrector(corrector_, paramvec, true);
     }
     else if (type_ == REDBAS_MAPPING)
     {
         GetRedBasCorrector(corrector_, G);
     }

     if ((direction == DesignElement::NO_DERIVATIVE))
     {
       GetModRedHomTensor(E, G,corrector_, notation);
     }
     else
     {
       Matrix<double> Gderiv(2,2);
       GetMappingGradient(Gderiv,direction);

       GetModRedHomTensor(E, G, Gderiv, corrector_, notation);
     }

}

void DesignMaterial::GetMappingGradient(Matrix<double>& G)
{
  //Here, we assume that the additional layer is on the north-east part of the domain and we assume that the ordering of the nodes for each element is as follows:

  //    4____________3
  //    |            |
  //    |            |
  //    |            |
  //    |            |
  //    |____________|
  //    1            2
  //
  //
  //An element contains the value of the mappings gx and gy on its south_west node


  G.Resize(2,2);
  G.Init();


  int south_west = 0;
  int south_east = 1;
  int north_east = 2;
  int north_west = 3;


  assert((type_== GREEDY_MAPPING) || (type_ == REDBAS_MAPPING));

  DesignSpace* space = domain->GetErsatzMaterial();

  DesignElement* de0_x = space->Find(current_elem->elemNum, DesignElement::G_MAP_X);
  DesignElement* de0_y = space->Find(current_elem->elemNum, DesignElement::G_MAP_Y);

  assert( (de0_x != NULL) && (de0_y!= NULL ));

  if ( (de0_x->vicinity->HasNeighbor(VicinityElement::X_P)) && ((de0_x->vicinity->HasNeighbor(VicinityElement::Y_P))) )
  {
     //Here we check that de0_y has the same vicinity pattern
    assert((de0_y->vicinity->HasNeighbor(VicinityElement::X_P)) && ((de0_y->vicinity->HasNeighbor(VicinityElement::Y_P))));

    DesignElement* depx_x = de0_x->vicinity->GetNeighbour(VicinityElement::X_P);
    DesignElement* depx_y = de0_y->vicinity->GetNeighbour(VicinityElement::X_P);

    DesignElement* depy_x = de0_x->vicinity->GetNeighbour(VicinityElement::Y_P);
    DesignElement* depy_y = de0_y->vicinity->GetNeighbour(VicinityElement::Y_P);

    if (depx_x->vicinity->HasNeighbor(VicinityElement::Y_P))
    {
        assert( (depx_y->vicinity->HasNeighbor(VicinityElement::Y_P)) && (depy_x->vicinity->HasNeighbor(VicinityElement::X_P)) && (depy_y->vicinity->HasNeighbor(VicinityElement::X_P)) );

        DesignElement* depxy_x = (depx_x->vicinity->GetNeighbour(VicinityElement::Y_P));
        DesignElement* depxy_y = (depx_y->vicinity->GetNeighbour(VicinityElement::Y_P));


        //I need to know the coordinates of the nodes of the cells I am working with
        Matrix<double>  coords; // we ignore the n times constructs

        StdVector<unsigned int> connect = current_elem->connect;
        // do not use updated coordinates up to now!!
        domain->GetGrid()->GetElemNodesCoord(coords, connect, false);

        double x_0 = coords(0,south_west);
        double y_0 = coords(1,south_west);

        double x_px = coords(0,south_east);
        double y_px = coords(1,south_east);

        double x_pxy = coords(0,north_east);
        double y_pxy = coords(1,north_east);

        double x_py = coords(0,north_west);
        double y_py = coords(1,north_west);

        //std::cout << "south west: x = " << x_0 << " y = " << y_0 << std::endl;
        //std::cout << "south east: x = " << x_px << " y = " << y_px << std::endl;
        //std::cout << "north west: x = " << x_py << " y = " << y_py << std::endl;
        //std::cout << "north east: x = " << x_pxy << " y = " << y_pxy << std::endl;


        double gx_0 = de0_x->GetDesign(DesignElement::SMART);
        double gy_0 = de0_y->GetDesign(DesignElement::SMART);

        double gx_px = depx_x->GetDesign(DesignElement::SMART);
        double gy_px = depx_y->GetDesign(DesignElement::SMART);


        double gx_py = depy_x->GetDesign(DesignElement::SMART);
        double gy_py = depy_y->GetDesign(DesignElement::SMART);

        double gx_pxy = depxy_x->GetDesign(DesignElement::SMART);
        double gy_pxy = depxy_y->GetDesign(DesignElement::SMART);

        // G(0,0) = d_x g_x
        // G(0,1) = d_y g_x
        // G(1,0) = d_x g_y
        // G(1,1) = d_y g_y

        G(0,0) = ( (gx_px - gx_0)*(x_px -x_0) + (gx_pxy - gx_px)*(x_pxy - x_px) + (gx_pxy - gx_py)*(x_pxy - x_py) + (gx_py - gx_0)*(x_py - x_0))/((x_px -x_0)*(x_px -x_0) + (x_pxy - x_px)*(x_pxy - x_px) + (x_pxy - x_py)*(x_pxy - x_py) + (x_py - x_0)*(x_py - x_0) );
        G(0,1) = ( (gx_px - gx_0)*(y_px -y_0) + (gx_pxy - gx_px)*(y_pxy - y_px) + (gx_pxy - gx_py)*(y_pxy - y_py) + (gx_py - gx_0)*(y_py - y_0))/((y_px -y_0)*(y_px -y_0) + (y_pxy - y_px)*(y_pxy - y_px) + (y_pxy - y_py)*(y_pxy - y_py) + (y_py - y_0)*(y_py - y_0) );
        G(1,0) = ( (gy_px - gy_0)*(x_px -x_0) + (gy_pxy - gy_px)*(x_pxy - x_px) + (gy_pxy - gy_py)*(x_pxy - x_py) + (gy_py - gy_0)*(x_py - x_0))/((x_px -x_0)*(x_px -x_0) + (x_pxy - x_px)*(x_pxy - x_px) + (x_pxy - x_py)*(x_pxy - x_py) + (x_py - x_0)*(x_py - x_0) );
        G(1,1) = ( (gy_px - gy_0)*(y_px -y_0) + (gy_pxy - gy_px)*(y_pxy - y_px) + (gy_pxy - gy_py)*(y_pxy - y_py) + (gy_py - gy_0)*(y_py - y_0))/((y_px -y_0)*(y_px -y_0) + (y_pxy - y_px)*(y_pxy - y_px) + (y_pxy - y_py)*(y_pxy - y_py) + (y_py - y_0)*(y_py - y_0) );


        //std::cout << "Elem = " << current_elem->elemNum << std::endl;
        //std::cout << "G = " << G << std::endl;
    }
    else
    {
      std::cout << "WARNING!!  We are computing the gradient of the mapping in cells of the ghost layer!!" << std::endl;
    }

  }
  else
  {
    std::cout << "WARNING!!  We are computing the gradient of the mapping in cells of the ghost layer!!" << std::endl;
  }

}


void DesignMaterial::GetMappingGradient(Matrix<double>& G, DesignElement::Type direction)
{
  //Here, we assume that the additional layer is on the north-east part of the domain and we assume that the ordering of the nodes for each element is as follows:

  //    4____________3
  //    |            |
  //    |            |
  //    |            |
  //    |            |
  //    |____________|
  //    1            2
  //
  //
  //An element contains the value of the mappings gx and gy on its south_west node


  G.Resize(2,2);
  G.Init();


  int south_west = 0;
  int south_east = 1;
  int north_east = 2;
  int north_west = 3;


  assert((type_== GREEDY_MAPPING) || (type_ == REDBAS_MAPPING));

  DesignSpace* space = domain->GetErsatzMaterial();

  DesignElement* de0_x = space->Find(current_elem->elemNum, DesignElement::G_MAP_X);
  DesignElement* de0_y = space->Find(current_elem->elemNum, DesignElement::G_MAP_Y);

  // otherwise de0_y is not used => error in release version
  de0_y = de0_y + 1;
  de0_y = de0_y - 1;

  //assert( (de0_x != NULL) && (de0_y!= NULL ));

  if ( (de0_x->vicinity->HasNeighbor(VicinityElement::X_P)) && ((de0_x->vicinity->HasNeighbor(VicinityElement::Y_P))) )
  {
     //Here we check that de0_y has the same vicinity pattern
  //  assert((de0_y->vicinity->HasNeighbor(VicinityElement::X_P)) && ((de0_y->vicinity->HasNeighbor(VicinityElement::Y_P))));

    DesignElement* depx_x = de0_x->vicinity->GetNeighbour(VicinityElement::X_P);
    //DesignElement* depx_y = de0_y->vicinity->GetNeighbour(VicinityElement::X_P);

    //DesignElement* depy_x = de0_x->vicinity->GetNeighbour(VicinityElement::Y_P);
    //DesignElement* depy_y = de0_y->vicinity->GetNeighbour(VicinityElement::Y_P);

    if (depx_x->vicinity->HasNeighbor(VicinityElement::Y_P))
    {
        //assert( (depx_y->vicinity->HasNeighbor(VicinityElement::Y_P)) && (depy_x->vicinity->HasNeighbor(VicinityElement::X_P)) && (depy_y->vicinity->HasNeighbor(VicinityElement::X_P)) );

        //We do not need these anymore
        //DesignElement* depxy_x = (depx_x->vicinity->GetNeighbour(VicinityElement::Y_P));
        //DesignElement* depxy_y = (depx_y->vicinity->GetNeighbour(VicinityElement::Y_P));


        //I need to know the coordinates of the nodes of the cells I am working with
        Matrix<double>  coords; // we ignore the n times constructs

        StdVector<unsigned int> connect = current_elem->connect;
        // do not use updated coordinates up to now!!
        domain->GetGrid()->GetElemNodesCoord(coords, connect, false);

        double x_0 = coords(0,south_west);
        double y_0 = coords(1,south_west);

        double x_px = coords(0,south_east);
        double y_px = coords(1,south_east);

        double x_pxy = coords(0,north_east);
        double y_pxy = coords(1,north_east);

        double x_py = coords(0,north_west);
        double y_py = coords(1,north_west);

        switch(direction)
        {

          case DesignElement::GX_0:
          {

              G(0,0) = ( -(x_px -x_0) + 0 + 0 -(x_py - x_0))/((x_px -x_0)*(x_px -x_0) + (x_pxy - x_px)*(x_pxy - x_px) + (x_pxy - x_py)*(x_pxy - x_py) + (x_py - x_0)*(x_py - x_0) );
              G(0,1) = ( -(y_px -y_0) + 0 + 0 -(y_py - y_0))/((y_px -y_0)*(y_px -y_0) + (y_pxy - y_px)*(y_pxy - y_px) + (y_pxy - y_py)*(y_pxy - y_py) + (y_py - y_0)*(y_py - y_0) );
              G(1,0) = 0;
              G(1,1) = 0;


              //std::cout << "Elem = " << current_elem->elemNum << std::endl;
              //std::cout << "GX_0 = " << G << std::endl;

                break;
            }


            case DesignElement::GX_PX:
            {

              G(0,0) = ( (x_px -x_0) -(x_pxy - x_px) + 0 + 0)/((x_px -x_0)*(x_px -x_0) + (x_pxy - x_px)*(x_pxy - x_px) + (x_pxy - x_py)*(x_pxy - x_py) + (x_py - x_0)*(x_py - x_0) );
              G(0,1) = ( (y_px -y_0) -(y_pxy - y_px) + 0 + 0)/((y_px -y_0)*(y_px -y_0) + (y_pxy - y_px)*(y_pxy - y_px) + (y_pxy - y_py)*(y_pxy - y_py) + (y_py - y_0)*(y_py - y_0) );
              G(1,0) = 0;
              G(1,1) = 0;

              //std::cout << "Elem = " << current_elem->elemNum << std::endl;
              //std::cout << "GX_PX = " << G << std::endl;

              break;
            }

            case DesignElement::GX_PY:
            {

              G(0,0) = ( 0+ 0 -(x_pxy - x_py) + (x_py - x_0))/((x_px -x_0)*(x_px -x_0) + (x_pxy - x_px)*(x_pxy - x_px) + (x_pxy - x_py)*(x_pxy - x_py) + (x_py - x_0)*(x_py - x_0) );
              G(0,1) = ( 0+ 0 -(y_pxy - y_py) + (y_py - y_0))/((y_px -y_0)*(y_px -y_0) + (y_pxy - y_px)*(y_pxy - y_px) + (y_pxy - y_py)*(y_pxy - y_py) + (y_py - y_0)*(y_py - y_0) );
              G(1,0) = 0;
              G(1,1) = 0;

              //std::cout << "Elem = " << current_elem->elemNum << std::endl;
              //std::cout << "GX_PY = " << G << std::endl;

              break;
            }

            case DesignElement::GX_PXY:
            {

              G(0,0) = ( 0 + (x_pxy - x_px) + (x_pxy - x_py) + 0)/((x_px -x_0)*(x_px -x_0) + (x_pxy - x_px)*(x_pxy - x_px) + (x_pxy - x_py)*(x_pxy - x_py) + (x_py - x_0)*(x_py - x_0) );
              G(0,1) = ( 0 + (y_pxy - y_px) + (y_pxy - y_py) +0)/((y_px -y_0)*(y_px -y_0) + (y_pxy - y_px)*(y_pxy - y_px) + (y_pxy - y_py)*(y_pxy - y_py) + (y_py - y_0)*(y_py - y_0) );
              G(1,0) = 0;
              G(1,1) = 0;


              //std::cout << "Elem = " << current_elem->elemNum << std::endl;
              //std::cout << "GX_PXY = " << G << std::endl;

              break;
            }


            case DesignElement::GY_0:
            {

              G(0,0) = 0;
              G(0,1) = 0;
              G(1,0) = ( -(x_px -x_0) + 0 + 0 -(x_py - x_0))/((x_px -x_0)*(x_px -x_0) + (x_pxy - x_px)*(x_pxy - x_px) + (x_pxy - x_py)*(x_pxy - x_py) + (x_py - x_0)*(x_py - x_0) );
              G(1,1) = ( -(y_px -y_0) + 0 + 0 + -(y_py - y_0))/((y_px -y_0)*(y_px -y_0) + (y_pxy - y_px)*(y_pxy - y_px) + (y_pxy - y_py)*(y_pxy - y_py) + (y_py - y_0)*(y_py - y_0) );

              //std::cout << "Elem = " << current_elem->elemNum << std::endl;
              //std::cout << "GY_0 = " << G << std::endl;

                break;
            }


            case DesignElement::GY_PX:
            {

              G(0,0) = 0;
              G(0,1) = 0;
              G(1,0) = ( (x_px -x_0) -(x_pxy - x_px) + 0 + 0)/((x_px -x_0)*(x_px -x_0) + (x_pxy - x_px)*(x_pxy - x_px) + (x_pxy - x_py)*(x_pxy - x_py) + (x_py - x_0)*(x_py - x_0) );
              G(1,1) = ( (y_px -y_0) -(y_pxy - y_px) + 0 + 0)/((y_px -y_0)*(y_px -y_0) + (y_pxy - y_px)*(y_pxy - y_px) + (y_pxy - y_py)*(y_pxy - y_py) + (y_py - y_0)*(y_py - y_0) );

              //std::cout << "Elem = " << current_elem->elemNum << std::endl;
              //std::cout << "GY_PX = " << G << std::endl;

              break;
            }

            case DesignElement::GY_PY:
            {

              G(0,0) = 0;
              G(0,1) = 0;
              G(1,0) = ( 0 + 0 -(x_pxy - x_py) + (x_py - x_0))/((x_px -x_0)*(x_px -x_0) + (x_pxy - x_px)*(x_pxy - x_px) + (x_pxy - x_py)*(x_pxy - x_py) + (x_py - x_0)*(x_py - x_0) );
              G(1,1) = ( 0 + 0 -(y_pxy - y_py) + (y_py - y_0))/((y_px -y_0)*(y_px -y_0) + (y_pxy - y_px)*(y_pxy - y_px) + (y_pxy - y_py)*(y_pxy - y_py) + (y_py - y_0)*(y_py - y_0) );

              //std::cout << "Elem = " << current_elem->elemNum << std::endl;
              //std::cout << "GY_PY = " << G << std::endl;

              break;
            }

            case DesignElement::GY_PXY:
            {
              G(0,0) = 0;
              G(0,1) = 0;
              G(1,0) = ( 0 + (x_pxy - x_px) + (x_pxy - x_py) + 0)/((x_px -x_0)*(x_px -x_0) + (x_pxy - x_px)*(x_pxy - x_px) + (x_pxy - x_py)*(x_pxy - x_py) + (x_py - x_0)*(x_py - x_0) );
              G(1,1) = ( 0 + (y_pxy - y_px) + (y_pxy - y_py) + 0)/((y_px -y_0)*(y_px -y_0) + (y_pxy - y_px)*(y_pxy - y_px) + (y_pxy - y_py)*(y_pxy - y_py) + (y_py - y_0)*(y_py - y_0) );

              //std::cout << "Elem = " << current_elem->elemNum << std::endl;
              //std::cout << "GY_PXY = " << G << std::endl;

              break;
            }
          default:
          //Zero matrix

              std::cout << "returning zero matrix" << std::endl;
            G(0,0) =0.0;
            G(0,1) =0.0;
            G(1,0) =0.0;
            G(1,1) =0.0;
            break;

        }


    }
    else
    {
      std::cout << "WARNING!!  We are computing the gradient of the mapping in cells of the ghost layer!!" << std::endl;
    }

  }
  else
  {
    std::cout << "WARNING!!  We are computing the gradient of the mapping in cells of the ghost layer!!" << std::endl;
  }

}




void DesignMaterial::GetModRedTensor(Matrix<double>& E, DesignElement::Type direction, Notation notation)
{

  //if type_ == REDBAS_FREE or type_ == GREEDY_FREE then we must have all_param_
  assert((type_ == GREEDY_FREE) || (type_ == REDBAS_FREE) || (type_ == REDBAS_PARAM) || (type_ == GREEDY_PARAM) );

  E.Resize(3,3);
  E.Init();

  Matrix<double> G(2,2);
  GetModRedGTensor(G,DesignElement::NO_DERIVATIVE, all_param_);

  Vector<double> paramvec(4);
  GetModRedParamVector(paramvec);

  double theta = paramvec[0];

  StdVector<Vector<double> > corrector_(3);
  //The corrector calculation is always done with Voigt notation
  if ((type_ == REDBAS_PARAM) | (type_== REDBAS_FREE)) GetRedBasCorrector(corrector_, G);
  else if ((type_ == GREEDY_PARAM) | (type_== GREEDY_FREE)) GetGreedyCorrector(corrector_, paramvec, all_param_);

  if ((direction == DesignElement::NO_DERIVATIVE) | ((!all_param_) & (direction == DesignElement::ROTANGLE)))
  {
    GetModRedHomTensor(E, G,corrector_, notation);
  }
  else
  {
    Matrix<double> Gderiv(2,2);
    GetModRedGTensor(Gderiv,direction, all_param_);
    GetModRedHomTensor(E, G, Gderiv, corrector_, notation);
  }
  if (!all_param_){
    RotateHMStiffnessTensor(E, PLANE, direction, -theta, notation);
    if (direction == DesignElement::ROTANGLE) E=-E;

  }


}


void DesignMaterial::GetModRedParamVector(Vector<double>& params)
{

  params.Resize(4,1);
  params.Init();

  if ((type_==GREEDY_PARAM) | (type_ == REDBAS_PARAM))
  {
    //Quad9FE fe;
    // std::cout << "GetModRedTensor" << std::endl;

    double theta = params_[DesignElement::ROTANGLE];
    double phi = params_[DesignElement::ROTANGLE2];
    double l1 = params_[DesignElement::SCALING1];
    double l2 = params_[DesignElement::SCALING2];

    params[0]=theta;
    params[1]=phi;
    params[2]=l1;
    params[3]=l2;

  }
  else if ((type_==GREEDY_FREE) | (type_ == REDBAS_FREE))
  {
      Matrix<double> G(2,2);
      GetModRedGTensor(G,DesignElement::NO_DERIVATIVE, true);

      GetSVDGTensorParameters(G, params);

  }
  else if ((type_ == GREEDY_MAPPING) || (type_ == REDBAS_MAPPING))
  {
    Matrix<double> G(2,2);
      GetMappingGradient(G);

      GetSVDGTensorParameters(G, params);
  }
}


/*void DesignMaterial::GetGreedyTensor(Matrix<double>& E, DesignElement::Type direction, Notation notation)
{

  E.Resize(3,3);
  E.Init();

  Matrix<double> params(4,1);
  GetGreedyParams(params);

  GetGreedyTensor(E,params,all_param_, direction,notation);

}*/



void DesignMaterial::GetGreedyCorrector(StdVector<Vector<double> >& corrector_, const Vector<double>& params, const bool& all_param)
{

  assert((type_== GREEDY_PARAM) || (type_ == GREEDY_FREE) || (type_ == GREEDY_MAPPING));

  double theta = params[0];
  double phi = params[1];
  double l1 = params[2];
  double l2 = params[3];

  if (l1 < lmin_) l1 = lmin_;
  else if (l1 > lmax_) l1 = lmax_;

  if (l2 < lmin_) l2 = lmin_;
  else if (l2 > lmax_) l2 = lmax_;

  corrector_ = StdVector<Vector<double> >(3);
  corrector_[0] = Vector<double>(3*dimension_);
  corrector_[1] = Vector<double>(3*dimension_);
  corrector_[2] = Vector<double>(3*dimension_);


  // corrector_[0] = Matrix<double>(dimension_, 1);
  // corrector_[1] = Matrix<double>(dimension_, 1);
  // corrector_[2] = Matrix<double>(dimension_, 1);

  //We begin by computing the coefficients of the three correctors in the reduced basis associated to the value of these parameters
  if (all_param)
  {

      for (UInt k=0; k<dimension_; k++)
      {
        Vector<double> coord1theta(2*Na_+1);
        Vector<double> coord2theta(2*Na_+1);
        Vector<double> coord3theta(2*Na_+1);
        Vector<double> coord1phi(2*Na_+1);
        Vector<double> coord2phi(2*Na_+1);
        Vector<double> coord3phi(2*Na_+1);
        Vector<double> coord1l1(Nl_+1);
        Vector<double> coord2l1(Nl_+1);
        Vector<double> coord3l1(Nl_+1);
        Vector<double> coord1l2(Nl_+1);
        Vector<double> coord2l2(Nl_+1);
        Vector<double> coord3l2(Nl_+1);

        matrices_param_[0].GetCol(coord1theta, k);
        matrices_param_[0].GetCol(coord2theta, dimension_+ k);
        matrices_param_[0].GetCol(coord3theta, 2*dimension_+ k);

        matrices_param_[1].GetCol(coord1phi, k);
        matrices_param_[1].GetCol(coord2phi, dimension_+ k);
        matrices_param_[1].GetCol(coord3phi, 2*dimension_+ k);

        matrices_param_[2].GetCol(coord1l1, k);
        matrices_param_[2].GetCol(coord2l1, dimension_+ k);
        matrices_param_[2].GetCol(coord3l1, 2*dimension_+ k);

        matrices_param_[3].GetCol(coord1l2, k);
        matrices_param_[3].GetCol(coord2l2, dimension_+ k);
        matrices_param_[3].GetCol(coord3l2, 2*dimension_+ k);

        //corrector_[0][k] = -AngleGreedyCalculus(coord1theta, theta)*AngleGreedyCalculus(coord1phi, phi)*ScalingGreedyCalculus(coord1l1,l1)*ScalingGreedyCalculus(coord1l2,l2);
        //corrector_[1][k] = -AngleGreedyCalculus(coord2theta, theta)*AngleGreedyCalculus(coord2phi, phi)*ScalingGreedyCalculus(coord2l1,l1)*ScalingGreedyCalculus(coord2l2,l2);
        //corrector_[2][k] = -AngleGreedyCalculus(coord3theta, theta)*AngleGreedyCalculus(coord3phi, phi)*ScalingGreedyCalculus(coord3l1,l1)*ScalingGreedyCalculus(coord3l2,l2);

        corrector_[0][k] = -AngleGreedyCalculus(coord1theta, theta)*AngleGreedyCalculus(coord1phi, phi)*ScalingGreedyCalculus(coord1l1,l1)*ScalingGreedyCalculus(coord1l2,l2);
        corrector_[1][dimension_ +k] = -AngleGreedyCalculus(coord2theta, theta)*AngleGreedyCalculus(coord2phi, phi)*ScalingGreedyCalculus(coord2l1,l1)*ScalingGreedyCalculus(coord2l2,l2);
        corrector_[2][2*dimension_+k] = -AngleGreedyCalculus(coord3theta, theta)*AngleGreedyCalculus(coord3phi, phi)*ScalingGreedyCalculus(coord3l1,l1)*ScalingGreedyCalculus(coord3l2,l2);


      }


  }
  else
  {
          for (UInt k=0; k<dimension_; k++)
          {
            Vector<double> coord1phi(2*Na_+1);
            Vector<double> coord2phi(2*Na_+1);
            Vector<double> coord3phi(2*Na_+1);
            Vector<double> coord1l1(Nl_+1);
            Vector<double> coord2l1(Nl_+1);
            Vector<double> coord3l1(Nl_+1);
            Vector<double> coord1l2(Nl_+1);
            Vector<double> coord2l2(Nl_+1);
            Vector<double> coord3l2(Nl_+1);

            matrices_param_[0].GetCol(coord1phi, k);
            matrices_param_[0].GetCol(coord2phi, dimension_+ k);
            matrices_param_[0].GetCol(coord3phi, 2*dimension_+ k);

            matrices_param_[1].GetCol(coord1l1, k);
            matrices_param_[1].GetCol(coord2l1, dimension_+ k);
            matrices_param_[1].GetCol(coord3l1, 2*dimension_+ k);

            matrices_param_[2].GetCol(coord1l2, k);
            matrices_param_[2].GetCol(coord2l2, dimension_+ k);
            matrices_param_[2].GetCol(coord3l2, 2*dimension_+ k);

            //corrector_[0][k] = -AngleGreedyCalculus(coord1phi, phi)*ScalingGreedyCalculus(coord1l1,l1)*ScalingGreedyCalculus(coord1l2,l2);
            //corrector_[1][k] = -AngleGreedyCalculus(coord2phi, phi)*ScalingGreedyCalculus(coord2l1,l1)*ScalingGreedyCalculus(coord2l2,l2);
            //corrector_[2][k] = -AngleGreedyCalculus(coord3phi, phi)*ScalingGreedyCalculus(coord3l1,l1)*ScalingGreedyCalculus(coord3l2,l2);

            corrector_[0][k] = -AngleGreedyCalculus(coord1phi, phi)*ScalingGreedyCalculus(coord1l1,l1)*ScalingGreedyCalculus(coord1l2,l2);
            corrector_[1][dimension_ +k] = -AngleGreedyCalculus(coord2phi, phi)*ScalingGreedyCalculus(coord2l1,l1)*ScalingGreedyCalculus(coord2l2,l2);
            corrector_[2][2*dimension_+k] = -AngleGreedyCalculus(coord3phi, phi)*ScalingGreedyCalculus(coord3l1,l1)*ScalingGreedyCalculus(coord3l2,l2);

          }

  }
}

void DesignMaterial::ApplyHomRectTensor(Matrix<double>& E,
    const Vector<double>& shape) const {
  E.Resize(3, 3);
  E.Init(); // for off-diagonal
  Vector<double> data;
  hom_rect_samples_.GetCol(data,
      DesignElement::TENSOR11 - DesignElement::TENSOR11);
  
  E[1 - 1][1 - 1] = shape * data;
  LOG_DBG(dm)<< "AHRT 11=" << E[1-1][1-1] << " data=" << data.ToString();
  hom_rect_samples_.GetCol(data,
      DesignElement::TENSOR12 - DesignElement::TENSOR11);
  E[1 - 1][2 - 1] = shape * data;
  E[2 - 1][1 - 1] = E[1 - 1][2 - 1];
  LOG_DBG(dm)<< "AHRT 12=" << E[1-1][2-1] << " data=" << data.ToString();
  hom_rect_samples_.GetCol(data,
      DesignElement::TENSOR22 - DesignElement::TENSOR11);
  E[2 - 1][2 - 1] = shape * data;
  LOG_DBG(dm)<< "AHRT 22=" << E[2-1][2-1] << " data=" << data.ToString();
  hom_rect_samples_.GetCol(data,
      DesignElement::TENSOR33 - DesignElement::TENSOR11);
  E[3 - 1][3 - 1] = shape * data;
  LOG_DBG(dm)<< "AHRT 33=" << E[3-1][3-1] << " data=" << data.ToString();
}

void DesignMaterial::ApplyHomRectC1Tensor(Matrix<double>& E, Vector<double>& p,
    DesignElement::Type direction, SubTensorType subTensor) const {
  // Method uses tricubic interpolation
  int m = hom_rect_a_.GetNumRows();
  int n = hom_rect_b_.GetNumRows();
  int o = hom_rect_c_.GetNumRows();
  double da = hom_rect_a_[1][0] - hom_rect_a_[0][0];
  double db = hom_rect_b_[1][0] - hom_rect_b_[0][0];
  double dc = 0.0;
  if (o > 1) {
    dc = hom_rect_c_[1][0] - hom_rect_c_[0][0];
  }
  int j, k, l(-1);
  j = GetInterpolationIndex(hom_rect_a_,p[0]);
  k = GetInterpolationIndex(hom_rect_b_,p[1]);
  if (o > 1) {
    l = GetInterpolationIndex(hom_rect_c_,p[2]);
  }

  if (subTensor == FULL) {
    E.Resize(6, 6);
    E.Init(); // for off-diagonal
    if (direction == DesignElement::NO_DERIVATIVE
        || direction == DesignElement::ROTANGLE || direction == DesignElement::ROTANGLEX || direction == DesignElement::ROTANGLEY) {
      E[1 - 1][1 - 1] = EvaluateC1Interpolation_3D(p, hom_rect_coeff11_, da,
          db, dc, j, k, l, m, n, o);
      E[1 - 1][2 - 1] = EvaluateC1Interpolation_3D(p, hom_rect_coeff12_, da,
          db, dc, j, k, l, m, n, o);
      E[1 - 1][3 - 1] = EvaluateC1Interpolation_3D(p, hom_rect_coeff13_, da,
          db, dc, j, k, l, m, n, o);
      E[2 - 1][3 - 1] = EvaluateC1Interpolation_3D(p, hom_rect_coeff23_, da,
          db, dc, j, k, l, m, n, o);
      E[2 - 1][2 - 1] = EvaluateC1Interpolation_3D(p, hom_rect_coeff22_, da,
          db, dc, j, k, l, m, n, o);
      E[3 - 1][3 - 1] = EvaluateC1Interpolation_3D(p, hom_rect_coeff33_, da,
          db, dc, j, k, l, m, n, o);
      E[4 - 1][4 - 1] = EvaluateC1Interpolation_3D(p, hom_rect_coeff44_, da,
          db, dc, j, k, l, m, n, o);
      E[5 - 1][5 - 1] = EvaluateC1Interpolation_3D(p, hom_rect_coeff55_, da,
          db, dc, j, k, l, m, n, o);
      E[6 - 1][6 - 1] = EvaluateC1Interpolation_3D(p, hom_rect_coeff66_, da,
          db, dc, j, k, l, m, n, o);
      E[2 - 1][1 - 1] = E[1 - 1][2 - 1];
      E[3 - 1][1 - 1] = E[1 - 1][3 - 1];
      E[3 - 1][2 - 1] = E[2 - 1][3 - 1];
      LOG_DBG(dm)<<"E11= "<<E[0][0]<<" E12= "<<E[0][1]<<" E22= "<< E[1][1]<<" E33= "<<E[2][2]<<" E23= "<<E[1][2]<<" E13= "<<E[0][2]<<" E44= "<<E[3][3]<<" E55= "<<E[4][4]<<" E66= "<<E[5][5];
    } else {
      E[1-1][1-1] = EvaluateC1Interpolation_Deriv_3D(p, hom_rect_coeff11_, da,db,dc,j,k,l,m,n,o,direction);
      E[1-1][2-1] = EvaluateC1Interpolation_Deriv_3D(p, hom_rect_coeff12_, da,db,dc,j,k,l,m,n,o,direction);
      E[1-1][3-1] = EvaluateC1Interpolation_Deriv_3D(p, hom_rect_coeff13_, da,db,dc,j,k,l,m,n,o,direction);
      E[2-1][3-1] = EvaluateC1Interpolation_Deriv_3D(p, hom_rect_coeff23_, da,db,dc,j,k,l,m,n,o,direction);
      E[2-1][2-1] = EvaluateC1Interpolation_Deriv_3D(p, hom_rect_coeff22_, da,db,dc,j,k,l,m,n,o,direction);
      E[3-1][3-1] = EvaluateC1Interpolation_Deriv_3D(p, hom_rect_coeff33_, da,db,dc,j,k,l,m,n,o,direction);
      E[4-1][4-1] = EvaluateC1Interpolation_Deriv_3D(p, hom_rect_coeff44_, da,db,dc,j,k,l,m,n,o,direction);
      E[5-1][5-1] = EvaluateC1Interpolation_Deriv_3D(p, hom_rect_coeff55_, da,db,dc,j,k,l,m,n,o,direction);
      E[6-1][6-1] = EvaluateC1Interpolation_Deriv_3D(p, hom_rect_coeff66_, da,db,dc,j,k,l,m,n,o,direction);
      E[2-1][1-1] = E[1-1][2-1];
      E[3-1][1-1] = E[1-1][3-1];
      E[3-1][2-1] = E[2-1][3-1];
      LOG_DBG(dm)<<"Derivative "<<((direction == DesignElement::STIFF1)?"1":((direction == DesignElement::STIFF2) ? "2":"3"))<<" E11= "<<E[0][0]<<" E12= "<<E[0][1]<<" E22= "<< E[1][1]<<" E33= "<<E[2][2]<<" E23= "<<E[1][2]<<" E13= "<<E[0][2]<<" E44= "<<E[3][3]<<" E55= "<<E[4][4]<<" E66= "<<E[5][5];
    }
  } else {
    E.Resize(3,3);
    E.Init(); // for off-diagonal
    if (direction == DesignElement::NO_DERIVATIVE || direction == DesignElement::ROTANGLE || direction == DesignElement::ROTANGLEX || direction == DesignElement::ROTANGLEY) {
      if (o == 0) { // no shearing
        E[1-1][1-1] = EvaluateC1Interpolation(p, hom_rect_coeff11_,da,db,j,k,m,n);
        E[1-1][2-1] = EvaluateC1Interpolation(p, hom_rect_coeff12_,da,db ,j,k,m,n);
        E[2-1][1-1] = E[1-1][2-1];
        E[2-1][2-1] = EvaluateC1Interpolation(p, hom_rect_coeff22_, da,db,j,k,m,n);
        E[3-1][3-1] = EvaluateC1Interpolation(p, hom_rect_coeff33_, da,db,j,k,m,n);
      } else { // shearing
        E[1-1][1-1] = EvaluateC1Interpolation_3D(p, hom_rect_coeff11_,da,db,dc,j,k,l,m,n,o);
        E[1-1][2-1] = EvaluateC1Interpolation_3D(p, hom_rect_coeff12_,da,db,dc,j,k,l,m,n,o);
        E[1-1][3-1] = EvaluateC1Interpolation_3D(p, hom_rect_coeff13_,da,db,dc,j,k,l,m,n,o);
        E[2-1][1-1] = E[1-1][2-1];
        E[2-1][2-1] = EvaluateC1Interpolation_3D(p, hom_rect_coeff22_, da,db,dc,j,k,l,m,n,o);
        E[2-1][3-1] = EvaluateC1Interpolation_3D(p, hom_rect_coeff23_, da,db,dc,j,k,l,m,n,o);
        E[3-1][1-1] = E[1-1][3-1];
        E[3-1][2-1] = E[2-1][3-1];
        E[3-1][3-1] = EvaluateC1Interpolation_3D(p, hom_rect_coeff33_, da,db,dc,j,k,l,m,n,o);
      }
      LOG_DBG(dm)<<p;
      LOG_DBG(dm)<<"E11= "<<E[0][0]<<" E12= "<<E[0][1]<<" E13= "<<E[0][2]<<" E22= "<< E[1][1]<<" E23= "<<E[1][2]<<" E33= "<<E[2][2];
    } else {
      if (o == 0) { // no shearing
        E[1-1][1-1] = EvaluateC1Interpolation_Deriv(p, hom_rect_coeff11_, da,db,j,k,m,n,direction);
        E[1-1][2-1] = EvaluateC1Interpolation_Deriv(p, hom_rect_coeff12_, da,db,j,k,m,n,direction);
        E[2-1][1-1] = E[1-1][2-1];
        E[2-1][2-1] = EvaluateC1Interpolation_Deriv(p, hom_rect_coeff22_, da,db,j,k,m,n,direction);
        E[3-1][3-1] = EvaluateC1Interpolation_Deriv(p, hom_rect_coeff33_, da,db,j,k,m,n,direction);
      } else { // shearing
        E[1-1][1-1] = EvaluateC1Interpolation_Deriv_3D(p, hom_rect_coeff11_, da,db,dc,j,k,l,m,n,o,direction);
        E[1-1][2-1] = EvaluateC1Interpolation_Deriv_3D(p, hom_rect_coeff12_, da,db,dc,j,k,l,m,n,o,direction);
        E[1-1][3-1] = EvaluateC1Interpolation_Deriv_3D(p, hom_rect_coeff13_, da,db,dc,j,k,l,m,n,o,direction);
        E[2-1][1-1] = E[1-1][2-1];
        E[2-1][2-1] = EvaluateC1Interpolation_Deriv_3D(p, hom_rect_coeff22_, da,db,dc,j,k,l,m,n,o,direction);
        E[2-1][3-1] = EvaluateC1Interpolation_Deriv_3D(p, hom_rect_coeff23_, da,db,dc,j,k,l,m,n,o,direction);
        E[3-1][1-1] = E[1-1][3-1];
        E[3-1][2-1] = E[2-1][3-1];
        E[3-1][3-1] = EvaluateC1Interpolation_Deriv_3D(p, hom_rect_coeff33_, da,db,dc,j,k,l,m,n,o,direction);
      }
      LOG_DBG(dm)<<p;
      LOG_DBG(dm)<<"Derivative "<<((direction == DesignElement::STIFF1)?"1":(direction == DesignElement::STIFF2)?"2":"3")
          <<" E11= "<<E[0][0]<<" E12= "<<E[0][1]<<" E13= "<<E[0][2]<<" E22= "<< E[1][1]<<" E23= "<<E[1][2]<<" E33= "<<E[2][2];
    }
  }

}

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
    bool dataIsHierarchised;
    Matrix<double> data;
    dataIsHierarchised = ReadDetailedStats(filename, data);
    // create regular grid
    sg::base::GridGenerator* gridGen = grid_->createGridGenerator();
    gridGen->regular(level_);
    delete gridGen;
    if (dataIsHierarchised) {
      FillSparseGridWithHierarchisedData(data);
    } else {
      FillSparseGridWithUnhierarchisedData(data);
    }
    return;
  }

  // new format
  unsigned int N, d, m;
  file >> N >> d >> m >> word;
  // standard assumptions
  assert((d == 2) || (d == 3));
  assert((m == 4) || (m == 6));
  // logically equivalent to "(shearIsDesign_) ==> ((d=3) and (m=6))"
  assert(!shearIsDesign_ || ((d == 3) && (m == 6)));
  bool hierarchised = (word == "hierarchised");
  file >> word;
  Notation notation = ((word == "voigt") ? VOIGT : HILL_MANDEL);
  // initialize coefficient vectors
  alpha1_ = boost::shared_ptr<sg::base::DataVector>(new sg::base::DataVector(N));
  alpha2_ = boost::shared_ptr<sg::base::DataVector>(new sg::base::DataVector(N));
  alpha3_ = boost::shared_ptr<sg::base::DataVector>(new sg::base::DataVector(N));
  alpha4_ = boost::shared_ptr<sg::base::DataVector>(new sg::base::DataVector(N));
  if (shearIsDesign_) {
    alpha5_ = boost::shared_ptr<sg::base::DataVector>(new sg::base::DataVector(N));
    alpha6_ = boost::shared_ptr<sg::base::DataVector>(new sg::base::DataVector(N));
  }
  std::vector<unsigned int> level(d, 0);
  std::vector<unsigned int> index(d, 0);
  sg::base::GridIndex grid_point(grid_->getStorage()->dim());
  level_ = 0; // we set that to the maximal level of the grid points
  double duck; // dummy
  unsigned int j = 0;
  for (unsigned int i = 0; i < N; i++) {
    if (shearIsDesign_) {
      // shearing angle should be optimized ==> Read ALL The Data!
      file >> level[0] >> level[1] >> level[2];
      file >> index[0] >> index[1] >> index[2];
      grid_point.set(0, level[0], index[0]);
      grid_point.set(1, level[1], index[1]);
      grid_point.set(2, level[2], index[2]);
      level_ = std::max(level_, std::max(level[0], std::max(level[1], level[2])));
    } else {
      // shearing angle should not be optimized ==> maybe we have to pick the right data
      if (d == 2) {
        file >> level[0] >> level[1];
        file >> index[0] >> index[1];
      } else {
        file >> level[0] >> level[1] >> level[2];
        file >> index[0] >> index[1] >> index[2];
        if (level[2] != 1) {
          // data is 3D, but this data point is not in the relevant z=0.5 plane ==> skip
          if (m == 4) {
            // Is there a "non-Ralph" way to do that?
            file >> duck >> duck >> duck >> duck;
          } else {
            file >> duck >> duck >> duck >> duck >> duck >> duck;
          }
          continue;
        }
      }
      grid_point.set(0, level[0], index[0]);
      grid_point.set(1, level[1], index[1]);
      level_ = std::max(level_, std::max(level[0], level[1]));
    }
    // add grid point
    grid_->getStorage()->insert(grid_point);
    if (shearIsDesign_) {
      // shearing angle should be optimized ==> Read ALL The Data!
      file >> (*alpha1_)[j] >> (*alpha2_)[j] >> (*alpha3_)[j] >> (*alpha4_)[j]
           >> (*alpha5_)[j] >> (*alpha6_)[j];
      if (notation == VOIGT) {
        (*alpha6_)[j] *= 2.0;
      }
    } else {
      // shearing angle should not be optimized ==> maybe we have to pick the right data
      if (m == 4) {
        file >> (*alpha1_)[j] >> (*alpha2_)[j] >> (*alpha3_)[j] >> (*alpha4_)[j];
      } else {
        file >> (*alpha1_)[j] >> (*alpha2_)[j] >> duck >> (*alpha3_)[j] >> duck >> (*alpha4_)[j];
      }
      if (notation == VOIGT) {
        (*alpha4_)[j] *= 2.0;
      }
    }
    j++;
  }
  LOG_DBG(dm) << "DM:ISG: level = " << level_ << "\n";
  file.close();
  if (!shearIsDesign_) {
    // coefficient vectors were too big, because we skipped grid points
    alpha1_->resize(grid_->getStorage()->size());
    alpha2_->resize(grid_->getStorage()->size());
    alpha3_->resize(grid_->getStorage()->size());
    alpha4_->resize(grid_->getStorage()->size());
  }
  // hierarchisation if needed
  if (!hierarchised) {
    if ((sgpp_basis_ == LINEAR) || (sgpp_basis_ == MODLINEAR)) {
      sg::base::OperationHierarchisation *hierOp =
          sg::op_factory::createOperationHierarchisation(*grid_);
      hierOp->doHierarchisation(*alpha1_);
      hierOp->doHierarchisation(*alpha2_);
      hierOp->doHierarchisation(*alpha3_);
      hierOp->doHierarchisation(*alpha4_);
      if (shearIsDesign_) {
        hierOp->doHierarchisation(*alpha5_);
        hierOp->doHierarchisation(*alpha6_);
      }
      delete hierOp;
    } else {
      std::vector<sg::base::DataVector *> alphas;
      alphas.push_back(alpha1_.get());
      alphas.push_back(alpha2_.get());
      alphas.push_back(alpha3_.get());
      alphas.push_back(alpha4_.get());
      if (shearIsDesign_) {
        alphas.push_back(alpha5_.get());
        alphas.push_back(alpha6_.get());
      }
      sg::opt::OperationMultipleHierarchisation *hierOp =
        sg::op_factory::createOperationMultipleHierarchisation(*grid_);
      hierOp->doHierarchisation(alphas);
      delete hierOp;
    }
  }
}

void DesignMaterial::FillSparseGridWithUnhierarchisedData(Matrix<double>& data) {
  sg::base::GridStorage* gridStorage = grid_->getStorage();

  // create coefficient vectors
  alpha1_ = boost::shared_ptr<sg::base::DataVector>(new sg::base::DataVector(gridStorage->size()));
  alpha2_ = boost::shared_ptr<sg::base::DataVector>(new sg::base::DataVector(gridStorage->size()));
  alpha3_ = boost::shared_ptr<sg::base::DataVector>(new sg::base::DataVector(gridStorage->size()));
  alpha4_ = boost::shared_ptr<sg::base::DataVector>(new sg::base::DataVector(gridStorage->size()));
  alpha5_ = boost::shared_ptr<sg::base::DataVector>(new sg::base::DataVector(gridStorage->size()));
  alpha6_ = boost::shared_ptr<sg::base::DataVector>(new sg::base::DataVector(gridStorage->size()));

  // put data values in coefficient vectors
  unsigned int dim1, dim2, dim3, index1, index2, index3, row;
  dim1 = catalogueSize_[0];
  dim2 = catalogueSize_[1];
  dim3 = catalogueSize_[2];
  sg::base::GridIndex* gp;
  unsigned int sz = gridStorage->size();
  sz = sz + 1;
  for (unsigned int i=0; i < gridStorage->size(); i++) {
    gp = gridStorage->get(i);
    if (catalogueSize_.GetSize() == 2) {
      index1 = gp->abs(0)*(dim1);
      index2 = gp->abs(1)*(dim2);
      row = (index1-1)*dim2 + index2 - 1;
    } else {
      index1 = gp->abs(0)*(dim1);
      index2 = gp->abs(1)*(dim2);
      if(shearIsDesign_) {
        index3 = gp->abs(2)*(dim3+1);
      } else {
        index3 = .5*(dim3+1);
      }
      row = (index1-1)*dim2*dim3 + (index2-1)*dim3 + index3 - 1;
    }
    alpha1_->set(i,data[row][0]);
    alpha2_->set(i,data[row][1]);
    if (shearIsDesign_) {
      alpha3_->set(i,data[row][2]);
      alpha4_->set(i,data[row][3]);
      alpha5_->set(i,data[row][4]);
      alpha6_->set(i,data[row][5]);
    } else {
      if (catalogueSize_.GetSize() == 2) {
        alpha3_->set(i,data[row][2]);
        alpha4_->set(i,data[row][3]);
      } else {
        alpha3_->set(i,data[row][3]);
        alpha4_->set(i,data[row][5]);
      }
    }
    LOG_DBG3(dm) << gp->abs(0) << " " << gp->abs(1) << " " << gp->abs(2) << " -> "
        << alpha1_->get(i) << " " << alpha2_->get(i) << " " << alpha3_->get(i) << " " << alpha4_->get(i);
  }
  // hierarchise data vectors
  if ((sgpp_basis_ == LINEAR) || (sgpp_basis_ == MODLINEAR)) {
    std::vector<sg::base::DataVector *> alphas;
    alphas.push_back(alpha1_.get());
    alphas.push_back(alpha2_.get());
    alphas.push_back(alpha3_.get());
    alphas.push_back(alpha4_.get());
    if (shearIsDesign_) {
      alphas.push_back(alpha5_.get());
      alphas.push_back(alpha6_.get());
    }
    sg::opt::OperationMultipleHierarchisation *hierOp =
        sg::op_factory::createOperationMultipleHierarchisation(*grid_);
    hierOp->doHierarchisation(alphas);
    delete hierOp;
  } else {
    sg::base::OperationHierarchisation *hierOp =
        sg::op_factory::createOperationHierarchisation(*grid_);
    hierOp->doHierarchisation(*alpha1_);
    hierOp->doHierarchisation(*alpha2_);
    hierOp->doHierarchisation(*alpha3_);
    hierOp->doHierarchisation(*alpha4_);
    if (shearIsDesign_) {
      hierOp->doHierarchisation(*alpha5_);
      hierOp->doHierarchisation(*alpha6_);
    }
    delete hierOp;
  }
}

void DesignMaterial::FillSparseGridWithHierarchisedData(Matrix<double>& data) {
  sg::base::GridStorage* gridStorage = grid_->getStorage();

  // create coefficient vectors
  alpha1_ = boost::shared_ptr<sg::base::DataVector>(new sg::base::DataVector(gridStorage->size()));
  alpha2_ = boost::shared_ptr<sg::base::DataVector>(new sg::base::DataVector(gridStorage->size()));
  alpha3_ = boost::shared_ptr<sg::base::DataVector>(new sg::base::DataVector(gridStorage->size()));
  alpha4_ = boost::shared_ptr<sg::base::DataVector>(new sg::base::DataVector(gridStorage->size()));
  alpha5_ = boost::shared_ptr<sg::base::DataVector>(new sg::base::DataVector(gridStorage->size()));
  alpha6_ = boost::shared_ptr<sg::base::DataVector>(new sg::base::DataVector(gridStorage->size()));

  // put data values in coefficient vectors
  sg::base::GridIndex* gp;
  unsigned int sz = gridStorage->size();
  sz = sz + 1;
  for (unsigned int i=0; i < gridStorage->size(); i++) {
    if (shearIsDesign_) {
      alpha1_->set(i,data[i][0]);
      alpha2_->set(i,data[i][1]);
      alpha3_->set(i,data[i][2]);
      alpha4_->set(i,data[i][3]);
      alpha5_->set(i,data[i][4]);
      alpha6_->set(i,data[i][5]);
    } else {
      alpha1_->set(i,data[i][0]);
      alpha2_->set(i,data[i][1]);
      if (catalogueSize_.GetSize() == 2) {
        alpha3_->set(i,data[i][2]);
        alpha4_->set(i,data[i][3]);
      } else {
        alpha3_->set(i,data[i][3]);
        alpha4_->set(i,data[i][5]);
      }
    }
    gp = gridStorage->get(i);
    LOG_DBG3(dm) << gp->abs(0) << " " << gp->abs(1) << " " << gp->abs(2) << " -> "
        << alpha1_->get(i) << " " << alpha2_->get(i) << " " << alpha3_->get(i) << " " << alpha4_->get(i);
  }

  // hierarchise data vectors
  if ((sgpp_basis_ == LINEAR) || (sgpp_basis_ == MODLINEAR)) {
    sg::base::OperationHierarchisation *hierOp =
        sg::op_factory::createOperationHierarchisation(*grid_);
    hierOp->doHierarchisation(*alpha1_);
    hierOp->doHierarchisation(*alpha2_);
    hierOp->doHierarchisation(*alpha3_);
    hierOp->doHierarchisation(*alpha4_);
    if (shearIsDesign_) {
      hierOp->doHierarchisation(*alpha5_);
      hierOp->doHierarchisation(*alpha6_);
    }
    delete hierOp;
  } else {
    std::vector<sg::base::DataVector *> alphas;
    alphas.push_back(alpha1_.get());
    alphas.push_back(alpha2_.get());
    alphas.push_back(alpha3_.get());
    alphas.push_back(alpha4_.get());
    if (shearIsDesign_) {
      alphas.push_back(alpha5_.get());
      alphas.push_back(alpha6_.get());
    }
    sg::opt::OperationMultipleHierarchisation *hierOp =
        sg::op_factory::createOperationMultipleHierarchisation(*grid_);
    hierOp->doHierarchisation(alphas);
    delete hierOp;
  }
}

void DesignMaterial::ApplyHomRectSGPPTensor(Matrix<double>& E, Vector<double>& p,
     DesignElement::Type direction, SubTensorType subTensor) const {
  // Method uses SGPP interpolation
  sg::base::DataVector point(p.GetSize());
  for (unsigned int i=0; i < p.GetSize(); i++) {
    point[i] = p[i];
  }
  sg::base::OperationEval* opEval = sg::op_factory::createOperationEval(*grid_);
  LOG_DBG(dm) << p;

  E.Resize(3,3);
  E.Init(); // for off-diagonal
  if (direction == DesignElement::NO_DERIVATIVE || direction == DesignElement::ROTANGLE || direction == DesignElement::ROTANGLEX || direction == DesignElement::ROTANGLEY) {
    if (!shearIsDesign_) { // no shearing
      E[1-1][1-1] = std::max(0.0, opEval->eval(*alpha1_, point));
      E[1-1][2-1] = std::max(0.0, opEval->eval(*alpha2_, point));
      E[2-1][1-1] = E[1-1][2-1];
      E[2-1][2-1] = std::max(0.0, opEval->eval(*alpha3_, point));
      E[3-1][3-1] = std::max(0.0, opEval->eval(*alpha4_, point));
    } else { // shearing
      E[1-1][1-1] = std::max(0.0, opEval->eval(*alpha1_, point));
      E[1-1][2-1] = std::max(0.0, opEval->eval(*alpha2_, point));
      E[1-1][3-1] = opEval->eval(*alpha3_, point);
      E[2-1][1-1] = E[1-1][2-1];
      E[2-1][2-1] = std::max(0.0, opEval->eval(*alpha4_, point));
      E[2-1][3-1] = opEval->eval(*alpha5_, point);
      E[3-1][1-1] = E[1-1][3-1];
      E[3-1][2-1] = E[2-1][3-1];
      E[3-1][3-1] = std::max(0.0, opEval->eval(*alpha6_, point));
    }
    LOG_DBG(dm)<<" E11= "<<E[0][0]<<" E12= "<<E[0][1]<<" E13= "<<E[0][2]<<" E22= "<< E[1][1]<<" E23= "<<E[1][2]<<" E33= "<<E[2][2];
  } else {
    if ((sgpp_basis_ == LINEAR) || (sgpp_basis_ == MODLINEAR)) {
      if (!shearIsDesign_) { // no shearing
        E[1-1][1-1] = EvaluateSGPPInterpolation_Deriv(opEval, *alpha1_, point, direction);
        E[1-1][2-1] = EvaluateSGPPInterpolation_Deriv(opEval, *alpha2_, point, direction);
        E[2-1][1-1] = E[1-1][2-1];
        E[2-1][2-1] = EvaluateSGPPInterpolation_Deriv(opEval, *alpha3_, point, direction);
        E[3-1][3-1] = EvaluateSGPPInterpolation_Deriv(opEval, *alpha4_, point, direction);
      } else { // shearing
        E[1-1][1-1] = EvaluateSGPPInterpolation_Deriv(opEval, *alpha1_, point, direction);
        E[1-1][2-1] = EvaluateSGPPInterpolation_Deriv(opEval, *alpha2_, point, direction);
        E[1-1][3-1] = EvaluateSGPPInterpolation_Deriv(opEval, *alpha3_, point, direction);
        E[2-1][1-1] = E[1-1][2-1];
        E[2-1][2-1] = EvaluateSGPPInterpolation_Deriv(opEval, *alpha4_, point, direction);
        E[2-1][3-1] = EvaluateSGPPInterpolation_Deriv(opEval, *alpha5_, point, direction);
        E[3-1][1-1] = E[1-1][3-1];
        E[3-1][2-1] = E[2-1][3-1];
        E[3-1][3-1] = EvaluateSGPPInterpolation_Deriv(opEval, *alpha6_, point, direction);
      }
    } else {
      sg::base::OperationNaiveEvalPartialDerivative* opEvalPartDeriv =
          sg::op_factory::createOperationNaiveEvalPartialDerivative(*grid_);
      if (!shearIsDesign_) { // no shearing
        E[1-1][1-1] = EvaluateSGPPInterpolation_Deriv_Exact(opEvalPartDeriv, *alpha1_, point, direction);
        E[1-1][2-1] = EvaluateSGPPInterpolation_Deriv_Exact(opEvalPartDeriv, *alpha2_, point, direction);
        E[2-1][1-1] = E[1-1][2-1];
        E[2-1][2-1] = EvaluateSGPPInterpolation_Deriv_Exact(opEvalPartDeriv, *alpha3_, point, direction);
        E[3-1][3-1] = EvaluateSGPPInterpolation_Deriv_Exact(opEvalPartDeriv, *alpha4_, point, direction);
      } else { // shearing
        E[1-1][1-1] = EvaluateSGPPInterpolation_Deriv_Exact(opEvalPartDeriv, *alpha1_, point, direction);
        E[1-1][2-1] = EvaluateSGPPInterpolation_Deriv_Exact(opEvalPartDeriv, *alpha2_, point, direction);
        E[1-1][3-1] = EvaluateSGPPInterpolation_Deriv_Exact(opEvalPartDeriv, *alpha3_, point, direction);
        E[2-1][1-1] = E[1-1][2-1];
        E[2-1][2-1] = EvaluateSGPPInterpolation_Deriv_Exact(opEvalPartDeriv, *alpha4_, point, direction);
        E[2-1][3-1] = EvaluateSGPPInterpolation_Deriv_Exact(opEvalPartDeriv, *alpha5_, point, direction);
        E[3-1][1-1] = E[1-1][3-1];
        E[3-1][2-1] = E[2-1][3-1];
        E[3-1][3-1] = EvaluateSGPPInterpolation_Deriv_Exact(opEvalPartDeriv, *alpha6_, point, direction);
      }
      delete opEvalPartDeriv;
    }
    LOG_DBG(dm)<<"Derivative "<<((direction == DesignElement::STIFF1)?"1":(direction == DesignElement::STIFF2)?"2":"3")
        <<" E11= "<<E[0][0]<<" E12= "<<E[0][1]<<" E13= "<<E[0][2]<<" E22= "<< E[1][1]<<" E23= "<<E[1][2]<<" E33= "<<E[2][2];
  }
  delete opEval;
}

void DesignMaterial::ApplyHomRectFullBsplineTensor(Matrix<double>& E, Vector<double>& p,
     DesignElement::Type direction, SubTensorType subTensor) const {
  E.Resize(3,3);
  E.Init(); // for off-diagonal
  const int margin = (bspline_degree_+1)/2 + 1;
  const double a = static_cast<double>((bspline_degree_+1)/2);
  sg::base::BsplineBasis<int, int> bspline(bspline_degree_);

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
        if (direction == DesignElement::NO_DERIVATIVE || direction == DesignElement::ROTANGLE || direction == DesignElement::ROTANGLEX || direction == DesignElement::ROTANGLEY) {
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
    // take positive part (if we're not calculating derivatives)
    if (direction == DesignElement::NO_DERIVATIVE || direction == DesignElement::ROTANGLE || direction == DesignElement::ROTANGLEX || direction == DesignElement::ROTANGLEY) {
      E[1-1][1-1] = std::max(0.0, E[1-1][1-1]);
      E[1-1][2-1] = std::max(0.0, E[1-1][2-1]);
      E[2-1][2-1] = std::max(0.0, E[2-1][2-1]);
      E[3-1][3-1] = std::max(0.0, E[3-1][3-1]);
    }
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
          if (direction == DesignElement::NO_DERIVATIVE || direction == DesignElement::ROTANGLE || direction == DesignElement::ROTANGLEX || direction == DesignElement::ROTANGLEY) {
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
    // take positive part (if we're not calculating derivatives)
    if (direction == DesignElement::NO_DERIVATIVE || direction == DesignElement::ROTANGLE || direction == DesignElement::ROTANGLEX || direction == DesignElement::ROTANGLEY) {
      E[1-1][1-1] = std::max(0.0, E[1-1][1-1]);
      E[1-1][2-1] = std::max(0.0, E[1-1][2-1]);
      E[2-1][2-1] = std::max(0.0, E[2-1][2-1]);
      E[3-1][3-1] = std::max(0.0, E[3-1][3-1]);
    }
    // symmetric entries
    E[2-1][1-1] = E[1-1][2-1];
    E[3-1][1-1] = E[1-1][3-1];
    E[3-1][2-1] = E[2-1][3-1];
  }
}

double DesignMaterial::EvaluateSGPPInterpolation_Deriv(sg::base::OperationEval* opEval,
  sg::base::DataVector alpha, sg::base::DataVector point, DesignElement::Type direction) const {
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
  sg::base::DataVector pointL = point;
  sg::base::DataVector pointU = point;
  pointL[dimension-1] -= h;
  pointU[dimension-1] += h;

  double valL = opEval->eval(alpha,pointL);
  double valU = opEval->eval(alpha,pointU);

  return (valU - valL) / (2*h);
}

inline double DesignMaterial::EvaluateSGPPInterpolation_Deriv_Exact(
  sg::base::OperationNaiveEvalPartialDerivative* opEvalPartDeriv,
  sg::base::DataVector& alpha, sg::base::DataVector& point, DesignElement::Type direction) const {
  unsigned int dimension;
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
  return opEvalPartDeriv->evalPartialDerivative(alpha, point, dimension-1);
}

#endif //USE_SGPP

int DesignMaterial::GetInterpolationIndex(Matrix<double> interval, double& point) const {
  PtrParamNode inf_warn = info->Get("optimization/header/designSpace");
  int nRows = interval.GetNumRows();
  double h = interval[1][0] - interval[0][0];
  int idx = -1;
  if (interval[0][0] <= point && point < interval[nRows - 1][0]) {
    idx = (int) ( (point - interval[0][0]) / h);
  } else if (point == interval[nRows - 1][0]) {
    idx = nRows - 2;
  } else if (point > interval[nRows - 1][0]) {
    idx = nRows - 2;
    point = 1.;
    if (point > 1.01) {
      inf_warn->Get(ParamNode::WARNING)->SetValue(
          "Interpolation of Hom_RectC1 tensor failed. Design Variable "
              + lexical_cast<string>(point) + " out of bounds ");
    }
  } else if (point < interval[0][0]) {
    idx = 0;
    point = interval[0][0];
    if (point < interval[0][0]-0.01) {
      inf_warn->Get(ParamNode::WARNING)->SetValue(
          "Interpolation of Hom_RectC1 tensor failed. Design Variable "
              + lexical_cast<string>(point) + " out of bounds ");
    }
  }
  assert(idx != -1);
  return idx;
}

double DesignMaterial::EvaluateC1Interpolation_3D(Vector<double>& p,
    const Matrix<double> & coeff, double & da, double & db, double & dc,
    int & j, int & k, int & l, int & m, int & n, int &o) const {
  LOG_DBG(dm)<<"p=["<<p[0]<<","<<p[1]<<", "<<p[2]<<"]";
  double t=(p[0]-hom_rect_a_[j][0])/da;
  double u =(p[1]-hom_rect_b_[k][0])/db;
  double v=(p[2]-hom_rect_c_[l][0])/dc;
  LOG_DBG(dm)<<"u = "<<u<<" t= "<<t<<" v= "<<v;
  double res = 0;
  for (int ii = 0;ii<4;ii++) {
    for (int jj=0;jj<4;jj++) {
      for (int kk=0;kk<4;kk++) {
        res += coeff[(n-1)*(o-1)*j+(o-1)*k+l][ii+4*jj+16*kk]*pow(t,ii)*pow(u,jj)*pow(v,kk);
      }
    }
  }
  LOG_DBG(dm) << "Result =" << res;
  return res;
}

double DesignMaterial::EvaluateC1Interpolation_Deriv_3D(Vector<double>& p,
    const Matrix<double> & coeff, double & da, double & db, double & dc,
    int & j, int & k, int & l, int & m, int & n, int & o,
    DesignElement::Type direction) const {
  double u = (p[0] - hom_rect_a_[j][0]) / (da);
  double t = (p[1] - hom_rect_b_[k][0]) / (db);
  double v = (p[2] - hom_rect_c_[l][0]) / (dc);
  LOG_DBG(dm)<<"Deriv: u = "<<u<<" t= "<<t<<" v= "<<v<<" j= "<<j<<" k= "<<k<<" l= "<<l;
  LOG_DBG(dm)<<"p_deriv: ["<<p[0]<<", "<<", "<<p[1]<<", "<<p[2];
  double deriv = 0;
  if (direction == DesignElement::STIFF1) {
    for (int ii = 1; ii < 4; ii++) {
      for (int jj = 0; jj < 4; jj++) {
        for (int kk = 0; kk < 4; kk++) {
          deriv += coeff[(n - 1) * (o - 1) * j + (o - 1) * k + l][ii + 4 * jj
              + 16 * kk] * ii * pow(u, ii - 1) * pow(t, jj) * pow(v, kk);
        }
      }
    }
    deriv /= da;
  }
  if (direction == DesignElement::STIFF2) {
    for (int ii = 0; ii < 4; ii++) {
      for (int jj = 1; jj < 4; jj++) {
        for (int kk = 0; kk < 4; kk++) {
          deriv += coeff[(n - 1) * (o - 1) * j + (o - 1) * k + l][ii + 4 * jj
              + 16 * kk] * jj * pow(u, ii) * pow(t, jj - 1) * pow(v, kk);
        }
      }
    }
    deriv /= db;
  }
  if (direction == DesignElement::STIFF3) {
    for (int ii = 0; ii < 4; ii++) {
      for (int jj = 0; jj < 4; jj++) {
        for (int kk = 1; kk < 4; kk++) {
          deriv += coeff[(n - 1) * (o - 1) * j + (o - 1) * k + l][ii + 4 * jj
              + 16 * kk] * kk * pow(u, ii) * pow(t, jj) * pow(v, kk - 1);
        }
      }
    }
    deriv /= dc;
  }
  LOG_DBG(dm)<< "Deriv Result =" << deriv;
  return deriv;
}

double DesignMaterial::EvaluateC1Interpolation(Vector<double>& p,
    const Matrix<double> & coeff, double & da, double & db, int & j, int & k,
    int & m, int & n) const {
  LOG_DBG(dm)<<"p=["<<p[0]<<","<<p[1]<<"]";
  double u =(p[1]-hom_rect_b_[k][0])/(db);
  double t=(p[0]-hom_rect_a_[j][0])/(da);
  LOG_DBG(dm)<<"u = "<<u<<" t= "<<t<<"\n";
  LOG_DBG(dm)<<"u = "<<u<<" t= "<<t;
  LOG_DBG(dm)<<"j = "<<j<<" k= "<<k;
  double res = 0;
  for (int i = 0;i<4;i++) {
    for (int l=0;l<4;l++) {
      res += coeff[(n-1)*j+k][(i)*4+l]*pow(t,i)*pow(u,l);
    }
  }
  LOG_DBG(dm) << "Result =" << res;
  return res;
}

double DesignMaterial::EvaluateC1Interpolation_Deriv(Vector<double>& p,
    const Matrix<double> & coeff, double & da, double & db, int & j, int & k,
    int & m, int & n, DesignElement::Type direction) const {
  double u = (p[1] - hom_rect_b_[k][0]) / (db);
  double t = (p[0] - hom_rect_a_[j][0]) / (da);
  LOG_DBG(dm)<<"Deriv: u = "<<u<<" t= "<<t<<"\n";

  double deriv = 0;
  if (direction == DesignElement::STIFF1) {
    for (int i = 1; i < 4; i++) {
      for (int l = 0; l < 4; l++) {
        deriv += coeff[(n - 1) * j + k][(i) * 4 + l] * i * pow(t, i - 1)
            * pow(u, l);
      }
    }
    deriv /= da;
  }
  if (direction == DesignElement::STIFF2) {
    for (int i = 0; i < 4; i++) {
      for (int l = 1; l < 4; l++) {
        deriv += coeff[(n - 1) * j + k][(i) * 4 + l] * l * pow(t, i)
            * pow(u, l - 1);
      }
    }
    deriv /= db;
  }
  LOG_DBG(dm)<< "Deriv Result =" << deriv;
  return deriv;
}

bool DesignMaterial::ReadDetailedStats(const char * filename, Matrix<double>& ret) {
  bool isHierarchised;
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
  Notation notation = !strVec[3].compare(0,strVec[3].size(),"voigt") ? VOIGT : HILL_MANDEL;

  // initialize matrix
  // if second entry starts with "L", we have hierarchised data
  if (strVec[1].compare(0,1,"L") == 0) {
    isHierarchised = true;
    unsigned int catalogueDimension = strVec[0].compare("3D") == 0 ? 3 : 2;
    level_ = boost::lexical_cast<int>(strVec[1].substr(1));
    for (unsigned int i=0; i<catalogueDimension; i++) {
      catalogueSize_.Push_back(2^level_);
    }
    nRows = boost::lexical_cast<int>(strVec[2]);
  } else {
    isHierarchised = false;
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

  if (shearIsDesign_) {
    assert(catalogueSize_.GetSize() == 3);
  }

  if (catalogueSize_.GetSize() == 2) {
    nCols = 4;
  } else {
    nCols = 6;
  }
  ret.Resize(nRows,nCols);
  ret.Init();
  LOG_DBG(dm)<<"ReadDetailedStats: nRows = "<<nRows<<", nCols = "<<nCols<<", level = "<<level_;

  for (unsigned int i=0; i<nRows; i++) {
    // read line and write data into matrix
    getline(file,line,'\n');
    strVec.Clear();
    SplitStringListWhitespace(line, strVec);
    for (unsigned int j=0; j<nCols; j++) {
      ret[i][j] = boost::lexical_cast<double>(strVec[j+catalogueSize_.GetSize()]);
    }
    if (notation == VOIGT) {
      ret[i][nCols-1] = ret[i][nCols-1] * 2.0;
    }
  }
  LOG_DBG3(dm)<<"Data: \n"<<ret;

  return isHierarchised;
}

void DesignMaterial::GetInterpolatedTensor(Matrix<double>& t,
    SubTensorType subTensor, DesignElement::Type direction, Notation notation){
  double a = (direction == DesignElement::INTERPOLATION) ? 1.0 : params_[DesignElement::INTERPOLATION];
  double ma = (direction == DesignElement::INTERPOLATION) ? -1.0 : 1.0-a;
  double E = params_[DesignElement::EMODUL];
  switch (subTensor) {
  case FULL:
    t.Resize(6, 6);
    t.Init();
//    SetOrthotropicTensor(t, subTensor, a*255.68181818+ma*294.03409091, a*99.43181818+ma*80.0, a*99.43181818+ma*114.34659091,
//        a*255.68181818+ma*166.19318182, a*99.43181818+ma*80.0, a*255.68181818+ma*294.03409091, a*78.125+ma*70.0, a*78.125+ma*70.0, a*78.125+ma*60.0);
//    SetOrthotropicTensor(t, subTensor, a*1.91761363636+ma*1.54548810664, a*0.745738636364+ma*0.529863106641,
//        a*0.745738636364+ma*0.622605363985, a*1.91761363636+ma*1.54548810664, a*0.745738636364+ma*0.622605363985,
//        a*1.91761363636+ma*2.87356321839, a*0.5859375+ma*0.57, a*0.5859375+ma*0.57, a*0.5859375+ma*0.5078125);
    SetOrthotropicTensor(t, subTensor, a*1.91761363636+ma*(1.54548810664-E/5), a*0.745738636364+ma*0.529863106641,
        a*0.745738636364+ma*0.622605363985, a*1.91761363636+ma*(1.54548810664-E/5), a*0.745738636364+ma*0.622605363985,
        a*1.91761363636+ma*(2.87356321839+E), a*0.5859375+ma*(0.57-E/10), a*0.5859375+ma*(0.57-E/10), a*0.5859375+ma*(0.5078125-E/10));
    break;
  case PLANE_STRAIN:
  case PLANE_STRESS:
    transIsoType_ = TRANSISO_XZ;
    SetTransIsoTensor(t, subTensor, a*1.91761363636+ma*1.54548810664, 0.0, 0.0, a*1.91761363636+ma*2.87356321839, a*0.745738636364+ma*0.529863106641, a*0.5859375+ma*0.5078125);
    break;
  default:
    throw Exception("subTensor not implemented yet");
  }

  if (type_ == D_INTERP_TENSOR || type_ == D_INTERP_TENSOR_ROT)
  {
    double dens = params_[DesignElement::DENSITY];
    if (direction == DesignElement::DENSITY)
    {
      if(penalty_ == 1.0)
        dens = 1.0;
      else
        dens = penalty_*std::pow(dens, penalty_-1.0);
    }
    else
    {
      dens = std::pow(dens, penalty_);
    }
    t *= dens;
  }
  if(type_ == D_INTERP_TENSOR_ROT){
    // for all rotated types, rotate the material tensor
    LOG_DBG3(dm) << "GetTransIsoMaterialTensor: tensor before rotation=" << t.ToString();
    RotateVoigtTensor(t, direction);
  }
}

void DesignMaterial::GetLaminatesTensor(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction, Notation notation){
  switch(subTensor){
  case PLANE_STRAIN:    //see Allaire: Shape optimization by the homogenization method, pp. 127 [(2.64),(2.65)]
  {
    t.Resize(3,3);
    t.Init();
    double eps = 5.0625e-4;
    double stiff1 = params_[DesignElement::STIFF1];
    double stiff2 = params_[DesignElement::STIFF2];
    double E = params_[DesignElement::EMODUL];
    double nu = params_[DesignElement::POISSON];
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
      ZeroTensor(t, subTensor);
      return;
    }
    break;
  }
  case PLANE_STRESS:      //see Bendsoe, Sigmund: Topology Optimization S. 166
  {
    double E33 = 0.02;
    double stiff1 = params_[DesignElement::STIFF1];
    double stiff2 = params_[DesignElement::STIFF2];
    double E = params_[DesignElement::EMODUL];
    double nu = params_[DesignElement::POISSON];
    double n = (stiff2 + stiff1 * stiff2 * (nu * nu - 1) - 1);
    switch (direction) {
    case DesignElement::NO_DERIVATIVE:
    case DesignElement::ROTANGLE:
    case DesignElement::DENSITY:
    {
      double E11 = -(E*stiff1)/n;
      double E22 = stiff2*E+stiff2*stiff2*nu*nu*E11;
      double E12 = stiff2*nu*E11;
      Set2dVoigtTensor(t, E11, E22, E33, 0.0, 0.0, E12);
      break;
    }
    case DesignElement::STIFF1: {
      double E11 = -(E * (stiff2 - 1)) / (n * n);
      double E22 = stiff2 * stiff2 * nu * nu * E11;
      double E12 = stiff2 * nu * E11;
      Set2dVoigtTensor(t, E11, E22, 0.0, 0.0, 0.0, E12);
      break;
    }
    case DesignElement::STIFF2: {
      double E11 = (E * stiff1 * (stiff1 * (nu * nu - 1) + 1)) / (n * n);
      double E22 = E - 2 * stiff2 * nu * nu * E * stiff1 / n
          + stiff2 * stiff2 * nu * nu * E11;
      double E12 = -nu * E * stiff1 / n + stiff2 * nu * E11;
      Set2dVoigtTensor(t, E11, E22, 0.0, 0.0, 0.0, E12);
    break;
    }
    default:
      ZeroTensor(t, subTensor);
      return;
    }
    break;
  }
  default:
    throw Exception("subTensor not implemented yet");
  }

  if (type_ == D_LAMINATES)
  {
    double dens = params_[DesignElement::DENSITY];
    if (direction == DesignElement::DENSITY)
    {
      if(penalty_ == 1.0)
        dens = 1.0;
      else
        dens = penalty_*std::pow(dens, penalty_-1.0);
    }
    else
    {
      dens = std::pow(dens, penalty_);
    }
    t *= dens;
  }

  double rotAngle = params_[DesignElement::ROTANGLE];
  RotateHMStiffnessTensor(t, subTensor, direction, rotAngle, notation);
  return;
}

void DesignMaterial::ZeroTensor(Matrix<double>& t, SubTensorType subTensor) {
  switch (subTensor) {
  case FULL:
    t.Resize(6, 6);
    LOG_DBG(dm)<<"Zero Tensor: "<<t.ToString(2);
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

void DesignMaterial::Set2dVoigtTensor(Matrix<double>& t, double t11, double t22,
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

void DesignMaterial::SetOrthotropicTensor(Matrix<double>& t,
    SubTensorType subTensor, double e11, double e12, double e13, double e22,
    double e23, double e33, double e44, double e55, double e66) {
  switch (subTensor) {
    case FULL:
      t.Resize(6, 6);
      t.Init();
      t[0][0] = e11;
      t[0][1] = e12;
      t[0][2] = e13;
      t[1][0] = e12;
      t[1][1] = e22;
      t[1][2] = e23;
      t[2][0] = e13;
      t[2][1] = e23;
      t[2][2] = e33;
      t[3][3] = e44;
      t[4][4] = e55;
      t[5][5] = e66;
      break;
    case PLANE_STRAIN:
    case PLANE_STRESS:
      SetTransIsoTensor(t, subTensor, e11, 0.0, 0.0, e22, e12, e66);
      break;
    default:
      throw Exception("subTensor not implemented yet");
    }
}

void DesignMaterial::SetTransIsoTensor(Matrix<double>& t,
    SubTensorType subTensor, double iD, double inD, double iG, double oD,
    double onD, double oG) {
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

void DesignMaterial::SetIsoTensor(Matrix<double>& t, SubTensorType subTensor,
    double D, double nd, double G) {
  SetTransIsoTensor(t, subTensor, D, nd, G, D, nd, G);
}


void DesignMaterial::RotateHMStiffnessTensor(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction, double a, Notation notation)
{
  // FIXME, DOCU needed
  // FIXME, rotation of Hill-Mandel-Tensors might be not a good idea
  // FIXME, this is specific to two dimensions!
  switch(subTensor){
  case PLANE_STRAIN:
  case PLANE_STRESS:
  case PLANE:
  {
    Matrix<double> theta(3,3);
    Matrix<double> help(3,3);
    const double sq2inv = 1/sqrt(2);
    theta.SetEntry(0,0, pow(cos(a),2));
    theta.SetEntry(0,1, pow(sin(a),2));
    theta.SetEntry(0,2, -sqrt(2)/2*sin(2*a));
    theta.SetEntry(1,0, theta(0,1));
    theta.SetEntry(1,1, theta(0,0));
    theta.SetEntry(1,2, -theta(0,2));
    theta.SetEntry(2,0, theta(1,2));
    theta.SetEntry(2,1, theta(0,2));
    theta.SetEntry(2,2, cos(2*a));
    t.Mult(theta, help);
    if(direction == DesignElement::ROTANGLE){
      Matrix<double> dtheta(3,3);
      dtheta.SetEntry(0,0, -sin(2*a));
      dtheta.SetEntry(0,1, -dtheta(0,0));
      dtheta.SetEntry(0,2, -sqrt(2)*cos(2*a));
      dtheta.SetEntry(1,0, dtheta(0,1));
      dtheta.SetEntry(1,1, dtheta(0,0));
      dtheta.SetEntry(1,2, -dtheta(0,2));
      dtheta.SetEntry(2,0, dtheta(1,2));
      dtheta.SetEntry(2,1, dtheta(0,2));
      dtheta.SetEntry(2,2, -2*sin(2*a));
      Matrix<double> dthetaTttheta(3,3);

      LOG_DBG3(dm) << "RHMST phi=" << a << " Theta=" << theta.ToString(2) << " dTheta=" << dtheta.ToString(2);

      dtheta.MultT(help, dthetaTttheta);
      t.Mult(dtheta, help);
      theta.MultT(help, dtheta);
      t = dthetaTttheta + dtheta;
      if(notation == VOIGT)
      {
        t(0,2)*=sq2inv;
        t(1,2)*=sq2inv;
        t(2,2)/=2;
        t(2,0)*=sq2inv;
        t(2,1)*=sq2inv;
      }
      return;
    }
    theta.MultT(help, t);
    if(notation == VOIGT)
    {
      t(0,2)*=sq2inv;
      t(1,2)*=sq2inv;
      t(2,2)/=2;
      t(2,0)*=sq2inv;
      t(2,1)*=sq2inv;
    }
    LOG_DBG3(dm) << "RHMST phi=" << a << " Theta=" << theta.ToString(2);

    return;
  }
  default:
    throw Exception("subTensor not implemented yet");
  }
}

void DesignMaterial::RotateVoigtTensor(Matrix<double>& t, DesignElement::Type direction){
  // rotation matrix is found in Dissertation of B. Schmidt: Topology Preserving Multi-Layer Shape and Material Optimization p. 62
  // and also found in Wikipedia Drehmatrix (german)
  // rotates the material by thetaz around the z-axis by thetay around the y-axis and by thetax around the x-axis in this given order
  // direction of rotation around an axis is positive (ccw), if the axis is pointing towards oneself
  // this is identical to BaseMaterial::RotateTensorByRotationAngles

  double thetax = 0.0, thetay = 0.0, thetaz = 0.0;
  if(dim == 3){
    thetax = params_[DesignElement::ROTANGLEX];
    thetay = params_[DesignElement::ROTANGLEY];
    thetaz = 0.0;
    if(params_.find(DesignElement::ROTANGLEZ) != params_.end()){
      thetaz = params_[DesignElement::ROTANGLEZ];
    }else if(transIsoType_ != TRANSISO_XY){
      throw Exception("The parameterization of the rotation is incompatible to the choice of the isotropic plane");
    }
  }else{ // dim == 2
    thetaz = params_[DesignElement::ROTANGLE];
  }

  double sthetax = sin(thetax);
  double cthetax = cos(thetax);
  double sthetay = sin(thetay);
  double cthetay = cos(thetay);
  double sthetaz = sin(thetaz);
  double cthetaz = cos(thetaz);
  Matrix<Double> R(dim, dim);
  SetRotationMatrix(R, sthetax, cthetax, sthetay, cthetay, sthetaz, cthetaz);
  LOG_DBG3(dm) << "Rotation matrix for tx=" << thetax << ", ty=" << thetay << ", tz=" << thetaz << " is " << R.ToString();
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
  LOG_DBG3(dm) << "Corresponding Q is " << Q.ToString();
  if(direction != DesignElement::ROTANGLEX && direction != DesignElement::ROTANGLEY && direction != DesignElement::ROTANGLEZ && direction != DesignElement::ROTANGLE){
    // calculate Q*t*Q' and store back to t. unfortunately MultT is the wrong way
    Matrix<Double> help(dimQ, dimQ);
    Q.Mult(t, help);
    Matrix<Double> QT(dimQ, dimQ);
    QT.Resize(dimQ, dimQ);
    Q.Transpose(QT);
    help.Mult(QT, t);
  }else{ // we need a derivative
    Matrix<Double> dR(dim, dim);
    SetRotationMatrix(dR, sthetax, cthetax, sthetay, cthetay, sthetaz, cthetaz, direction); // this now produces the derivative
    LOG_DBG3(dm) << "Corresponding dR is " << dR.ToString();

    Matrix<Double> dQ(dimQ, dimQ);
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

    if(dim == 3){
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
    LOG_DBG3(dm) << "Corresponding dQ is " << dQ.ToString();

    // we now, have to calculate dQ*t*Q' + Q*t*dQ'
    Matrix<Double> help(dimQ, dimQ);
    Matrix<Double> dQT(dimQ, dimQ);
    dQT.Resize(dimQ, dimQ);
    dQ.Transpose(dQT);
    dQ.Mult(t, help);
    Matrix<Double> QT(dimQ, dimQ);
    QT.Resize(dimQ, dimQ);
    Q.Transpose(QT);
    help.Mult(QT, dQ); // dQ is no longer needed, we overwrite it
    Q.Mult(t, help);
    help.Mult(dQT, t); // here, we overwrite t
    t.Add(1.0, dQ);    // and add the rest
    //FIXME: this section is ugly and should be fixed if expression templates work reliably
  }

}
void DesignMaterial::SetRotationMatrix(Matrix<double>& R, double sthetax, double cthetax, double sthetay, double cthetay, double sthetaz, double cthetaz, DesignElement::Type direction){
  // rotation matrix is found in Dissertation of B. Schmidt: Topology Preserving Multi-Layer Shape and Material Optimization p. 62
  // and also found in Wikipedia Drehmatrix (german)
  // rotates the material by thetaz around the z-axis by thetay around the y-axis and by thetax around the x-axis in this given order
  // direction of rotation around an axis is positive (ccw), if the axis is pointing towards oneself
  R.Resize(dim, dim);
  double ndx = 1.0, ndy = 1.0, ndz = 1.0;
  // for the derivative, we replace the correspond sin by cos and cos by -sin and set nd to 0.0
  switch(direction){
  case DesignElement::ROTANGLEX:
    ndx = sthetax;
    sthetax = cthetax;
    cthetax = -ndx;
    ndx = 0.0;
    break;
  case DesignElement::ROTANGLEY:
    ndy = sthetay;
    sthetay = cthetay;
    cthetay = -ndy;
    ndy = 0.0;
    break;
  case DesignElement::ROTANGLEZ:
  case DesignElement::ROTANGLE:
    ndz = sthetaz;
    sthetaz = cthetaz;
    cthetaz = -ndz;
    ndz = 0.0;
    break;
  default:
    break;
  }
  R[0][0] =  ndx * cthetay * cthetaz;
  R[0][1] = -ndx * cthetay * sthetaz;
  R[1][0] =  cthetax * ndy * sthetaz + sthetax * sthetay * cthetaz;
  R[1][1] =  cthetax * ndy * cthetaz - sthetax * sthetay * sthetaz;
  if(dim == 3){
    R[0][2] =  ndx * sthetay * ndz;
    R[1][2] = -sthetax * cthetay * ndz;
    R[2][0] =  sthetax * ndy * sthetaz - cthetax * sthetay * cthetaz;
    R[2][1] =  sthetax * ndy * cthetaz + cthetax * sthetay * sthetaz;
    R[2][2] =  cthetax * cthetay * ndz;
  }
}

void DesignMaterial::RotatePiezoCouplingTensor(Matrix<double>& E, double phi, DesignElement::Type direction)
{
  // R(phi) * [e] * Q(phi)^T
  // derivative: dR(phi)/dphi * ([e] * Q(phi)^T) + R(phi) * ([e] * dQ(phi)/dphi)^T

  // Note, that we use VOIGT rotation matrix Q here, while in RotateHMStiffnessTensor Hill-Mandel is used. Also the QT here is QT^T of the HM rotation!

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

    LOG_DBG3(dm) << "RPCT phi=" << phi << " R=" << R.ToString(2) << " QT=" << QT.ToString(2);
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
    LOG_DBG3(dm) << "RPCT phi=" << phi << " R=" << R.ToString(2) << " dR=" << dR.ToString(2) << " QT=" << QT.ToString(2) << " dQT=" << dQT.ToString(2);
    return;
  }
}


void DesignMaterial::RotateElecTensor(Matrix<double>& E, double phi, DesignElement::Type direction)
{
  // R(phi) * [e] * R(phi)^T
  // derivative: dR(phi)/dphi * ([e] * R(phi)^T) + R(phi) * ([e] * dR(phi)/dphi)^T

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


double DesignMaterial::GetTransIsoMass(double iD, double iG, double oD, double oG){
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

double DesignMaterial::GetIsoMass(double D, double G) {
  return (GetTransIsoMass(D, G, D, G));
}

void DesignMaterial::GetMaterialTensor(Matrix<double>& t, SubTensorType subTensor, DesignElement::Type direction, Notation notation) {
  assert(!(notation == HILL_MANDEL && type_ != FMO && type_ != LAMINATES && type_ != D_LAMINATES && type_ != HOM_RECT && type_ != D_HOM_RECT && type_ != HOM_RECT_C1 && type_ !=  DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC && type_ != DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED && type_ != ORTHOTROPIC && type_ != DENSITY_TIMES_ROT_PA12 && type_ != REDBAS_PARAM && type_ != REDBAS_FREE && type_ != GREEDY_PARAM && type_ != GREEDY_FREE && type_ != GREEDY_MAPPING));

  switch (type_) {
  case FMO:
    GetElasticFMOTensor(t, direction, notation);
    break;
  case ORTHOTROPIC:
  case DENSITY_TIMES_ORTHOTROPIC:
    GetOrthotropicMaterialTensor(t, subTensor, direction, notation);
    break;
  case ISOTROPIC:
    GetIsoMaterialTensor(t, subTensor, direction);
    break;
  case LAME_ISOTROPIC: // LAME_ISOTROPIC
    GetLameMaterialTensor(t, subTensor, direction);
    break;
  case TRANSVERSAL_ISOTROPIC:
  case TRANSVERSAL_ISOTROPIC_BOXED:
  case DENSITY_TIMES_TRANSVERSAL_ISOTROPIC:
  case DENSITY_TIMES_TRANSVERSAL_ISOTROPIC_BOXED:
  case DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC:
  case DENSITY_TIMES_ROT_TRANSVERSAL_ISOTROPIC_BOXED:
  case DENSITY_TIMES_ROT_PA12:
    GetTransIsoMaterialTensor(t, subTensor, direction, notation);
    break;
  case DENSITY_TIMES_2D_TENSOR:
  case DENSITY_TIMES_2D_TENSOR_CONSTANT_TRACE:
  case DENSITY_TIMES_ROTATED_2D_TENSOR:
    GetDensityTimes2dTensorTensor(t, subTensor, direction);
    break;
  case LAMINATES:
  case D_LAMINATES:
    GetLaminatesTensor(t, subTensor, direction, notation);
    break;
  case HOM_RECT:
  case D_HOM_RECT:
  case HOM_RECT_C1:
    GetHomRectTensor(t, subTensor, direction, notation);
    break;
  case GREEDY_PARAM:
  case GREEDY_FREE:
  case REDBAS_PARAM:
  case REDBAS_FREE:
    GetModRedTensor(t, direction, notation);
    break;
  case GREEDY_MAPPING:
  case REDBAS_MAPPING:
    GetMappingTensor(t, direction, notation);
    break;
  case D_INTERP_TENSOR:
  case D_INTERP_TENSOR_ROT:
    GetInterpolatedTensor(t, subTensor, direction, notation);
    break;
  default: // case default
    throw Exception("DesignMaterial Type not implemented yet");
  }
}


void DesignMaterial::GetElecTensor(Matrix<double>& E, DesignElement::Type direction)
{

  // only 2D!
  bool set = direction == DesignElement::NO_DERIVATIVE || direction == DesignElement::ROTANGLE;

  double e11 = set ? params_[DesignElement::DIELEC_11] : 0;
  double e22 = set ? params_[DesignElement::DIELEC_22] : 0;
  double e12 = set ? params_[DesignElement::DIELEC_12] : 0;
  double rotAngle = set ? params_[DesignElement::ROTANGLE] : 0;

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

  LOG_DBG2(dm) << "GET: E before rotation = " << E.ToString(2);
  RotateElecTensor(E, rotAngle, direction);
  LOG_DBG2(dm) << "GET: E after rotation =  " << E.ToString(2);
}


void DesignMaterial::GetPiezoCouplingTensor(Matrix<double>& E, DesignElement::Type direction)
{
  // only 2D!
  bool set = direction == DesignElement::NO_DERIVATIVE || direction == DesignElement::ROTANGLE;
  double e11 = set ? params_[DesignElement::PIEZO_11] : 0;
  double e12 = set ? params_[DesignElement::PIEZO_12] : 0;
  double e13 = set ? params_[DesignElement::PIEZO_13] : 0;
  double e21 = set ? params_[DesignElement::PIEZO_21] : 0;
  double e22 = set ? params_[DesignElement::PIEZO_22] : 0;
  double e23 = set ? params_[DesignElement::PIEZO_23] : 0;
  double rotAngle = set ? params_[DesignElement::ROTANGLE] : 0;
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

  LOG_DBG2(dm) << "GPCT: E before rotation = " << E.ToString(2) << " ra=" << rotAngle << " d=" << DesignElement::type.ToString(direction);
  RotatePiezoCouplingTensor(E, rotAngle, direction);
  LOG_DBG2(dm) << "GPCT: E after rotation =  " << E.ToString(2);

}

double DesignMaterial::GetMaterialMass(DesignElement::Type direction) {
  if (massIsDesign_) {
    switch (direction) {
    case DesignElement::MASS:
      return (massFactor_);
    case DesignElement::NO_DERIVATIVE:
      return (params_[DesignElement::MASS] * massFactor_);
    default:
      return (0.0);
    }
  } else {
    switch (type_) {
    case ISOTROPIC:
      return (GetIsoMaterialMass(direction) * massFactor_);
    case LAME_ISOTROPIC: // LAME_ISOTROPIC
      return (GetLameMaterialMass(direction) * massFactor_);
    case TRANSVERSAL_ISOTROPIC:
    case TRANSVERSAL_ISOTROPIC_BOXED:
      return (GetTransIsoMaterialMass(direction) * massFactor_);
    case DENSITY_TIMES_ROT_PA12:
      return (GetDensityTimesTensorMass(direction) * massFactor_);
    default: // case default
      throw Exception("DesignMaterial Type not implemented yet");
    }
  }
}

bool DesignMaterial::GetMaterialDamping(double& alpha, double& beta,
    DesignElement::Type direction) {
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
      alpha = params_[DesignElement::DAMPINGALPHA];
      beta = params_[DesignElement::DAMPINGBETA];
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

void DesignMaterial::DumpParams() {
  std::map<DesignElement::Type, double>::iterator iter;

  for (iter = params_.begin(); iter != params_.end(); ++iter)
    std::cout << "params[" << DesignElement::type.ToString(iter->first)
        << "] = " << iter->second << std::endl;
}

void DesignMaterial::SetEnums() {
  type.SetName("DesignMaterial::Type");
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
  type.Add(REDBAS_PARAM, "reducedBasis-param");
  type.Add(REDBAS_FREE, "reducedBasis-free");
  type.Add(GREEDY_PARAM, "greedy-param");
  type.Add(GREEDY_FREE, "greedy-free");
  type.Add(GREEDY_MAPPING, "greedy-mapping");
  type.Add(D_HOM_RECT, "density-times-hom-rect");
  type.Add(HOM_RECT_C1, "hom-rect-C1");
  type.Add(D_INTERP_TENSOR, "density-times-interpolated-tensor");
  type.Add(D_INTERP_TENSOR_ROT, "density-times-rotated-interpolated-tensor");
  type.Add(REDBAS_MAPPING, "reducedBasis-mapping");
  transIsoType.SetName("DesignMaterial::TransIsoType");
  transIsoType.Add(TRANSISO_XY, "xy");
  transIsoType.Add(TRANSISO_YZ, "yz");
  transIsoType.Add(TRANSISO_XZ, "xz");

  notation.SetName("DesignMaterial::Notation");
  notation.Add(VOIGT, "voigt");
  notation.Add(HILL_MANDEL, "hill_mandel");
  notation.Add(HILL_MANDEL_NO_DENSITY, "hill_mandel_no_density");
}

void DesignMaterial::GetSVDGTensorParameters(const Matrix<double>& G, Vector<double>& paramvec)
{
   /* G= R(theta)*(l1 0 \\ 0l2) R(phi) paramvec = (theta, phi, l1, l2) */
  //Check first that G is a 2 by 2 matrix
  assert((G.GetNumRows() == 2) & (G.GetNumCols() == 2));

  if((G.GetNumRows() != 2) || (G.GetNumCols() != 2)){
      throw Exception("GetSVDGTensor only works for 2*2 matrices");
  }


  paramvec.Resize(4);
  paramvec.Init();

  //help = G'*G
  //help2 = G*G'
  Matrix<double> help(2,2);
  Matrix<double> transpose(2,2);
  G.Transpose(transpose);
  Matrix<double> help2(2,2);
  G.MultT(G,help);
  G.Mult(transpose,help2);

  double theta = 0.5*std::atan2(help2(0,1)+help2(1,0), help2(0,0)-help2(1,1));

  double phi = -0.5*std::atan2(help(0,1)+help(1,0), help(0,0)-help(1,1));

  Matrix<double> Rthetainv(2,2);
  Matrix<double> Rphiinv(2,2);
  Rthetainv(0,0) = cos(theta);
  Rthetainv(1,1) = cos(theta);
  Rthetainv(0,1) = -sin(theta);
  Rthetainv(1,0) = sin(theta);

  Rphiinv(0,0) = cos(phi);
  Rphiinv(1,1) = cos(phi);
  Rphiinv(0,1) = -sin(phi);
  Rphiinv(1,0) = sin(phi);

  Matrix<double> t1(2,2);
  Matrix<double> t2(2,2);

  G.Mult(Rphiinv, t1);
  Rthetainv.Mult(t1,t2);

  double l1 = t2(0,0);
  double l2 = t2(1,1);

  if (l1 <0)
  {
    l1 = -l1;
    l2 = -l2;
    if (phi< 0) phi = PI +phi;
    else phi = phi - PI;
  }

  if (phi <0) phi = phi +2*PI;
  if (theta<0 ) theta = theta + 2*PI;

  paramvec[0] = theta;
  paramvec[1] = phi;
  paramvec[2] = l1;
  paramvec[3] = l2;

}

double DesignMaterial::AngleGreedyCalculus(const Vector<double>& coeffs, const double& angle)
{
   assert(coeffs.size() == 2*Na_+1);

   double s=0;
   for (UInt k= 0; k<2*Na_+1; k++ )
   {
     int indk = k - Na_;
     if (indk<0)
     {
       s += coeffs[k]*sin(-indk*angle);
     }
     else if (indk>0)
     {
       s += coeffs[k]*cos(indk*angle);
     }
     else
     {
       s += coeffs[k];
     }
   }
   return s;

}

double DesignMaterial::ScalingGreedyCalculus(const Vector<double>& coeffs, const double& l)
{
   assert(coeffs.size() == Nl_+1);
   assert((l>=lmin_) & (l<=lmax_));
   double dl = (lmax_ -lmin_)/Nl_;

   int ind = static_cast <int> (std::floor((l-lmin_)/dl));

   if (ind < 0) ind=0;
   if (ind > int(Nl_)-1) ind = Nl_-1;
   double r = l-lmin_-ind*dl;

   return (1.0-r)*coeffs[ind] + r*coeffs[ind+1];

}
