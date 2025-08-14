#include <numeric>
#include <unordered_set>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/XmlReader.hh"
#include "DataInOut/ProgramOptions.hh"
#include "MatVec/Matrix.hh"
#include "MatVec/Vector.hh"
#include "Optimization/Design/SplineBoxDesign.hh"
#include "Optimization/Excitation.hh"
#include "Utils/tools.hh"

namespace CoupledField {

DEFINE_LOG(SBD, "SplineBoxDesign")

Enum<SplineBoxDesign::Interpolation> SplineBoxDesign::interpolation;
Enum<SplineBoxDesign::AnalyticFunc> SplineBoxDesign::analyticFunc;

SplineBoxDesign::SplineBoxDesign(StdVector<RegionIdType>& regionIds, PtrParamNode pn, ErsatzMaterial::Method method)
: FeaturedDesign(regionIds, pn, method)
{
  setup_timer_->Start();

  interpolation.SetName("SplineBoxDesign::Interpolation");
  interpolation.Add(NONE, "none");
  interpolation.Add(CUBIC, "cubic");

  analyticFunc.SetName("SplineBoxDesign::AnalyticFunc");
  analyticFunc.Add(FILE, "File");
  analyticFunc.Add(SUM_OF_SINE, "SumOfSine");
  analyticFunc.Add(MAX_SINE, "MaxSine");
  analyticFunc.Add(SINE_X, "SineX");
  analyticFunc.Add(SINE_Y, "SineY");
  analyticFunc.Add(SINE_Z, "SineZ");

  this->degree_ = pn->Get("splineBox/degree")->As<unsigned int>();

  this->num_cp_.Resize(3);
  this->num_cp_[0] = pn->Get("splineBox/cp_x")->As<unsigned int>();
  this->num_cp_[1] = pn->Get("splineBox/cp_y")->As<unsigned int>();
  this->num_cp_[2] = this->dim_ == 3 ? pn->Get("splineBox/cp_z")->As<unsigned int>() : 1;
  this->total_num_cp_ = num_cp_.Product();
  LOG_DBG(SBD) << "SplineBoxDesign: number of cp = " << this->num_cp_.ToString()
      << ", total number of cp = " << total_num_cp_;

  bounding_box_.Resize(this->dim_, 2);
  if(!pn->Has("splineBox/boundingBox")) {
    bounding_box_ = domain->GetGrid()->CalcGridBoundingBox();
  } else {
    PtrParamNode pt = pn->Get("splineBox/boundingBox");
    bounding_box_[0][0] = pt->Get("lowerXLimit")->As<double>();
    bounding_box_[0][1] = pt->Get("upperXLimit")->As<double>();
    bounding_box_[1][0] = pt->Get("lowerYLimit")->As<double>();
    bounding_box_[1][1] = pt->Get("upperYLimit")->As<double>();
    if(this->dim_ == 3) {
      bounding_box_[2][0] = pt->Get("lowerZLimit")->As<double>();
      bounding_box_[2][1] = pt->Get("upperZLimit")->As<double>();
    }
  }

  // Generate bsplines and get greville abscissae
  bsplines_.Resize(this->dim_);
  StdVector<StdVector<double>> GAbscissae(this->dim_);
  for(unsigned int i = 0; i < this->dim_; ++i) {
    BSpline* bspline = new BSpline(degree_, num_cp_[i]+degree_+1, bounding_box_[i][0], bounding_box_[i][1]);
    bsplines_[i] = bspline;
    GAbscissae[i] = bspline->GrevilleAbscissae();
  }

  // set initial control_points to greville abscissae
  this->initial_control_points_.Reserve(this->total_num_cp_);
  for(unsigned int i = 0; i < this->total_num_cp_; ++i) {
    StdVector<int> sub(this->dim_);
    Ind2Sub(this->num_cp_, i, sub);
    Point p;
    for(unsigned int j = 0; j < this->dim_; ++j) {
      p[j] = GAbscissae[j][sub[j]];
    }
    this->initial_control_points_.Push_back(p);
  }
  SetControlPoints(initial_control_points_, false);

  if(pn->Has("splineBox/feature")) {
    forward_ = false;
    if(pn->Get("splineBox/feature/type")->As<string>() == "file") {
      this->analyticFunc_ = AnalyticFunc::FILE;
      if(!pn->Has("splineBox/feature/file")) {
        throw Exception("Feature density file has to be given");
      } else {
        cover_box_.Resize(this->dim_, 2);
        if(!pn->Has("splineBox/feature/file/coverBox")) {
          cover_box_[0][0] = 0.0;
          cover_box_[0][1] = 1.0;
          cover_box_[1][0] = 0.0;
          cover_box_[1][1] = 1.0;
          if(this->dim_ == 3) {
            cover_box_[2][0] = 0.0;
            cover_box_[2][1] = 1.0;
          }
        } else {
          PtrParamNode pt = pn->Get("splineBox/feature/file/coverBox");
          cover_box_[0][0] = pt->Get("lowerXLimit")->As<double>();
          cover_box_[0][1] = pt->Get("upperXLimit")->As<double>();
          cover_box_[1][0] = pt->Get("lowerYLimit")->As<double>();
          cover_box_[1][1] = pt->Get("upperYLimit")->As<double>();
          if(this->dim_ == 3) {
            cover_box_[2][0] = pt->Get("lowerZLimit")->As<double>();
            cover_box_[2][1] = pt->Get("upperZLimit")->As<double>();
          }
        }

        string set = pn->Get("splineBox/feature/file/set")->As<string>();
        ReadFeature(pn->Get("splineBox/feature/file/path")->As<string>(), set);

        this->periodic_ = pn->Get("splineBox/feature/file/periodic")->As<bool>();
        this->feature_scale_ = pn->Get("splineBox/feature/file/scale")->As<double>();

        string interpol = pn->Get("splineBox/feature/file/interpolation/type")->As<string>();
        this->interpolation_ = interpolation.Parse(interpol);
        this->beta_ = pn->Get("splineBox/feature/file/interpolation/beta")->As<double>();
        assert(this->beta_ < 500);
        LOG_DBG(SBD) << "SBD: interpolation = " << interpol << ", beta = " << this->beta_;
        InterpolateFeature();
      }
    } else {
      // type is "analytic"
      string func = pn->Get("splineBox/feature/function/name")->As<string>();
      this->analyticFunc_ = analyticFunc.Parse(func);
      this->feature_scale_ = pn->Get("splineBox/feature/function/scale")->As<double>();
      this->periodic_ = true;
      this->beta_ = 100;
    }
  } else {
    forward_ = true;
  }

  // this are coordinates of nodes in FE mesh
  StdVector<Vector<double>> points = GetPointsForBasis();

  GenerateBasis(points);

  SetupDesign(pn->Get("splineBox"));

  SetupOptParam();

  this->numInt_.Init(this, pn->Get("splineBox"), info_->Get("splineBox/numInt"));

  StdVector<Point> offset(total_num_cp_);
  for(unsigned int i = 0; i < offset.GetSize(); ++i) {
    for(unsigned int j = 0; j < dim_; ++j) {
      offset[i][j] = shape_param_[i*dim_ + j]->GetDesign(BaseDesignElement::Access::PLAIN);
    }
  }
  SetControlPoints(offset, true);
  LOG_DBG(SBD) << "SplineBoxDesign: initial control points: ";
  for(unsigned int i = 0; i < this->total_num_cp_; ++i) {
    LOG_DBG(SBD) << this->control_points_[i].GetCoordVector().ToString();
  }

  if(forward_) {
    Matrix<double> new_coords;
    EvalAll(new_coords);
    UpdateFEMesh(new_coords); // changes FE mesh
  } else {
    MapFeatureToDensity();
  }

  setup_timer_->Stop();
}

/* <splinebox degree="3" cpx="4" cpy="5" cpz="3" fixed_boundary="true">
     <boundingBox lowerXLimit="0" upperXLimit="1" .../>   optional
     <allControlPoints lower="0" upper="1"/>              oder
     <cp index="0" dof="X|Y|Z" lower="0" upper="1" initial="0.5"/>
   </splinebox> */
void SplineBoxDesign::SetupDesign(PtrParamNode pn)
{
  // shape, i.e. distortion of control points
  unsigned int nshapeparams = total_num_cp_ * dim_;
  double l = -std::numeric_limits<double>::infinity();
  double u =  std::numeric_limits<double>::infinity();
  double v = 0.0;

  if(pn->Has("allControlPoints")){
    l = pn->Get("allControlPoints")->Get("lower")->As<double>();
    u = pn->Get("allControlPoints")->Get("upper")->As<double>();
  }

  /** This are the design variables, which is distortion and rotation.
   *  distortion variables are offset to initial control_points
   *  (i.e. deformation). Contains offset for all control_points.
   *  Last three variables are rotation variables for X, Y and Z axis.
   *  Rotation will always be first Z, then Y, then X.
   *  In 2D the angles for Y and X are 0. */
  shape_param_.Reserve(nshapeparams);

  for(unsigned int i = 0; i < nshapeparams; ++i) {
    ShapeParamElement* de = new ShapeParamElement(BaseDesignElement::Type::CP, shape_param_.GetSize());
    de->SetLowerBound(l);
    de->SetUpperBound(u);
    de->SetDesign(v);
    shape_param_.Push_back(de);
  }

  ParamNodeList cp = pn->GetList("controlpoint");
  if(!cp.IsEmpty()) {
    unsigned int ind;
    for(unsigned int i = 0; i < cp.GetSize(); ++i) {
      unsigned int index = cp[i]->Get("index")->As<unsigned int>();
      string dof = cp[i]->Get("dof")->As<string>();
      if(cp[i]->Has("lower")) {
        l = cp[i]->Get("lower")->As<double>();
      }
      if(cp[i]->Has("upper")) {
        u = cp[i]->Get("upper")->As<double>();
      }
      if(cp[i]->Has("initial")) {
        v = cp[i]->Get("initial")->As<double>();
      }

      if(dof == "x") {
        ind = index*dim_;
      } else if(dof == "y") {
        ind = index*dim_ + 1;
      } else {
        ind = index*dim_ + 2;
      }

      ShapeParamElement* de = shape_param_[ind];
      de->SetLowerBound(l);
      de->SetUpperBound(u);
      de->SetDesign(v);
      LOG_DBG3(SBD) << "SD: shape_param_[" << ind << "]= " << v;
    }
  }

  if(pn->Has("fixedBoundary")) {
    fixed_boundary_ = pn->Get("fixedBoundary")->As<bool>();
  }

  SetupMapping();
}

void SplineBoxDesign::SetupOptParam()
{
  unsigned int nshapeparams = total_num_cp_ * dim_;
  assert(shape_param_.GetSize() == nshapeparams);

  is_opt_.Resize(shape_param_.GetSize(), false);
  opt_shape_param_.Reserve(shape_param_.GetSize());
  // distortion
  for(unsigned int i = 0; i < nshapeparams; ++i) {
    StdVector<int> sub(dim_);
    // intentional integer division to map dof -> cp_index -> cp_sub
    Ind2Sub(num_cp_, i/dim_, sub);
    LOG_DBG3(SBD) << "SOP: i= " << i << " sub= " << sub.ToString() << " i/dim= " << i/dim_;
    bool isOnBoundary = sub[0] == 0 || sub[1] == 0 || (dim_ == 3 ? sub[2] == 0 : false )
        || sub[0] == (int)num_cp_[0]-1 || sub[1] == (int)num_cp_[1]-1 || ((dim_ == 3) ? sub[2] == (int)num_cp_[2]-1 : false);
    if(!(fixed_boundary_ && isOnBoundary)) {
      opt_shape_param_.Push_back(shape_param_[i]);
      is_opt_[i] = true;
    }
  }

  opt_shape_param_.Trim();

  for(unsigned int i = 0, n = opt_shape_param_.GetSize(); i < n; i++) {
    opt_shape_param_[i]->SetOptIndex(i);
  }

  LOG_DBG2(SBD)<< "SOP: design variables =" << opt_shape_param_.GetSize() << ", variables =" << shape_param_.GetSize() << ", feature design variables =" << data.GetSize();
}

void SplineBoxDesign::PostInit(int objectives, int constraints)
{
  FeaturedDesign::PostInit(objectives, constraints);

  setup_timer_->Start();

  WriteBoxToVTK();

  setup_timer_->Stop();
}

void SplineBoxDesign::CheckPlausibility()
{
  FeaturedDesign::CheckPlausibility();

  if(num_cp_[0] < degree_)
    throw Exception("Number of control points has to be at least degree+1");
  if(num_cp_[1] < degree_)
    throw Exception("Number of control points has to be at least degree+1");
  if(dim_ == 3 && num_cp_[2] < degree_)
    throw Exception("Number of control points has to be at least degree+1");
}

int SplineBoxDesign::ReadDesignFromExtern(const double* space_in, bool setAndWriteCurrent)
{
  int old_design = design_id;

  FeaturedDesign::ReadDesignFromExtern(space_in, setAndWriteCurrent);

  if(design_id != old_design) {
    StdVector<Point> offset(total_num_cp_);
    unsigned int k = 0;
    for(unsigned int i = 0; i < offset.GetSize(); ++i) {
      for(unsigned int j = 0; j < dim_; ++j) {
        offset[i][j] = is_opt_[i*dim_ + j] ? opt_shape_param_[k++]->GetDesign(BaseDesignElement::Access::PLAIN) : 0.0;
      }
    }
    SetControlPoints(offset, true);

    if(forward_) {
      Matrix<double> new_coords;
      EvalAll(new_coords);
      UpdateFEMesh(new_coords); // changes FE mesh
    }
  }

  if(!forward_ && mapped_design_ != design_id)
    MapFeatureToDensity();

  WriteBoxToVTK();

  LOG_DBG(SBD) << "RDFE: design_id -> " << design_id;
  return design_id;
}

void SplineBoxDesign::EvalAll(Matrix<double>& out)
{
  Matrix<double> box_mtx = GetValuesOfBasisFunctions();

  // "Casting" StdVector<Point> to Matrix<double>
  Matrix<double> cp;
  cp.Resize(total_num_cp_, dim_);
  cp.Init();
  for(unsigned int i = 0; i < total_num_cp_; ++i) {
    Vector<double> p = control_points_[i].GetCoordVector();
    for(unsigned int j = 0; j < dim_; ++j) {
      cp[i][j] = p[j];
    }
  }
//  LOG_DBG3(SBD) << "EA: control points: \n" << cp.ToString(0);

  // matrix * curr_control_points = deformed nodes
  // dimensions: num_nodes_inside x total_num_cp_ * total_num_cp_ x dim
  out.Resize(box_mtx.GetNumRows(), dim_);
  box_mtx.Mult(cp, out);
//  LOG_DBG3(SBD) << "EA: deformed nodes: \n" << out.ToString(0);
}

Vector<double> SplineBoxDesign::Eval(Vector<double> point)
{
  StdVector<int> sub(dim_);

  Matrix<double> cp;
  cp.Resize(total_num_cp_, dim_);
  cp.Init();
  for(unsigned int i = 0; i < total_num_cp_; ++i) {
    Vector<double> p = control_points_[i].GetCoordVector();
    for(unsigned int j = 0; j < dim_; ++j) {
      cp[i][j] = p[j];
    }
  }

  Matrix<double> new_coords;
  // matrix * curr_control_points = deformed nodes
  // dimensions: 1 x total_num_cp_ * total_num_cp_ x dim
  new_coords.Resize(1, dim_);

  Matrix<double> mtx = GetValuesOfBasisFunctions(point);
  mtx.Mult(cp, new_coords);
  LOG_DBG3(SBD) << "Eval: basis=" << mtx.ToString() << "\ncp=" << cp.ToString() << "\nnew_coords=" << new_coords.ToString();

  // convert Matrix<double> of size 1 x dim to Vector<double>
  Vector<double> out(dim_);
  new_coords.GetRow(out, 0);
  LOG_DBG3(SBD) << "Eval: point= " << point.ToString() << " out= " << out.ToString();
  return out;
}

StdVector<Vector<double>> SplineBoxDesign::EvalDerivative(Vector<double> point)
{
  StdVector<Vector<double>> out(dim_);

  Matrix<double> mtx = GetValuesOfBasisFunctions(point);

  for(unsigned int d = 0; d < dim_; ++d) {
    mtx.GetRow(out[d],0);
  }
  return out;
}


void SplineBoxDesign::SetupVirtualShapeElementMap(Function* f, StdVector<Function::Local::Identifier>& vem, Function::Local::Locality locality)
{
  // FIXME Todo copy from SetupVirtualMultiShapeElementMap
  assert(f != NULL);
  assert(f->IsLocal(f->GetType()));
  // we shall be called by Local::PostInit() therefore local shall exist
  assert(f->GetLocal() != NULL);

  assert(locality == Function::Local::NEXT || locality == Function::Local::PREV_NEXT_AND_REVERSE || locality == Function::Local::PREV_NEXT || locality == Function::Local::NEXT_AND_REVERSE || locality == Function::Local::PREV_NEXT_AND_REVERSE);

  // a lot copy&paste from Function::SetupVirtualElementMap()
  bool prev = locality == Function::Local::PREV_NEXT_AND_REVERSE || locality == Function::Local::PREV_NEXT;
  // next is always true!
  bool two_signs = locality == Function::Local::NEXT_AND_REVERSE || locality == Function::Local::PREV_NEXT_AND_REVERSE;

  assert(f->GetType() == Function::CONES);
  assert(f->GetDesignType() == DesignElement::SPLINE_BOX);

  // we don't set Function::Local::element_dimension_
  // we wont't use the full space as the individual shape_ are not connected
  unsigned int nRows = 0.0;
  Vector<unsigned int> ncp(dim_);
  for(unsigned int d = 0; d < dim_; ++d) {
    ncp = num_cp_;
    --ncp[d];
    nRows += ncp.Product();
  }
  nRows *= std::pow(2, dim_-1);
  vem.Reserve(nRows);

  for(int e = 0; e < (int)total_num_cp_; ++e) {
    int psz = shape_param_.GetSize();

    ShapeParamElement* bde = shape_param_[e*dim_];
    assert(f->GetDesignType() == bde->GetType());
    assert(bde->dof_ == ShapeParamElement::Dof::X);

    int prev_idx = (e-1)*dim_;
    int next_idx = (e+1)*dim_;
    BaseDesignElement* prev_de = prev ? prev_idx < 0 ? NULL : shape_param_[prev_idx] : NULL;
    BaseDesignElement* next_de =        next_idx > psz ? NULL : shape_param_[next_idx];

    LOG_DBG3(SBD) << "SVSEM po=" << (prev_de != NULL ? (int) prev_de->GetOptIndex() : -1) << " eo=" << e << " no=" << (next_de != NULL ? (int) next_de->GetOptIndex() : -1)
                       << " p="  << (prev_de != NULL ? (int) prev_de->GetIndex() : -1)    << " e=" << bde->GetIndex() << " n=" << (next_de != NULL ? (int) next_de->GetIndex() : -1);

    // we add constraints only if control point is subject to optimization
    if(!(is_opt_[e] || is_opt_[e+1] || is_opt_[e+2]) && !(1))
      continue;

    vem.Push_back(Function::Local::Identifier(bde, prev_de, next_de, 1));
    vem.Push_back(Function::Local::Identifier(bde, prev_de, next_de, 2));
    vem.Push_back(Function::Local::Identifier(bde, prev_de, next_de, 3));
    vem.Push_back(Function::Local::Identifier(bde, prev_de, next_de, 4));
  }

  LOG_DBG(SBD) << "SVSEM final f=" << f->ToString() << " loc=" << locality << " ts=" << two_signs << " prev=" << prev << " -> vem=" << vem.GetSize();
}

void SplineBoxDesign::SetupVirtualMultiShapeElementMap(Function* f, StdVector<Function::Local::Identifier>& vem, Function::Local::Locality locality)
{
  assert(f != NULL);
  assert(f->IsLocal(f->GetType()));
  // we shall be called by Local::PostInit() therefore local shall exist
  assert(f->GetLocal() != NULL);

  assert(locality == Function::Local::MULT_DESIGNS_PREV_NEXT_AND_REVERSE || locality == Function::Local::MULT_DESIGNS_PREV_NEXT || locality == Function::Local::MULT_DESIGNS_NEXT_AND_REVERSE || locality == Function::Local::MULT_DESIGNS_NEXT);

  // a lot copy&paste from SetupVirtualShapeElementMap()
  bool prev = locality == Function::Local::MULT_DESIGNS_PREV_NEXT_AND_REVERSE || locality == Function::Local::MULT_DESIGNS_PREV_NEXT;

  assert(f->GetType() == Function::CONES);
  assert(f->GetDesignType() == DesignElement::SPLINE_BOX || f->GetDesignType() == DesignElement::CP);

  // we don't set Function::Local::element_dimension_
  // we wont't use the full space as the individual shape_ are not connected
  unsigned int nRows = 0.0;
  Vector<unsigned int> ncp(dim_);
  for(unsigned int d = 0; d < dim_; ++d) {
    ncp = num_cp_;
    --ncp[d];
    nRows += ncp.Product();
  }
  nRows *= std::pow(2, dim_-1);
  vem.Reserve(nRows);

  StdVector<BaseDesignElement*> opt_buddies; // to be reused temporary vector
  StdVector<BaseDesignElement*> all_buddies; // to be reused temporary vector
  int nsp = total_num_cp_ * dim_;
  BaseDesignElement* bde = nullptr;
  BaseDesignElement *bdex = nullptr, *prev_dex = nullptr, *next_dex = nullptr;
  BaseDesignElement *bdey = nullptr, *prev_dey = nullptr, *next_dey = nullptr;
  BaseDesignElement *bdez = nullptr, *prev_dez = nullptr, *next_dez = nullptr;

  // traverse control points only and do the corresponding coordinates implicitly
  for(int e = 0; e < (int)total_num_cp_; ++e) {
    // we add a constraint only if the involved control point is subject to optimization
    if(is_opt_[e*dim_] || is_opt_[e*dim_+1] || (dim_ == 3 ? is_opt_[e*dim_+2] : false)) {
      // neighbors in x direction
      int curr_idx =  e     *dim_;
      int prev_idx = (e - 1)*dim_;
      int next_idx = (e + 1)*dim_;
      // In x direction we check, if the control points are in the same row.
      // Imagine the following simple two by two cp grid:
      //   2 -- 3
      //   |    |
      //   0 -- 1
      // Without the row check, the predecessor of "2" would be "1".
      // With the check we have curr_row = 1 != 0 = prev_row.
      int curr_row = curr_idx / (num_cp_[0] * dim_); // intentional integer division
      int prev_row = prev_idx / (num_cp_[0] * dim_); // intentional integer division
      int next_row = next_idx / (num_cp_[0] * dim_); // intentional integer division
      bdex = shape_param_[curr_idx];
      assert(bdex->GetType() == DesignElement::CP);
      prev_dex = prev ? (prev_idx > -1 && prev_row == curr_row ? shape_param_[prev_idx] : NULL) : NULL;
      next_dex =        next_idx < nsp && next_row == curr_row ? shape_param_[next_idx] : NULL;

      curr_idx =  e     *dim_+1;
      prev_idx = (e - 1)*dim_+1;
      next_idx = (e + 1)*dim_+1;
      bdey = shape_param_[curr_idx];
      prev_dey = prev ? (prev_idx > -1 && prev_row == curr_row ? shape_param_[prev_idx] : NULL) : NULL;
      next_dey =        next_idx < nsp && next_row == curr_row ? shape_param_[next_idx] : NULL;

      if(dim_ == 3) {
        curr_idx =  e     *dim_+2;
        prev_idx = (e - 1)*dim_+2;
        next_idx = (e + 1)*dim_+2;
        bdez = shape_param_[curr_idx];
        prev_dez = prev ? (prev_idx > -1 && prev_row == curr_row ? shape_param_[prev_idx] : NULL) : NULL;
        next_dez =        next_idx < nsp && next_row == curr_row ? shape_param_[next_idx] : NULL;
      }


      if(prev_dex) {
        opt_buddies.Clear(true);
        if(is_opt_[bdex->GetIndex()]) {
          bde = bdex;
          if(is_opt_[bdey->GetIndex()])
            opt_buddies.Push_back(bdey);
          if(dim_ == 3 && is_opt_[bdez->GetIndex()])
            opt_buddies.Push_back(bdez);
        } else if(is_opt_[bdey->GetIndex()]) {
          bde = bdey;
          if(dim_ == 3 && is_opt_[bdez->GetIndex()])
            opt_buddies.Push_back(bdez);
        } else {
          bde = bdez;
        }
        if(is_opt_[prev_dex->GetIndex()])
          opt_buddies.Push_back(prev_dex);
        if(is_opt_[prev_dey->GetIndex()])
          opt_buddies.Push_back(prev_dey);
        if(dim_ == 3 && is_opt_[prev_dez->GetIndex()])
          opt_buddies.Push_back(prev_dez);

        all_buddies.Clear(true);
        all_buddies.Push_back(bdex);
        all_buddies.Push_back(bdey);
        if(dim_ == 3)
          all_buddies.Push_back(bdez);
        all_buddies.Push_back(prev_dex);
        all_buddies.Push_back(prev_dey);
        if(dim_ == 3)
          all_buddies.Push_back(prev_dez);
        // 2D: cases 1 and 2 describe cone in x direction
        vem.Push_back(Function::Local::Identifier(bde, opt_buddies, all_buddies, -1));
        vem.Push_back(Function::Local::Identifier(bde, opt_buddies, all_buddies, -2));
        if(dim_ == 3) {
          // 3D: cases 1,2,3 and 4 describe cone in x direction
          vem.Push_back(Function::Local::Identifier(bde, opt_buddies, all_buddies, -3));
          vem.Push_back(Function::Local::Identifier(bde, opt_buddies, all_buddies, -4));
        }
        LOG_DBG3(SBD) << "SVMSEM: x opt_buddies=";
        for(auto b : opt_buddies)
          LOG_DBG3(SBD) << b->GetIndex();
        LOG_DBG3(SBD) << "SVMSEM: x all_buddies=";
        for(auto b : all_buddies)
          LOG_DBG3(SBD) << b->GetIndex();
      }

      if(next_dex) {
        opt_buddies.Clear(true);
        if(is_opt_[bdex->GetIndex()]) {
          bde = bdex;
          if(is_opt_[bdey->GetIndex()])
            opt_buddies.Push_back(bdey);
          if(dim_ == 3 && is_opt_[bdez->GetIndex()])
            opt_buddies.Push_back(bdez);
        } else if(is_opt_[bdey->GetIndex()]) {
          bde = bdey;
          if(dim_ == 3 && is_opt_[bdez->GetIndex()])
            opt_buddies.Push_back(bdez);
        } else {
          bde = bdez;
        }
        if(is_opt_[next_dex->GetIndex()])
          opt_buddies.Push_back(next_dex);
        if(is_opt_[next_dey->GetIndex()])
          opt_buddies.Push_back(next_dey);
        if(dim_ == 3 && is_opt_[next_dez->GetIndex()])
          opt_buddies.Push_back(next_dez);

        all_buddies.Clear(true);
        all_buddies.Push_back(bdex);
        all_buddies.Push_back(bdey);
        if(dim_ == 3)
          all_buddies.Push_back(bdez);
        all_buddies.Push_back(next_dex);
        all_buddies.Push_back(next_dey);
        if(dim_ == 3)
          all_buddies.Push_back(next_dez);
        // 2D: cases 1 and 2 describe cone in x direction
        vem.Push_back(Function::Local::Identifier(bde, opt_buddies, all_buddies, 1));
        vem.Push_back(Function::Local::Identifier(bde, opt_buddies, all_buddies, 2));
        if(dim_ == 3) {
          // 3D: cases 1,2,3 and 4 describe cone in x direction
          vem.Push_back(Function::Local::Identifier(bde, opt_buddies, all_buddies, 3));
          vem.Push_back(Function::Local::Identifier(bde, opt_buddies, all_buddies, 4));
        }
        LOG_DBG3(SBD) << "SVMSEM: x opt_buddies=";
        for(auto b : opt_buddies)
          LOG_DBG3(SBD) << b->GetIndex();
        LOG_DBG3(SBD) << "SVMSEM: x all_buddies=";
        for(auto b : all_buddies)
          LOG_DBG3(SBD) << b->GetIndex();
      }

      LOG_DBG3(SBD) << "SVMSEM: x p="  << (prev_dex != NULL ? (int) prev_dex->GetIndex() : -1)
          << " e=" << bde->GetIndex() << " n=" << (next_dex != NULL ? (int) next_dex->GetIndex() : -1);

      // neighbors in y direction
      curr_idx =  e              *dim_;
      prev_idx = (e - num_cp_[0])*dim_;
      next_idx = (e + num_cp_[0])*dim_;
      // See explanation above. The check "prev_row == curr_row" matters only in 3D
      // and will always return true in 2D. Those 'rows' are slices in the x-y-plane.
      curr_row = curr_idx / (num_cp_[0] * num_cp_[1] * dim_); // intentional integer division
      prev_row = prev_idx / (num_cp_[0] * num_cp_[1] * dim_); // intentional integer division
      next_row = next_idx / (num_cp_[0] * num_cp_[1] * dim_); // intentional integer division
      bdex = shape_param_[curr_idx];
      prev_dex = prev ? (prev_idx > -1 && prev_row == curr_row ? shape_param_[prev_idx] : NULL) : NULL;
      next_dex =        next_idx < nsp && next_row == curr_row ? shape_param_[next_idx] : NULL;

      curr_idx =  e              *dim_+1;
      prev_idx = (e - num_cp_[0])*dim_+1;
      next_idx = (e + num_cp_[0])*dim_+1;
      bdey = shape_param_[curr_idx];
      prev_dey = prev ? (prev_idx > -1 && prev_row == curr_row ? shape_param_[prev_idx] : NULL) : NULL;
      next_dey =        next_idx < nsp && next_row == curr_row ? shape_param_[next_idx] : NULL;

      if(dim_ == 3) {
        curr_idx =  e              *dim_+2;
        prev_idx = (e - num_cp_[0])*dim_+2;
        next_idx = (e + num_cp_[0])*dim_+2;
        bdez = shape_param_[curr_idx];
        prev_dez = prev ? (prev_idx > -1 && prev_row == curr_row ? shape_param_[prev_idx] : NULL) : NULL;
        next_dez =        next_idx < nsp && next_row == curr_row ? shape_param_[next_idx] : NULL;
      }


      if(prev_dex) {
        opt_buddies.Clear(true);
        if(is_opt_[bdex->GetIndex()]) {
          bde = bdex;
          if(is_opt_[bdey->GetIndex()])
            opt_buddies.Push_back(bdey);
          if(dim_ == 3 && is_opt_[bdez->GetIndex()])
            opt_buddies.Push_back(bdez);
        } else if(is_opt_[bdey->GetIndex()]) {
          bde = bdey;
          if(dim_ == 3 && is_opt_[bdez->GetIndex()])
            opt_buddies.Push_back(bdez);
        } else {
          bde = bdez;
        }
        if(is_opt_[prev_dex->GetIndex()])
          opt_buddies.Push_back(prev_dex);
        if(is_opt_[prev_dey->GetIndex()])
          opt_buddies.Push_back(prev_dey);
        if(dim_ == 3 && is_opt_[prev_dez->GetIndex()])
          opt_buddies.Push_back(prev_dez);

        all_buddies.Clear(true);
        all_buddies.Push_back(bdex);
        all_buddies.Push_back(bdey);
        if(dim_ == 3)
          all_buddies.Push_back(bdez);
        all_buddies.Push_back(prev_dex);
        all_buddies.Push_back(prev_dey);
        if(dim_ == 3)
          all_buddies.Push_back(prev_dez);
        // 2D: cases 3 and 4 describe cone in y direction
        if(dim_ == 2) {
          vem.Push_back(Function::Local::Identifier(bde, opt_buddies, all_buddies, -3));
          vem.Push_back(Function::Local::Identifier(bde, opt_buddies, all_buddies, -4));
        } else {
          // 3D: cases 5,6,7 and 8 describe cone in y direction
          vem.Push_back(Function::Local::Identifier(bde, opt_buddies, all_buddies, -5));
          vem.Push_back(Function::Local::Identifier(bde, opt_buddies, all_buddies, -6));
          vem.Push_back(Function::Local::Identifier(bde, opt_buddies, all_buddies, -7));
          vem.Push_back(Function::Local::Identifier(bde, opt_buddies, all_buddies, -8));
        }
        LOG_DBG3(SBD) << "SVMSEM: y opt_buddies=";
        for(auto b : opt_buddies)
          LOG_DBG3(SBD) << b->GetIndex();
        LOG_DBG3(SBD) << "SVMSEM: y all_buddies=";
        for(auto b : all_buddies)
          LOG_DBG3(SBD) << b->GetIndex();
      }

      if(next_dex) {
        opt_buddies.Clear(true);
        if(is_opt_[bdex->GetIndex()]) {
          bde = bdex;
          if(is_opt_[bdey->GetIndex()])
            opt_buddies.Push_back(bdey);
          if(dim_ == 3 && is_opt_[bdez->GetIndex()])
            opt_buddies.Push_back(bdez);
        } else if(is_opt_[bdey->GetIndex()]) {
          bde = bdey;
          if(dim_ == 3 && is_opt_[bdez->GetIndex()])
            opt_buddies.Push_back(bdez);
        } else {
          bde = bdez;
        }
        if(is_opt_[next_dex->GetIndex()])
          opt_buddies.Push_back(next_dex);
        if(is_opt_[next_dey->GetIndex()])
          opt_buddies.Push_back(next_dey);
        if(dim_ == 3 && is_opt_[next_dez->GetIndex()])
          opt_buddies.Push_back(next_dez);

        all_buddies.Clear(true);
        all_buddies.Push_back(bdex);
        all_buddies.Push_back(bdey);
        if(dim_ == 3)
          all_buddies.Push_back(bdez);
        all_buddies.Push_back(next_dex);
        all_buddies.Push_back(next_dey);
        if(dim_ == 3)
          all_buddies.Push_back(next_dez);
        // 2D: cases 3 and 4 describe cone in y direction
        if(dim_ == 2) {
          vem.Push_back(Function::Local::Identifier(bde, opt_buddies, all_buddies, 3));
          vem.Push_back(Function::Local::Identifier(bde, opt_buddies, all_buddies, 4));
        } else {
          // 3D: cases 5,6,7 and 8 describe cone in y direction
          vem.Push_back(Function::Local::Identifier(bde, opt_buddies, all_buddies, 5));
          vem.Push_back(Function::Local::Identifier(bde, opt_buddies, all_buddies, 6));
          vem.Push_back(Function::Local::Identifier(bde, opt_buddies, all_buddies, 7));
          vem.Push_back(Function::Local::Identifier(bde, opt_buddies, all_buddies, 8));
        }
        LOG_DBG3(SBD) << "SVMSEM: y opt_buddies=";
        for(auto b : opt_buddies)
          LOG_DBG3(SBD) << b->GetIndex();
        LOG_DBG3(SBD) << "SVMSEM: y all_buddies=";
        for(auto b : all_buddies)
          LOG_DBG3(SBD) << b->GetIndex();
      }
      LOG_DBG3(SBD) << "SVMSEM: y p="  << (prev_dex != NULL ? (int) prev_dex->GetIndex() : -1)
        << " e=" << bde->GetIndex() << " n=" << (next_dex != NULL ? (int) next_dex->GetIndex() : -1);

      if(dim_ == 3) {
        // neighbors in z direction
        curr_idx =  e                         *dim_;
        prev_idx = (e - num_cp_[0]*num_cp_[1])*dim_;
        next_idx = (e + num_cp_[0]*num_cp_[1])*dim_;
        // No need for a row check
        bdex = shape_param_[curr_idx];
        prev_dex = prev ? (prev_idx > -1 ? shape_param_[prev_idx] : NULL) : NULL;
        next_dex =        next_idx < nsp ? shape_param_[next_idx] : NULL;

        curr_idx =  e                         *dim_+1;
        prev_idx = (e - num_cp_[0]*num_cp_[1])*dim_+1;
        next_idx = (e + num_cp_[0]*num_cp_[1])*dim_+1;
        bdey = shape_param_[curr_idx];
        prev_dey = prev ? (prev_idx > -1 ? shape_param_[prev_idx] : NULL) : NULL;
        next_dey =        next_idx < nsp ? shape_param_[next_idx] : NULL;

        curr_idx =  e                         *dim_+2;
        prev_idx = (e - num_cp_[0]*num_cp_[1])*dim_+2;
        next_idx = (e + num_cp_[0]*num_cp_[1])*dim_+2;
        bdez = shape_param_[curr_idx];
        prev_dez = prev ? (prev_idx > -1 ? shape_param_[prev_idx] : NULL) : NULL;
        next_dez =        next_idx < nsp ? shape_param_[next_idx] : NULL;


        if(prev_dex) {
          opt_buddies.Clear(true);
          if(is_opt_[bdex->GetIndex()]) {
            bde = bdex;
            if(is_opt_[bdey->GetIndex()])
              opt_buddies.Push_back(bdey);
            if(dim_ == 3 && is_opt_[bdez->GetIndex()])
              opt_buddies.Push_back(bdez);
          } else if(is_opt_[bdey->GetIndex()]) {
            bde = bdey;
            if(dim_ == 3 && is_opt_[bdez->GetIndex()])
              opt_buddies.Push_back(bdez);
          } else {
            bde = bdez;
          }
          if(is_opt_[prev_dex->GetIndex()])
            opt_buddies.Push_back(prev_dex);
          if(is_opt_[prev_dey->GetIndex()])
            opt_buddies.Push_back(prev_dey);
          if(dim_ == 3 && is_opt_[prev_dez->GetIndex()])
            opt_buddies.Push_back(prev_dez);

          all_buddies.Clear(true);
          all_buddies.Push_back(bdex);
          all_buddies.Push_back(bdey);
          all_buddies.Push_back(bdez);
          all_buddies.Push_back(prev_dex);
          all_buddies.Push_back(prev_dey);
          all_buddies.Push_back(prev_dez);
          // 3D: cases 9,10,11 and 12 describe cone in y direction
          vem.Push_back(Function::Local::Identifier(bde, opt_buddies, all_buddies, -9));
          vem.Push_back(Function::Local::Identifier(bde, opt_buddies, all_buddies, -10));
          vem.Push_back(Function::Local::Identifier(bde, opt_buddies, all_buddies, -11));
          vem.Push_back(Function::Local::Identifier(bde, opt_buddies, all_buddies, -12));
          LOG_DBG3(SBD) << "SVMSEM: z opt_buddies=";
          for(auto b : opt_buddies)
            LOG_DBG3(SBD) << b->GetIndex();
          LOG_DBG3(SBD) << "SVMSEM: z all_buddies=";
          for(auto b : all_buddies)
            LOG_DBG3(SBD) << b->GetIndex();
        }

        if(next_dex) {
          opt_buddies.Clear(true);
          if(is_opt_[bdex->GetIndex()]) {
            bde = bdex;
            if(is_opt_[bdey->GetIndex()])
              opt_buddies.Push_back(bdey);
            if(dim_ == 3 && is_opt_[bdez->GetIndex()])
              opt_buddies.Push_back(bdez);
          } else if(is_opt_[bdey->GetIndex()]) {
            bde = bdey;
            if(dim_ == 3 && is_opt_[bdez->GetIndex()])
              opt_buddies.Push_back(bdez);
          } else {
            bde = bdez;
          }
          if(is_opt_[next_dex->GetIndex()])
            opt_buddies.Push_back(next_dex);
          if(is_opt_[next_dey->GetIndex()])
            opt_buddies.Push_back(next_dey);
          if(dim_ == 3 && is_opt_[next_dez->GetIndex()])
            opt_buddies.Push_back(next_dez);

          all_buddies.Clear(true);
          all_buddies.Push_back(bdex);
          all_buddies.Push_back(bdey);
          all_buddies.Push_back(bdez);
          all_buddies.Push_back(next_dex);
          all_buddies.Push_back(next_dey);
          all_buddies.Push_back(next_dez);
          // 3D: cases 9,10,11 and 12 describe cone in y direction
          vem.Push_back(Function::Local::Identifier(bde, opt_buddies, all_buddies, 9));
          vem.Push_back(Function::Local::Identifier(bde, opt_buddies, all_buddies, 10));
          vem.Push_back(Function::Local::Identifier(bde, opt_buddies, all_buddies, 11));
          vem.Push_back(Function::Local::Identifier(bde, opt_buddies, all_buddies, 12));
          LOG_DBG3(SBD) << "SVMSEM: z opt_buddies=";
          for(auto b : opt_buddies)
            LOG_DBG3(SBD) << b->GetIndex();
          LOG_DBG3(SBD) << "SVMSEM: z all_buddies=";
          for(auto b : all_buddies)
            LOG_DBG3(SBD) << b->GetIndex();
        }
        LOG_DBG3(SBD) << "SVMSEM: z p="  << (prev_dex != NULL ? (int) prev_dex->GetIndex() : -1)
          << " e=" << bde->GetIndex() << " n=" << (next_dex != NULL ? (int) next_dex->GetIndex() : -1);
      }
    }
  }
  vem.Trim();
  LOG_DBG(SBD) << "SVMSEM: final f=" << f->ToString() << " loc=" << locality << " prev=" << prev << " -> vem=" << vem.GetSize();
}

void SplineBoxDesign::MapFeatureToDensity()
{
  assert(map.GetSize() == domain->GetGrid()->GetNumElems());
  assert(map.GetSize() == data.GetSize());

  mapping_timer_->Start();

  LOG_DBG(SBD) << "MSTD: di=" << design_id;

  int res_idx_r = GetSpecialResultIndex(DesignElement::DEFAULT, DesignElement::SPLINE_BOX_INT_ORDER);
  int res_idx_c = GetSpecialResultIndex(DesignElement::DEFAULT, DesignElement::SPLINE_BOX_INT_CORNER);

  // prepare all corner values
  EvalAllCornerValues();

  // statistics - be sure to handle them correctly when running parallel.
  // Set them to NumInt::int_cells_* after the parallel loop, omp just does not allow the class attributes
  int cells_cnt = 0;
  int cells_order_sum = 0;

  Grid* grid = domain->GetGrid();

  #pragma omp parallel num_threads(CFS_NUM_THREADS)
  {
    StdVector<double> local_ip(dim_); // the current ip within the element

    // loop over FE grid elements
    for(unsigned int r = 0; r < map.GetSize(); ++r) {
      Item& item = map[r];
      DesignElement* de = item.elemval;

      Vector<int> order(item.min_corner_value.GetSize());

      // check what we need to integrate. The order is the maximal order of all relevant shapes
      // when we integrate we integrate all shapes at the same integration points.
      // If this could be relaxed we might be able to save quite some time!!
      int max_order = item.GetOrder(order, numInt_); // sets order. 0, 1 or >= 2

      if(res_idx_r >= 0)
        de->specialResult[res_idx_r] = max_order;

      if(res_idx_c >= 0)
        de->specialResult[res_idx_c] = item.MaxDiffCornerValue();

      double rho = -1.0;

      if(max_order == 0)
        rho = 0.0;

      LOG_DBG2(SBD) << "MFTD: de=" << de->elem->elemNum << ", max order=" << max_order << ", order=" << order.ToString() << ", rho=" << rho;

      // the number of integration points per dimension is actually max_order
      int num_ip = max_order;

      if(rho == -1.0) {
        ++cells_cnt;

        rho = 0.0; // such that we can sum the ip to it
        for(int ip_x = 0; ip_x < num_ip; ++ip_x) {
          for(int ip_y = 0; ip_y < num_ip; ++ip_y) {
            for(int ip_z = 0; ip_z < (dim_ == 2 ? 1 : num_ip); ++ip_z) {
              ++cells_order_sum;
              // get local coordinates of ip and integration weights
              double weight = FeaturedDesign::Item::SetIPGetWeight(this, local_ip, ip_x, ip_y, ip_z, numInt_.max_order_);
              // get global coordinates of lower left node of current element
              StdVector<unsigned int> nodes(dim_ == 2 ? 4 : 8);
              grid->GetElemNodes(nodes, r+1);
              Matrix<double> coords;
              grid->GetElemNodesCoord(coords, nodes, false);

              // Get element edge length
              double sze, tmp;
              LagrangeElemShapeMap sm(domain->GetGrid());
              sm.SetElem(de->elem, false);
              sm.GetMaxMinEdgeLength(sze, tmp);

              // add local ip to global node.
              Vector<double> ip(dim_);
              for(unsigned int i = 0; i < dim_; ++i) {
                ip[i] = coords[i][0] + local_ip[i] * sze;
              }

              double ip_rho = GetDensityAtCoord( Eval(ip) ); // the final value for one integration point.
              assert(ip_rho >= 0.0 && ip_rho <= 1.01);
              rho += weight * ip_rho;
              LOG_DBG3(SBD) << "MFTD: de=" << de->GetIndex() << " p=" << de->GetLocation()->ToString()
                  << " local_ip=" << local_ip.ToString() << " ip=" << ip.ToString() << " Eval(ip)=" << Eval(ip).ToString()
                  << " ip_rho=" << ip_rho;
            } // end ip_z
          } // end ip_y
        } // end ip_x
      } // end real integration
      assert(rho >= 0.0 && rho <= 1.01);
      de->SetDesign(de->GetLowerBound() + (de->GetUpperBound() - de->GetLowerBound()) * rho); // we assume 0 <= rho <= 1
      LOG_DBG2(SBD) << "MFTD: el=" << de->elem->elemNum << " avg=" << de->GetPlainDesignValue()
                    << " delta=" << (de->GetPlainDesignValue() - de->GetLowerBound());
      assert(de->GetPlainDesignValue() >= de->GetLowerBound() - 1e-10);
    } // end loop over density elements
  } // end of omp parallel section

  numInt_.int_cells_cnt_ = cells_cnt;
  numInt_.int_cells_order_sum_ = cells_order_sum;
  numInt_.ToInfo(info_->Get("splineBox/numInt"));
  mapped_design_ = design_id;
  mapping_timer_->Stop();
}

void SplineBoxDesign::MapFeatureGradient(Function* f)
{
  assert(design_id == mapped_design_); // we need the Item setting from MapShapeDesign for the current design!
  assert(!(!f->IsObjective() && dynamic_cast<Condition*>(f)->IsLocalCondition())); // it makes no sense for a local condition!!

  gradient_timer_->Start();

  // Optimization::EvalObjectiveConstraints() triggers MapFeatureGradient() via WriteGradientToExtern() but rho::constraintGrad might not be set yet
  // From DesignSpace and ErsatzMaterial we have d_function/d_rho.
  // By chain rule we get d_function/d_shape by summing up for each shape_var the corresponding d_function/d_rho.
  // Note that we do this based on each integration point!
  //
  // !! Also all simp gradients are not computed before the first WriteGradientToExtern() which triggers this MapShapeGradient

  // shall we store max_grad for output? -1 if not

  int res_idx_dip_rho_x = GetSpecialResultIndex(DesignElement::DENSITY, DesignElement::SPLINE_BOX_GRAD_X, DesignElement::SP_CP);
  int res_idx_dip_rho_y = GetSpecialResultIndex(DesignElement::DENSITY, DesignElement::SPLINE_BOX_GRAD_Y, DesignElement::SP_CP);
  int res_idx_dip_rho_z = GetSpecialResultIndex(DesignElement::DENSITY, DesignElement::SPLINE_BOX_GRAD_Z, DesignElement::SP_CP);

  // to speed up performance and allow parallelization we have to not call BaseDesignElement::AddGradient()
  // within each integration point but have a flat vector which is added after the map loop
  Vector<double> shape_f_grad(opt_shape_param_.GetSize());
  shape_f_grad.Init(0.0);

  Grid* grid = domain->GetGrid();

  //#pragma omp parallel num_threads(CFS_NUM_THREADS)
  {
    // this are thread private constructs to be reused over the for loop iterations

    // as we have no OpenMP 4.5 but 3.1 we cannot use array reduction or reduction functions :(
    Vector<double> private_shape_f_grad(shape_f_grad); // the local array is initialized with 0
    assert(private_shape_f_grad.Max() == 0.0);

    // within the element coordinates we perform the integration
    Vector<unsigned int> idx(dim_);
    StdVector<double>    local_ip(dim_); // the current ip within the element

    Vector<int> order(map[0].min_corner_value.GetSize());

    // loop over FE grid elements -> why not parallel?!
    for(unsigned int r = 0; r < map.GetSize(); r++)
    {
      Item& item = map[r];
      DesignElement* de = item.elemval;
      Vector<double> log_dip_rho(dim_, 0.0);

      int max_order = item.GetOrder(order, numInt_); // sets order. 0, 1 or >= 2

      LOG_DBG2(SBD) << "MSG: f=" << f->ToString() << " de=" << de->elem->elemNum
          << " rho=" << de->GetPlainDesignValue() << " maxorder=" << max_order << " order=" << order.ToString();

      // if max_order = 0=void there is no gradient to add.
      if(max_order >= 1) {
        // case were we have a non-zero gradient!
        double de_plain_f_grad = de->GetPlainGradient(f);
        assert(!std::isnan(de_plain_f_grad));

        // @see SplineBoxDesign::MapFeatureToDensity()
        int num_ip = max_order;

        for(int ip_x = 0; ip_x < num_ip; ++ip_x) {
          for(int ip_y = 0; ip_y < num_ip; ++ip_y) {
            for(int ip_z = 0; ip_z < (dim_ == 2 ? 1 : num_ip); ++ip_z) {
              // get local coordinates of ip and integration weights
              double weight = FeaturedDesign::Item::SetIPGetWeight(this, local_ip, ip_x, ip_y, ip_z, numInt_.max_order_);
              // get global coordinates of lower left node of current element
              StdVector<unsigned int> nodes(dim_ == 2 ? 4 : 8);
              grid->GetElemNodes(nodes, r+1);
              Matrix<double> coords;
              grid->GetElemNodesCoord(coords, nodes, false);

              // Get element edge length
              double sze, tmp;
              LagrangeElemShapeMap sm(domain->GetGrid());
              sm.SetElem(de->elem,false);
              sm.GetMaxMinEdgeLength(sze, tmp);

              // add local ip to global node. no casting available ...
              Vector<double> ip(dim_);
              for(unsigned int i = 0; i < dim_; ++i) {
                ip[i] = coords[i][0] + local_ip[i] * sze;
              }

              // gradient of feature density w.r.t. integration point
              Vector<double> dip_rho = GetDensityDerivativeAtCoord( Eval(ip) );
              assert(!dip_rho.ContainsNaN() && !dip_rho.ContainsInf());

              Vector<double> dip_rho_normalized;
              dip_rho_normalized = dip_rho * (de->GetUpperBound() - de->GetLowerBound()) * weight;
              log_dip_rho += dip_rho_normalized;

              // derivative of integration point w.r.t. control point
              assert(shape_param_.GetSize() == total_num_cp_*dim_);
              Vector<double> dcp_rho(opt_shape_param_.GetSize());
              Matrix<double> mtx = GetValuesOfBasisFunctions(ip);
              unsigned int k = 0;
              for(unsigned int cp = 0; cp < total_num_cp_; ++cp) {
                for(unsigned int d = 0; d < dim_; ++d) {
                  if(is_opt_[cp*dim_+d]) {
                    dcp_rho[k++] = dip_rho_normalized[d] * mtx[0][cp];
                  }
                }
              }

              private_shape_f_grad += dcp_rho * de_plain_f_grad;
              LOG_DBG3(SBD) << "MSG: de=" << de->GetIndex() << " p=" << de->GetLocation()->ToString()
                  << " local_ip=" << local_ip.ToString() << " ip=" << ip.ToString() << " Eval(ip)=" << Eval(ip).ToString()
                  << " dip_rho=" << dip_rho.ToString() << " dcp_rho=" << dcp_rho.ToString() << " fg=drho_f=" << de_plain_f_grad;
            } // end ip_z
          } // end ip_y
        } // end ip_x
      } // normalize by integration points.
      LOG_DBG3(SBD) << "MSG: el=" << de->elem->elemNum << " rho=" << de->GetPlainDesignValue() << " dip_rho=" << log_dip_rho.ToString();

      if(res_idx_dip_rho_x >= 0)
        de->specialResult[res_idx_dip_rho_x] = log_dip_rho[0];
      if(res_idx_dip_rho_y >= 0)
        de->specialResult[res_idx_dip_rho_y] = log_dip_rho[1];
      if(res_idx_dip_rho_z >= 0 && dim_ == 3)
        de->specialResult[res_idx_dip_rho_z] = log_dip_rho[2];
    } // end loop over density elements

    // now gather private_shape_f_grad to shape_f_grad as we may not use array reduction yet
    // https://stackoverflow.com/questions/20413995/reducing-on-array-in-openmp
    #pragma omp critical
    {
      // we are in the omp parallel section with n threads. The critical only of the n threads at one time
      // but this is executed for all.
      shape_f_grad.Add(private_shape_f_grad);
    }
  } // end of omp parallel

  // write back shape_f_grad
  LOG_DBG2(SBD) << "MSG: f=" << f->ToString() << " sfg=" << shape_f_grad.ToString();

  assert(shape_f_grad.GetSize() == opt_shape_param_.GetSize());
  for(unsigned int i = 0; i < shape_f_grad.GetSize(); i++) {
    // if it is not an opt variable DesignElement::*Gradient has size zero. Note negative gradients!
    if(shape_f_grad[i] != 0) {
      opt_shape_param_[i]->AddGradient(f, shape_f_grad[i]);
    }
  }
  gradient_timer_->Stop();
}

void SplineBoxDesign::EvalAllCornerValues()
{
  // for a mesh we evaluate all points and set each value for all adjacent Item.
  // Hence the value is repeated almost 4 times in 2D and 8 times in 3D
  Grid* grid = domain->GetGrid();

  Vector<double> glob;
  glob.Resize(grid->GetNumNodes());

  // get transformed mesh coordinates (forward_) or control points respectively
  Matrix<double> new_coords;
  EvalAll(new_coords);

//  LOG_DBG3(SBD) << "EACV: new_coords: " << new_coords.ToString(0);

  unsigned int nodes_per_elem = dim_ == 2 ? 4 : 8;
  StdVector<unsigned int> nodes(nodes_per_elem);
  Vector<double> vals(nodes_per_elem);
  for(unsigned int e = 0; e < grid->GetNumElems(); ++e) {
    // get nodes of element
    // node and element numbers are 1-based!
    grid->GetElemNodes(nodes, e+1);
    // get values at nodes
    for(unsigned int i = 0; i < nodes_per_elem; ++i) {
      if(local_index_[nodes[i]-1] != -1.0) {
        Vector<double> p(dim_);
        new_coords.GetRow(p, local_index_[nodes[i]-1]);
        vals[i] = GetDensityAtCoord(p);
      } else {
        // Todo
        vals[i] = 0.0;
        throw Exception("Not implemented");
      }
    }

    map[e].min_corner_value[0] = vals.Min();
    map[e].max_corner_value[0] = vals.Max();
    LOG_DBG2(SBD) << "EACV: e= " << e << " min= " << map[e].min_corner_value[0] << " max= " << map[e].max_corner_value[0];
  }
}

double SplineBoxDesign::GetDensityAtCoord(Vector<double> point) const
{
  double val = 0.0;
//  if(!IsInside(point)) {
//    throw Exception("point outside of spline box");
//    return 0.0; // never reached
//  } else {
  if(analyticFunc_ == AnalyticFunc::FILE) {
    StdVector<double> relative_coords(dim_);
    for(unsigned int d = 0; d < dim_; ++d) {
      // Assume that bounding_box of initial spline box is equal to bounding box of density, i.e.
      // spline box covers density. Thus relative_coords are relative to (constant) density bounding box.
      relative_coords[d] = (point[d] - bounding_box_[d][0]) / (bounding_box_[d][1] - bounding_box_[d][0]);
      if(periodic_) {
        while(relative_coords[d] < 0) {
          ++relative_coords[d];
        }
        while(relative_coords[d] > 1) {
          --relative_coords[d];
        }
      }
    }
    if(interpolation_ == Interpolation::CUBIC) {
      if(dim_ == 2) {
        val = bicubicInterpolator_->EvaluateFunc(relative_coords[0], relative_coords[1]);
      } else {
        val = tricubicInterpolator_->EvaluateFunc(relative_coords[0], relative_coords[1], relative_coords[2]);
      }
      // Interpolated density might be outside (0,1).
      val = SmoothMin(SmoothMax(val, 0.0, beta_), 1.0, beta_);
    } else {
      StdVector<int> sub(dim_);
      for(unsigned int d = 0; d < dim_; ++d) {
        sub[d] = std::min((unsigned int) (relative_coords[d] * density_resolution_[d]), density_resolution_[d]-1);
      }
      unsigned int idx;
      Sub2Ind(density_resolution_, sub, idx);
      val = density_[idx];
    }
  } else {
    double relative_coords;
    StdVector<double> sub(dim_);
    Vector<double> values(dim_);
    for(unsigned int d = 0; d < dim_; ++d) {
      relative_coords = (point[d] - bounding_box_[d][0]) / (bounding_box_[d][1] - bounding_box_[d][0]);
      // shift sine by pi/2 -> could use cosine instead
      sub[d] = (relative_coords / feature_scale_ + 1./4.) * 2. * M_PI;
      values[d] = std::sin(sub[d]) / 2. + .5;
    }

    switch(analyticFunc_) {
    case AnalyticFunc::SUM_OF_SINE:
      val = (values[0] + values[1] + (dim_ == 3 ? values[2] : 0.0)) / dim_;
      break;
    case AnalyticFunc::MAX_SINE:
      val = SmoothMax(values, beta_);
      break;
    case AnalyticFunc::SINE_X:
      val = values[0];
      break;
    case AnalyticFunc::SINE_Y:
      val = values[1];
      break;
    case AnalyticFunc::SINE_Z:
      assert(dim_ == 3);
      val = values[2];
      break;
    default:
      throw Exception("Function not implemented");
    }
  }
  LOG_DBG2(SBD) << "EAC: p=" << point.ToString() << " -> f=" << val;
  return val;
}

Vector<double> SplineBoxDesign::GetDensityDerivativeAtCoord(Vector<double> point) const
{
  Vector<double> out(dim_, 0.0);

//  if(!IsInside(point)) {
//    throw Exception("point outside of spline box");
//    return out; // never reached
//  } else {
  if(analyticFunc_ == AnalyticFunc::FILE) {
    StdVector<double> relative_coords(dim_);
    for(unsigned int d = 0; d < dim_; ++d) {
      // @see GetDensityAtCoord
      relative_coords[d] = (point[d] - bounding_box_[d][0]) / (bounding_box_[d][1] - bounding_box_[d][0]);
      if(periodic_) {
        while(relative_coords[d] < 0) {
          ++relative_coords[d];
        }
        while(relative_coords[d] > 1) {
          --relative_coords[d];
        }
      }
    }
    if(interpolation_ == Interpolation::CUBIC) {
      double val;
      if(dim_ == 2) {
        out = bicubicInterpolator_->EvaluatePrime(relative_coords[0], relative_coords[1]);
        val = bicubicInterpolator_->EvaluateFunc(relative_coords[0], relative_coords[1]);
      } else {
        out = tricubicInterpolator_->EvaluatePrime(relative_coords[0], relative_coords[1], relative_coords[2]);
        val = tricubicInterpolator_->EvaluateFunc(relative_coords[0], relative_coords[1], relative_coords[2]);
      }
      // chain rule
      out *= DerivSmoothMax(val, 0.0, beta_, -1);
      out *= DerivSmoothMin(SmoothMax(val, 0.0, beta_), 1.0, beta_, -1);
      for(unsigned int d = 0; d < dim_; ++d) {
        // derivative normalized w.r.t. size of spline box
        out[d] *= 1/(bounding_box_[d][1]-bounding_box_[d][0]);
      }
    } else {
      StdVector<int> sub(dim_);
      for(unsigned int d = 0; d < dim_; ++d) {
        sub[d] = std::min((unsigned int) (relative_coords[d] * density_resolution_[d]), density_resolution_[d]-1);
      }
      unsigned int idx;
      Sub2Ind(density_resolution_, sub, idx);

      for(unsigned int d = 0; d < dim_; ++d) {
        out[d] = density_derivative_[d][idx];
      }
    }
    //LOG_DBG(SBD) << "EAC: point= " << point.ToString() << " relative_coords= " << sub.ToString();
  } else {
    double relative_coords;
    StdVector<double> sub(dim_);
    Vector<double> values(dim_);
    for(unsigned int d = 0; d < dim_; ++d) {
      relative_coords = (point[d] - bounding_box_[d][0]) / (bounding_box_[d][1] - bounding_box_[d][0]);
      sub[d] = (relative_coords / feature_scale_ + 1./4.) * 2. * M_PI;
      values[d] = std::sin(sub[d]) / 2. + .5;
    }
    for(unsigned int d = 0; d < dim_; ++d) {
      switch(analyticFunc_) {
      case AnalyticFunc::SUM_OF_SINE:
        out[d] = 1.0 / dim_;
        break;
      case AnalyticFunc::MAX_SINE:
        out[d] = DerivSmoothMax(values, beta_, d);
        break;
      case AnalyticFunc::SINE_X:
        out[d] = d == 0 ? 1.0 : 0.0;
        break;
      case AnalyticFunc::SINE_Y:
        out[d] = d == 1 ? 1.0 : 0.0;
        break;
      case AnalyticFunc::SINE_Z:
        out[d] = d == 2 ? 1.0 : 0.0;
        break;
      default:
        throw Exception("Function not implemented");
      }
      out[d] *= std::cos(sub[d]) / (bounding_box_[d][1] - bounding_box_[d][0]) / feature_scale_ * M_PI;
    }
  }
  LOG_DBG2(SBD) << "EAC: p=" << point.ToString() << " -> out=" << out.ToString();
  return out;
}

void SplineBoxDesign::ReadDensityXml(PtrParamNode set, double& lower_violation, double& upper_violation)
{
  // Todo
  MapFeatureToDensity();
}

Matrix<double> SplineBoxDesign::GetInjectivityMatrix()
{
  double tana = std::tan(45./180. * M_PI);
  unsigned int nRows = 0.0;
  Vector<unsigned int> ncp(dim_);
  for(unsigned int d = 0; d < dim_; ++d) {
    ncp = num_cp_;
    --ncp[d];
    nRows += ncp.Product();
  }
  nRows *= std::pow(2, dim_-1);

  Matrix<double> mtx(nRows, num_cp_.Product()*dim_);

  unsigned int idx = 0;
  if(dim_ == 2) {
    // finite differences in x direction
    for(unsigned int y = 0; y < num_cp_[1]; ++y) {
      for(unsigned int x = 0; x < num_cp_[0] - 1; ++x) {
        unsigned int index1 = y * num_cp_[0] + x;
        unsigned int index2 = y * num_cp_[0] + x + 1;
        index1 *= dim_;
        index2 *= dim_;

        // a_1 * x >= y -> - a_1 * x + y <= 0 -> - a_1 * x_1 + a_1 * x_0 + y_1 - y_0 <= 0
        mtx[idx][index1    ] =  tana;
        mtx[idx][index1 + 1] =    -1;
        mtx[idx][index2    ] = -tana;
        mtx[idx][index2 + 1] =     1;
        ++idx;
        // a_2 * x >= -y -> - a_2 * x - y <= 0 -> - a_2 * x_1 + a_2 * x_0 - y_1 + y_0 <= 0
        mtx[idx][index1    ] =  tana;
        mtx[idx][index1 + 1] =     1;
        mtx[idx][index2    ] = -tana;
        mtx[idx][index2 + 1] =    -1;
        ++idx;
      }
    }
    // finite differences in y direction
    for(unsigned int y = 0; y < num_cp_[1] - 1; ++y) {
      for(unsigned int x = 0; x < num_cp_[0]; ++x) {
        unsigned int index1 = y * num_cp_[0] + x;
        unsigned int index2 = y * num_cp_[0] + x + num_cp_[0];
        index1 *= dim_;
        index2 *= dim_;

        // b_1 * y >= x -> x - b_1 * y <= 0 -> x_1 - x_0 - b_1 * y_1 + b_1 * y_0 <= 0
        mtx[idx][index1    ] =    -1;
        mtx[idx][index1 + 1] =  tana;
        mtx[idx][index2    ] =     1;
        mtx[idx][index2 + 1] = -tana;
        ++idx;
        // b_2 * y >= - x -> - x - b_2 * y <= 0 -> - x_1 + x_0 - b_2 * y_1 + b_2 * y_0 <= 0
        mtx[idx][index1    ] =     1;
        mtx[idx][index1 + 1] =  tana;
        mtx[idx][index2    ] =    -1;
        mtx[idx][index2 + 1] = -tana;
        ++idx;
      }
    }
  } else {
    // finite differences in x direction
    for(unsigned int z = 0; z < num_cp_[2]; ++z) {
      for(unsigned int y = 0; y < num_cp_[1]; ++y) {
        for(unsigned int x = 0; x < num_cp_[0] - 1; ++x) {
          unsigned int index1 = z * num_cp_[1] * num_cp_[0] + y * num_cp_[0] + x;
          unsigned int index2 = z * num_cp_[1] * num_cp_[0] + y * num_cp_[0] + x+1;
          index1 *= dim_;
          index2 *= dim_;

          // x >= y + z -> - x + y + z <= 0 -> - x_1 + x_0 + y_1 - y_0 + z_1 - z_0 <= 0
          mtx[idx][index1    ] =  1;
          mtx[idx][index1 + 1] = -1;
          mtx[idx][index1 + 2] = -1;
          mtx[idx][index2    ] = -1;
          mtx[idx][index2 + 1] =  1;
          mtx[idx][index2 + 2] =  1;
          ++idx;
          // x >= y - z -> - x + y - z <= 0 -> - x_1 + x_0 + y_1 - y_0 - z_1 + z_0 <= 0
          mtx[idx][index1    ] =  1;
          mtx[idx][index1 + 1] = -1;
          mtx[idx][index1 + 2] =  1;
          mtx[idx][index2    ] = -1;
          mtx[idx][index2 + 1] =  1;
          mtx[idx][index2 + 2] = -1;
          ++idx;
          // x >= - y + z -> - x - y + z <= 0 -> - x_1 + x_0 - y_1 + y_0 + z_1 - z_0 <= 0
          mtx[idx][index1    ] =  1;
          mtx[idx][index1 + 1] =  1;
          mtx[idx][index1 + 2] = -1;
          mtx[idx][index2    ] = -1;
          mtx[idx][index2 + 1] = -1;
          mtx[idx][index2 + 2] =  1;
          ++idx;
          // x >= - y - z -> - x - y - z <= 0 -> - x_1 + x_0 - y_1 + y_0 - z_1 + z_0 <= 0
          mtx[idx][index1    ] =  1;
          mtx[idx][index1 + 1] =  1;
          mtx[idx][index1 + 2] =  1;
          mtx[idx][index2    ] = -1;
          mtx[idx][index2 + 1] = -1;
          mtx[idx][index2 + 2] = -1;
          ++idx;
        }
      }
    }
    // finite differences in y direction
    for(unsigned int z = 0; z < num_cp_[2]; ++z) {
      for(unsigned int y = 0; y < num_cp_[1] - 1; ++y) {
        for(unsigned int x = 0; x < num_cp_[0]; ++x) {
          unsigned int index1 = z * num_cp_[1] * num_cp_[0] + y * num_cp_[0] + x;
          unsigned int index2 = z * num_cp_[1] * num_cp_[0] + (y+1) * num_cp_[0] + x;
          index1 *= dim_;
          index2 *= dim_;

          // y >= x - z -> x - y - z <= 0 -> x_1 - x_0 - y_1 + y_0 - z_1 + z_0 <= 0
          mtx[idx][index1    ] = -1;
          mtx[idx][index1 + 1] =  1;
          mtx[idx][index1 + 2] =  1;
          mtx[idx][index2    ] =  1;
          mtx[idx][index2 + 1] = -1;
          mtx[idx][index2 + 2] = -1;
          ++idx;
          // y >= x + z -> x - y + z <= 0 -> x_1 - x_0 - y_1 + y_0 + z_1 - z_0 <= 0
          mtx[idx][index1    ] = -1;
          mtx[idx][index1 + 1] =  1;
          mtx[idx][index1 + 2] = -1;
          mtx[idx][index2    ] =  1;
          mtx[idx][index2 + 1] = -1;
          mtx[idx][index2 + 2] =  1;
          ++idx;
          // y >= - x + z -> - x - y + z <= 0 -> - x_1 + x_0 - y_1 + y_0 + z_1 - z_0 <= 0
          mtx[idx][index1    ] =  1;
          mtx[idx][index1 + 1] =  1;
          mtx[idx][index1 + 2] = -1;
          mtx[idx][index2    ] = -1;
          mtx[idx][index2 + 1] = -1;
          mtx[idx][index2 + 2] =  1;
          ++idx;
          // y >= - x - z -> - x - y - z <= 0 -> - x_1 + x_0 - y_1 + y_0 - z_1 + z_0 <= 0
          mtx[idx][index1    ] =  1;
          mtx[idx][index1 + 1] =  1;
          mtx[idx][index1 + 2] =  1;
          mtx[idx][index2    ] = -1;
          mtx[idx][index2 + 1] = -1;
          mtx[idx][index2 + 2] = -1;
          ++idx;
        }
      }
    }
    // finite differences in z direction
    for(unsigned int z = 0; z < num_cp_[2] - 1; ++z) {
      for(unsigned int y = 0; y < num_cp_[1]; ++y) {
        for(unsigned int x = 0; x < num_cp_[0]; ++x) {
          unsigned int index1 = z * num_cp_[1] * num_cp_[0] + y * num_cp_[0] + x;
          unsigned int index2 = (z+1) * num_cp_[1] * num_cp_[0] + y * num_cp_[0] + x;
          index1 *= dim_;
          index2 *= dim_;

          // x >= y + z -> - x + y + z <= 0 -> - x_1 + x_0 + y_1 - y_0 + z_1 - z_0 <= 0
          mtx[idx][index1    ] = -1;
          mtx[idx][index1 + 1] =  1;
          mtx[idx][index1 + 2] =  1;
          mtx[idx][index2    ] =  1;
          mtx[idx][index2 + 1] = -1;
          mtx[idx][index2 + 2] = -1;
          ++idx;
          // z >= - x + y -> - x + y - z <= 0 -> - x_1 + x_0 + y_1 - y_0 - z_1 + z_0 <= 0
          mtx[idx][index1    ] =  1;
          mtx[idx][index1 + 1] = -1;
          mtx[idx][index1 + 2] =  1;
          mtx[idx][index2    ] = -1;
          mtx[idx][index2 + 1] =  1;
          mtx[idx][index2 + 2] = -1;
          ++idx;
          // z >= x + y -> x + y - z <= 0 -> x_1 - x_0 + y_1 - y_0 - z_1 + z_0 <= 0
          mtx[idx][index1    ] = -1;
          mtx[idx][index1 + 1] = -1;
          mtx[idx][index1 + 2] =  1;
          mtx[idx][index2    ] =  1;
          mtx[idx][index2 + 1] =  1;
          mtx[idx][index2 + 2] = -1;
          ++idx;
          // z >= - x - y -> - x - y - z <= 0 -> - x_1 + x_0 - y_1 + y_0 - z_1 + z_0 <= 0
          mtx[idx][index1    ] =  1;
          mtx[idx][index1 + 1] =  1;
          mtx[idx][index1 + 2] =  1;
          mtx[idx][index2    ] = -1;
          mtx[idx][index2 + 1] = -1;
          mtx[idx][index2 + 2] = -1;
          ++idx;
        }
      }
    }
  }

  return mtx;
}

void SplineBoxDesign::ToInfo(ErsatzMaterial* em)
{
  AuxDesign::ToInfo(em);

  PtrParamNode sb = info_->Get("splineBox");
  PtrParamNode msh = sb->Get("mesh");
  numInt_.ToInfo(sb->Get("numInt"));
  PtrParamNode base = info_->Get("designVariables");
  for(unsigned int i = 0; i < shape_param_.GetSize(); i++) {
    PtrParamNode sBP = base->Get("splineBoxParam", ParamNode::APPEND);
    ShapeParamElement* de = shape_param_[i];
    sBP->Get("lower")->SetValue(de->GetLowerBound());
    sBP->Get("upper")->SetValue(de->GetUpperBound());
    sBP->Get("design")->SetValue(de->GetPlainDesignValue());
  }
}

void SplineBoxDesign::WriteBoxToVTK()
{
  unsigned int iteration = this->opt->GetCurrentIteration();
  // add leading zeros
  double nzeros = std::ceil(std::log10(domain->GetOptimization()->GetMaxIterations())) - std::to_string(iteration).length();
  std::string iter = std::string(nzeros > 0 ? nzeros : 0, '0') + std::to_string(iteration);

  string filename = "results_hdf5/" + progOpts->GetSimName() + ".spb_" + iter + ".vtk";

  std::ofstream out;
  out.open(filename);

  string header = "";
  header += "# vtk DataFile Version 2.6\n";
  header += "splinebox file\n";
  header += "ASCII\n";
  header += "DATASET UNSTRUCTURED_GRID\n\n";
  out << header;

  // (control) points
  out << "POINTS " + std::to_string(total_num_cp_) + " float\n";
  out.setf(std::ios::scientific); // change to scientific notation
  for(unsigned int i = 0; i < total_num_cp_; ++i) {
    for(unsigned int d = 0; d < dim_; ++d) {
      out << control_points_[i][d] << " ";
    }
    if(dim_ == 2) {
      out << "0.0";
    }
    out << "\n";
  }
  out << "\n";
  // revert to normal notation
  out.unsetf(std::ios::fixed | std::ios::scientific);

  // number of lines in x dimension
  unsigned int nLines = 0;
  nLines += (num_cp_[0]-1)*num_cp_[1]*num_cp_[2];
  nLines += num_cp_[0]*(num_cp_[1]-1)*num_cp_[2];
  nLines += num_cp_[0]*num_cp_[1]*(num_cp_[2]-1);

  // lines (i.e. "cells" in vtk language)
  out << "CELLS " + std::to_string(nLines) + " " + std::to_string(3*nLines) + "\n";
  for(unsigned int nz = 0; nz < num_cp_[2]; ++nz) {
    for(unsigned int ny = 0; ny < num_cp_[1]; ++ny) {
      for(unsigned int nx = 0; nx < num_cp_[0]; ++nx) {
        unsigned int node = nz*num_cp_[1]*num_cp_[0] + ny*num_cp_[0] + nx;
        // line in x dimension
        if(nx < num_cp_[0]-1) {
          out << "2 " << std::to_string(node) << " " << std::to_string(node+1) << "\n";
        }
        // line in y dimension
        if(ny < num_cp_[1]-1) {
          out << "2 " << std::to_string(node) << " " << std::to_string(node + num_cp_[0]) << "\n";
        }
        // line in z dimension
        if(nz < num_cp_[2]-1) {
          out << "2 " << std::to_string(node) << " " << std::to_string(node + num_cp_[1]*num_cp_[0]) << "\n";
        }
      }
    }
  }
  out << "\n";

  out << "CELL_TYPES " + std::to_string(nLines) + "\n";
  for(unsigned int i = 0; i < nLines; ++i) {
    out << "3\n";
  }
  out.close();
}

void SplineBoxDesign::SetControlPoint(int idx, Point coords, bool add)
{
  assert(idx < static_cast<int>(total_num_cp_));

  if(add == false) {
    control_points_[idx] = coords;
  } else {
    control_points_[idx] = initial_control_points_[idx] + coords;
  }
}

void SplineBoxDesign::SetControlPoints(StdVector<Point> coords, bool add)
{
  assert(coords.GetSize() == total_num_cp_);
  for(auto p : coords)
    LOG_DBG3(SBD) << "SCP: setting control point to " << p.ToString();

  if(add == false) {
    control_points_ = coords;
  } else {
    for(unsigned int i = 0; i < coords.GetSize(); ++i) {
      control_points_[i] = initial_control_points_[i] + coords[i];
    }
  }
}

void SplineBoxDesign::UpdateFEMesh(Matrix<double> new_coords)
{
  assert(forward_);
  Grid* grid = domain->GetGrid();
  unsigned int num_nodes = grid->GetNumNodes();
  for(UInt i = 0; i < num_nodes; ++i) {
    if(local_index_[i] != -1.0) {
      Vector<Double> p;
      new_coords.GetRow(p, local_index_[i]);
      grid->SetNodeCoordinate(i+1, p); //1-based
    }
  }
}

void SplineBoxDesign::ReadFeature(string file_in, string key)
{

  std::cout << "++ Reading feature ... \"" << file_in << "\"" << std::flush;

  PtrParamNode xml = XmlReader::ParseFile(file_in);

  density_resolution_.Resize(dim_);
  density_resolution_[0] = xml->Get("header/mesh/x")->As<unsigned int>();
  density_resolution_[1] = xml->Get("header/mesh/y")->As<unsigned int>();
  if(dim_ == 3 && xml->Has("header/mesh/z")) {
    density_resolution_[2] = xml->Get("header/mesh/z")->As<unsigned int>();
  }

  // find the proper design set. This is either 'first', 'last' or the * in <set id="*"> ...
  PtrParamNode set;
  if (key == "first")
    set = xml->GetList("set")[0];
  if (key == "last")
    set = xml->GetList("set").Last();
  if (set == NULL)
    set = xml->GetByVal("set", "id", key);
  if (xml->Count("set") == 0)
    throw Exception("There are no design sets in the ersatz material file");

  ParamNodeList elems = set->GetList("element"); // we be 0 for shape map

  const unsigned int elsize = elems.GetSize();

  density_.Resize(elsize);
  for (unsigned int e = 0; e < elsize; ++e) {
    unsigned int nr = elems[e]->Get("nr")->As<unsigned int>(); // 1-based
    string name = "design";
    if (elems[e]->Has(name))
      density_[nr-1] = elems[e]->Get(name)->As<double>();
  }

  density_derivative_.Resize(dim_);
  // finite differences
  unsigned int h = 1; // in elements of density field discretization
  unsigned int num_density_elem = density_resolution_.Product();
  for(unsigned int d = 0; d < dim_; ++d) {
    assert(h < density_resolution_[d]);
    density_derivative_[d].Resize(num_density_elem);
    for(unsigned int e = 0; e < num_density_elem; ++e) {
      unsigned int left, right;
      // at upper boundary do backward FD else forward
      left = (e+h) % density_resolution_[d] < h ? e : e+h;
      right = (e+h) % density_resolution_[d] < h ? e-h : e;
      // derivative normalized w.r.t. density resolution
      density_derivative_[d][e] = (density_[left] - density_[right]) / h;
      // derivative normalized w.r.t. size of density
      // TODO normalize w.r.t. cover_box_
      density_derivative_[d][e] *= density_resolution_[d];
      // derivative normalized w.r.t. size of FE mesh
      density_derivative_[d][e] *= 1/(bounding_box_[d][1]-bounding_box_[d][0]);
    }
  }
  std::cout << " -> resolution: " << density_resolution_.ToString(TS_PLAIN,"x") << "\n" << std::flush;
}

void SplineBoxDesign::InterpolateFeature()
{
  if(interpolation_ == Interpolation::CUBIC) {
    StdVector<int> sub(dim_);

    // vectors with "sampling points"
    StdVector<double> x(density_resolution_[0]);
    std::iota(x.Begin(), x.End(), 0);
    for(StdVector<double>::iterator it = x.Begin(); it != x.End(); ++it)
      *it /= (density_resolution_[0]-1);

    StdVector<double> y(density_resolution_[1]);
    std::iota(y.Begin(), y.End(), 0);
    for(StdVector<double>::iterator it = y.Begin(); it != y.End(); ++it)
      *it /= (density_resolution_[0]-1);

    LOG_DBG3(SBD) << "IF: x=" << x.ToString();
    LOG_DBG3(SBD) << "IF: y=" << y.ToString();

    if(dim_ == 2) {
      bicubicInterpolator_ = new BiCubicInterpolate(density_, x, y, true);
      bicubicInterpolator_->CalcApproximation();
    } else {
      StdVector<double> z(density_resolution_[2]);
      std::iota(z.begin(), z.end(), 0);
      for(StdVector<double>::iterator it = z.Begin(); it != z.End(); ++it)
        *it /= (density_resolution_[0]-1);

      LOG_DBG3(SBD) << "IF: z=" << z.ToString();

      tricubicInterpolator_ = new TriCubicInterpolate(density_, x, y, z, true);
      tricubicInterpolator_->CalcApproximation();
    }
  }
}

// Todo: Return also coordinates of integration points. Those will then be
//       subject to basis generation and can be evaluated very fast in MapFeatureToDensity.
StdVector<Vector<Double>> SplineBoxDesign::GetPointsForBasis()
{
  Grid* grid = domain->GetGrid();
  assert(dim_ == grid->GetDim());
  unsigned int num_nodes = grid->GetNumNodes();

  // Vector for all node coordinates for each dimension
  StdVector<Vector<Double>> points;
  points.Reserve(num_nodes);

  // local_index is true (i.e. > 0) if node is inside box, else false (i.e. = -1)
  local_index_.Resize(num_nodes);
  local_index_.Init(-1);
  unsigned int index = 0;
  for(unsigned int i = 0; i < num_nodes; ++i) {
    Vector<Double> node(dim_);
    grid->GetNodeCoordinate(node, i+1, false); // 1-based!
    if(IsInside(node)) {
      local_index_[i] = index;
      points.Push_back(node);
      ++index;
    }
  }
  points.Trim();
  LOG_DBG3(SBD) << "GPFB: local_index: " << local_index_.ToString();
  return points;
}

void SplineBoxDesign::GenerateBasis(StdVector<Vector<Double>> points)
{
  unsigned int num_points = points.GetSize();

  // Get unique coordinates in each dimension.
  // This will save computation time when evaluating the bsplines
  // (especially for regular grids). After evaluating the new coordinates
  // with the unique values, we will reconstruct the original coordinate
  // vectors by means of the index array 'unique_inverse', i.e.
  // unique_coords[i][unique_inverse] = nodes[i];
  unique_coords_.Resize(dim_);
  unique_inverse_.Resize(dim_);
  basis_.Resize(dim_);
  for(unsigned int d = 0; d < dim_; ++d) {
    // Fastest way of unique in C++ is using unordered_set
    std::unordered_set<double> set;
    for(unsigned int j = 0; j < num_points; ++j) {
      set.insert(points[j][d]); // value is only inserted if not already present
    }
    // convert back to StdVector
    StdVector<double> unique_coordsd;
    unique_coordsd.Reserve(set.size());
    for(std::unordered_set<double>::iterator it = set.begin(); it != set.end(); ++it) {
      unique_coordsd.Push_back(*it);
    }
    std::sort(unique_coordsd.Begin(), unique_coordsd.End()); // might be unnecessary
    unique_coords_[d] = &unique_coordsd;

    // Find inverse indices
    unique_inverse_[d].Resize(num_points);
    for(unsigned int j = 0; j < num_points; ++j) {
      for(unsigned int index = 0; index < unique_coordsd.GetSize(); ++index) {
        if(unique_coordsd[index] == points[j][d]) {
          unique_inverse_[d][j] = index;
          break;
        }
      }
    }

    // evaluate basis to get spline coefficients
    basis_[d] = bsplines_[d]->Eval(&unique_coordsd);

    LOG_DBG3(SBD) << "GB: unique_coords[" << d << "]: \n" << unique_coordsd.ToString();
    LOG_DBG3(SBD) << "GB: unique_inverse_[" << d << "]: \n" << unique_inverse_[d].ToString();
    LOG_DBG3(SBD) << "GB: basis_[" << d << "]: \n" << basis_[d].ToString();
  }
}

bool SplineBoxDesign::IsInside(Vector<double> point) const
{
  // check if point inside bounding_box
  bool inside = true;
  for(unsigned int d = 0; d < dim_; ++d) {
    inside = inside && (bounding_box_[d][0] <= point[d] && point[d] <= bounding_box_[d][1]);
  }
  return inside;
}

Matrix<double> SplineBoxDesign::GetValuesOfBasisFunctions()
{
  assert(basis_.GetSize() > 0);
  assert(basis_[0].GetNumCols() > 0);

  // Get number of nodes local_index splinebox
  unsigned int num_nodes_inside = 0;
  for(unsigned int i = 0; i < local_index_.GetSize(); ++i) {
    if(local_index_[i] >= 0) {
      num_nodes_inside += 1;
    }
  }
  LOG_DBG2(SBD) << "num_nodes_inside: " << num_nodes_inside;

  StdVector<int> sub(dim_);
  Matrix<double> mtx;
  mtx.Resize(num_nodes_inside, total_num_cp_);
  mtx.InitValue(1.0);

  // box is tensor product of 1D splines
  for(unsigned int d = 0; d < dim_; ++d) {
    Matrix<double> mtx_temp;
    mtx_temp.Resize(num_nodes_inside, total_num_cp_);
    for(unsigned int i = 0; i < num_nodes_inside; ++i) {
      for(unsigned int j = 0; j < total_num_cp_; ++j) {
        Ind2Sub(num_cp_, j, sub);
//        LOG_DBG3(SBD) << "GBM: d: " << d << " u_i[d][i]: " << unique_inverse_[d][i] << " sub[d]: " << sub[d];
        mtx_temp[i][j] = basis_[d][unique_inverse_[d][i]][sub[d]];
      }
    }
    mtx = mtx.EntryMult(mtx_temp);
  }

//  LOG_DBG3(SBD) << "GBM: basis matrix: \n" << mtx.ToString();
  return mtx;
}

Matrix<double> SplineBoxDesign::GetValuesOfBasisFunctions(Vector<double> point)
{
  StdVector<int> sub(dim_);
  Matrix<double> mtx;
  mtx.Resize(1, total_num_cp_);
  mtx.InitValue(1.0);

  for(unsigned int d = 0; d < dim_; ++d) {
    StdVector<double> t(1);
    t[0] = point[d];

    // values of all basis functions in dimension d at t
    Matrix<double> basis_t = bsplines_[d]->Eval(&t);

    // box is tensor product of 1D splines
    Matrix<double> mtx_temp;
    mtx_temp.Resize(1, total_num_cp_);
    for(unsigned int j = 0; j < total_num_cp_; ++j) {
      Ind2Sub(num_cp_, j, sub);
      mtx_temp[0][j] = basis_t[0][sub[d]];
    }
    mtx = mtx.EntryMult(mtx_temp);
  }
  return mtx;
}

} // end of namespace
