#include <sstream>
#include <iomanip>
#include <unordered_set>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"
#include "MatVec/Matrix.hh"
#include "MatVec/Vector.hh"
#include "Optimization/Design/SplineBoxDesign.hh"

namespace CoupledField {

DECLARE_LOG(SplBox)
DEFINE_LOG(SplBox, "SplineBox")

//StdVector<double>* bounding_box
SplineBox::SplineBox(StdVector<RegionIdType>& regionIds, PtrParamNode pn, ErsatzMaterial::Method method)
: AuxDesign(regionIds, pn, method)
{
  dim = domain->GetGrid()->GetDim();
  degree = pn->Get("splinebox/degree")->As<unsigned int>();

  num_cp.Reserve(dim);
  num_cp[0] = pn->Get("splinebox/cpx")->As<unsigned int>();
  num_cp[1] = pn->Get("splinebox/cpy")->As<unsigned int>();
  if(dim == 3) {
    num_cp[2] = pn->Get("splinebox/cpz")->As<unsigned int>();
  }

  int prod = 1;
  for(StdVector<unsigned int>::iterator it = num_cp.Begin(); it != num_cp.End(); ++it) {
    prod *= *it;
  }
  total_num_cp = prod;
  LOG_DBG(SplBox) << "SplineBox: number of cp, total number of cp: "
      << num_cp.ToString() << ", " << total_num_cp;

  Matrix<double> bounding_box;
  bounding_box.Resize(dim, 2);
  if(!pn->Has("splinebox/boundingBox")) {
    Grid* grid = domain->GetGrid();
    bounding_box = grid->CalcGridBoundingBox();
  } else {
    bounding_box[0][0] = pn->Get("splinebox/boundingBox/lowerXLimit")->As<double>();
    bounding_box[0][1] = pn->Get("splinebox/boundingBox/upperXLimit")->As<double>();
    bounding_box[1][0] = pn->Get("splinebox/boundingBox/lowerYLimit")->As<double>();
    bounding_box[1][1] = pn->Get("splinebox/boundingBox/upperYLimit")->As<double>();
    if(dim == 3) {
      bounding_box[2][0] = pn->Get("splinebox/boundingBox/lowerZLimit")->As<double>();
      bounding_box[2][1] = pn->Get("splinebox/boundingBox/upperZLimit")->As<double>();
    }
  }

  // Generate bsplines and get greville abscissae
  StdVector<BSpline*> bsplines(dim);
  StdVector<StdVector<double>> GAbscissae;
  for(unsigned int i = 0; i < dim; ++i) {
    BSpline* bspline = new BSpline(degree, num_cp[i]+degree+1);
    bsplines.Push_back(bspline);
    GAbscissae.Push_back(bspline->GrevilleAbscissae());
  }

  // set control_points to greville abscissae
  control_points.Reserve(total_num_cp);
  for(unsigned int i = 0; i < total_num_cp; ++i) {
    StdVector<int> sub(dim);
    ind2sub(num_cp, i, sub);
    Point p;
    for(unsigned int j = 0; j < dim; ++j) {
      p[j] = GAbscissae[j][sub[j]];
    }
    control_points.Push_back(p);
  }
  LOG_DBG(SplBox) << "SplineBox: initial control points: ";
  for(unsigned int i = 0; i < total_num_cp; ++i) {
    LOG_DBG(SplBox) << control_points[i].GetCoordVector().ToString();
  }

  GenerateBasis(bounding_box, bsplines);
}


/* <splinebox degree="3" cpx="4" cpy="5" cpz="3">
     <boundingBox lowerXLimit="0" upperXLimit="1" .../>   optional
     <allControlPoints lower="0" upper="1"/>              oder
     <cp index="0" dof="X|Y|Z" lower="0" upper="1" initial="0.5"/>
   </splinebox> */
void SplineBox::SetupDesign(PtrParamNode pn) {
  unsigned int nshapeparams = total_num_cp * dim;
  double l = -1.0;
  double u = 1.0;
  double v = 0.0;
  scaling_ = 1.0;

  shape_param.Reserve(nshapeparams);

  if(pn->Has("allControlPoints")){
    l = pn->Get("allControlPoints")->Get("lower")->As<double>();
    u = pn->Get("allControlPoints")->Get("upper")->As<double>();
    scaling_ = pn->Get("allControlPoints")->Get("scaling")->As<double>();
    for(unsigned int i = 0; i < nshapeparams; ++i) {
      BaseDesignElement &de = shape_param[i];
      de.SetLowerBound(l);
      de.SetUpperBound(u);
      de.SetDesign(v);
    }
  }

  ParamNodeList cp = pn->GetList("cp");
  if(!cp.IsEmpty()) {
    unsigned int ind;
    for(unsigned int i = 0; i < cp.GetSize(); ++i) {
      unsigned int index = cp[i]->Get("index")->As<unsigned int>();
      std::string dof = cp[i]->Get("dof")->As<std::string>();
      if(cp[i]->Has("lower")) {
        l = cp[i]->Get("lower")->As<double>();
      }
      if(cp[i]->Has("upper")) {
        u = cp[i]->Get("upper")->As<double>();
      }
      if(cp[i]->Has("initial")) {
        v = cp[i]->Get("initial")->As<double>();
      }

      if(dof == "X") {
        ind = index*dim;
      } else if(dof == "Y") {
        ind = index*dim + 1;
      } else {
        ind = index*dim + 2;
      }

      BaseDesignElement &de = shape_param[ind];
      de.SetLowerBound(l);
      de.SetUpperBound(u);
      de.SetDesign(v);
    }
  }
}

int SplineBox::ReadDesignFromExtern(const double* space_in) {
  assert(!std::isnan(scaling_));
  int old_design = design_id;

  bool new_design = false;

  for(unsigned int i = 0; i < shape_param.GetSize(); ++i) {
    double v = space_in[i] * scaling_;
    assert(!std::isnan(v));
    if(!new_design && v != shape_param[i].GetPlainDesignValue()) {
      new_design = true;
    }

    shape_param[i].SetDesign(v);

    //todo: apply periodic bcs

    LOG_DBG3(SplBox) << "RDFE: i=" << i << ", " << shape_param[i].ToString() << " -> " << v;
  }

  // append aux design, might also change design_id
  AuxDesign::ReadDesignFromExtern(space_in);

  if(new_design && design_id <= old_design) {
    StdVector<Point> offset(total_num_cp);
    for(unsigned int i = 0; i < offset.GetSize(); ++i) {
      for(unsigned int j = 0; j < dim; ++j) {
        offset[i][j] = shape_param[i*dim + j].GetDesign();
      }
    }
    SetControlPoints(offset, true);
    UpdateCoordinates();
    ++design_id;
  }

  LOG_DBG(SplBox) << "RDFE: design_id -> " << design_id;
  return design_id;
}

int SplineBox::WriteDesignToExtern(double* space_out, bool scaling) const {
  double rscaling = scaling ? 1.0 / scaling_ : 1.0;

  for(unsigned int i=0; i < shape_param.GetSize(); ++i) {
    space_out[i] = shape_param[i].GetPlainDesignValue() * rscaling;
    LOG_DBG3(SplBox) << "WDTE: out[" << i << "]=" << space_out[i];
  }

  AuxDesign::WriteDesignToExtern(space_out, scaling);

  LOG_DBG(SplBox) << "WDTE: di -> " << design_id;
  return design_id;
}

void SplineBox::WriteGradientToExtern(StdVector<double>& out, DesignElement::ValueSpecifier vs, DesignElement::Access access, Function* f, bool scaling) {

}

void SplineBox::SetControlPoint(int idx, Point coords, bool add) {
  assert(idx < static_cast<int>(total_num_cp));

  if(add == false) {
    control_points[idx] = coords;
  } else {
    control_points[idx] += coords;
  }
}

void SplineBox::SetControlPoints(StdVector<Point> coords, bool add) {
  assert(coords.GetSize() == total_num_cp);

  if(add == false) {
    control_points = coords;
  } else {
    for(unsigned int i = 0; i < coords.GetSize(); ++i) {
      control_points[i] += coords[i];
    }
  }
  LOG_DBG(SplBox) << "SCP: control points: ";
  for(unsigned int i = 0; i < total_num_cp; ++i) {
    LOG_DBG(SplBox) << control_points[i].GetCoordVector().ToString();
  }
}

void SplineBox::UpdateCoordinates() {
  Matrix<double> box_mtx = GetBoxMatrix();

  // "Casting" StdVector<Point> control_points to Matrix<double>
  Matrix<double> cp;
  cp.Resize(total_num_cp, 3);
  cp.Init();
  for(unsigned int i = 0; i < total_num_cp; ++i) {
    Vector<double> p = control_points[i].GetCoordVector();
    for(unsigned int j = 0; j < dim; ++j) {
      cp[i][j] = p[j];
    }
  }
  LOG_DBG(SplBox) << "Eval: control points: \n" << cp.ToString(0);
  Matrix<double> coords;
  // deformed nodes = matrix * control_points
  box_mtx.Mult(cp, coords);
  LOG_DBG3(SplBox) << "Eval: deformed FE nodes: \n" << coords.ToString(0);

  Grid* grid = domain->GetGrid();
  assert(dim == grid->GetDim());
  unsigned int num_nodes = grid->GetNumNodes();
  unsigned int j = 0;
  for(UInt i = 0; i < num_nodes; ++i) {
    if(inside[i] == true) {
      Vector<Double> p;
      coords.GetRow(p, j);
      grid->SetNodeCoordinate(i, p);
      ++j;
    }
  }
}

void SplineBox::EvalDerivative() {
}

inline void SplineBox::sub2ind(StdVector<unsigned int> size, StdVector<int> sub, unsigned int &ind) {
  assert(size.GetSize() == sub.GetSize());

  // cumulative product
  StdVector<int> cumprod = StdVector<int>(size.GetSize());
  cumprod[0] = size[0];
  for(unsigned int i = 1; i < size.GetSize(); ++i) {
    cumprod[i] = cumprod[i-1] * size[i];
  }

  int idx = sub[0] + 1; // zero based
  for(unsigned int i = 1; i < sub.GetSize(); ++i) {
    idx += sub[i] * cumprod[i-1];
  }
  ind = idx - 1; // zero based

  LOG_DBG3(SplBox) << "S2I: sub, ind: " << sub.ToString() << ", " << ind;
}

inline void SplineBox::ind2sub(StdVector<unsigned int> size, unsigned int ind, StdVector<int> &sub) {
  assert(size.GetSize() == sub.GetSize());

  // cumulative product
  StdVector<int> cumprod = StdVector<int>(size.GetSize());
  cumprod[0] = size[0];
  for(unsigned int i = 1; i < size.GetSize(); ++i) {
    cumprod[i] = cumprod[i-1] * size[i];
  }

  ind += 1; //zero based
  for(unsigned int i = sub.GetSize()-1; i > 0; i--) {
    int v = (ind-1)%cumprod[i-1] + 1;
    sub[i] = (ind - v)/cumprod[i-1];
    ind = v;
  }
  sub[0] = ind;

  LOG_DBG3(SplBox) << "I2S: ind, sub: " << ind << ", " << sub.ToString();
}

void SplineBox::GenerateBasis(Matrix<double> bounding_box,StdVector<BSpline*> bsplines) {
  Grid* grid = domain->GetGrid();
  assert(dim == grid->GetDim());
  unsigned int num_nodes = grid->GetNumNodes();

  // Vector for all node coordinates for each dimension
  Matrix<double> nodes;
  nodes.Resize(dim, num_nodes);
  nodes.Init();

  // inside is true if node is inside box, else false
  inside.Resize(grid->GetNumNodes());
  inside.Init(true);
  for(unsigned int i = 0; i < num_nodes; ++i) {
    Vector<Double> node(dim);
    grid->GetNodeCoordinate(node, i, false);
    for(unsigned int j = 0; j < dim; ++j) {
      nodes[j][i] = node[j]; // fill nodes
      // check if inside bounding_box
      inside[i] = inside[i] && (bounding_box[j][0] <= node[j] && node[j] <= bounding_box[j][1]);
    }
  }

  LOG_DBG3(SplBox) << "GB: nodes: \n" << nodes.ToString(0);
  LOG_DBG3(SplBox) << "GB: inside: " << inside.ToString();

  // Get unique coordinates in each dimension.
  // This will save computation time when evaluating the bsplines
  // (especially for regular grids). After evaluating the new coordinates
  // with the unique values, we will reconstruct the original coordinate
  // vectors by means of the index array 'unique_inverse', i.e.
  // unique_coords[i][unique_inverse] = nodes[i];
  unique_coords.Resize(dim);
  unique_inverse.Resize(dim);
  basis.Resize(dim);
  for(unsigned int i = 0; i < dim; ++i) {
    // Fastest way of unique in C++ is using unordered_set
    std::unordered_set<double> s;
    for(unsigned int j = 0; j < num_nodes; ++j) {
      if(inside[j]) { // we are only interested in nodes inside the splinebox
        s.insert(nodes[i][j]); // value is only inserted if not already present
      }
    }
    StdVector<double> unique_coordsd;
    unique_coordsd.Reserve(num_nodes);
    for(std::unordered_set<double>::iterator it = s.begin(); it != s.end(); ++it) {
      unique_coordsd.Push_back(*it);
    }
    unique_coordsd.Trim();
    std::sort(unique_coordsd.Begin(), unique_coordsd.End()); //todo: might be unnecessary
    unique_coords.Push_back(&unique_coordsd);

    // Find inverse indices
    unique_inverse[i]->Reserve(num_nodes);
    for(unsigned int j = 0; j < num_nodes; ++j) {
      for(unsigned int index = 0; index <= unique_coordsd.GetSize(); i++) {
        if(unique_coordsd[index] == nodes[i][j]) {
          (*unique_inverse[i])[j] = index;
          break;
        }
      }
    }

    // evaluate basis to get spline coefficients
    basis.Push_back(bsplines[i]->EvalBasis(&unique_coordsd));

    LOG_DBG3(SplBox) << "GB: unique_coords[" << dim << "]: \n" << unique_coordsd.ToString();
    LOG_DBG3(SplBox) << "GB: unique_inverse[" << dim << "]: \n" << unique_inverse[i]->ToString();
    LOG_DBG3(SplBox) << "GB: basis[" << dim << "]: \n" << basis[i].ToString();
  }
}

Matrix<double> SplineBox::GetBoxMatrix() {
  // Get number of nodes inside splinebox
  unsigned int num_nodes_inside = 0;
  for(StdVector<bool>::iterator it = inside.Begin(); it != inside.End(); ++it) {
    num_nodes_inside += *it;
  }
  LOG_DBG(SplBox) << "num_nodes_inside: " << num_nodes_inside;

  StdVector<int> sub(dim);
  Matrix<double> mtx;
  mtx.Resize(num_nodes_inside, total_num_cp);
  mtx.InitValue(1.0);

  for(unsigned int d = 0; d < dim; ++d) {
    Matrix<double> mtx_temp;
    mtx_temp.Resize(num_nodes_inside, total_num_cp);
    for(unsigned int i = 0; i < num_nodes_inside; ++i) {
      for(unsigned int j = 0; j < total_num_cp; ++j) {
        ind2sub(num_cp, j, sub);
        mtx_temp[i][j] = basis[d][(*unique_inverse[d])[i]][sub[d]];
      }
    }
    mtx.Mult(mtx_temp, mtx);
  }

  LOG_DBG3(SplBox) << "GBM: box matrix: \n" << mtx.ToString(2);
  return mtx;
}

} // end of namespace
