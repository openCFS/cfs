#include "DataInOut/Logging/LogConfigurator.hh"
#include "MatVec/Matrix.hh"
#include "Utils/BSpline.hh"

namespace CoupledField {

DEFINE_LOG(BSpl, "BSpline")

DEFINE_LOG(BSplC, "BSplineCurve")

BSpline::BSpline(const unsigned int degree, const unsigned int numKnots, double min_knot, double max_knot) {
  this->degree_ = degree;
  this->numKnots_ = numKnots;
  MakeKnotVector(min_knot, max_knot);
}

Matrix<double> BSpline::Eval(StdVector<double>* t) {
  if(t == NULL) {
    StdVector<double> tmp(100);
    for(unsigned int i = 0; i < 100; ++i)
      tmp[i] = i*0.01;

    tmp[99] -= 5*std::numeric_limits<double>::epsilon();

    return BasisFuncDg(degree_, &tmp, knots_, numKnots_ - degree_ - 1);
  }

  for(StdVector<double>::iterator it = t->Begin(); it != t->End(); ++it) {
    if(*it == knots_.Last())
      *it -= 5*std::numeric_limits<double>::epsilon();
  }
  return BasisFuncDg(degree_, t, knots_, numKnots_ - degree_ - 1);
}

Matrix<double> BSpline::Eval(double t) {
  StdVector<double> temp(1);
  temp[0] = t;
  return Eval(&temp);
}

void BSpline::MakeKnotVector(double min_knot, double max_knot) {
  unsigned int numInnerKnots = numKnots_ - 2 * degree_;
  double delta = (max_knot-min_knot)/(numInnerKnots - 1);

  StdVector<double> knots;
  knots.Reserve(numKnots_);
  for(unsigned int i = 0; i < degree_; ++i) {
    knots.Push_back(min_knot);
  }
  for(unsigned int i = 0; i < numInnerKnots; ++i) {
    knots.Push_back(min_knot + i*delta);
  }
  for(unsigned int i = 0; i < degree_; ++i) {
    knots.Push_back(max_knot);
  }
  this->knots_ = knots;
  LOG_DBG(BSpl) << "MKV: knots: " << knots.ToString();
}

StdVector<double> BSpline::GrevilleAbscissae() {
  unsigned int numGA = numKnots_ - degree_ - 1;
  StdVector<double> GAbscissae = StdVector<double>(numGA);

  StdVector<double> knots = GetKnots();
  for(unsigned int j = 0; j < numGA; ++j) {
    for(unsigned int k = 0; k < degree_; ++k) {
      GAbscissae[j] += knots[j+1+k]/degree_;
    }
  }

  LOG_DBG(BSpl) << "GA: Greville abscissae: " << GAbscissae.ToString();
  return GAbscissae;
}

Matrix<double> BSpline::BasisFuncDg0(StdVector<double>* t, StdVector<double> knots, unsigned int numIntervals) {

  Matrix<double> bfun;
  bfun.Resize(t->GetSize(), numIntervals);
  bfun.Init();

  // Recursive BSpline definition by Cox-de Boor:
  // B_i,0(t) = 1 if knots_i <= t < knots_i+1 else 0
  for(unsigned int i = 0; i < numIntervals; ++i) {
    for(unsigned int j = 0; j < t->GetSize(); ++j) {
      if(knots[i] <= (*t)[j] && (*t)[j] < knots[i+1]) {
        bfun[j][i] = 1;
      }
    }
  }

  LOG_DBG(BSpl) << "BFD0: bfun: \n" << bfun.ToString();
  return bfun;
}

Matrix<double> BSpline::BasisFuncDg(unsigned int degree, StdVector<double>* t, StdVector<double> knots, unsigned int numIntervals) {
  Matrix<double> bfun;
  bfun.Resize(t->GetSize(), numIntervals);
  bfun.Init();

  // Recursive BSpline definition by Cox-de Boor:
  // B_i,p(t) = (t - knots_i)/(knots_i+p - knots_i) B_i,p-1(t) + (knots_i+p+1 - t)/(knots_i+p+1 - knots_i+1) B_i+1,p-1(t)
  if(degree >= 1) {
    Matrix<double> bfun_dgm1 = BasisFuncDg(degree-1, t, knots, numIntervals + 1);
    double c1, c2;
    for(unsigned int i = 0; i < numIntervals; ++i) {
      c2 = knots[i + degree] - knots[i];
      if(std::abs(c2) > std::numeric_limits<double>::epsilon()) {
        for(unsigned int j = 0; j < t->GetSize(); ++j) {
          c1 = (*t)[j] - knots[i];
          bfun[j][i] = c1 / c2 * bfun_dgm1[j][i];
        }
      }

      c2 = knots[i + degree + 1] - knots[i+1];
      if(std::abs(c2) > std::numeric_limits<double>::epsilon()) {
        for(unsigned int j = 0; j < t->GetSize(); ++j) {
          c1 = knots[i + degree + 1] - (*t)[j];
          bfun[j][i] += c1 / c2 * bfun_dgm1[j][i+1];
        }
      }

    }
  } else {
    bfun = BSpline::BasisFuncDg0(t, knots, numIntervals);
  }

  LOG_DBG(BSpl) << "BFD: bfun: \n" << bfun.ToString();
  return bfun;
}


BSplineCurve::BSplineCurve(unsigned int degree, unsigned int numControlPoints) {
  numCP_ = numControlPoints;
  bspline_ = new BSpline(degree, numControlPoints + degree + 1);

  StdVector<double> ga = bspline_->GrevilleAbscissae();
  control_points_.Reserve(ga.GetSize());
  for(unsigned int i = 0; i < ga.GetSize(); ++i) {
    Point p(0,0,0);
    p[0] = ga[i];
    control_points_.Push_back(p);
  }
}

void BSplineCurve::SetControlPoint(int index, Point coords) {
  assert(index < (int) control_points_.GetSize());
  control_points_[index] = coords;
}

void BSplineCurve::SetControlPoints(StdVector<Point> coords) {
  assert(numCP_ == coords.GetSize());
  control_points_ = coords;
}

Matrix<double> BSplineCurve::Eval(StdVector<double>* t) {
  assert(!control_points_.IsEmpty());

  Matrix<double> basis = bspline_->Eval(t);
  assert (basis.GetNumCols() == control_points_.GetSize());
  assert (basis.GetNumRows() == t->GetSize());

  unsigned int dim = 3; // Point is always 3D

  Matrix<double> curve_coors;
  curve_coors.Resize(basis.GetNumRows(), dim);
  curve_coors.Init();

  for(unsigned int i = 0; i < dim; ++i) {
    for(unsigned int j = 0; j < basis.GetNumRows(); ++j) {
      curve_coors[j][i] = basis[j][0] * control_points_[0].data[i];
      for(unsigned int k = 1; k < control_points_.GetSize(); ++k) {
        curve_coors[j][i] += basis[j][k] * control_points_[k].data[i];
      }
    }
  }

  LOG_DBG(BSplC) << "Eval: curve corrdinates: \n" << curve_coors.ToString();
  return curve_coors;
}

} // end of namespace
