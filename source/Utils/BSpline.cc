#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"
#include "MatVec/Matrix.hh"
#include "Utils/BSpline.hh"

namespace CoupledField {

DECLARE_LOG(BSpl)
DEFINE_LOG(BSpl, "BSpline")

DECLARE_LOG(BSplC)
DEFINE_LOG(BSplC, "BSplineCurve")

BSpline::BSpline(const unsigned int degree, const unsigned int numKnots, StdVector<double> knot_range) {
  this->degree = degree;
  this->numKnots = numKnots;
  MakeKnotVector(knot_range);
};

Matrix<double> BSpline::EvalBasis(StdVector<double>* t) {
  if(t == NULL) {
    StdVector<double> t;
    t.Reserve(100);
    for(unsigned int i = 0; i < 100; ++i) {
      t[i] = i*0.01;
    }
    t[99] -= 5*std::numeric_limits<double>::epsilon();

    return BasisFuncDg(degree, &t, knots, numKnots - degree - 1);
  }

  for(StdVector<double>::iterator it = t->Begin(); it != t->End(); ++it) {
    if(*it == knots.Last()) {
      *it -= 5*std::numeric_limits<double>::epsilon();
    }
  }
  return BasisFuncDg(degree, t, knots, numKnots - degree - 1);
}

void BSpline::MakeKnotVector(StdVector<double> knot_range) {
  double min_knot = 0.0;
  double max_knot = 1.0;
  if(!knot_range.IsEmpty()) {
    min_knot = knot_range[0];
    max_knot = knot_range[1];
  }
  unsigned int numInnerKnots = numKnots - 2 * degree;
  double delta = (max_knot-min_knot)/(numInnerKnots - 1);

  StdVector<double> knots;
  knots.Reserve(numKnots);
  for(unsigned int i = 0; i < degree; ++i) {
    knots.Push_back(min_knot);
  }
  for(unsigned int i = 0; i < numInnerKnots; ++i) {
    knots.Push_back(min_knot + i*delta);
  }
  for(unsigned int i = 0; i < degree; ++i) {
    knots.Push_back(max_knot);
  }
  this->knots = knots;
  LOG_DBG(BSpl) << "MKV: knots: " << knots.ToString();
}

StdVector<double> BSpline::GrevilleAbscissae() {
  unsigned int numGA = numKnots - degree - 1;
  StdVector<double> GAbscissae = StdVector<double>(numGA);

  StdVector<double> knots = GetKnots();
  for(unsigned int j = 0; j < numGA; ++j) {
    for(unsigned int k = 0; k < degree; ++k) {
      GAbscissae[j] += knots[j+1+k]/degree;
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
      if(knots[i] <= (*t)[j] and (*t)[j] < knots[i+1]) {
        bfun[j][i] = 1;
      }
    }
  }

  LOG_DBG(BSpl) << "BFD0: bfun: \n" << bfun.ToString(2);
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

  LOG_DBG(BSpl) << "BFD: bfun: \n" << bfun.ToString(2);
  return bfun;
}


BSplineCurve::BSplineCurve(unsigned int degree, unsigned int numControlPoints) {
  numCP = numControlPoints;
  bspline = new BSpline(degree, numControlPoints + degree + 1);

  StdVector<double> ga = bspline->GrevilleAbscissae();
  control_points.Reserve(ga.GetSize());
  for(unsigned int i = 0; i < ga.GetSize(); ++i) {
    Point p(0,0,0);
    p[0] = ga[i];
    control_points.Push_back(p);
  }
}

void BSplineCurve::SetControlPoint(int index, Point coords) {
  assert (index < control_points.GetSize());
  control_points[index] = coords;
}

void BSplineCurve::SetControlPoints(StdVector<Point> coords) {
  assert (numCP == coords.GetSize());
  control_points = coords;
}

Matrix<double> BSplineCurve::Eval(StdVector<double>* t) {
  assert(!control_points.IsEmpty());

  Matrix<double> basis = bspline->EvalBasis(t);
  assert (basis.GetNumCols() == control_points.GetSize());
  assert (basis.GetNumRows() == t->GetSize());

  unsigned int dim = 3; // Point is always 3D

  Matrix<double> curve_coors;
  curve_coors.Resize(basis.GetNumRows(), dim);
  curve_coors.Init();

  for(unsigned int i = 0; i < dim; ++i) {
    for(unsigned int j = 0; j < basis.GetNumRows(); ++j) {
      curve_coors[j][i] = basis[j][0] * control_points[0].data[i];
      for(unsigned int k = 1; k < control_points.GetSize(); ++k) {
        curve_coors[j][i] += basis[j][k] * control_points[k].data[i];
      }
    }
  }

  LOG_DBG(BSplC) << "Eval: curve corrdinates: \n" << curve_coors.ToString(2);
  return curve_coors;
}

} // end of namespace
