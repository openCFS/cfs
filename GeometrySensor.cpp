/***************************************************************************
    File        : GeometrySensor.cpp
    Description : 

 ---------------------------------------------------------------------------
    Begin       : Thu Apr 25 2002
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#include "GeometrySensor.h"

using namespace std;
namespace grd {

class SetNormals {
public:
  SetNormals () {}
  ~SetNormals () {}

  // Methods
  void operator() (Element* t)
  {
    TriSkeleton* triSk = ((Triangle*)t)->getSkeleton();
    if (triSk == 0)
    {
      curvature.push_back(1.0);
    }
    else if (!triSk->isRefined()) {
      curvature.push_back(1.0);
    }
    else {
      int i,ch;
      double len0,len1;
      double normal[3];
      double chnorm[4][3];
      double v1[3],v2[3];

      for (i = 0; i < 3; i++)
      {
        v1[i] = triSk->pos[1][i] - triSk->pos[0][i];
        v2[i] = triSk->pos[2][i] - triSk->pos[0][i];
      }
      crossProduct(normal,v1,v2);
      for (ch = 0; ch < 4; ch++)
      {
        for (i = 0; i < 3; i++)
        {
          v1[i] = triSk->child[ch]->pos[1][i] - triSk->child[ch]->pos[0][i];
          v2[i] = triSk->child[ch]->pos[2][i] - triSk->child[ch]->pos[0][i];
        }
        crossProduct(chnorm[ch],v1,v2);
      }

      len0 = grd::length(normal);
      double curv = 1.0;
      double tmp;
      for (ch = 0; ch < 4; ch++)
      {
        len1 = grd::length(chnorm[ch]);
        tmp = grd::innerProduct(normal,chnorm[ch]);
        tmp /= (len0*len1);
        tmp = fabs(tmp);
        curv = (curv < tmp) ? curv : tmp;
      }

      curvature.push_back(curv);
    }

  }

  float getMin() {
    int sz = curvature.size();
    if (sz > 0)
    {
      double min = curvature[0];
      for (int i = 0; i < sz; i++)
        min = (min < curvature[i]) ? min : curvature[i];
      return (float) min;
    }
    else
      return (float) 0.0;
  }
  float getMax() {
    int sz = curvature.size();
    if (sz > 0)
    {
      double max = curvature[0];
      for (int i = 0; i < sz; i++)
        max = (max > curvature[i]) ? max : curvature[i];
      return (float) max;
    }
    else
      return (float) 0.0;
  }

  vector<float>* getCurvature() { return &curvature; }
private:
  vector<float> curvature;
};



class MarkElement {
public:
  MarkElement(vector<float>* curv,float min,float max,float thr) {
    curvature = curv;
    minCurv = min;
    maxCurv = max;
    threshold = thr;
    limes = threshold; // maxCurv*threshold + (1.0 - threshold)*minCurv;
    counter = 0;
  }
  ~MarkElement() {}

  // Method
  void operator() (Element* t) {
    float value = (*curvature)[counter];
    if (value < limes)
      t->markForRefinement();
    counter++;
  }

private:
  int counter;
  float minCurv;
  float maxCurv;
  float threshold;
  float limes;
  vector<float>* curvature;
};


//
//
//
void
GeometrySensor::markForRefinement(MultilevelGrid& grid)
{
  // for each triangle mark
  SetNormals sn;
  grid.forEachUnrefinedBoundaryElement(sn);

  // Mark elements
  float minCurv = sn.getMin();
  float maxCurv = sn.getMax();
  vector<float>* curvs = sn.getCurvature();
  MarkElement  me(curvs,minCurv,maxCurv,0.99);
  grid.forEachUnrefinedBoundaryElement(me);
}

} // namespace