#ifndef GBEDGECOLLAPSE_HH
#define GBEDGECOLLAPSE_HH

#include<GbVec3.hh>


struct GbHalfEdgeCollapse {
  int i;                           // origin
  int j;                           // new 
  std::vector<int> n;              // neighbours
  int b;                           // is it a boundary
  GbVec3<float> pos;               // position
  int del;
  float energy;                    // for flow vertex infos
  float density;
  GbVec3<float> momentum;
};



struct GbEdgeCollapse {
  int i;                           // origin
  int j;                           // new 
  std::vector<int> n;              // neighbours
  int b;                           // is it a boundary
  GbVec3<float> pos;               // position of new vertex
  GbVec3<float> pos2;              // position of old vertex
};


#endif // GBEDGECOLLAPSE_HH
