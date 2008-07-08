#include <iostream>
#include <cmath>
#include <iterator>

#include "SphericalShellGenerator.hh"

namespace CoupledField
{
  SphericalShellGenerator::SphericalShellGenerator() :
    numLevels_(4),
    numInnerPoints_(0),
    innerRadius_(0.58),
    outerRadius_(0.6),
    numLayers_(2),
    maxNumElemNodes_(8)
  {
  }

  int SphericalShellGenerator::CalcHashVal(double x, double y, double z) 
  {
    double hashVal = (x+3)*1e8 + (y+3)*1e5 + (z+3)*1e2;

    //    std::cout << "CalcHashVal " << x << " " << y << " " << z << " " << (int) hashVal << std::endl;
  
    return (int) hashVal;
  }

  void SphericalShellGenerator::RefineQuad(int elemIdx,
                  std::vector<UInt>& connectOut,
                  std::vector<UInt>& regionsOut,
                  std::vector<UInt>& elemTypesOut) 
  {
    double x1,y1,z1;
    double x2,y2,z2;
    double x3,y3,z3;
    double x4,y4,z4;
    double x5,y5,z5;
    double x6,y6,z6;
    double x7,y7,z7;
    double x8,y8,z8;
    double x9,y9,z9;
    int hashVal;
    int idx1, idx2, idx3, idx4, idx5, idx6, idx7, idx8, idx9;

    idx1 = connect_[elemIdx * maxNumElemNodes_+0];  
    x1 = coords_[(idx1 - 1) * 3+0];
    y1 = coords_[(idx1 - 1) * 3+1];
    z1 = coords_[(idx1 - 1) * 3+2];

    idx2 = connect_[elemIdx * maxNumElemNodes_+1];  
    x2 = coords_[(idx2 - 1) * 3+0];
    y2 = coords_[(idx2 - 1) * 3+1];
    z2 = coords_[(idx2 - 1) * 3+2];

    idx3 = connect_[elemIdx * maxNumElemNodes_+2];  
    x3 = coords_[(idx3 - 1) * 3+0];
    y3 = coords_[(idx3 - 1) * 3+1];
    z3 = coords_[(idx3 - 1) * 3+2];

    idx4 = connect_[elemIdx * maxNumElemNodes_+3];  
    x4 = coords_[(idx4 - 1) * 3+0];
    y4 = coords_[(idx4 - 1) * 3+1];
    z4 = coords_[(idx4 - 1) * 3+2];


    x5 = x1 + (x2 - x1) / 2.0;
    y5 = y1 + (y2 - y1) / 2.0;
    z5 = z1 + (z2 - z1) / 2.0;

    hashVal = CalcHashVal(x5, y5, z5);

    if(hashMap_.find(hashVal) == hashMap_.end())
    {
      coords_.push_back(x5);
      coords_.push_back(y5);
      coords_.push_back(z5);

      hashMap_[hashVal] = coords_.size() / 3;
    }
  
    idx5 = hashMap_[hashVal];



    x6 = x2 + (x3 - x2) / 2.0;
    y6 = y2 + (y3 - y2) / 2.0;
    z6 = z2 + (z3 - z2) / 2.0;

    hashVal = CalcHashVal(x6, y6, z6);

    if(hashMap_.find(hashVal) == hashMap_.end())
    {
      coords_.push_back(x6);
      coords_.push_back(y6);
      coords_.push_back(z6);

      hashMap_[hashVal] = coords_.size() / 3;
    }

    idx6 = hashMap_[hashVal];
  


    x7 = x3 + (x4 - x3) / 2.0;
    y7 = y3 + (y4 - y3) / 2.0;
    z7 = z3 + (z4 - z3) / 2.0;

    hashVal = CalcHashVal(x7, y7, z7);

    if(hashMap_.find(hashVal) == hashMap_.end()) 
    {
      coords_.push_back(x7);
      coords_.push_back(y7);
      coords_.push_back(z7);

      hashMap_[hashVal] = coords_.size() / 3;
    }
  
    idx7 = hashMap_[hashVal];
  

    x8 = x4 + (x1 - x4) / 2.0;
    y8 = y4 + (y1 - y4) / 2.0;
    z8 = z4 + (z1 - z4) / 2.0;

    hashVal = CalcHashVal(x8, y8, z8);

    if(hashMap_.find(hashVal) == hashMap_.end())
    {
      coords_.push_back(x8);
      coords_.push_back(y8);
      coords_.push_back(z8);

      hashMap_[hashVal] = coords_.size() / 3;
    }
  
    idx8 = hashMap_[hashVal];



    x9 = x2 + (x4 - x2) / 2.0;
    y9 = y2 + (y4 - y2) / 2.0;
    z9 = z2 + (z4 - z2) / 2.0;

    hashVal = CalcHashVal(x9, y9, z9);

    if(hashMap_.find(hashVal) == hashMap_.end())
    {
      coords_.push_back(x9);
      coords_.push_back(y9);
      coords_.push_back(z9);

      hashMap_[hashVal] = coords_.size() / 3;
    }
  
    idx9 = hashMap_[hashVal];
  
    UInt regionIdx = regions_[elemIdx];
    UInt baseIdx;

    baseIdx = 4 * elemIdx + 0;
    connectOut[baseIdx * maxNumElemNodes_+0] = idx1;
    connectOut[baseIdx * maxNumElemNodes_+1] = idx5;
    connectOut[baseIdx * maxNumElemNodes_+2] = idx9;
    connectOut[baseIdx * maxNumElemNodes_+3] = idx8;
    regionsOut[baseIdx] = regionIdx;
    elemTypesOut[baseIdx] = ET_QUAD4;
  
    baseIdx = 4 * elemIdx + 1;
    connectOut[baseIdx * maxNumElemNodes_+0] = idx5;
    connectOut[baseIdx * maxNumElemNodes_+1] = idx2;
    connectOut[baseIdx * maxNumElemNodes_+2] = idx6;
    connectOut[baseIdx * maxNumElemNodes_+3] = idx9;
    regionsOut[baseIdx] = regionIdx;
    elemTypesOut[baseIdx] = ET_QUAD4;
 
    baseIdx = 4 * elemIdx + 2;
    connectOut[baseIdx * maxNumElemNodes_+0] = idx9;
    connectOut[baseIdx * maxNumElemNodes_+1] = idx6;
    connectOut[baseIdx * maxNumElemNodes_+2] = idx3;
    connectOut[baseIdx * maxNumElemNodes_+3] = idx7;
    regionsOut[baseIdx] = regionIdx;
    elemTypesOut[baseIdx] = ET_QUAD4;

    baseIdx = 4 * elemIdx + 3;
    connectOut[baseIdx * maxNumElemNodes_+0] = idx8;
    connectOut[baseIdx * maxNumElemNodes_+1] = idx9;
    connectOut[baseIdx * maxNumElemNodes_+2] = idx7;
    connectOut[baseIdx * maxNumElemNodes_+3] = idx4;
    regionsOut[baseIdx] = regionIdx;
    elemTypesOut[baseIdx] = ET_QUAD4;

  }

  void SphericalShellGenerator::TriangulateQuads(std::vector<UInt>& connectIn,
                                                 std::vector<UInt>& connectOut,
                                                 std::vector<UInt>& regionsIn,
                                                 std::vector<UInt>& regionsOut)
  {
    UInt p1, p2, p3, p4;
    UInt maxNumElemNodes = 4;
    UInt numQuads = connectIn.size() / maxNumElemNodes;
    connectOut.clear();
    regionsOut.clear();
    
    for(UInt i=0; i<numQuads; i++)
    {
      p1 = connectIn[i*maxNumElemNodes+0];
      p2 = connectIn[i*maxNumElemNodes+1];
      p3 = connectIn[i*maxNumElemNodes+2];
      p4 = connectIn[i*maxNumElemNodes+3];
      
      connectOut.push_back(p1);
      connectOut.push_back(p2);
      connectOut.push_back(p3);
      regionsOut.push_back(regionsIn[i]);
      
      connectOut.push_back(p3);
      connectOut.push_back(p4);
      connectOut.push_back(p1);
      regionsOut.push_back(regionsIn[i]);
    }
  }
  

  void SphericalShellGenerator::ProjectCoordsToSphere(std::vector<double>& coordsIn,
                             std::vector<double>& coordsOut,
                             double radius) 
  {
    std::vector<double>::iterator it, end;
    std::vector<double>::iterator it2;
    
    coordsOut.resize(coordsIn.size());
  
    it = coordsIn.begin();
    end = coordsIn.end();

    it2 = coordsOut.begin();
  
    for( ; it != end; it+=3, it2+=3) 
    {
      double x,y,z;
      double factor;
      x = *it;
      y = *(it+1);
      z = *(it+2);

      double length = sqrt(x*x + y*y + z*z);
      factor = (1/length) * radius;
    
      x *= factor;
      y *= factor;
      z *= factor;

      // std::cout << "length: " << length << " " << x << " " << y << " " << z << std::endl;
    
      *it2 = x;
      *(it2+1) = y;
      *(it2+2) = z;
    }
  }

  void SphericalShellGenerator::CreateOuterQuads()
  {
    UInt numElems = regions_.size();
    UInt numRegions = 6;
    UInt nodeOffset = numInnerPoints_ * numLayers_;
    UInt p1, p2, p3, p4;
    UInt region;

    regionNames_.push_back("outer_bottom");
    regionNames_.push_back("outer_front");
    regionNames_.push_back("outer_right");
    regionNames_.push_back("outer_back");
    regionNames_.push_back("outer_left");
    regionNames_.push_back("outer_top");
    
    connect_.resize(numElems * 2 * maxNumElemNodes_);
    
    for(UInt i=0; i<numElems; i++) 
    {
      p1 = connect_[i*maxNumElemNodes_+0];
      p2 = connect_[i*maxNumElemNodes_+1];
      p3 = connect_[i*maxNumElemNodes_+2];
      p4 = connect_[i*maxNumElemNodes_+3];

      region = regions_[i];
      
      connect_[(numElems+i)*maxNumElemNodes_+0] = p1 + nodeOffset;
      connect_[(numElems+i)*maxNumElemNodes_+1] = p2 + nodeOffset;
      connect_[(numElems+i)*maxNumElemNodes_+2] = p3 + nodeOffset;
      connect_[(numElems+i)*maxNumElemNodes_+3] = p4 + nodeOffset;
      
      regions_.push_back(region + numRegions);

      elemTypes_.push_back(ET_QUAD4);
    }
  }

  void SphericalShellGenerator::CreateHexas()
  {
    UInt numInnerElems = regions_.size() / 2;
    UInt numRegions = 12;
    UInt elemOffset;
    UInt elem;
    UInt nodeOffsetBottom;
    UInt nodeOffsetTop;
    UInt p1, p2, p3, p4;
    UInt region;

    regionNames_.push_back("hexas_bottom");
    regionNames_.push_back("hexas_front");
    regionNames_.push_back("hexas_right");
    regionNames_.push_back("hexas_back");
    regionNames_.push_back("hexas_left");
    regionNames_.push_back("hexas_top");
    
    connect_.resize(numInnerElems * (2 + numLayers_) * maxNumElemNodes_);
    elemOffset = numInnerElems * 2;

    elem = elemOffset;
    for(UInt i=0; i<numInnerElems; i++) 
    {
      p1 = connect_[i*maxNumElemNodes_+0];
      p2 = connect_[i*maxNumElemNodes_+1];
      p3 = connect_[i*maxNumElemNodes_+2];
      p4 = connect_[i*maxNumElemNodes_+3];
      region = regions_[i];

      for(UInt layer=0; layer<numLayers_; layer++) 
      {
        nodeOffsetBottom = numInnerPoints_ * layer;
        nodeOffsetTop = numInnerPoints_ * (layer+1);

        connect_[elem*maxNumElemNodes_+0] = p1 + nodeOffsetBottom;
        connect_[elem*maxNumElemNodes_+1] = p2 + nodeOffsetBottom;
        connect_[elem*maxNumElemNodes_+2] = p3 + nodeOffsetBottom;
        connect_[elem*maxNumElemNodes_+3] = p4 + nodeOffsetBottom;

        connect_[elem*maxNumElemNodes_+4] = p1 + nodeOffsetTop;
        connect_[elem*maxNumElemNodes_+5] = p2 + nodeOffsetTop;
        connect_[elem*maxNumElemNodes_+6] = p3 + nodeOffsetTop;
        connect_[elem*maxNumElemNodes_+7] = p4 + nodeOffsetTop;
        
        regions_.push_back(region + numRegions);
        
        elemTypes_.push_back(ET_HEXA8);

        elem++;
      }
    }
  }
  
  void SphericalShellGenerator::CreateBaseMesh()
  {
    coords_.clear();
    connect_.clear();
    hashMap_.clear();
    regions_.clear();
    elemTypes_.clear();
    
    // Unten
    coords_.push_back(-1);
    coords_.push_back(-1);
    coords_.push_back(-1);

    coords_.push_back(1);
    coords_.push_back(-1);
    coords_.push_back(-1);

    coords_.push_back(1);
    coords_.push_back(1);
    coords_.push_back(-1);
  
    coords_.push_back(-1);
    coords_.push_back(1);
    coords_.push_back(-1);

    // oben
    coords_.push_back(-1);
    coords_.push_back(-1);
    coords_.push_back(1);

    coords_.push_back(1);
    coords_.push_back(-1);
    coords_.push_back(1);

    coords_.push_back(1);
    coords_.push_back(1);
    coords_.push_back(1);
  
    coords_.push_back(-1);
    coords_.push_back(1);
    coords_.push_back(1);

    for(int i=0, n=coords_.size()/3; i<n; i++)
    {
      double x = coords_[i*3+0];
      double y = coords_[i*3+1];
      double z = coords_[i*3+2];
    
      unsigned int hashVal = CalcHashVal(x, y, z);
    
      hashMap_[hashVal] = i+1;
    }

    connect_.resize(6 * maxNumElemNodes_);

    connect_[0 * maxNumElemNodes_ + 0] = 1;
    connect_[0 * maxNumElemNodes_ + 1] = 2;
    connect_[0 * maxNumElemNodes_ + 2] = 3;
    connect_[0 * maxNumElemNodes_ + 3] = 4;

    regions_.push_back(0);
    elemTypes_.push_back(ET_QUAD4);

    connect_[1 * maxNumElemNodes_ + 0] = 1;
    connect_[1 * maxNumElemNodes_ + 1] = 2;
    connect_[1 * maxNumElemNodes_ + 2] = 6;
    connect_[1 * maxNumElemNodes_ + 3] = 5;

    regions_.push_back(1);
    elemTypes_.push_back(ET_QUAD4);
  
    connect_[2 * maxNumElemNodes_ + 0] = 2;
    connect_[2 * maxNumElemNodes_ + 1] = 3;
    connect_[2 * maxNumElemNodes_ + 2] = 7;
    connect_[2 * maxNumElemNodes_ + 3] = 6;

    regions_.push_back(2);
    elemTypes_.push_back(ET_QUAD4);

    connect_[3 * maxNumElemNodes_ + 0] = 3;
    connect_[3 * maxNumElemNodes_ + 1] = 4;
    connect_[3 * maxNumElemNodes_ + 2] = 8;
    connect_[3 * maxNumElemNodes_ + 3] = 7;

    regions_.push_back(3);
    elemTypes_.push_back(ET_QUAD4);

    connect_[4 * maxNumElemNodes_ + 0] = 4;
    connect_[4 * maxNumElemNodes_ + 1] = 1;
    connect_[4 * maxNumElemNodes_ + 2] = 5;
    connect_[4 * maxNumElemNodes_ + 3] = 8;

    regions_.push_back(4);
    elemTypes_.push_back(ET_QUAD4);

    connect_[5 * maxNumElemNodes_ + 0] = 5;
    connect_[5 * maxNumElemNodes_ + 1] = 6;
    connect_[5 * maxNumElemNodes_ + 2] = 7;
    connect_[5 * maxNumElemNodes_ + 3] = 8;

    regions_.push_back(5);
    elemTypes_.push_back(ET_QUAD4);

    regionNames_.push_back("inner_bottom");
    regionNames_.push_back("inner_front");
    regionNames_.push_back("inner_right");
    regionNames_.push_back("inner_back");
    regionNames_.push_back("inner_left");
    regionNames_.push_back("inner_top");
  }

  void SphericalShellGenerator::RefineInnerQuads()
  {
    std::vector<UInt> connect2;
    std::vector<UInt> regions2;
    std::vector<UInt> elemTypes2;
    UInt numElems;
    
    for(UInt level=0; level < numLevels_; level++) 
    {
      numElems = elemTypes_.size();
      connect2.resize(numElems * 4 * maxNumElemNodes_);
      std::fill(connect2.begin(), connect2.end(), 0);
      regions2.resize(numElems * 4);
      std::fill(regions2.begin(), regions2.end(), 0);
      elemTypes2.resize(numElems * 4);
      std::fill(elemTypes2.begin(), elemTypes2.end(), 0);

      for(UInt i=0; i<numElems; i++) 
      {
        RefineQuad(i, connect2, regions2, elemTypes2);
      }

      connect_ = connect2;
      regions_ = regions2;
      elemTypes_ = elemTypes2;
    }
  }

  void SphericalShellGenerator::GenSphericalShell(std::vector<double>& coords,
                         std::vector<UInt>& connect,
                         std::vector<UInt>& elemTypes,
                         UInt& maxNumElemNodes,
                         std::map<std::string, std::vector<UInt> >& regionElems)
  {
    coords.clear();
    connect.clear();
    elemTypes.clear();
  

    CreateBaseMesh();
    
    RefineInnerQuads();
    
#if 0
    TriangulateQuads(connect, connect2, regions, regions2);
    maxNumElemNodes = 3;
    elemType = ET_TRIA3;
    
    connect = connect2;
    connect2.clear();

    regions = regions2;
    regions2.clear();
#endif

    numInnerPoints_ = coords_.size() / 3;
    
    Double deltaRadius = (outerRadius_ - innerRadius_) / numLayers_;
    std::vector<Double> myCoords, myCoords2;
    
    myCoords = coords_;
    coords_.clear();
    for(UInt layer=0; layer <= numLayers_; layer++) 
    {
      ProjectCoordsToSphere(myCoords, myCoords2, innerRadius_ + layer * deltaRadius);
      std::copy(myCoords2.begin(), myCoords2.end(), std::back_inserter(coords_));
    }

    CreateOuterQuads();

    CreateHexas();

    UInt numElems = connect_.size() / maxNumElemNodes_;
    std::cout << "NumElems: " << numElems << std::endl;
    std::cout << "NumInnerNodes: " << numInnerPoints_ << std::endl;

    //    UInt k = 1;
    //    UInt numFaceElems = numElems/6;
    //    std::map<std::string, std::vector<UInt> >::iterator it;
    
    for(UInt i=0; i<numElems; i++) 
    {
      UInt regionIdx = regions_[i];
      regionElems[regionNames_[regionIdx]].push_back(i+1);
    }
    
    elemTypes = elemTypes_;
    connect = connect_;
    coords = coords_;
    maxNumElemNodes = maxNumElemNodes_;


    /*    for(int i=0; i<numFaceElems; i++)
    {
      regionElems["inner_bottom"].push_back(k);
      k++;
    }

    for(int i=0; i<numFaceElems; i++)
    {
      regionElems["inner_front"].push_back(k);
      k++;
    }

    for(int i=0; i<numFaceElems; i++)
    {
      regionElems["inner_right"].push_back(k);
      k++;
    }

    for(int i=0; i<numFaceElems; i++)
    {
      regionElems["inner_back"].push_back(k);
      k++;
    }
  
    for(int i=0; i<numFaceElems; i++)
    {
      regionElems["inner_left"].push_back(k);
      k++;
    }

    for(int i=0; i<numFaceElems; i++)
    {
      regionElems["inner_top"].push_back(k);
      k++;
    }
    */
  }
}
