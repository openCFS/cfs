#include <iostream>
#include <cmath>
#include <iterator>

#include "SphericalShellGenerator.hh"

namespace CoupledField
{
  SphericalShellGenerator::SphericalShellGenerator() :
    numLevels_(3),
    numInnerPoints_(0),
    numInnerElems_(0),
    innerRadius_(0.55),
    outerRadius_(0.6),
    numLayers_(2),
    maxNumElemNodes_(27),
    octantRegions_(false),
    makeQuadratic_(true),
    makeSuperQuadratic_(false)
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
                                           std::vector<UInt>& elemTypesOut,
                                           bool generateElems)
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
  
    if(!generateElems)
    {
      quadConnectFull_[elemIdx*13+0] = idx1;
      quadConnectFull_[elemIdx*13+1] = idx2;
      quadConnectFull_[elemIdx*13+2] = idx3;
      quadConnectFull_[elemIdx*13+3] = idx4;
      quadConnectFull_[elemIdx*13+4] = idx5;
      quadConnectFull_[elemIdx*13+5] = idx6;
      quadConnectFull_[elemIdx*13+6] = idx7;
      quadConnectFull_[elemIdx*13+7] = idx8;
      quadConnectFull_[elemIdx*13+8] = idx9;

      Double x,y,z;
      x = x1 + (x9 - x1) / 2.0;
      y = y1 + (y9 - y1) / 2.0;
      z = z1 + (z9 - z1) / 2.0;
      
      hashVal = CalcHashVal(x, y, z);
      
      if(hashMap_.find(hashVal) == hashMap_.end())
      {
        coords_.push_back(x);
        coords_.push_back(y);
        coords_.push_back(z);
        
        hashMap_[hashVal] = coords_.size() / 3;
      }

      UInt idx10 = hashMap_[hashVal];

      x = x2 + (x9 - x2) / 2.0;
      y = y2 + (y9 - y2) / 2.0;
      z = z2 + (z9 - z2) / 2.0;
      
      hashVal = CalcHashVal(x, y, z);
      
      if(hashMap_.find(hashVal) == hashMap_.end())
      {
        coords_.push_back(x);
        coords_.push_back(y);
        coords_.push_back(z);
        
        hashMap_[hashVal] = coords_.size() / 3;
      }

      UInt idx11 = hashMap_[hashVal];

      x = x3 + (x9 - x3) / 2.0;
      y = y3 + (y9 - y3) / 2.0;
      z = z3 + (z9 - z3) / 2.0;
      
      hashVal = CalcHashVal(x, y, z);
      
      if(hashMap_.find(hashVal) == hashMap_.end())
      {
        coords_.push_back(x);
        coords_.push_back(y);
        coords_.push_back(z);
        
        hashMap_[hashVal] = coords_.size() / 3;
      }

      UInt idx12 = hashMap_[hashVal];

      x = x4 + (x9 - x4) / 2.0;
      y = y4 + (y9 - y4) / 2.0;
      z = z4 + (z9 - z4) / 2.0;
      
      hashVal = CalcHashVal(x, y, z);
      
      if(hashMap_.find(hashVal) == hashMap_.end())
      {
        coords_.push_back(x);
        coords_.push_back(y);
        coords_.push_back(z);
        
        hashMap_[hashVal] = coords_.size() / 3;
      }
      
      UInt idx13 = hashMap_[hashVal];

      quadConnectFull_[elemIdx*13+9] = idx10;
      quadConnectFull_[elemIdx*13+10] = idx11;
      quadConnectFull_[elemIdx*13+11] = idx12;
      quadConnectFull_[elemIdx*13+12] = idx13;
      
      return;
    }
    
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
    UInt numRegions;
    UInt nodeOffset = numInnerPoints_ * 4 * numLayers_;
    UInt p1, p2, p3, p4;
    UInt p1_outer, p2_outer, p3_outer, p4_outer;
    UInt region;
    UInt baseIdx;

    numRegions=regionBaseNames_.size();
    
    for(UInt i=0; i < numRegions; i++)
    {
      std::stringstream sstr;
      
      sstr << "outer_" << regionBaseNames_[i];
      regionNames_.push_back(sstr.str());
    }

    connect_.resize(numElems * 2 * maxNumElemNodes_);
    
    for(UInt i=0; i<numElems; i++) 
    {
      baseIdx = i*maxNumElemNodes_;
      p1 = connect_[baseIdx + 0];
      p2 = connect_[baseIdx + 1];
      p3 = connect_[baseIdx + 2];
      p4 = connect_[baseIdx + 3];

      region = regions_[i];
      
      baseIdx = (numElems+i)*maxNumElemNodes_;
      p1_outer = p1 + nodeOffset;
      p2_outer = p2 + nodeOffset;
      p3_outer = p3 + nodeOffset;;
      p4_outer = p4 + nodeOffset;;
      
      connect_[baseIdx + 0] = p1_outer;
      connect_[baseIdx + 1] = p2_outer;
      connect_[baseIdx + 2] = p3_outer;
      connect_[baseIdx + 3] = p4_outer;
      
      regions_.push_back(region + numRegions);

      elemTypes_.push_back(ET_QUAD4);
    }

    nodeGroups_["outer_face_bottom_pole"].push_back(25 + nodeOffset);
    nodeGroups_["outer_face_top_pole"].push_back(26 + nodeOffset);
    nodeGroups_["outer_face_front_pole"].push_back(21 + nodeOffset);
    nodeGroups_["outer_face_right_pole"].push_back(22 + nodeOffset);
    nodeGroups_["outer_face_back_pole"].push_back(23 + nodeOffset);
    nodeGroups_["outer_face_left_pole"].push_back(24 + nodeOffset);

    nodeGroups_["outer_edge_nine_pole"].push_back(17 + nodeOffset);
    nodeGroups_["outer_edge_ten_pole"].push_back(18 + nodeOffset);
    nodeGroups_["outer_edge_eleven_pole"].push_back(19 + nodeOffset);
    nodeGroups_["outer_edge_twelve_pole"].push_back(20 + nodeOffset);

    nodeGroups_["outer_corner_one"].push_back(1 + nodeOffset);
    nodeGroups_["outer_corner_two"].push_back(2 + nodeOffset);
    nodeGroups_["outer_corner_three"].push_back(3 + nodeOffset);
    nodeGroups_["outer_corner_four"].push_back(4 + nodeOffset);
    nodeGroups_["outer_corner_five"].push_back(5 + nodeOffset);
    nodeGroups_["outer_corner_six"].push_back(6 + nodeOffset);
    nodeGroups_["outer_corner_seven"].push_back(7 + nodeOffset);
    nodeGroups_["outer_corner_eight"].push_back(8 + nodeOffset);

  }

  void SphericalShellGenerator::CreateHexas()
  {
    UInt numRegions;
    UInt elemOffset;
    UInt elem;
    UInt nodeOffsetBottom;
    UInt nodeOffsetMiddle;
    UInt nodeOffsetTop;
    UInt p1, p2, p3, p4, p5, p6, p7, p8, p9;
    UInt region;

    numRegions = regionBaseNames_.size();
    
    for(UInt i=0; i < numRegions; i++)
    {
      std::stringstream sstr;
      
      sstr << "volelems_" << regionBaseNames_[i];
      regionNames_.push_back(sstr.str());
    }
    
    numRegions *= 2;

    connect_.resize(numInnerElems_ * (2 + numLayers_) * maxNumElemNodes_);
    elemOffset = numInnerElems_ * 2;

    elem = elemOffset;
    for(UInt i=0; i<numInnerElems_; i++) 
    {
      p1 = quadConnectFull_[i*13+0];
      p2 = quadConnectFull_[i*13+1];
      p3 = quadConnectFull_[i*13+2];
      p4 = quadConnectFull_[i*13+3];

      if(makeQuadratic_)
      {
        p5 = quadConnectFull_[i*13+4];
        p6 = quadConnectFull_[i*13+5];
        p7 = quadConnectFull_[i*13+6];
        p8 = quadConnectFull_[i*13+7];

        if(makeSuperQuadratic_)
        {
          p9 = quadConnectFull_[i*13+8];
        }
      }
      
      region = regions_[i];

      for(UInt layer=0; layer<numLayers_; layer++) 
      {
        nodeOffsetBottom = numInnerPoints_ * 4 * layer;
        nodeOffsetMiddle = numInnerPoints_ * 4 * layer + 2 * numInnerPoints_;
        nodeOffsetTop = numInnerPoints_ * 4 * (layer+1);

        connect_[elem*maxNumElemNodes_+0] = p1 + nodeOffsetBottom;
        connect_[elem*maxNumElemNodes_+1] = p2 + nodeOffsetBottom;
        connect_[elem*maxNumElemNodes_+2] = p3 + nodeOffsetBottom;
        connect_[elem*maxNumElemNodes_+3] = p4 + nodeOffsetBottom;

        connect_[elem*maxNumElemNodes_+4] = p1 + nodeOffsetTop;
        connect_[elem*maxNumElemNodes_+5] = p2 + nodeOffsetTop;
        connect_[elem*maxNumElemNodes_+6] = p3 + nodeOffsetTop;
        connect_[elem*maxNumElemNodes_+7] = p4 + nodeOffsetTop;
        
        if(makeQuadratic_)
        {
          connect_[elem*maxNumElemNodes_+8] = p5 + nodeOffsetBottom;
          connect_[elem*maxNumElemNodes_+9] = p6 + nodeOffsetBottom;
          connect_[elem*maxNumElemNodes_+10] = p7 + nodeOffsetBottom;
          connect_[elem*maxNumElemNodes_+11] = p8 + nodeOffsetBottom;

          connect_[elem*maxNumElemNodes_+12] = p5 + nodeOffsetTop;
          connect_[elem*maxNumElemNodes_+13] = p6 + nodeOffsetTop;
          connect_[elem*maxNumElemNodes_+14] = p7 + nodeOffsetTop;
          connect_[elem*maxNumElemNodes_+15] = p8 + nodeOffsetTop;

          connect_[elem*maxNumElemNodes_+16] = p1 + nodeOffsetMiddle;
          connect_[elem*maxNumElemNodes_+17] = p2 + nodeOffsetMiddle;
          connect_[elem*maxNumElemNodes_+18] = p3 + nodeOffsetMiddle;
          connect_[elem*maxNumElemNodes_+19] = p4 + nodeOffsetMiddle;

          if(makeSuperQuadratic_)
          {
            connect_[elem*maxNumElemNodes_+20] = p5 + nodeOffsetMiddle;
            connect_[elem*maxNumElemNodes_+21] = p6 + nodeOffsetMiddle;
            connect_[elem*maxNumElemNodes_+22] = p7 + nodeOffsetMiddle;
            connect_[elem*maxNumElemNodes_+23] = p8 + nodeOffsetMiddle;

            connect_[elem*maxNumElemNodes_+24] = p9 + nodeOffsetBottom;
            connect_[elem*maxNumElemNodes_+25] = p9 + nodeOffsetTop;
            connect_[elem*maxNumElemNodes_+26] = p9 + nodeOffsetMiddle;
            elemTypes_.push_back(ET_HEXA27);
          }
          else
          {
            elemTypes_.push_back(ET_HEXA20);
          }
        }
        else 
        {
          elemTypes_.push_back(ET_HEXA8);
        }
        
        regions_.push_back(region + numRegions);
        elem++;
      }
    }
  }

  void SphericalShellGenerator::CreatePyras()
  {
    UInt numRegions;
    UInt elemOffset;
    UInt elem;
    UInt nodeOffsetBottom;
    UInt nodeOffsetTop;
    UInt nodeOffsetTip;
    UInt nodeOffsetLower;
    UInt nodeOffsetUpper;
    UInt p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13;
    UInt region;
    UInt numElems;
    
    numRegions = regionBaseNames_.size();
    
    for(UInt i=0; i < numRegions; i++)
    {
      std::stringstream sstr;
      
      sstr << "volelems_" << regionBaseNames_[i];
      regionNames_.push_back(sstr.str());
    }

    numRegions *= 2;

    numElems = numInnerElems_ * 2; // inner and outer elems
    numElems += numInnerElems_ * (6 * numLayers_); // 6 pyras inner element and layer

    connect_.resize(numElems * maxNumElemNodes_);
    elemOffset = numInnerElems_ * 2;

    elem = elemOffset;
    for(UInt i=0; i<numInnerElems_; i++) 
    {
      p1 = quadConnectFull_[i*13+0];
      p2 = quadConnectFull_[i*13+1];
      p3 = quadConnectFull_[i*13+2];
      p4 = quadConnectFull_[i*13+3];
      p5 = quadConnectFull_[i*13+8];

      if(makeQuadratic_)
      {
        p6 = quadConnectFull_[i*13+4];
        p7 = quadConnectFull_[i*13+5];
        p8 = quadConnectFull_[i*13+6];
        p9 = quadConnectFull_[i*13+7];
        
        p10 = quadConnectFull_[i*13+9];
        p11 = quadConnectFull_[i*13+10];
        p12 = quadConnectFull_[i*13+11];
        p13 = quadConnectFull_[i*13+12];
      }
      
      region = regions_[i];

      if(p1 == 0 ||
         p2 == 0 ||
         p3 == 0 ||
         p4 == 0 ||
         p5 == 0)
        EXCEPTION("Some point is zero! " << i);
      

      for(UInt layer=0; layer<numLayers_; layer++) 
      {
        nodeOffsetBottom = numInnerPoints_ * 4 * layer;
        nodeOffsetTop = numInnerPoints_ * 4 * (layer+1);
        nodeOffsetTip = numInnerPoints_ * 4 * layer + 2 * numInnerPoints_;
        if(makeQuadratic_)
        {
          nodeOffsetLower = numInnerPoints_ * 4 * layer + 1 * numInnerPoints_;
          nodeOffsetUpper = numInnerPoints_ * 4 * layer + 3 * numInnerPoints_;
        }

        connect_[elem*maxNumElemNodes_+0] = p1 + nodeOffsetBottom;
        connect_[elem*maxNumElemNodes_+1] = p2 + nodeOffsetBottom;
        connect_[elem*maxNumElemNodes_+2] = p3 + nodeOffsetBottom;
        connect_[elem*maxNumElemNodes_+3] = p4 + nodeOffsetBottom;
        connect_[elem*maxNumElemNodes_+4] = p5 + nodeOffsetTip;

        if(makeQuadratic_)
        {
          connect_[elem*maxNumElemNodes_+5] = p6 + nodeOffsetBottom;
          connect_[elem*maxNumElemNodes_+6] = p7 + nodeOffsetBottom;
          connect_[elem*maxNumElemNodes_+7] = p8 + nodeOffsetBottom;
          connect_[elem*maxNumElemNodes_+8] = p9 + nodeOffsetBottom;

          connect_[elem*maxNumElemNodes_+9] = p10 + nodeOffsetLower;
          connect_[elem*maxNumElemNodes_+10] = p11 + nodeOffsetLower;
          connect_[elem*maxNumElemNodes_+11] = p12 + nodeOffsetLower;
          connect_[elem*maxNumElemNodes_+12] = p13 + nodeOffsetLower;

          elemTypes_.push_back(ET_PYRA13);
        }
        else
        {
          elemTypes_.push_back(ET_PYRA5);
        }
        
        regions_.push_back(region + numRegions);
        elem++;
        
        // Reihenfolge geändert!
        connect_[elem*maxNumElemNodes_+0] = p2 + nodeOffsetBottom;
        connect_[elem*maxNumElemNodes_+1] = p1 + nodeOffsetBottom;
        connect_[elem*maxNumElemNodes_+2] = p1 + nodeOffsetTop;
        connect_[elem*maxNumElemNodes_+3] = p2 + nodeOffsetTop;
        connect_[elem*maxNumElemNodes_+4] = p5 + nodeOffsetTip;

        if(makeQuadratic_)
        {
          connect_[elem*maxNumElemNodes_+5] = p6 + nodeOffsetBottom;
          connect_[elem*maxNumElemNodes_+6] = p2 + nodeOffsetTip;
          connect_[elem*maxNumElemNodes_+7] = p6 + nodeOffsetTop;
          connect_[elem*maxNumElemNodes_+8] = p1 + nodeOffsetTip;

          connect_[elem*maxNumElemNodes_+9] = p10 + nodeOffsetLower;
          connect_[elem*maxNumElemNodes_+10] = p11 + nodeOffsetLower;
          connect_[elem*maxNumElemNodes_+11] = p11 + nodeOffsetUpper;
          connect_[elem*maxNumElemNodes_+12] = p10 + nodeOffsetUpper;

          elemTypes_.push_back(ET_PYRA13);
        }
        else
        {
          elemTypes_.push_back(ET_PYRA5);
        }

        regions_.push_back(region + numRegions);
        elem++;

        // Reihenfolge geändert!
        connect_[elem*maxNumElemNodes_+0] = p3 + nodeOffsetBottom;
        connect_[elem*maxNumElemNodes_+1] = p2 + nodeOffsetBottom;
        connect_[elem*maxNumElemNodes_+2] = p2 + nodeOffsetTop;
        connect_[elem*maxNumElemNodes_+3] = p3 + nodeOffsetTop;
        connect_[elem*maxNumElemNodes_+4] = p5 + nodeOffsetTip;

        if(makeQuadratic_)
        {
          connect_[elem*maxNumElemNodes_+5] = p7 + nodeOffsetBottom;
          connect_[elem*maxNumElemNodes_+6] = p3 + nodeOffsetTip;
          connect_[elem*maxNumElemNodes_+7] = p7 + nodeOffsetTop;
          connect_[elem*maxNumElemNodes_+8] = p2 + nodeOffsetTip;

          connect_[elem*maxNumElemNodes_+9] = p11 + nodeOffsetLower;
          connect_[elem*maxNumElemNodes_+10] = p12 + nodeOffsetLower;
          connect_[elem*maxNumElemNodes_+11] = p12 + nodeOffsetUpper;
          connect_[elem*maxNumElemNodes_+12] = p11 + nodeOffsetUpper;

          elemTypes_.push_back(ET_PYRA13);
        }
        else
        {
          elemTypes_.push_back(ET_PYRA5);
        }

        regions_.push_back(region + numRegions);
        elem++;
        
        // Reihenfolge geändert!
        connect_[elem*maxNumElemNodes_+0] = p4 + nodeOffsetBottom;
        connect_[elem*maxNumElemNodes_+1] = p3 + nodeOffsetBottom;
        connect_[elem*maxNumElemNodes_+2] = p3 + nodeOffsetTop;
        connect_[elem*maxNumElemNodes_+3] = p4 + nodeOffsetTop;
        connect_[elem*maxNumElemNodes_+4] = p5 + nodeOffsetTip;

        if(makeQuadratic_)
        {
          connect_[elem*maxNumElemNodes_+5] = p8 + nodeOffsetBottom;
          connect_[elem*maxNumElemNodes_+6] = p4 + nodeOffsetTip;
          connect_[elem*maxNumElemNodes_+7] = p8 + nodeOffsetTop;
          connect_[elem*maxNumElemNodes_+8] = p3 + nodeOffsetTip;

          connect_[elem*maxNumElemNodes_+9] = p12 + nodeOffsetLower;
          connect_[elem*maxNumElemNodes_+10] = p13 + nodeOffsetLower;
          connect_[elem*maxNumElemNodes_+11] = p13 + nodeOffsetUpper;
          connect_[elem*maxNumElemNodes_+12] = p12 + nodeOffsetUpper;

          elemTypes_.push_back(ET_PYRA13);
        }
        else
        {
          elemTypes_.push_back(ET_PYRA5);
        }

        regions_.push_back(region + numRegions);
        elem++;

        // Reihenfolge geändert!
        connect_[elem*maxNumElemNodes_+0] = p1 + nodeOffsetBottom;
        connect_[elem*maxNumElemNodes_+1] = p4 + nodeOffsetBottom;
        connect_[elem*maxNumElemNodes_+2] = p4 + nodeOffsetTop;
        connect_[elem*maxNumElemNodes_+3] = p1 + nodeOffsetTop;
        connect_[elem*maxNumElemNodes_+4] = p5 + nodeOffsetTip;

        if(makeQuadratic_)
        {
          connect_[elem*maxNumElemNodes_+5] = p9 + nodeOffsetBottom;
          connect_[elem*maxNumElemNodes_+6] = p1 + nodeOffsetTip;
          connect_[elem*maxNumElemNodes_+7] = p9 + nodeOffsetTop;
          connect_[elem*maxNumElemNodes_+8] = p4 + nodeOffsetTip;

          connect_[elem*maxNumElemNodes_+9] = p13 + nodeOffsetLower;
          connect_[elem*maxNumElemNodes_+10] = p10 + nodeOffsetLower;
          connect_[elem*maxNumElemNodes_+11] = p10 + nodeOffsetUpper;
          connect_[elem*maxNumElemNodes_+12] = p13 + nodeOffsetUpper;

          elemTypes_.push_back(ET_PYRA13);
        }
        else
        {
          elemTypes_.push_back(ET_PYRA5);
        }

        regions_.push_back(region + numRegions);
        elem++;

        // Reihenfolge geändert!
        connect_[elem*maxNumElemNodes_+0] = p4 + nodeOffsetTop;
        connect_[elem*maxNumElemNodes_+1] = p3 + nodeOffsetTop;
        connect_[elem*maxNumElemNodes_+2] = p2 + nodeOffsetTop;
        connect_[elem*maxNumElemNodes_+3] = p1 + nodeOffsetTop;
        connect_[elem*maxNumElemNodes_+4] = p5 + nodeOffsetTip;

        if(makeQuadratic_)
        {
          connect_[elem*maxNumElemNodes_+5] = p6 + nodeOffsetTop;
          connect_[elem*maxNumElemNodes_+6] = p7 + nodeOffsetTop;
          connect_[elem*maxNumElemNodes_+7] = p8 + nodeOffsetTop;
          connect_[elem*maxNumElemNodes_+8] = p9 + nodeOffsetTop;

          connect_[elem*maxNumElemNodes_+9] = p10 + nodeOffsetUpper;
          connect_[elem*maxNumElemNodes_+10] = p11 + nodeOffsetUpper;
          connect_[elem*maxNumElemNodes_+11] = p12 + nodeOffsetUpper;
          connect_[elem*maxNumElemNodes_+12] = p13 + nodeOffsetUpper;

          elemTypes_.push_back(ET_PYRA13);
        }
        else
        {
          elemTypes_.push_back(ET_PYRA5);
        }

        regions_.push_back(region + numRegions);
        elem++;

      }
    }
  }
  
  void SphericalShellGenerator::MakeQuadraticQuads() 
  {
    UInt connectOffset=0;
    UInt quadConnectOffset=0;
    bool outer;
    UInt nodeOffsetOuter = numInnerPoints_ * 4 * numLayers_;
    
    for(UInt i=0, n=numInnerElems_*2; i<n; i++) 
    {
      quadConnectOffset = (i % numInnerElems_) * 13;
      outer = i / numInnerElems_;

      connect_[connectOffset + 4] = quadConnectFull_[quadConnectOffset + 4];
      connect_[connectOffset + 5] = quadConnectFull_[quadConnectOffset + 5];
      connect_[connectOffset + 6] = quadConnectFull_[quadConnectOffset + 6];
      connect_[connectOffset + 7] = quadConnectFull_[quadConnectOffset + 7];

      if(outer)
      {
        connect_[connectOffset + 4] += nodeOffsetOuter;
        connect_[connectOffset + 5] += nodeOffsetOuter;
        connect_[connectOffset + 6] += nodeOffsetOuter;
        connect_[connectOffset + 7] += nodeOffsetOuter;
      }
      
      if(makeSuperQuadratic_) 
      {
        connect_[connectOffset + 8] = quadConnectFull_[quadConnectOffset + 8];

        if(outer)
          connect_[connectOffset + 8] += nodeOffsetOuter;

        elemTypes_[i] = ET_QUAD9;
      }
      else 
      {
        elemTypes_[i] = ET_QUAD8;
      }
      
      connectOffset += maxNumElemNodes_;
    }
  }

  void SphericalShellGenerator::CreateBaseMesh()
  {
    coords_.clear();
    connect_.clear();
    hashMap_.clear();
    regions_.clear();
    elemTypes_.clear();
    
    // Bottom
    // Node 1
    coords_.push_back(-1);
    coords_.push_back(-1);
    coords_.push_back(-1);

    // Node 2
    coords_.push_back(1);
    coords_.push_back(-1);
    coords_.push_back(-1);

    // Node 3
    coords_.push_back(1);
    coords_.push_back(1);
    coords_.push_back(-1);
  
    // Node 4
    coords_.push_back(-1);
    coords_.push_back(1);
    coords_.push_back(-1);

    // Top
    // Node 5
    coords_.push_back(-1);
    coords_.push_back(-1);
    coords_.push_back(1);

    // Node 6
    coords_.push_back(1);
    coords_.push_back(-1);
    coords_.push_back(1);

    // Node 7
    coords_.push_back(1);
    coords_.push_back(1);
    coords_.push_back(1);
  
    // Node 8
    coords_.push_back(-1);
    coords_.push_back(1);
    coords_.push_back(1);

    // Bottom middle
    // Node 9
    coords_.push_back(0);
    coords_.push_back(-1);
    coords_.push_back(-1);

    // Node 10
    coords_.push_back(1);
    coords_.push_back(0);
    coords_.push_back(-1);

    // Node 11
    coords_.push_back(0);
    coords_.push_back(1);
    coords_.push_back(-1);
  
    // Node 12
    coords_.push_back(-1);
    coords_.push_back(0);
    coords_.push_back(-1);

    // Top middle
    // Node 13
    coords_.push_back(0);
    coords_.push_back(-1);
    coords_.push_back(1);

    // Node 14
    coords_.push_back(1);
    coords_.push_back(0);
    coords_.push_back(1);

    // Node 15
    coords_.push_back(0);
    coords_.push_back(1);
    coords_.push_back(1);
  
    // Node 16
    coords_.push_back(-1);
    coords_.push_back(0);
    coords_.push_back(1);

    // Middle
    // Node 17
    coords_.push_back(-1);
    coords_.push_back(-1);
    coords_.push_back(0);

    // Node 18
    coords_.push_back(1);
    coords_.push_back(-1);
    coords_.push_back(0);

    // Node 19
    coords_.push_back(1);
    coords_.push_back(1);
    coords_.push_back(0);
  
    // Node 20
    coords_.push_back(-1);
    coords_.push_back(1);
    coords_.push_back(0);


    // Middle middle
    // Node 21
    coords_.push_back(0);
    coords_.push_back(-1);
    coords_.push_back(0);

    // Node 22
    coords_.push_back(1);
    coords_.push_back(0);
    coords_.push_back(0);

    // Node 23
    coords_.push_back(0);
    coords_.push_back(1);
    coords_.push_back(0);
  
    // Node 24
    coords_.push_back(-1);
    coords_.push_back(0);
    coords_.push_back(0);

    // Node 25
    coords_.push_back(0);
    coords_.push_back(0);
    coords_.push_back(-1);
  
    // Node 26
    coords_.push_back(0);
    coords_.push_back(0);
    coords_.push_back(1);

    for(int i=0, n=coords_.size()/3; i<n; i++)
    {
      double x = coords_[i*3+0];
      double y = coords_[i*3+1];
      double z = coords_[i*3+2];
    
      unsigned int hashVal = CalcHashVal(x, y, z);
    
      hashMap_[hashVal] = i+1;
    }

    connect_.resize(24 * maxNumElemNodes_);
    UInt elem = 0;
    
    // Bottom
    connect_[elem * maxNumElemNodes_ + 0] = 12;
    connect_[elem * maxNumElemNodes_ + 1] = 25;
    connect_[elem * maxNumElemNodes_ + 2] = 9;
    connect_[elem * maxNumElemNodes_ + 3] = 1;

    elemTypes_.push_back(ET_QUAD4);
    elem++;

    connect_[elem * maxNumElemNodes_ + 0] = 25;
    connect_[elem * maxNumElemNodes_ + 1] = 10;
    connect_[elem * maxNumElemNodes_ + 2] = 2;
    connect_[elem * maxNumElemNodes_ + 3] = 9;

    elemTypes_.push_back(ET_QUAD4);
    elem++;

    connect_[elem * maxNumElemNodes_ + 0] = 11;
    connect_[elem * maxNumElemNodes_ + 1] = 3;
    connect_[elem * maxNumElemNodes_ + 2] = 10;
    connect_[elem * maxNumElemNodes_ + 3] = 25;

    elemTypes_.push_back(ET_QUAD4);
    elem++;

    connect_[elem * maxNumElemNodes_ + 0] = 4;
    connect_[elem * maxNumElemNodes_ + 1] = 11;
    connect_[elem * maxNumElemNodes_ + 2] = 25;
    connect_[elem * maxNumElemNodes_ + 3] = 12;

    elemTypes_.push_back(ET_QUAD4);
    elem++;


    // Front
    connect_[elem * maxNumElemNodes_ + 0] = 1;
    connect_[elem * maxNumElemNodes_ + 1] = 9;
    connect_[elem * maxNumElemNodes_ + 2] = 21;
    connect_[elem * maxNumElemNodes_ + 3] = 17;

    elemTypes_.push_back(ET_QUAD4);
    elem++;

    connect_[elem * maxNumElemNodes_ + 0] = 9;
    connect_[elem * maxNumElemNodes_ + 1] = 2;
    connect_[elem * maxNumElemNodes_ + 2] = 18;
    connect_[elem * maxNumElemNodes_ + 3] = 21;

    elemTypes_.push_back(ET_QUAD4);
    elem++;

    connect_[elem * maxNumElemNodes_ + 0] = 21;
    connect_[elem * maxNumElemNodes_ + 1] = 18;
    connect_[elem * maxNumElemNodes_ + 2] = 6;
    connect_[elem * maxNumElemNodes_ + 3] = 13;

    elemTypes_.push_back(ET_QUAD4);
    elem++;

    connect_[elem * maxNumElemNodes_ + 0] = 17;
    connect_[elem * maxNumElemNodes_ + 1] = 21;
    connect_[elem * maxNumElemNodes_ + 2] = 13;
    connect_[elem * maxNumElemNodes_ + 3] = 5;

    elemTypes_.push_back(ET_QUAD4);
    elem++;


    // Right 
    connect_[elem * maxNumElemNodes_ + 0] = 2;
    connect_[elem * maxNumElemNodes_ + 1] = 10;
    connect_[elem * maxNumElemNodes_ + 2] = 22;
    connect_[elem * maxNumElemNodes_ + 3] = 18;

    elemTypes_.push_back(ET_QUAD4);
    elem++;

    connect_[elem * maxNumElemNodes_ + 0] = 10;
    connect_[elem * maxNumElemNodes_ + 1] = 3;
    connect_[elem * maxNumElemNodes_ + 2] = 19;
    connect_[elem * maxNumElemNodes_ + 3] = 22;

    elemTypes_.push_back(ET_QUAD4);
    elem++;

    connect_[elem * maxNumElemNodes_ + 0] = 22;
    connect_[elem * maxNumElemNodes_ + 1] = 19;
    connect_[elem * maxNumElemNodes_ + 2] = 7;
    connect_[elem * maxNumElemNodes_ + 3] = 14;

    elemTypes_.push_back(ET_QUAD4);
    elem++;

    connect_[elem * maxNumElemNodes_ + 0] = 18;
    connect_[elem * maxNumElemNodes_ + 1] = 22;
    connect_[elem * maxNumElemNodes_ + 2] = 14;
    connect_[elem * maxNumElemNodes_ + 3] = 6;

    elemTypes_.push_back(ET_QUAD4);
    elem++;

    // Back
    connect_[elem * maxNumElemNodes_ + 0] = 3;
    connect_[elem * maxNumElemNodes_ + 1] = 11;
    connect_[elem * maxNumElemNodes_ + 2] = 23;
    connect_[elem * maxNumElemNodes_ + 3] = 19;

    elemTypes_.push_back(ET_QUAD4);
    elem++;

    connect_[elem * maxNumElemNodes_ + 0] = 11;
    connect_[elem * maxNumElemNodes_ + 1] = 4;
    connect_[elem * maxNumElemNodes_ + 2] = 20;
    connect_[elem * maxNumElemNodes_ + 3] = 23;

    elemTypes_.push_back(ET_QUAD4);
    elem++;

    connect_[elem * maxNumElemNodes_ + 0] = 23;
    connect_[elem * maxNumElemNodes_ + 1] = 20;
    connect_[elem * maxNumElemNodes_ + 2] = 8;
    connect_[elem * maxNumElemNodes_ + 3] = 15;

    elemTypes_.push_back(ET_QUAD4);
    elem++;

    connect_[elem * maxNumElemNodes_ + 0] = 19;
    connect_[elem * maxNumElemNodes_ + 1] = 23;
    connect_[elem * maxNumElemNodes_ + 2] = 15;
    connect_[elem * maxNumElemNodes_ + 3] = 7;

    elemTypes_.push_back(ET_QUAD4);
    elem++;


    //Left
    connect_[elem * maxNumElemNodes_ + 0] = 4;
    connect_[elem * maxNumElemNodes_ + 1] = 12;
    connect_[elem * maxNumElemNodes_ + 2] = 24;
    connect_[elem * maxNumElemNodes_ + 3] = 20;

    elemTypes_.push_back(ET_QUAD4);
    elem++;

    connect_[elem * maxNumElemNodes_ + 0] = 12;
    connect_[elem * maxNumElemNodes_ + 1] = 1;
    connect_[elem * maxNumElemNodes_ + 2] = 17;
    connect_[elem * maxNumElemNodes_ + 3] = 24;

    elemTypes_.push_back(ET_QUAD4);
    elem++;

    connect_[elem * maxNumElemNodes_ + 0] = 24;
    connect_[elem * maxNumElemNodes_ + 1] = 17;
    connect_[elem * maxNumElemNodes_ + 2] = 5;
    connect_[elem * maxNumElemNodes_ + 3] = 16;

    elemTypes_.push_back(ET_QUAD4);
    elem++;

    connect_[elem * maxNumElemNodes_ + 0] = 20;
    connect_[elem * maxNumElemNodes_ + 1] = 24;
    connect_[elem * maxNumElemNodes_ + 2] = 16;
    connect_[elem * maxNumElemNodes_ + 3] = 8;

    elemTypes_.push_back(ET_QUAD4);
    elem++;


    // Top
    connect_[elem * maxNumElemNodes_ + 0] = 5;
    connect_[elem * maxNumElemNodes_ + 1] = 13;
    connect_[elem * maxNumElemNodes_ + 2] = 26;
    connect_[elem * maxNumElemNodes_ + 3] = 16;

    elemTypes_.push_back(ET_QUAD4);
    elem++;

    connect_[elem * maxNumElemNodes_ + 0] = 13;
    connect_[elem * maxNumElemNodes_ + 1] = 6;
    connect_[elem * maxNumElemNodes_ + 2] = 14;
    connect_[elem * maxNumElemNodes_ + 3] = 26;

    elemTypes_.push_back(ET_QUAD4);
    elem++;

    connect_[elem * maxNumElemNodes_ + 0] = 26;
    connect_[elem * maxNumElemNodes_ + 1] = 14;
    connect_[elem * maxNumElemNodes_ + 2] = 7;
    connect_[elem * maxNumElemNodes_ + 3] = 15;

    elemTypes_.push_back(ET_QUAD4);
    elem++;

    connect_[elem * maxNumElemNodes_ + 0] = 16;
    connect_[elem * maxNumElemNodes_ + 1] = 26;
    connect_[elem * maxNumElemNodes_ + 2] = 15;
    connect_[elem * maxNumElemNodes_ + 3] = 8;

    elemTypes_.push_back(ET_QUAD4);
    elem++;

    if(octantRegions_) 
    {
      regionBaseNames_.push_back("octant_one");
      regionBaseNames_.push_back("octant_two");
      regionBaseNames_.push_back("octant_three");
      regionBaseNames_.push_back("octant_four");
      regionBaseNames_.push_back("octant_five");
      regionBaseNames_.push_back("octant_six");
      regionBaseNames_.push_back("octant_seven");
      regionBaseNames_.push_back("octant_eight");

      for(UInt i=0, numRegions=regionBaseNames_.size(); i < numRegions; i++)
      {
        std::stringstream sstr;
      
        sstr << "inner_" << regionBaseNames_[i];
        regionNames_.push_back(sstr.str());
      }

      // Bottom
      regions_.push_back(0);
      regions_.push_back(1);
      regions_.push_back(2);
      regions_.push_back(3);

      // Front
      regions_.push_back(0);
      regions_.push_back(1);
      regions_.push_back(5);
      regions_.push_back(4);

      // Right
      regions_.push_back(1);
      regions_.push_back(2);
      regions_.push_back(6);
      regions_.push_back(5);

      // Back
      regions_.push_back(2);
      regions_.push_back(3);
      regions_.push_back(7);
      regions_.push_back(6);

      // Left
      regions_.push_back(3);
      regions_.push_back(0);
      regions_.push_back(4);
      regions_.push_back(7);

      // Top
      regions_.push_back(4);
      regions_.push_back(5);
      regions_.push_back(6);
      regions_.push_back(7);
    }
    else
    {
      regionBaseNames_.push_back("bottom");
      regionBaseNames_.push_back("front");
      regionBaseNames_.push_back("right");
      regionBaseNames_.push_back("back");
      regionBaseNames_.push_back("left");
      regionBaseNames_.push_back("top");

      for(UInt i=0, numRegions=regionBaseNames_.size(); i < numRegions; i++)
      {
        std::stringstream sstr;
      
        sstr << "inner_" << regionBaseNames_[i];
        regionNames_.push_back(sstr.str());

        for(UInt j=0; j<4; j++)
        {
          regions_.push_back(i);
        }
      }
    }

    nodeGroups_["inner_face_bottom_pole"].push_back(25);
    nodeGroups_["inner_face_top_pole"].push_back(26);
    nodeGroups_["inner_face_front_pole"].push_back(21);
    nodeGroups_["inner_face_right_pole"].push_back(22);
    nodeGroups_["inner_face_back_pole"].push_back(23);
    nodeGroups_["inner_face_left_pole"].push_back(24);

    nodeGroups_["inner_edge_nine_pole"].push_back(17);
    nodeGroups_["inner_edge_ten_pole"].push_back(18);
    nodeGroups_["inner_edge_eleven_pole"].push_back(19);
    nodeGroups_["inner_edge_twelve_pole"].push_back(20);

    nodeGroups_["inner_corner_one"].push_back(1);
    nodeGroups_["inner_corner_two"].push_back(2);
    nodeGroups_["inner_corner_three"].push_back(3);
    nodeGroups_["inner_corner_four"].push_back(4);
    nodeGroups_["inner_corner_five"].push_back(5);
    nodeGroups_["inner_corner_six"].push_back(6);
    nodeGroups_["inner_corner_seven"].push_back(7);
    nodeGroups_["inner_corner_eight"].push_back(8);

    numInnerElems_ = regions_.size();
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
        RefineQuad(i, connect2, regions2, elemTypes2, true);
      }

      connect_ = connect2;
      regions_ = regions2;
      elemTypes_ = elemTypes2;
    }

    numInnerElems_ = regions_.size();
    quadConnectFull_.resize(numInnerElems_ * 13);

    for(UInt i=0; i<numInnerElems_; i++) 
    {
      RefineQuad(i, connect2, regions2, elemTypes2, false);
    }

    numInnerPoints_ = coords_.size() / 3;

    std::cout << "Number of quadconnectfull: " << (quadConnectFull_.size() / 13) << std::endl;

  }

  void SphericalShellGenerator::GenSphericalShell(std::vector<double>& coords,
                                                  std::vector<UInt>& connect,
                                                  std::vector<UInt>& elemTypes,
                                                  UInt& maxNumElemNodes,
                                                  std::map<std::string, std::vector<UInt> >& regionElems,
                                                  std::map<std::string, std::vector<UInt> >& nodeGroups)
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

    Double deltaRadius = (outerRadius_ - innerRadius_) / ( 4 * numLayers_ );
    std::vector<Double> myCoords, myCoords2;
    
    myCoords = coords_;
    coords_.clear();
    for(UInt layer=0, n=numLayers_*4; layer <= n; layer++) 
    {
      ProjectCoordsToSphere(myCoords, myCoords2, innerRadius_ + layer * deltaRadius);
      std::copy(myCoords2.begin(), myCoords2.end(), std::back_inserter(coords_));
    }

    CreateOuterQuads();

    CreateHexas();

    // CreatePyras();

    if(makeQuadratic_)
      MakeQuadraticQuads();

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
    nodeGroups = nodeGroups_;
    

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
