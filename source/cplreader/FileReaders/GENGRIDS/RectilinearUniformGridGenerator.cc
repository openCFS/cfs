// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <cmath>
#include <iterator>

#include "Domain/elem.hh"
#include "RectilinearUniformGridGenerator.hh"

namespace CoupledField
{
  RectilinearUniformGridGenerator::RectilinearUniformGridGenerator() :
    maxNumElemNodes_(27)
  {
  }
  
  void RectilinearUniformGridGenerator::GenUniformGrid(std::vector<double>& coords,
                                                       std::vector<UInt>& connect,
                                                       std::vector<UInt>& elemTypes,
                                                       UInt& maxNumElemNodes,
                                                       std::map<std::string, std::vector<UInt> >& regionElems,
                                                       std::map<std::string, std::vector<UInt> >& nodeGroups,
                                                       HexaType hexaType)
  {
    coords.clear();
    connect.clear();
    elemTypes.clear();
    
    std::vector<Double> xCoords;
    std::vector<Double> yCoords;
    std::vector<Double> zCoords;
    std::vector<UInt> xNElems;
    std::vector<UInt> yNElems;
    std::vector<UInt> zNElems;
    
    /*
    xCoords.push_back(0.0);
    xCoords.push_back(10.0);
    //    xCoords.push_back(11.0);
    
    //    yCoords.push_back(-0.1);
    yCoords.push_back(0.0);
    yCoords.push_back(1.0);
    
    zCoords.push_back(0.0);
    zCoords.push_back(0.5);
    //    zCoords.push_back(0.7);
    
    xNElems.push_back(1);
    xNElems.push_back(2);
    
    yNElems.push_back(2);
    yNElems.push_back(2);
    
    zNElems.push_back(3);
    zNElems.push_back(1);
    */

    /*
    // Cube 3D
    xCoords.push_back(0);
    xNElems.push_back(10);
    xCoords.push_back(1);

    yCoords.push_back(0);
    yNElems.push_back(10);
    yCoords.push_back(1);

    zCoords.push_back(0);
    zNElems.push_back(10);
    zCoords.push_back(1);
    */

    Double xCoord = 0;
    xCoords.push_back(xCoord);
    xNElems.push_back(5); // 10
    xCoord += 1.0e-1;
    xCoords.push_back(xCoord);
    xNElems.push_back(2); // 5
    xCoord += 4.5e-2;
    xCoords.push_back(xCoord);
    xNElems.push_back(12); // 25
    xCoord += 2.5e-1;
    xCoords.push_back(xCoord);
    xNElems.push_back(6); // 12
    xCoord += 4.5e-2;
    xCoords.push_back(xCoord);
    xNElems.push_back(1);
    xCoord += 1.0e-2;
    xCoords.push_back(xCoord);

    Double yCoord = 0;
    yCoords.push_back(yCoord);
    yNElems.push_back(7); // 15
    yCoord += 1.525e-1;
    yCoords.push_back(yCoord);
    yNElems.push_back(6); // 12
    yCoord += 4.5e-2;
    yCoords.push_back(yCoord);
    yNElems.push_back(7); // 15
    yCoord += 1.525e-1;
    yCoords.push_back(yCoord);

    zCoords.push_back(0.0);
    zNElems.push_back(2);
    zCoords.push_back(1.0e-3);
    zNElems.push_back(2);
    zCoords.push_back(2.0e-3);

    CreateMesh(xCoords, yCoords, zCoords,
               xNElems, yNElems, zNElems, hexaType);
    
    UInt numElems = connect_.size() / maxNumElemNodes_;
    std::cout << "NumElems: " << numElems << std::endl;
    
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
    
  }

  void RectilinearUniformGridGenerator::CreateMesh(const std::vector<Double>& xCoords,
                                                   const std::vector<Double>& yCoords,
                                                   const std::vector<Double>& zCoords,
                                                   const std::vector<UInt>& xNElems,
                                                   const std::vector<UInt>& yNElems,
                                                   const std::vector<UInt>& zNElems,
                                                   HexaType hexaType)
  {
    coords_.clear();
    connect_.clear();
    regions_.clear();
    elemTypes_.clear();

    std::vector<Double> xSampleCoords;
    std::vector<Double> ySampleCoords;
    std::vector<Double> zSampleCoords;

    std::vector<UInt> xNodeOff;
    std::vector<UInt> yNodeOff;
    std::vector<UInt> zNodeOff;
    
    UInt nX=xCoords.size();
    UInt nY=yCoords.size();
    UInt nZ=zCoords.size();

    UInt nRegionX = nX - 1;
    UInt nRegionY = nY - 1;
    UInt nRegionZ = nZ - 1;

    // Generate x-node coords
    UInt xIdx=0;
    UInt nOffY = 1;
    xNodeOff.push_back(0);
    xSampleCoords.push_back(xCoords[xIdx]);
    for(xIdx=1; xIdx < nX; xIdx++) {
      Double dx = xCoords[xIdx] - xCoords[xIdx-1];      
      dx /= 2.0 * (Double) xNElems[xIdx-1];
      for(UInt xNE=1; xNE <= xNElems[xIdx-1] * 2; xNE++) 
      {
        xSampleCoords.push_back(xCoords[xIdx-1] + xNE * dx);
      }
      xNodeOff.push_back((*xNodeOff.rbegin()) + xNElems[xIdx-1] * 2);
      nOffY += xNElems[xIdx-1] * 2;
    }
    nOffY *=2;

    // Generate y-node coords
    UInt yIdx=0;
    UInt nOffZ = 1;
    yNodeOff.push_back(0);
    ySampleCoords.push_back(yCoords[yIdx]);
    for(yIdx=1; yIdx < nY; yIdx++) {
      Double dy = yCoords[yIdx] - yCoords[yIdx-1];      
      dy /= 2.0 * (Double) yNElems[yIdx-1];
      for(UInt yNE=1; yNE <= yNElems[yIdx-1] * 2; yNE++) 
      {
        ySampleCoords.push_back(yCoords[yIdx-1] + yNE * dy);
      }
      yNodeOff.push_back((*yNodeOff.rbegin()) + yNElems[yIdx-1] * 2);
      nOffZ += yNElems[yIdx-1] * 2;
    }
    nOffZ *= nOffY;

    // Generate z-node coords
    UInt zIdx=0;
    zNodeOff.push_back(0);
    zSampleCoords.push_back(zCoords[zIdx]);
    for(zIdx=1; zIdx < nZ; zIdx++) {
      Double dz = zCoords[zIdx] - zCoords[zIdx-1];      
      dz /= 2.0 * (Double) zNElems[zIdx-1];
      for(UInt zNE=1; zNE <= zNElems[zIdx-1] * 2; zNE++) 
      {
        zSampleCoords.push_back(zCoords[zIdx-1] + zNE * dz);
      }
      zNodeOff.push_back((*zNodeOff.rbegin()) + zNElems[zIdx-1] * 2);
    }
    
    for(UInt zIdx=0; zIdx < zSampleCoords.size(); zIdx++) {
      for(UInt yIdx=0; yIdx < ySampleCoords.size(); yIdx++) {
        for(UInt xIdx=0; xIdx < xSampleCoords.size(); xIdx++) {
          coords_.push_back(xSampleCoords[xIdx]);
          coords_.push_back(ySampleCoords[yIdx]);
          coords_.push_back(zSampleCoords[zIdx]);
        }
      }
    }
    
    // Generate elems
    std::cout << "nOffY " << nOffY << std::endl;
    std::cout << "nOffZ " << nOffZ << std::endl;
    

    UInt elem = 0;
    std::stringstream sstr;
    std::string xyzIdx;
    std::vector<UInt> myConnect(maxNumElemNodes_);
    std::vector<UInt> hex27Con(maxNumElemNodes_);
    UInt regionIdx = 0;


    for(UInt zIdx=0; zIdx < nRegionZ; zIdx++) 
    {
      std::cout << " zNodeOff[zIdx] " <<  zNodeOff[zIdx] << std::endl;
      UInt zBase = zNodeOff[zIdx] * (nOffZ / 2);
      std::cout << " zBase " <<  zBase << std::endl;
      for(UInt yIdx=0; yIdx < nRegionY; yIdx++) 
      {
        std::cout << " yNodeOff[yIdx] " <<  yNodeOff[yIdx] << std::endl;
        UInt yBase = yNodeOff[yIdx] * (nOffY / 2);
        std::cout << " yBase " <<  yBase << std::endl;
        for(UInt xIdx=0; xIdx < nRegionX; xIdx++) 
        {
          UInt xBase = zBase + yBase + xNodeOff[xIdx] + 1;
          
          sstr.str(""); sstr.clear();
          sstr << (xIdx+1) << "_" << (yIdx+1) << "_" << (zIdx+1);
          xyzIdx = sstr.str();

          sstr.str(""); sstr.clear(); sstr << "volume_" << xyzIdx;
          regionNames_.push_back(sstr.str());

          sstr.str(""); sstr.clear(); sstr << "surf_bottom_" << xyzIdx;
          regionNames_.push_back(sstr.str());

          sstr.str(""); sstr.clear(); sstr << "surf_top_" << xyzIdx;
          regionNames_.push_back(sstr.str());

          sstr.str(""); sstr.clear(); sstr << "surf_left_" << xyzIdx;
          regionNames_.push_back(sstr.str());

          sstr.str(""); sstr.clear(); sstr << "surf_right_" << xyzIdx;
          regionNames_.push_back(sstr.str());

          sstr.str(""); sstr.clear(); sstr << "surf_front_" << xyzIdx;
          regionNames_.push_back(sstr.str());

          sstr.str(""); sstr.clear(); sstr << "surf_back_" << xyzIdx;
          regionNames_.push_back(sstr.str());

          for(UInt xNE=0; xNE < xNElems[xIdx]; xNE++) 
          {
            for(UInt yNE=0; yNE < yNElems[yIdx]; yNE++) 
            {
              for(UInt zNE=0; zNE < zNElems[zIdx]; zNE++) 
              {

                std::fill(myConnect.begin(), myConnect.end(), 0);
                // bottom
                myConnect[0] = xBase + zNE*nOffZ + yNE*nOffY + xNE * 2;
                myConnect[1] = myConnect[0] + 2;
                myConnect[2] = xBase + zNE*nOffZ + yNE*nOffY + nOffY + (xNE + 1) * 2;
                myConnect[3] = myConnect[2] - 2;
                //top
                myConnect[4] = myConnect[0] + nOffZ;
                myConnect[5] = myConnect[1] + nOffZ;
                myConnect[6] = myConnect[2] + nOffZ;
                myConnect[7] = myConnect[3] + nOffZ;
                // bottom mid
                myConnect[8] = myConnect[0] + 1;
                myConnect[11] = myConnect[0] + nOffY / 2;
                myConnect[9] = myConnect[11] + 2;
                myConnect[10] = myConnect[3] + 1;
                // top mid
                myConnect[12] = myConnect[4] + 1;
                myConnect[15] = myConnect[4] + nOffY / 2;
                myConnect[13] = myConnect[15] + 2;
                myConnect[14] = myConnect[7] + 1;
                // middle
                myConnect[16] = myConnect[0] + nOffZ / 2;
                myConnect[17] = myConnect[1] + nOffZ / 2;
                myConnect[18] = myConnect[2] + nOffZ / 2;
                myConnect[19] = myConnect[3] + nOffZ / 2;
                // middle mid
                myConnect[20] = myConnect[8] + nOffZ / 2;
                myConnect[21] = myConnect[9] + nOffZ / 2;
                myConnect[22] = myConnect[10] + nOffZ / 2;
                myConnect[23] = myConnect[11] + nOffZ / 2;
                // vertical
                myConnect[24] = myConnect[11] + 1;
                myConnect[25] = myConnect[15] + 1;
                myConnect[26] = myConnect[23] + 1;

                std::copy(myConnect.begin(), myConnect.end(), hex27Con.begin());
                std::copy(myConnect.begin(), myConnect.end(), std::back_inserter(connect_));
                regions_.push_back(regionIdx+0);
                switch(hexaType) 
                {
                case LINEAR:
                  elemTypes_.push_back(Elem::HEXA8);
                  break;
                case SERENDIPITY:
                  elemTypes_.push_back(Elem::HEXA20);
                  break;
                case LAGRANGE:
                  elemTypes_.push_back(Elem::HEXA27);
                  break;
                }
                elem++;


                if(zNE == 0) 
                {  
                  // Bottom
                  std::fill(myConnect.begin(), myConnect.end(), 0);
                  myConnect[0] = hex27Con[0];
                  myConnect[1] = hex27Con[1];
                  myConnect[2] = hex27Con[2];
                  myConnect[3] = hex27Con[3];

                  myConnect[4] = hex27Con[8];
                  myConnect[5] = hex27Con[9];
                  myConnect[6] = hex27Con[10];
                  myConnect[7] = hex27Con[11];

                  myConnect[8] = hex27Con[24];

                  std::copy(myConnect.begin(), myConnect.end(), std::back_inserter(connect_));
                  regions_.push_back(regionIdx+1);
                  switch(hexaType) 
                  {
                  case LINEAR:
                    elemTypes_.push_back(Elem::QUAD4);
                    break;
                  case SERENDIPITY:
                    elemTypes_.push_back(Elem::QUAD8);
                    break;
                  case LAGRANGE:
                    elemTypes_.push_back(Elem::QUAD9);
                    break;
                  }
                  elem++;
                }
                
                if(zNE == zNElems[zIdx]-1) 
                {  
                  // Top
                  std::fill(myConnect.begin(), myConnect.end(), 0);
                  myConnect[0] = hex27Con[4];
                  myConnect[1] = hex27Con[5];
                  myConnect[2] = hex27Con[6];
                  myConnect[3] = hex27Con[7];

                  myConnect[4] = hex27Con[12];
                  myConnect[5] = hex27Con[13];
                  myConnect[6] = hex27Con[14];
                  myConnect[7] = hex27Con[15];

                  myConnect[8] = hex27Con[25];

                  std::copy(myConnect.begin(), myConnect.end(), std::back_inserter(connect_));
                  regions_.push_back(regionIdx+2);
                  switch(hexaType) 
                  {
                  case LINEAR:
                    elemTypes_.push_back(Elem::QUAD4);
                    break;
                  case SERENDIPITY:
                    elemTypes_.push_back(Elem::QUAD8);
                    break;
                  case LAGRANGE:
                    elemTypes_.push_back(Elem::QUAD9);
                    break;
                  }
                  elem++;
                }
                
                if(xNE == 0) 
                {  
                  // Left
                  std::fill(myConnect.begin(), myConnect.end(), 0);
                  myConnect[0] = hex27Con[0];
                  myConnect[1] = hex27Con[3];
                  myConnect[2] = hex27Con[7];
                  myConnect[3] = hex27Con[4];

                  myConnect[4] = hex27Con[11];
                  myConnect[5] = hex27Con[19];
                  myConnect[6] = hex27Con[15];
                  myConnect[7] = hex27Con[16];

                  myConnect[8] = hex27Con[23];

                  std::copy(myConnect.begin(), myConnect.end(), std::back_inserter(connect_));
                  regions_.push_back(regionIdx+3);
                  switch(hexaType) 
                  {
                  case LINEAR:
                    elemTypes_.push_back(Elem::QUAD4);
                    break;
                  case SERENDIPITY:
                    elemTypes_.push_back(Elem::QUAD8);
                    break;
                  case LAGRANGE:
                    elemTypes_.push_back(Elem::QUAD9);
                    break;
                  }
                  elem++;
                }

                if(xNE == xNElems[xIdx]-1) 
                {  
                  // Right
                  std::fill(myConnect.begin(), myConnect.end(), 0);
                  myConnect[0] = hex27Con[1];
                  myConnect[1] = hex27Con[2];
                  myConnect[2] = hex27Con[6];
                  myConnect[3] = hex27Con[5];

                  myConnect[4] = hex27Con[9];
                  myConnect[5] = hex27Con[18];
                  myConnect[6] = hex27Con[13];
                  myConnect[7] = hex27Con[17];

                  myConnect[8] = hex27Con[21];

                  std::copy(myConnect.begin(), myConnect.end(), std::back_inserter(connect_));
                  regions_.push_back(regionIdx+4);
                  switch(hexaType) 
                  {
                  case LINEAR:
                    elemTypes_.push_back(Elem::QUAD4);
                    break;
                  case SERENDIPITY:
                    elemTypes_.push_back(Elem::QUAD8);
                    break;
                  case LAGRANGE:
                    elemTypes_.push_back(Elem::QUAD9);
                    break;
                  }
                  elem++;
                }

                if(yNE == 0) 
                {  
                  // Front
                  std::fill(myConnect.begin(), myConnect.end(), 0);
                  myConnect[0] = hex27Con[0];
                  myConnect[1] = hex27Con[1];
                  myConnect[2] = hex27Con[5];
                  myConnect[3] = hex27Con[4];

                  myConnect[4] = hex27Con[8];
                  myConnect[5] = hex27Con[17];
                  myConnect[6] = hex27Con[12];
                  myConnect[7] = hex27Con[16];

                  myConnect[8] = hex27Con[20];

                  std::copy(myConnect.begin(), myConnect.end(), std::back_inserter(connect_));
                  regions_.push_back(regionIdx+5);
                  switch(hexaType) 
                  {
                  case LINEAR:
                    elemTypes_.push_back(Elem::QUAD4);
                    break;
                  case SERENDIPITY:
                    elemTypes_.push_back(Elem::QUAD8);
                    break;
                  case LAGRANGE:
                    elemTypes_.push_back(Elem::QUAD9);
                    break;
                  }
                  elem++;
                }

                if(yNE == yNElems[yIdx]-1) 
                {  
                  // Back
                  std::fill(myConnect.begin(), myConnect.end(), 0);                
                  myConnect[0] = hex27Con[3];
                  myConnect[1] = hex27Con[2];
                  myConnect[2] = hex27Con[6];
                  myConnect[3] = hex27Con[7];

                  myConnect[4] = hex27Con[10];
                  myConnect[5] = hex27Con[18];
                  myConnect[6] = hex27Con[14];
                  myConnect[7] = hex27Con[19];

                  myConnect[8] = hex27Con[22];

                  std::copy(myConnect.begin(), myConnect.end(), std::back_inserter(connect_));
                  regions_.push_back(regionIdx+6);
                  switch(hexaType) 
                  {
                  case LINEAR:
                    elemTypes_.push_back(Elem::QUAD4);
                    break;
                  case SERENDIPITY:
                    elemTypes_.push_back(Elem::QUAD8);
                    break;
                  case LAGRANGE:
                    elemTypes_.push_back(Elem::QUAD9);
                    break;
                  }
                  elem++;
                }
              }
            }
          }

          regionIdx += 7;
          
        }
      }
    }
  }
  
} // end of namespace
