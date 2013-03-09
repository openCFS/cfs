/*
 * FieldViewFile.cpp
 *
 *  Created on: 27.09.2011
 *      Author: willmitzer
 */

#include "FieldViewFile.h"
#include "fv_reader_tags.h"

#include <iostream>
#include <fstream>
#include <cstring>

#define TRACE 0
#define DEBUG 0
#define VERBOSE 1
#define NODEVERBOSE 0
#define ELEMENTVERBOSE 0
#define VARVERBOSE 0

namespace CoupledField {

FieldViewFile::FieldViewFile() {

  nodeCoords_ = NULL;
  nodeVariables_ = NULL;
  numResultFacesStd_.resize(0);
  numResultFacesPoly_.resize(0);

}

void FieldViewFile::open(std::string fileName) {

  // error checking
  if( fileName.empty() ) {
    std::cout << "ERROR: empty file name" << std::endl;
    return;
  }

  // create the ifstream
  fileStream_.open(fileName.c_str(), std::ios::in | std::ios::binary);

  // check if its good for IO
  if( !fileStream_.good() ) {
    std::cout << "ERROR: opening file '" << fileName.c_str() << "' failed" << std::endl;
    return;
  }

  // save the file name
  fileName_ = fileName.c_str();

  // read the header and find out the file type and number of grids
  readHeader();

  // depending on whether we are interested in mesh or in value information, skip other stuff

  // grids
  readGrids(false, false);

  // EOF reached, close ifstream
  if( fileStream_.is_open() ) {
    fileStream_.close();
  }

}

void FieldViewFile::openGrid(std::string fileName) {

  // error checking
  if( fileName.empty() ) {
    std::cout << "ERROR: empty file name" << std::endl;
    return;
  }

  // create the ifstream
  fileStream_.open(fileName.c_str(), std::ios::in | std::ios::binary);

  // check if its good for IO
  if( !fileStream_.good() ) {
    std::cout << "ERROR: opening file '" << fileName.c_str() << "' failed" << std::endl;
    return;
  }

  // save the file name
  fileName_ = fileName.c_str();

  // read the header of the file with basic information
  readHeader();

  // depending on whether we are interested in mesh or in value information, skip other stuff

  // grids
  readGrids(false, true);

  // EOF reached, close ifstream
  if( fileStream_.is_open() ) {
    fileStream_.close();
  }

}

void FieldViewFile::openVars(std::string fileName) {

  // error checking
  if( fileName.empty() ) {
    std::cout << "ERROR: empty file name" << std::endl;
    return;
  }

  // create the ifstream
  fileStream_.open(fileName.c_str(), std::ios::in | std::ios::binary);

  // check if its good for IO
  if( !fileStream_.good() ) {
    std::cout << "ERROR: opening file '" << fileName.c_str() << "' failed" << std::endl;
    return;
  }

  // save the file name
  fileName_ = fileName.c_str();

  // read the header and find out the file type and number of grids
  readHeader();

  // depending on whether we are interested in mesh or in value information, skip other stuff

  // grids
  readGrids(true, false);

  // EOF reached, close ifstream
  if( fileStream_.is_open() ) {
    fileStream_.close();
  }

}

FieldViewFile::~FieldViewFile() {

  // close the ifstream
  if( fileStream_.is_open() ) {
    fileStream_.close();
  }

  // free the node coordinate arrays
  if( nodeCoords_ != NULL ) {

    for( int i = 0; i < numGrids_; i++ ) {
      delete[] nodeCoords_[i];
    }
    delete[] nodeCoords_;

  }

  // free the variables arrays
  if( nodeVariables_ != NULL ) {

    for( int i = 0; i < numGrids_; i++ ) {
      delete[] nodeVariables_[i];
    }
    delete[] nodeVariables_;

  }

}

float* FieldViewFile::getNodalVars(int gridIdx) {
  return nodeVariables_[gridIdx];
}

std::vector<std::string> FieldViewFile::getNodalVarNames() {
  return nodalVarNames_;
}

int FieldViewFile::getNumGrids() {

  return numGrids_;

}

int FieldViewFile::getNumNodes(int gridIdx) {

  return numNodes_[gridIdx];

}

int FieldViewFile::getNumElements(int gridIdx) {

  int numElements = numTetElems_[gridIdx] + numPyrElems_[gridIdx] + numPriElems_[gridIdx] + numHexElems_[gridIdx];

  return numElements;

}

float* FieldViewFile::getNodeCoordinates(int gridIdx) {

  if( nodeCoords_ != NULL ) {
    return nodeCoords_[gridIdx];
  } else {
    return NULL;
  }

}

std::string FieldViewFile::getFileName() {
  return fileName_;
}

int FieldViewFile::getNumTetElems(int gridIdx) {
  return numTetElems_[gridIdx];
}

int FieldViewFile::getNumHexElems(int gridIdx) {
  return numHexElems_[gridIdx];
}

int FieldViewFile::getNumPriElems(int gridIdx) {
  return numPriElems_[gridIdx];
}

int FieldViewFile::getNumPyrElems(int gridIdx) {
  return numPyrElems_[gridIdx];
}

int* FieldViewFile::getTetElems(int gridIdx) {
  return tetElemNodes_[gridIdx];
}

int* FieldViewFile::getHexElems(int gridIdx) {
  return hexElemNodes_[gridIdx];
}

int* FieldViewFile::getPriElems(int gridIdx) {
  return priElemNodes_[gridIdx];
}

int* FieldViewFile::getPyrElems(int gridIdx) {
  return pyrElemNodes_[gridIdx];
}

void FieldViewFile::readBoundaryVariables(int gridIdx) {

  if( TRACE ) std::cout << "TRACE: FieldViewFile::readBoundaryVariables()" << std::endl;

  if ( fileType_ == FV_GRIDS_FILE ) {
    // does not contain variables sections
    return;
  }

  // buffer for ifstream.read()
  char cbuf[80];
  // reinterpret_cast target for ints
  int *ibuf;

  // read section header
  fileStream_.read(cbuf, sizeof(int));
  ibuf = reinterpret_cast<int*>(cbuf);

  if( ibuf[0] != FV_BNDRY_VARS ) {
    std::cout << "ERROR: expected FV_BNDRY_VARS section, found " << ibuf[0] << std::endl;
    return;
  }

  // boundary information is ignored

  // number of floats in this section
  if( numResultFacesStd_.empty() ) {
    std::cout << "ERROR: cannot read FV_BNDRY_VARS section: no mesh information available" << std::endl;
    return;
  }

  int offset = numResultFacesStd_.at(gridIdx) * numBndryVars_ * sizeof(float);
  if( DEBUG ) std::cout << "int offset = " << numResultFacesStd_.at(gridIdx) << " * " << numBndryVars_ << " * " << sizeof(float) << std::endl;

  if( DEBUG ) std::cout << "FV_BNDRY_VARS: size of " << offset << " bytes" << std::endl;

  fileStream_.seekg(offset, std::ios_base::cur);

  // read section header
  fileStream_.read(cbuf, sizeof(int));
  ibuf = reinterpret_cast<int*>(cbuf);

  if( ibuf[0] != FV_ARB_POLY_BNDRY_VARS ) {
    std::cout << "ERROR: expected FV_ARB_POLY_BNDRY_VARS section, found " << ibuf[0] << std::endl;
    isValidFile_ = false;
    return;
  }

  // boundary information is ignored

  // number of floats in this section
  if( numResultFacesPoly_.empty() ) {
    std::cout << "ERROR: cannot read FV_ARB_POLY_BNDRY_VARS section: no mesh information available" << std::endl;
    return;
  }

  offset = numResultFacesPoly_.at(gridIdx) * numBndryVars_ * sizeof(float);
  if( DEBUG ) std::cout << "FV_ARB_POLY_BNDRY_VARS: size of " << offset << " bytes" << std::endl;

  // advance file pointer
  fileStream_.seekg(offset, std::ios_base::cur);

}

void FieldViewFile::readVariables(int gridIdx, bool isIgnored) {

  if( TRACE ) std::cout << "TRACE: FieldViewFile::readVariables()" << std::endl;

  if ( fileType_ == FV_GRIDS_FILE ) {
    // does not contain variables sections
    return;
  }

  // buffer for ifstream.read()
  char cbuf[80];
  // reinterpret_cast target for ints
  int *ibuf;
  // reinterpret_cast target for floats
  float *fbuf;

  // read section header
  fileStream_.read(cbuf, sizeof(int));
  ibuf = reinterpret_cast<int*>(cbuf);

  if( ibuf[0] != FV_VARIABLES ) {
    std::cout << "ERROR: expected FV_VARIABLES section, found " << ibuf[0] << std::endl;
    return;
  }

  if( DEBUG ) std::cout << "Section variables (FV_VARIABLES)" << std::endl;

  if( isIgnored ) {

    // data should be ignored, just modify the pointer
    if( DEBUG ) std::cout << "Ignoring variables" << std::endl;
    // for each entry, there are numNodes_[gridIdx] * numNodalVars_ floats
    int offset = sizeof(float) * numNodes_[gridIdx] * numNodalVars_;
    fileStream_.seekg(offset, std::ios_base::cur);

  } else {

    // number of float values to read
    int numFloats = numNodes_[gridIdx] * numNodalVars_;

    // allocate memory
    nodeVariables_[gridIdx] = new float[numFloats];

    // read data
    char *vbuf = new char[sizeof(float) * numFloats];
    fileStream_.read(vbuf, sizeof(float) * numFloats);
    fbuf = reinterpret_cast<float*>(vbuf);

    // copy array
    for( int i = 0; i < numFloats; i++ ) {
      nodeVariables_[gridIdx][i] = fbuf[i];
    }

  }

  if( VARVERBOSE ) {
    for( int nodeIdx = 0; nodeIdx < numNodes_[gridIdx]; nodeIdx++ ) {
      std::cout << "Node " << nodeIdx << ": ";
      for( int i = 0; i < numNodalVars_; i++ ) {
        std::cout << nodeVariables_[gridIdx][(i*numNodes_[gridIdx]+nodeIdx)];
        if( i < numNodalVars_-1 ) std::cout << " ";
      }
      std::cout << std::endl;

    }

  }

  if( DEBUG ) std::cout << "Section variables (FV_VARIABLES) finished" << std::endl;

}

void FieldViewFile::readElements(int gridIdx, bool isIgnored) {

  if( TRACE ) std::cout << "TRACE: FieldViewFile::readElements()" << std::endl;

  if ( fileType_ == FV_RESULTS_FILE ) {
    // does not contain elements sections
    std::cout << "readElements() was called for FV_RESULTS_FILE" << std::endl;
    return;
  }

  // buffer for ifstream.read()
  char cbuf[80];
  // reinterpret_cast target for ints
  int *ibuf;

  // whether the elements sections are finished
  bool isFinished = false;

  while( isFinished == false ) {

    fileStream_.read(cbuf, sizeof(int));
    ibuf = reinterpret_cast<int*>(cbuf);

    // check if end of file has been reached
    if( fileStream_.eof() ) {
      if( VERBOSE ) std::cout << "End of file has been reached" << std::endl;
      isFinished = true;
      continue;
    }

    switch( ibuf[0] ) {

    case FV_ELEMENTS: {

      if( DEBUG ) std::cout << "Section standard 3-D elements (FV_ELEMENTS)" << std::endl;

      // number of elements in this section
      fileStream_.read(cbuf, sizeof(int) * 4);
      ibuf = reinterpret_cast<int*>(cbuf);
      int numTet = ibuf[0];
      int numHex = ibuf[1];
      int numPri = ibuf[2];
      int numPyr = ibuf[3];
      int numAll = ibuf[0] + ibuf[1] + ibuf[2] + ibuf[3];

      if( VERBOSE ) {
        std::cout << "Number of tetrahedron elements: " << numTet << std::endl;
        std::cout << "Number of hexahedron elements: " << numHex << std::endl;
        std::cout << "Number of prism elements: " << numPri << std::endl;
        std::cout << "Number of pyramidal elements: " << numPyr << std::endl;
        std::cout << "Number of elements: " << numAll << std::endl;
      }

      // initialize the arrays for the nodes

      // only if there actually are elements of this type
      if( numTet != 0 ) {

        if( numTetElems_[gridIdx] == 0 ) {

          // first element section for this grid
          tetElemNodes_[gridIdx] = new int[4*numTet];

        } else {

          // not the first element section for this grid, resize
          int *temp = new int[4*(numTetElems_[gridIdx] + numTet)];

          for( int i = 0; i < numTetElems_[gridIdx]; i++ ) {
            temp[i] = tetElemNodes_[gridIdx][i];
          }

          // delete old array, set new one
          delete[] tetElemNodes_[gridIdx];
          tetElemNodes_[gridIdx] = temp;

        }

      }

      if( numHex != 0 ) {

        if( numHexElems_[gridIdx] == 0 ) {

          // first element section for this grid
          hexElemNodes_[gridIdx] = new int[8*numHex];

        } else {

          // not the first element section for this grid, resize
          int *temp = new int[8*(numHexElems_[gridIdx] + numHex)];

          for( int i = 0; i < numHexElems_[gridIdx]; i++ ) {
            temp[i] = hexElemNodes_[gridIdx][i];
          }

          // delete old array, set new one
          delete[] hexElemNodes_[gridIdx];
          hexElemNodes_[gridIdx] = temp;

        }

      }

      if( numPri != 0 ) {

        if( numPriElems_[gridIdx] == 0 ) {

          // first element section for this grid
          priElemNodes_[gridIdx] = new int[6*numPri];

        } else {

          // not the first element section for this grid, resize
          int *temp = new int[6*(numPriElems_[gridIdx] + numPri)];

          for( int i = 0; i < numPriElems_[gridIdx]; i++ ) {
            temp[i] = priElemNodes_[gridIdx][i];
          }

          // delete old array, set new one
          delete[] priElemNodes_[gridIdx];
          priElemNodes_[gridIdx] = temp;

        }

      }

      if( numPyr != 0 ) {

        if( numPyrElems_[gridIdx] == 0 ) {

          // first element section for this grid
          pyrElemNodes_[gridIdx] = new int[5*numPyr];

        } else {

          // not the first element section for this grid, resize
          int *temp = new int[5*(numPyrElems_[gridIdx] + numPyr)];

          for( int i = 0; i < numPyrElems_[gridIdx]; i++ ) {
            temp[i] = pyrElemNodes_[gridIdx][i];
          }

          // delete old array, set new one
          delete[] pyrElemNodes_[gridIdx];
          pyrElemNodes_[gridIdx] = temp;

        }

      }

      // walk the elements
      for( int i = 0; i < numAll; i++ ) {

        // determine element type; 4-byte word
        fileStream_.read(cbuf, 4);
        ibuf = reinterpret_cast<int*>(cbuf);
        int elemType = ibuf[0] >> 18;

        // node numbers from .fvuns file are stored i.e. 1...numNodes, NOT 0...numNodes-1

        switch( elemType ) {

        case 1:

          // tetrahedron: 4 nodes

          fileStream_.read(cbuf, 4*sizeof(int));
          ibuf = reinterpret_cast<int*>(cbuf);

          // store the nodes
          for( int j = 0; j < 4; j++ ) {
            tetElemNodes_[gridIdx][(4*numTetElems_[gridIdx])+j] = ibuf[j];
          }

          if( ELEMENTVERBOSE ) {
            std::cout << "Element " << i << ": ";
            std::cout << ibuf[0] << " " << ibuf[1] << " ";
            std::cout << ibuf[2] << " " << ibuf[3] << " ";
            std::cout << "0 0 0 0";
            std::cout << std::endl;
          }

          // increase the grid elem counter
          numTetElems_[gridIdx]++;

          break;

        case 4:

          // hexahedron: 8 nodes

          fileStream_.read(cbuf, 8*sizeof(int));
          ibuf = reinterpret_cast<int*>(cbuf);

          // store the nodes
          for( int j = 0; j < 8; j++ ) {
            hexElemNodes_[gridIdx][(8*numHexElems_[gridIdx])+j] = ibuf[j];
          }

          if( ELEMENTVERBOSE ) {
            std::cout << "Element " << i << ": ";
            std::cout << ibuf[0] << " " << ibuf[1] << " ";
            std::cout << ibuf[2] << " " << ibuf[3] << " ";
            std::cout << ibuf[4] << " " << ibuf[5] << " ";
            std::cout << ibuf[6] << " " << ibuf[7];
            std::cout << std::endl;
          }

          // increase the grid elem counter
          numHexElems_[gridIdx]++;

          break;

        case 3:

          // prism: 6 nodes

          fileStream_.read(cbuf, 6*sizeof(int));
          ibuf = reinterpret_cast<int*>(cbuf);

          // store the nodes
          for( int j = 0; j < 6; j++ ) {
            priElemNodes_[gridIdx][(6*numPriElems_[gridIdx])+j] = ibuf[j];
          }

          if( ELEMENTVERBOSE ) {
            std::cout << "Element " << i << ": ";
            std::cout << ibuf[0] << " " << ibuf[1] << " ";
            std::cout << ibuf[2] << " " << ibuf[3] << " ";
            std::cout << ibuf[4] << " " << ibuf[5] << " ";
            std::cout << "0 0";
            std::cout << std::endl;
          }

          // increase the grid elem counter
          numPriElems_[gridIdx]++;

          break;

        case 2:

          // pyramid: 5 nodes

          fileStream_.read(cbuf, 5*sizeof(int));
          ibuf = reinterpret_cast<int*>(cbuf);

          // store the nodes
          for( int j = 0; j < 5; j++ ) {
            pyrElemNodes_[gridIdx][(5*numPyrElems_[gridIdx])+j] = ibuf[j];
          }

          if( ELEMENTVERBOSE ) {
            std::cout << "Element " << i << ": ";
            std::cout << ibuf[0] << " " << ibuf[1] << " ";
            std::cout << ibuf[2] << " " << ibuf[3] << " ";
            std::cout << ibuf[4] << " 0 0 0";
            std::cout << std::endl;
          }

          // increase the grid elem counter
          numPyrElems_[gridIdx]++;

          break;

        default:

          std::cout << "ERROR: unrecognized element: " << elemType << std::endl;

          break;

        }

      }

      if( DEBUG ) std::cout << "Section standard 3-D elements (FV_ELEMENTS) finished" << std::endl;

      break; }

    case FV_ARB_POLY_ELEMENTS:

      if( DEBUG ) std::cout << "Section arbitrary polyhedron elements (FV_ARB_POLY_ELEMENTS)" << std::endl;

      std::cout << "ERROR: FV_ARB_POLY_ELEMENTS section" << std::endl;

      {

        fileStream_.read(cbuf, 1*sizeof(int));
        ibuf = reinterpret_cast<int*>(cbuf);
        int numElems = ibuf[0];

        for( int i = 0; i < numElems; i++ ) {

          fileStream_.read(cbuf, 3*sizeof(int));
          ibuf = reinterpret_cast<int*>(cbuf);
          int numFaces = ibuf[0];

          for( int j = 0; j < numFaces; j++ ) {
            fileStream_.read(cbuf, 2*sizeof(int));
            ibuf = reinterpret_cast<int*>(cbuf);
            int numVertices = ibuf[1];

            fileStream_.read(cbuf, (numVertices+1)*sizeof(int));
            ibuf = reinterpret_cast<int*>(cbuf);
            int numHangNodes = ibuf[numVertices];

            if( numHangNodes > 0 ) {
              fileStream_.read(cbuf, numHangNodes*sizeof(int));
            }

          }

        }

      }

      if( DEBUG ) std::cout << "Section arbitrary polyhedron elements (FV_ARB_POLY_ELEMENTS) finished" << std::endl;

      break;

    case FV_VARIABLES:

      // elements section is finished
      isFinished = true;

      // rewind ifstream for next section
      fileStream_.seekg(-1 * sizeof(int), std::ios_base::cur);

      break;

    case FV_NODES:

      // elements section is finished, new grid
      isFinished = true;

      // rewind ifstream for next section
      fileStream_.seekg(-1 * sizeof(int), std::ios_base::cur);

      break;

    default:

      std::cout << "ERROR: Expected FV_ELEMENTS/FV_ARB_POLY_ELEMENTS/FV_VARIABLES/FV_NODES, found " << ibuf[0] << std::endl;
      return;

      break;

    }

  }

}

void FieldViewFile::readFaces(int gridIdx) {

  if( TRACE ) std::cout << "TRACE: FieldViewFile::readFacesSection()" << std::endl;

  if ( fileType_ == FV_RESULTS_FILE ) {
    // does not contain faces sections
    std::cout << "readFacesSection() was called for FV_RESULTS_FILE" << std::endl;
    return;
  }

  // buffer for ifstream.read()
  char cbuf[80];
  // reinterpret_cast target for ints
  int *ibuf;
  // whether the faces sections are finished
  bool isFinished = false;
  // number of faces per section
  int numFaces;

  // init some vectors that have the size of grids
  numResultFacesStd_.resize(numGrids_);
  numResultFacesStd_.at(gridIdx) = 0;
  numResultFacesPoly_.resize(numGrids_);
  numResultFacesPoly_.at(gridIdx) = 0;

  while( isFinished == false ) {

    fileStream_.read(cbuf, sizeof(int));
    ibuf = reinterpret_cast<int*>(cbuf);

    switch( ibuf[0] ) {

    case FV_FACES:

      // number of faces in this section
      fileStream_.read(cbuf, sizeof(int) * 2);
      ibuf = reinterpret_cast<int*>(cbuf);
      numFaces = ibuf[1];
      if( DEBUG ) std::cout << numFaces << " faces belonging to " << bndryNames_.at(ibuf[0]-1) << std::endl;

      // if the boundary the faces belong to has results, boundary variables will follow
      if( bndryResultFlag_.at(ibuf[0]-1) == 1 ) {
        numResultFacesStd_.at(gridIdx) += numFaces;
        if( DEBUG ) std::cout << "numResultFacesStd_.at(" << gridIdx << ") = " << numResultFacesStd_.at(gridIdx) << std::endl;
      }

      // skip the information
      fileStream_.seekg(numFaces * sizeof(int) * 4, std::ios_base::cur);

      break;

    case FV_ARB_POLY_FACES:

      // number of faces in this section
      fileStream_.read(cbuf, sizeof(int) * 2);
      ibuf = reinterpret_cast<int*>(cbuf);
      numFaces = ibuf[1];
      if( DEBUG ) std::cout << numFaces << " faces belonging to " << bndryNames_.at(ibuf[0]-1) << std::endl;

      // if the boundary the faces belong to has results, boundary variables will follow
      if( bndryResultFlag_.at(ibuf[0]-1) == 1 ) {
        numResultFacesPoly_.at(gridIdx) += numFaces;
      }

      // skip each face
      for( int i = 0; i < numFaces; i++ ) {

        fileStream_.read(cbuf, sizeof(int));
        ibuf = reinterpret_cast<int*>(cbuf);
        fileStream_.seekg(ibuf[0] * sizeof(int), std::ios_base::cur);

      }

      break;

    case FV_ELEMENTS:

      // faces section is finished
      isFinished = true;

      // rewind ifstream for next section
      fileStream_.seekg(-1 * sizeof(int), std::ios_base::cur);

      break;

    default:

      std::cout << "ERROR: Expected FV_FACES/FV_ARB_POLY_FACES/FV_ELEMENTS, found " << ibuf[0] << std::endl;
      return;

      break;

    }

  }

}

void FieldViewFile::readCoordinates(int gridIdx, bool isIgnored) {

  if( TRACE ) std::cout << "TRACE: FieldViewFile::readCoordinates()" << std::endl;

  if( fileType_ == FV_RESULTS_FILE ) {

    // this file does not contain face information
    if( DEBUG ) std::cout << "FV_RESULTS_FILE, readCoordinates() returns" << std::endl;
    return;

  }

  // buffer for ifstream.read()
  char *cbuf;
  // reinterpret_cast target for floats
  float *fbuf;

  if( isIgnored == true ) {

    // data should be ignored, just modify the pointer
    if( DEBUG ) std::cout << "Ignoring node coordinates" << std::endl;
    // for each entry, there are numNodes_[gridIdx] * 3 floats
    int offset = sizeof(float)* numNodes_[gridIdx] * 3;
    fileStream_.seekg(offset, std::ios_base::cur);

  } else {

    // allocate memory
    nodeCoords_[gridIdx] = new float[numNodes_[gridIdx] * 3];

    // read data
    cbuf = new char[sizeof(float)* numNodes_[gridIdx] * 3];
    fileStream_.read(cbuf, sizeof(float)* numNodes_[gridIdx] * 3);
    fbuf = reinterpret_cast<float*>(cbuf);

    // copy array
    for( int i = 0; i < numNodes_[gridIdx] * 3; i++ ) {
      nodeCoords_[gridIdx][i] = fbuf[i];
    }

    if( NODEVERBOSE ) {

      for( int nodeIdx = 0; nodeIdx < numNodes_[gridIdx]; nodeIdx++ ) {
        std::cout << "Node " << nodeIdx << ": " << nodeCoords_[gridIdx][numNodes_[gridIdx]*0+nodeIdx] << ", ";
        std::cout << nodeCoords_[gridIdx][numNodes_[gridIdx]*1+nodeIdx] << ", " << nodeCoords_[gridIdx][numNodes_[gridIdx]*2+nodeIdx] << std::endl;
      }
    }

  }

}

void FieldViewFile::readGrids(bool ignoreGrid, bool ignoreVars) {

  if( TRACE ) std::cout << "TRACE: FieldViewFile::readGrids()" << std::endl;

  // buffer for ifstream.read()
  char cbuf[80];
  // reinterpret_cast target for ints
  int *ibuf;

  // init the array for the grid nodes
  numNodes_ = new int[numGrids_];

  // init the array for the node coordinates if necessary
  if( ignoreGrid == false ) {
    nodeCoords_ = new float*[numGrids_];
  }

  // init the arrays for the elements if necessary
  if( ignoreGrid == false ) {
    numTetElems_ = new int[numGrids_];
    std::fill(numTetElems_, numTetElems_ + numGrids_, 0);
    numHexElems_ = new int[numGrids_];
    std::fill(numHexElems_, numHexElems_ + numGrids_, 0);
    numPriElems_ = new int[numGrids_];
    std::fill(numPriElems_, numPriElems_ + numGrids_, 0);
    numPyrElems_ = new int[numGrids_];
    std::fill(numPyrElems_, numPyrElems_ + numGrids_, 0);
    tetElemNodes_ = new int*[numGrids_];
    hexElemNodes_ = new int*[numGrids_];
    priElemNodes_ = new int*[numGrids_];
    pyrElemNodes_ = new int*[numGrids_];
  }

  // init the variables array
  if( ignoreVars == false ) {
    nodeVariables_ = new float*[numGrids_];
  }

  for( int currentGrid = 0; currentGrid < numGrids_; currentGrid++ ) {

    if( VERBOSE ) std::cout << "-----------------------------------------------" << std::endl;
    if( VERBOSE ) std::cout << "Now processing grid " << currentGrid+1 << "..." << std::endl;

    // FV_NODES section
    fileStream_.read(cbuf, sizeof(int));
    ibuf = reinterpret_cast<int*>(cbuf);

    // check if section is correct
    if( ibuf[0] != FV_NODES ) {

      std::cout << "ERROR: could not read FV_NODES section" << std::endl;

      if( DEBUG ) std::cout << "Expected " << FV_NODES << ", found " << ibuf[0] << std::endl;

      isValidFile_ = false;
      return;

    }

    if( DEBUG ) std::cout << "Section FV_NODES" << std::endl;

    // read number of nodes in this grid
    fileStream_.read(cbuf, sizeof(int));
    ibuf = reinterpret_cast<int*>(cbuf);
    numNodes_[currentGrid] = ibuf[0];

    if( VERBOSE ) std::cout << "Number of nodes: " << numNodes_[currentGrid] << std::endl;

    readCoordinates(currentGrid, ignoreGrid);

    if( DEBUG ) std::cout << "Section FV_NODES complete" << std::endl;

    if ( fileType_ != FV_RESULTS_FILE ) {

      // FV_FACES / FV_ARB_POLY_FACES
      if( DEBUG ) std::cout << "Section face information (FV_FACES/FV_ARB_POLY_FACES)" << std::endl;
      readFaces(currentGrid);
      if( DEBUG ) std::cout << "Section face information (FV_FACES/FV_ARB_POLY_FACES) finished" << std::endl;

      // FV_ELEMENTS / FV_ARB_POLY_ELEMENTS
      if( DEBUG ) std::cout << "Section element information (FV_ELEMENTS/FV_ARB_POLY_ELEMENTS)" << std::endl;
      readElements(currentGrid, ignoreGrid);
      if( DEBUG ) std::cout << "Section element information (FV_ELEMENTS/FV_ARB_POLY_ELEMENTS) finished" << std::endl;

    }

    // check if end of file has been reached
    if( fileStream_.eof() ) return;

    if( fileType_ != FV_GRIDS_FILE ) {

      // FV_VARIABLES
      readVariables(currentGrid, ignoreVars);

      // FV_BNDRY_VARS
      readBoundaryVariables(currentGrid);

    }

    // check if end of file has been reached
    if( fileStream_.eof() ) return;

  }

}

void FieldViewFile::readHeader() {


  if( VERBOSE ) std::cout << std::endl << "-----------------------------------------------" << std::endl;

  if( TRACE ) std::cout << "TRACE: FieldViewFile::readHeader()" << std::endl;

  // buffer for ifstream.read()
  char cbuf[80];
  // reinterpret_cast target for ints
  int *ibuf;
  // reinterpret_cast target for floats
  float *fbuf;

  // FV_MAGIC (bit pattern to identify the file as FIELDVIEW file)
  fileStream_.read(cbuf, sizeof(int));
  ibuf = reinterpret_cast<int*>(cbuf);

  if( ibuf[0] != FV_MAGIC ) {

    std::cout << "ERROR: failed to read FV_MAGIC" << std::endl;

    if( DEBUG ) {
      std::cout << "ERROR: Expected " << FV_MAGIC << ", found " << ibuf[0] << std::endl;
    }

    isValidFile_ = false;
    return;

  } else {

    if( DEBUG ) {
      std::cout << "FV_MAGIC found: " << ibuf[0] << std::endl;
    }

  }

  // FIELDVIEW string
  fileStream_.read(cbuf, sizeof(char)*80);
  cbuf[9] = '\0';

  if( strncmp(cbuf, "FIELDVIEW", 9) != 0 ) {

    std::cout << "ERROR: failed to read 'FIELDVIEW' string" << std::endl;

    if( DEBUG ) {
      std::cout << "ERROR: Expected 'FIELDVIEW', found '" << cbuf << "'" << std::endl;
    }

    isValidFile_ = false;
    return;

  } else {

    if( DEBUG ) {
      std::cout << "'FIELDVIEW' string found" << std::endl;
    }

  }

  // file version numbers
  fileStream_.read(cbuf, sizeof(int)*2);
  ibuf = reinterpret_cast<int*>(cbuf);

  if( ibuf[0] != 3 || ibuf[1] != 0 ) {

    std::cout << "ERROR: failed to read version numbers" << std::endl;

    if( DEBUG ) {
      std::cout << "ERROR: Expected 3 0, found " << ibuf[0] << " " << ibuf[1] << std::endl;
    }

    isValidFile_ = false;
    return;

  } else {

    if( DEBUG ) {
      std::cout << "Version numbers found: " << ibuf[0] << " " << ibuf[1] << std::endl;
    }

  }


  // file type
  fileStream_.read(cbuf, sizeof(int));
  ibuf = reinterpret_cast<int*>(cbuf);
  fileType_ = ibuf[0];

  if( fileType_ == FV_GRIDS_FILE ) {

    std::cout << "File is of type FV_GRIDS_FILE" << std::endl;

  } else if( fileType_ == FV_RESULTS_FILE ) {

    std::cout << "File is of type FV_RESULTS_FILE" << std::endl;

  } else if( fileType_ == FV_COMBINED_FILE ) {

    std::cout << "File is of type FV_COMBINED_FILE" << std::endl;

  } else {

    std::cout << "Error: could not recognize file type" << std::endl;

    if( DEBUG ) {
      std::cout << "File type field contains: " << fileType_ << std::endl;
    }

    isValidFile_ = false;
    return;

  }

  // reserved field containing 0
  fileStream_.read(cbuf, sizeof(int));
  ibuf = reinterpret_cast<int*>(cbuf);

  if( ibuf[0] != 0 ) {

    std::cout << "ERROR: failed to read reserved field" << std::endl;

    if( DEBUG ) {
      std::cout << "Expected 0, found " << ibuf[0] << std::endl;
    }

    isValidFile_ = false;
    return;

  } else {

    if( DEBUG ) {
      std::cout << "Found reserved field: " << ibuf[0] << std::endl;
    }

  }

  // solution time and 3 constants
  // only solution time is stored, no use for the other three
  if( fileType_ != FV_GRIDS_FILE ) {

    fileStream_.read(cbuf, sizeof(float)*4);
    fbuf = reinterpret_cast<float*>(cbuf);
    solutionTime_ = fbuf[0];

    if( VERBOSE ) std::cout << "Time: " << solutionTime_ << std::endl;

  }

  // number of grids
  fileStream_.read(cbuf, sizeof(int));
  ibuf = reinterpret_cast<int*>(cbuf);

  if( ibuf[0] < 1 ) {

    std::cout << "ERROR: failed to read number of grids" << std::endl;

    if( DEBUG ) {
      std::cout << "Found value of " << ibuf[0] << std::endl;
    }

    isValidFile_ = false;
    return;

  } else {

    numGrids_ = ibuf[0];
    std::cout << "Number of grids: " << numGrids_ << std::endl;

  }

  // if control reaches this point, it is a valid file
  isValidFile_ = true;

  if( fileType_ != FV_RESULTS_FILE ) {

    // number of face types
    int numFaceTypes;

    // number of face types
    fileStream_.read(cbuf, sizeof(int));
    ibuf = reinterpret_cast<int*>(cbuf);
    numFaceTypes = ibuf[0];

    if( DEBUG ) std::cout << "Number of boundaries: " << numFaceTypes << std::endl;

    // resize vector for result flags
    bndryResultFlag_.resize(numFaceTypes);
    bndryNames_.resize(numFaceTypes);

    // walk the face types
    for( int i = 0; i < numFaceTypes; i++ ) {

      // read the face type data
      fileStream_.read(cbuf, sizeof(int)*2);
      ibuf = reinterpret_cast<int*>(cbuf);

      // store the result flag
      bndryResultFlag_.at(i) = ibuf[0];

      fileStream_.read(cbuf, sizeof(char)*80);
      bndryNames_.at(i) = std::string(cbuf);
      if( DEBUG ) std::cout << cbuf << " (results: " << bndryResultFlag_.at(i) << ")" << std::endl;

    }

  }

  // nodal and boundary variables
  if( fileType_ != FV_GRIDS_FILE ) {

    // number of nodal variables
    fileStream_.read(cbuf, sizeof(int));
    ibuf = reinterpret_cast<int*>(cbuf);
    numNodalVars_ = ibuf[0];
    nodalVarNames_.resize(numNodalVars_);

    if( VERBOSE ) std::cout << "Number of nodal variables: " << numNodalVars_ << std::endl;

    // walk the variables
    for( int i = 0; i < numNodalVars_; i++ ) {

      fileStream_.read(cbuf, sizeof(char)*80);
      nodalVarNames_.at(i) = std::string(cbuf);
      if( VERBOSE ) std::cout << cbuf << std::endl;

    }

    // number of nodal variables
    fileStream_.read(cbuf, sizeof(int));
    ibuf = reinterpret_cast<int*>(cbuf);
    numBndryVars_ = ibuf[0];

    if( DEBUG ) std::cout << "Number of boundary variables: " << numBndryVars_ << std::endl;

    // walk the variables
    for( int i = 0; i < numBndryVars_; i++ ) {

      fileStream_.read(cbuf, sizeof(char)*80);
      if( DEBUG ) std::cout << cbuf << std::endl;

    }

  }

}

} /* namespace CoupledField */


