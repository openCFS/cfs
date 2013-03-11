/*
 * FieldViewFile.h
 *
 *  Created on: 27.09.2011
 *      Author: willmitzer
 */

#ifndef FIELDVIEWFILE_H_
#define FIELDVIEWFILE_H_

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

namespace CoupledField {

class FieldViewFile {
public:

  FieldViewFile();
  ~FieldViewFile();
  void open(std::string fileName);
  void openGrid(std::string fileName);
  void openVars(std::string fileName);

  int getNumGrids();
  int getNumNodes(int gridIndex);
  int getNumElements(int gridIndex);
  float* getNodeCoordinates(int gridIndex);
  std::string getFileName();
  int getNumTetElems(int gridIdx);
  int getNumHexElems(int gridIdx);
  int getNumPriElems(int gridIdx);
  int getNumPyrElems(int gridIdx);
  int* getTetElems(int gridIdx);
  int* getHexElems(int gridIdx);
  int* getPriElems(int gridIdx);
  int* getPyrElems(int gridIdx);
  std::vector<std::string> getNodalVarNames();
  float* getNodalVars(int gridIdx);

private:

  // name of fvuns file
  std::string fileName_;

  // associated ifstream
  std::ifstream fileStream_;

  // indicates whether the file is of valid format
  bool isValidFile_;

  // file data

  // FVUNS format knows three different file types (see fv_reader_tags.h)
  int fileType_;
  // number of grids in file (one partition per grid in cplreader)
  int numGrids_;
  float solutionTime_;

  // variable data

  // number of vars for each node
  int numNodalVars_;
  // number of vars for each boundary node
  int numBndryVars_;
  std::vector<std::string> nodalVarNames_;
  // size is numGrids_
  float **nodeVariables_;

  // node data

  // size is numGrids_
  float **nodeCoords_;
  int *numNodes_;

  // element data

  // size is numGrids_
  int *numTetElems_;
  int *numHexElems_;
  int *numPriElems_;
  int *numPyrElems_;
  int **tetElemNodes_;
  int **hexElemNodes_;
  int **priElemNodes_;
  int **pyrElemNodes_;

  // boundary data

  // contains the results flag for boundary types (boundaries)
  std::vector<int> bndryResultFlag_;
  std::vector<std::string> bndryNames_;

  // number of faces per grid which have results
  // size is numGrids_
  std::vector<int> numResultFacesStd_;
  std::vector<int> numResultFacesPoly_;

  // functions that read part of the file
  void readHeader();
  void readGrids(bool ignoreGrid, bool ignoreVars);
  void readCoordinates(int gridIdx, bool ignoreGrid);
  void readFaces(int gridIdx);
  void readElements(int gridIdx, bool ignoreGrid);
  void readVariables(int gridIdx, bool ignoreVars);
  void readBoundaryVariables(int gridIdx);

};

} /* namespace CoupledField */

#endif /* FIELDVIEWFILE_H_ */
