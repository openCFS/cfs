/** This is the reader for HDF5 files from CFS++. */
#ifndef __vtkCFSReader_h
#define __vtkCFSReader_h


#include "vtkSetGet.h"
#include "vtkObjectBase.h"
#include "vtkWin32Header.h"
#include "vtkMultiBlockDataSetAlgorithm.h"
#include "vtkUnstructuredGridAlgorithm.h"
#include "vtkCellType.h"
#include "vtkMultiBlockDataSet.h"

#include "vtkIOCFSReaderModule.h" // For export macro

#include "hdf5Reader.h"

class vtkDoubleArray;
class vtkCell;
class VTKIOCFSREADER_EXPORT vtkCFSReader : public vtkMultiBlockDataSetAlgorithm
{
public:
  static vtkCFSReader *New();
  vtkTypeMacro(vtkCFSReader,vtkMultiBlockDataSetAlgorithm)
  void PrintSelf(ostream& os, vtkIndent indent);

  // Check, if the current file is really a CFS hdf5 file.
  
#if 0
  // This method gets called in order to determine the "right" reader in case
  // several readers can open files with the same extension.
  virtual int CanReadFile(const char* fname);
#endif
  
  // Description:
  // Specify the file name of the CFS data file to read.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  // Description:
  // Get the multisequence range. Filled during RequestInformation.
  vtkGetVector2Macro(MultiSequenceRange, int);

  // Description:
  // Set/Get the current multisequence step
  void SetMultiSequenceStep(int step);
  //vtkSetMacro(MultiSequenceStep, int);
  vtkGetMacro(MultiSequenceStep, int);
    
  // Description:
  // Set/Get the current timestep index
  virtual void SetTimeStep(unsigned int step);
  vtkGetMacro(TimeStep, unsigned int);

  // Description:
  // Get current frequency / time value
  virtual const char* GetTimeFreqValStr(); 

  // Description:
  // Add dimensions to array names
  virtual void SetAddDimensionsToArrayNames(int);
  vtkGetMacro(AddDimensionsToArrayNames, int);
  vtkBooleanMacro(AddDimensionsToArrayNames, int);

    // Description:
  // Add dimensions to array names
  virtual void SetHarmonicDataAsModeShape(int);
//  vtkSetMacro(HarmonicDataAsModeShape, int);
  vtkGetMacro(HarmonicDataAsModeShape, int);
  vtkBooleanMacro(HarmonicDataAsModeShape, int);

  // Description:
  // Get the timesteprange. Filled during RequestInformation.
  vtkGetVector2Macro(TimeStepNumRange, int);
 
  virtual double *GetTimeStepRange ()
  {
    vtkDebugMacro(<< this->GetClassName() << " (" << this << "): returning TimeStepRange pointer "
                  << this->TimeStepValuesRange);
    
    return this->TimeStepValuesRange;
  }
  
  virtual void GetTimeStepRange (double &_arg1, double &_arg2)
  {
    _arg1 = this->TimeStepValuesRange[0];
    _arg2 = this->TimeStepValuesRange[1];
    vtkDebugMacro(<< this->GetClassName() << " (" << this << "): returning TimeStepRange = ("
                  << _arg1 << "," << _arg2 << ")"); 
  };
  
  virtual void GetTimeStepRange (double _arg[2])
  {
    this->GetTimeStepRange (_arg[0], _arg[1]);
  } 

  // Description:
  // Type of complex result treating
  // 0 = REAL/IMAG
  // 1 = AMPL/PHASE
  vtkSetMacro(ComplexMode, int);
  vtkGetMacro(ComplexMode, int);
  
  //! switch if missing results should get filled with 0
  //! 0 = omit empty regions (only partial results available)
  //! 1 = fill empty missing results with 0-valued vector
  void SetFillMissingResults(int);
  vtkGetMacro(FillMissingResults, int);
  vtkBooleanMacro(FillMissingResults, int)

  // Description:
  // Get the number of cell arrays available in the input.
  int GetNumberOfRegionArrays(void);

  // Description:
  // Get/Set whether the cell array with the given name is to
  // be read.
  int GetRegionArrayStatus(const char *name);
  void SetRegionArrayStatus(const char *name, int status);

  // Description:
  // Get the name of the  cell array with the given index in
  // the input.
  const char *GetRegionArrayName(int index);

  // Description:
  // Get the number of named node arrays available in the input.
  int GetNumberOfNamedNodeArrays(void);

  // Description:
  // Get/Set whether the named node array with the given name is to
  // be read.
  int GetNamedNodeArrayStatus(const char *name);
  void SetNamedNodeArrayStatus(const char *name, int status);

  // Description:
  // Get the name of the  named node array with the given index in
  // the input.
  const char *GetNamedNodeArrayName(int index);

  // Description:
  // Get the number of named element arrays available in the input.
  int GetNumberOfNamedElemArrays(void);

  // Description:
  // Get/Set whether the named elem array with the given name is to
  // be read.
  int GetNamedElemArrayStatus(const char *name);
  void SetNamedElemArrayStatus(const char *name, int status);

  // Description:
  // Get the name of the  named elem array with the given index in
  // the input.
  const char *GetNamedElemArrayName(int index);
  
protected:
  vtkCFSReader();
  ~vtkCFSReader();
  int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  char *FileName;

private:
  
  // Functions or variables, which should not be Python-wrapped go between BTX and ETX
  //BTX
  
  //! hdf5 reader
  H5CFS::Hdf5Reader in_;
  
  //! map cfs/hdf5 element type to vtk native ones
  VTKCellType GetCellIdType( H5CFS::ElemType type);
  
  
  //! save vector into vtkDoubleArray
  vtkDoubleArray* SaveToArray( const std::vector<double>& vals,
                               unsigned int numDofs,
                               unsigned int numEntities,
                               H5CFS::EntryType entryType,
                               const std::string& name );

  //! Pretty print frequency including suffix (Hz, kHz, MHz...)
  void PrettyPrintFreq(double freq, std::string& str);
  
  //! Pretty print frequency including suffix (ns, ms, s)
  void PrettyPrintTime(double time, std::string& str);
		       
  //! vector containing all region names
  std::vector<std::string> regionNames_;

  //! vector containing all named nodes
  std::vector<std::string> namedElemNames_;

  //! vector containing all named elems
  std::vector<std::string> namedNodeNames_;
  
  //! map (region, globalNodeNumber)->(regionLocalNodeNumber)
  std::vector<std::vector< unsigned int> > nodeMap_;
  
  //! time steps of current multisequence step
  std::vector<double> stepVals;
  
  //! vector with step values having a result
  std::vector<unsigned int> stepNums; 
  
  //! multiblock dataset, which contains for each
  //! region a initialized grid
  vtkMultiBlockDataSet * mbDataSet;
  
  //! reduced multiblock dataset, which contains
  //! only blocks for currently active regions
  vtkMultiBlockDataSet * mbActiveDataSet;
  
  //! flag array indicating active regions
  std::map<std::string, int> regionSwitch;

  //! flag array indicating named nodes to read
  std::map<std::string, int> namedNodeSwitch;

  //! flag array indicating named elems to read
  std::map<std::string, int> namedElemSwitch;
  
  //! flag indicating change of multisequence step
  bool msStepChanged;
  
  //ETX
  
  //! current multisequence step
  unsigned int MultiSequenceStep;
  
  //! analysis type of current multisequence step
  H5CFS::AnalysisType analysisType;
  
  //! current time step (for discrete time steps)
  unsigned int TimeStep;
  
  //! current time step / frequency value pretty printed
  double timeFreqVal;
  
  //! pretty formatted time value / frequency value
  std::string timeFreqValStr;
  
  //! The time value which the pipeline requests
  double requestedTimeValue;
  
  //! complex mode
  //! 0 = real-imag
  //! 1 = ampl-phase
  unsigned int ComplexMode;
  
  //! switch if missing results should get filled with 0
  //! 0 = omit empty regions (only partial results available)
  //! 1 = fill empty missing results with 0-valued vector
  unsigned int FillMissingResults; 
  
  // add dimensions to array names
  int AddDimensionsToArrayNames;

  // Interpret harmonic data as transient mode shape
  int HarmonicDataAsModeShape;
  
  // Are we dealing with harmonic data?
  bool harmonicData;
  
  //! number of time/freq steps in current ms step
  unsigned int NumberOfTimeSteps;
  
  //! array with first/last time/freq step value
  int TimeStepNumRange[2];

  //! array with first/last time/freq step value
  double TimeStepValuesRange[2];
  
  //! array with first/last ms step number
  int MultiSequenceRange[2];
  
  //! flag if hdf5 file is already read in
  bool isInitialized;

  //! flag if hdf5 infos (region names, named node/elems) have been read in
  bool hdf5InfoRead;
  
  //! flag, if region settings were modified
  bool activeRegionsChanged;
  
  //! flag, to reset the data value arrays (e.g. after
  //! switching the complex mode etc.)
  bool resetDataArrays;
  
  
  void ReadFile(vtkMultiBlockDataSet *output);
  void ReadNodes(vtkMultiBlockDataSet *output);
  void ReadCells(vtkMultiBlockDataSet *output);
  
  void ReadNodeCellData(vtkMultiBlockDataSet *output, bool isNode);

  void UpdateActiveRegions(); 
  
  vtkCFSReader(const vtkCFSReader&);  // Not implemented.
  void operator=(const vtkCFSReader&);  // Not implemented.
};

#endif
