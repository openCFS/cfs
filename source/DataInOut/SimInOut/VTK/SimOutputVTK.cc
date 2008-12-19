// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdio.h>
#include <complex>
#include <ctime>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/exception.hpp>
namespace fs = boost::filesystem;

#include <def_cfs_stats.hh>

#include <DataInOut/SimInOut/VTK/SimOutputVTK.hh>
#include <DataInOut/simInput.hh>
#include <DataInOut/Logging/cfslog.hh>
#include <DataInOut/ParamHandling/ParamNode.hh>
#include <Domain/entityList.hh>


// ------- include VTK stuff
// this stuff from OLAS prevents VTK to be properly included! :(
#undef New
#undef DeleteArray
// VTK assuses this within it's includes
using namespace std;
#include <vtkCellType.h>
#include <vtkXMLDataElement.h>
#include <vtkPointSet.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkXMLMultiBlockDataWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkZLibDataCompressor.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkXMLDataElement.h>
#include <vtkXMLUtilities.h>


using namespace CoupledField;

// declare logging stream
DECLARE_LOG(vtk)
DEFINE_LOG(vtk, "SimOutputVTK")

SimOutputVTK::SimOutputVTK(const std::string& fileName, ParamNode * outputNode ) :  SimOutput(fileName, outputNode)
{
  // Initialize variables
  formatName_ = "VTK";
  capabilities_.insert( MESH );
  capabilities_.insert( MESH_RESULTS );

  stepNumOffset_ = 0;
  stepValOffset_ = 0.0;
  dirName_ = ".";
  fileName_ = fileName;

  writeSingleFiles_ = true;
}

SimOutputVTK::~SimOutputVTK() 
{
  for(unsigned int i = 0; i < data_.GetSize(); i++)
    data_[i].grid->Delete();
}

void SimOutputVTK::Init( Grid* ptGrid, bool printGridOnly ) 
{
  ptGrid_ = ptGrid;

  // init own structure
  ptGrid_->GetRegionIds(regions);

  // allocate our storage data
  data_.Resize(regions.GetSize());
  
  // a region represents an unstructred grid with vtk cells = cfs elements
  for(unsigned int i = 0; i < data_.GetSize(); i++)
  {
    assert(regions[i] == (int) i);
    RegionData& data = data_[i];
    data.regionId = i;
    CreateMesh(data);
  }
}


void SimOutputVTK::CreateMesh(RegionData& data)
{
  data.grid = vtkUnstructuredGrid::New();

  // TODO: what to do with the name? 
  //std::string name = ptGrid_->RegionIdToName(regionId);
  
  // handle nodes/points
  StdVector<unsigned int> nodes;
  ptGrid_->GetNodesByRegion(nodes, data.regionId);
  
  std::map<unsigned int, unsigned int> n2p;

  // store the nodes in vtk points 
  vtkPoints* points = vtkPoints::New();
  // allocates
  points->SetNumberOfPoints(nodes.GetSize());
  
  // within the vtk grid we always start with 0
  Point p;
  for(unsigned int i = 0, n = nodes.GetSize(); i < n; i++)
  {
    unsigned int node_num = nodes[i];
    n2p[node_num] =  i; // store mapping
    
    ptGrid_->GetNodeCoordinate(p, node_num);
    points->SetPoint(i, p[0], p[1], p[2]);
  }
  
  data.grid->SetPoints(points);
  // according to the VTK example we can delete the points
  points->Delete();
  
  // now add the elements
  StdVector<Elem*> elems;
  ptGrid_->GetElems(elems,data.regionId);
  
  // reserve as we will add cell by cell
  data.grid->Allocate(elems.GetSize());
  
  // traverse element
  for(unsigned int i = 0, n = elems.GetSize(); i < n; i++)
  {
    Elem* elem = elems[i];
    
    // store here the nodes
    vtkIdType nodeIds[27]; // maximal size
    // set to the node number
    for(unsigned int i = 0, n = elem->connect.GetSize(); i < n; i++)
    {
      std::map<unsigned int, unsigned int>::iterator iter = n2p.find(elem->connect[i]);
      assert(iter != n2p.end());
      nodeIds[i] = iter->second; // vtk point id is 0 based within a grid (we assume!)
    }
    
    VTKCellType ct = (VTKCellType) GetCellIdType(elem->ptElem->feType());
    data.grid->InsertNextCell(ct, elem->connect.GetSize(), nodeIds);
  }
}



void SimOutputVTK::BeginMultiSequenceStep( UInt step,
    BasePDE::AnalysisType type,
    UInt numSteps )
{
  currAnalysis_ = type;
  currMsStep_ = step;
}

void SimOutputVTK::RegisterResult( shared_ptr<BaseResult> sol,
    UInt saveBegin, UInt saveInc,
    UInt saveEnd,
    bool isHistory )
{

  ResultInfo & actDof = *(sol->GetResultInfo());

  LOG_DBG(vtk) << "Registering output '" << actDof.resultName 
  << "' with saveBegin " << saveBegin
  << ", saveInc " << saveInc 
  << ", saveEnd " << saveEnd;
}


//! Begin single analysis step
void SimOutputVTK::BeginStep( UInt stepNum, Double stepVal )
{
  for(unsigned int i = 0; i < data_.GetSize(); i++)
  {
    data_[i].results.Clear();
    data_[i].grid->GetPointData()->Reset();
    data_[i].grid->GetCellData()->Reset();
  }

  actStep_ = stepNum;
  actStepVal_ = stepVal;
  if( currAnalysis_ == BasePDE::TRANSIENT ||
      currAnalysis_ == BasePDE::STATIC  ) { 
    actStep_ += stepNumOffset_;
    actStepVal_ += stepValOffset_;
  }

}


//! Add result to current step
void SimOutputVTK::AddResult( shared_ptr<BaseResult> sol )
{
  LOG_DBG(vtk) << "AddResult: " << sol->ToString() << " with name " << sol->GetEntityList()->GetName();
  
  // retrieve the region id from the name
  RegionIdType regionId = ptGrid_->RegionNameToId(sol->GetEntityList()->GetName());
  data_[regionId].results.Push_back(sol);
}


//! End single analysis step
void SimOutputVTK::FinishStep( )
{
  // iterate over all regions
  for(unsigned int r = 0; r < data_.GetSize(); r++)
  {
    RegionData& data = data_[r];
    
    // iterate over all solutions within this region
    for(unsigned int s = 0; s < data.results.GetSize(); s++)
    {
      shared_ptr<BaseResult> sol = data.results[s];

      // get result info object and results for current result type
      ResultInfo& actInfo = *(sol->GetResultInfo());

      if(!ValidateNodesAndElements(actInfo)) continue;    

      LOG_DBG(vtk) << "Writing result '" << actInfo.resultName << "'";

      if(sol->GetEntryType() == EntryType::DOUBLE)
        WriteData(data, sol, actStepVal_, false);
      else
        WriteData(data, sol, actStepVal_, true);
    }
  }

  if(writeSingleFiles_) 
    WriteSingleXMlFile();
}

void SimOutputVTK::WriteSingleXMlFile()
{
  // finally write the file
  vtkXMLWriter* writer = CreateXMlWriter();
  
  // determine filename
  ostringstream os;
  os << fileName_ << "_";
  os.width(5);
  os.fill('0');
  os << singleFiles.GetSize(); // in the end we push "us" back to this index
  os << "." << writer->GetDefaultFileExtension();
  string fn = "simoutput_vtk/" + os.str();
  writer->SetFileName(fn.c_str());
  LOG_DBG(vtk) << "WriteSingleXMlFile: name=" << writer->GetFileName();
  
  // set modes
  if(myParam_->Get("format")->AsString() == "ascii")
  {
    writer->SetDataModeToAscii();
  }
  else
  {
    writer->SetDataModeToAppended(); 
    
    if(myParam_->Get("format")->AsString() == "base64")
      writer->EncodeAppendedDataOn();
    else
     writer->EncodeAppendedDataOff(); // raw-binary
    
    // in contrast to the docu we seem to have a compressor as default
    LOG_DBG(vtk) << "default compressor: " << writer->GetCompressor()->GetClassName();
    assert(writer->GetCompressor()->IsA("vtkZLibDataCompressor"));
    
    vtkZLibDataCompressor* zipper = dynamic_cast<vtkZLibDataCompressor*>(writer->GetCompressor());

    int cl = myParam_->Get("compression")->AsInt(); // -1 is off
    LOG_DBG(vtk) << "default compression level: " <<  zipper->GetCompressionLevel() << " in [" 
                 << zipper->GetCompressionLevelMinValue() << " .. " << zipper->GetCompressionLevelMaxValue() << "]";
    zipper->SetCompressionLevel(cl);
    LOG_DBG(vtk) << "set compression level: " << zipper->GetCompressionLevel();
  }

  // write the stuff
  int res = writer->Write();
  if(res != 1)
  {
    string msg = "Could not write single VTK file '" + fn + "'";
    Warning(msg.c_str());
  }

  // save what we did
  SingleVTUFile item;
  item.name = os.str();
  item.num_label = actStepVal_;
  item.id = singleFiles.GetSize();
  item.an_step  = actStep_;
  item.ms_step = actMSStep_;
  
  singleFiles.Push_back(item);

  // clean
  writer->Delete();
}

vtkXMLWriter* SimOutputVTK::CreateXMlWriter()
{
  if(data_.GetSize() > 1)
  {
    vtkXMLMultiBlockDataWriter* writer = vtkXMLMultiBlockDataWriter::New();

    // here we hold out unstructured grids
    vtkMultiBlockDataSet* mb = vtkMultiBlockDataSet::New();
    mb->SetNumberOfBlocks(data_.GetSize());
    for(unsigned int i = 0; i < data_.GetSize(); i++)
    {
      mb->SetBlock(i, data_[i].grid);
    }

    // transfer data  
    writer->SetInput(mb);
    mb->Delete();
    return writer;
  }
  else
  {
    vtkXMLUnstructuredGridWriter* writer = vtkXMLUnstructuredGridWriter::New();
    writer->SetInput(data_[0].grid);
    return writer;
  }
}


void SimOutputVTK::WriteData(RegionData& data, shared_ptr<BaseResult> sol, double actStepVal, bool harmonic)
{
  ResultInfo& info = *(sol->GetResultInfo());  

  string unit = info.unit.size() > 0 ? "_(" + info.unit + ")" : "";
  
  // complex data is split into two vtk results. Either real and imag or ampl/phase.
  // also the imag part might be ignored.
  if(!harmonic)
  {
    // pure real 
    string label = info.resultName + unit;
    WriteSingleVTKData<double>(data, sol, actStepVal, label, harmonic, true);
    return;
  }
  // we are harmonic!
  if(IgnoreImaginaryPart(info))
  {
    string label = info.resultName + unit;
    WriteSingleVTKData<complex<double> >(data, sol, actStepVal, label, harmonic, true);
    return;
  }

  if(info.complexFormat == REAL_IMAG)
  {
    string label = info.resultName + "_real" + unit;
    WriteSingleVTKData<complex<double> >(data, sol, actStepVal, label, harmonic, true);

    label = info.resultName + "_imag" + unit;
    WriteSingleVTKData<complex<double> >(data, sol, actStepVal, label, harmonic, false);
    return;
  }
  else
  {
    string label = info.resultName + "_ampl" + unit;
    WriteSingleVTKData<complex<double> >(data, sol, actStepVal, label, harmonic, true);

    label = info.resultName + "_phase";
    WriteSingleVTKData<complex<double> >(data, sol, actStepVal, label, harmonic, false);
    return;
  }
}

template<class TYPE>
void SimOutputVTK::WriteSingleVTKData(RegionData& data, shared_ptr<BaseResult> sol, double actStepVal, string& label, bool harmonic, bool real)
{
  ResultInfo& actInfo = *(sol->GetResultInfo());  

   // our data Array
  vtkDoubleArray* out = vtkDoubleArray::New();
  // meta data
  out->SetName(label.c_str());
 
  TransferData<TYPE>(actInfo, out, sol->GetCFSVector(), harmonic, real); 
  
  // did we hande node/point data or element/cell data?
  // nodes or elements
  bool nodes = actInfo.definedOn == ResultInfo::NODE || actInfo.definedOn == ResultInfo::PFEM;
  if(nodes)
    data.grid->GetPointData()->AddArray(out);
  else
    data.grid->GetCellData()->AddArray(out);
  
  out->Delete();
}


template<class TYPE>
void SimOutputVTK::TransferData(ResultInfo& info, vtkDoubleArray* out, CFSVector* sol_values, bool harmonic, bool real)
{
  // extract data -> the ordering is ok!
  Vector<TYPE>& values = dynamic_cast<Vector<TYPE>& >(*sol_values);

  
  unsigned int dofs = info.dofNames.GetSize(); 

  // dimension of stuff to be written
  unsigned int out_dofs = info.resultType == MECH_DISPLACEMENT ? 3 : dofs;

  unsigned int entities = values.GetSize() / dofs; 
  
  // tuple size 
  out->SetNumberOfComponents(out_dofs);
  out->SetNumberOfTuples(entities);
  
  // data
  StdVector<double> entity;
  entity.Resize(out_dofs);
  
  // transfer data
  for(unsigned int e = 0; e < entities; e++)
  {
    for(unsigned int d = 0; d < dofs; d++)
    {
      if(!harmonic)
      {
        entity[d] = ((complex<double>)values[e * dofs + d]).real();
      }
      else
      {
        if(info.complexFormat == REAL_IMAG)
        {
          if(real)
            entity[d] = ((complex<double>)values[e * dofs + d]).real();
          else
            entity[d] = ((complex<double>)values[e * dofs + d]).imag();
        }
        else
        {
          if(real)
            entity[d] = abs(values[e * dofs + d]);
          else
            entity[d] = CPhase((complex<double>) values[e * dofs + d]);
        }
      }
    }

    out->SetTuple(e, entity.GetPointer());
  }
}

bool SimOutputVTK::IgnoreImaginaryPart(ResultInfo& info)
{
  switch(info.resultType)
  {
  case MECH_PSEUDO_DENSITY:
  case ELEC_PSEUDO_POLARIZATION:
    return true;
  default:
    return false;
  }
}


//! End multisequence step
void SimOutputVTK::FinishMultiSequenceStep( ) 
{
  // set offset for step value and number to last values
  if( currAnalysis_ == BasePDE::TRANSIENT ||
      currAnalysis_ == BasePDE::STATIC ) {
    stepNumOffset_ = actStep_;
    stepValOffset_ = actStepVal_;
  }

}


//! Finalize the output
void SimOutputVTK::Finalize() 
{
  if(writeSingleFiles_)
  {
    // initialize VTK stuff
    vtkXMLDataElement* pvd = vtkXMLDataElement::New();
    pvd->SetName("VTKFile");
    pvd->SetAttribute("type", "Collection");
    pvd->SetAttribute("version", "0.1");

    vtkXMLDataElement* col = vtkXMLDataElement::New();
    col->SetName("Collection");
    pvd->AddNestedElement(col);

    for(unsigned int i = 0; i < singleFiles.GetSize(); i++)
    {
      SingleVTUFile& item = singleFiles[i];
      
      vtkXMLDataElement* file = vtkXMLDataElement::New();
      file->SetName("DataSet");
      file->SetDoubleAttribute("timestep", item.num_label);
      file->SetAttribute("name", "this is a name");
      file->SetAttribute("group", "");
      file->SetAttribute("part", "0");
      file->SetAttribute("file", item.name.c_str());
      
      col->AddNestedElement(file);
      file->Delete();
    }
    
    string name = "simoutput_vtk/" + fileName_ + ".pvd";
    int res = vtkXMLUtilities::WriteElementToFile(pvd, name.c_str());
    if(res != 1)
    {
      string txt =  "Could not write VTK 'mother file' '" + name + "'";
      Warning(txt.c_str());
    }
    
    col->Delete();
    pvd->Delete();
  }
}

int SimOutputVTK::GetCellIdType(FEType cfsType) 
{
  switch( cfsType ) 
  {
  //  case ET_UNDEF:
  case ET_POINT:    return VTK_VERTEX;
  case ET_LINE2:    return VTK_LINE;
  case ET_LINE3:    return VTK_QUADRATIC_EDGE;
  case ET_TRIA3:    return VTK_TRIANGLE;
  case ET_TRIA6:    return VTK_QUADRATIC_TRIANGLE;
  case ET_QUAD4:    return VTK_QUAD;
  case ET_QUAD8:    return VTK_QUADRATIC_QUAD;
  //case ET_QUAD9:
  case ET_TET4:     return VTK_TETRA;
  case ET_TET10:    return VTK_QUADRATIC_TETRA;
  case ET_HEXA8:    return VTK_HEXAHEDRON;
  case ET_HEXA20:   return VTK_QUADRATIC_HEXAHEDRON;
  //case ET_HEXA27:
  case ET_PYRA5:    return VTK_PYRAMID;
  case ET_PYRA13:   return VTK_QUADRATIC_PYRAMID;
  case ET_WEDGE6:   return VTK_WEDGE;
  case ET_WEDGE15:  return VTK_QUADRATIC_WEDGE;
  default:  EXCEPTION("no VTK expression for CFS element " << cfsType);
  }
}

