/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkOFFReader.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkOFFReader.h"

#include "vtkCellData.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkUnsignedIntArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkTriangle.h"
#include "vtkQuad.h"
#include "vtkPolygon.h"
#include "vtkLookupTable.h"
// #include "Rendering/Core/vtkColorTransferFunction.h"

#include <sstream>
#include <fstream>

vtkStandardNewMacro(vtkOFFReader);

// Description:
// Instantiate object with NULL filename.
vtkOFFReader::vtkOFFReader()
{
  this->FileName = NULL;

  this->SetNumberOfInputPorts(0);
}

vtkOFFReader::~vtkOFFReader()
{
  if (this->FileName)
  {
    delete [] this->FileName;
    this->FileName = NULL;
  }
}

/*--------------------------------------------------------

This is only partial support for the OFF format, which is 
quite complicated. To find a full specification,
search the net for "OFF format", eg.:

http://people.sc.fsu.edu/~burkardt/data/off/off.html

We support only vertices and faces composed of 3 vertices.

The above is not true any more. We now support also 
arbitrary flat polygons with vertex and face colors.
But still just ASCII .off files are supported, but this
should be sufficient for supporting the output meshes from
CGAL.

---------------------------------------------------------*/

// a replacement for isspace()
int is_whitespace(char c)
{
  if ( c==' ' || c=='\t' || c=='\n' || c=='\r' || c=='\v' || c=='\f')
    return 1;
  else
    return 0;
}

int vtkOFFReader::RequestData(
                              vtkInformation *vtkNotUsed(request),
    vtkInformationVector **vtkNotUsed(inputVector),
                                      vtkInformationVector *outputVector)
{
  // get the info object
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the ouptut
  vtkPolyData *output = vtkPolyData::SafeDownCast(
      outInfo->Get(vtkDataObject::DATA_OBJECT()));

  if (!this->FileName) 
  {
    vtkErrorMacro(<< "A FileName must be specified.");
    return 0;
  }
    
  std::ifstream in(this->FileName);
      
  if (!in) 
  {
    vtkErrorMacro(<< "File " << this->FileName << " not found");
    return 0;
  }

  vtkDebugMacro(<<"Reading file");

  // intialise some structures to store the file contents in
  vtkPoints *points = vtkPoints::New(); 
  vtkCellArray *polys = vtkCellArray::New();
  
  // -- work through the file line by line, assigning into the above structures as appropriate --

    float xyz[3];

  std::string line;

    getline(in, line); //throw the 1st line away, it always says "OFF"
  
  //Read until we get to the first line that is not a comment or a blank line.
  //This line second states the number of vertices, faces, and edges. 
  while(getline(in, line))
  {
    //std::cout << line << std::endl;
    if(line.size() == 0)
      continue;
    if(line[0] == '#')
      continue;
    
    //if we get to here, this is the info line
    break;
    
  }
  //std::cout << "Found info line." << std::endl;
  //At this point, the line were are interested in is stored in 'line'
  //We are only interested in vertices and faces.
  std::stringstream ss;
  ss << line;
  unsigned int NumberOfVertices, NumberOfFaces;
  ss >> NumberOfVertices >> NumberOfFaces;
  //std::cout << "Vertices: " << NumberOfVertices << " Faces: " << NumberOfFaces << std::endl;
  //Now we can start to read the vertices
  unsigned int VertexCounter = 0;

  vtkFloatArray *vscalars = vtkFloatArray::New();
  vscalars->SetNumberOfComponents(1);
  vscalars->SetName("VertexColor");
  // output->GetPointData()->AddArray(vscalars);
  vtkLookupTable* vlut = vtkLookupTable::New();
  vscalars->SetLookupTable(vlut);

  vtkUnsignedCharArray *vcolors = vtkUnsignedCharArray::New();
  vcolors->SetNumberOfComponents(4);
  vcolors->SetName("VertexColor");
  output->GetPointData()->SetScalars(vcolors);

  vlut->SetTable(vcolors);

  while(getline(in, line) && VertexCounter < NumberOfVertices)
  {
    if(line.size() == 0) //skip blank lines (they should only occur before the vertices start)
      continue;

    std::stringstream ssVertex;
    ssVertex << line;
    float x,y,z;
    if( !(ssVertex >> x) ) continue;
    ssVertex >> y >> z;
    points->InsertNextPoint(x, y, z);

    double r=1.0;
    double g=1.0;
    double b=1.0;
    double a=1.0;
    if( !(ssVertex >> r) ) r=1.0;
    if( !(ssVertex >> g) ) g=1.0;
    if( !(ssVertex >> b) ) b=1.0;
    if( !(ssVertex >> a) ) a=1.0;

    vscalars->InsertNextValue(VertexCounter);
    vcolors->InsertNextValue(static_cast<unsigned char>(255.0*r));
    vcolors->InsertNextValue(static_cast<unsigned char>(255.0*g));
    vcolors->InsertNextValue(static_cast<unsigned char>(255.0*b));
    vcolors->InsertNextValue(static_cast<unsigned char>(255.0*a));

    //std::cout << "adding vertex: " << x << " " << y << " " << z << std::endl;
    VertexCounter++;
  } // (end of vertex while loop)

  unsigned int FaceCounter = 0;

  vtkFloatArray *fscalars = vtkFloatArray::New();
  fscalars->SetNumberOfComponents(1);
  fscalars->SetName("FaceColor");
  output->GetCellData()->AddArray(fscalars);
  vtkLookupTable* flut = vtkLookupTable::New();
  fscalars->SetLookupTable(flut);

  vtkUnsignedCharArray *fcolors = vtkUnsignedCharArray::New();
  fcolors->SetNumberOfComponents(4);
  fcolors->SetName("FaceColor");
  //  output->GetCellData()->SetScalars(fcolors);

  flut->SetTable(fcolors);
  
  //read faces
  do
  {
    std::stringstream ssFace;
    ssFace << line;
    unsigned int NumFaceVerts;
    unsigned int Vertices[1024];
    if( !(ssFace >> NumFaceVerts) ) continue;
    double coords[3];
    double* colors = NULL;

    // Read node ids of face
    for(unsigned int i=0; i<NumFaceVerts; i++){
        ssFace >> Vertices[i];
    }

    // Read color of face
    double r=1.0;
    double g=1.0;
    double b=1.0;
    double a=1.0;
    if( !(ssFace >> r) ) r=1.0;
    if( !(ssFace >> g) ) g=1.0;
    if( !(ssFace >> b) ) b=1.0;
    if( !(ssFace >> a) ) a=1.0;
    fcolors->InsertNextValue(static_cast<unsigned char>(255.0*r));
    fcolors->InsertNextValue(static_cast<unsigned char>(255.0*g));
    fcolors->InsertNextValue(static_cast<unsigned char>(255.0*b));
    fcolors->InsertNextValue(static_cast<unsigned char>(255.0*a));

    //    std::cout << "Face " << FaceCounter << " " << Vertices[0] << " " << Vertices[1] << " " << Vertices[2] << " "<< Vertices[3] << std::endl;
    //    std::cout << "Face " << FaceCounter << " " << r << " " << g << " " << b << " "<< a << std::endl;

    // Create faces in vtkPolyData.
    switch(NumFaceVerts) 
    {
    case 3:
      {
        vtkTriangle* triangle = vtkTriangle::New();
        
        triangle->GetPointIds()->SetId(0, Vertices[0]);
        triangle->GetPointIds()->SetId(1, Vertices[1]);
        triangle->GetPointIds()->SetId(2, Vertices[2]);
        
        polys->InsertNextCell(triangle);
        
        triangle->Delete();
      }
      break;
    case 4:
      {
        vtkQuad* quad = vtkQuad::New();
      
        quad->GetPointIds()->SetId(0, Vertices[0]);
        quad->GetPointIds()->SetId(1, Vertices[1]);
        quad->GetPointIds()->SetId(2, Vertices[2]);
        quad->GetPointIds()->SetId(3, Vertices[3]);
        
        polys->InsertNextCell(quad);
        
        quad->Delete();
      }
      break;
    default:
      {
        vtkPolygon* poly = vtkPolygon::New();
      
        poly->GetPointIds()->SetNumberOfIds(NumFaceVerts);
        for(unsigned i=0; i<NumFaceVerts; i++) 
        {
          poly->GetPointIds()->SetId(i, Vertices[i]);
        }
        
        polys->InsertNextCell(poly);
        
        poly->Delete();
      }
      break;      
    }
    
    fscalars->InsertNextValue((1.0*FaceCounter)/(NumberOfFaces-1));

    FaceCounter++;

    //    std::cout << "FaceCounter: " << FaceCounter << std::endl;
    
  }while(getline(in, line) && FaceCounter < NumberOfFaces);

  //  std::cout << "JIPI" << std::endl;

  for(unsigned i=0; i<VertexCounter; i++)
    vscalars->SetTuple1(i, vscalars->GetTuple1(i) / (VertexCounter-1));
  
  vlut->SetNumberOfTableValues(VertexCounter);
  vlut->SetRange(0, 1);
  //vlut->SetValueRange(0, VertexCounter-1);
  //  vlut->Build();
  vlut->SetRampToLinear();
  //output->GetPointData()->SetActiveScalars("VertexColor");

  flut->SetNumberOfTableValues(NumberOfFaces);
  flut->SetRange(0, 1);
  //flut->SetValueRange(0, NumberOfFaces-1);

  //  std::cout << "Index for 3.0: " << flut->GetIndex(3.0) << std::endl;
  //  flut->Build();

  
  // we have finished with the file
  in.close(); 

  output->SetPoints(points);
  output->SetPolys(polys);

  points->Delete();
  polys->Delete();
  
  return 1;
}

void vtkOFFReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "File Name: " 
      << (this->FileName ? this->FileName : "(none)") << "\n";

}
