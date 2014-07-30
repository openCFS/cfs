/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkOFFReader.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkOBJReader - read Wavefront .obj files
// .SECTION Description
// vtkOBJReader is a source object that reads Wavefront .obj
// files. The output of this source object is polygonal data.
// .SECTION See Also
// 

#ifndef __vtkOFFReader_h
#define __vtkOFFReader_h

#include "vtkIOOFFReaderModule.h" // For export macro
#include "vtkPolyDataAlgorithm.h"

class VTKIOOFFREADER_EXPORT vtkOFFReader : public vtkPolyDataAlgorithm 
{
public:
  static vtkOFFReader *New();
  vtkTypeMacro(vtkOFFReader,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Specify file name of .off file.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

protected:
  vtkOFFReader();
  ~vtkOFFReader();
  
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  char *FileName;
private:
  vtkOFFReader(const vtkOFFReader&);  // Not implemented.
  void operator=(const vtkOFFReader&);  // Not implemented.
};

#endif
