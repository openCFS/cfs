#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

#include "outGMV.hh"

namespace CoupledField
{

WriteResultsGMV :: WriteResultsGMV(const Char * const filename)
: WriteResults(filename)
{
#ifdef TRACE
  (*trace) << "entering WriteResultsGMV :: WriteResultsGMV" << std::endl;
#endif

 namedir_=new Char[30];

 Char S[50];
 strcpy(namedir_,filename);
 strcat(namedir_,"_gmv");
 sprintf(S,"mkdir -p %s",namedir_);

 system(S);

 OpenFile(0);

 if (!output)
   Error(" File for output results in .gmv-format ", __FILE__, __LINE__);

 ptgrid=NULL;
 currstep_=0;

}

WriteResultsGMV ::~WriteResultsGMV()
{
#ifdef TRACE
  (*trace) << "entering WriteResultsGMV::~ WriteResultsGMV" << std::endl;
#endif

 // write keyword
 (*output) << "endgmv " << std::endl;

 delete [] namedir_;
}

void WriteResultsGMV :: WriteHeader()
{
 (*output) << "gmvinput" << " ascii" << std::endl;
}

void WriteResultsGMV :: WriteNodes(const Integer alevel)
{
  Integer level=alevel;

  // write keyword
 (*output) << "nodev ";

 //get and write number of nodes on the level
 Integer numnodes=ptgrid->GetMaxnumnodes(level);
 (*output) << numnodes << std::endl;

 Integer dim=ptgrid->GetDim();

 //get and write coodinates of nodes
 Integer i;

 if (dim==2)
{
 Point2D point;

 // write x,y,z-coordinate
  for (i=0; i<numnodes; i++)
    {
      ptgrid->GetCoordinateNode(i,level,point);
      (*output) << " " << point.x << " " << point.y << " " << 0 << std::endl;
    }
}
  else
{
 Point3D point;

 //
 if (!ptgrid)
    Error("ptgrid is not initialized", __FILE__,__LINE__);

 // write x,y,z-coordinate
  for (i=0; i<numnodes; i++)
    {
      ptgrid->GetCoordinateNode(i,level,point);
      (*output) << " " << point.x << " " << point.y << " " << point.z << std::endl;
    }
}
}

void WriteResultsGMV:: WriteCells(const Integer alevel) 
{
  Integer level=alevel;

// write keyword
 (*output) << "cells ";

 //
 if (!ptgrid)
    Error("ptgrid is not initialized", __FILE__,__LINE__);

// read information about number of elements 
  Integer numelem; 
  numelem=ptgrid->GetMaxnumElem(level);

  (*output) << numelem << std::endl;

  Vector<Integer> connect;

  Integer dim=ptgrid->GetDim();

  Integer i;
  for (i=0; i<numelem; i++)
   {

     ptgrid->GetConnection(connect, i, level);

   if (dim==2)
{
     switch (connect.size())
      {
        case 3: 
                (*output) << "tri 3" << std::endl;
                break;
        case 4:
                (*output) << "quad 4" << std::endl;
                 break;
        default:
          Error("This type of element is not implemented", __FILE__, __LINE__);
      }
}
   else
{
     switch (connect.size())
      {
        case 4:
                (*output) << "tet 4" << std::endl;
                break;
        case 8:
                (*output) << "hex 8" << std::endl;
                 break;
        default:
          Error("This type of element is not implemented", __FILE__, __LINE__);
      }
}

     Integer j;
     for (j=0; j< connect.size(); j++)
       (*output) << " " << connect[j] ;

     (*output) << std::endl;
   }
}

void WriteResultsGMV::WriteVariable(const Vector<Double> var, const std::string name, const Integer type)
{
  (*output) << "variable" << std::endl;

  (*output) << name << " " << type << std::endl;

  Integer i;
  for (i=0; i<var.size(); i++)
    (*output) << var[i] << " ";

    (*output) << std::endl;

  (*output) << "endvars" << std::endl;
}

void WriteResultsGMV::WriteGrid(const Integer level)
{

 WriteHeader();
 WriteNodes(level);
 WriteCells(level);
}

void WriteResultsGMV::WriteSolution(const Vector<Double> & sol, const Integer step, const Double time, const std::string title)
{
#ifdef TRACE
 (*trace) << " entering WriteResultsGMV::WriteSolution " << std::endl;
#endif

  Integer i;
  if (NeedHistory_)
       for (i=0; i<nodeshist_.size(); i++)
       { if (sol.size()<=nodeshist_[i])
         Error("Please, check history-nodes in config-file.",__FILE__,__LINE__);
         if (lastsavetime[i] != time )
          AddInHistory(time,sol[nodeshist_[i]],i);
       }

  Integer type=1; // 0 - for cell 
                  // 1 - for node
                  // 2 - for face data

     if (step!=currstep_)
{
   (*output) << "endgmv " << std::endl;
   delete output;

   OpenFile(step);

   WriteGrid(ptgrid->GetLastLevel());
}

  WriteVariable(sol,title,type);
  (*output) << "probtime " << time << std::endl;

  currstep_=step;
}

void WriteResultsGMV::OpenFile(const Integer num)
{
   Char * name=new Char[30];
   Char * aux=new Char[2];
   sprintf(aux,"%i",num);

   strcpy(name,namedir_);
   strcat(name,"/");
   strcat(name,namefile_);
   if (num/10 < 1) strcat(name,".gmv00");
     else if (num/100 < 1) strcat(name,".gmv0");
       else strcat(name,".gmv");

   strcat(name,aux);

   output=new std::ofstream(name);

   delete [] name;
   delete [] aux;
}

void WriteResultsGMV::Init(Grid * aptgrid)
{
ptgrid=aptgrid;
}

} // end of namespace
