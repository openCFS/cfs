#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

#include "outGMV.hh"

namespace CoupledField
{

template<class Dim>
WriteResultsGMV<Dim> :: WriteResultsGMV(const Char * const filename)
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

template<class Dim>
WriteResultsGMV<Dim> ::~WriteResultsGMV()
{
#ifdef TRACE
  (*trace) << "entering WriteResultsGMV::~ WriteResultsGMV" << std::endl;
#endif

 // write keyword
 (*output) << "endgmv " << std::endl;

 delete output; 
 delete [] namedir_;
}

template<class Dim>
void WriteResultsGMV<Dim> :: WriteHeader()
{
 (*output) << "gmvinput" << " ascii" << std::endl;
}

template<>
void WriteResultsGMV<Point2D> :: WriteNodes(const Integer alevel)
{
  Integer level=alevel;

/*
 // write keyword
 (*output) << "nodes ";

 //get and write number of nodes on the level 
 Integer numnodes=ptgrid->GetMaxnumnodes(level); 
 (*output) << numnodes << std::endl;

 //get and write coodinates of nodes
 Integer i;
 Point2D point;

 // 
 if (!ptgrid)
    Error("ptgrid is not initialized", __FILE__,__LINE__);

 // write x-coordinate
  for (i=0; i<numnodes; i++)
    { 
      ptgrid->GetCoordinateNode(i,level,point);
      (*output) << point.x << " ";
    }
  (*output) << std::endl;

 // write y-coordinate
  for (i=0; i<numnodes; i++)
    {
      ptgrid->GetCoordinateNode(i,level,point);
      (*output) << point.y << " ";
    }
   (*output) << std::endl;

 // write z-coordinate
   for (i=0; i<numnodes; i++)
    {
      (*output) << 0 << " ";
    }
   (*output) << std::endl;
*/

   // alternative
    // write keyword
 (*output) << "nodev ";

 //get and write number of nodes on the level
 Integer numnodes=ptgrid->GetMaxnumnodes(level);
 (*output) << numnodes << std::endl;

 //get and write coodinates of nodes
 Integer i;
 Point2D point;

 // write x,y,z-coordinate
  for (i=0; i<numnodes; i++)
    {
      ptgrid->GetCoordinateNode(i,level,point);
      (*output) << " " << point.x << " " << point.y << " " << 0 << std::endl;
    }
}

template<>
void WriteResultsGMV<Point3D> :: WriteNodes(const Integer alevel)
{
  Integer level=alevel;

/*
 // write keyword
 (*output) << "nodes ";

 //get and write number of nodes on the level
 Integer numnodes=ptgrid->GetMaxnumnodes(level);
 (*output) << numnodes << std::endl;

 //get and write coodinates of nodes
 Integer i;
 Point2D point;

 // write x-coordinate
  for (i=0; i<numnodes; i++)
    {
      ptgrid->GetCoordinateNode(i,level,point);
      (*output) << point.x << " ";
    }
  (*output) << std::endl;

 // write y-coordinate
  for (i=0; i<numnodes; i++)
    {
      ptgrid->GetCoordinateNode(i,level,point);
      (*output) << point.y << " ";
    }
   (*output) << std::endl;

 // write z-coordinate
   for (i=0; i<numnodes; i++)
    {
      (*output) << 0 << " ";
    }
   (*output) << std::endl;
*/

   // alternative
    // write keyword
 (*output) << "nodev ";

 //get and write number of nodes on the level
 Integer numnodes=ptgrid->GetMaxnumnodes(level);
 (*output) << numnodes << std::endl;

 //get and write coodinates of nodes
 Integer i;
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

template<>
void WriteResultsGMV<Point2D>:: WriteCells(const Integer alevel) 
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

  Integer i;
  for (i=0; i<numelem; i++)
   {

     ptgrid->GetConnection(connect, i, level);

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

     Integer j;
     for (j=0; j< connect.size(); j++)
       (*output) << " " << connect[j] ;

     (*output) << std::endl;
   }
}

template<>
void WriteResultsGMV<Point3D>:: WriteCells(const Integer alevel)
{
  Integer level=alevel;

// write keyword
 (*output) << "cells ";

 //
 if (!ptgrid)
    Error("ptgrid is not initialized", __FILE__,__LINE__);

// read information about number of elements and number of nodes per element
  Integer numelem;

  numelem=ptgrid->GetMaxnumElem(level);
  (*output) << numelem << std::endl;

  Vector<Integer> connect;

  Integer i;
  for (i=0; i<numelem; i++)
   {
     ptgrid->GetConnection(connect, i, level);

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

     Integer j;
     for (j=0; j<connect.size(); j++)
       (*output) << " " << connect[j] ;

     (*output) << std::endl;
   }
}

template<class Dim>
void WriteResultsGMV<Dim>::WriteVariable(const Vector<Double> var, const std::string name, const Integer type)
{
  (*output) << "variable" << std::endl;

  (*output) << name << " " << type << std::endl;

  Integer i;
  for (i=0; i<var.size(); i++)
    (*output) << var[i] << " ";

    (*output) << std::endl;

  (*output) << "endvars" << std::endl;
}

template<class Dim>
void WriteResultsGMV<Dim>::WriteGrid(const Integer level)
{

 WriteHeader();
 WriteNodes(level);
 WriteCells(level);
}

template<class Dim>
void WriteResultsGMV<Dim>::WriteSolution(const Vector<Double> & sol, const Integer step, const Double time, const std::string title)
{
#ifdef TRACE
 (*trace) << " entering WriteResultsGMV<Dim>::WriteSolution " << std::endl;
#endif

  if (sol.size()<=history_node_) Error("Please, check history-nodes in config-file.",__FILE__,__LINE__);

  Integer i;
  if (NeedHistory_)
       for (i=0; i<nodeshist_.size(); i++)
          AddInHistory(time,sol[nodeshist_[i]],i);

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

template<class Dim>
void WriteResultsGMV<Dim>::OpenFile(const Integer num)
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

template<class Dim>
void WriteResultsGMV<Dim>::Init(Grid * aptgrid)
{
ptgrid=aptgrid;
}

} // end of namespace
