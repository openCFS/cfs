#include <iostream>
#include <fstream>
#include <string>

#include "outGMV.hh"

namespace CoupledField
{

template<class Dim>
WriteResultsGMV<Dim> :: WriteResultsGMV(const Char * filename, Grid<Dim> * aptgrid)
: WriteResults<Dim>()
{
#ifdef TRACE
  (*trace) << "entering WriteResultsGMV :: WriteResultsGMV" << std::endl;
#endif

 Char * help=new Char[20];
 strcpy(help,filename);
 output=new std::ofstream(strcat(help,".gmv"));

 if (!output)
   Error(" File for output results in .gmv-format ", __FILE__, __LINE__); 

 delete [] help;

 ptgrid=aptgrid; 

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

// read information about number of elements and number of nodes per element
  Integer numelem, elemsize; 

  numelem=ptgrid->GetMaxnumElem(level);
  elemsize=ptgrid->GetNumNodesPerElem(0,level);  //// !!! Attention !!!

  (*output) << numelem << std::endl;

  Integer * connect;

  Integer i;
  for (i=0; i<numelem; i++)
   {

     elemsize=ptgrid->GetNumNodesPerElem(i,level);
     connect=new Integer[elemsize];

     ptgrid->GetConnection(connect, level, i, elemsize);

     switch (elemsize)
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
     for (j=0; j<elemsize; j++)
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

// read information about number of elements and number of nodes per element
  Integer numelem, elemsize;

  numelem=ptgrid->GetMaxnumElem(level);
  elemsize=ptgrid->GetNumNodesPerElem(0,level);  //// !!! Attention !!!

  (*output) << numelem << std::endl;

  Integer * connect;

  Integer i;
  for (i=0; i<numelem; i++)
   {

     elemsize=ptgrid->GetNumNodesPerElem(i,level);
     connect=new Integer[elemsize];

     ptgrid->GetConnection(connect, level, i, elemsize);

     switch (elemsize)
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
     for (j=0; j<elemsize; j++)
       (*output) << " " << connect[j] ;

     (*output) << std::endl;
   }
}

template<class Dim>
void WriteResultsGMV<Dim>::WriteGrid(const Integer level)
{

 WriteHeader();
 WriteNodes(level);
 WriteCells(level);
}

template<class Dim>
void WriteResultsGMV<Dim>::WriteSolution(const Vector<Double> & sol, const Integer step, const Integer time)
{
;
}

template<class Dim>
void WriteResultsGMV<Dim>::WriteFirstDerSolution(const Vector<Double> & sol, const Integer step, const Integer time)
{
;
}

template<class Dim>
void WriteResultsGMV<Dim>::WriteSecondDerSolution(const Vector<Double> & sol, const Integer step, const Integer time)
{
;
}

template<class Dim>
void WriteResultsGMV<Dim>::Init(Grid<Dim> * aptgrid)
{
ptgrid=aptgrid;
}

} // end of namespace
