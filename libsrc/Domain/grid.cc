#include <iostream.h>
#include <fstream.h>
#include <string>

#include <general_head.hh>
#include <utils_head.hh>

#include <datainout_head.hh>
#include "grid.hh"

namespace CoupledField
{

template<class Dim>
 Grid<Dim> :: Grid(FileType * const aptFileType)
{
  Integer i,ii;
  
#ifdef TRACE
  (*trace) << "entering Grid::Grid" << endl;
#endif
  numlevel = 0;
// ----------------------------- Initialize gh
  Integer dataHelp[1];
  aptFileType->ReadGeneralAnalChoice(dataHelp,FileType::numnode,
                                  FileType::endGAnal);

  gh[0].maxnumnode=dataHelp[0];

  gh[0].ptCoordinate=new Dim[gh[0].maxnumnode];
  aptFileType->ReadCoordinate(gh[0].ptCoordinate, gh[0].maxnumnode);
 
//  if (InfoPrint) { counter1=0;
//                   if (counter1==0)  PrintCoordinate(0, infofile);
//                   counter1++;}  

  Integer data[3];
  aptFileType->ReadGeneralElemChoice(0,data, FileType::numelem,
                   FileType::ielemtyp, FileType::maxnode,
                   FileType::endGElem);

   gh[0].maxnumelem=data[0];

   Integer NumNodeperElem=data[2];

   sizeConnectElem=NumNodeperElem*gh[0].maxnumelem;
   gh[0].Connect=new Integer[sizeConnectElem];   

   Integer sizeFullInfoElem=4*gh[0].maxnumelem;
   gh[0].Info=new Integer[sizeFullInfoElem];  
  
   gh[0].fp=new Integer[gh[0].maxnumelem];

   Integer start=0;
   Integer * help=new Integer[20];
   Integer ihelp=0;   

   aptFileType->ReadElemConnectionGH(gh[0].maxnumelem, gh[0].Connect, NumNodeperElem, 0);
   for (i=0; i<gh[0].maxnumelem; i++) 
   {
      gh[0].Info[start+0]=i;  // global element number
      gh[0].Info[start+1]=0; // element level of last touch
      gh[0].Info[start+2]=ihelp;// start position of connection
      gh[0].Info[start+3]=999; // address of pointer to Element
   
      ihelp+=NumNodeperElem;
      gh[0].fp[i]=start;
      start+=4;
   }

  if (InfoPrint && gh[0].maxnumelem<100)
 {  
  static Integer counter=0;
    if (counter==0)
  { 
    PrintCoordinate(0, infofile);
    for (i=0; i<gh[0].maxnumelem; i++ )
    {
       PrintInfoElem(0, i, infofile);
    }
  }
    counter++;
 }
 
}


/// Deconstructor  
template<class Dim>
 Grid<Dim>:: ~Grid()
{
#ifdef TRACE
  (*trace) << "entering Grid::~Grid" << endl;
#endif

  delete [] gh[0].ptCoordinate;
  delete [] gh[0].Connect;
  delete [] gh[0].Info;
  delete [] gh[0].fp;
  
}

template<class Dim>
void Grid<Dim> :: GetCoordOfNodesElem(const Integer i, const Integer l, 
                                      Dim * ptCoordElem)
{
#ifdef TRACE
  (*trace) << "entering Grid:: GetCoordOfNodesElem" << endl;
#endif

 Integer n,k; 
 Integer stpos=gh[l].Info[gh[l].fp[i]+2]; 

 // Define number nodes per element 7777777777777 To get from Iterator
  if (i==gh[l].maxnumelem-1)
                    n=sizeConnectElem -  gh[l].Info[gh[l].fp[i]+2];
  else
            n=gh[l].Info[gh[l].fp[i+1]+2]- gh[l].Info[gh[l].fp[i]+2];
 
  for (k=0; k < n; k++)
       ptCoordElem[k]=gh[l].ptCoordinate[gh[l].Connect[stpos+k]-1];
}

template<class Dim>
Integer Grid<Dim> :: GetNumNodesPerElem(const Integer i, const Integer l)
{
  Integer n;
  if (i==gh[l].maxnumelem-1)
                    n=sizeConnectElem -  gh[l].Info[gh[l].fp[i]+2];
  else
            n=gh[l].Info[gh[l].fp[i+1]+2]- gh[l].Info[gh[l].fp[i]+2];
  return n;
}

template <class Dim>
inline void Grid<Dim> :: GetConnection(Integer * p, const Integer l, const Integer i,
                                const Integer n)
{
#ifdef TRACE
  (*trace) << "entering Grid:: GetConnection" << endl;
#endif
 
 Integer ii;
 Integer stpos=gh[l].Info[gh[l].fp[i]+2];
 
 Integer fpConnect=gh[l].Info[gh[l].fp[i]+2];

 for (ii=0; ii<n; ii++) 
  p[ii]=gh[l].Connect[fpConnect+ii] ;

}

template<class Dim>
void Grid<Dim> :: PrintCoordinate(const Integer level,ostream * out) const
{
#ifdef TRACE
  (*trace) << "entering Grid::PrintCoordinate" << endl;
#endif

  (*out) << "# coordinates of grid with level " << level << endl; 

  for (Integer i=0; i<gh[level].maxnumnode; i++)
    {
      (*out) << i+1 <<"." ; 
      PrintPoint(gh[level].ptCoordinate[i],out);
      (*out) << endl;  
    }
}

template<class Dim>
void Grid<Dim> :: PrintInfoElem(const Integer l, const Integer i, ostream * out) const
{
#ifdef TRACE
  (*trace) << "entering Grid::PrintInfoElem" << endl;
#endif
 
  (*out) << i+1 <<
  ".\t Level:" <<  gh[l].Info[gh[l].fp[i]+1]
  << "\t ElemType:" <<  gh[l].Info[gh[l].fp[i]+3]
  << endl;

  (*out) << " Connection: ";
  Integer n,ii;
  if (i==gh[l].maxnumelem-1)
     n=sizeConnectElem -
    gh[l].Info[gh[l].fp[i]+2];
  else 
  n=gh[l].Info[gh[l].fp[i+1]+2]-
    gh[l].Info[gh[l].fp[i]+2];
  for (ii=0; ii<n; ii++) (*out) <<
  gh[l].Connect[gh[l].Info[gh[l].fp[i]+2]+ii] << " "; 
  (*out) << endl;
}
} // end namespace
