#include <iostream>
#include <fstream>

#include "grid_cfs.hh"
#include "elements_header.hh"

namespace CoupledField
{

template<class Dim>
 GridCFS<Dim> :: GridCFS(FileType * const aptFileType)
{
#ifdef TRACE
  (*trace) << "entering GridCFS::GridCFS" << std::endl;
#endif

  numlevel = 0;
  InFile = aptFileType;
}

template<class Dim>
void GridCFS<Dim> :: Read()
{
#ifdef TRACE
  (*trace) << "entering GridCFS::GridCFSRead" << std::endl;
#endif

// ----------------------------- Initialize gh
  Integer i;
  Integer dataHelp[1];
  InFile->ReadGeneralAnalChoice(dataHelp,FileType::numnode,
                                  FileType::endGAnal);

  gh[0].maxnumnode=dataHelp[0];

  gh[0].ptCoordinate=new Dim[gh[0].maxnumnode];
  InFile->ReadCoordinate(gh[0].ptCoordinate, gh[0].maxnumnode);
 
  Integer data[3];
  InFile->ReadGeneralElemChoice(0,data, FileType::numelem,
                   FileType::ielemtyp, FileType::maxnode,
                   FileType::endGElem);

   gh[0].maxnumelem=data[0];

//#################### V etom meste budet check na 3 tochki
   Integer NumNodeperElem=data[2];
   ptArrayElem_=new BaseElem*[gh[0].maxnumelem+1];   
   
   ptQ_=new Quad1(GaussOrder2);
   ptTr_=new Triangle1(GaussOrder3);

   for (i=0; i<gh[0].maxnumelem; i++) 
   { 
       switch(NumNodeperElem)
     {
       case 3:
          ptArrayElem_[i]=ptTr_;
          break;

        case 4:
          ptArrayElem_[i]=ptQ_;
          break;

        default:
           Error("Number of nodes per element is strange",__FILE__,__LINE__);
     }
   }
   ptArrayElem_[gh[0].maxnumelem]=NULL;

   sizeConnectElem=NumNodeperElem*gh[0].maxnumelem;
   gh[0].Connect=new Integer[sizeConnectElem];   

   Integer sizeFullInfoElem=4*gh[0].maxnumelem;
   gh[0].Info=new Integer[sizeFullInfoElem];  
  
   gh[0].fp=new Integer[gh[0].maxnumelem];

   Integer start=0;
   Integer * help=new Integer[20];
   Integer ihelp=0;   

   InFile->ReadElemConnectionGH(gh[0].maxnumelem, gh[0].Connect, NumNodeperElem, 0);
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

}


/// Deconstructor  
template<class Dim>
 GridCFS<Dim>:: ~GridCFS()
{
#ifdef TRACE
  (*trace) << "entering GridCFS::~GridCFS" << std::endl;
#endif
  delete [] gh[0].ptCoordinate;
  delete [] gh[0].Connect;
  delete [] gh[0].Info;
  delete [] gh[0].fp;

  if (ptQ_) delete ptQ_;
  if (ptTr_) delete ptTr_;

   if (ptArrayElem_) delete [] ptArrayElem_;
  
}

template<class Dim>
void GridCFS<Dim> :: GetCoordOfNodesElem(const Integer i, const Integer l, const Integer numnodes,  Dim * ptCoordElem)   
{
#ifdef TRACE
  (*trace) << "entering GridCFS:: GetCoordOfNodesElem" << std::endl;
#endif

 Integer k; 
 Integer stpos=gh[l].Info[gh[l].fp[i]+2]; 

 // Define number nodes per element 7777777777777 To get from Iterator
//  if (i==gh[l].maxnumelem-1)
//                    n=sizeConnectElem -  gh[l].Info[gh[l].fp[i]+2];
//  else
//            n=gh[l].Info[gh[l].fp[i+1]+2]- gh[l].Info[gh[l].fp[i]+2];
 
  for (k=0; k < numnodes; k++)
       ptCoordElem[k]=gh[l].ptCoordinate[gh[l].Connect[stpos+k]-1];
}

template<class Dim>
Integer GridCFS<Dim> :: GetNumNodesPerElem(const Integer i, const Integer l)
{
  Integer n;
  if (i==gh[l].maxnumelem-1)
                    n=sizeConnectElem -  gh[l].Info[gh[l].fp[i]+2];
  else
            n=gh[l].Info[gh[l].fp[i+1]+2]- gh[l].Info[gh[l].fp[i]+2];
  return n;
}

template <class Dim>
inline void GridCFS<Dim> :: GetConnection(Integer * p, const Integer l, const Integer i,
                                const Integer n)
{
 
 Integer ii;
// Integer stpos=gh[l].Info[gh[l].fp[i]+2];
 
 Integer fpConnect=gh[l].Info[gh[l].fp[i]+2];

 for (ii=0; ii<n; ii++) 
  p[ii]=gh[l].Connect[fpConnect+ii] ;

}

template<class Dim>
void GridCFS<Dim> :: PrintCoordinate(const Integer level,std::ostream * out) const
{
#ifdef TRACE
  (*trace) << "entering GridCFS::PrintCoordinate" << std::endl;
#endif

  (*out) << "# coordinates of grid with level " << level << std::endl; 

  for (Integer i=0; i<gh[level].maxnumnode; i++)
    {
      (*out) << i+1 <<"." ; 
      PrintPoint(gh[level].ptCoordinate[i],out);
      (*out) << std::endl;  
    }
}

template<class Dim>
void GridCFS<Dim> :: PrintInfoElem(const Integer l, const Integer i, std::ostream * out) const
{
#ifdef TRACE
  (*trace) << "entering GridCFS::PrintInfoElem" << std::endl;
#endif
 
  (*out) << i+1 <<
  ".\t Level:" <<  gh[l].Info[gh[l].fp[i]+1]
  << "\t ElemType:" <<  gh[l].Info[gh[l].fp[i]+3]
  << std::endl;

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
  (*out) << std::endl;
}

template<class Dim>
void GridCFS<Dim> :: GetCoordinateNode(const Integer inode, const Integer numlevel, Dim & rfPoint) 
{
#ifdef TRACE
  (*trace) << "entering GridCFS::GetCoordinateNode" << std::endl;
#endif
  rfPoint=gh[numlevel].ptCoordinate[inode];
}

} // end namespace
