#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include "outUnverg.hh"

namespace CoupledField
{

template<class Dim>
WriteResultsUnverg<Dim> :: WriteResultsUnverg(const Char * const filename)
:WriteResults<Dim>()
{
#ifdef TRACE
  (*trace) << "entering WriteResultsUnverg :: WriteResultsUnverg" << std::endl;
#endif

  Char * help = new Char[20];
  strcpy(help,filename);
  output=new std::ofstream(strcat(help,".unverg"));

  delete [] help;
}

template<class Dim>
WriteResultsUnverg<Dim> ::~ WriteResultsUnverg()
{
#ifdef TRACE
  (*trace) << "entering WriteResultsUnverg ::~ WriteResultsUnverg" << std::endl;
#endif
  ;
}

template<class Dim>
void WriteResultsUnverg<Dim> :: WriteGrid(const Integer level)
{
 Dataset666(level);
 Dataset781(level);
 Dataset780(level);
}


template<>
void  WriteResultsUnverg<Point2D>::Dataset666(const Integer level)
{
 //
 if (!ptgrid)
    Error("ptgrid is not initialized", __FILE__,__LINE__);

 (*output)<< std::setw(6) << -1 << std::endl << std::setw(6) << -666 << std::endl ;

 Integer maxnumnodes=ptgrid-> GetMaxnumnodes(level);
 Integer maxnumelem=ptgrid-> GetMaxnumElem(level);

 (*output)<< std::setw(10) << 1 << std::setw(10) << 1 << std::setw(10) << 2 << std::endl << std::setw(10) << maxnumnodes << std::setw(10) << maxnumelem << std::endl;

 (*output) << std::setw(6) << -1 << std::endl;
}

template<>
void  WriteResultsUnverg<Point3D>::Dataset666(const Integer level)
{

 (*output)<< std::setw(6) << -1 << std::endl << std::setw(6) << -666 << std::endl ;

 Integer maxnumnodes=ptgrid-> GetMaxnumnodes(level);
 Integer maxnumelem=ptgrid-> GetMaxnumElem(level);

 (*output)<< std::setw(10) << 1 << std::setw(10) << 1 << std::setw(10) << 3 << std::endl << std::setw(10) << maxnumnodes << std::setw(10) << maxnumelem << std::endl;

 (*output) << std::setw(6) << -1 << std::endl;
}

template< >
void  WriteResultsUnverg<Point2D>::Dataset781(const Integer level)
{
  //
 if (!ptgrid)
    Error("ptgrid is not initialized", __FILE__,__LINE__);

 (*output) << std::setw(6) << -1 << std::endl << std::setw(6) << 781 << std::endl;

 Integer maxnumnodes=ptgrid-> GetMaxnumnodes(level);

 (*output).setf(std::ios::scientific);
 (*output).precision(16);

 for (Integer i=0; i<maxnumnodes; i++)
   {
     (*output) << std::setw(10) << i+1 << std::setw(10) << 0 << std::setw(10) << 0 << std::setw(10) << 11 << std::endl;

     (*output).setf(std::ios::uppercase);

     Point2D Point;

     ptgrid->GetCoordinateNode(i,level,Point);

     (*output) << "   " << 0.0 ;

     PrintPoint(Point,output);

     (*output) << std::endl;           
   }
 
 (*output) << std::setw(6) << -1 << std::endl;
}

template< >
void  WriteResultsUnverg<Point3D>::Dataset781(const Integer level)
{
  //
 if (!ptgrid)
    Error("ptgrid is not initialized", __FILE__,__LINE__);

 (*output) << std::setw(6) << -1 << std::endl << std::setw(6) << 781 << std::endl;

 Integer maxnumnodes=ptgrid-> GetMaxnumnodes(level);

 (*output).setf(std::ios::scientific);
 (*output).precision(16);

 for (Integer i=0; i<maxnumnodes; i++)
   {
     (*output) << std::setw(10) << i+1 << std::setw(10) << 0 << std::setw(10) << 0 << std::setw(10) << 11 << std::endl;

     (*output).setf(std::ios::uppercase);

     Point3D Point;

     ptgrid->GetCoordinateNode(i,level,Point);

     PrintPoint(Point,output);

     (*output) << std::endl;
   }

 (*output) << std::setw(6) << -1 << std::endl;
}

template<>
void  WriteResultsUnverg<Point2D>::Dataset780(const Integer level)
{
  //
 if (!ptgrid)
    Error("ptgrid is not initialized", __FILE__,__LINE__);

 (*output) << std::setw(6) << -1 << std::endl << std::setw(6) << 780 << std::endl;

 Integer maxnumnodes=ptgrid-> GetMaxnumnodes(level);
 Integer maxnumelem=ptgrid-> GetMaxnumElem(level);

 Integer * ConnectElem=NULL;
 Integer numnodesPerElem=0;

 for (Integer i=0; i<maxnumelem; i++)
   {
     numnodesPerElem=ptgrid->GetNumNodesPerElem(i,level); 

     (*output) << std::setw(10) << i+1 << std::setw(10);

     switch(numnodesPerElem)
     {
       case 3: (*output) << 91 ; break;
       case 4: (*output) << 94 ; break;
       case 6: (*output) << 92; break;
       case 8: (*output) << 95; break;
       default: Error("Please, put element type according to unverg-format for this number of nodes per element", __FILE__,__LINE__);
     }

     (*output) << std::setw(10) << 2 << std::setw(10) << 2 << std::setw(10) << 1 << std::setw(10) << 1 << std::setw(10) << 1 << std::setw(10) << numnodesPerElem << std::endl;

     ConnectElem=new Integer[numnodesPerElem];      
     ptgrid->GetConnection(ConnectElem, level, i, numnodesPerElem);

     for (Integer ii=0; ii < numnodesPerElem; ii++) 
       { 
	 (*output).width(10);
	 (*output) << ConnectElem[ii];
       }

     delete [] ConnectElem;
 
     (*output) << std::endl;
   }
 (*output) << std::setw(6) << -1 << std::endl;
}

template<>
void  WriteResultsUnverg<Point3D>::Dataset780(const Integer level)
{
  //
 if (!ptgrid)
    Error("ptgrid is not initialized", __FILE__,__LINE__);

 (*output) << std::setw(6) << -1 << std::endl << std::setw(6) << 780 << std::endl;

 Integer maxnumnodes=ptgrid-> GetMaxnumnodes(level);
 Integer maxnumelem=ptgrid-> GetMaxnumElem(level);

 Integer * ConnectElem=NULL;
 Integer numnodesPerElem=0;

 for (Integer i=0; i<maxnumelem; i++)
   {
     numnodesPerElem=ptgrid->GetNumNodesPerElem(i,level);

     (*output) << std::setw(10) << i+1 << std::setw(10);

     switch(numnodesPerElem)
     {
       case 4: (*output) << 111 ; break;
       case 6: (*output) << 112; break;
       case 8: (*output) << 115; break;
       case 15: (*output) << 113; break;
       case 20: (*output) << 116; break;
       default: Error("Please, put element type according to unverg-format for this number of nodes per element", __FILE__,__LINE__);
     }

     (*output) << std::setw(10) << 11 << std::setw(10) << 1 << std::setw(10) << 1 << std::setw(10) << 1 << std::setw(10) << 1 << std::setw(10) << numnodesPerElem << std::endl;

     ConnectElem=new Integer[numnodesPerElem];
     ptgrid->GetConnection(ConnectElem, level, i, numnodesPerElem);

     for (Integer ii=0; ii < numnodesPerElem; ii++)
       {
         (*output).width(10);
         (*output) << ConnectElem[ii];
       }

     delete [] ConnectElem;

     (*output) << std::endl;
   }
 (*output) << std::setw(6) << -1 << std::endl;
}

template<class Dim>
void  WriteResultsUnverg<Dim>::Dataset55(const std::string & title, const Vector<Double> & x, const Integer step, const Double time)
{
  //
 if (!ptgrid)
    Error("ptgrid is not initialized", __FILE__,__LINE__);

 (*output) << std::setw(6) << -1 << std::endl << std::setw(6) << 55 << std::endl;

 (*output).setf(std::ios::scientific);
 (*output).precision(6);
 (*output).setf(std::ios::uppercase);

 (*output) << " " << title << ", step" << std::setw(6) << step <<
             " time   " << time << std::endl;  
 (*output) << std::endl << std::endl << std::endl << std::endl;
 (*output) << std::setw(10) << 1 << std::setw(10) << 4 << std::setw(10) << 1 << std::setw(10) << 0
           << std::setw(10) << 2 << std::setw(10) << 1 << std::endl;
 (*output) << std::setw(10) << 2 << std::setw(10) << 1 << std::setw(10) << 1 << std::setw(10) <<
              step << std::endl;
 (*output) << " " << time << std::endl;       

 Integer i,n;
 n=x.size();  
 for (i=0; i<n; i++)
   (*output) << std::setw(10) << i+1 << std::endl << " " << x[i] << std::endl;
     
 (*output) << std::setw(6) << -1 << std::endl;
}  

template<class Dim>
void  WriteResultsUnverg<Dim>::Dataset56()
{
 (*output) << std::setw(6) << -1 << std::endl << std::setw(6) << 56 << std::endl;
 
 (*output) << std::setw(6) << -1;
}

template<class Dim>
void  WriteResultsUnverg<Dim>::Init(Grid<Dim> * aptgrid)
{
 ptgrid=aptgrid;
}

template<class Dim>
void  WriteResultsUnverg<Dim>::WriteSolution(const Vector<Double> & sol, const Integer step, const Double time)
{
 Dataset55(" fluid potential", sol, step+1, time);

/*  switch(typesol)
 {
  case fluid:
   Dataset55(" fluid potential", sol, step+1, time);
   break;
  case temperature:
   Dataset55(" temperature", sol, step+1, time);
   break;
  default:
   Error("Unknown type of results", __FILE__,__LINE__);
 }
*/

}

template<class Dim>
void  WriteResultsUnverg<Dim>::WriteFirstDerSolution(const Vector<Double> & sol, const Integer step, const Double time)
{
 Dataset55(" fluid potential, 1st deriv., ", sol, step+1,time);
}

template<class Dim>
void  WriteResultsUnverg<Dim>::WriteSecondDerSolution(const Vector<Double> & sol, const Integer step, const Double time)
{
 Dataset55(" fluid potential, 2nd deriv., ", sol, step+1, time);
}

} // end of namespace
