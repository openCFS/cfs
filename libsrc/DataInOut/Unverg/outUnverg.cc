#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include "outUnverg.hh"

namespace CoupledField
{

template<class Dim>
OutResultUnverg<Dim> :: OutResultUnverg(const Char * const filename)
{
#ifdef TRACE
  (*trace) << "entering OutResultUnverg :: OutResultUnverg" << std::endl;
#endif

  Char * help = new Char[20];
  strcpy(help,filename);
  output=new std::ofstream(strcat(help,".unverg"));

}

template<class Dim>
OutResultUnverg<Dim> ::~ OutResultUnverg()
{
#ifdef TRACE
  (*trace) << "entering OutResultUnverg ::~ OutResultUnverg" << std::endl;
#endif

  delete output;
}

template<class Dim>
void OutResultUnverg<Dim> :: Create(Grid<Dim> * ptgrid, const Integer level)
{
 Dataset666(ptgrid,level);
 Dataset781(ptgrid,level);
 Dataset780(ptgrid,level);
 ;
}


template<class Dim>
void  OutResultUnverg<Dim>::Dataset666(Grid<Dim> * ptgrid, const Integer level)
{

 (*output)<< std::setw(6) << -1 << std::endl << std::setw(6) << -666 << std::endl ;

 Integer maxnumnodes=ptgrid-> GetMaxnumnodes(level);
 Integer maxnumelem=ptgrid-> GetMaxnumElem(level);

 (*output)<< std::setw(10) << 1 << std::setw(10) << 1 ;

// !!!!!! Pay attention that there is switch

 if (ptCoordinate->is2D()) (*output) << std::setw(10) << 2 << std::endl; 
 else (*output) << std::setw(10) << 3 << std::endl; 

 (*output) << std::setw(10) << maxnumnodes;
 (*output) << std::setw(10) << maxnumelem << std::endl;
 
 (*output) << std::setw(6) << -1 << std::endl;
}


template<class Dim>
void  OutResultUnverg<Dim>::Dataset781(Grid<Dim> * ptgrid, const Integer level)
{
 (*output) << std::setw(6) << -1 << std::endl << std::setw(6) << 781 << std::endl;

 Integer maxnumnodes=ptgrid-> GetMaxnumnodes(level);

 (*output).setf(std::ios::scientific);
 (*output).precision(16);

// ptCoordinate=ptgrid->GetptCoordinate(level); 

 for (Integer i=0; i<maxnumnodes; i++)
   {
     (*output) << std::setw(10) << i+1 << std::setw(10) << 0 << std::setw(10) << 0 << std::setw(10) << 11 << std::endl;

     (*output).setf(std::ios::uppercase);

     Dim Point;

     ptgrid->GetCoordinateNode(i,level,Point);

     (*output) << "   " << 0.0 ;

     PrintPoint(Point,output);

     (*output) << std::endl;           
   }
 
 (*output) << std::setw(6) << -1 << std::endl;
}

template<class Dim>
void  OutResultUnverg<Dim>::Dataset780(Grid<Dim> * ptgrid, const Integer level)
{
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

     numnodesPerElem=ptgrid->GetNumNodesPerElem(i,level);

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
void  OutResultUnverg<Dim>::Dataset55(const std::string & title, const Vector<Double> & x, const Integer step, const Double time)
{
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
void  OutResultUnverg<Dim>::Dataset56()
{
 (*output) << std::setw(6) << -1 << std::endl << std::setw(6) << 56 << std::endl;
 
 (*output) << std::setw(6) << -1;
}

} // end of namespace
