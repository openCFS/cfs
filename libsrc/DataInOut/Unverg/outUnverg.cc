#include <iostream.h>
#include <fstream.h>
#include <string>
#include <iomanip.h>

#include <general_head.hh>
#include <utils_head.hh>
#include <domain_head.hh>
#include "outUnverg.hh"

namespace CoupledField
{

//template<class Point2D>
OutResultUnverg :: OutResultUnverg(const Char * const filename)
{
#ifdef TRACE
  (*trace) << "entering OutResultUnverg :: OutResultUnverg" << endl;
#endif

  Char * help = new Char[20];
  strcpy(help,filename);
  output=new ofstream(strcat(help,".unverg"));

}

//template<class Point2D>
OutResultUnverg ::~ OutResultUnverg()
{
#ifdef TRACE
  (*trace) << "entering OutResultUnverg ::~ OutResultUnverg" << endl;
#endif

  delete output;
}

//template<class Point2D>
void OutResultUnverg :: Create(Grid<Point2D> * ptgrid, const Integer level)
{
 Dataset666(ptgrid,level);
 Dataset781(ptgrid,level);
 Dataset780(ptgrid,level);
// Dataset55(ptgrid,level);
// Dataset56(ptgrid,level);
 ;
}


//template<class Point2D>
void  OutResultUnverg::Dataset666(Grid<Point2D> * ptgrid, const Integer level)
{

 (*output)<< setw(6) << -1 << endl << setw(6) << -666 << endl ;

 Integer maxnumnodes=ptgrid-> GetMaxnumnodes(level);
 Integer maxnumelem=ptgrid-> GetMaxnumElem(level);

 (*output)<< setw(10) << 1 << setw(10) << 1 ;

 if (ptCoordinate->is2D()) (*output) << setw(10) << 2 << endl; 
 else (*output) << setw(10) << 3 << endl; 

 (*output) << setw(10) << maxnumnodes;
 (*output) << setw(10) << maxnumelem << endl;
 
 (*output) << setw(6) << -1 << endl;
}


//template<class Point2D>
void  OutResultUnverg::Dataset781(Grid<Point2D> * ptgrid, const Integer level)
{
 (*output) << setw(6) << -1 << endl << setw(6) << 781 << endl;

 Integer maxnumnodes=ptgrid-> GetMaxnumnodes(level);

 (*output).setf(ios::scientific);
 (*output).precision(16);

 ptCoordinate=ptgrid->GetptCoordinate(level); 

 for (Integer i=0; i<maxnumnodes; i++)
   {
     (*output) << setw(10) << i+1 << setw(10) << 0 << setw(10) << 0 << setw(10) << 11 << endl;

     (*output).setf(ios::uppercase);
     if (ptCoordinate->is2D()) 
       {
	 (*output) << "   " << 0.0 ; 
	 PrintPoint(ptCoordinate[i],output);
       }
     else  PrintPoint(ptCoordinate[i],output);


     (*output) << endl;           
   }
 
 (*output) << setw(6) << -1 << endl;
}

//template<class Point2D>
void  OutResultUnverg::Dataset780(Grid<Point2D> * ptgrid, const Integer level)
{
 (*output) << setw(6) << -1 << endl << setw(6) << 780 << endl;

 Integer maxnumnodes=ptgrid-> GetMaxnumnodes(level);
 Integer maxnumelem=ptgrid-> GetMaxnumElem(level);

 Integer * ConnectElem=NULL;
 Integer numnodesPerElem=0;

 for (Integer i=0; i<maxnumelem; i++)
   {
     numnodesPerElem=ptgrid->GetNumNodesPerElem(i,level); 

     (*output) << setw(10) << i+1 << setw(10) << 94 << setw(10) << 2 << setw(10) << 2 << setw(10) << 1 << setw(10) << 1 << setw(10) << 1 << setw(10) << numnodesPerElem << endl;

     numnodesPerElem=ptgrid->GetNumNodesPerElem(i,level);

     ConnectElem=new Integer[numnodesPerElem];      
     ptgrid->GetConnection(ConnectElem, level, i, numnodesPerElem);

     for (Integer ii=0; ii < numnodesPerElem; ii++) 
       { 
	 (*output).width(10);
	 (*output) << ConnectElem[ii];
       }

     delete [] ConnectElem;
 
     (*output) << endl;
   }
 (*output) << setw(6) << -1 << endl;
}

//template<class Point2D>
void  OutResultUnverg::Dataset55(const string & title, const Vector<Double> & x, const Integer step, const Double time)
{
 (*output) << setw(6) << -1 << endl << setw(6) << 55 << endl;

 (*output).setf(ios::scientific);
 (*output).precision(6);
 (*output).setf(ios::uppercase);

 (*output) << " " << title << ", step" << setw(6) << step <<
             " time   " << time << endl;  
 (*output) << endl << endl << endl << endl;
 (*output) << setw(10) << 1 << setw(10) << 4 << setw(10) << 1 << setw(10) << 0
           << setw(10) << 2 << setw(10) << 1 << endl;
 (*output) << setw(10) << 2 << setw(10) << 1 << setw(10) << 1 << setw(10) <<
              step << endl;
 (*output) << " " << time << endl;       

 Integer i,n;
 n=x.size();  
 for (i=0; i<n; i++)
   (*output) << setw(10) << i+1 << endl << " " << x[i] << endl;
     
 (*output) << setw(6) << -1 << endl;
}  

//template<class Point2D>
void  OutResultUnverg::Dataset56()
{
 (*output) << setw(6) << -1 << endl << setw(6) << 56 << endl;
 
 (*output) << setw(6) << -1;
}

} // end of namespace
