#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include "outUnverg.hh"

namespace CoupledField
{

WriteResultsUnverg :: WriteResultsUnverg(const Char * const filename)
:WriteResults(filename)
{
#ifdef TRACE
  (*trace) << "entering WriteResultsUnverg :: WriteResultsUnverg" << std::endl;
#endif

  output=new std::ofstream(strcat(namefile_,".unverg"));
}

WriteResultsUnverg ::~WriteResultsUnverg()
{
#ifdef TRACE
 (*trace) << "entering WriteResultsUnverg ::~ WriteResultsUnverg" << std::endl;
#endif
  
 delete output;
}

void WriteResultsUnverg :: WriteGrid(const Integer level)
{
#ifdef TRACE
  (*trace) << " entering WriteResultsUnverg :: WriteGrid " << std::endl;
#endif
 if (!output) Error(" File for output results is not initialized");
 Dataset666(level);
 Dataset781(level);
 Dataset780(level);
}


void  WriteResultsUnverg::Dataset666(const Integer level)
{
 //
 if (!ptgrid)
    Error("ptgrid is not initialized", __FILE__,__LINE__);

 (*output)<< std::setw(6) << -1 << std::endl << std::setw(6) << -666 << std::endl ;

 Integer dim=ptgrid->GetDim();
 Integer maxnumnodes=ptgrid-> GetMaxnumnodes(level);
 Integer maxnumelem=ptgrid-> GetMaxnumElem(level);

 (*output)<< std::setw(10) << 1 << std::setw(10) << 1 << std::setw(10) << dim << std::endl << std::setw(10) << maxnumnodes << std::setw(10) << maxnumelem << std::endl;

 (*output) << std::setw(6) << -1 << std::endl;
}

void  WriteResultsUnverg::Dataset781(const Integer level)
{
  //
 if (!ptgrid)
    Error("ptgrid is not initialized", __FILE__,__LINE__);

 (*output) << std::setw(6) << -1 << std::endl << std::setw(6) << 781 << std::endl;

 Integer dim=ptgrid->GetDim();
 Integer maxnumnodes=ptgrid->GetMaxnumnodes(level);
 std::cout << " maxnumnodes " << maxnumnodes << std::endl;

 (*output).setf(std::ios::scientific);
 (*output).precision(16);

 for (Integer i=0; i<maxnumnodes; i++)
   {
     (*output) << std::setw(10) << i+1 << std::setw(10) << 0 << std::setw(10) << 0 << std::setw(10) << 11 << std::endl;

     (*output).setf(std::ios::uppercase);

     if (dim==2)
     { Point2D Point;
        ptgrid->GetCoordinateNode(i,level,Point);

        (*output) << "   " << 0.0 ;

         PrintPoint(Point,output);

         (*output) << std::endl;
      }
      else
     { Point3D Point;

     ptgrid->GetCoordinateNode(i,level,Point);

     PrintPoint(Point,output);

     (*output) << std::endl;  
     }         
   }
 
 (*output) << std::setw(6) << -1 << std::endl;
}

void  WriteResultsUnverg::Dataset780(const Integer level)
{
  //
 if (!ptgrid)
    Error("ptgrid is not initialized", __FILE__,__LINE__);

 (*output) << std::setw(6) << -1 << std::endl << std::setw(6) << 780 << std::endl;
 Integer dim=ptgrid->GetDim();

 Integer maxnumelem=ptgrid->GetMaxnumElem(level);

 std::cout << " maxnumelem " << maxnumelem << std::endl;

 Vector<Integer> connect;

 for (Integer i=0; i<maxnumelem; i++)
   {
     (*output) << std::setw(10) << i+1 << std::setw(10);

     ptgrid->GetConnection(connect, i, level);

     if (dim==2)
{     switch(connect.size())
     {
       case 3: (*output) << 91 ; break;
       case 4: (*output) << 94 ; break;
       case 6: (*output) << 92; break;
       case 8: (*output) << 95; break;
       default: Error("Please, put element type according to unverg-format for this number of nodes per element", __FILE__,__LINE__);
     }

     (*output) << std::setw(10) << 2 << std::setw(10) << 2 << std::setw(10) << 1 << std::setw(10) << 1 << std::setw(10) << 1 << std::setw(10) << connect.size()
<< std::endl;
}
   else
{
     switch(connect.size())
     {
       case 4: (*output) << 111 ; break;
       case 6: (*output) << 112; break;
       case 8: (*output) << 115; break;
       case 15: (*output) << 113; break;
       case 20: (*output) << 116; break;
       default: Error("Please, put element type according to unverg-format for this number of nodes per element", __FILE__,__LINE__);
     }

     (*output) << std::setw(10) << 11 << std::setw(10) << 1 << std::setw(10) <<
1 << std::setw(10) << 1 << std::setw(10) << 1 << std::setw(10) << connect.size() << std::endl;
}

     for (Integer ii=0; ii < connect.size(); ii++) 
       { 
	 (*output).width(10);
	 (*output) << connect[ii];
       }

     (*output) << std::endl;
   }
 (*output) << std::setw(6) << -1 << std::endl;
}

void  WriteResultsUnverg::Dataset55(const std::string & title, const Vector<Double> & x, const Integer step, const Double time)
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

void  WriteResultsUnverg::Dataset56()
{
 (*output) << std::setw(6) << -1 << std::endl << std::setw(6) << 56 << std::endl;
 
 (*output) << std::setw(6) << -1;
}

void  WriteResultsUnverg::Init(Grid * aptgrid)
{
 ptgrid=aptgrid;
}

void  WriteResultsUnverg::WriteSolution(const Vector<Double> & sol, const Integer step, const Double time, const std::string title)
{

 Integer i;
 if (NeedHistory_) 
   for (i=0; i< nodeshist_.size(); i++)
    {
      if (sol.size() <= nodeshist_[i])
        Error("Please, check history-nodes in config-file.",__FILE__,__LINE__);
      if (time != lastsavetime[i])
          AddInHistory(time,sol[nodeshist_[i]],i); 
    }
 Dataset55(title, sol, step+1, time);
}

} // end of namespace
