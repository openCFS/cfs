#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include "outUnverg.hh"

namespace CoupledField
{

WriteResultsUnverg :: WriteResultsUnverg(const Char * const filename, Boolean withHistory)
:WriteResults(filename,withHistory)
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

 (*output).setf(std::ios::scientific);
 (*output).precision(16);

 for (Integer i=0; i<maxnumnodes; i++)
   {
     (*output) << std::setw(10) << i+1 << std::setw(10) << 0 << std::setw(10) << 0 << std::setw(10) << 11 << std::endl;

     (*output).setf(std::ios::uppercase);

	if (dim==2)
{
 	Point<2> Point;
        ptgrid->GetCoordinateNode(i,level,Point);

        (*output) << "   " << 0.0 ;
         PrintPoint(Point,output);
         (*output) << std::endl;
}
      else
{
	Point<3> Point;
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

void  WriteResultsUnverg::Dataset55(const std::string & title, const Vector<Double> & x, const Integer step, const Double time, const Integer nrDofs)
{
  const Integer nrBlankLinesIn55 = 8;
  
  //
  if (!ptgrid)
     Error("ptgrid is not initialized", __FILE__,__LINE__);

 (*output) << std::setw(6) << -1 << std::endl << std::setw(6) << 55 << std::endl;

 (*output).setf(std::ios::scientific);
 (*output).precision(6);
 (*output).setf(std::ios::uppercase);

 Integer valsPerNode = 1;
 if (nrDofs > 1)
   valsPerNode = 3;
 

 (*output) << " " << title << " step" << std::setw(6) << step <<
             " time   " << time << std::endl;  
 for(Integer i=0; i<nrBlankLinesIn55; i++)
   (*output) << std::endl;
 
 (*output) << std::setw(10) << 1 << std::setw(10) << 4 << std::setw(10) << 1 << std::setw(10) << 0
           << std::setw(10) << 2 << std::setw(10) << valsPerNode << std::endl;
 (*output) << std::setw(10) << 2 << std::setw(10) << 1 << std::setw(10) << 1 << std::setw(10) <<
              step << std::endl;
 (*output) << " " << time << std::endl;       

 Integer i,j,n;
 n=x.size()/nrDofs;  
 for (i=0; i<n; i++)
   {
     (*output) << std::setw(10) << i+1 << std::endl << " ";

     // in the universal file eihter one or three results datas must exist
     if (nrDofs == 1)
       (*output) << 0.0 << " " << 0.0;

     // in the universal file eihter one or three results datas must exist
     if (nrDofs == 2)
       (*output) << 0.0;

     for (j=0; j<nrDofs; j++)
       (*output) << " " << x[i*nrDofs + j];
     
     (*output) << std::endl;
   }    
 (*output) << std::setw(6) << -1 << std::endl;
}  



void  WriteResultsUnverg::Write56Header(const std::string & title, const Integer step, const Double time)
{
  const Integer nrBlankLinesIn56 = 8;
  //
  if (!ptgrid)
     Error("ptgrid is not initialized", __FILE__,__LINE__);

 (*output) << std::setw(6) << -1 << std::endl << std::setw(6) << 56 << std::endl;

 (*output).setf(std::ios::scientific);
 (*output).precision(6);
 (*output).setf(std::ios::uppercase);

 (*output) << " " << title << ", step" << std::setw(6) << step <<
             " time   " << time << std::endl;  

 for(Integer i=0; i<nrBlankLinesIn56; i++)
   (*output) << std::endl;

 (*output) << std::setw(10) << 1 << std::setw(10) << 4 << std::setw(10) << 1 << std::setw(10) << 0
           << std::setw(10) << 2 << std::setw(10) << 3 << std::endl;
 (*output) << std::setw(10) << 2 << std::setw(10) << 1 << std::setw(10) << 1 << std::setw(10) <<
              step << std::endl;
 (*output) << " " << time << std::endl;     
}



// writes the header of a 56 element dataset
void  WriteResultsUnverg::Dataset56(const std::string & title, const Vector<Double> & x, const Integer step, const Double time)
{
  Write56Header(title, step+1, time);
  

  Integer i,n;
  n=x.size();  
  for (i=0; i<n; i++)
    (*output) << std::setw(10) << i+1 << std::setw(10) << 3 << std::endl << "   " << 0.0 << "    " << 0.0 << "   " << x[i] << std::endl;
  
  (*output) << std::setw(6) << -1 << std::endl;
}  




// writes the header of a 56 element dataset
void  WriteResultsUnverg::Dataset56(const std::string & title, const Matrix<Double> & x, const Integer step, const Double time)
{
  Write56Header(title, step+1, time);

  const Integer  nrRows = x.size_row();  
  for (Integer i=0; i<nrRows; i++)
    (*output) << std::setw(10) << i+1 << std::setw(10) << 3 << std::endl 
	      << "   " << x[i][0] << "    " << x[i][1] << "   " << x[i][2] << std::endl;
  
  (*output) << std::setw(6) << -1 << std::endl;
}  






void  WriteResultsUnverg::Init(Grid * aptgrid)
{
 ptgrid=aptgrid;
}

void  WriteResultsUnverg::WriteSolution(const Vector<Double> & sol, const Integer step, const Double time, const std::string title, const Integer nrDofs)
{

 Integer i,j;
 if (NeedHistory_) 
   for (i=0; i< nodeshist_.size(); i++)
     {
      if (sol.size() <= nodeshist_[i])
        Error("Please, check history-nodes in config-file.",__FILE__,__LINE__);
      if (lastsavetime[i] != time )
	if (nrDofs > 1)	
	  {
	    std::vector<Double> solVec;
	    solVec.resize(nrDofs);
	    for (j=0; j<nrDofs; j++)
	      solVec[j] = sol[(nodeshist_[i]-1) * nrDofs + j];
	    
	    AddVecInHistory(time, solVec, i);
	  }
	else
	  // sol[node-1] since internal node starts at zero!        
	  AddInHistory(time,sol[nodeshist_[i]-1],i);

    }
 Dataset55(title, sol, step+1, time, nrDofs);
}

void  WriteResultsUnverg::WriteDataOnCell(const Vector<Double> & sol, const Integer step, const Double time, const std::string title)
{
  std::string aux;
  if (title == "elecfield") {
    aux="electric field";
    Dataset56(aux, sol, step+1, time);
  } else
    if (title == "magnetic field density") 
      Dataset56(title, sol, step+1, time);

  else 
    Warning("This cell-data is not printed, since this type is not supported by Capapost",__FILE__,__LINE__);
}



void  WriteResultsUnverg::WriteDataOnCell(const Matrix<Double> & sol, const Integer step, const Double time, const std::string title)
{
  if (title == "magnetic field density") 
    Dataset56(title, sol, step+1, time);
  else 
    Warning("This cell-data is not printed, since this type is not supported!",__FILE__,__LINE__);
}

} // end of namespace
