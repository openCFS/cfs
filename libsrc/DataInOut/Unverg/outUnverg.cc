#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include "outUnverg.hh"

namespace CoupledField
{

WriteResultsUnverg :: WriteResultsUnverg(const Char * const filename, Boolean withHistory, FileType * const aInFile)
:WriteResults(filename, withHistory, aInFile)
{
#ifdef TRACE
  (*trace) << "entering WriteResultsUnverg :: WriteResultsUnverg" << std::endl;
#endif

  output = NULL;
  if (!withHistory)
    output=new std::ofstream(strcat(namefile_,".unv"));
}

WriteResultsUnverg ::~WriteResultsUnverg()
{
#ifdef TRACE
 (*trace) << "entering WriteResultsUnverg :: ~WriteResultsUnverg" << std::endl;
#endif
 if (output)
   {
     output->close();
     //   delete output;
   }
}

void WriteResultsUnverg :: WriteGrid(const Integer level)
{
#ifdef TRACE
  (*trace) << " entering WriteResultsUnverg :: WriteGrid " << std::endl;
#endif

 if (!NeedHistory_)
   {    
     if (!output) 
       Error(" File for output results is not initialized");

     Dataset666(level);
     Dataset781(level);
     Dataset780(level);
   }
 
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

  StdVector<Integer> connect;
  StdVector<Elem*> elemssd;
  Integer elmsgrp=1;

  StdVector<std::string>* subdoms;
  subdoms=ptgrid->GetAllSDs();
  Integer i, j, k;
  k = 0;
  for (i=0; i<subdoms->GetSize(); i++)
    {
      ptgrid->GetElemSD(elemssd,(*subdoms)[i],level);

      for (j=0; j < elemssd.GetSize(); j++)
	{  
	  k++; 
	  connect=elemssd[j]->connect;

	  (*output) << std::setw(10) << elemssd[j]->elemNum << std::setw(10);

	  if (dim==2)
	    {     switch(connect.GetSize())
	      {
	      case 3: (*output) << 91 ; break;
	      case 4: (*output) << 94 ; break;
	      case 6: (*output) << 92; break;
	      case 8: (*output) << 95; break;
	      default: Error("Please, put element type according to unverg-format for this number of nodes per element", __FILE__,__LINE__);
	      }

	    (*output) << std::setw(10) << 2 << std::setw(10) << 2 << std::setw(10) << 1 << std::setw(10) << 1 << std::setw(10) << elmsgrp << std::setw(10) << connect.GetSize() << std::endl;
	    }
	  else
	    {
	      switch(connect.GetSize())
		{
		case 4: (*output) << 111 ; break;
		case 6: (*output) << 112; break;
		case 8: (*output) << 115; break;
		case 15: (*output) << 113; break;
		case 20: (*output) << 116; break;
		default: Error("Please, put element type according to unverg-format for this number of nodes per element", __FILE__,__LINE__);
		}

	      (*output) << std::setw(10) << 11 << std::setw(10) << 1 << std::setw(10) << 1
			<< std::setw(10) << 1 << std::setw(10) << elmsgrp << std::setw(10) 
			<< connect.GetSize() << std::endl;
	    }

	  if (dim == 2 && (connect.GetSize() == 6 || connect.GetSize() == 8))
	    {
	      //quadratic elements
	      Integer offset = Integer(connect.GetSize()/2);
	      for (Integer ii=0; ii < offset; ii++)
		{
		 (*output).width(10);
		 (*output) << connect[ii];
		 (*output).width(10);
		 (*output) << connect[ii+offset]; 
		}

	    }
	  else
	    {
	      for (Integer ii=0; ii < connect.GetSize(); ii++) 
		{ 
		  (*output).width(10);
		  (*output) << connect[ii];
		}
	    }

	  (*output) << std::endl;
	} // over elements of group
      elmsgrp++;
    } // over groups
  (*output) << std::setw(6) << -1 << std::endl;
}

void  WriteResultsUnverg::Dataset55(const std::string & title, 
				    const Vector<Double> & x, 
				    const Integer step, 
				    const Double time, 
				    const Integer nrNodes,
				    const Integer nrDofs)
{
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
 (*output) << std::endl << std::endl << std::endl << std::endl;
 (*output) << std::setw(10) << 1 << std::setw(10) << 4 << std::setw(10) << 1 << std::setw(10) << 0
           << std::setw(10) << 2 << std::setw(10) << valsPerNode << std::endl;
 (*output) << std::setw(10) << 2 << std::setw(10) << 1 << std::setw(10) << 1 << std::setw(10) <<
              step << std::endl;
 (*output) << " " << time << std::endl;       

 Integer i,j,n;
 n=nrNodes;;  
 for (i=0; i<n; i++)
   {
     (*output) << std::setw(10) << i+1 << std::endl;
     
     // in the universal file either one or three results datas must exist
     if (nrDofs == 2)
       (*output) << 0.0;

     for (j=0; j<nrDofs; j++)
       {
	 //std::cerr << "trying to write " << i << ", " << j << std::endl;
	 (*output) << std::setw(14) << x[i*nrDofs +j];
       }
     
     (*output) << std::endl;
   }    
 (*output) << std::setw(6) << -1 << std::endl;
}  

void  WriteResultsUnverg::Dataset56(const std::string & title, 
				    const Vector<Double> & x, 
				    const Integer step, 
				    const Double time, 
				    const Integer numElems,
				    const Integer nrDofs)
{
  
   if (!ptgrid)
      Error("ptgrid is not initialized", __FILE__,__LINE__);

  (*output) << std::setw(6) << -1 << std::endl << std::setw(6) << 56 << std::endl;

  (*output).setf(std::ios::scientific);
  (*output).precision(6);
  (*output).setf(std::ios::uppercase);

  Integer valsPerNode = 1;
  if (nrDofs > 1)
    valsPerNode = 3;

  (*output) << " " << title << ", step" << std::setw(6) << step <<
              " time   " << time << std::endl;  
  (*output) << std::endl << std::endl << std::endl << std::endl;
  (*output) << std::setw(10) << 1 << std::setw(10) << 4 << std::setw(10) << 2 << std::setw(10) << 0
            << std::setw(10) << 2 << std::setw(10) << valsPerNode << std::endl;
  (*output) << std::setw(10) << 2 << std::setw(10) << 1 << std::setw(10) << 1 << std::setw(10) <<
               step << std::endl;
  (*output) << " " << time << std::endl;       

  Integer i,j,n;
  n=numElems;  

  // for 2-dimensional solution, the plane has to be rotated
  if (nrDofs == 2)
    {
    for (i=0; i<n; i++)
 	{
	  (*output) << std::setw(10) << i+1 << std::setw(10) << 3 << std::endl;

	  (*output) << std::setw(13) << 0.0 << std::setw(13) << x[i*nrDofs];
	  (*output)<< std::setw(13) << x[i*nrDofs+1] << std::endl;
	}
    } 
  else
    {
      for (i=0; i<n; i++)
 	{
	  (*output) << std::setw(10) << i+1 << std::setw(10) << 3 << std::endl;
	  for (j=0; j<nrDofs; j++)
	    (*output) << std::setw(14) << x[i*nrDofs + j];
	  
	  (*output) << std::endl;
 	}
    }

 (*output) << std::setw(6) << -1 << std::endl;
}  

void  WriteResultsUnverg::Init(Grid * aptgrid)
{
 ptgrid=aptgrid;
}

void  WriteResultsUnverg::WriteNodeSolution(const NodeStoreSol<Double> & sol, 
					    const Integer step, 
					    const Double time, 
					    const std::string title)
{

 Integer i,j;
 Integer nrDofs = 1;
 Double help;
 
 Vector<Double> globalSolution;

 // Transform local nodal solution to global one
 // WARNING: Level for refinemet is hardcoded to 1
 sol.TransformNodeSolution(globalSolution,ptgrid,1);

 Integer numNodes =  ptgrid->GetMaxnumnodes(1);
 
 if (NeedHistory_) 
   for (i=0; i< nodeshist_.GetSize(); i++)
     {
      if (sol.GetDof() * sol.GetNumNodes() <= nodeshist_[i])
        Error("Please, check history-nodes in config-file.",__FILE__,__LINE__);
      if (lastsavetime[i] != time )
	if (nrDofs > 1)	
	  {
	    Vector<Double> solVec;
	    solVec.Resize(nrDofs);
	    for (j=0; j<nrDofs; j++)
	      sol.Get(nodeshist_[i]-1,j,solVec[j]);

	    AddVecInHistory(time, solVec, i);
	  }
	else
	  {
	    sol.Get(nodeshist_[i]-1,0,help);
	    AddInHistory(time,help,i);
	  }
     }
 else
   Dataset55(title, globalSolution, step, time, numNodes ,sol.GetDof());
}


void  WriteResultsUnverg::WriteElemSolution(const ElemStoreSol<Double>& sol, const Integer step, const Double time, const std::string title)
{
#ifdef TRACE
  (*trace) << " entering WriteResultsUnverg::WriteElemSolution " << std::endl;
#endif

  Vector<Double> globalSolution;
  Integer numElems =  ptgrid->GetMaxnumElem(1);
  // Transform local nodal solution to global one
  // WARNING: Level for refinemet is hardcoded to 1
  sol.TransformElemSolution(globalSolution,ptgrid,1);
  
   if (!NeedHistory_)
     Dataset56(title, globalSolution, step, time, numElems, sol.GetDof());
}



// void  WriteResultsUnverg::WriteDataOnCell(const Vector<Double> & sol, const Integer step, const Double time, const std::string title, const Integer nrDofs)
// {
//   std::string aux;
//   if (title == "elecfield") {
//     aux="electric field";
//     Dataset56(aux, sol, step+1, time);
//   }
//   else 
//     Dataset56(title, sol, step+1, time, nrDofs);
//     //Warning("This cell-data is not printed, since this type is not supported by Capapost",__FILE__,__LINE__);
// }

} // end of namespace
