#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include "outUnverg.hh"

namespace CoupledField
{

WriteResultsUnverg :: WriteResultsUnverg(const Char * const filename,
					 FileType * const aInFile)
:WriteResults(filename, aInFile)
{
  ENTER_FCN( "WriteResultsUnverg::WriteResultsUnverg" );

  std::string name = namefile_ + ".unv";
  output = NULL;
  output=new std::ofstream(name.c_str());
 }

WriteResultsUnverg ::~WriteResultsUnverg()
{
  ENTER_FCN( "WriteResultsUnverg::~WriteResultsUnverg" );
 if (output)
   {
     output->close();
     //   delete output;
   }
}

void WriteResultsUnverg :: WriteGrid(const Integer level)
{
  ENTER_FCN( "WriteResultsUnverg::WriteGrid" );

     if (!output) 
       Error(" File for output results is not initialized");

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

void  WriteResultsUnverg::NodeElemDataTransient(const Integer dataSetNr,
						const std::string & title, 
						const Vector<Double> & x, 
						const Integer step, 
						const Double time, 
						const Integer nrNodes,
						const Integer nrDofs)
{
  //
  if (!ptgrid)
     Error("ptgrid is not initialized", __FILE__,__LINE__);

 (*output) << std::setw(6) << -1 << std::endl 
	   << std::setw(6) << dataSetNr << std::endl;

 (*output).setf(std::ios::scientific);
 (*output).precision(6);
 (*output).setf(std::ios::uppercase);

 // Standard for scalar values
 Integer dataCharac = 1;
 Integer valsPerNode = 1;
 
 // needed for undocumented value of
 // Dataset 55/56 in record8 , field4
 Integer specDataCharac = 0;

 // Vector type
 if (nrDofs > 1 && nrDofs <= 3) 
   {
   valsPerNode = 3;
   dataCharac = 2;
 } 
 // Tensor
 else if(nrDofs == 6) 
   { 
     valsPerNode = 6;
     specDataCharac = 2;
     dataCharac = 4; // symmetric tensor
//      dataCharac = 3; // vector with 6 components
//      dataCharac = 5; // unsymmetric tensor
   }
     
 
 (*output) << " " << title << " step" << std::setw(6) << step <<
             " time   " << time << std::endl;  
 (*output) << std::endl << std::endl << std::endl << std::endl;
 (*output) << std::setw(10) << 1 << std::setw(10) << 4 << std::setw(10) 
	   << dataCharac  << std::setw(10) << specDataCharac
           << std::setw(10) << 2 << std::setw(10) << valsPerNode << std::endl;
 (*output) << std::setw(10) << 2 << std::setw(10) << 1 << std::setw(10) << 1 
	   << std::setw(10) <<
              step << std::endl;
 (*output) << " " << time << std::endl;       

 Integer i,j,n;
 n=nrNodes;;  
 for (i=0; i<n; i++)
   {
     
     (*output) << std::setw(10) << i+1;
     if (dataSetNr == 56)
       (*output) << std::setw(10) << valsPerNode;

     (*output) << std::endl;
     
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

void WriteResultsUnverg::NodeElemDataHarmonic(const Integer dataSetNr,
					      const std::string & title, 
					      const Vector<Complex> & x, 
					      const Integer step,
					      const Double frequency, 
					      const ComplexFormat format,
					      const Integer nrNodes,
					      const Integer nrDofs)
{
  
  Integer dataCharact = 1;
  if (!ptgrid)
    Error("ptgrid is not initialized", __FILE__,__LINE__);
  
  (*output) << std::setw(6) << -1 << std::endl 
	    << std::setw(6) << dataSetNr << std::endl;
  
  (*output).setf(std::ios::scientific);
  (*output).precision(6);
  (*output).setf(std::ios::uppercase);
  
  Integer valsPerNode = 1;
  if (nrDofs > 1)
    {
      dataCharact = 2;
      valsPerNode = 3;
    }
 

  if (format == REAL_IMAG)
    {
      // write out realpart
      (*output) << " " << title << " cw realpart" << std::setw(6) <<" frequency   " << frequency << std::endl;  
      (*output) << std::endl << std::endl << std::endl << std::endl;
      (*output) << std::setw(10) << 1 << std::setw(10) << 5 << std::setw(10) << dataCharact << std::setw(10) << 0
		<< std::setw(10) << 2 << std::setw(10) << valsPerNode << std::endl;
      (*output) << std::setw(10) << 2 << std::setw(10) << 1 << std::setw(10) << -2 << std::setw(10) <<
	step << std::endl;
      (*output) << " " << frequency << std::endl;       
      
      Integer i,j,n;
      n=nrNodes;
      for (i=0; i<n; i++)
	{
	  (*output) << std::setw(10) << i+1 << std::endl;
	  
	  // in the universal file either one or three results datas must exist
	  if (nrDofs == 2)
	    (*output) << 0.0;
	  
	  for (j=0; j<nrDofs; j++)
	    {
	      //std::cerr << "trying to write " << i << ", " << j << std::endl;
	      (*output) << std::setw(14) << real(x[i*nrDofs +j]);
	    }
	  
	  (*output) << std::endl;
	}    
      (*output) << std::setw(6) << -1 << std::endl;

      // write out imag part
      (*output) << std::setw(6) << -1 << std::endl << std::setw(6) << 55 << std::endl;
      (*output) << " " << title << " cw imagpart" << std::setw(6) <<" frequency   " << frequency << std::endl;  
      (*output) << std::endl << std::endl << std::endl << std::endl;
      (*output) << std::setw(10) << 1 << std::setw(10) << 5 << std::setw(10) << dataCharact << std::setw(10) << 0
		<< std::setw(10) << 2 << std::setw(10) << valsPerNode << std::endl;
      (*output) << std::setw(10) << 2 << std::setw(10) << 1 << std::setw(10) << -12 << std::setw(10) <<
	step << std::endl;
      (*output) << " " << frequency << std::endl;       
      
      for (i=0; i<n; i++)
	{
	  (*output) << std::setw(10) << i+1 << std::endl;
	  
	  // in the universal file either one or three results datas must exist
	  if (nrDofs == 2)
	    (*output) << 0.0;
	  
	  for (j=0; j<nrDofs; j++)
	    {
	      //std::cerr << "trying to write " << i << ", " << j << std::endl;
	      (*output) << std::setw(14) << imag(x[i*nrDofs +j]);
	    }
	  
	  (*output) << std::endl;
	}    
      (*output) << std::setw(6) << -1 << std::endl;
    }
  
  else if (format == AMPLITUDE_PHASE) {
    // write out amplitude
    (*output) << " " << title << " cw amplitude" << std::setw(6) <<" frequency   " << frequency << std::endl;  
    (*output) << std::endl << std::endl << std::endl << std::endl;
    (*output) << std::setw(10) << 1 << std::setw(10) << 5 << std::setw(10) << dataCharact << std::setw(10) << 0
	      << std::setw(10) << 2 << std::setw(10) << valsPerNode << std::endl;
    (*output) << std::setw(10) << 2 << std::setw(10) << 1 << std::setw(10) << -1 << std::setw(10) <<
      step << std::endl;
    (*output) << " " << frequency << std::endl;       
    
    Integer i,j,n;
    n=nrNodes;
    for (i=0; i<n; i++)
      {
	(*output) << std::setw(10) << i+1 << std::endl;
	
	  // in the universal file either one or three results datas must exist
	if (nrDofs == 2)
	  (*output) << 0.0;
	
	for (j=0; j<nrDofs; j++)
	  {
	    //std::cerr << "trying to write " << i << ", " << j << std::endl;
	    (*output) << std::setw(14) << abs(x[i*nrDofs +j]);
	  }
	
	(*output) << std::endl;
      }    
    (*output) << std::setw(6) << -1 << std::endl;
    
    // write out phase
    (*output) << std::setw(6) << -1 << std::endl << std::setw(6) << 55 << std::endl;
    (*output) << " " << title << " cw phase" << std::setw(6) <<" frequency   " << frequency << std::endl;  
    (*output) << std::endl << std::endl << std::endl << std::endl;
    (*output) << std::setw(10) << 1 << std::setw(10) << 5 << std::setw(10) << dataCharact << std::setw(10) << 0
	      << std::setw(10) << 2 << std::setw(10) << valsPerNode << std::endl;
    (*output) << std::setw(10) << 2 << std::setw(10) << 1 << std::setw(10) << -11 << std::setw(10) <<
      step << std::endl;
    (*output) << " " << frequency << std::endl;       
    
    for (i=0; i<n; i++)
      {
	(*output) << std::setw(10) << i+1 << std::endl;
	
	// in the universal file either one or three results datas must exist
	if (nrDofs == 2)
	  (*output) << 0.0;
	
	for (j=0; j<nrDofs; j++)
	  {
	    if (abs(imag(x[i*nrDofs +j])) > 1e-16)
	      (*output) << std::setw(14) << arg(x[i*nrDofs +j])*180/PI;
	    else 
	      (*output) << std::setw(14) << 0.0;
	  }
	
	(*output) << std::endl;
      }    
    (*output) << std::setw(6) << -1 << std::endl;
  }
  
}

// void  WriteResultsUnverg::Dataset56_Transient(const std::string & title, 
// 					      const Vector<Double> & x, 
// 					      const Integer step, 
// 					      const Double time, 
// 					      const Integer numElems,
// 					      const Integer nrDofs)
// {
  
//    if (!ptgrid)
//       Error("ptgrid is not initialized", __FILE__,__LINE__);

//   (*output) << std::setw(6) << -1 << std::endl << std::setw(6) << 56 << std::endl;

//   (*output).setf(std::ios::scientific);
//   (*output).precision(6);
//   (*output).setf(std::ios::uppercase);

//   Integer valsPerNode = nrDofs;
//   if (nrDofs ==2)
//     valsPerNode = 3;

//   //just for stresses
//   Integer stress1 =2;
//   Integer stress2 = 0;
//   if (nrDofs==6) {
//     stress1 = 4;
//     stress2 = 2;
//   }

//   (*output) << " " << title << ", step" << std::setw(6) << step <<
//               " time   " << time << std::endl;  
//   (*output) << std::endl << std::endl << std::endl << std::endl;
//   (*output) << std::setw(10) << 1 << std::setw(10) << 4 << std::setw(10) << stress1 << std::setw(10) << stress2
//             << std::setw(10) << 2 << std::setw(10) << valsPerNode << std::endl;
//   (*output) << std::setw(10) << 2 << std::setw(10) << 1 << std::setw(10) << 1 << std::setw(10) <<
//                step << std::endl;
//   (*output) << " " << time << std::endl;       

//   Integer i,j,n;
//   n=numElems;  

//   // for 2-dimensional solution, the plane has to be rotated
//   if (nrDofs == 2)
//     {
//     for (i=0; i<n; i++)
//  	{
// 	  (*output) << std::setw(10) << i+1 << std::setw(10) << 3 << std::endl;

// 	  (*output) << std::setw(13) << 0.0 << std::setw(13) << x[i*nrDofs];
// 	  (*output)<< std::setw(13) << x[i*nrDofs+1] << std::endl;
// 	}
//     } 
//   else
//     {
//       for (i=0; i<n; i++)
//  	{
// 	  (*output) << std::setw(10) << i+1 << std::setw(10) <<  nrDofs << std::endl;
// 	  for (j=0; j<nrDofs; j++)
// 	    (*output) << std::setw(14) << x[i*nrDofs + j];
	  
// 	  (*output) << std::endl;
//  	}
//     }

//  (*output) << std::setw(6) << -1 << std::endl;
// }  


// void WriteResultsUnverg::Dataset56_Harmonic(const std::string & title, 
// 					    const Vector<Complex> & x, 
// 					    const Integer step,
// 					    const Double frequency, 
// 					    const ComplexFormat format, 
// 					    const Integer numElems,
// 					    const Integer nrDofs)
// {

//   Error("WriteResultsUnverg::Dataset56_Harmonic: Not implemented yet",
// 	__FILE__, __LINE__);
// }


void  WriteResultsUnverg::Init(Grid * aptgrid, BCs * aptbcs)
{
  ptgrid=aptgrid;
  ptBCs_ = aptbcs;
}

void  WriteResultsUnverg::WriteNodeSolutionTransient(const NodeStoreSol<Double> & sol, 
						     const Integer step, 
						     const Double time)
{
  
  ENTER_FCN( "WriteResultsUnverg::WriteNodeSolutionTransient" );
    
  Vector<Double> globalSolution;
  StdVector<SolutionType> solTypes;
  Integer numNodes =  ptgrid->GetMaxnumnodes(1);
  std::string title;
  sol.GetSolutionTypes(solTypes);

 for (Integer iSol=0; iSol<solTypes.GetSize(); iSol++)
   {
     sol.GetGlobalSolVector(solTypes[iSol],globalSolution);
     title = SolutionTypeToString(solTypes[iSol]);
     
     NodeElemDataTransient(55,title, globalSolution, step, 
			   time, numNodes ,sol.GetDof(solTypes[iSol]));
   }
}

void  WriteResultsUnverg::WriteElemSolutionTransient(const ElemStoreSol<Double>& sol, 
						     const Integer step, 
						     const Double time)
{
  ENTER_FCN( "WriteResultsUnverg::WriteElemSolutionTransient" );

  Vector<Double> globalSolution;
  StdVector<SolutionType> solTypes;
  std::string title;
  Integer numElems =  ptgrid->GetMaxnumElem(0);  
  
  sol.GetSolutionTypes(solTypes);
  sol.TransformElemSolution(globalSolution,ptgrid,0);
  title = SolutionTypeToString(solTypes[0]);
  NodeElemDataTransient(56,title, globalSolution, step, 
			time, numElems, sol.GetDof());
}

void  WriteResultsUnverg::WriteNodeSolutionHarmonic(const NodeStoreSol<Complex> & sol, 
						    const Integer step,
						    const Double frequency, 
						    const ComplexFormat format)
{

  ENTER_FCN( "WriteResultsUnverg::WriteNodeSolutionHarmonic" );
  
  
  Integer i,j;
  Integer nrDofs = 1;
  Double help;
 
  Vector<Complex> globalSolution;
  
  StdVector<SolutionType> solTypes;
  sol.GetSolutionTypes(solTypes);
  
  Integer numNodes =  ptgrid->GetMaxnumnodes(1);
  std::string title;

  for (Integer iSol=0; iSol<solTypes.GetSize(); iSol++)
    {
      sol.GetGlobalSolVector(solTypes[iSol],globalSolution);
      title = SolutionTypeToString(solTypes[iSol]);
      NodeElemDataHarmonic(55, title, globalSolution, step, frequency, 
			   format, numNodes ,sol.GetDof(solTypes[iSol]));
    }
  
  
}


void  WriteResultsUnverg::WriteElemSolutionHarmonic(const ElemStoreSol<Complex>& sol, 
						    const Integer step,
						    const Double frequency, 
						    const ComplexFormat format)
{
  ENTER_FCN( "WriteResultsUnverg::WriteElemSolutionHarmonic" );
 Vector<Complex> globalSolution;
  StdVector<SolutionType> solTypes;
  Integer numElems =  ptgrid->GetMaxnumElem(0);  

  std::string title;

  sol.GetSolutionTypes(solTypes);
  sol.TransformElemSolution(globalSolution,ptgrid,0);
  title = SolutionTypeToString(solTypes[0]);  

  NodeElemDataHarmonic(55, title, globalSolution, step, frequency, 
			   format, numElems ,sol.GetDof(solTypes[0]));
}

std::string WriteResultsUnverg::SolutionTypeToString(const SolutionType type) const
{
  ENTER_FCN( "WriteResultsUnverg::SolutionTypeToString" );

  switch (type)
    {
    case MECH_DISPLACEMENT:
      return "displacement";
      break;
    case MECH_ACCELERATION:
      return "acceleration";
      break;
    case MECH_VELOCITY:
      return "velocity";
      break;
    case MECH_FORCE:
      break;
    case MECH_STRESS:
      return "stress";
      break;
    case MECH_STRAIN:
      Error("Not implemented", __FILE__, __LINE__);
      break;
    case ELEC_POTENTIAL:
      return "electric potential";
      break;
    case ELEC_FIELD_INTENSITY:
      return "electric field";
      break;
    case ELEC_FORCE_VWP: 
      Error("Not implemented", __FILE__, __LINE__);
      break;
    case ELEC_INTERFACE_FORCE:
      Error("Not implemented", __FILE__, __LINE__);
      break; 
    case ELEC_CHARGE:
      return "electric charge";
      break;
    case ELEC_FLUX_DENSITY:
      Error("Not implemented", __FILE__, __LINE__);
      break; 
    case ELEC_ENERGY:
      Error("Not implemented", __FILE__, __LINE__);
    case SMOOTH_DISPLACEMENT:
      return "displacement";
      break;
    case ACOU_POTENTIAL:
      return "fluid potential";
      break;
    case ACOU_FORCE:
      Error("Not implemented", __FILE__, __LINE__);
      break;
    case ACOU_POTENTIAL_DERIV_1:
      return "fluid potential, 1st deriv.";
      break;
    case ACOU_POTENTIAL_DERIV_2:
      return "fluid potential, 1st deriv.";
      break;
    case MAG_POTENTIAL:
      return "mag. vector potential";
      break;
    case MAG_FLUX_DENSITY:
      return "mag. flux density";
      break;
    case MAG_EDDY_CURRENT:
      return "eddy current";
      break;
    case MAG_FORCE_VWP:
      Error("Not implemented", __FILE__, __LINE__);
      break;
    case MAG_FORCE_LORENTZ:
      Error("Not implemented", __FILE__, __LINE__);
      break;
    case MAG_ENERGY:
      Error("Not implemented", __FILE__, __LINE__);
      break;
    default:
      Error( "Wrong type of solution or 'SolutionType2String' not implemented for\
this type of solution", __FILE__, __LINE__);
    }
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
