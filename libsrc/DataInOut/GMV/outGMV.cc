#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

#include "outGMV.hh"
#include "DataInOut/conffile.hh"

namespace CoupledField
{

WriteResultsGMV :: WriteResultsGMV(const Char * const filename, Boolean withHistory, FileType * const aInFile)
: WriteResults(filename,withHistory, aInFile)
{
#ifdef TRACE
  (*trace) << "entering WriteResultsGMV :: WriteResultsGMV" << std::endl;
#endif

 namedir_=new Char[30];

 Char S[50];

 strcpy(namedir_,filename);

 strcat(namedir_,"_gmv");
 sprintf(S,"mkdir -p %s",namedir_);

 system(S);

 currstep_=0;
 output=NULL;

 ptgrid=NULL;


 ascii_=FALSE;

 std::string typedata;
 if (conf->ifget("format_data_output",typedata)) {
   if (typedata == "binary")
     ascii_=FALSE;
   else ascii_=TRUE;
 }
 
 std::string flag;
 if (conf->ifget("fixed_grid",flag)) 
   {
   if (flag == "no")
     {
       fixedgrid_=FALSE;
       OpenFile(0);   
     }
   else
     { 
       fixedgrid_=TRUE;
       OpenFile(-1);

     }
   }
 if (!output)
    Error(" File for output results in .gmv-format ", __FILE__, __LINE__);

}

WriteResultsGMV ::~WriteResultsGMV()
{
#ifdef TRACE
  (*trace) << "entering WriteResultsGMV::~ WriteResultsGMV" << std::endl;
#endif

 // write keyword
 (*output) << "endgmv  ";
 if (ascii_) (*output) << std::endl;

 delete [] namedir_;

 delete output;
}

void WriteResultsGMV :: WriteHeader()
{
  (*output) << "gmvinput";
 if (ascii_) (*output) << " ascii" << std::endl;
 else (*output) << "ieeei4r8";

}

void WriteResultsGMV :: WriteNodes(const Integer alevel)
{
  Integer level=alevel;

  // write keyword
 (*output) << "nodev   ";

 //get and write number of nodes on the level
 Integer numnodes=ptgrid->GetMaxnumnodes(level);
 if (ascii_)
   (*output) << numnodes << std::endl;
 else
   output->write((char*)&numnodes,sizeof(Integer));

 Integer dim=ptgrid->GetDim();

 //get and write coodinates of nodes
 Integer i;

if (dim==2) 
  {
    Point<2> point;

 // write x,y,z-coordinate
  for (i=0; i<numnodes; i++)
    {
      ptgrid->GetCoordinateNode(i,level,point);

      if (ascii_)
      (*output) << " " << point[0] << " " << point[1] << " " << 0 << std::endl;
      else {
        Double z=0.;
        output->write((char*)&point[0],sizeof(Double));
        output->write((char*)&point[1],sizeof(Double));
        output->write((char*)&z,sizeof(Double));
      }
	
    }
}
else {
    Point<3> point;

 // write x,y,z-coordinate
    for (i=0; i<numnodes; i++)
      {
	ptgrid->GetCoordinateNode(i,level,point);
	
	if (ascii_)
	  (*output) << " " << point[0] << " " << point[1] << " " << point[2] << std::endl;
	else {
	  output->write((char*)&point[0],sizeof(Double));
	  output->write((char*)&point[1],sizeof(Double));
	  output->write((char*)&point[2],sizeof(Double));
	}
	
      }
}
}

void WriteResultsGMV::WriteCells(const Integer alevel) 
{
#ifdef TRACE
  (*trace) << " entering WriteResultsGMV::WriteCells \n";
#endif

  Integer level=alevel;

  // write keyword
  (*output) << "cells   ";

  //
  if (!ptgrid)
    Error("ptgrid is not initialized", __FILE__,__LINE__);

  // read information about number of elements 
  Integer numelem; 
  numelem=ptgrid->GetMaxnumElem(level);

  if (ascii_)
    (*output) << numelem << std::endl;
  else
    output->write((char*)&numelem,sizeof(Integer));

  Vector<Integer> connect;

  Integer dim=ptgrid->GetDim();
  
  Integer i;
  for (i=0; i<numelem; i++)
    {

      ptgrid->GetConnection(connect, i, level);

      if (dim==2)
	{
	  switch (connect.size())
	    {
// 	    case 2:
// 	      if (ascii_)
// 		(*output) << "line 2" << std::endl;
// 	      else {
// 		(*output) << "line    ";
// 		Integer nn=2;
// 		output->write((char*)&nn,sizeof(Integer));
// 	      }
	    case 3: 
	      if (ascii_)
		(*output) << "tri 3" << std::endl;
	      else {
		(*output) << "tri     ";
		Integer nn=3;
		output->write((char*)&nn,sizeof(Integer));
	      }
	      break;
	    case 4:
	      if (ascii_)
		(*output) << "quad 4" << std::endl;
	      else {
		(*output) << "quad    ";
		Integer nn=4;
		output->write((char*)&nn,sizeof(Integer));
	      }
	      break;
	    case 6: 
	      if (ascii_)
		(*output) << "6tri 6" << std::endl;
	      else {
		(*output) << "6tri    ";
		Integer nn=6;
		output->write((char*)&nn,sizeof(Integer));
	      }
	      break;
	    case 8:
	      if (ascii_)
		(*output) << "8quad 8" << std::endl;
	      else {
		(*output) << "8quad   ";
		Integer nn=8;
		output->write((char*)&nn,sizeof(Integer));
	      }
	      break;
	    default:
	      Error("This type of element is not implemented", __FILE__, __LINE__);
	    }
	}
      else
	{
	  switch (connect.size())
	    {
	    case 4:
	      if (ascii_)
		(*output) << "tet 4" << std::endl;
	      else {
		(*output) << "tet     ";
		Integer nn=4;
		output->write((char*)&nn,sizeof(Integer));
	      }
	      break;
	    case 5:
	      if (ascii_)
		(*output) << "ppyrmd5 5" << std::endl;
	      else {
		(*output) << "ppyrmd5 ";
		Integer nn=5;
		output->write((char*)&nn,sizeof(Integer));
	      }
	      break;
	    case 6:
	      if (ascii_)
		(*output) << "pprism6 6" << std::endl;
	      else {
		(*output) << "pprism6 ";
		Integer nn=6;
		output->write((char*)&nn,sizeof(Integer));
	      }
	      break;
	    case 8:
	      if (ascii_)
		(*output) << "hex 8" << std::endl;
	      else 
		{
		  (*output) << "hex     ";
		  Integer nn=8;
		  output->write((char*)&nn,sizeof(Integer));
		}
	      break;
	    default:
	      std::cout << connect.size() << std::endl;
	      Error("This type of element is not implemented", __FILE__, __LINE__);
	    }
	}

      if (ascii_) 
	{
	  Integer j;
	  for (j=0; j< connect.size(); j++)
	    (*output) << " " << connect[j] ;
	  (*output) << std::endl;
	}
      else 
	{
	  Integer * ptcon=connect.get();
	  Integer len=connect.size();
	  output->write((char*)ptcon,len * sizeof(Integer));
	}

    }
}

void WriteResultsGMV::WriteNodeVariable(const Vector<Double> var, const std::string name, const Integer dataType)
{
  (*output) << "variable";
  if (ascii_) (*output) << std::endl;


  if (ascii_)
    (*output) << name << " " << dataType << std::endl;
  else 
    {
      Char * str=new Char[8];  
      to8Char(name,str);
      (*output) << str;
      output->write((Char*)&dataType,sizeof(Integer));
      delete [] str;
    }
  
  if (ascii_) 
    {
      Integer i;
      for (i=0; i<var.size(); i++)
	(*output) << var[i] << " ";
      (*output) << "\n endvars \n" ;
    }
  else 
    {
      Double * ptvar=var.get();
      Integer len=var.size();
      output->write((char*)ptvar,len * sizeof(Double));
      (*output) << "endvars ";
    }

}

  // only for 3D
void WriteResultsGMV::WriteVelocity(const Vector<Double> *  var, const std::string name, const Integer dataType)
{

  (*output) << "velocity";
  if (ascii_)
    (*output) << " " << dataType << std::endl;
  else {
     output->write((char*)&dataType,sizeof(Integer));
  }

  if (ascii_) {
  Integer i,j;
  for (i=0; i<3; i++) {
    for (j=0; j<var[i].size(); j++)
      (*output) << var[i][j] << " ";
    (*output) << std::endl;
  }
  }
  else {
    Integer i;
    for(i=0; i<3; i++) {
     Double * ptvar=var[i].get();
     Integer len=var[i].size();
     output->write((char*)ptvar,len * sizeof(Double));
    }
  }
}



void WriteResultsGMV::WriteGrid(const Integer level)
{
  WriteHeader();
  WriteNodes(level);
  WriteCells(level);
  if (fixedgrid_)
    OpenFile(0);
}





void WriteResultsGMV::WriteNodeSolution(const Array<Double>& sol, const Integer step, const Double time, const std::string title)
{
#ifdef TRACE
 (*trace) << " entering WriteResultsGMV::WriteNodeSolution " << std::endl;
#endif

  Integer i,j;
  if (NeedHistory_)
    for (i=0; i<nodeshist_.size(); i++) {
      {
	if (nodeshist_[i] > sol.size())
	  Error("Nr. of history-node(s) is too high --> not in Solution! ",__FILE__,__LINE__);

	if (sol.dim() * sol.size() <=nodeshist_[i])
        Error("Please, check history-nodes in config-file.",__FILE__,__LINE__);
	//     if (lastsavetime[i] != time )
	if (sol.dim() > 1)	
	  {
	    std::vector<Double> solVec;
	    solVec.resize(sol.dim());
	    for (j=0; j<sol.dim(); j++)
	      solVec[j] = sol[j][(nodeshist_[i]-1)];
	    
	    AddVecInHistory(time, solVec, i);
	  }
	else
	  AddInHistory(time,sol[0][nodeshist_[i]-1],i);
      }
      
    }

  Integer type=1; // 0 - for cell 
                  // 1 - for node
                  // 2 - for face data

  static Integer flagGrid=0;
  if (flagGrid==0 && fixedgrid_) {
    WriteGrid(ptgrid->GetLastLevel());   
    flagGrid++;
  }
  
  if (step!=currstep_)
    {
    OpenFile(step);
    if (!fixedgrid_)
      {
	WriteGrid(ptgrid->GetLastLevel());   
      }
    }
  
    if (fixedgrid_)
      {
	WriteHeader();

	Char * name=new Char[80];
	strcpy(name,"");
	strcat(name,namefile_);
	strcat(name,"_GRID.gmv");
	
	if (ascii_)
	  (*output) << "nodev fromfile \"" << name <<"\""<< std::endl;
	else 
	  (*output) << "nodev   fromfile\"" << name <<"\"";

	
	if (ascii_)
	  (*output) << "cells fromfile \"" << name <<"\""<< std::endl;
	else 
	  (*output) << "cells   fromfile\"" << name <<"\"";

	delete [] name;
      }


      for (i=0; i< sol.dim(); i++)
	{
	  char nrStr[10];
	  sprintf(nrStr,"%i",i+1);
	  std::string sumString = title + nrStr;
	  
	  WriteNodeVariable(sol[i], sumString , type);
	}
      
      
      if (ascii_)
	(*output) << "probtime " << time << std::endl;
      else {
	(*output) << "probtime";
	output->write((char*)&time,sizeof(Double));
      }
      //}
  

  currstep_=step;
}


void WriteResultsGMV::WriteElemSolution(const Array<Double>& data, const Integer step, const Double time, const std::string title)
{
#ifdef TRACE
  (*trace) << " entering WriteResultsGMV::WriteElemSolution " << std::endl;
#endif  

 Integer type=0; // 0 - for cell 
                 // 1 - for node
                 // 2 - for face data
 Integer i = 0;

   if (step!=(currstep_)) {
     Error("You should write solution of this step before printing some cell data",__FILE__,__LINE__);
   }

  for (i=0; i<data.dim(); i++)
    {
      char nrStr[10];
      sprintf(nrStr,"%i",i+1);
      std::string sumString = title + nrStr;
      
      WriteNodeVariable(data[i], sumString , type);
    }
 
  if (ascii_)
    (*output) << "probtime " << time << std::endl;
  else {
    (*output) << "probtime";
    output->write((char*)&time,sizeof(Double));
  }
  
}

// 

// void WriteResultsGMV::WriteDataOnCell(const Vector<Double>&sol,const Integer step, const Double time, const std::string title)
// {
//   Integer type=0; // 0 - for cell 
//                   // 1 - for node
//                   // 2 - for face data

//   if (step!=currstep_) {
//     Error("You should write solution of this step before printing some cell data",__FILE__,__LINE__);
//   }

//   WriteVariable(sol,title,type);

//   if (ascii_)
//     (*output) << "probtime " << time << std::endl;
//   else {
//     (*output) << "probtime";
//     output->write((char*)&time,sizeof(Double));
//   }

// }

// void WriteResultsGMV::WriteVecDataOnCell(const Vector<Double>*vec,const Integer step, const Double time, const std::string title)
// {
//   Integer type=0; // 0 - for cell-centered 
//                   // 1 - for node-centered
//                   // 2 - for face-centered data

//   if (step!=currstep_) {
//     Error("You should write solution of this step before printing some cell data",__FILE__,__LINE__);
//   }

//   WriteVelocity(vec,title,type);

//   if (ascii_)
//     (*output) << "probtime " << time << std::endl;
//   else {
//     (*output) << "probtime";
//     output->write((char*)&time,sizeof(Double));
//   }
  
// }

// void WriteResultsGMV::OpenFile(const Integer num)
// {
//    Char * name=new Char[80];
//    Char * aux=new Char[2];
//    sprintf(aux,"%i",num);

//    strcpy(name,namedir_);
//    strcat(name,"/");
//    strcat(name,namefile_);
//    if (num/10 < 1) strcat(name,".gmv00");
//      else if (num/100 < 1) strcat(name,".gmv0");
//        else strcat(name,".gmv");

//    strcat(name,aux);
//    if (output) {
//      if (ascii_)
//        (*output) << "endgmv " << std::endl; 
//      else
//        (*output) << "endgmv  ";

//      delete output;
//    }

//    // check what kind of data for input
//    std::string typedata;

//    if (ascii_)
//      output=new std::ofstream(name);
//    else
//      output=new std::ofstream(name,std::ofstream::binary);

//    delete [] name;
//    delete [] aux;
// }

void WriteResultsGMV::OpenFile(const Integer num)
{
   Char * name=new Char[80];
   Char * aux=new Char[2];
   sprintf(aux,"%i",num);

   strcpy(name,namedir_);
   strcat(name,"/");
   strcat(name,namefile_);
   if (num==-1)
     {
       strcat(name,"_GRID.gmv");
     }
   else
     {
       if (num/10 < 1) strcat(name,".gmv00");
       else if (num/100 < 1) strcat(name,".gmv0");
       else strcat(name,".gmv");
       strcat(name,aux);
     }
   
   if (output) {
   if (num!=-1)
     {
       if (ascii_)
       (*output) << "endgmv " << std::endl; 
       else
	 (*output) << "endgmv  ";
     }
   
     delete output;
   }



   // check what kind of data for input
   std::string typedata;

   if (ascii_)
     output=new std::ofstream(name);
   else
     output=new std::ofstream(name,std::ofstream::binary);

   delete [] name;
   delete [] aux;
}

void WriteResultsGMV::Init(Grid * aptgrid)
{
ptgrid=aptgrid;
}

void WriteResultsGMV::to8Char(const std::string name, char * result)
{
  std::string aux;
  Integer i;
 
  if (name.size()!= 8) {
      aux="        ";
      for (i=0; i<8; i++)
	aux[i]=name[i];
  }
  else aux=name;

  strcpy(result,aux.c_str());

}


} // end of namespace
