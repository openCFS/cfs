#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

#include "outGMV.hh"
#include "conffile.hh"

namespace CoupledField
{

WriteResultsGMV :: WriteResultsGMV(const Char * const filename, Boolean withHistory)
: WriteResults(filename,withHistory)
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

 output=NULL;
 OpenFile(0);

 if (!output)
   Error(" File for output results in .gmv-format ", __FILE__, __LINE__);

 ptgrid=NULL;
 currstep_=0;

 ascii_=FALSE;

 std::string typedata;
 if (conf->ifget("format_data_output",typedata)) {
   if (typedata == "binary")
     ascii_=FALSE;
   else ascii_=TRUE;
 }
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

if (dim==2) {
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
      case 2:
	if (ascii_)
	(*output) << "line 2" << std::endl;
	else {
	  (*output) << "line    ";
	   Integer nn=2;
	   output->write((char*)&nn,sizeof(Integer));
	}
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
        case 8:
	  if (ascii_)
	  (*output) << "hex 8" << std::endl;
	  else {
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

if (ascii_) {
     Integer j;
     for (j=0; j< connect.size(); j++)
       (*output) << " " << connect[j] ;

     (*output) << std::endl;
}
else {
  Integer * ptcon=connect.get();
  Integer len=connect.size();
  output->write((char*)ptcon,len * sizeof(Integer));
}

   }
}

void WriteResultsGMV::WriteVariable(const Vector<Double> var, const std::string name, const Integer dataType)
{
  (*output) << "variable";
  if (ascii_) (*output) << std::endl;


  if (ascii_)
  (*output) << name << " " << dataType << std::endl;
  else {
    Char * str=new Char[8];  
    to8Char(name,str);
    (*output) << str;
    output->write((Char*)&dataType,sizeof(Integer));
    delete [] str;
  }
  

  if (ascii_) {
  Integer i;
  for (i=0; i<var.size(); i++)
    (*output) << var[i] << " ";
  (*output) << "\n endvars \n" ;
  }
  else {
     Double * ptvar=var.get();
     Integer len=var.size();
     output->write((char*)ptvar,len * sizeof(Double));
     (*output) << "endvars ";
  }

}

  // only for 3D
void WriteResultsGMV::WriteVelocity(const Vector<Double>* var, const std::string name, const Integer dataType)
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
}

void WriteResultsGMV::WriteSolution(const Vector<Double> & sol, const Integer step, const Double time, const std::string title)
{
#ifdef TRACE
 (*trace) << " entering WriteResultsGMV::WriteSolution " << std::endl;
#endif

  Integer i;
  if (NeedHistory_)
    for (i=0; i<nodeshist_.size(); i++) {
      if (sol.size()<=nodeshist_[i])
        Error("Please, check history-nodes in config-file.",__FILE__,__LINE__);
      if (lastsavetime[i] != time )
	AddInHistory(time,sol[nodeshist_[i]],i);
    }

  Integer type=1; // 0 - for cell 
                  // 1 - for node
                  // 2 - for face data

  if (step!=currstep_) {
    OpenFile(step);

    WriteGrid(ptgrid->GetLastLevel());
  }

  WriteVariable(sol,title,type);

  if (ascii_)
    (*output) << "probtime " << time << std::endl;
  else {
    (*output) << "probtime";
    output->write((char*)&time,sizeof(Double));
  }

  currstep_=step;
}

void WriteResultsGMV::WriteDataOnCell(const Vector<Double>&sol,const Integer step, const Double time, const std::string title)
{
  Integer type=0; // 0 - for cell 
                  // 1 - for node
                  // 2 - for face data

  if (step!=currstep_) {
    Error("You should write solution of this step before printing some cell data",__FILE__,__LINE__);
  }

  WriteVariable(sol,title,type);

  if (ascii_)
    (*output) << "probtime " << time << std::endl;
  else {
    (*output) << "probtime";
    output->write((char*)&time,sizeof(Double));
  }

}

void WriteResultsGMV::WriteVecDataOnCell(const Vector<Double>*vec,const Integer step, const Double time, const std::string title)
{
  Integer type=0; // 0 - for cell-centered 
                  // 1 - for node-centered
                  // 2 - for face-centered data

  if (step!=currstep_) {
    Error("You should write solution of this step before printing some cell data",__FILE__,__LINE__);
  }

  WriteVelocity(vec,title,type);

  if (ascii_)
    (*output) << "probtime " << time << std::endl;
  else {
    (*output) << "probtime";
    output->write((char*)&time,sizeof(Double));
  }
  
}

void WriteResultsGMV::OpenFile(const Integer num)
{
   Char * name=new Char[30];
   Char * aux=new Char[2];
   sprintf(aux,"%i",num);

   strcpy(name,namedir_);
   strcat(name,"/");
   strcat(name,namefile_);
   if (num/10 < 1) strcat(name,".gmv00");
     else if (num/100 < 1) strcat(name,".gmv0");
       else strcat(name,".gmv");

   strcat(name,aux);

   if (output) {
     if (ascii_)
       (*output) << "endgmv " << std::endl; 
     else
       (*output) << "endgmv  ";

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
