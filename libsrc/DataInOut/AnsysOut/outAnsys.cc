#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include "outAnsys.hh"
#include "DataInOut/conffile.hh"

namespace CoupledField
{

WriteResultsAnsys :: WriteResultsAnsys(const Char * const filename, Boolean withHistory)
: WriteResults(filename,withHistory)
{
#ifdef TRACE
  (*trace) << "entering WriteResultsAnsys :: WriteResultsAnsys" << std::endl;
#endif

 namedir_=new Char[30];

 Char S[50];
 strcpy(namedir_,filename);
 strcat(namedir_,"_ansys");
 sprintf(S,"mkdir -p %s",namedir_);

 system(S);

 ptgrid=NULL;
 currnum_=0;

 out_nodes_=NULL;
 out_cells_=NULL;

 out_res_=NULL;

 OpenFile(0);
}

WriteResultsAnsys ::~WriteResultsAnsys()
{
#ifdef TRACE
  (*trace) << "entering WriteResultsAnsys::~ WriteResultsAnsys" << std::endl;
#endif

 delete [] namedir_;

 delete out_nodes_;
 delete out_cells_;
 if (out_res_) delete out_res_;
}

void WriteResultsAnsys :: WriteNodes(const Integer alevel)
{
  Integer level=alevel;

  // write keyword
  (*out_nodes_) << std::setw(5) << -777 << std::endl;

 //get and write number of nodes on the level
  Integer numnodes=ptgrid->GetMaxnumnodes(level);

  Integer dim=ptgrid->GetDim();

  (*out_nodes_).setf(std::ios::scientific);
  (*out_nodes_).precision(13);

  //get and write coodinates of nodes
  Integer i;
  if (dim==2) {
    // write x,y,z-coordinate
    for (i=0; i<numnodes; i++)
      {
	Point<2> point;
	ptgrid->GetCoordinateNode(i,level,point);

	(*out_nodes_).setf(std::ios::uppercase);

	(*out_nodes_) << std::setw(10) << i+1 << "    ";
	(*out_nodes_) << std::setw(20) << point[0]
		      << std::setw(20) << point[1]
		      << std::setw(20) << .0 << std::endl;
      
         (*out_nodes_) << std::setw(14) << " " << 
	   std::setw(20) << .0 << 
           std::setw(20) << .0 <<
           std::setw(20) << .0 << std::endl; 
      }
  }
  
  if (dim==3) {
    // write x,y,z-coordinate
    for (i=0; i<numnodes; i++)
      {
	Point<3> point;
	ptgrid->GetCoordinateNode(i,level,point);

	(*out_nodes_).setf(std::ios::uppercase);

	(*out_nodes_) << std::setw(10) << i+1 << "    ";
	(*out_nodes_) << std::setw(20) << point[0]
	               << std::setw(20) << point[1]
		      << std::setw(20) << point[2] << std::endl;

	(*out_nodes_) << std::setw(14) << " " << std::setw(20) << .0 <<
	  std::setw(20) << .0 <<
	  std::setw(20) << .0 << std::endl;
      }
}
}

void WriteResultsAnsys::WriteCells(const Integer alevel) 
{
#ifdef TRACE
  (*trace) << " entering WriteResultsAnsys::WriteCells \n";
#endif

  Integer level=alevel;

  if (!ptgrid)
    Error("ptgrid is not initialized", __FILE__,__LINE__);

  // read information about number of elements 
  Integer numelem; 
  numelem=ptgrid->GetMaxnumElem(level);

  Vector<Integer> connect;
  Integer dim=ptgrid->GetDim();

  Integer f=8; // format of output

  Integer i,j;
  for (i=0; i<numelem; i++)
   {
     ptgrid->GetConnection(connect, i, level);

     if (dim==2)
       {
	 switch (connect.GetSize())
	   {
	   case 2: // line
	     (*out_cells_) << std::setw(f) << connect[0] << 
	                      std::setw(f) << connect[1];
	     for (j=0; j<6; j++)
	       (*out_cells_) << std::setw(f) << connect[1];
	     break;
	   case 3: // triangle
	      (*out_cells_) << std::setw(f) << connect[0] 
		            << std::setw(f) << connect[1]
                            << std::setw(f) << connect[2];
	     for (j=0; j<5; j++)
	       (*out_cells_) << std::setw(f) << connect[2];
	     break;
	   case 4: // quad
	     for (j=0; j<4; j++)
	       (*out_cells_) << std::setw(f) << connect[j];
	     for (j=0; j<4; j++)
	       (*out_cells_) << std::setw(f) << connect[3];
	     break;
	   default:
	     Error("This type of element is not implemented", __FILE__, __LINE__);
	   }
       }
     else
{
  switch (connect.GetSize())
    {
    case 4: // tet
      for (j=0; j<2; j++)
	(*out_cells_) << std::setw(f) << connect[j];
       for (j=0; j<2; j++)
      (*out_cells_) << std::setw(f) << connect[2] ;
      for (j=0; j<4; j++)
	(*out_cells_) << std::setw(f) << connect[3];
	     break;
	  break;
    case 8: // brick
      for (j=0; j<8; j++)
	(*out_cells_) << std::setw(f) << connect[j]; 
	  break;
      default:
	std::cout << " connection size for the element is :" << connect.GetSize() << std::endl;
          Error("This type of element is not implemented", __FILE__, __LINE__);
      }
}

     for (j=0; j<4; j++)
     (*out_cells_) << std::setw(f)  << 1 ;
     (*out_cells_) << std::setw(f)  << 0 
		   << std::setw(f)  << i+1 << std::endl;	  
   }
}

void WriteResultsAnsys::WriteVariableNode(const Vector<Double> var)
{
  Integer size=var.GetSize();
  Integer i;
  for (i=0; i< size; i++)
    (*out_res_) << "dnsol,   " << i+1 << ",u,x,"
		 << var[i] << "," << std::endl;
}

void WriteResultsAnsys::WriteVariableCell(const Vector<Double> var, const std::string title)
{

  OpenFileRes(title);

  (*out_res_).setf(std::ios::scientific);
  (*out_res_).setf(std::ios::uppercase);
  (*out_res_).precision(7);

  Integer size=var.GetSize();
  Integer i;
  for (i=0; i< size; i++)
    (*out_res_) << "desol," << std::setw(10) << i+1 << ",all,s,x,"
		<< std::setw(15) << var[i] << std::endl;

}

void WriteResultsAnsys::WriteVelocity(const Vector<Double>* var, const std::string name, const Integer dataType)
{

  Warning("sorry, but in output file ansys-format for vector data is not implemented yet");

}

void WriteResultsAnsys::WriteGrid(const Integer level)
{
 WriteNodes(level);
 WriteCells(level);
}

void WriteResultsAnsys::WriteNodeSolution(const Vector<Double> & sol, const Integer step, const Double time, const std::string title)
{
#ifdef TRACE
 (*trace) << " entering WriteResultsAnsys::WriteSolution " << std::endl;
#endif

 // diese funktion ist nicht verarbeitet

  Integer i;
  if (NeedHistory_)
    for (i=0; i<nodeshist_.size(); i++) {
      if (sol.GetSize()<=nodeshist_[i])
        Error("Please, check history-nodes in config-file.",__FILE__,__LINE__);
      if (lastsavetime[i] != time )
	AddInHistory(time,sol[nodeshist_[i]],i);
    }

  //  Integer type=1; // 0 - for cell 
                      // 1 - for node
                      // 2 - for face data

  if (step!=currnum_) {
    OpenFile(step);

    WriteGrid(ptgrid->GetLastLevel());
  }

  WriteVariableNode(sol);

  currnum_=step;
}

void WriteResultsAnsys::WriteElemSolution(const Vector<Double>&sol,const Integer step, const Double time, const std::string title)
{

  // cut: only for printing my info about mesh
  WriteVariableCell(sol,title);

}

void WriteResultsAnsys::WriteVecDataOnCell(const Vector<Double>*vec,const Integer step, const Double time, const std::string title)
{
 
   Warning("sorry, but in output file ansys-format for vector data is not implemented yet");

}

void WriteResultsAnsys::OpenFile(const Integer num)
{
  currnum_=num;

   Char * name=new Char[30];
   Char * aux=new Char[2];
   sprintf(aux,"%i",num);

   strcpy(name,namedir_);
   strcat(name,"/");
   strcat(name,namefile_);

   Char * nameElmFile=new Char[30];
   strcpy(nameElmFile,name);
   strcat(nameElmFile,".elem");
   strcat(nameElmFile,aux);

   Char * nameNdFile=new Char[30];
   strcpy(nameNdFile,name);
   strcat(nameNdFile,".node");
   strcat(nameNdFile,aux);

   // check what kind of data for input
   std::string typedata;

   out_nodes_=new std::ofstream(nameNdFile);
   out_cells_=new std::ofstream(nameElmFile);

   delete [] name;
   delete [] aux;
   delete [] nameElmFile;
   delete [] nameNdFile;
   
  
}

void WriteResultsAnsys::OpenFileRes(const std::string title)
{
  if (out_res_) delete out_res_;

   Char * name=new Char[50];
   Char * aux=new Char[2];
   sprintf(aux,"%i",currnum_);

   strcpy(name,namedir_);
   strcat(name,"/");
   strcat(name,namefile_);

   Char * nameFile=new Char[30];
   strcpy(nameFile,name);
   strcat(nameFile,".");
   strcat(nameFile,title.c_str());
    if (currnum_/10 < 1) strcat(name,"00");
     else if (currnum_/100 < 1) strcat(name,"0");
   strcat(nameFile,aux);

   out_res_=new std::ofstream(nameFile);

   delete [] name;
   delete [] aux;
   delete [] nameFile;

}

void WriteResultsAnsys::Init(Grid * aptgrid)
{
ptgrid=aptgrid;
}

} // end of namespace
