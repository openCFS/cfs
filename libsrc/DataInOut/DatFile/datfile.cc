#include <string>
#include <fstream>
#include <stdarg.h>
#include <iostream>
#include <iomanip>

#include "datfile.hh"

namespace CoupledField
{

// ----------------- Constructor of Datfile --------------------------------
   DatFile :: DatFile( const Char * const afilename)
  : FileType(afilename)
{
#ifdef TRACE
  (*trace) << "entering DatFile::DatFile" << std::endl;
#endif

  infile.open(strcat(filename,".dat"));
  if (!infile) {std::cerr << "ERROR(" << __FILE__ << " " << __LINE__ <<
                         ") Can't open " << filename << std::endl;
                exit(1);}

  infile.seekg(0, std::ios::end);
  pos_end=infile.tellg();
}
  
// --------------------------- Deconstructor --------------------------------


DatFile :: ~DatFile()
{
#ifdef TRACE
  (*trace) << "entering DatFile::~DatFile" << std::endl;
#endif

infile.close();

}
// -------------------------ReadTitle---------------------------------------

Boolean DatFile :: ReadTitle(std::string& title)
{
 
#ifdef TRACE
  (*trace) << "entering DatFile:: ReadTitle" << std::endl;
#endif
  infile.seekg(0, std::ios::beg);   
  std::getline(infile,title,'\n');   
  std::getline(infile,title,'\n');         
  Integer isize_title=title.size();
  if (title.empty()) {return FALSE;}
  for (Integer i=0; i < isize_title; i++)
  { if (title[i] != ' ')
         { return TRUE; }
  }
  return FALSE;
}

// ---------------------- Print Title ------------------------------------
 void  DatFile :: PrintTitle (std::ostream * out)
{
 std::string title;
  Boolean indicator_of_title;
  indicator_of_title = ReadTitle(title);
   if (indicator_of_title)
    (*out) << "Title of run process is " << title << std::endl;
  else (*out) << " Title of run process isn't specified " << std::endl;
}  
// ------------------ Read analysis data from General section -------------

void DatFile :: ReadGeneralAnal(Integer * dataGAnal)
{
#ifdef TRACE
  (*trace) << "entering DatFile::ReadGeneralAnal" << std::endl;
#endif 
  std::string::size_type pos=0;
  TakePos("soltype",pos);
  infile.seekg(pos, std::ios::beg);
 
  std::string buf;
  infile >> dataGAnal[soltype] >> dataGAnal[analtype] >> dataGAnal[numnode] >>
            dataGAnal[numgroup] >>  dataGAnal[restart];
  infile.ignore(100, '\n');
  std::getline(infile, buf, '\n');

  infile >> dataGAnal[inactdofs] >> dataGAnal[circuit];
  infile.ignore(100, '\n');
  std::getline(infile, buf, '\n');

  std::getline(infile, buf, '\n');
  pos=buf.find("masstype");
  if (pos == std::string::npos ) {Integer size;
                             size = buf.size() - CountCharInStr(buf,' ');
                             char hold[100];
                             Peel(buf,hold);
                             dataGAnal[deactDf]=TransformInNameDf(hold);
                             std::getline(infile, buf,'\n');
                             }
  infile >> dataGAnal[masstype] >> dataGAnal[damptype] >> dataGAnal[nooptimiz];
}
// ------------------- Read analysis data from General section by choice -----

void DatFile:: ReadGeneralAnalChoice(Integer * dataGAnalCh, 
                enum nameGAn first ...)
{
#ifdef TRACE
  (*trace) << "entering DatFile::ReadGeneralAnalChoice" << std::endl;
#endif
 va_list argumentList;
 va_start(argumentList, first);
 enum nameGAn runner=first;
 Integer dataGAnal[11];
 ReadGeneralAnal(dataGAnal);
 Integer count=0;
 do {
      dataGAnalCh[count]=dataGAnal[runner];
      count++;
      runner=va_arg(argumentList, enum nameGAn);
      }
  while (runner!=endGAnal);
 va_end(argumentList);
}

// ---------------- Read output option: save first, second derivative vp
 void DatFile :: ReadOutputOptions(Boolean & Der1, Boolean & Der2)
{
#ifdef TRACE
  (*trace) << "entering DatFile::ReadOutputOptions" << std::endl;
#endif
  std::string buf;
  std::string::size_type pos=0;
  std::string seek="save first deriv. of nodal degrees of freedom"; 
  Der1=FALSE; Der2=FALSE;
  Char c=' ';

  Integer j;
  for (j=0; j <2; j++) 
{ 
  TakePos(seek,pos); 
  while (c!='\n')
{
  while (c==' ') { infile.get(c);}

  if (c=='\n') break;

  infile.seekg(-1, std::ios::cur);
  infile >> buf;

  if (vp==TransformInNameDf(buf.c_str())) { if (j==0) Der1=TRUE;
                                              else Der2=TRUE; }
  c=' ';
}
   seek="save secnd. deriv. of nodal degrees of freedom"; 
   c=' ';
}
}

// ------------------------- Read output specification data ----------------

 void DatFile :: ReadOutputData(Integer * dataOut)
{
#ifdef TRACE
  (*trace) << "entering DatFile::ReadOutputData" << std::endl;
#endif
  std::string buf;
  std::string::size_type pos=0;
  Integer n;

  TakePos("formsave",pos);
  infile.seekg(pos, std::ios::beg);
  infile >> dataOut[formsave] >> dataOut[nsavenode] >> dataOut[savemax];
  infile.ignore(100, '\n');
  
  Integer k=3;
 for (Integer i=0; i < 3; i++)
{
  std::getline(infile, buf, '\n');
  pos=infile.tellg();
  n=CountElemInString(pos);
  infile.seekg(pos, std::ios::beg);
  if (n==0) {
           dataOut[k+i]=non;
           infile.ignore(100,'\n'); }
  if (n==1) { 
                 std::string help_str;
                 infile >> help_str;
                 dataOut[k+i]=TransformInNameDf(help_str.c_str());
                infile.ignore(100,'\n');
               }
  if (n>1)  { 
             dataOut[k+i]=Cryptogr(n); 
             infile.ignore(100,'\n');  
             } 
}
}

// -------------- Read output specification data for dofs in array ----------

 void DatFile :: ReadOutputDataArr(Integer * dataOutArr, Integer sizeArr,
                                   std::string & seekStr)
{  
#ifdef TRACE
  (*trace) << "entering DatFile::ReadOutputDataArr" << std::endl;
#endif
  std::string buf;
  std::string::size_type pos=0;
  TakePos(seekStr,pos);
  infile.seekg(pos, std::ios::beg);
  Char c=' ';
  for (Integer i=0; i < sizeArr; i++) {
                                        while (c==' ') infile.get(c);
                                        infile.seekg(-1, std::ios::cur);
                                        infile >> buf;
                                        c=' ';
                               dataOutArr[i]=TransformInNameDf(buf.c_str());
                                      }
}  

// ------------------------ Read Nodal Coordinates 3D -------------------------

void DatFile :: ReadCoordinate(Point3D * const InitNodalCo,                                                     const Integer maxnumNod)
{
#ifdef TRACE
  (*trace) << "entering DatFile::ReadCoordinate 3D" << std::endl;
#endif
 Integer ii;
 std::string::size_type pos=0;
 TakePos("coordinate",pos);
 infile.seekg(pos,std::ios::beg);
 for (Integer i=0; i < maxnumNod; i++)
  {
   infile >> ii >>InitNodalCo[i].x >>  InitNodalCo[i].y >>  InitNodalCo[i].z; 
  }
}

// ------------------------ Read Nodal Coordinates 2D------------------------
 
void DatFile :: ReadCoordinate(Point2D * const InitNodalCo,
                                const Integer maxnumNod)
{
#ifdef TRACE
  (*trace) << "entering DatFile::ReadCoordinate 2D" << std::endl;
#endif
 Integer ii;
 Double dummy;
 std::string::size_type pos=0;
 TakePos("coordinate",pos);
 infile.seekg(pos,std::ios::beg);
 for (Integer i=0; i < maxnumNod; i++)
  {
   infile >> ii >> dummy >> InitNodalCo[i].x >>  InitNodalCo[i].y ;
  }
}

// ------------------Read General Info about element -------------------------
void DatFile :: ReadGeneralElem(Integer * dataGElem, const Integer num)
{
#ifdef TRACE
  (*trace) << "entering DatFile::ReadGeneralElem" << std::endl;
#endif 

  std::string::size_type pos=0;

  Integer i;
  for (i=0; i<=num; i++)
     TakePos("group definition",pos); 

  infile.seekg(pos, std::ios::beg);
  std::string buf;
  std::getline(infile,buf,'\n' );
 
  infile >> dataGElem[numelem] >> dataGElem[ielemtyp] >> dataGElem[isubtype] >>
       dataGElem[ielemsave] >>  dataGElem[maxnode] >>  dataGElem[nonlinear] >>
       dataGElem[form1] >>  dataGElem[form2];
 
}
// ------------------- Read general data about elements section by choice ----
 
void DatFile:: ReadGeneralElemChoice(const Integer num, Integer * dataGElemCh,
                enum nameGElem first ...)
{
#ifdef TRACE
  (*trace) << "entering DatFile::ReadGeneralElemCh" << std::endl;
#endif

 va_list argumentList;
 va_start(argumentList, first);
 enum nameGElem runner=first;
 Integer dataGElem[8];
 ReadGeneralElem(dataGElem, num);
 Integer count=0;
 do {
      dataGElemCh[count]=dataGElem[runner];
      count++;
      runner=va_arg(argumentList, enum nameGElem);
      }
  while (runner!=endGElem);
 va_end(argumentList);
}

// ------------------ Read number of time functions ---------------------------

void DatFile :: ReadNumTimeFunc(Integer & maxnumTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering DatFile:: ReadNumTimeFunc" << std::endl;
#endif

   std::string::size_type pos=0;
   TakePos("timefunctions",pos);
   if (pos==pos_end) {std::cerr << "ERROR time function is absent in your dat file. Don't create object TFuncDat." << std::endl; 
                      exit(1);}
   infile.seekg(pos,std::ios::beg); 
   Char psign;
   infile >> psign; 
   if (psign == '#') maxnumTimeFunc=1;
   else  { 
           infile.seekg(pos,std::ios::beg);
           infile >> maxnumTimeFunc ; 
         }
}
// ----------------- Read Max Value for each time function ----------------

 void DatFile:: ReadInfoTimeFunc ( Integer * maxvalTimeFunc,                                                 const Integer maxnumTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering DatFile:: ReadInfoTimeFunc" << std::endl;
#endif
 Integer num;
 std::string buffer;
 std::string::size_type pos=0;
 TakePos("timefunctions",pos);
 if (pos==pos_end) {std::cerr << "ERROR" << std::endl; exit(1);}
 infile.seekg(pos, std::ios::beg);
 for (Integer i=0; i < maxnumTimeFunc; i++)
 {
       pos=std::string::npos;
       while (pos == std::string::npos) 
       {
        std::getline(infile,buffer,'\n');
        pos=buffer.find("nfunc");
       }
   infile >> num >> maxvalTimeFunc[i];
 } 
}

// ------------------ Read time function --------------------------------

void DatFile :: ReadTimeFunc(Double * const timeTimeFunc,
                              Double * const valTimeFunc,
                              const Integer numTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering DatFile::ReadTimeFunc" << std::endl;
#endif 
  Integer num=0, nval;
  std::string buf;
  infile.seekg(0,std::ios::beg);
  std::string::size_type pos=std::string::npos;
  while (num != numTimeFunc)
{ while (pos == std::string::npos) 
    {
     std::getline(infile, buf, '\n');
     pos=buf.find("nfunc");
     }
    infile >> num >> nval;
  }
  for (int i=0; i< nval; i++)
  {
  infile >> timeTimeFunc[i] >> valTimeFunc[i];
  }
}

// -------------------pre Read transient analysis data -----------------------
void DatFile :: preReadTransAnal(Integer & soltype, Integer & statickey)
{
#ifdef TRACE
  (*trace) << "entering DatFile::preReadTransAnal" << std::endl;
#endif
  std::string::size_type pos=0;
  TakePos("ibemcheck",pos);
  infile.seekg(pos, std::ios::beg);
  infile >> soltype;

  TakePos("statkey",pos);
  infile.seekg(pos, std::ios::beg);

  infile >> statickey;
}

// ------------------- Read transient analysis data for direct solver ---------

void DatFile :: ReadTransAnalDir(Integer * dataTAnalDir, Double & bemtolDir)
{
#ifdef TRACE
  (*trace) << "entering DatFile::ReadTransAnalDir" << std::endl;
#endif
  std::string::size_type pos=0;
  TakePos("ibemcheck",pos);
  infile.seekg(pos, std::ios::beg);

 std::string buf;
 Integer dum;
 infile >> dum >> dataTAnalDir[impexpDir] >> dataTAnalDir[ncorrectDir] >>
                  dataTAnalDir[bemcheckDir];
 infile.ignore(100, '\n');
 std::getline(infile, buf,'\n');

 infile >> bemtolDir;
}

// ------------------- Read transient analysis data for iterative solver ------

void  DatFile :: ReadTransAnalIter(Integer * dataTAnalIt, Double & bemtolIt,
                                   Double & convtolIt)
{
#ifdef TRACE
  (*trace) << "entering DatFile::ReadTransAnalDir" << std::endl;
#endif
 std::string buf;
 Integer dum;
 std::string::size_type pos=0;
 TakePos("ibemcheck",pos);
 infile.seekg(pos, std::ios::beg);

 infile >> dum >> dataTAnalIt[itersoltype] >> dataTAnalIt[ncorrectIt] >>
                  dataTAnalIt[bemcheckIt] >>  dataTAnalIt[iprecon] >>
                   dataTAnalIt[iterlog] >>  dataTAnalIt[vecgmres];
 infile.ignore(100, '\n');
 std::getline(infile, buf,'\n');

 infile >> bemtolIt >> convtolIt;
}

// --------------------- Read transient analysis data , nonlinear parameters--

void DatFile :: ReadTransAnalNonl(Integer * dataTAnalNonl,
                                   Double * dataTAnalNonlTol)
{
#ifdef TRACE
  (*trace) << "entering DatFile::ReadTransAnalNonl" << std::endl;
#endif

 infile.seekg(0, std::ios::beg);
 std::string buf;
 std::string::size_type pos=0;
 TakePos("nonlinear paramters",pos);
 infile.seekg(pos, std::ios::beg);
 std::getline(infile, buf,'\n');
 infile >> dataTAnalNonl[nlinalg] >> dataTAnalNonl[nlcrit] >>
           dataTAnalNonl[nform] >>  dataTAnalNonl[nlmaxit] >>
           dataTAnalNonl[updgeo];
 
 infile.ignore(100, '\n');
 std::getline(infile, buf,'\n');
 infile >> dataTAnalNonlTol[0] >>  dataTAnalNonlTol[1] >>
           dataTAnalNonlTol[2] >>  dataTAnalNonlTol[3];
}

// ---------------- Read transient analysis data for transient behavior

void DatFile :: ReadTransAnalTran(Integer * dataTAnalTran,
                                   Double * dataTAnalTranReal)
{
#ifdef TRACE
  (*trace) << "entering DatFile::ReadTransAnalTran" << std::endl;
#endif
 
  std::string::size_type pos=0;
  TakePos("numsteps",pos);
  infile.seekg(pos, std::ios::beg);

  infile >> dataTAnalTran[numstep] >> dataTAnalTran[isavebegin] >>
           dataTAnalTran[isaveend] >>  dataTAnalTran[isaveinc] >>
           dataTAnalTranReal[deltatTr];
  std::string buf;
  infile.ignore(100, '\n');
  std::getline(infile, buf,'\n');

 infile >> dataTAnalTranReal[alpha] >>  dataTAnalTranReal[beta] >>
           dataTAnalTranReal[gamma_h] >>  dataTAnalTranReal[gamma_p];
}

// ----------------- Read transient analysis data for static behavior 

void DatFile :: ReadTransAnalStat(Integer & numstepStat,
                                   Double & deltatStat )
{
#ifdef TRACE
  (*trace) << "entering DatFile::ReadTransAnalStat" << std::endl;
#endif
   std::string::size_type pos=0;
   TakePos("numsteps",pos);
   infile.seekg(pos, std::ios::beg);
   infile >> numstepStat >> deltatStat; 
}

// ----------------- Read eigenvalue analisys data ----------------------

void DatFile :: ReadEigenvalAnal(Integer * dataEAnal,
                                 Double * dataEAnalReal )
{
#ifdef TRACE
  (*trace) << "entering DatFile::ReadEigenvalAnal" << std::endl;
#endif
 infile.seekg(0, std::ios::beg);
 std::string buf;
 std::string::size_type pos=std::string::npos;
 while (pos == std::string::npos)
 {
  std::getline(infile,buf,'\n');
  pos=buf.find("neigval");       // ?!??????????????????????????????
 }
 infile >> dataEAnal[method] >>  dataEAnal[neigval] >> dataEAnal[maxiter];
 infile.ignore(100, '\n');
 std::getline(infile, buf,'\n');
 infile >> dataEAnalReal[flower] >> dataEAnalReal[fupper] >>
           dataEAnalReal[fshift] >>  dataEAnalReal[convtol];
}

// ----------------- Read harmonic analysis data ------------------------

void DatFile :: ReadHarmAnal(Integer * dataHarm,
                                 Double * dataHarmReal )
{
#ifdef TRACE
  (*trace) << "entering DatFile::ReadHarmAnal" << std::endl;
#endif

  std::string::size_type pos=0;
  TakePos("lower",pos);        // ???????????????????????????????
  infile.seekg(pos, std::ios::beg);

   infile >> dataHarmReal[cwlower] >> dataHarmReal[cwupper] >>
           dataHarm[cwnumf] >> dataHarm[cwspace] >> dataHarm[cwout];
}


// ----------------- Read general info about boundary conditions -------------

void DatFile :: ReadGeneralBound (Integer * dataGBound, Double * dataGBoundD)
{
#ifdef TRACE
  (*trace) << "entering DatFile::ReadBoundGen" << std::endl;
#endif 
  std::string buf;
  std::string::size_type pos=0;
  TakePos("numdofs",pos);
  infile.seekg(pos, std::ios::beg);

  infile >> dataGBound[numdofs] >> dataGBound[numconstr] >> 
            dataGBound[numrestr] >> dataGBound[numloads] >>
            dataGBound[resistors] >> dataGBound[numspring] >>
            dataGBound[bembdry] >> dataGBound[numflux];
  infile.ignore(100, '\n');
  std::getline(infile, buf, '\n');

  infile >> dataGBound[numrad];
  
   TakePos("numpress",pos);
   infile.seekg(pos, std::ios::beg);

   infile >> dataGBound[numpress] >> dataGBound[ncurrdens];
   infile.ignore(100, '\n');
   std::getline(infile, buf, '\n');
   infile >> dataGBoundD[regler_offset] >>  dataGBoundD[regler_tol] >> 
             dataGBound[regl_itmx] >> dataGBound[reglr_use];
}
// ---------------- Read number of nodes at which there is Dirichlet BC
void DatFile::ReadNumNodesforDirichletBC(Integer & n)
{
  #ifdef TRACE
  (*trace) << "entering DatFile::ReadNumNodesforDirichletBC" << std::endl;
#endif
  std::string::size_type pos=0;
  TakePos("numrestr",pos);
  infile.seekg(pos, std::ios::beg);

  Integer dummy;
  infile >> dummy >> dummy >> n;

}

// ---------------- Read array of nodes at which there is Dirichlet BC

void DatFile::ReadDirichletBC(Integer * nodes)
{
  #ifdef TRACE
  (*trace) << "entering DatFile::ReadDirichletBC" << std::endl;
#endif
  Integer n;
  ReadNumNodesforDirichletBC(n);     

  std::string::size_type pos=0;
  TakePos("restraints",pos);
  infile.seekg(pos, std::ios::beg);

  Integer i;
  for (i=0; i <n; i++)
    {
         infile >> nodes[i];
         infile.ignore(100,'\n');
    }
}

void DatFile::ReadDirichletBC(Vector<Integer> & nodes)
{
  #ifdef TRACE
  (*trace) << "entering DatFile::ReadDirichletBC" << std::endl;
#endif
  Integer n;

  std::string::size_type pos=0;
  TakePos("numrestr",pos);
  infile.seekg(pos, std::ios::beg);
 
  Integer dummy;
  infile >> dummy >> dummy >> n;

  nodes.Allocate(n);
 
  TakePos("restraints",pos);
  infile.seekg(pos, std::ios::beg);

//  std::string buf;
//  std::getline(infile, buf, '\n');
 
  Integer i;
  for (i=0; i <n; i++)
    {
         infile >> nodes[i];
         infile.ignore(100,'\n');
    }
}

// ---------------- Read boundary conditions: nodal dof definitions ----------
void DatFile:: ReadBoundDof(Integer ** dataBDof, Integer numberdofs, 
                           Integer maxrecord)
{
#ifdef TRACE
  (*trace) << "entering DatFile::ReadBoundDof" << std::endl;
#endif
  std::string::size_type pos=0;
  TakePos("dof-settings",pos);
  infile.seekg(pos, std::ios::beg); 
  std::string str;
  for (Integer i=0; i < numberdofs; i++) {
    infile >> dataBDof[i][0];
     for (Integer ii=1; ii < maxrecord; ii++) {
     infile >> str;
    if (str.empty() || str[0]=='#') dataBDof[i][ii]=non;
    else
    dataBDof[i][ii]=TransformInNameDf(str.c_str());
    }
    infile.ignore(100,'\n');  
   }
}

// ---------------- Read boundary conditions: loads ----------
void DatFile:: ReadBoundLoads(Integer ** dataBLoads, Integer number,
                           Double * factor)
{
#ifdef TRACE
  (*trace) << "entering DatFile::ReadBoundLoads" << std::endl;
#endif
  std::string::size_type pos=0;
  TakePos(" loads",pos);
  infile.seekg(pos, std::ios::beg);
  std::string str;
  for (Integer i=0; i < number; i++) {
    infile >> dataBLoads[i][0];
     infile >> str;
    dataBLoads[i][1]=TransformInNameDf(str.c_str());
     infile >> dataBLoads[i][2];
     infile >> factor[i];
    infile.ignore(100,'\n');
   }
}
// ---------------- Read boundary conditions: constraint ----------
void DatFile:: ReadBoundConstr(Integer ** dataBDof, Integer numberdofs,
                           Double * factor)
{
#ifdef TRACE
  (*trace) << "entering DatFile::ReadBoundConstr" << std::endl;
#endif
  std::string::size_type pos=0;
  TakePos("constraints",pos);
  infile.seekg(pos, std::ios::beg);
  std::string str;
  for (Integer i=0; i < numberdofs; i++) {
    infile >> dataBDof[i][0]; 
     infile >> str;
    dataBDof[i][1]=TransformInNameDf(str.c_str());
     infile >> dataBDof[i][2];
     infile >> str;
    dataBDof[i][3]=TransformInNameDf(str.c_str());
     infile >> factor[i]; 
    infile.ignore(100,'\n');
   }
}

// ---------------- Read boundary conditions: restraint conditions  ----------
void DatFile:: ReadBoundRestr(Integer ** dataBRestr, Integer numberRestr, 
                              Double * factorRestr)
{
#ifdef TRACE
  (*trace) << "entering DatFile::ReadBoundRestr" << std::endl;
#endif
  std::string::size_type pos=0;
  TakePos("restraints",pos);
  infile.seekg(pos, std::ios::beg);
  std::string str;
  for (Integer i=0; i < numberRestr; i++) {
    infile >> dataBRestr[i][0] >> str >> dataBRestr[i][2] >> factorRestr[i];
    dataBRestr[i][1]=TransformInNameDf(str.c_str());
    infile.ignore(100,'\n');
   }
}

// --------------- Read parameters about boundary condition by choice -------- 
void DatFile:: ReadGeneralBoundChoice(Integer * dataGBoundCh,
                enum nameBound first ...)
{
#ifdef TRACE
  (*trace) << "entering DatFile::ReadBoundGenChoice" << std::endl;
#endif
 va_list argumentList;
 va_start(argumentList, first);
 enum nameBound runner=first;
 Integer dataGBound[13];
 Double dataGBoundD[2];
 ReadGeneralBound(dataGBound,dataGBoundD);
 Integer count=0;
 do {
      dataGBoundCh[count]=dataGBound[runner];
      count++;
      runner=va_arg(argumentList, enum nameBound);
      }
  while (runner!=endBound);
 va_end(argumentList);
}

// ----------------- Read material number and calculation expession for elem

void DatFile:: ReadElemMatCalc(Integer & matnumber, std::string & calc_expr)
{
#ifdef TRACE
  (*trace) << "entering DatFile::ReadElemMatCalc" << std::endl;
#endif
  std::string::size_type pos=0;
  TakePos("matnum",pos);
  infile.seekg(pos, std::ios::beg);

  infile >> matnumber;
  infile >> calc_expr;
  if (calc_expr.empty() || calc_expr[0] == '#') calc_expr="non";

}

// ----------------- Read element thickness -------------------------------

void DatFile:: ReadThickness(Double & thickness)
{
#ifdef TRACE
  (*trace) << "entering DatFile::ReadThickness" << std::endl;
#endif
  std::string::size_type pos=0;
  TakePos("thickness",pos);
  infile.seekg(pos, std::ios::beg);
 
  infile >> thickness;
}

// ------------------- Read connection for element in 3d ------------------
void DatFile:: ReadElemRecord3d(Integer elemnum, Integer * connect,
                     Integer maxnode, const Integer numelemgr)
{
#ifdef TRACE
  (*trace) << "entering DatFile::ReadElemRecord3d" << std::endl;
#endif
  Integer maxrecords,i;
  if (maxnode<9) {maxrecords=1;}
  else if (maxnode < 17) {maxrecords=2;}
  else {maxrecords=3;}

  std::string::size_type pos=0;
  for (i=0; i<=numelemgr; i++)
     TakePos("element definition",pos);
  infile.seekg(pos, std::ios::beg);
  
  Integer help;
  std::string buf;
  infile >> help;
  infile.ignore(100,'\n');
  while (help!=elemnum)
  { for (i=0; i< maxrecords; i++)
          std::getline(infile,buf,'\n');
    infile >> help;
    infile.ignore(100,'\n');
  }

  switch(maxrecords)
  {
  case 1: for (i=0; i<maxnode; i++) infile >> connect[i];
          break;
  case 2: for (i=0; i<8; i++) infile >> connect[i];
          infile.ignore(100,'\n');
          for (i=8; i<maxnode; i++) infile >> connect[i];
          break;
  case 3: for (i=0; i<8; i++) infile >> connect[i];
          infile.ignore(100,'\n'); 
          for (i=8; i<16; i++) infile >> connect[i];
          infile.ignore(100,'\n'); 
          for (i=16; i<maxnode; i++) infile >> connect[i];
          break;
  }

} // end of func ReadElemRecord3d

// ------------------- Read connection for element in 3d ------------------
void DatFile:: ReadElem(const Integer maxelem, Integer ** connect,
                     const Integer maxnode, const Integer numelemgr)
{
#ifdef TRACE
  (*trace) << "entering DatFile::ReadElem" << std::endl;
#endif
  Integer maxrecords,i;
  std::string buf;
  if (maxnode<9) {maxrecords=1;}
  else if (maxnode < 17) {maxrecords=2;}
  else {maxrecords=3;}
 
  std::string::size_type pos=0;
  for (i=0; i<=numelemgr; i++)
     TakePos("element definition",pos);
  infile.seekg(pos, std::ios::beg);

  Integer j;
  std::getline(infile, buf, '\n'); 
  for (j=0; j<maxelem; j++) { 
  switch(maxrecords)
  {
  case 1: for (i=0; i<maxnode; i++) infile >> connect[j][i];
          break;
  case 2: for (i=0; i<8; i++) infile >> connect[j][i];
          infile.ignore(100,'\n');
          for (i=8; i<maxnode; i++) infile >> connect[j][i];
          break;
  case 3: for (i=0; i<8; i++) infile >> connect[j][i];
          infile.ignore(100,'\n');
          for (i=8; i<16; i++) infile >> connect[j][i];
          infile.ignore(100,'\n');
          for (i=16; i<maxnode; i++) infile >> connect[j][i];
          break;
  }
  infile.ignore(100,'\n');
  std::getline(infile, buf, '\n');
  }
 
} // end of func ReadElemRecord3d

// --------------- Read connection for element in gridhierarchy----------
void DatFile:: ReadElemConnectionGH(const Integer maxelem, Integer * connect,
            const Integer maxnode, const Integer numelemgr)
/// number of groups is not implemented
{
#ifdef TRACE
  (*trace) << "entering DatFile::ReadElemConnectionGH" << std::endl;
#endif
  Integer maxrecords,i;
  std::string buf;
  if (maxnode<9) {maxrecords=1;}
  else if (maxnode < 17) {maxrecords=2;}
  else {maxrecords=3;}
 
  std::string::size_type pos=0;
  for (i=0; i<=numelemgr; i++)
     TakePos("element definition",pos);
  infile.seekg(pos, std::ios::beg);
 
  Integer j, counter=0;
  std::getline(infile, buf, '\n');
  for (j=0; j<maxelem; j++)
 {
  switch(maxrecords)
  {
  case 1: for (i=0; i<maxnode; i++, counter++)  infile >> connect[counter]; 
          break;
  case 2: for (i=0; i<8; i++, counter ++) infile >> connect[counter];
          infile.ignore(100,'\n');
          for (i=8; i<maxnode; i++, counter++) infile >> connect[counter];
          break;
  case 3: for (i=0; i<8; i++, counter++) infile >> connect[counter];
          infile.ignore(100,'\n');
          for (i=8; i<16; i++, counter++) infile >> connect[counter];
          infile.ignore(100,'\n');
          for (i=16; i<maxnode; i++, counter++) infile >> connect[counter];
          break;
  }
  infile.ignore(100,'\n');
  std::getline(infile, buf, '\n');
  }
} // end of func ReadElemRecord3d

// ------------------- Read connection for element  ------------------
void DatFile:: ReadElemConnect(Integer elemnum, Integer * connect,
                     Integer maxnode)
{
#ifdef TRACE
  (*trace) << "entering DatFile::ReadElemConnect" << std::endl;
#endif
  Integer maxrecords,i;
  if (maxnode<9) {maxrecords=1;}
  else if (maxnode < 17) {maxrecords=2;}
  else {maxrecords=3;}
 
  std::string::size_type pos=0;
  TakePos("element definition",pos);
  infile.seekg(pos, std::ios::beg);
 
  Integer help;
  std::string buf;
  infile >> help;
  infile.ignore(100,'\n');
  while (help!=elemnum)
  { for (i=0; i< maxrecords; i++)
          std::getline(infile,buf,'\n');
    infile >> help;
    infile.ignore(100,'\n');
  }
 
  switch(maxrecords)
  {
  case 1: for (i=0; i<maxnode; i++) infile >> connect[i];
          break;
  case 2: for (i=0; i<8; i++) infile >> connect[i];
          infile.ignore(100,'\n');
          for (i=8; i<maxnode; i++) infile >> connect[i];
          break;
  case 3: for (i=0; i<8; i++) infile >> connect[i];
          infile.ignore(100,'\n');
          for (i=8; i<16; i++) infile >> connect[i];
          infile.ignore(100,'\n');
          for (i=16; i<maxnode; i++) infile >> connect[i];
          break;
  }
 
} // end of func ReadElemConnect
// --------------------- Read number of records for dof-settings ----------

void DatFile::ReadMaxRecord(Integer & maxrecord)
{
#ifdef TRACE 
 (*trace)<<" entering DatFile::ReadMaxRecord" << std::endl;
#endif  
  std::string::size_type pos=0;
  TakePos("dof-settings",pos);
  infile.seekg(pos, std::ios::beg);
  std::string buf;
  std::getline(infile,buf,'\n');
  pos=infile.tellg();
  maxrecord=CountElemInString(pos);
}

// ------------------ For driver.cc --------------------------------
void DatFile::ReadIntegrationParam(Double & alpha, Double & beta, Double & gamma)
{
#ifdef TRACE 
 (*trace)<<" entering DatFile::ReadIntegrationParam" << std::endl;
#endif 

  std::string::size_type pos=0;
  TakePos("alpha",pos);
  infile.seekg(pos, std::ios::beg);

  infile >> alpha >> beta >> gamma;
}

void DatFile :: ReadNumStepsAndTimeSteps(Integer & numsteps, Double & dt)
{
#ifdef TRACE 
 (*trace)<<" entering DatFile::ReadTimeSteps" << std::endl;
#endif 
  std::string::size_type pos=0;
  TakePos("delta-t",pos, "deltat");
  infile.seekg(pos, std::ios::beg);

  Integer dummy;
  infile >> numsteps >>  dummy >>  dummy >>  dummy >> dt;
}


// --------------------- Take position in file -----------------------------

void DatFile::TakePos(const std::string seekexp, std::string::size_type & pos, const std::string reservexp)
{ 
  infile.seekg(pos, std::ios::beg);
  std::string buf, buf1;
  std::string::size_type pos1=pos;
  pos=std::string::npos;

  while ( pos == std::string::npos & !infile.eof() )
  { 
    std::getline(infile, buf, '\n');
    pos=buf.find(seekexp);
  }

  pos=infile.tellg();

  if (pos==pos_end & reservexp!="") 
  {
      infile.seekg(0, std::ios::beg);
      pos=std::string::npos;



      while ( pos == std::string::npos & pos<=pos_end )
      { 
        std::getline(infile, buf1, '\n');
        pos=buf.find(reservexp);
      }


      pos=infile.tellg();
  }

  if (pos==pos_end) {std::cerr << "ERROR: (" << __FILE__ <<" "<< __LINE__ 
               << ") Cannot find string: " << seekexp ;
                    if (reservexp!="") std::cerr <<" and " << reservexp ;
          std::cerr << " in your dat file.\n\t\t Please, change your dat file."<< std::endl;
                      exit(1);}
}

// ------ Take position in file with saving in buf previous std::string -----
 
void DatFile::TakePos(const std::string seekexp, std::string::size_type & pos, std::string & buf_prev, const std::string reservexp)
{
  infile.seekg(pos, std::ios::beg);
  std::string buf;
  std::string::size_type pos1=pos;
  pos=std::string::npos;
 
  while ( pos == std::string::npos & !infile.eof() )
  { std::getline(infile, buf, '\n');
    pos=buf.find(seekexp);
    buf_prev=buf;
  }
  pos=infile.tellg();
  
  if (pos==pos_end & reservexp!="")
  {
      infile.seekg(pos1, std::ios::beg);
      pos=std::string::npos;
 
      while ( pos == std::string::npos & !infile.eof() )
      { std::getline(infile, buf, '\n');
        pos=buf.find(reservexp);
        buf_prev=buf;
      }
      pos=infile.tellg();
  }

  if (pos==pos_end) {std::cerr << "ERROR: (" << __FILE__ <<" "<< __LINE__
                     << ") Cannot find string: " << seekexp ;
                    if (reservexp!="") std::cerr<<" and "<<reservexp; 
     std::cerr <<" in your dat file.\n\t\t Please, change your dat file."<< std::endl;
                      exit(1);}
}

// ---------------------- Return TRUE if there is seekexp in file, else FALSE
Boolean DatFile::IsThere(const std::string seekexp) 
{
  std::string buf;
  std::string::size_type pos=std::string::npos;

  infile.seekg(0, std::ios::beg);
   while ( pos == std::string::npos & !infile.eof() )
  { std::getline(infile, buf, '\n');
    pos=buf.find(seekexp);
  }
   pos=infile.tellg();
   if (pos==pos_end) return FALSE;
   else return TRUE;
} 
  
/// Print Point3D
//void DatFile::PrintPoint(const Point3D point, std::ostream * out) const 
//{
// (*out)<<".\t " << point.x << "\t " <<  point.y << "\t " <<  point.z ;
//}
/// Print Point2D
//void DatFile::PrintPoint(const Point2D point, std::ostream * out) const
//{
// (*out)<<" \t " << point.x << "\t " << point.y ;
//}

// ---------------------- Print element of enum nameDf for Int argv----------
  std::string DatFile :: PrintnameDf( Integer el)
{ std::string result; 
 switch (el)
  { case 0: result="non"; break;
    case 1: result="dx"; break;
    case 2: result="dy"; break;
    case 3: result="dz"; break;
    case 4: result="ep"; break;
    case 5: result="vp"; break;
    case 6: result="chg"; break; 
    case 7: result="ax"; break;
    case 8: result="ay"; break;
    case 9: result="az"; break;
    case 10: result="esp"; break;
    case 11: result="anx"; break;
    case 12: result="any"; break;
    case 13: result="anz"; break;
    case 14: result="tp"; break;
    case 15: result="tnp"; break;
    case 16: result="enp"; break;
    case 17: result="vnp"; break;
    case 18: result="mue"; break;
    case 19: result="sig"; break;
  }
  return result;
}
 
//  ---------------------- Print element of enum nameDf for enum argv--------
  std::string DatFile :: PrintnameDf( enum nameDf el)
{ std::string result;
 switch (el)
  { case dx: result="dx"; break;
    case dy: result="dy"; break;
    case dz: result="dz"; break;
    case ep: result="ep"; break;
    case vp: result="vp"; break;
    case chg: result="chg"; break;
    case ax: result="ax"; break;
    case ay: result="ay"; break;
    case az: result="az"; break;
    case esp: result="esp"; break;
    case anx: result="anx"; break;
    case any: result="any"; break;
    case anz: result="anz"; break;
    case tp: result="tp"; break;
    case tnp: result="tnp"; break;
    case enp: result="enp"; break;
    case vnp: result="vnp"; break;
    case mue: result="mue"; break;
    case sig: result="sig"; break;
    case non: result="non";
  }
  return result;
}

// -------------- Transform std::string in element of type nameDf ---------------

   enum DatFile::nameDf DatFile :: TransformInNameDf(const Char * el)
{ enum DatFile::nameDf result;
  if (!strcmp("dx", el)) result=dx;
  if (!strcmp("dy", el)) result=dy;
  if (!strcmp("dz", el)) result=dz;
  if (!strcmp("ep", el)) result=ep;  
  if (!strcmp("vp", el)) result=vp;
  if (!strcmp("chg", el)) result=chg;
  if (!strcmp("ax", el)) result=ax;
  if (!strcmp("ay", el)) result=ay;
  if (!strcmp("az", el)) result=az;
  if (!strcmp("esp", el)) result=esp;
  if (!strcmp("anx", el)) result=anx;
  if (!strcmp("any", el)) result=any;
  if (!strcmp("anz", el)) result=anz;
  if (!strcmp("tp", el)) result=tp;
  if (!strcmp("tnp", el)) result=tnp;
  if (!strcmp("enp", el)) result=enp;
  if (!strcmp("vnp", el)) result=vnp;
  if (!strcmp("mue", el)) result=mue;
  if (!strcmp("sig", el)) result=sig;
  if (!strcmp("non", el)) result=non;
  return result;
}

// -------------------- Preprocess for transform std::string in enum ---------

  void DatFile::Peel(std::string & buf, Char * init)
{ 
  Integer k=0, size;
  size = buf.size() - CountCharInStr(buf,' '); 
  Char * result = new char[size+1]; 
  for (Integer i=0; i < strlen(buf.c_str()); i++)
  { if (buf.c_str()[i]!=' ') {result[k]=buf.c_str()[i]; k++;} }
  strcpy(init,result);
  delete result;
}
// ----------- Code number --------------------------------------------------
   Integer DatFile::FindDigits(const Integer x)
{  return (x ? 1+DatFile::FindDigits(x/10):0);
}
   Integer DatFile::Step(Integer n, Integer a)
{   Integer b=a;
    for (Integer i=1; i < n; i++) b*=a;
    return b;}

   Integer DatFile::Cryptogr(const Integer x) 
 { Integer n=DatFile::FindDigits(x); 
   return (x + 9*DatFile::Step(n,10));
 }

   Integer DatFile::Encode(const Integer x)
{ if ( x < 20) return 1;
  Integer n=DatFile::FindDigits(x);
  return (x-9*DatFile::Step(n-1,10));
}

// ------------ Count number of character a in std::string buf  -----------------

  Integer DatFile::CountCharInStr(std::string & buf, const Char a)
{  
   Integer count=0;
   for (Integer i=0; i < buf.size(); i++)
   { if (buf.c_str()[i]==a) count++;}
   return count;
}
 
// --------- Count number of words  in std::string. It is needfull for ReadOut -----

  Integer DatFile::CountElemInString(std::string::size_type pos)
{ 
  Char a,c;
  Integer count;
  infile.seekg(pos, std::ios::beg);
  do {infile.get(c);
     } while(c == ' ' &&  c!= '\n');
  if (c=='\n') {count=0;}
  else { count=1; a=c;
       while (c!= '\n') {
           if (a ==' ' && c!=' ') count++;
           a=c;
           infile.get(c); }
      }
   return count;
}
}
