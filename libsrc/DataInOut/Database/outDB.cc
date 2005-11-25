#include "outDB.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"

namespace CoupledField
{

  WriteResultsDatabase::WriteResultsDatabase(const Char * const filename)
                                             
    :WriteResults(filename) {
    ENTER_FCN( "WriteResultsDatabase::WriteResultsDatabase" );
  }

  WriteResultsDatabase::~WriteResultsDatabase()
  {
    ENTER_FCN("WriteResultsDatabase::~WriteResultsDatabase");
    Db_.Close();
  }

  void WriteResultsDatabase::Init (Grid *aptgrid)
  {
    ENTER_FCN("WriteResultsDatabase::Init");
   
    InitHistoryFiles();
    
    std::string hostName, userName, passwd, databaseName;
    Integer port;

    ptGrid_ = aptgrid;

    // get connection parameters
    params->Get("hostName", hostName, "database");
    params->Get("port", port, "database");
    params->Get("databaseName", databaseName, "database");
    params->Get("userName", userName, "database");
    params->Get("password", passwd, "database");

    Db_.Connect (hostName, port, userName,
                 passwd, databaseName);

    WriteBasisData();
  }

  void WriteResultsDatabase::WriteBasisData()
  {
    ENTER_FCN("WriteResultsDatabase::WriteBasisData");
    WriteConfFile();              // save .conf-file in database

    if (!ptGrid_)
      Error("ptGrid_ is not initialized", __FILE__,__LINE__);

    Integer dim         = ptGrid_->GetDim();

    dbLineData d("Calculation");
    d.Set("idx","0");
    d.Set("stamp","NULL",FALSE);
    d.Set("inputparam_idx",InputParamIdx_);
    d.Set("solution_type","1");
    d.Set("analysis_type","1");
    d.Set("model_dimension",dim);
    CalculationIdx_ = Db_.InsertAndGetIndex(d);
  }

  void WriteResultsDatabase::WriteGrid ()
  {
    ENTER_FCN("WriteResultsDatabase::WriteGrid");

    if (!Db_.IsConnected)
      Error(" Database is not initialized",__FILE__,__LINE__);
  
    dbLineData d("Result");
    d.Set("idx","0");
    d.Set("calculation_idx",CalculationIdx_);
    ResultIdx_ = Db_.InsertAndGetIndex(d);
    WriteNodeCoordinates();
    WriteElementNodes();

  }

  long int WriteResultsDatabase::WriteNodeCoordinates()
  {
    ENTER_FCN("WriteResultsDatabase::WriteNodeCoordinates");
    if (!ptGrid_)
      Error("ptGrid_ is not initialized", __FILE__,__LINE__);
 
    Integer dim     = ptGrid_->GetDim();
    Integer maxnumnodes = ptGrid_->GetNumNodes();  

    dbLineData d("Node");
    d.Set("idx","0");
    d.Set("result_idx",ResultIdx_);
    long int idx = Db_.InsertAndGetIndex(d);

    Db_.Lock("Node_coordinates");
    Db_.SetMultipleTuple(maxnumnodes);
    for (Integer i=0; i<maxnumnodes; i++)
      {
        d.Clear();
        d.SetTableName("Node_coordinates");
        d.Set("node_label",(i+1));
        if (dim==2)
          {
            Point<2> Point;
            ptGrid_->GetNodeCoordinate(Point,i);
            d.Set("x_coord","0");
            d.Set("y_coord",Point[0]);
            d.Set("z_coord",Point[1]);
            d.Set("node_idx",idx);
          }
        if (dim==3)
          {
            Point<3> Point;
            ptGrid_->GetNodeCoordinate(Point,i);
            d.Set("x_coord",Point[0]);
            d.Set("y_coord",Point[1]);
            d.Set("z_coord",Point[2]);
            d.Set("node_idx",idx);
          }
        Db_.Insert(d);
      }
    Db_.Unlock();
  }


  long int WriteResultsDatabase::WriteElementNodes()
  {
    ENTER_FCN("WriteResultsDatabase::WriteElementNodes");
    if (!ptGrid_)
      Error("ptGrid_ is not initialized", __FILE__,__LINE__);

    Integer dim = ptGrid_->GetDim();

    StdVector<Integer> connect;
    StdVector<Elem*> elemssd;
    Integer elmsgrp=1;

    StdVector<RegionIdType> subdoms;
    ptGrid_->GetVolRegionIds(subdoms);
    Integer i, j, k;
    k = 0;
    Integer elemlabel, elemtypegeo, elemtypephys, subtype, elemgrpno, nofnodes;
    for (i=0; i<subdoms.GetSize(); i++)
      {
        ptGrid_->GetVolElems(elemssd,subdoms[i]);

        for (j=0; j < elemssd.GetSize(); j++)
          {  
            k++; 
            connect=elemssd[j]->connect;

            elemlabel = elemssd[j]->elemNum;

            if (dim==2)
              {     
                switch(connect.GetSize())
                  {
                  case 3: elemtypegeo = 91; break;
                  case 4: elemtypegeo = 94; break;
                  case 6: elemtypegeo = 92; break;
                  case 8: elemtypegeo = 95; break;
                  default: Error("Please, put element type according to unverg-format for this number of nodes per element", 
                                 __FILE__,__LINE__);
                  }

                elemtypephys = 2;
                subtype      = 2;
                elemgrpno    = elmsgrp;
                nofnodes     = connect.GetSize();
              }

            else
              {
                switch(connect.GetSize())
                  {
                  case 4: elemtypegeo = 111; break;
                  case 6: elemtypegeo = 112; break;
                  case 8: elemtypegeo = 115; break;
                  case 15: elemtypegeo = 113; break;
                  case 20: elemtypegeo = 116; break;
                  default: Error("Please, put element type according to unverg-format for this number of nodes per element", 
                                 __FILE__,__LINE__);
                  }
              
                elemtypephys = 11;
                subtype      = 1;
                elemgrpno    = elmsgrp;
                nofnodes     = connect.GetSize();

              }

            dbLineData d("Element");
            d.Set("idx","0");
            d.Set("elem_label",elemlabel);
            d.Set("elem_type_geo",elemtypegeo);
            d.Set("subtype",subtype);
            d.Set("elem_grp_no",elemgrpno);
            d.Set("result_idx",ResultIdx_); 
            unsigned long int idx = Db_.InsertAndGetIndex(d);

            Db_.Lock("Element_nodes");        
            if (dim == 2 && (connect.GetSize() == 6 || connect.GetSize() == 8))
              {
                //quadratic elements
                Integer offset = Integer(connect.GetSize()/2);
                Db_.SetMultipleTuple(2*offset);
                for (Integer ii=0; ii < offset; ii++)
                  {
                    d.Clear(); d.SetTableName("Element_nodes");
                    d.Set("element_idx",idx);
                    d.Set("node",connect[ii]);
                    d.Set("localidx",ii*2);
                    Db_.Insert(d);
                    d.Clear(); d.SetTableName("Element_nodes");
                    d.Set("element_idx",idx);
                    d.Set("node",connect[ii+offset]);
                    d.Set("localidx",ii*2+1);
                    Db_.Insert(d);
                  }

              }
            else
              {
                Db_.SetMultipleTuple(connect.GetSize());
                for (Integer ii=0; ii < connect.GetSize(); ii++) 
                  { 
                    d.Clear(); d.SetTableName("Element_nodes"); 
                    d.Set("element_idx",idx);
                    d.Set("node",connect[ii]);
                    d.Set("localidx",ii);
                    Db_.Insert(d);
                  }
              }
            Db_.Unlock();

          } // over elements of group
        elmsgrp++;
      } // over groups
  
  
  }


  void WriteResultsDatabase::WriteNodalResult(const std::string & title, 
                                              const Vector<Double> & x, 
                                              const Integer step, 
                                              const Double time, 
                                              const Integer nrNodes,
                                              const Integer nrDofs)
  {
    ENTER_FCN("WriteResultsDatabase::WriteNodalResult");
    if (!ptGrid_)
      Error("ptGrid_ is not initialized", __FILE__,__LINE__);

    Integer valsPerNode = 1;
    if (nrDofs > 1)
      valsPerNode = 3;

    dbLineData d("");
    Integer transIdx=0,normIdx=0,freqIdx=0,Idx=0;
    if (1)                // Transient analysis
      {
        d.SetTableName("Nodal_result_trans");
        d.Set("idx","0"); 
        d.Set("result_type",title); 
        d.Set("time_step",step); 
        d.Set("time_value",time);
        transIdx = Db_.InsertAndGetIndex(d);
      }

    d.Clear(); d.SetTableName("Nodal_result");
    d.Set("idx","0"); d.Set("analysis_type", 4);
    d.Set("nodal_result_trans_idx",transIdx);
    d.Set("nodal_result_norm_idx", normIdx);
    d.Set("nodal_result_freq_idx", freqIdx);
    if (valsPerNode==1)
      d.Set("data_characteristic",1); // data characteristics = scalar
    else if ((valsPerNode>1)&&(valsPerNode<=3))
      d.Set("data_characteristic",2);  // data characteristics = vector with 3 components
    else
      d.Set("data_characteristic",3); // data  characteristics = vector with 6 components
    d.Set("data_type",2);
    d.Set("values_per_node",valsPerNode);
    d.Set("result_idx",ResultIdx_);
    Idx = Db_.InsertAndGetIndex(d);

    Integer i,j,n;
    n = nrNodes;
    Db_.Lock("Nodal_result_value");
    Db_.SetMultipleTuple(n*nrDofs);
    for (i=0; i<n;i++)
      {

        for (j=0; j<nrDofs; j++)
          {
            d.Clear(); d.SetTableName("Nodal_result_value");
            d.Set("nodal_result_idx", Idx);
            d.Set("result", x[i*nrDofs+j]);
            d.Set("node_no", (i+1));
            d.Set("dof", (j+1));
            Db_.Insert(d);
          }
      }
    Db_.Unlock();
  }

  void WriteResultsDatabase::WriteNodeSolutionTransient (const NodeStoreSol<Double>&sol, 
                                                         const Integer step, 
                                                         const Double time)
  {
    ENTER_FCN("WriteResultsDatabase::WriteNodeSolutionTransient");

    /* Integer i,j;
       Integer nrDofs = 1;
       Double help; */

    /* if (NeedHistory_) 
       {
       for (i=0; i< nodeshist_.size(); i++)
       {
       if (sol.GetDof() * sol.GetNumNodes() <= nodeshist_[i])
       Error("Please, check history-nodes.",__FILE__,__LINE__);
       if (lastsavetime[i] != time )
       {
       if (nrDofs > 1)  
       {
       std::vector<Double> solVec;
       solVec.resize(nrDofs);
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
       }
       }
       else*/  
    Vector<Double> globalSolution;
    StdVector<SolutionType> solTypes;
    sol.GetSolutionTypes(solTypes);
    std::string title;
    Integer numNodes =  ptGrid_->GetNumNodes();

    for (Integer iSol=0; iSol<solTypes.GetSize(); iSol++)
      {
        sol.GetGlobalSolVector(solTypes[iSol],globalSolution);
        title = SolutionTypeToString(solTypes[iSol]);
        WriteNodalResult(title, globalSolution, step, time,
                         numNodes,sol.GetDof(solTypes[iSol]));
      }

  }

  /*void WriteResultsDatabase::AddInHistory(const Double time, 
    const Double val, 
    int nodeidx)
    {
    ENTER_FCN("WriteResultsDatabase::AddInHistory");
    if (nodeidx>=nodehistoryindex.size())
    {
    dbLineData d("NodeHistory");
    d.Set("idx","0");
    d.Set("node_no",nodeshist_[nodeidx]);
    d.Set("result_idx",ResultIdx_);
    int idx = Db_.InsertAndGetIndex(d);
    nodehistoryindex.push_back(idx);
    }
    lastsavetime[nodeidx] = time;
    dbLineData v("NodeHistoryValue");
    v.Set("nodehistory_idx",nodehistoryindex[nodeidx]);
    v.Set("dof",0);
    v.Set("time",time);
    v.Set("value",val);
    Db_.Insert(v);
    }*/

  /*void WriteResultsDatabase::AddVecInHistory(const Double time, const std::vector<double> val, const int nodeidx)
    {
    ENTER_FCN("WriteResultsDatabase::AddVecInHistory");
    if (nodeidx>=nodehistoryindex.size())
    {
    dbLineData d("NodeHistory");
    d.Set("idx","0");
    d.Set("node_no",nodeshist_[nodeidx]);
    int idx = Db_.InsertAndGetIndex(d);
    nodehistoryindex.push_back(idx);
    }
    lastsavetime[nodeidx] = time;
    dbLineData v("NodeHistoryValue");
    for (int i=0; i<val.size(); i++)
    {
    v.Clear();
    v.SetTableName("NodeHistoryValue");
    v.Set("nodehistory_idx",nodehistoryindex[nodeidx]);
    v.Set("dof",i);
    v.Set("time",time);
    v.Set("value",val[i]);
    Db_.Insert(v);
    }
    }*/

  void WriteResultsDatabase::WriteElemSolutionTransient (const ElemStoreSol<Double>&sol, 
                                                         const Integer step, 
                                                         const Double time)
  {
    ENTER_FCN("WriteResultsDatabase::WriteElemSolution");
    Vector<Double> globalSolution;
    StdVector<SolutionType> solTypes;
    std::string title;
    Integer numElems =  ptGrid_->GetNumVolElems();  
  
    sol.GetSolutionTypes(solTypes);
    sol.TransformElemSolution(globalSolution,ptGrid_);
    title = SolutionTypeToString(solTypes[0]);
    WriteElementResult(title, globalSolution, step, time,numElems, sol.GetDof());
  }


  void WriteResultsDatabase::WriteConfFile()
  {
    ENTER_FCN("WriteResultsDatabase::WriteConfFile");

    std::string filename = namefile_ + ".xml";

    std::ifstream istr;
    istr.open(filename.c_str(),std::ios_base::in);
    if (!istr)
      {
        Warning ("Cannot open ConfigFile for reading",__FILE__,__LINE__);
        return;
      }  
    std::stringstream configstream;
    char ch;
    while (!istr.eof())
      {
        istr.get(ch);
        configstream << ch;
      }
    istr.close();

    struct stat fileinfo;
    unsigned long int seconds;
    if (stat(filename.c_str(),&fileinfo)==0)
      {
        seconds = (unsigned long int) fileinfo.st_mtime;
      }
    else
      {
        time_t tm;
        time (&tm);
        seconds = (unsigned long int) tm;
      }

    // Check if file is already in database
    std::stringstream wherestr;
    wherestr<<"(filename='"<<filename;
    wherestr<<"') AND (date_modified=FROM_UNIXTIME("<<seconds<<"))";
    Db_.SelectFrom("idx,filename,date_modified,no_references",
                   "InputParam",wherestr.str());
    dbMatrix mat;
    Db_.FetchFields(mat);
    if (mat.getNoOfRow()==1)
      {
        int ref; 
        dbColumn *nOfRef;
        nOfRef = mat["no_references"];
        nOfRef->get(ref,0);
        std::stringstream Set;
        Set<<"no_references="<<(ref+1);
        int idx;
        dbColumn *pIdx;
        pIdx = mat["idx"];
        pIdx->get(idx,0);
        std::stringstream Where;
        Where<<"idx="<<idx;
        Db_.Update("InputParam",Set.str(),Where.str());
        InputParamIdx_ = idx;
        return;
      }
    else if (mat.getNoOfRow()>1)
      Warning("Input file exists more than once in database...");
  

    std::string configstring;
    configstring = Db_.EscapeString(configstream.str());
    dbLineData d("InputParam");
    d.Set("idx","0");
    d.Set("filename",filename);
    std::stringstream moddatestr;
    moddatestr<<"FROM_UNIXTIME("<<seconds<<")";
    d.Set("date_modified", moddatestr.str(),FALSE);
    d.Set("file",configstring);
    InputParamIdx_ = Db_.InsertAndGetIndex(d);
  }

  void WriteResultsDatabase::WriteElementResult(const std::string &title, 
                                                const Vector<Double> & x, 
                                                const Integer step, 
                                                const Double time, 
                                                const Integer nrElems,
                                                const Integer nrDofs)
  {
    ENTER_FCN("WriteResultsDatabase::WriteElementResult");
    if (!ptGrid_)
      Error("ptGrid_ is not initialized", __FILE__,__LINE__);

    Integer valsPerNode = 1;
    if (nrDofs > 1)
      valsPerNode = 3;

    dbLineData d("");
    Integer transProbIdx=0,normProbIdx=0,freqProbIdx=0, transParamIdx=0,normParamIdx=0,freqParamIdx=0, Idx=0;
    if (1)                // Transient analysis
      {
        d.Clear(); 
        d.SetTableName("Element_result_trans");
        d.Set("idx","0"); 
        d.Set("result_type",title); 
        d.Set("time_step", time); 
        d.Set("time_value",time);
        transProbIdx = Db_.InsertAndGetIndex(d);
      }

    d.Clear(); 
    d.SetTableName("Element_result");
    d.Set("idx","0");
    d.Set("analysis_type",4);
    d.Set("element_result_trans_idx", transProbIdx);
    d.Set("element_result_norm_idx", normProbIdx);
    d.Set("element_result_freq_idx", freqProbIdx);

    if (valsPerNode==1)
      d.Set("data_characteristic",1);
    else if ((valsPerNode>1)&&(valsPerNode<=3))
      d.Set("data_characteristic",2);
    else
      d.Set("data_characteristic",3);
    d.Set("data_type",2);
    d.Set("values_per_node", valsPerNode);
    d.Set("result_idx",ResultIdx_);
    Idx = Db_.InsertAndGetIndex(d);

    Integer i,j,n;
    n = nrElems;
    Db_.Lock("Element_result_value");
    Db_.SetMultipleTuple(n*nrDofs);
    for (i=0; i<n;i++)
      {
        d.SetTableName("Element_result_value");
        for (j=0; j<nrDofs; j++)
          {
            d.Clear();
            d.Set("element_result_idx",Idx);
            d.Set("result", x[i*nrDofs+j]);
            d.Set("elem_no", (i+1));
            d.Set("dof", (j+1));
            Db_.Insert(d);
          }
      }
    Db_.Unlock();
  }


  void WriteResultsDatabase::InitHistoryFiles() {

    ENTER_FCN( "WriteResultsDatabase::InitHistoryFiles" );

    (*warning) << "Writing a history is not implemented for "
               << "WriteResultsDatabase";
    Warning( __FILE__, __LINE__ );
  }


  std::string
  WriteResultsDatabase::SolutionTypeToString(const SolutionType type) const {

    ENTER_FCN( "WriteResultsDatabase::SolutionTypeToString" );

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
        return "fluid potential, 2nd deriv.";
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
      case HEAT_TEMPERATURE:
        return "temperature";
        break;
      default:
        (*error) << "Wrong type of solution or 'SolutionType2String' not "
                 << "implemented for this type of solution";
        Error( __FILE__, __LINE__ );
      }
  }

} // End of namespace
