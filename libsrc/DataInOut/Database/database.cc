#include "database.hh"

namespace CoupledField
{

  Database::Database()
  {
    ENTER_FCN("Database::Database");
    Conn_ = mysql_init(NULL);
    if (Conn_==NULL)
      Error("Database::Database() mysql_init failed",__FILE__,__LINE__);
    IsConnected = FALSE;
    Port_ = 0;
    ResultStored_ = FALSE;
    FieldsStored_ = FALSE;
    TupleNo_      = 1;
    CurTuple      = 0;
    PendTupleNo_  = 0;
  }

  int Database::Connect (std::string hostname, 
                         unsigned int port, 
                         std::string user, 
                         std::string passwd, 
                         std::string db)
  {
    ENTER_FCN("Database::Connect(..)");
    Hostname_ = hostname;
    User_     = user;
    Passwd_   = passwd;
    Db_       = db;
    Port_=port;

    return Connect();
  }

  int Database::Connect ()
  {
    ENTER_FCN("Database::Connect()");
    if (mysql_real_connect(Conn_, Hostname_.c_str(), User_.c_str(), Passwd_.c_str(), 
                           Db_.c_str(), Port_, NULL, 0)==NULL)
      {
        std::string err("Can't establish connection to ");
        err += Db_;
        Error(err.c_str(), __FILE__, __LINE__);
        return (-1);
      }
#ifdef DEBUG
    (*debug) << std::endl;
    (*debug) << "Connected to database " << Db_ << std::endl;     
    (*debug) << "hostName     = " << Hostname_ << std::endl;
    (*debug) << "port         = " << Port_ << std::endl;
    (*debug) << "databaseName = " << Db_ << std::endl;
    (*debug) << "userName     = " << User_ << std::endl;
    (*debug) << std::endl;
#endif
    IsConnected = TRUE;
    return (0);
  }

  void Database::Close ()
  {
    ENTER_FCN("Database::Close()");
    if (IsConnected)
      mysql_close(Conn_);
  }

  int Database::Query (std::string statement)
  {
    ENTER_FCN("Database::Query(string)");
    if (IsConnected==FALSE)
      exit (-1);
#ifdef DEBUG
    (*debug) << "Database::Query : " << statement << std::endl;   
#endif
    FreeResult();
    if (mysql_query(Conn_, statement.c_str()))
      {
        Error("Query failed",__FILE__,__LINE__);
        return (-1);
      }
    ResultStored_=FALSE;
    return (0);
  }


  Boolean Database::Lock (std::string tablename)
  {
    std::stringstream qstream;
    qstream<<"LOCK TABLE "<<tablename<<" WRITE";
    Query(qstream.str());
  }

  Boolean Database::Unlock ()
  {
    Query("UNLOCK TABLES");
  }


  // Not for use with tables with binary data or special characters like binary zero
  Boolean Database::GetTable(dbMatrix &matrix, std::string tablename)
  {
    ENTER_FCN("Database::GetTable");
    SelectFrom("*",tablename,"");

    FetchFields(matrix);
    return TRUE;
  }

  unsigned int Database::GetNoOfResultFields()
  {
    ENTER_FCN("Database::GetNoOfResultFields");
    if (!StoreResults())
      return 0;
    return mysql_num_fields(Result_);
  }

  int Database::GetNoOfResultLines()
  {
    ENTER_FCN("Database::GetNoOfResultLines");
    if (!StoreResults())
      return 0;
    return mysql_num_rows(Result_);
  }

  // Not for use with tables with binary data or special characters like binary zero
  Boolean Database::FetchFields(dbMatrix &matrix)
  {
    ENTER_FCN("Database::FetchFields");
    if (!StoreResults())
      return FALSE;
    unsigned int numfields = GetNoOfResultFields();
    MYSQL_FIELD *fields;
    fields = mysql_fetch_fields(Result_);
    for (int i=0; i<numfields; i++)
      { // Iterate through table, field by field
        switch (fields[i].type)
          {
          case FIELD_TYPE_TINY:
          case FIELD_TYPE_SHORT:
          case FIELD_TYPE_LONG:
          case FIELD_TYPE_INT24:
          case FIELD_TYPE_LONGLONG:
            matrix.appendIntColumn(std::string(fields[i].name));  
            break;
          case FIELD_TYPE_FLOAT:
          case FIELD_TYPE_DOUBLE:
            matrix.appendDoubleColumn(std::string(fields[i].name));
            break;
          case FIELD_TYPE_DATETIME:
          case FIELD_TYPE_BLOB:          // use only for text fields, not for real blobs
          case FIELD_TYPE_STRING:
            matrix.appendStringColumn(std::string(fields[i].name));
            break;
          default:
#ifdef DEBUG2        
            (*debug)<< "Try to add column of unknown type. Append it as std::string instead"
                    <<std::endl;
#endif
            matrix.appendStringColumn(std::string(fields[i].name));
          }  
      }
    MYSQL_ROW row;                         // line index
    int rowIdx=0;
    while ((row=mysql_fetch_row(Result_)))  // line by line
      {
        for (int i=0; i<numfields; i++)      // i=fieldindex
          {
            switch (matrix.getFieldType(i))
              {
              case DBFIELDTYPEINT:
                matrix[i]->append(atoi(row[i]));
                break;
              case DBFIELDTYPEDOUBLE:
                matrix[i]->append(atof(row[i]));
                break;
              case DBFIELDTYPESTRING:
                matrix[i]->append(row[i]);
                break;
              }
          }
        rowIdx++;
      }
    //#ifdef DEBUG
    //  matrix.printMatrixDebug();
    //#endif
    return TRUE;
  }


  int Database::Insert (dbLineData &d) 
  {
    ENTER_FCN("Database::Insert");

#ifdef DEBUG2
    (*debug)<<"Table            = "<<d.tablename<<std::endl;
    (*debug)<<"InsertTableName_ = "<<InsertTableName_<<std::endl;
    (*debug)<<"InsertField_     = "<<InsertField_<<std::endl;
    (*debug)<<"InsertValues_    = "<<InsertValues_.str()<<std::endl;
    (*debug)<<"CurTuple         = "<<CurTuple<<std::endl;
    (*debug)<<"TupleNo_         = "<<TupleNo_<<std::endl;
    (*debug)<<"PendTupleNo_     = "<<PendTupleNo_<<std::endl;
    (*debug)<<std::endl;
#endif

    if (CurTuple==0)
      {
        InsertField_ = d.GetFieldNames();
        InsertTableName_ = d.tablename;
      }
    else
      {
        if (InsertField_!=d.GetFieldNames())
          Error("Try to use multiple tupel for INSERT, but uses different InsertField_.",
                __FILE__,__LINE__);
        if (InsertTableName_!=d.tablename)
          Error("Try to use multiple tuple for INSERT, but uses different tables.",
                __FILE__,__LINE__);
      }

    InsertValues_<<"(";
    InsertValues_<<d.GetValues();
    InsertValues_<<")";

    CurTuple++;

    if (CurTuple<TupleNo_)
      {
        InsertValues_<<",";
        return (0);
      }

    std::stringstream querystr;
    querystr<<"INSERT INTO "<<d.tablename<<" (";
    querystr<<InsertField_;
    querystr<<") VALUES ";
    querystr<<InsertValues_.str();
    Integer qresult=Query(querystr.str());
    CurTuple=0;
    InsertField_="";
    ResetMultipleTuple();
    InsertValues_.str("");
    return qresult;
  }  

  std::string Database::EscapeString(std::string source)
  {  
    ENTER_FCN("Database::EscapeString");
    unsigned long res;
    unsigned long position;
    // Worst case: escape char for every char + terminating zero
    char szOut[2*SIZEOFPARAMBLOCK+1];     
    char szIn[SIZEOFPARAMBLOCK];
    std::string to;
    if (source.length()<SIZEOFPARAMBLOCK)
      {
        mysql_real_escape_string(Conn_,szOut,source.c_str(),source.length());
        to.assign(szOut);
      }
    else
      {
        position=0;
        to.erase();
        std::string substring;
        while (position<source.length())
          {
            substring = source.substr(position,SIZEOFPARAMBLOCK);
            strcpy(szIn, substring.c_str());
            mysql_real_escape_string(Conn_,szOut,szIn,strlen(szIn));
            to += szOut;
            position += SIZEOFPARAMBLOCK;
          }
      }
    return to;
  }

  Boolean Database::SelectFrom(std::string fields, 
                               std::string table, 
                               std::string where)
  {
    ENTER_FCN("Database::SelectFrom");
    std::stringstream statement;
    statement<<"SELECT "<<fields<<" FROM "<<table;
    if (where.size()>0)
      statement<<" WHERE "<<where;
    if (Query(statement.str())==-1)
      return FALSE;
    return TRUE;
  }

  Boolean Database::Update (std::string table, 
                            std::string set, 
                            std::string where)
  {
    std::stringstream querystr;
    querystr<<"UPDATE "<<table<<" SET "<<set<<" WHERE "<<where;
    if (Query(querystr.str())==-1)
      return FALSE;
    return TRUE;
  }

  int Database::InsertAndGetIndex (dbLineData &d) 
  {
    ENTER_FCN("Database::InsertAndGetIndex");
    Insert (d);
    return mysql_insert_id(Conn_);
  }

  Boolean Database::StoreResults()
  {
    ENTER_FCN("Database::StoreResults");
    if (ResultStored_)
      {
        return TRUE;
      }
    Result_ = mysql_store_result(Conn_);
    if (Result_==(MYSQL_RES*)NULL)
      {
        Warning ("mysql_store_result-error",__FILE__,__LINE__);
        return FALSE;
      }
    ResultStored_=TRUE;
    return TRUE;
  }

  Boolean Database::StoreFields()
  {
    ENTER_FCN("Database::StoreFields");
    if (FieldsStored_)
      return TRUE;
    if (!StoreResults())
      return FALSE;
    Fields_ = mysql_fetch_fields(Result_);
    FieldsStored_=TRUE;
    return TRUE;
  }

  void Database::FreeResult()
  {
    ENTER_FCN("Database::FreeResult");
    if (ResultStored_)
      mysql_free_result(Result_);
  }

  void Database::SetMultipleTuple(Integer n)
  {
    ENTER_FCN("Database::SetMultipleTuple");
    if (n>MAXTUPLEPERSTATEMENT)
      {
        TupleNo_    = MAXTUPLEPERSTATEMENT;
        PendTupleNo_= n-MAXTUPLEPERSTATEMENT;
      }
    else
      {
        TupleNo_    = n;
        PendTupleNo_= 0;
      }
  }

  void Database::ResetMultipleTuple()
  {
    ENTER_FCN("Database::ResetMultipleTuple");
    if (PendTupleNo_==0)
      {
        TupleNo_ = 1;
      }
    else if (PendTupleNo_<MAXTUPLEPERSTATEMENT)
      {
        TupleNo_     = PendTupleNo_;
        PendTupleNo_ = 0;
      }
    else
      {
        TupleNo_    = MAXTUPLEPERSTATEMENT;
        PendTupleNo_= PendTupleNo_ - MAXTUPLEPERSTATEMENT;
      }
  }

  Database::~Database()
  {
    ENTER_FCN("Database::~Database");
    FreeResult();
  }

} // end of namespace

