#include "database.hh"

namespace CoupledField
{

Database::Database()
{
  ENTER_FCN("Database::Database");
  Conn_ = mysql_init(NULL);
  if (Conn_==NULL)
  {
#ifdef DEBUG
  (*debug) << "Database::Database() mysql_init failed" << std::endl;
#endif
  }
  IsConnected = FALSE;
  Port_ = 0;
  ResultStored_ = FALSE;
  FieldsStored_ = FALSE;
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
  (*debug)<<"DB:"<<Db_<<" on "<<Hostname_<<std::endl;
  if (mysql_real_connect(Conn_, Hostname_.c_str(), User_.c_str(), Passwd_.c_str(), Db_.c_str(), Port_, NULL, 0)==NULL)
  {
    std::string err("Can't establish connection to");
    err += Db_;
    Error(err.c_str(), __FILE__, __LINE__);
    return (-1);
  }
#ifdef DEBUG
  (*debug) << "Connected to database " << Db_ << std::endl;	
  (*debug) << "hostName = " << Hostname_ << std::endl;
  (*debug) << "port = " << Port_ << std::endl;
  (*debug) << "databaseName = " << Db_ << std::endl;
  (*debug) << "userName = " << User_ << std::endl;
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
#ifdef DEBUG2
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
#ifdef DEBUG
//    (*debug)<<"numfields = "<<numfields<<std::endl;
#endif 
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
      case FIELD_TYPE_BLOB:          // use only for text fields, not for real blobs
      case FIELD_TYPE_STRING:
        matrix.appendStringColumn(std::string(fields[i].name));
        break;
      default:
#ifdef DEBUG        
        (*debug)<< "Try to add column of unknown type. Append it as std::string instead"<<std::endl;
#endif
        matrix.appendStringColumn(std::string(fields[i].name));
    }  
  }
  MYSQL_ROW row;			 // line index
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
// int size = d.Size();
 std::stringstream querystr;
 querystr<<"INSERT INTO "<<d.tablename<<" (";
/* int i;
 for (i=0; i<size; i++)
 {
   querystr<<d.GetFieldName(i);
   if (i!=size-1)
     querystr<<",";
 }*/
 querystr<<d.GetFieldNames();
 querystr<<") values (";
/* for (i=0; i<size; i++)
 {
   querystr<<d.GetValue(i);
   if (i!=size-1)
     querystr<<","; 
 }
 querystr<<")";*/
 querystr<<d.GetValues();
 querystr<<")";
 return Query(querystr.str());
}  

std::string Database::EscapeString(std::string source)
{  
 ENTER_FCN("Database::EscapeString");
  unsigned long res;
  unsigned long position;
  char szOut[2*SIZEOFPARAMBLOCK+1];	// Worst case: escape char for every char + terminating zero
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

/*Boolean Database::Update (dbLineData &d, dbLineData &where)
{
  int size = d.Size();
  std::stringstream querystr;
#ifdef DEBUG
  if (d.tablename!=where.tablename)
    (*debug)<<"Database::Update : Table names differ!!"<<std::endl;
#endif
  querystr<<"UPDATE "<<d.tablename<<" SET ";

  char *pStr = */

/*  for (int i=0; i<size; i++)
  {
    querystr<<d.GetFieldName(i)<<"="<<d.GetValue(i);
    if (i!=size-1)
      querystr<<",";
  }
  int wsize = where.Size();
  if (wsize>0)
  {
    querystr<<" WHERE ";
    for (int i=0; i<wsize; i++)
    {
      querystr<<where.GetFieldName(i)<<"="<<where.GetValue(i);
      if (i!=wsize-1)
        querystr<<" AND ";
    }
  }*/
/*  if (Query(querystr.str())==-1)
    return FALSE;*/
/*  return TRUE;
}*/

int Database::InsertAndGetIndex (dbLineData &d) //std::string table, std::vector<std::string> field, std::vector<std::string> value)
{
  ENTER_FCN("Database::InsertAndGetIndex");
  Insert (d);
//   Insert (table, field, value);
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
#ifdef DEBUG
  (*debug)<<"Error store the result"<<std::endl;
#endif
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

Database::~Database()
{
  ENTER_FCN("Database::~Database");
  FreeResult();
}

} // end of namespace

