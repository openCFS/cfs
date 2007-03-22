// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

// we need some character conversion routines
#include <string>

#include "BaseParamHandler.hh"
#include "PlainXMLParamHandler.hh"
#include "General/environment.hh"
#include "DataInOut/WriteInfo.hh"

namespace CoupledField {


  // Shortcut for standard type
  typedef std::string::size_type sType;


  // **************************************************************************
  //   Constructors and Destructors
  // **************************************************************************

  // ===============
  //   Constructor
  // ===============
  PlainXMLParamHandler::PlainXMLParamHandler(const char *fname ) {

    ENTER_FCN( "PlainParamHandler::PlainParamHandler" );
    std::string errmsg;   
    // check that file is close
    if  (inFile.is_open())
      inFile.close();

    // clear all flags for the file
    inFile.clear();

    // open input-file
    inFile.open(fname);
    
    // send error-msg, if we can't open input-file
    if (!inFile)
      {
        errmsg  = "The file '";
        errmsg += fname;
        errmsg += "' could not be opened!";
        Info->Error(errmsg.c_str(),__FILE__,__LINE__);
      }
    
    // determine the final position in file
    inFile.seekg(0,std::ios::end);
    pos_end_ = inFile.tellg();
  }


  // ==============
  //   Destructor
  // ==============
  PlainXMLParamHandler::~PlainXMLParamHandler()
  {
    ENTER_FCN( "PlainParamHandler::PlainParamHandler" );

    // close the input-file
    inFile.close();
  }


  // **************************************************************************
  //   Public Methods
  // **************************************************************************


  // ================================================================
  //   Return as string the value of a parameter matching a keyword
  // ================================================================
  void PlainXMLParamHandler::Get( const std::string keyword,
                                  std::string &value,
                                  const std::string section,
                                  const std::string subsection )
  {
    ENTER_FCN( "PlainParamHandler::Get (string)" );

    std::cout<<"section: "<<section<<std::endl;
    std::cout<<"subsection: "<<subsection<<std::endl;   

    // find all elements
    StdVector<std::string> list;
    GetList( keyword, list, section, subsection );

    if (list.GetSize() > 1) {
      std::string errmsg = "There is more than 1 match for key in .xml file.";
      errmsg += "Use GetList for this case.";
      Info->Error( errmsg, __FILE__, __LINE__ );
    }

    // copy to the value
    value = list[0];

  }


  // =================================================================
  //   Return as integer the value of a parameter matching a keyword
  // =================================================================
  void PlainXMLParamHandler::Get( const std::string keyword, int &value,
                                  const std::string section,
                                  const std::string subsection )
  {
    ENTER_FCN( "PlainParamHandler::Get (integer)" );

    // find all elements
    StdVector<std::string> list;
    GetList( keyword, list, section, subsection );

    if (list.GetSize() > 1) {
      std::string errmsg = "There is more than 1 match for key in .xml file.";
      errmsg += "Use GetList for this case.";
      Info->Error( errmsg, __FILE__, __LINE__ );
    }

    value = atoi(list[0].c_str());
  }


  // ================================================================
  //   Return as double the value of a parameter matching a keyword
  // ================================================================
  void PlainXMLParamHandler::Get( const std::string keyword, double &value,
                                  const std::string section,
                                  const std::string subsection )
  {
    ENTER_FCN( "PlainParamHandler::Get (double)" );

    // find all elements
    StdVector<std::string> list;
    GetList( keyword, list, section, subsection );

    if (list.GetSize() > 1) {
      std::string errmsg = "There is more than 1 match for key in .xml file.";
      errmsg += "Use GetList for this case.";
      Info->Error( errmsg, __FILE__, __LINE__ );
    }

    value = atof(list[0].c_str());
  }

  void PlainXMLParamHandler::CGet( const std::string key,
                                   std::string &value,
                                   const std::string attribute,
                                   const std::string aValue,
                                   Integer applyToElem,
                                   const std::string section,
                                   const std::string subsection )
  {
    ENTER_FCN( "PlainParamHandler::CGet (std::string)" );
    exit(-1);
  }

  void  PlainXMLParamHandler::CGet( const std::string key,
                                    Double &value,
                                    const std::string attribute,
                                    const std::string aValue,
                                    Integer applyToElem,
                                    const std::string section,
                                    const std::string subsection )
  {
    ENTER_FCN( "PlainParamHandler::CGet (double)" );
    exit(-1);
  }

  void  PlainXMLParamHandler::CGet( const std::string key,
                                    Integer &value,
                                    const std::string attribute,
                                    const std::string aValue,
                                    Integer applyToElem,
                                    const std::string section,
                                    const std::string subsection )
  {
    ENTER_FCN( "PlainParamHandler::CGet (Integer)" );
    exit(-1);
  }
  
  // ====================================================
  //   Return list of strings values matching a keyword
  // ====================================================
  void PlainXMLParamHandler::GetList( const std::string key,
                                      StdVector<std::string> &list,
                                      const std::string section,
                                      const std::string subsection )
  {
    ENTER_FCN( "PlainParamHandler::GetList" );

    if (!list.IsEmpty())
      list.Clear();

    StdVector<sType> s_elems; // start pos for the element
    StdVector<sType> e_elems; // end pos for the element
    StdVector<sType> s_attrs; // start pos for the attribute
    StdVector<sType> e_attrs; // end position for the attribute
    StdVector<std::string> keys;

    // indicator that there are elements
    bool match_elems = false;
    // indicator that there are attributes
    bool match_attrs  = false;

    // ++++++++++++++++++++++
    // 1. Look for the keywords in the elements tag
    // ++++++++++++++++++++++

    // form array with keys
    if ( section != "" ) {
      keys.Push_back( section );
    }

    if ( subsection != "" ) {
      keys.Push_back( subsection );
    }
    keys.Push_back(key);

    // array with start position of the section
    StdVector<sType> s_section;
    // array with end position of the section
    StdVector<sType> e_section;

    // find position of all elements
    FindPosElems(keys,1,s_elems,e_elems,s_section,e_section);

    // check that there are elements
    if (s_elems.GetSize())
      match_elems = true;

    // ++++++++++++++++++++
    // 2. Look for the keywords as attributes of the element
    // ++++++++++++++++++++

    // look for the attributes, if we didn't find the elements
    if (!match_elems)       
      {
        
        // find start and end positions of all attributes
        FindPosAttrs(keys,1,s_elems,e_elems,s_section,e_section);  

        // check that there are attributes
        if (s_elems.GetSize()) {
          match_attrs = true;
        }
        else {
          std::string errmsg = "There are not elements nor attributes '";
          errmsg += key + "' in section '" + section + "' and in subsection '";
          errmsg += subsection;
          //  Info->Error( errmsg, __FILE__,__LINE__ );
          Info->Warning( errmsg, __FILE__,__LINE__ );
        }
      }

    // ++++++++++++++++++++
    // 3. Cut the values from the elements or attrs
    // ++++++++++++++++++++
    
    // null all flags of the inFile
    inFile.clear();
    
    // read values of elements
    if (match_elems)
      ReadValuesElem(list,s_elems,e_elems);
    
    // read values of the attributes
    if (match_attrs)
      ReadAttrsElem(list,s_elems,e_elems);
    
  }


  // ===================================================
  //   Return list of Double values matching a keyword
  // ===================================================
  void PlainXMLParamHandler::GetList( const std::string key,
                                      StdVector<Double> &list,
                                      const std::string section,
                                      const std::string subsection )
  {
    ENTER_FCN( "PlainParamHandler::GetList" );

    if (!list.IsEmpty())
      list.Clear();

    // First assemble matches as list of strings
    StdVector<std::string> matches;

    StdVector<sType> s_elems; // start pos for the element
    StdVector<sType> e_elems; // end pos for the element
    StdVector<sType> s_attrs; // start pos for the attribute
    StdVector<sType> e_attrs; // end position for the attribute
    StdVector<std::string> keys;

    // indicator that there are elements
    bool match_elems = false;
    // indicator that there are attributes
    bool match_attrs  = false;

    // ++++++++++++++++++++++
    // 1. Look for the keywords in the elements tag
    // ++++++++++++++++++++++

    // form array with keys
    if ( section != "" ) {
      keys.Push_back( section );
    }
    if ( subsection != "" ) {
      keys.Push_back( subsection );
    }
    keys.Push_back(key);

    // array with start position of the section
    StdVector<sType> s_section;
    // array with end position of the section
    StdVector<sType> e_section;

    // find position of all elements
    FindPosElems(keys,1,s_elems,e_elems,s_section,e_section);

    // check that there are elements
    if (s_elems.GetSize())
      match_elems = true;

    // ++++++++++++++++++++
    // 2. Look for the keywords as attributes of the element
    // ++++++++++++++++++++

    // look for the attributes, if we didn't find the elements
    if (!match_elems)       
      {
        
        // find start and end positions of all attributes
        FindPosAttrs(keys,1,s_elems,e_elems,s_section,e_section);  

        // check that there are attributes
        if (s_elems.GetSize()) {
          match_attrs = true;
        }
        else {
          std::string errmsg = "There are not elements or attributes '";
          errmsg += key + "' in section '" + section + "' and in subsection '";
          errmsg += subsection;
          Info->Warning( errmsg, __FILE__,__LINE__ );
        }
      }

    // ++++++++++++++++++++
    // 3. Cut the values from the elements or attrs
    // ++++++++++++++++++++
    
    // null all flags of the inFile
    inFile.clear();
    
    // read values of elements
    if (match_elems)
      ReadValuesElem(matches,s_elems,e_elems);
    
    // read values of the attributes
    if (match_attrs)
      ReadAttrsElem(matches,s_elems,e_elems);

    // convert matches to Double
    for ( UInt k = 0; k < matches.GetSize(); k++ ) {
      list.Push_back( (Double)atof( matches[k].c_str() ) );
    }
  }

  void  PlainXMLParamHandler::CGetList( const std::string key,
                                        StdVector<std::string> &list,
                                        const std::string attribute,
                                        const std::string value,
                                        Integer applyToElem,
                                        const std::string section,
                                        const std::string subsection )
  {
    
    ENTER_FCN( "PlainParamHandler::CGetList" );
    exit(-1);
    
  }

  void  PlainXMLParamHandler::GetCoilList( StdVector<std::string> &list,
                                           const std::string pde,
                                           const std::string sequenceTag
                                           )
  {
    
    ENTER_FCN( "PlainParamHandler::GetCoilList" );
    exit(-1);
    
  }

  void  PlainXMLParamHandler::GetCoilType( std::string &coilType,
                                           const std::string coilName,
                                           const std::string pde )
  {
    
    ENTER_FCN( "PlainParamHandler::GetCoilType" );
    exit(-1);
    
  }
  
  // =====================================
  //   Return a list of the defined PDEs
  // =====================================
  void PlainXMLParamHandler::GetPDEList( StdVector<std::string> &list )
  {

    ENTER_FCN( "PlainParamHandler::GetPDEList" );

    // Check if vector is empty and erase its entries,
    // if it is not empty
    if ( list.IsEmpty() != true )
      list.Clear();
    
    sType pos0, pos1, pos;
     
    // Find initial position of PDE section 
    pos0 = getposElem("pdeList");
    
    // Find end position of PDE section
    pos1 = getposElem("/pdeList");

    // Get the names of the PDEs and store them in vector
    inFile.seekg(pos0,std::ios::beg);

    std::string pdename;
    std::string buf, buf1;
    pos = 0;

    // read lines in the input-file upto final tag of element PDE_list
    std::getline(inFile,buf);
    do {

      // check, whether buf contains tag
      if (is_tag(buf)) {
        
        peel(buf,pdename);
        
        pos = inFile.tellg();
        buf1 = '/'+ pdename;
        
        pos = getposElem(buf1,pos);       
        inFile.seekg(pos,std::ios::beg);
        inFile.ignore(100,'\n');
        
        list.Push_back(pdename);
      }

      // read next line
      std::getline(inFile,buf);
      pos = inFile.tellg();
      
    } 
    while (pos<pos1);
    
    
    if (list.GetSize()==0)
      Info->Error( "We didn't find any PDE in XML parameter file",
                   __FILE__, __LINE__ );

  }


  // =========================================
  //   Return a list of iterative coupled PDEs
  // =========================================
  void PlainXMLParamHandler::GetIterCoupledPDEList( StdVector<std::string> &list,
                                                    const std::string sequenceTag)
                                        
  {
    ENTER_FCN( "PlainXMLParamHandler::GetIterCoupledPDEList" );
    Error (" Implement me!", __FILE__, __LINE__);
  }
  
  void PlainXMLParamHandler::GetDirectCouplingList( StdVector<std::string> &list,
                                                    const std::string sequenceTag)
  {
    ENTER_FCN( "PlainXMLParamHandler::GetDirectCouplingList" );
    Error (" Implement me!", __FILE__, __LINE__);
  }
  

  // ========================================
  //   Query on/off status of a flag/switch
  // ========================================
  bool PlainXMLParamHandler::IsSet( const std::string key,
                                       const std::string section,
                                       const std::string subsection )
  {
    // assume that flag/switch can be only element

    ENTER_FCN( "PlainXMLParamHandler::IsSet" );
    
    // option 'on' for the flag/switch
    std::string yes="yes";

    return HasValue(key,yes,section,subsection);

  }


  // =================================================
  //   Query whether a parameter has a certain value
  // =================================================
  bool PlainXMLParamHandler::HasValue( const std::string key,
                                          const std::string value,
                                          const std::string section,
                                          const std::string subsection )
  {
    ENTER_FCN( "PlainParamHandler::HasList" );

    StdVector<std::string> list;

    GetList(key,list,section,subsection);

    if (list.GetSize()>1)
      Info->Error("There is more than one element with key "+key+
                  " in subsection " +subsection+ " in section "+section,
                  __FILE__,__LINE__);

    if (list.GetSize() == 0)
      return false;
   

    if (list[0] == value)
      return true;
    else
      return false;
  }


  // ===================================================================
  //   To find an initial position of an element with the name keyword  
  // ===================================================================
  sType PlainXMLParamHandler::getposElem( const std::string keyword,
                                          const sType startpos )
  {

    // have begun a search from the position startpos
    inFile.clear();
    inFile.seekg(startpos,std::ios::beg);


    // read a line from the file and check, whether it contains a required word
    sType help = 0;
    sType pos = std::string::npos;
    std::string buf;

    while ( pos == std::string::npos && !inFile.eof() )
      {
        help=inFile.tellg();
        std::getline(inFile, buf, '\n');

        pos = buf.find(keyword);        
      }
   
    // return std::string::npos, if have found nothing
    if (pos>=pos_end_) 
      return std::string::npos;
   
    // return position of the next line after the line with keyword
    UInt size=buf.size();
    return help+size;
   
  }


  // =================================================
  //   To cut PDE's name from the string
  // =================================================
  void PlainXMLParamHandler::peel( const std::string unpeeled,
                                   std::string & peeled)
  {
    ENTER_FCN( "PlainXMLParamHandler::peel" );

    peeled.clear();

    UInt i;
    for (i=0; i<unpeeled.size(); i++) {
      if ((unpeeled.c_str())[i] == '<')
        break;
    }
      
    i++;
    Char ctmp = unpeeled[i];
  
    while ( ctmp != '>') {
      peeled.push_back(ctmp);
      i++;
      ctmp = unpeeled[i];
    }

  }


  // =================================================
  //   To find final and initial positions of 'sought for' elements 
  // =================================================
  void PlainXMLParamHandler::FindPosElems(StdVector<std::string> keys,
                                          int level,
                                          StdVector<sType> & s_elems,
                                          StdVector<sType> & e_elems,
                                          StdVector<sType> s_section,
                                          StdVector<sType> e_section)
  {
    ENTER_FCN( "PlainXMLParamHandler::FindPosElems" );
       
    if (level == 0 || level>keys.GetSize()) 
      Info->Error("The level is strange, check data",__FILE__,__LINE__);
    
    // get start and end position of the elements
    getElems(keys[level-1], s_elems, e_elems, s_section, e_section);
    
    // recursion
    if (level != keys.GetSize())
      {

        s_section = s_elems;
        e_section = e_elems;

        s_elems.Clear();
        e_elems.Clear();

        FindPosElems(keys,level+1,s_elems,e_elems,s_section,e_section);
      }
    
  }


  // =================================================
  //   To find final and initial positions of 'sought for' elements 
  // =================================================
  void PlainXMLParamHandler::FindPosAttrs(StdVector<std::string> keys,
                                          int level,
                                          StdVector<sType> & s_elems,
                                          StdVector<sType> & e_elems,
                                          StdVector<sType> s_section,
                                          StdVector<sType> e_section)
  {
    ENTER_FCN( "PlainXMLParamHandler::FindPosAttrs" );

    //     if (level == 0 || level>keys.GetSize()) 
    //       Info->Error("The level is strange, check data",__FILE__,__LINE__);

    // get start and end position of the elements
    getElems(keys[level-1], s_elems, e_elems, s_section, e_section);
  
    // do recursion for the section,subsection,etc.
    //orignal if ( (level+1) != keys.GetSize())
    if ( (level+1) < keys.GetSize())
      {
        s_section = s_elems;
        e_section = e_elems;

        s_elems.Clear();
        e_elems.Clear();

        FindPosAttrs(keys,level+1,s_elems,e_elems, s_section, e_section);
      }
    else // last keyword is attribute, so we use fnc getAttr
      if ( (level+1) == keys.GetSize() || keys.GetSize()==1) {

        s_section = s_elems;
        e_section = s_elems;

        s_elems.Clear();
        e_elems.Clear();
        if (keys.GetSize()==1)
          getAttr(keys[0],"",s_elems,e_elems, s_section,
                  e_section);
        else
          getAttr(keys[level],keys[level-1],s_elems,e_elems, s_section,
                  e_section);
      }
        
  }

  
  // =================================================
  //   To find final and initial positions of 'sought for' elements 
  // =================================================
  void PlainXMLParamHandler::getElems(const std::string key,
                                      StdVector<sType> & s_elems,
                                      StdVector<sType> & e_elems,
                                      StdVector<sType> s_section,
                                      StdVector<sType> e_section,
                                      sType start,
                                      const int type)
  {
    ENTER_FCN( "PlainXMLParamHandler::getElems" );
    
    bool testt=false;

    // have begun a search from the position startpos
    inFile.clear();
    inFile.seekg(0,std::ios::beg);

    UInt i;
    // read a line from the file and check, whether it contains the key
    sType       pos=std::string::npos,
      posE=std::string::npos;
    
    pos = findPos('<'+key, start);

    // indicator whether our found key is enclosed in section, subsection, etc.
    bool put=true;

    // we've found the key
    if (pos < pos_end_)
      {

        // check that the key is enclosed
        for (i=0; i<s_section.GetSize(); i++)
          if ((pos < s_section[i]) || (pos > e_section[i]))
            put = false;
        
        // find the end position
        posE = findPos('/'+key,pos);
        
        if (posE == std::string::npos)
          {
            pos = posE;
            put = false;
            testt=true;
          }

        if (put)
          s_elems.Push_back(pos);

        //      Info->Error("We can't find the end-tag for key '" + key + "'\n",
        //                  __FILE__,__LINE__);
        
        if (put)
          e_elems.Push_back(posE);

        if (!testt)
          getElems(key,s_elems,e_elems,s_section,e_section,posE+key.size());
      }

  }


  // ==================================================
  //   To find final and initial positions of 'sought
  //   for' element of the second type 
  // ==================================================
  void PlainXMLParamHandler::getAttr(const std::string attr,
                                     const std::string element,
                                     StdVector<sType> & s_elems,
                                     StdVector<sType> & e_elems,
                                     StdVector<sType> s_section,
                                     StdVector<sType> e_section,
                                     sType start)
  {
    ENTER_FCN( "PlainXMLParamHandler::getAttr" );
    
    // have begun a search from the position startpos
    inFile.clear();
    inFile.seekg(0,std::ios::beg);
    
    // read a line from the file and check, whether it contains the key
    sType       pos=std::string::npos,
      posEEl=std::string::npos,
      posE = std::string::npos;
    
    std::string endAttr, startAttr;
    startAttr.push_back('=');
    startAttr.push_back('"');
    endAttr.push_back('"');
    
    // find the element
    pos = findPos('<'+element, start);

    bool next=true;
    if (pos < pos_end_) {
      
      // check that element is enclosed in sections
      int i;
      for (i=0; i<s_section.GetSize(); i++)
        if ((pos < s_section[i]) || (pos > e_section[i]))
          next = false;
      
      // find end position of the element
      posEEl = findEndPosElementType2(pos+element.size()+1);

      if (posEEl == std::string::npos)
        Info->Error("We can't find the end of the element",__FILE__,__LINE__);

      // find position of the attribute
      pos = findPos(attr, pos);

      // check that the attribute is enclosed in the element
      if (posEEl > pos)
        next = true;
      else
        next = false;
      
      if (next) {
        // find the end of the attribute
        posE = findPos(startAttr,pos);
        posE = findPos(endAttr,posE);
          
        // put the found attribute
        s_elems.Push_back(pos);
        e_elems.Push_back(posE);
      }
      
      // find the next attribute
      getAttr(attr,element,s_elems,e_elems,s_section,e_section,posEEl);
      
    }

  }


  // =================================================
  //   To find position of the keyword
  // =================================================
  sType PlainXMLParamHandler::findPos(const std::string key,
                                      sType start)
  {
    ENTER_FCN( "PlainXMLParamHandler::findPos" );
    
    // have begun a search from the position startpos
    inFile.clear();
    inFile.seekg(start,std::ios::beg);
    
    sType       help,  pos=std::string::npos;
    std::string                  buf;

    // read a line from the file and check, whether it contains the key
    while ( pos == std::string::npos && !inFile.eof() ) {
      
      help=inFile.tellg(); 
      std::getline(inFile, buf, '\n');

      pos = buf.find(key);      
    }

    if (pos == std::string::npos)
      return pos;
    else
      return help+pos;

  }


  // =================================================
  //   To read values of the element
  // =================================================
  void PlainXMLParamHandler::ReadValuesElem(StdVector<std::string> & list,
                                            StdVector<sType> s_elems,
                                            StdVector<sType> e_elems)
  {
    ENTER_FCN( "PlainXMLParamHandler::ReadValuesElem" );

    char tchar;

    // resize list for the resultes
    list.Resize(s_elems.GetSize());
    list.Init();

    // loop over list with given initial and final positions of the elements
    UInt i;
    for (i=0; i<s_elems.GetSize(); i++) {
      
      inFile.seekg(s_elems[i],std::ios::beg);

      // read out the value under letters
      do {
        inFile >> tchar;
        
        if (tchar == '>')
          while (tchar != '<') {
            inFile >> tchar;
            if (tchar != ' ' && tchar != '<')
              list[i].push_back(tchar);     
          } ;
      } while (inFile.tellg() < e_elems[i]);
      

    }
 
  }


  // =================================================
  //   To read attributes of the element
  // =================================================
  void PlainXMLParamHandler::ReadAttrsElem(StdVector<std::string> & list,
                                           StdVector<sType> s_elems,
                                           StdVector<sType> e_elems)
  {
    ENTER_FCN( "PlainXMLParamHandler::ReadAttrsElem" );

    char tchar;
    
    list.Resize(s_elems.GetSize());
    list.Init();
    
    UInt i;
    for (i=0; i<s_elems.GetSize(); i++) {
      
      inFile.seekg(e_elems[i],std::ios::beg);
      
      do {
        inFile >> tchar;
        
        if (tchar == '"') {
          inFile >> tchar;
          while (tchar != '"') {
            if (tchar != ' ')
              list[i].push_back(tchar);     
            inFile >> tchar;
          };
        }
      } while (inFile.tellg() < e_elems[i]);

    }

  }

  
  // =========================================
  //   Read the last position of the element
  // =========================================
  sType PlainXMLParamHandler::findEndPosElementType2(sType start)
  {
    ENTER_FCN( "PlainXMLParamHandler::findEndPosElementType2" );

    // to count subelements in the element
    int counter = 0;
    int counter_end = 0;

    inFile.seekg(start,std::ios::beg);

    char tmp;
    do {
      inFile >> tmp;
      if (tmp == '<')
        counter++;
      if (tmp == '/')
        if (tmp == '>')
          counter_end++;

      if (counter == counter_end)
        return inFile.tellg();
    } while (inFile.tellg() < pos_end_);
    
    // we are in trouble
    Info->Error( "Not as many '>' as '<'!", __FILE__, __LINE__ );

    // get rid of compiler warning
    return 0;
  }


  // =================================================
  //   To test, whether this line contains tag
  // =================================================
  bool PlainXMLParamHandler::is_tag( const std::string buf)
  {

    ENTER_FCN( "PlainXMLParamHandler::is_tag" );

    sType pos, pos1;
    
    pos = buf.find('<');
    pos1 = buf.find('>');

    if (pos != std::string::npos && pos1 != std::string::npos)
      return true;
    else
      return false;

  }





  void PlainXMLParamHandler::GetList (const StdVector< std::string > &keyVec,
                                      const StdVector< std::string > &attrVec,
                                      const StdVector< std::string > &valVec, StdVector< std::string > &list)
  {

    ENTER_FCN( "PlainXMLParamHandler::GetList(StdVector<std::string>,StdVector<std::string>,StdVector<std::string>,StdVector<std::string>)" );

    std::string keyword;
    //   StdVector<Double> list;
    std::string section;
    std::string subsection;
      

    if (keyVec.GetSize()==2)
      {
        section = keyVec[0];
        keyword = keyVec[1];
      }
    if (keyVec.GetSize()==3)
      {
        section = keyVec[0];
        subsection = keyVec[1];
        keyword = keyVec[2];
      }
    if (keyVec.GetSize()>=4)
      {
        if(valVec.GetSize()==3 && valVec[0]==""&& valVec[1]=="")
          {
            section = keyVec[2];
            subsection = keyVec[3];
            keyword = valVec[2];
          }
        else
          {
            section = keyVec[1];
            subsection = keyVec[2];
            keyword = keyVec[3];            
          }
          
      }
  
    GetList( keyword, list, section, subsection );


  }
    

  void    PlainXMLParamHandler::GetList (const StdVector< std::string > &keyVec,
                                         const StdVector< std::string > &attrVec,
                                         const StdVector< std::string > &valVec, StdVector< Double > &list)
  {

    ENTER_FCN( "PlainXMLParamHandler::GetList(StdVector<std::string>,StdVector<std::string>,StdVector<std::string>,StdVector<Double>)" );

    std::string keyword;
    //     StdVector<Double> list;
    std::string section;
    std::string subsection;
      
    section = keyVec[0];
    if (keyVec.GetSize()==2)
      {
        keyword = keyVec[1];
      }
    if (keyVec.GetSize()==3)
      {
        subsection = keyVec[1];
        keyword = keyVec[2];
      }
    if (keyVec.GetSize()>=4)
      {
        section = keyVec[1];
        subsection = keyVec[2];
        keyword = keyVec[3];
      }

    GetList( keyword, list, section, subsection );


  }

  void PlainXMLParamHandler::GetList (const StdVector< std::string > &keyVec,
                                      const StdVector< std::string > &attrVec,
                                      const StdVector< std::string > &valVec, StdVector< Integer > &list)
  {
    Error("This GetList is not yet implemented for PlainXMLParamHandler",__FILE__,__LINE__);
  }

  void    PlainXMLParamHandler::Get (const StdVector< std::string > &keyVec,
                                     const StdVector< std::string > &attrVec,
                                     const StdVector< std::string > &valVec, std::string &value)
  {
    Error("This Get is not yet implemented for PlainXMLParamHandler",__FILE__,__LINE__);
  }

  void    PlainXMLParamHandler::Get (const StdVector< std::string > &keyVec,
                                     const StdVector< std::string > &attrVec,
                                     const StdVector< std::string > &valVec, Double &value)
  {

    ENTER_FCN( "PlainXMLParamHandler::Get(StdVector<std::string>,StdVector<std::string>,StdVector<std::string>,Double)" );

    std::string keyword;
    StdVector<Double> list;
    std::string section;
    std::string subsection;
      
    section = keyVec[0];
    if (keyVec.GetSize()==2)
      {
        keyword = keyVec[1];
      }
    if (keyVec.GetSize()==3)
      {
        subsection = keyVec[1];
        keyword = keyVec[2];
      }
  
    GetList( keyword, list, section, subsection );

    if (list.GetSize() > 1) {
      std::string errmsg = "There is more than 1 match for key in .xml file.";
      errmsg += "Use GetList for this case.";
      Info->Error( errmsg, __FILE__, __LINE__ );
    }

    // copy to the value
    value = list[0];      
     
  }


  void    PlainXMLParamHandler::Get (const StdVector< std::string > &keyVec,
                                     const StdVector< std::string > &attrVec, 
                                     const StdVector< std::string > &valVec, Integer &value)
  {


    ENTER_FCN( "PlainXMLParamHandler::Get(StdVector<std::string>,StdVector<std::string>,StdVector<std::string>,Integer)" );

    std::string keyword;
    StdVector<Double> list;
    std::string section;
    std::string subsection;
      
    section = keyVec[0];
    if (keyVec.GetSize()==2)
      {
        keyword = keyVec[1];
      }
    if (keyVec.GetSize()==3)
      {
        subsection = keyVec[1];
        keyword = keyVec[2];
      }
  
    GetList( keyword, list, section, subsection );

    if (list.GetSize() > 1) {
      std::string errmsg = "There is more than 1 match for key in .xml file.";
      errmsg += "Use GetList for this case.";
      Info->Error( errmsg, __FILE__, __LINE__ );
    }

    // copy to the value
    value = (int)list[0];         
    
  }

   
  void PlainXMLParamHandler::
  GetDim1xDim2Tensor( const StdVector<std::string> &keyVec,
                      const StdVector<std::string> &attrVec,
                      const StdVector<std::string> &valVec,
                      const unsigned int &dim1,
                      const unsigned int &dim2,
                      Matrix<Double> &matr ) {
    Error( "GetDim1xDim2Tensor not implemented for PlainParamHandler",
       __FILE__, __LINE__ );
  }
    
  
  void PlainXMLParamHandler::
  GetDim1xDim2Tensor( const StdVector<std::string> &keyVec,
                             const StdVector<std::string> &attrVec,
                             const StdVector<std::string> &valVec,
                             const unsigned int &dim1,
                             const unsigned int &dim2,
                      Matrix<Integer> &matr ) {
    Error( "GetDim1xDim2Tensor not implemented for PlainParamHandler",
       __FILE__, __LINE__ );
  }

   
  void PlainXMLParamHandler::
  GetDim1xDim2Tensor( const StdVector<std::string> &keyVec,
                      const StdVector<std::string> &attrVec,
                      const StdVector<std::string> &valVec,
                      const unsigned int &dim1,
                      const unsigned int &dim2,
                      Matrix<UInt> &matr ) {
    Error( "GetDim1xDim2Tensor not implemented for PlainParamHandler",
       __FILE__, __LINE__ );
  }

  void PlainXMLParamHandler::
  GetVector( const StdVector<std::string> &keyVec,
             const StdVector<std::string> &attrVec,
             const StdVector<std::string> &valVec,
             StdVector<std::string> &strVec ) {
    Error( " GetVector not implemented for PlainParamHandler",
       __FILE__, __LINE__ );
  }
             

  void PlainXMLParamHandler::GetList( const std::string key, StdVector<Integer> &list,
                                      const std::string section ,
                                      const std::string subsection )   
  {
    Error("This GetList is not yet implemented for PlainXMLParamHandler",__FILE__,__LINE__);
  }

  void PlainXMLParamHandler::Get( const StdVector<std::string> &keyVec,
                                  std::string &value ) 
  {

    ENTER_FCN( "PlainXMLParamHandler::Get(StdVector<std::string>,std::string)" );

    std::string keyword;
    StdVector<std::string> list;
    std::string section;
    std::string subsection;
      

    if (keyVec.GetSize()==2)
      {
        section = keyVec[0];
        keyword = keyVec[1];
      }
    if (keyVec.GetSize()==3)
      {
        section = keyVec[0];
        subsection = keyVec[1];
        keyword = keyVec[2];
      }
    if (keyVec.GetSize()>=4)
      {
        section = keyVec[1];
        subsection = keyVec[2];
        keyword = keyVec[3];
      }
  
    GetList( keyword, list, section, subsection );
    if (list.GetSize() >= 1)
      value=list[0];
    else
      {
        std::string errmsg = "There are not elements nor attributes '";
        errmsg += keyword + "' in section '" + section + "' and in subsection '";
        errmsg += subsection;
        Info->Error( errmsg, __FILE__,__LINE__ );
      }
    


  }


  void PlainXMLParamHandler::Get( const StdVector<std::string> &keyVec,
                                  Double &value )
  {
    Error("This Get is not yet implemented for PlainXMLParamHandler",__FILE__,__LINE__);
  }
  void PlainXMLParamHandler::Get( const StdVector<std::string> &keyVec,
                                  Integer &value )
  {
    Error("This Get is not yet implemented for PlainXMLParamHandler",__FILE__,__LINE__);
  }

}
