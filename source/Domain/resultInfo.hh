// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_RESULTINFO_HH
#define FILE_CFS_RESULTINFO_HH

#include "Utils/StdVector.hh"
#include "General/environment.hh"

namespace CoupledField {

  //! Forward class declaration
  //class AnsatzFct;
  class EntityList;

  //! This class describes the resultType
  struct ResultInfo {
    
  public:

    //! Typedef for the unknown entities where the result is defined on
    typedef enum{ NODE, EDGE, FACE, ELEMENT, SURF_ELEM, REGION, 
                  SURF_REGION, NODELIST, COIL, FREE } EntityUnknownType;
    
    //! Typedef describing the entryType of the result
    typedef enum { UNKNOWN, SCALAR, VECTOR, TENSOR, STRING } EntryType;
    

    //! Friend declaration for operator==
    friend bool operator==(ResultInfo& a, ResultInfo& b );

    //! Constructor
    ResultInfo();

    // =======================================================================
    // A C C E S S   F U N C T I O N S
    // =======================================================================
    
    //! Returns for a given dof name the related index position (1-based)

    //! This method returns for a given dof name (u,ep,...)
    //! the according index in the dof vector.
    //! \param dof (in) name of a coordinate component
    //! \return index of the coordinate component
    UInt GetDofIndex( const std::string & dof ) const;
    

    //! Returns for a given dof index the according name
    
    //! This method returns for a given dof index (1,2,3)
    //! the according name (x,y,z,rad,...).
    //! \param dof (in) index of the dof component
    //! \return name of the coordinate component
    std::string GetDofName( const UInt dof ) const;
   
    //! Copy operator
    ResultInfo& operator=( const ResultInfo& data );


    // =======================================================================
    // D A T A    M E M B E R S 
    // =======================================================================

    //! Physical type of result
    SolutionType resultType;

    //! General string representation of result
    std::string resultName;

    //! Number of degrees of freedoms
    StdVector<std::string> dofNames;

    //! Unit of the result
    std::string unit;

    //! Format for complex entries
    ComplexFormat complexFormat;
    
    //! Type of entries (scalar, vectorial, tensor)
    EntryType entryType;

    //! Type of entity the unkowns are defined on
    EntityUnknownType definedOn;

    /////! Type of approximation used for the result
    ///shared_ptr<AnsatzFct> fctType;

    /** Gives back a debug summery of the result info */
    std::string ToString() const; 
 
    //! Conversion from EntityUnknownType to string
    static void Enum2String( EntityUnknownType in, 
                             std::string& out );

    //! Conversion from EntityUnknownType to string
    static void String2Enum( const std::string& in,
                             EntityUnknownType& out );
  };

  //! Comparison operator
  bool operator==( const ResultInfo& a, const ResultInfo&b );

  //! Negative comparison operator
  bool operator!=( const ResultInfo& a, const ResultInfo&b );

  //! external operator for comparing two resultInfo s
  bool operator<( const ResultInfo& a, const ResultInfo& b );
}

#endif
