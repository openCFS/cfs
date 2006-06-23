#ifndef FILE_CFS_RESULTDOF_HH
#define FILE_CFS_RESULTDOF_HH

#include "General/environment.hh"
#include "Domain/entityList.hh"
#include "Utils/StdVector.hh"
#include "Domain/ansatzFct.hh"

namespace CoupledField {

  //! Forward class declaration
  class AnsatzFct;
  class EntityList;

  struct ResultDof {
    
  public:

    //! Typedef for the unknown entities where the result is defined on
    typedef enum{NODE, EDGE, SURFACE, ELEMENT, REGION, 
                 NODELIST, FREE} EntityUnknownType;
    

    //! Friend declaration for operator==
    friend bool operator==(ResultDof& a, ResultDof& b );

    //! Constructor
    ResultDof();

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
    ResultDof& operator=( const ResultDof& data );


    // =======================================================================
    // D A T A    M E M B E R S 
    // =======================================================================

    //! Physical type of result
    SolutionType resultType;

    //! Number of degrees of freedoms
    StdVector<std::string> dofNames;

    //! Type of entity the unkowns are defined on
    EntityUnknownType definedOn;

    //! Type of approximation used for the result
    shared_ptr<AnsatzFct> fctType;

  };

  //! Comparison operator
  bool operator==( const ResultDof& a, const ResultDof&b );

  //! Negative comparison operator
  bool operator!=( const ResultDof& a, const ResultDof&b );

  //! external operator for comparing two resultDofs
  bool operator<( const ResultDof& a, const ResultDof& b );
}

#endif
