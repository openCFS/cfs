#ifndef FILE_CFS_RESULTINFO_HH
#define FILE_CFS_RESULTINFO_HH

#include <string>

#include "Utils/StdVector.hh"
#include "General/Environment.hh"
#include "FeBasis/FeFunctions.hh"

namespace CoupledField {

  //! Forward class declaration
  //class AnsatzFct;
  class EntityList;
  class BaseFeFunction;

  //! This class describes the resultType
  struct ResultInfo {
    
  public:

    //! Typedef for the unknown entities where the result is defined on
    typedef enum{ NODE, EDGE, FACE, ELEMENT, SURF_ELEM, REGION, REGION_AVERAGE,
                  SURF_REGION, NODELIST, COIL, FREE } EntityUnknownType;
    static Enum<EntityUnknownType> EntityUnknownTypeEnum_;
    
    //! Typedef describing the entryType of the result
    typedef enum { UNKNOWN = 0, SCALAR = 1, VECTOR = 2, TENSOR = 3, STRING = 4 } EntryType;
    static Enum<EntryType> EntryTypeEnum_;

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

    //! Returns the size of the DoF vector
    UInt GetDimDof() {
      return dofNames.GetSize();
    }
   
    //! Copy operator
    ResultInfo& operator=( const ResultInfo& data );


    /** Set the FeFunction of the current result **/
    void SetFeFunction(shared_ptr<BaseFeFunction > feFunct){
      this->feFct_ = feFunct;
    }

    /** Get the FeFunction of the current result **/
    weak_ptr<BaseFeFunction>  GetFeFunction(){
      return feFct_;
    }

    /** Maps a SolutionType to its canonical EntityUnknownType (definedOn).
     * This function returns the standard entity type where a given solution is defined.
     * For unmapped types, returns ResultInfo::FREE as a fallback.
     * @param solType the solution type
     * @return the EntityUnknownType where this solution is typically defined
     */
    ResultInfo::EntityUnknownType MapSolTypeToDefinedOn(SolutionType solType);

    // =======================================================================
    // D A T A    M E M B E R S 
    // =======================================================================

    //! Physical type of result
    SolutionType resultType;

    //! General string representation of result
    std::string resultName;

    /** Helper to set dof names for vector. Sets to x, y (,z) or r, z */
    void SetVectorDOFs(UInt dim, bool is_axi);

    /** If an extended 2D (2.5D) formulation is used, vectors must have 3 components */
    void SetVectorDOFs(UInt dim, bool is_axi, bool is2p5);

    //! Number of degrees of freedoms
    StdVector<std::string> dofNames;

    //! Unit of the result
    std::string unit;

    //! Format for complex entries
    ComplexFormat complexFormat;
    
    //! Type of entries (scalar, vectorial, tensor)
    EntryType entryType;

    //! Type of entity the unknowns are defined on
    EntityUnknownType definedOn;

    //! Is this result provided from the optimization part
    bool fromOptimization;
    
    //! this is true if the result is not depending on time/frequency
    bool isStatic;

    /** Gives back a debug summary of the result info */
    std::string ToString() const; 
 
    //! Conversion from EntityUnknownType to string
    static void Enum2String( EntityUnknownType in, 
                             std::string& out );

    //! Conversion from EntityUnknownType to string
    static void String2Enum( const std::string& in,
                             EntityUnknownType& out );

  protected:

    weak_ptr<BaseFeFunction> feFct_;

  };

  //! Comparison operator
  bool operator==( const ResultInfo& a, const ResultInfo&b );

  //! Negative comparison operator
  bool operator!=( const ResultInfo& a, const ResultInfo&b );

  //! external operator for comparing two resultInfo s
  bool operator<( const ResultInfo& a, const ResultInfo& b );
}

#endif
