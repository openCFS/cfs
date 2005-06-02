#ifndef FILE_SUPERBLOCKEQN_2004
#define FILE_SUPERBLOCKEQN_2004

#include "nodeEQN.hh"

namespace CoupledField
{

  //! Class for numbering equations in a superblock way.

  class SuperBlockEQN : public NodeEQN 
  {
  public:
  
    //! Constructor
    SuperBlockEQN(Grid * aptGrid, 
                  StdVector<RegionIdType> & asubdoms, 
                  UInt dofsPerNode);
  
    virtual ~SuperBlockEQN();
  
    //! Calculate the mapping after Dirichlet and
    //! constraint nodes were set
    void CalcMapping();
  
    //! Print the mapping nodes <-> EQNs
    void Print(std::ostream & out) const;
  
    //! 
    UInt GetNumMechEQNs()
    {return numMechEQNs_;}

    //!
    UInt GetNumElecEQNs()
    {return numElecEQNs_;}

    // -----------------------------------------------------------------------
    // Functions for mapping node numbers <-> EQN numbers
    // -----------------------------------------------------------------------
  
    //! Map node number and dof to according equation number
    void Node2EQN(const UInt nodeNr, 
                  const UInt dof,
                  Integer & eqnNr,
                  UInt & eqnDof) const;
  
    //! Map node number to according equation number(s)
    void Node2EQN(const UInt nodeNr, StdVector<Integer> &eqns) const;
  
    //! Map vector of node numbers to according
    //! vector of equiation numbers
    void Node2EQN(const StdVector<UInt> &nodeNr,
                  StdVector<Integer> &eqnNr) const;

    //! Maps the equation numbers according to the reordering
    void ReorderMapping(Integer *order) {
      Error("ReorderMapping not working for SuperBlockEQN");
    }
  

  private:

    //! number of mechanic eqns
    UInt numMechEQNs_;

    //! number of electric eqns
    UInt numElecEQNs_;

    //! Contains the equation numbers
    Matrix<Integer> mechNode2EQN_;
  
    //! Contains the according column numbers
    StdVector<Integer> elecNode2EQN_;

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class SuperBlockEQN
  //! 
  //! \purpose 
  //!  This class is hardcoded for use with piezoPDE,
  //! since in summer there (hopefully) will be a general
  //! mechanism to handel direct-coupled PDEs, so that
  //! this class will superfluous.
  //! 
  //! \collab 
  //! 
  //! \implement 
  //! 
  //! \status 
  //! Deprecated. Only used for the singlepde PiezoPDE, which is also 
  //! deprecated
  //! 
  //! \unused 
  //! 
  //! \improve
  //! 

#endif
  
}  // end of namespace

#endif // FILE_SUPERBLOCKEQN
