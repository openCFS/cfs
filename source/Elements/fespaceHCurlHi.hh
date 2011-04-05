// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
#ifndef FILE_CFS_FESPACE_HCURL_HI_HH
#define FILE_CFS_FESPACE_HCURL_HI_HH

#include "Elements/fespaceH1.hh"
#include <boost/array.hpp>

namespace CoupledField {


//! Finite element space for hierarchical H1 elements
class FeSpaceHCurlHi : public FeSpaceH1 {

public:

  //! Constructor
  FeSpaceHCurlHi(ParamNode* aNode);

  //! Destructor
  ~FeSpaceHCurlHi();
  
  //! Initialize class (read order etc.)
  void Init();
  //! Set Solution strategy and solution step
    void SetStrategy( SolStrategyType strategy, UInt step );

  //! Return pointer to reference element
  virtual BaseFE* GetFe( const EntityIterator ent );

  //! Return pointer to reference element (by element number)
  virtual BaseFE* GetFe( UInt elemNum );

  //! @see FeSpace::GetNumEntityOrder
  UInt GetEntityOrder( UInt elemNum, BaseFE::EntityType type, 
                       UInt entityNum, UInt comp = 1 );

  //! @see FeSpace::GetMaxEntityOrder
  UInt GetMaxEntityOrder( UInt elemNum, BaseFE::EntityType type, 
                          UInt entityNum );

  //! Map equations i.e. initialize object
  virtual void Finalize();

  //! Set type of map
  virtual void SetMapType(MappingType mapT);

protected:

  //! Specialized version of the method of the base class
  //! We number the lowest order dofs consecutively 
  virtual void CreateVirtualNodes();
  
  //! Adjust orders of edges / faces according to min/max rule
  void AdjustEntityOrder();

  //! Virtual node numbers for lowest order edge functions
  boost::array<UInt,2> nedelecNodeRange_;
  
  //! Virtual node numbers for higher order edge functions
  boost::array<UInt,2> edgeHiNodeRange_;
  
  //! Virtual node numbers for higher order face functions
  boost::array<UInt,2> faceNodeRange_;
  
  //! Virtual node numbers for higher order inner functions
  boost::array<UInt,2> innerNodeRange_;
  
  
  //! Map containing the polynomial order for every edge
  std::map<UInt, StdVector<UInt> > edgeOrder_;

private:
};
} // end of namespace
#endif // header guard
