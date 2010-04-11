// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
#ifndef FILE_CFS_FESPACE_H1_HI_HH
#define FILE_CFS_FESPACE_H1_HI_HH

#include "Elements/fespaceH1.hh"
#include "Elements/H1ElemsHi.hh"

namespace CoupledField {


//! Finite element space for hierarchical H1 elements
class FeSpaceH1Hi : public FeSpaceH1 {

  public:

    //! Constructor
    FeSpaceH1Hi();

    //! Destructor
    ~FeSpaceH1Hi();

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
      
    //! Adjust orders of edges / faces according to min/max rule
    void AdjustEntityOrder();
    
    //! Map containing the polynomial order for every edge
    std::map<UInt, StdVector<UInt> > edgeOrder_;
    
    // For the polynomial degrees of faces / interior we have to think
    //! Map containing the polynomial order for every face
//    std::map<UInt, UInt> faceOrder_;
    
    //! Map containing the polynomial order for every element interior
//    std::map<UInt, UInt> faceOrder_;

  private:
};
} // end of namespace
#endif // header guard