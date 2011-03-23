// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
#ifndef FILE_CFS_FESPACE_H1_LAGRANGE_HH
#define FILE_CFS_FESPACE_H1_LAGRANGE_HH

// =====================================================================================
// 
//       Filename:  fespaceH1Lagrange.hh
// 
//    Description:  This implements the finite element space for Lagrange polinomial based
//                  approximations.
//                  This space is capabale to treat the mapping of equation numbers for 
//                  elements of first and second order as well as arbitrary order tesnsor 
//                  product elements.
//                  Virtual Node Array: Simplify the implementation in the referenece 
//                                      by taking element edge and face orientations
//                                      into account
//                  
// 
//        Version:  0.1
//        Created:  01/21/2010 10:31:09 AM
//       Revision:  none
// 
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
// 
// =====================================================================================

#include "Elements/fespaceH1.hh"
#include "Elements/H1ElemsLagExpl.hh"
#include "Elements/H1ElemsLagVar.hh"


namespace CoupledField {

class FeSpaceH1Lagrange : public FeSpaceH1 {

  public:

    //! Constructor
    FeSpaceH1Lagrange();

    //! Destructor
    ~FeSpaceH1Lagrange();

    //! Return pointer to reference element
    virtual BaseFE* GetFe( const EntityIterator ent );
    
    //! Return pointer to reference element (by element number)
    virtual BaseFE* GetFe( UInt elemNum );
    
    ////! Return equation numbers
    //virtual void GetEqns( StdVector<Integer>& eqns, const EntityIterator ent ); 

    ////! Return equation numbers for a specific dof
    //virtual void GetEqns( StdVector<Integer>& eqns, const EntityIterator ent
    //                      , UInt dof );
    
    //! @see FeSpace::GetEntityOrder
    UInt GetEntityOrder( UInt elemNum, BaseFE::EntityType type, 
                           UInt entityNum, UInt comp = 1 );
    
    //! @see FeSpace:: GetMaxEntityOrder
    UInt GetMaxEntityOrder( UInt elemNum, BaseFE::EntityType type, 
                            UInt entityNum );
    
    //! Map equations i.e. initialize object
    virtual void Finalize();

    //! Set type of map
    virtual void SetMapType(MappingType mapT);

  protected:
      
    //! This array stores all nodes, including the virtual ones, which are available in the feSpace
    //StdVector<UInt> nodes;


  private:
};
}
#endif //
