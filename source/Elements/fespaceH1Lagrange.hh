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
#include "Elements/H1Elems.hh"
#include "Elements/H1LagrangeVar.hh"


namespace CoupledField {

class FeSpaceH1Lagrange : public FeSpaceH1 {

  public:

    //! Constructor
    FeSpaceH1Lagrange();

    //! Destructor
    ~FeSpaceH1Lagrange();

    ////! Return equation numbers
    //virtual void GetEqns( StdVector<Integer>& eqns, const EntityIterator ent ); 

    ////! Return equation numbers for a specific dof
    //virtual void GetEqns( StdVector<Integer>& eqns, const EntityIterator ent
    //                      , UInt dof ); 

    //! Map equations i.e. intialize object
    virtual void Finalize();

    //! Reorder the equation Map (just for comptibility)
    virtual void ReorderEqnMap( StdVector<UInt> newOrder );

    //! print the equation map
    virtual void PrintEqnMap();

    virtual void SetMapType(MappingType mapT);

  protected:
      
    //! Map BC Equation NUmbers
    virtual void MapBCs();

    //! This array stores all nodes, including the virtual ones, which are available in the feSpace
    StdVector<UInt> nodes;

    virtual void CreateVirtualNodes();

  private:
};
}
