#ifndef FILE_BASEOPERATOR_2003
#define FILE_BASEOPERATOR_2003

#include "General/environment.hh"
#include "Utils/StdVector.hh"
#include "PDE/eqnMap.hh"

namespace CoupledField
{

  class Grid;
  class StdPDE;

  //! Base class for operators working on elements

  //! Class BaseOperator is base class from which different kinds of 
  //! operators are derived

  class BaseOperator
  {
  public:

    //! Constructor
    BaseOperator(Grid * ptGrid, 
                 StdPDE * ptPDE,  
                 shared_ptr<EqnMap> eqnMap,
                 bool isaxi=false,
                 bool coordUpdate_ = false );
  
    //! Destructor
    virtual ~BaseOperator() = 0;

  protected:

    Grid * ptGrid_;     //!< pointer to grid
    StdPDE * ptPDE_;   //!< pointer to PDE
    shared_ptr<EqnMap> eqnMap_;
    bool isaxi_;
    bool coordUpdate_;
  };

} // end of namespace

#endif
