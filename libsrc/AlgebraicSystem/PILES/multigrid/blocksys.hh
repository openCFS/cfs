#ifndef FILE_BLOCKSYSTEM_CLA
#define FILE_BLOCKSYSTEM_CLA

namespace CoupledField
{

//! Block structure for algebraic system
class BlockSystem : public AlgebraicSystem
{
public:
  //! Constructor
  /*! \param anumsys number of systems
      \param anumgraph number of needed graphs
  */
  BlockSystem(Integer anumsys, Integer anumgraph);

  //! Destructor
  virtual ~BlockSystem();

  ///
  virtual void CalcPrecond();

  ///
  virtual void CalcPrecond(Integer newprecond, Integer mat);

  ///
  virtual void Solve();

  ///
  virtual void Solve(Double * f, Double * u);

  //! Set the solver parameter for each system in the bloch structure
  //! to predefined values
  virtual void CreateParameter();
};

}

#endif // FILE_BLOCKSYSTEM_CLA
