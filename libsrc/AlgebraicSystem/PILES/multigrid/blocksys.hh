#ifndef FILE_BLOCKSYSTEM_PILES
#define FILE_BLOCKSYSTEM_PILES

namespace CoupledField
{

class BlockSystem : public AlgebraicSystem
{
public:
  ///
  BlockSystem(Integer anumsys, Integer anumgraph);

  ///
  virtual ~BlockSystem();

  ///
  virtual void CalcPrecond();

  ///
  virtual void CalcPrecond(Integer newprecond, Integer mat);

  ///
  virtual void Solve();

  ///
  virtual void Solve(Double * f, Double * u);

  ///
  virtual void CreateParameter();
};

}

#endif // FILE_BLOCKSYSTEM_PILES
