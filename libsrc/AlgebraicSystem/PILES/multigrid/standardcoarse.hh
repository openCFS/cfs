#ifndef FILE_STANDARDCOARSE_PILES
#define FILE_STANDARDCOARSE_PILES

namespace CoupledField
{

class StandardCoarse : public BaseCoarse
{
public:
  ///
  StandardCoarse(BaseMatrix * amat, Integer amaxnh);

  ///
  virtual ~StandardCoarse();

  ///
  virtual void SetNeighbour(Double alpha, Double epsmat);

  ///
  virtual void CalcTopology();

  ///
  virtual void SetSize();
};

}

#endif // FILE_STANDARDCOARSE_PILES
