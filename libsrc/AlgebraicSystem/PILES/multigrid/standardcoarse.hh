#ifndef FILE_STANDARDCOARSE_CLA
#define FILE_STANDARDCOARSE_CLA

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

#endif // FILE_STANDARDCOARSE_CLA
