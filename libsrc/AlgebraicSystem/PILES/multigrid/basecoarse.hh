#ifndef FILE_BASECOARSE_CLA
#define FILE_BASECOARSE_CLA

namespace CoupledField
{

class BaseCoarse
{
public:
  ///
  BaseCoarse(Integer maxnh);

  ///
  virtual ~BaseCoarse();

  ///
  Integer GetSysCSize() {return csize_sys;};

  ///
  Integer GetAuxCSize() {return csize_aux;};

  ///
  Integer GetSysCNNE() {return cnne_sys;};

  ///
  Integer GetAuxCNNE() {return cnne_aux;};

  ///
  virtual void SetNeighbour(Double alpha, Double epsmat) = 0;

  ///
  virtual void SetSize() = 0;

  ///
  void Calc(Integer * rsw);

  ///
  void Print() const;

  ///
  virtual void CalcTopology() = 0;

  BaseTopology * GetTopology() const {return topology;};

protected:
  ///
  Integer fsize_sys, fsize_aux, fnne_sys, fnne_aux;

  ///
  Integer csize_sys, csize_aux, cnne_sys, cnne_aux;

  ///
  Integer startpoint, maxnh;

  ///
  BaseMatrix * mat;

  ///
  BaseTopology * topology;

};

}

#endif // FILE_BASECOARSE_CLA
