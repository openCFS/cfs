#ifndef FILE_BASETOPOLOGY_CLA
#define FILE_BASETOPOLOGY_CLA

namespace CoupledField
{

class BaseTopology
{
public:
  ///
  BaseTopology(Integer asize, Integer anne, Integer * pos, Integer * start);

  ///
  virtual ~BaseTopology();

  ///
  Boolean MemberS1(Integer p, Integer q);

  ///
  void UpdateS1(Integer p, Integer q);

  ///
  void UpdateS2(Integer p, Integer q);

  ///
  void SetDirichlet(Integer p) {ss1[p-1] = 0; ss2[p-1] = 0;};

  ///
  Integer GetNHSize(Integer i) const {return (snh[i]-snh[i-1]);};

  ///
  Integer GetS1Size(Integer i) const {return ss1[i-1];};

  ///
  Integer GetS2Size(Integer i) const {return ss2[i-1];};

  ///
  Integer GetNHNode(Integer i, Integer j) const {return nh[snh[i-1]+j-1];};

  ///
  Integer GetS1Node(Integer i, Integer j) const {return s1[snh[i-1]+j-1];};

  ///
  Integer GetS2Node(Integer i, Integer j) const {return s2[snh[i-1]+j-1];};

  ///
  void SetCoarseGridPoint(Integer p) {csize++; cnn[p-1] = csize;}

  ///
  void SetFineGridPoint(Integer p);

  ///
  virtual void CalcCoarseGraph() = 0;

  ///
  Integer GetNextCoarseGridPoint(Integer p, Integer maxnh) const;

  ///
  Integer GetNextStartPoint() const;

  ///
  Integer GetStrongNumCoarse(Integer p) const;

  ///
  Integer GetNumCoarse(Integer p) const;

  ///
  Integer GetNumFine(Integer p) const;

  ///
  Integer GetCoarseSize() const {return csize;};

  ///
  Integer GetCoarseNNE() const {return cnne;};

  ///
  Integer MaxNumNeighbour() const;

  ///
  Integer GetStrongCoarseNeighbour(Integer p, Integer * cnh);

  ///
  Integer GetCoarseNeighbour(Integer p, Integer * cnh);

  ///
  Integer GetFineNeighbour(Integer p, Integer * fnh);

  ///
  Integer GetCNN(Integer p) const {return cnn[p-1];};

  ///
  Integer GetCoarseRowSize(Integer p) const {return nrs[p-1];};

  ///
  Integer * GetCoarsePos(Integer p) {return &crs[(p-1)*offset];};

  ///
  void SetRSW(Integer * rsw);

  ///
  void Insert(Integer p, Integer q);

  ///
  void Print() const;

protected:
  ///
  Integer size, nne, csize, cnne, *cnn, * nh, * s1, *s2, * snh, * ss1, * ss2, *active;

  ///
  Integer * crs, * nrs, offset;
  
};

class LocalSet
{
public:
  ///
  LocalSet(Integer asize, Integer adof=1);

  ///
  virtual ~LocalSet();

  ///
  void Init();

  ///
  Integer GetSize() const {return acts;};

  ///
  void Insert(Integer q);

  ///
  void Insert(Integer q, Double val);

  ///
  void Insert(Integer q, DenseMatrix * val);

  ///
  void Sort();

  ///
  Integer Pos(Integer q) const {return p[q-1];};

  ///
  Double Get(Integer q) const {return v[q-1];};

  ///
  DenseMatrix * GetBlock(Integer q) const 
    {*dm = &v[(q-1)*dof*dof]; return dm;};

  ///
  void Print() const;

private:
  ///
  Double * v;

  ///
  Integer * p;

  ///
  Integer size, acts, dof;

  ///
  DenseMatrix * dm;
};

}

#endif // FILE_BASETOPOLOGY_CLA
