#ifndef FILE_STANDARDTOPOLOGY_CLA
#define FILE_STANDARDTOPOLOGY_CLA

namespace CoupledField
{

class StandardTopology : public BaseTopology
{
public:
  ///
  StandardTopology(Integer asize, Integer anne, Integer * pos, Integer * start);

  ///
  virtual ~StandardTopology();

  ///
  virtual void CalcCoarseGraph();
  
};

}

#endif // FILE_STANDARDTOPOLOGY_CLA
