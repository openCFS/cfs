#ifndef FILE_STANDARDTOPOLOGY
#define FILE_STANDARDTOPOLOGY

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

#endif // FILE_STANDARDTOPOLOGY
