#ifndef FILE_BEMMESH_PILES
#define FILE_BEMMESH_PILES

namespace CoupledField
{

class BEMMesh : public BaseMesh
{
public:
  ///
  BEMMesh(char * ainputfile);

  ///
  virtual ~BEMMesh();

  ///
  virtual void Read();
};

}

#endif // FILE_BEMMESH_PILES
