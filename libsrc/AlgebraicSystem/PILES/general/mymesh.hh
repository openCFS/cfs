#ifndef FILE_MYMESH_PILES
#define FILE_MYMESH_PILES

namespace CoupledField
{

class MyMesh : public BaseMesh
{
public:
  ///
  MyMesh(char * ainputfile);

  ///
  virtual ~MyMesh();

  ///
  virtual void Read();
};

}

#endif // FILE_MYMESH_PILES
