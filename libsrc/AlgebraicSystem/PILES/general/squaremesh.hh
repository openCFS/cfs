#ifndef FILE_SQUAREMESH_PILES
#define FILE_SQUAREMESH_PILES

namespace CoupledField
{

class SquareMesh : public BaseMesh
{
public:
  ///
  SquareMesh(char * ainputfile);

  ///
  virtual ~SquareMesh();

  ///
  virtual void Read();
};

}

#endif // FILE_SQUAREMESH_PILES
