#ifndef FILE_NEUROFEMMESH_PILES
#define FILE_NEUROFEMMESH_PILES

namespace CoupledField
{

class NeuroFEMMesh : public BaseMesh
{
public:
  ///
  NeuroFEMMesh(char * ainputfile);

  ///
  virtual ~NeuroFEMMesh();

  ///
  virtual void Read();
};

}

#endif // FILE_NEUROFEMMESH_PILES
