#ifndef FILE_NETGENMESH_PILES
#define FILE_NETGENMESH_PILES

namespace CoupledField
{

class NETGENMesh : public BaseMesh
{
public:
  ///
  NETGENMesh(char * ainputfile);

  ///
  virtual ~NETGENMesh();

  ///
  virtual void Read();
};

}

#endif // FILE_NETGENMESH_PILES
