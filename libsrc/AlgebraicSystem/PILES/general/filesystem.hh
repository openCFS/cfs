#ifndef FILE_FILESYSTEM_PILES
#define FILE_FILESYSTEM_PILES

namespace CoupledField
{

class FileSystem
{
public:
  ///
  FileSystem(char * ainputfile);

  ///
  virtual ~FileSystem();

private:
  ///
  char * filename;

};

}

#endif // FILE_FILESYSTEM_PILES
