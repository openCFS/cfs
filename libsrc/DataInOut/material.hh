#ifndef FILE_READMATERIAL_2001
#define FILE_READMATERIAL_2001

#include "tools.hh"
 
namespace CoupledField
{

  //! Class for working with material data
  /*!
    This class contains methods for the reading information about material from the 
    input a file with material's data.
  */
 
class Material
{
public:

  //! Constructor
  Material(const Char * const name);

  //! Deconstructor
  ~Material();

  //! get density and compressity for acoustic equation from material data-file
   /*!
	\param density out: value of density from the material file
	\param compress out: value of compressity from the material file
	\param matnum material number
	\param keyword name of material in the  material file
  */
  void ReadDensityAndCompressity(Double & density, Double & compress, const Integer matnum, const std::string keyword);

  //! get dielectric term for an electrostatic equation from the material file
  /*!
	\param dielectr out: value of dielectric term from the material file
	\param matnum in: material number
  */
  void ReadDielectricTerms(Double & dielectr,const Integer matnum);

private:

  //! ifstream which is associated with the material file
  std::ifstream infile;
  
};
} // end of namespace

#endif
