#ifndef FILE_READMATERIAL_2001
#define FILE_READMATERIAL_2001
 
namespace CoupledField
{

  //! Class for working with material data
  /*!
    This class contains methods for reading information about material from special
    file, which contains data about it. In future we plan to implement here also methods for the interpolation / approximation of nonlinear material behaviour, transformations of material tensors as well as the handling of inhomogeneous material data.
  */
 
class Material
{
public:

  //! Constructor
  Material(const Char * const name);

  //! Deconstructor
  ~Material();

  //! Get density and compressity
  void ReadDensityAndCompress(Double & density, Double & compress);

private:

  //!
  ifstream infile;
  
};
} // end of namespace

#endif
