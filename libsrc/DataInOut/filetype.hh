#ifndef FILE_FILETYPE_2001
#define FILE_FILETYPE_2001

#include <Utils/vector.hh>
#include <list>
#include <vector>
#include <string>

#ifdef ADAPTGRID
#include "Element.h"  
#endif

namespace CoupledField
{

struct Elem;

  //! base class for reading initial data
  /*! 
    Base class for reading mesh-data. Functions of the class are virtual in order to handle the different types of input files.
  */
class FileType
{

public:

  //! constructor with name of mesh-file
  FileType(const Char * const afilename);

  //! deconstructor
  virtual ~FileType();

  //! read maximum number of nodes in the mesh
  /*!
	\param maxnumnodes maximum number of nodes
  */
  virtual void ReadMaxnumnodes(Integer & maxnumnodes)=0;  

  //! read coordinates of nodes of 3d-mesh 
  /*!
	\param coordinates_node returned pointer to array with data
	\param maxnumnodes should be provided number of nodes in mesh
   */
  virtual void ReadCoordinate(Point<3> * const coordinates_node,                                                     const Integer maxnumnodes)=0;

  //! read coordinates of nodes of 2d-mesh
  /*!
        \param coordinates_node returned pointer to array with data
        \param maxnumnodes should be provided number of nodes in mesh
   */
  virtual void ReadCoordinate(Point<2> * const coordinates_node,
			      const Integer maxnumnodes)=0;

  //! read boundary conditions
  /*!
	\param bcs  out: returned pointer to list with global number of nodes to which boundary conditions are applied. 
	\param colors in: vector with colors of nodes, which are requested
  */
  virtual void ReadBCs(std::list<Integer> * bcs, std::vector<std::string> colors)=0;  

  //! read information about elements of the mesh
  /*!
	\param elems out: pointer to vector with elements
	\param sd vector with color of subdomains, for which elements are read
  */
  virtual void ReadEl(std::vector<Elem*> * elems, const std::vector<std::string>
sd)
 { Error(" not implemented",__FILE__,__LINE__);}

  //! read 1D element. we call it directly when we set BCs
  /*!
   \param allelems out: pointer to vector with 1D-elements
   \param sd color of subdomains, for which elements are read
  */
  virtual void ReadEl1d(std::vector<Elem*> * allelems, const std::vector<std::string> sd)
 { Error(" not implemented",__FILE__,__LINE__);}  

 //! read 2d - elements from the mesh-file
    /*!
   \param allelems out: pointer to vector with 3D-elements
   \param sd color of subdomains, for which elements are read
  */
  virtual void ReadEl2d(std::vector<Elem*> * allelems, const std::vector<std::string> sd)
 { Error(" not implemented",__FILE__,__LINE__);}  

//! read 3d - elements from the mesh-file
    /*!
   \param allelems out: pointer to vector with 3D-elements
   \param sd color of subdomains, for which elements are read
  */
  virtual void ReadEl3d(std::vector<Elem*> * allelems, const std::vector<std::string> sd)
 { Error(" not implemented",__FILE__,__LINE__);} 

  //! read 3d -elements from the mesh-file and extractes the data for the conf-file
   /*!
   \param sd color of subdomains, for which elements are read
  */
   virtual void ReadEl3dConf(std::vector<std::string> &sd)
 { Error(" not implemented",__FILE__,__LINE__);} 

  //! read 2d -elements from the mesh-file and extractes the data for the conf-file
   /*!
   \param sd color of subdomains, for which elements are read
  */
   virtual void ReadEl2dConf(std::vector<std::string> &sd)
  { Error(" not implemented",__FILE__,__LINE__);} 

   //! read 1d -elements from the mesh-file and extractes the data for the conf-file
   /*!
   \param sd color of subdomains, for which elements are read
  */
   virtual void ReadEl1dConf(std::vector<std::string> &sd)
  { Error(" not implemented",__FILE__,__LINE__);} 

  //! read BCs from the mesh-file and extractes the data for the conf-file
   /*!
   \param sd color of subdomains, for which elements are read
  */
  virtual void ReadBCsConf(std::vector<std::string> &sd)
  { Error(" not implemented",__FILE__,__LINE__);} 

#ifdef ADAPTGRID
  //! read the mesh from mesh-file for Grid_RG
  /*!
	\param elems out: vector with elements
	\param vertex out: vector with vertices
	\param sd in: vector with color of subdomains which is put in Grid_RG
  */
  virtual void ReadGrid_RG(std::vector<grd::Element*> & elems, std::vector<grd::Vertex*> * vertex, const std::vector<std::string> sd)
   { Error(" not implemented",__FILE__,__LINE__);}

  //! read the mesh from mesh-file for Grid_RG
  /*!
        \param elems out: vector with elements
        \param vertex out: vector with vertices
        \param sd in: vector with color of subdomains which is put in Grid_RG
  */
  virtual void ReadBCs_GridRG(std::vector<Integer> & idBCs,std::vector<Integer> &colorBCs)
   { Error(" not implemented",__FILE__,__LINE__);}
#endif
 
  //! return dimension of the mesh
  virtual Integer ReadDim()
  { Error(" not implemented",__FILE__,__LINE__);} 

  //! returns the number of 3D elements
  virtual Integer GetNum3DElems()
  { Error(" not implemented",__FILE__,__LINE__);} 

  //! returns the number of 2D elements
  virtual Integer GetNum2DElems()
  { Error(" not implemented",__FILE__,__LINE__);} 

  //! returns the number of 1D elements
  virtual Integer GetNum1DElems()
  { Error(" not implemented",__FILE__,__LINE__);} 

  //! retuns the number of specified boundary conditions
  virtual Integer GetNumBCs()
  { Error(" not implemented",__FILE__,__LINE__);} 

protected:

  //! name of input file
  Char * filename;

};

}


#endif // FILE_FILETYPE
