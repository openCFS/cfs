#ifndef FILE_ANSYSFILE_2002
#define FILE_ANSYSFILE_2002

#include "DataInOut/filetype.hh"
#include "Domain/grid.hh"

namespace CoupledField
{

  //! base class for reading initial data
  /*! 
    Class, that is derived from class FileType for reading mesh-input data, which is produced by Ansys interface. 
  */

class AnsysFile: virtual public FileType 
{

public:

  //! constructor with name of mesh-file
  AnsysFile(const Char * const afilename);

  //! deconstructor
  virtual ~AnsysFile();

   //! read maximum number of nodes in the mesh
  /*!
	\param maxnumnodes maximum number of nodes
  */
  virtual void ReadMaxnumnodes(Integer & maxnumnodes);

   //! read number of save nodes in the mesh
  /*!
	\param nrNodes number of save nodes
  */
  virtual void ReadNumSaveNodes(Integer & nrNodes);
  
  //! read coordinates of nodes of 3d-mesh 
  /*!
	\param coordinates_node returned pointer to array with data
	\param maxnumnodes should be provided number of nodes in mesh
   */
  virtual void ReadCoordinate(Point<3> * const coordinates_node,
                                  const Integer maxnumnodes);

  //! read coordinates of nodes of 2d-mesh
  /*!
        \param coordinates_node returned pointer to array with data
        \param maxnumnodes should be provided number of nodes in mesh
   */
  virtual void ReadCoordinate(Point<2> * const coordinates_node,
			      const Integer maxnumnodes);

# //! read information about elements of the mesh
  /*!
    \param elems out: pointer to vector with elements for each subdomain
    \param orderedElems out: vector with pointers to elements, ordered
    by element numbers
    \param sd vector with color of subdomains, for which elements are read
  */
  void ReadEl(StdVector<Elem*> * elems, 
	      StdVector<Elem*> & orderedElems,
	      const StdVector<std::string> sd); 
  
  //! read 1D-elements. we cause it directly when we set BCs
  /*!
   \param allelems out: pointer to vector with 1D-elements
   \param orderedElems out: vector with pointers to elements, ordered
    by element numbers
   \param sd color of subdomains, for which elements are read
  */
  void ReadEl1d(StdVector<Elem*> * allelems, 
		StdVector<Elem*> & orderedElems,
		const StdVector<std::string> sd);

  //! read 2d - elements from the mesh-file
    /*!
   \param allelems out: pointer to vector with 2D-elements
   \param orderedElems out: vector with pointers to elements, ordered
    by element numbers
   \param sd color of subdomains, for which elements are read
  */
  void ReadEl2d(StdVector<Elem*> * allelems, 
		StdVector<Elem*> & orderedElems,
		const StdVector<std::string> sd);

  //! read 3d -elements from the mesh-file
   /*!
   \param allelems out: pointer to vector with 3D-elements
   \param orderedElems out: vector with pointers to elements, ordered
   by element numbers
   \param sd color of subdomains, for which elements are read
  */
  void ReadEl3d(StdVector<Elem*> * allelems, 
		StdVector<Elem*> & orderedElems,
		const StdVector<std::string> sd);
  
  //! read 3d -elements from the mesh-file and extractes the data for the conf-file
   /*!
   \param sd color of subdomains, for which elements are read
  */
  void ReadEl3dConf(StdVector<std::string> &sd);

  //! read 2d -elements from the mesh-file and extractes the data for the conf-file
   /*!
   \param sd color of subdomains, for which elements are read
  */
  void ReadEl2dConf(StdVector<std::string> &sd);

  //! read 1d -elements from the mesh-file and extractes the data for the conf-file
   /*!
   \param sd color of subdomains, for which elements are read
  */
  void ReadEl1dConf(StdVector<std::string> &sd);

  //! read BCs from the mesh-file and extractes the data for the conf-file
   /*!
   \param sd color of subdomains, for which elements are read
  */
  void ReadBCsConf(StdVector<std::string> &sd);

  //! return dimension of the mesh
  Integer ReadDim();

  //! returns the number of 3D elements
  Integer GetNum3DElems();

  //! returns the number of 2D elements
  Integer GetNum2DElems();

  //! returns the number of 1D elements
  Integer GetNum1DElems();

  //! retuns the number of specified boundary conditions
  Integer GetNumBCs();

  //! retuns the number of specified boundary conditions
  Integer GetNumSaveNodes();

   //! read the mesh from mesh-file for Grid_RG
  /*!
        \param bcs out: vector with global number of nodes which are applied to boundary condition
        \param levels in: vector with color of nodes
  */
  void ReadBCs(std::list<Integer> * bcs, const StdVector<std::string> levels);


 //! read the save nodes
  /*!
        \param saveNodes out: vector with global number of nodes
        \param level in: name of nodes
  */
  virtual void ReadSaveNodes(StdVector<Integer> & saveNodes , const std::string level);
  
  //! read only levels (names) of save nodes
  /*! \param levels out: list with names of save node levels  */
  void ReadLevelOfSaveNodes(StdVector<std::string>& levels);
  

  //!
#ifdef ADAPTGRID
 //! read the mesh from mesh-file for Grid_RG
  /*!
	\param elems out: vector with elements
	\param vertex out: vector with vertices
	\param sd in: vector with color of subdomains which is put in Grid_RG
  */
  void ReadGrid_RG(StdVector<grd::Element*> & elems, StdVector<grd::Vertex*> * vertex, const StdVector<std::string> sd);

 //! read the mesh from mesh-file for Grid_RG
  /*!
        \param elems out: vector with elements
        \param vertex out: vector with vertices
        \param sd in: vector with color of subdomains which is put in Grid_RG
  */
  void ReadBCs_GridRG(StdVector<Integer> & idBCs,StdVector<Integer> &colorBCs);
#endif

protected:
  //! ifstream of input mesh-file
  std::ifstream infile;
  
private:

  //! end position in input mesh-file
  std::string::size_type pos_end;

  //! dimension of problem
  Integer dim_;

  //! number of elems
  Integer maxNumElems_;

  //! maximum number of elements read in so far
  Integer actMaxElemNum_;

  //! number of nodes
  Integer maxNumNodes_;

  //! get a sinlge integer in a save way
  Integer GetInteger(std::string seekexp);

  //! tests the next line for emptyness
  //! \param actPos last position of file pointer before next line
  Boolean IsNextLineEmpty(std::string::size_type actPos);
  
  //! get position after line with seekexp and comments lines
  void getPosLine(const std::string seekexp, std::string::size_type & pos);

  //! get position in line
  void getPosition(const std::string seekexp, std::string::size_type & pos);

  //! read number of nodes for boundary condition  from INFO section of the mesh-file
  void ReadMaxnumnodesbc(Integer & nbc);

  //!  read maximum number of elements from INFO section of the mesh-file
  void ReadMaxnumelem(Integer & , const std::string keyword);

  //! transform type of elem in pointer to base class BaseFE
  BaseFE * Type2ptElem(const Integer itype);

#ifdef ADAPTGRID
  //! read 2d elements from the input mesh-file for the Grid_RG 
  void ReadEl4AdaptGrid2d(StdVector<grd::Element*> & elems, StdVector<grd::Vertex*> * vertices,  const StdVector<std::string> sd);  

  //! read 3d elements from the input mesh-file for the Grid_RG 
  void ReadEl4AdaptGrid3d(StdVector<grd::Element*> & elems, StdVector<grd::Vertex*> * vertices,  const StdVector<std::string> sd);

  //! for each element of the mesh in format of Grid_RG set value: color of subdomain
  /*!
    \param namesd name of the subdomain
    \param sd vector with colors of all subdomains
  */
  void SetNumSD(grd::Element * ptEl, const std::string namesd,  const StdVector<std::string> sd);
#endif

};

}

#endif // FILE_ANSYSFILE
