/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOMESH_HH
#define GOMESH_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <vector>
#include <iostream>
#include <typeinfo>
#include <algorithm>

#include "GoVertex.hh"
#include "GoGeometryElement.hh"
#include "GoTetrahedronElement.hh"
#include "GbMeshFunctions.hh"
#include "GbBox3.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

//! class GoMesh
/*!
  Base class for 2d or 3d meshes. 
 */
class GoMesh 
{
public:
  // Constructor
  GoMesh();
  virtual ~GoMesh();

  // Add objects to the mesh
  virtual INLINE void addElement(GoGeometryElement<float> *face, int level=0);
  virtual INLINE void addVertex(GoVertex<float> *vert, int level=0);

  // Query contents
  virtual GbBool has(GbContents f);
  virtual GbContents content();

  // Query sizes
  virtual INLINE int getNumElements(int level=0) const;
  virtual INLINE int getNumVertices(int level=0) const;
  virtual INLINE int getNumLevels() const;
  virtual GbBox3<float> getAxisAlignedBox ();
  virtual GbBox3<float> getOrientedBox ();
  virtual GbBool getOrientedBox (GbVertexFlag f, GbBox3<float>& rkBox);

  // Query topology
  virtual GoVertex<float> *nextVertex(GoVertex<float> *v) = 0;
  virtual GoGeometryElement<float>   *nextElement(GoVertex<float> *v) = 0;
  virtual INLINE GoVertex<float> *getVertex(int i, int level=0) const;
  virtual INLINE GoGeometryElement<float>   *getElement(int i, int level=0) const;

  // Actions on the mesh
  virtual INLINE void faceNormals() const;
  virtual INLINE void vertexNormals() const;
  virtual void setupNeighbours() = 0;
  virtual void switchNormals() = 0;
  virtual int  consistency();
  virtual void shrink();
  virtual void renumberVertices();
  virtual void removeUnusedVertices();
  virtual void statistic();

  virtual std::vector<GoGeometryElement<float> *> *halfEdgeCollapse(GoGeometryElement<float> *, int, int, int &, std::vector<GoVertex<float> *> *);
  // virtual void halfEdgeCollapse(int, int);

  virtual std::vector<GoGeometryElement<float> *> * getAllNeighbourElements(GoVertex<float>* vertex);


  // some types to make code more readable
  typedef std::vector<GoGeometryElement<float> *>   ElementVector;
  typedef std::vector<GoVertex<float> *> VertexVector;

  // convenient iterators over mesh members
  template <class T> INLINE void
  forAllElements(const T &t, int level=0) const {
    if (level<elements_.size())
      std::for_each(elements_[level]->begin(), elements_[level]->end(), t);
    else
      errormsg("cannot iterate over elements level "<<level<<" - valid: [0;"<<elements_.size()-1<<"]");
  }

  template <class T> INLINE void
  forAllElementsOnLevel(ElementVector* eV,const T &t) const {
    std::for_each(eV->begin(), eV->end(), t);
  }

  template <class T> INLINE void
  forAllVertices(const T &t, int level=0) const {
    if (level<vertices_.size()) {
//      for (int i=0; i<=level; i++)
	std::for_each(vertices_[level]->begin(), vertices_[level]->end(), t);
    }
    else
      errormsg("cannot iterate over vertices level "<<level<<" - valid: [0;"<<vertices_.size()-1<<"]");
  }

  template <class T> INLINE void
  forAllVerticesOnLevel(VertexVector* vV,const T &t) const {
    std::for_each(vV->begin(), vV->end(), t);
  }

  // This operator pretty prints info about the mesh
//  friend std::ostream& operator<<(std::ostream&, const GoMesh&);

protected:

  GbContents contents_;

  std::vector<VertexVector *>  vertices_;
  std::vector<ElementVector *> elements_;

private:

  // Inhibit assignment to a mesh
  GoMesh& operator=(const GoMesh& rhs);
};

#ifndef OUTLINE
#include "GoMesh.in"
#endif

#endif // GOMESH_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
|
| Revision 1.16  2001/06/15 08:55:57  prkipfer
| added compiler hints
|
| Revision 1.15  2001/05/03 08:05:51  uflabsik
| functionality added
|
| Revision 1.14  2001/01/02 14:29:55  prkipfer
| new bounding box computation now uses new GbBox3 class
|
| Revision 1.13  2000/12/12 10:03:49  prkipfer
| moved consistency check to GoMesh
|
| Revision 1.12  2000/11/07 17:27:54  prkipfer
| introduced external polymorphism for vertices
|
| Revision 1.11  2000/11/07 09:56:18  uflabsik
| docu added
|
| Revision 1.10  2000/10/31 16:29:06  uflabsik
| new version
|
| Revision 1.9  2000/10/23 09:35:19  prkipfer
| moved bounding box calculation to GoMesh
|
| Revision 1.8  2000/09/19 07:33:59  uflabsik
| modified
|
| Revision 1.7  2000/09/07 16:55:50  prkipfer
| added subdivision
|
| Revision 1.6  2000/08/28 09:25:25  prkipfer
| added localized subdivision mechanism
|
| Revision 1.5  2000/08/22 17:15:23  prkipfer
| added possibility to query mesh contents
|
| Revision 1.4  2000/07/21 10:40:47  prkipfer
| dropped g++ support
|
| Revision 1.3  2000/07/20 13:32:41  prkipfer
| complete rework: now KCC is the standard compiler for Linux
|
| Revision 1.2  2000/07/03 16:15:12  prkipfer
| extended IO subsystem, created class registry, added tetrahedron
|
| Revision 1.1  2000/06/16 13:32:21  prkipfer
| started mesh class hierarchy
|
|
+---------------------------------------------------------------------*/
