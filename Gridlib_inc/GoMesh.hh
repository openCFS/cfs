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

#include <algorithm>

#include "GbTypes.hh"
#include "GbBitArray.hh"
#include "GoVertex.hh"
#include "GoGeometryElement.hh"
#include "GoEdge.hh"
#include "GbMeshFunctions.hh"
#include "GbBox3.hh"
#include "GbEdgeCollapse.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

// some types to make code more readable
typedef std::vector<GoGeometryElement<float> *> ElementVector;
typedef std::vector<GoVertex<float> *> VertexVector;
typedef std::vector<GoEdge<float> *> EdgeVector;

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
  virtual INLINE void addEdge(GoEdge<float> *edge, int level=0);

  // Query contents
  virtual GbBool has(GbContents f) const;
  virtual GbContents content();

  // Query sizes
  virtual INLINE int getNumElements(int level=0) const;
  virtual INLINE int getNumVertices(int level=0) const;
  virtual INLINE int getNumEdges(int level=0) const;
  virtual INLINE int getNumLevels() const;

  // Query geometry and topology
  virtual GoVertex<float> *         nextVertex(GoVertex<float> *v) = 0;
  virtual GoGeometryElement<float> *nextElement(GoVertex<float> *v) = 0;
  virtual GoEdge<float> *           nextEdge(GoVertex<float> *v) = 0;

  virtual INLINE GoVertex<float> *         getVertex(int i, int level=0) const;
  virtual INLINE GoGeometryElement<float> *getElement(int i, int level=0) const;
  virtual INLINE GoEdge<float> *           getEdge(int i, int level=0) const;
  virtual INLINE GoEdge<float> *           getEdge(GoVertex<float>* v0, GoVertex<float>* v1, int level);

  virtual void getNeighboursAtEdge(int i0, int i1, int level, ElementVector *tlist);
  virtual void getMarkedNeighboursAtEdge(int i0, int i1, int level, GbGeoStatusFlag flag, ElementVector *tlist);

  virtual void getNeighboursAtVertex(int i, int level, ElementVector *tlist);
  virtual void getMarkedNeighboursAtVertex(int i, int level, GbGeoStatusFlag flag, ElementVector *tlist);

  virtual void getVerticesAtVertex(int i, int level, VertexVector *vlist);
  virtual void getMarkedVerticesAtVertex(int i, int level, GbVertexFlag flag, VertexVector *vlist);

  // Dimensions
  virtual GbBox3<float> getAxisAlignedBox (int level=0);
  virtual GbBox3<float> getOrientedBox (int level=0);
  virtual GbBool        getOrientedBox (GbVertexFlag f, GbBox3<float>& rkBox, int level=0);

  // Actions on the mesh
  virtual INLINE void faceNormals() const;
  virtual INLINE void vertexNormals() const;
  virtual void setupNeighbours();
  virtual void setupNeighbours(int level);
  virtual void createEdges();
  virtual void createEdges(int level);
  virtual void setupEdges();
  virtual void switchNormals() = 0;
  virtual int  consistency() const;
  virtual void shrink();

  // Helpers to restore consistent mesh state
  virtual void renumberVertices();
  virtual void renumberVerticesOnLevel(int l);
  virtual void renumberElements();
  virtual void renumberElementsOnLevel(int l);
  virtual void renumberEdges();
  virtual void renumberEdgesOnLevel(int l);
  virtual void updateNeighbourPositions();
  virtual void updateNeighbourPositionsOnLevel(int l) = 0;

  // Utilities
  virtual void removeUnusedVertices();
  virtual void removeUnusedVerticesOnLevel(int l);
//  virtual void removeUnusedEdges();
//  virtual void removeUnusedEdgesOnLevel(int l);
  virtual void statistic() const;

  virtual void smoothUmbrella(int level) {};
  virtual void smoothCurvatureFlow(int level) {};
  virtual void smoothMultiLevelUmbrella(int level) {};


  virtual bool mayCollapseEdge(GoGeometryElement<float> *, int, int);
  virtual std::vector<GoGeometryElement<float> *> *halfEdgeCollapse(GoGeometryElement<float> *, int, int, int &, std::vector<GoVertex<float> *> *);
  virtual void halfEdgeCollapse(int, int);
  virtual int vertexSplit( GbHalfEdgeCollapse &); 


  // convenient iterators over mesh members
  // elements
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

  // vertices
  template <class T> INLINE void
  forAllVertices(const T &t, int level=0) const {
    if (level<vertices_.size())
      std::for_each(vertices_[level]->begin(), vertices_[level]->end(), t);
    else
      errormsg("cannot iterate over vertices level "<<level<<" - valid: [0;"<<vertices_.size()-1<<"]");
  }

  template <class T> INLINE void
  forAllVerticesOnLevel(VertexVector* vV,const T &t) const {
    std::for_each(vV->begin(), vV->end(), t);
  }

  // edges
  template <class T> INLINE void
  forAllEdges(const T &t, int level=0) const {
    if (level<edges_.size())
      std::for_each(edges_[level]->begin(), edges_[level]->end(), t);
    else
      errormsg("cannot iterate over edges level "<<level<<" - valid: [0;"<<edges_.size()-1<<"]");
  }

  template <class T> INLINE void
  forAllEdgesOnLevel(EdgeVector* eV,const T &t) const {
    std::for_each(eV->begin(), eV->end(), t);
  }

  // This operator pretty prints info about the mesh
//  friend std::ostream& operator<<(std::ostream&, const GoMesh&);

protected:

  // internal recursion routines
  virtual void getNeighboursAtEdgeRec(GoGeometryElement<float> *el,
				      GoVertex<float>* v0, GoVertex<float>* v1,
				      ElementVector *tlist);
  virtual void getNeighboursAtVertexRec(GoGeometryElement<float>* firstEl,
					GoVertex<float>* v0,
					ElementVector *tlist);
  virtual void getMarkedNeighboursAtEdgeRec(GoGeometryElement<float>* firstEl,
					    GoVertex<float>* v0, GoVertex<float>* v1,
					    GbGeoStatusFlag flag,
					    ElementVector *tlist);
  virtual void getMarkedNeighboursAtVertexRec(GoGeometryElement<float>* firstEl,
					      GoVertex<float>* v0,
					      GbGeoStatusFlag flag,
					      ElementVector *tlist);

  struct SetupEdges {

    INLINE void swap(EdgeVector &edges, int a, int b) {
      if (a!=b) {
	GoEdge<float> *bufferEdge = edges[a];
	edges[a] = edges[b];
	edges[b] = bufferEdge;
      }
    }

    void operator()(GoGeometryElement<float> *t) {
      // debugmsg("setup Edges for Tetrahedron Id: "<<t->getId());
      int numEdges = t->getNumEdges();
      EdgeVector edges;
      edges.reserve(numEdges);
      int *edgeList;
      int edgeListLen = t->getEdgeList(&edgeList);

      // get edges
      for (int i=0; i<numEdges; ++i) {
	edges[i] = t->getEdge(i);
	// debugmsg("edge "<<i<<" = "<<edges[i]);
      }

      // search each connecting edge
      for (int i=0; i<numEdges-1; ++i) {
	int otherEdge = i;
	GoVertex<float> *t0 = t->getVertex(edgeList[i*2]);
	GoVertex<float> *t1 = t->getVertex(edgeList[i*2+1]);
	while (otherEdge < numEdges) {
	  GoVertex<float> *e0 = edges[otherEdge]->getVertex(0);
	  GoVertex<float> *e1 = edges[otherEdge]->getVertex(1);
	  if ( (e0 == t0) && (e1 == t1) ) break;
	  if ( (e0 == t1) && (e1 == t0) ) break;
	  ++otherEdge;
	}
	if (otherEdge == numEdges) {
	  errormsg("element "<<t->getId()<<": edge not found ?");
	  return;
	}
	swap(edges,i,otherEdge);
      }
    
      // now all edges are at the right position in the edge array

      // write edges in right order
      for(int i=0; i<numEdges; ++i)
	t->setEdge(i, edges[i]);
    }
  };

  // caches for edges and faces

  struct EdgeCacheEntry {
    int idx;
    GoGeometryElement<float> *el;
    GoGeometryElement<float> *el2;
    EdgeCacheEntry *next;
  };

  virtual void updateEdgeCache(int level=0);
  virtual void deleteEdgeCache(int level=0);

  class FaceCacheEntry {
  public:
    virtual ~FaceCacheEntry() {}
    GoGeometryElement<float> *f;
    GoGeometryElement<float> *f2;
    FaceCacheEntry *next;
  };
  
  class FaceCacheTriEntry
    : public FaceCacheEntry
  {
  public:
    virtual ~FaceCacheTriEntry() {}
    int i,j;
  };
 
  class FaceCacheQuadEntry
    : public FaceCacheEntry
  {
  public:
    virtual ~FaceCacheQuadEntry() {}
    int i,j,k;
  };
 
  virtual void updateFaceCache(int level=0) = 0;
  virtual void deleteFaceCache(int level=0);

  // helper routines to construct caches
  void searchTable(int ii, int jj, GoGeometryElement<float>* el, EdgeCacheEntry **table);
  void searchTable(int ii, int jj, int kk, GoGeometryElement<float>* el, FaceCacheEntry **table);
  void searchTable(int ii, int jj, int kk, int ll, GoGeometryElement<float>* el, FaceCacheEntry **table);
  INLINE void initTempFlags(unsigned int num);

  // symbols for contained object types
  GbContents contents_;

  // finally here are the contained objects
  std::vector<VertexVector *>  vertices_;
  std::vector<ElementVector *> elements_;
  std::vector<EdgeVector *>    edges_;
  std::vector<EdgeCacheEntry **> edgeCache_;
  std::vector<FaceCacheEntry **> faceCache_;
  std::vector<int> edgeCacheSize_;
  std::vector<int> faceCacheSize_;
  GbBitArray * tempFlags_;
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
| Revision 1.3  2002/04/03 11:17:08  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.21  2002/03/18 09:57:34  prkipfer
| refactored element structure
|
| Revision 1.20  2001/12/18 13:28:51  prkipfer
| implemented general edge setup method
|
| Revision 1.19  2001/12/11 12:44:59  prkipfer
| fixes for KCC compiler on new PCs
|
| Revision 1.18  2001/09/12 11:53:04  prkipfer
| introduced adaptive tet subdivision
|
| Revision 1.17  2001/08/03 13:04:03  uflabsik
| *** empty log message ***
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
