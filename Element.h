// - C++ -
/***************************************************************************
    File        : Element.h
    Description : 

 ---------------------------------------------------------------------------
    Begin       : Fri Dec 7 2001
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#ifndef ELEMENT_H
#define ELEMENT_H


// Grid
#include "Vertex.h"
#include "Edge.h"
#include "Jacobian.h"


namespace grd {

#define GRD_TRIANGLE    1
#define GRD_QUADRANGLE  2
#define GRD_TETRAHEDRON 3
#define GRD_OCTAHEDRON  4
#define GRD_HEXAHEDRON  5
#define GRD_PRISM       6
#define GRD_PYRAMID     7



// Adaptive Refinement
class GridLevel;
// Conforming closure
class ConformingClosure;


class Element {
public:
  // Constructors
  Element() {
    parent = 0;
    child  = 0;
    // all new elements are regular
    // due to virtual refinement
    tag    = 4; //ELEMENT_TYPE_REGULAR
    value  = 0;
    iso    = false;
    minVal = 0.0;
    maxVal = 0.0;
  }

  // Destructor
  virtual ~Element() { }

  // Element state
  inline bool isRegular();
  inline bool isRefined();
  inline bool isIrregular();

  // Element marks manipulation
  inline void markForRefinement();
  inline void markForCoarsening();
  inline void markParentForRefinement();
  inline void markParentForCoarsening();
  inline bool isMarked();
  inline bool isMarkedForRefinement();
  inline bool isMarkedForCoarsening();


  // Element refinement
  virtual void refine(GridLevel& gridLevel) = 0;
  virtual void coarsen(GridLevel& gridLevel) = 0;

  inline  void  close();
  virtual void  close(ConformingClosure& closure) = 0;
  virtual void  close(int refP,ConformingClosure& closure) = 0;

  virtual int getRefinementPattern() = 0;

  inline Element* getParent();
  inline Element* getChild(int ch);

  virtual int getNoOfChildren() = 0;

  virtual bool queryChildrenRefinement() = 0;
  virtual bool queryEdgeRefinement() = 0;
  virtual bool queryChildrenEdgeRefinement() = 0;
  virtual bool queryChildrenMarks() = 0;

  virtual int  getCommonFace(Element* theElem) = 0;
  virtual void getRefElements(int f,Vertex** theVertex,Edge** theEdg,Element** theChild,int* theChildFace) = 0;

  inline void resetDescendant();

  // Topology
  virtual int getNoOfVertices() = 0;
  virtual int getNoOfEdges() = 0;
  virtual int getNoOfFaces() = 0;

  virtual void setVertex(int vertex,Vertex* theVertex) = 0;
  virtual void setEdge(int edge,Edge* theEdge) = 0;
  virtual void setNeighbor(int face,Element* theNeighbor) = 0;

  virtual Vertex*  getVertex(int vertex) = 0;
  virtual Edge*    getEdge(int edge) = 0;
  virtual Element* getNeighbor(int face) = 0;

  virtual void setVertices(Vertex** vt) = 0;
  virtual void getVertices(Vertex** vt) = 0;

  virtual void setEdges(Edge** edg) = 0;
  virtual void getEdges(Edge** edg) = 0;

  // Topology
  virtual int type() = 0;
  virtual int noOfVertexAtFace(int fc) = 0;
  virtual Vertex* getVertexAtFace(int fc,int vt) = 0;
  virtual Vertex* getVertexAtEdge(int ed,int vt) = 0;
  virtual Edge*   getEdgeAtFace(int fc,int ed)   = 0;

//  // Geometry
//  virtual void getGlobalCoords(const double* lc,double* coord) = 0;
//  virtual void getJacobian(Jacobian& Jac) = 0;

  // Element internal value
  void setValue(int val);
  int  getValue();

  // Element refinement
  inline void setRefined();
  inline void setRegular();
  inline void setIrregular();
  inline void setParent(Element* theParent);

  // For Iso-surface
  void setMinVal(float mv) { minVal = mv; }
  void setMaxVal(float mv) { maxVal = mv; }
  float getMinVal() { return minVal; }
  float getMaxVal() { return maxVal; }

  // iso-surface
  void setIsoRefined()   { iso = true;  }
  void unsetIsoRefined() { iso = false; }
  bool isIsoRefined()    { return iso;  }

protected:
  // tree structure
  Element*  parent;
  Element** child;

private:
  // mark and element type; bit wise
  char tag;
  // internal element information
  char value;

  // iso-surface
  bool iso;
  // only for iso-surfaces
  float minVal;
  float maxVal;
};




// ____________________________________________________________
//
//	
//	         inline functions
//
//



// Note: element states and marks
// -------------------
// the variable tag is used bit wise
// Element Marks:
//   element can be marked for refinement or coarsening
// Element States:
//  (i)   an element can be regular, if it is not regularly or
//        irregularly refined. This is the default
//  (ii)  an element can be irregular, if one of its neighbors
//        is refined
//  (iii) an element can be refined, if it is regularly refined
//
// Bit usage:
// bit correspondence:	00 00 00 00
//	first two bits for element mark
//		00 00 00 00	not marked
// 		00 00 00 01	marked for refinement
//		00 00 00 10 	marked for coarsenig
// 	third, fourth and fifth bits for element state
//		00 00 01 xx	regular, default
//    00 00 1x xx     irregular
//		00 01 xx xx	regular refined


#define ELEMENT_MARK_REFINE    1
#define ELEMENT_MARK_COARSEN   2
#define ELEMENT_TYPE_REGULAR   4
#define ELEMENT_TYPE_IRREGULAR 8
#define ELEMENT_TYPE_REFINED   16


// Description
//
inline bool
Element::isRegular()
{
  if (tag & ELEMENT_TYPE_REGULAR)
    return true;
  else
    return false;
  
}


// Description
// 1 if regular refined, 0 if not or irregular refined
inline bool
Element::isRefined()
{
  if (tag & ELEMENT_TYPE_REFINED)
    return true;
  else
    return false;
}


// Description
//
inline bool
Element::isIrregular()
{
  if (tag & ELEMENT_TYPE_IRREGULAR)
    return true;
  else
    return false;
}


// Description
//
inline void
Element::markForRefinement()
{
  if (tag & ELEMENT_MARK_COARSEN)
    tag = (char)(tag ^ ELEMENT_MARK_COARSEN);
  tag = (char)(tag | ELEMENT_MARK_REFINE);
}


// Description
//
inline void
Element::markForCoarsening()
{
  if (tag & ELEMENT_MARK_REFINE)
    tag = (char)(tag ^ ELEMENT_MARK_REFINE);
  tag = (char)(tag | ELEMENT_MARK_COARSEN);
}


// Description
//
inline void
Element::markParentForRefinement()
{
  if (parent)
    parent->markForRefinement();
}


// Description
//
inline void
Element::markParentForCoarsening()
{
  if (parent)
    parent->markForCoarsening();
}


// Description
//
inline bool
Element::isMarked()
{
  if (tag & ELEMENT_MARK_REFINE)
    return true;
  else if (tag & ELEMENT_MARK_COARSEN)
    return true;
  else
    return false;
}


// Description
//
inline bool
Element::isMarkedForRefinement()
{
  if (tag & ELEMENT_MARK_REFINE)
    return true;
  else
    return false;
}


// Description
//
inline bool
Element::isMarkedForCoarsening()
{
 if (tag & ELEMENT_MARK_COARSEN)
   return true;
 else
   return false;
}


// Description
// for virtual close set element to irregular
inline void
Element::close()
{
  tag = ELEMENT_TYPE_IRREGULAR;
}


// Description
//
inline Element*
Element::getParent()
{
  return parent;
}


// Description
//
inline Element*
Element::getChild(int ch)
{
  return child[ch];
}


// Description
// for removal purposes
inline void
Element::resetDescendant()
{
  child = 0;
  tag = ELEMENT_TYPE_REGULAR;
}


//
inline void
Element::setValue(int val)
{
  value = (char)val;
}


// Description
//
inline int
Element::getValue()
{
  return value;
}



//-----------------------------------------------------
//
// Protected function for derived classes only
//
//-----------------------------------------------------

// Description
// if refined set all marks to 0 and set the
// return bit on
inline void
Element::setRefined()
{
  tag = ELEMENT_TYPE_REFINED;
}


// Description
//
inline void
Element::setRegular()
{
  tag = ELEMENT_TYPE_REGULAR;
  child = 0;
}


// Description
//
inline void
Element::setIrregular()
{
  tag = ELEMENT_TYPE_IRREGULAR;
}

// Description
//
inline void
Element::setParent(Element* theParent)
{
  parent = theParent;
}

} // namespace grd

#endif // ELEMENT_H
