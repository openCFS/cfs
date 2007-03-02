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
#include "Vertex.hh"
#include "Edge.hh"
#include "Jacobian.hh"


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


class Element {
public:
  // Constructors
  Element() : parent(0),child(0),chcount(0),tag(0),value(0) {}

  // Destructor
  virtual ~Element() { }

  // Element type, use the first three bits
  inline bool isRegular() const;
  inline bool isRefined() const;
  inline bool isIrregular() const;
  inline bool isClosure() const;

  inline void setRegular();
  inline void setRefined();
  inline void setIrregular();
  inline void setClosure();

  // Element marks, uses the fourth and fifth bits
  inline bool isMarkedForRefinement() const;
  inline bool isMarkedForCoarsening() const;

  inline void markForRefinement();
  inline void markForCoarsening();
  inline void markParentForRefinement();
  inline void markChildrenForCoarsening();
  inline void unmarkForRefinement();
  inline void unmarkForCoarsening();
  inline void unsetAllRefinementMarks();
  inline void unsetChildrenRefinementMarks();

  // Reference for DB operations, use sixth bit
  inline void ref();
  inline void unref();
  inline bool queryRef();

  // Mark for general purpose, seventh bit
  inline void mark();
  inline void unmark();
  inline bool queryMark();

  // Set element marks to default, i.e. 0
  inline void resetAllMarks();

  // Element refinement
  virtual void refine(GridLevel& gridLevel)  = 0;
  virtual void coarsen(GridLevel& gridLevel) = 0;
  virtual void close(GridLevel& gridLevel)   = 0;

  inline void close();

  virtual int getRefinementPattern() const = 0;

  inline Element* getParent()      const;
  inline void setParent(Element* theParent);

  inline int      getNoOfChildren() const;
  inline Element* getChild(int ch)  const;
  inline void     resetPointerToChildren();

  virtual bool queryEdgeRefinement()         const = 0;
  virtual bool queryChildrenEdgeRefinement() const = 0;

  inline bool queryChildrenRefinement()     const;
  inline bool queryChildrenMarks()          const;

  virtual int  getCommonFace(Element* theElem) const = 0;
  virtual void getRefElements(int f,Vertex** theVertex,Edge** theEdg,Element** theChild,int* theChildFace) const = 0;

  // Topology
  virtual int getNoOfVertices() const = 0;
  virtual int getNoOfEdges() const = 0;
  virtual int getNoOfFaces() const = 0;

  virtual void setVertex(int vertex,Vertex* theVertex) = 0;
  virtual void setEdge(int edge,Edge* theEdge) = 0;
  virtual void setNeighbor(int face,Element* theNeighbor) = 0;

  virtual Vertex*  getVertex(int vertex) const = 0;
  virtual Edge*    getEdge(int edge) const = 0;
  virtual Element* getNeighbor(int face) const = 0;

  virtual void setVertices(Vertex* const* vt) = 0;
  virtual void getVertices(Vertex** vt) const = 0;

  virtual void setEdges(Edge** edg) = 0;
  virtual void getEdges(Edge** edg) const = 0;

  // Topology
  virtual int     type() const = 0;
  virtual int     noOfVerticesAtFace(int fc) const = 0;
  virtual int     topologicalVertexAtFace(int fc,int vt) const = 0;
  virtual Vertex* getVertexAtFace(int fc,int vt) const = 0;
  virtual int     noOfEdgesAtFace(int fc) const = 0;
  virtual Vertex* getVertexAtEdge(int ed,int vt) const = 0;
  virtual Edge*   getEdgeAtFace(int fc,int ed)   const = 0;

  // Element internal value
  void setValue(int val);
  int  getValue() const;

protected:
  // tree structure
  Element*  parent;
  Element** child;
  unsigned char chcount;

private:
  // mark and element type; bit wise
  unsigned char tag;
  // internal element information
  unsigned char value;
};




// ____________________________________________________________ //
//                                                              //
//                                                              //
//         INLINE FUNCTIONS                                     //
//                                                              //
// ____________________________________________________________ //



// Note: element states and marks
// -------------------
// the variable tag is used bit wise
// Element Type (state):
//  (i)   an element can be regular, if it is not regularly or
//        irregularly refined. This is the default
//  (ii)  an element can be refined, if it is regularly refined
//  (iii) an element can be irregular, if one or more of its neighbors
//        is regularly refined
//  (iv)  an element can be closure, i.e. for conforming closure
//        result of an irregular refinement
// Element Marks:
//   element can be marked for refinement or coarsening
// Element Reference
//   element can be referenced and unreferenced
// Element general marks
//   element can be marked for a general purpose
//
// Bit correspondence:  00 00 00 00
//  first two bits for element type
//    00 00 00 00 regular
//    00 00 00 01 refined
//    00 00 00 10 irregular
//    00 00 01 00 closure
//  fourth and fifth bits for element refinement
//    00 00 0x xx nor marked
//    00 00 1x xx marked for refinement
//    00 01 0x xx marked for coarsening
//  sixth bit for reference: element is referenced
//    00 0x xx xx not reference
//    00 1x xx xx reference
//  sevenh bit element is mark for some general purpose
//   00 xx xx xx not marked
//   01 xx xx xx marked

#define ELEMENT_TYPE_REGULAR   0
#define ELEMENT_TYPE_REFINED   1
#define ELEMENT_TYPE_IRREGULAR 2
#define ELEMENT_TYPE_CLOSURE   4

#define ELEMENT_TYPE_NEG       7

#define ELEMENT_MARK_REFINE    8
#define ELEMENT_MARK_COARSEN   16

#define ELEMENT_MARK_REFERENCE 32
#define ELEMENT_MARK_GENERAL   64



// Implementation
inline bool
Element::isRegular() const
{
  if (tag&ELEMENT_TYPE_NEG)
    return false;
  else
    return true;
}

inline bool
Element::isRefined() const
{
  return (tag&ELEMENT_TYPE_REFINED);
}

inline bool
Element::isIrregular() const
{
  return (tag&ELEMENT_TYPE_IRREGULAR);
}

inline bool
Element::isClosure() const
{
  return (tag&ELEMENT_TYPE_CLOSURE);
}


inline void
Element::setRegular()
{
  tag |= ELEMENT_TYPE_REGULAR;
  if (tag&ELEMENT_TYPE_REFINED)
    tag ^=ELEMENT_TYPE_REFINED;
  if (tag&ELEMENT_TYPE_IRREGULAR)
    tag ^=ELEMENT_TYPE_IRREGULAR;
  if (tag&ELEMENT_TYPE_CLOSURE)
    tag ^=ELEMENT_TYPE_CLOSURE;
}

inline void
Element::setRefined()
{
  tag |= ELEMENT_TYPE_REFINED;
  if (tag&ELEMENT_TYPE_REGULAR)
    tag ^=ELEMENT_TYPE_REGULAR;
  if (tag&ELEMENT_TYPE_IRREGULAR)
    tag ^=ELEMENT_TYPE_IRREGULAR;
  if (tag&ELEMENT_TYPE_CLOSURE)
    tag ^=ELEMENT_TYPE_CLOSURE;
}

inline void
Element::setIrregular()
{
  tag |=ELEMENT_TYPE_IRREGULAR;
  if (tag&ELEMENT_TYPE_REGULAR)
    tag ^=ELEMENT_TYPE_REGULAR;
  if (tag&ELEMENT_TYPE_REFINED)
    tag ^= ELEMENT_TYPE_REFINED;
  if (tag&ELEMENT_TYPE_CLOSURE)
    tag ^=ELEMENT_TYPE_CLOSURE;
}

inline void
Element::setClosure()
{
  tag |=ELEMENT_TYPE_CLOSURE;
  if (tag&ELEMENT_TYPE_REGULAR)
    tag ^=ELEMENT_TYPE_REGULAR;
  if (tag&ELEMENT_TYPE_REFINED)
    tag ^= ELEMENT_TYPE_REFINED;
  if (tag&ELEMENT_TYPE_IRREGULAR)
    tag ^=ELEMENT_TYPE_IRREGULAR;
}


inline bool
Element::isMarkedForRefinement() const
{
  return (tag&ELEMENT_MARK_REFINE);
}


inline bool
Element::isMarkedForCoarsening() const
{
  return (tag&ELEMENT_MARK_COARSEN);
}


inline void
Element::markForRefinement()
{
  tag |= ELEMENT_MARK_REFINE;
  if (tag&ELEMENT_MARK_COARSEN)
    tag ^= ELEMENT_MARK_COARSEN;
}

inline void
Element::markForCoarsening()
{
  tag |= ELEMENT_MARK_COARSEN;
  if (tag&ELEMENT_MARK_REFINE)
    tag ^= ELEMENT_MARK_REFINE;
}


inline void
Element::markParentForRefinement()
{
  if (parent) parent->markForRefinement();
}

inline void
Element::markChildrenForCoarsening()
{
  for (int nn = 0; nn < chcount; nn++)
  {
    child[nn]->markForCoarsening();
  }
}

inline void
Element::unmarkForRefinement()
{
  if (tag&ELEMENT_MARK_REFINE)
    tag ^= ELEMENT_MARK_REFINE;
}

inline void
Element::unmarkForCoarsening()
{
  if (tag&ELEMENT_MARK_COARSEN)
    tag ^= ELEMENT_MARK_COARSEN;
}


inline void
Element::unsetAllRefinementMarks()
{
  // unmark for refinement
  if (tag&ELEMENT_MARK_REFINE)
    tag ^= ELEMENT_MARK_REFINE;
  // unmark for coarsening
  if (tag&ELEMENT_MARK_COARSEN)
    tag ^= ELEMENT_MARK_COARSEN;
}

inline void
Element::unsetChildrenRefinementMarks()
{
  for (int i = 0; i < static_cast<int>(chcount); i++)
    child[i]->unsetAllRefinementMarks();
}

inline void
Element::ref()
{
  tag |= ELEMENT_MARK_REFERENCE;
}

inline void
Element::unref()
{
  if (tag&ELEMENT_MARK_REFERENCE)
    tag ^= ELEMENT_MARK_REFERENCE;
}

inline bool
Element::queryRef()
{
  return (tag&ELEMENT_MARK_REFERENCE);
}


inline void
Element::mark()
{
  tag |= ELEMENT_MARK_GENERAL;
}

inline void
Element::unmark()
{
  if (tag&ELEMENT_MARK_GENERAL)
    tag ^= ELEMENT_MARK_GENERAL;
}

inline bool
Element::queryMark()
{
  return (tag&ELEMENT_MARK_GENERAL);
}

inline void
Element::resetAllMarks()
{
  tag = 0;
}


inline void
Element::close()
{
  tag ^=ELEMENT_TYPE_IRREGULAR;
  if (tag&ELEMENT_TYPE_REGULAR)
    tag ^= ELEMENT_TYPE_REGULAR;
  if (tag&ELEMENT_TYPE_REFINED)
    tag ^= ELEMENT_TYPE_REFINED;
  if (tag&ELEMENT_TYPE_CLOSURE)
    tag ^=ELEMENT_TYPE_CLOSURE;
}

//=============================================================================
//  ELEMENT REFINEMENT
//=============================================================================


// Description
//
inline Element*
Element::getParent() const
{
  return parent;
}

// Description
//
inline void
Element::setParent(Element* theParent)
{
  parent = theParent;
}

// Description
//
inline int
Element::getNoOfChildren() const
{
  return (static_cast<int>(chcount));
}
// Description
//
inline Element*
Element::getChild(int ch) const
{
  return child[ch];
}


// Description
// for removal purposes
inline void
Element::resetPointerToChildren()
{
  if (isRefined() || isIrregular())
  {
    if (child)
    {
      // for (int nn = 0; nn < static_cast<int>(chcount); nn++)
      // {
      //   child[nn]->setParent(0);
      // }
      // delete the pointer
      delete [] child;
    }
  }

  // set no. of children
  chcount = 0;
  // set pointer to zero
  child = 0;
}



// Description
//
inline bool
Element::queryChildrenRefinement() const
{
  if(isRefined())
  {
    for (int i = 0; i < static_cast<int>(chcount); i++)
    {
      if (child[i]->isRefined())
        return true;
    }
  }
  return false;
}


// Description
//
inline bool
Element::queryChildrenMarks() const
{
  int i;
  int counter = 0;
  int noOfChildren = static_cast<int>(chcount);
  int minNoOfMarkedChildren = 1; //noOfChildren*refinementPolicy;

  if (isRegular())
    return false;

  for (i = 0; i < noOfChildren; i++)
  {
    if (child[i]->isMarkedForCoarsening())
    {
      bool edgFlag = true;
      for (int nn = 0; nn < child[i]->getNoOfEdges(); nn++)
      {
        Edge* edg = child[i]->getEdge(nn);
        if (edg->isRefined())
        {
          edgFlag = false;
          break;
        }
      }
      if (edgFlag)
        counter++;
    }
  }
  // set back children makrs
  if (counter < minNoOfMarkedChildren)
  {
    for (int i = 0; i < noOfChildren; i++)
      child[i]->unmarkForCoarsening();
    return false;
  }
  else
    return true;
}


//=============================================================================
//  GENERAL ELEMENT VALUE
//=============================================================================

// Description
//
inline void
Element::setValue(int val)
{
  value = (char)val;
}


// Description
//
inline int
Element::getValue() const
{
  return value;
}




} // namespace grd

#endif // ELEMENT_H
