/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GBEXTERNALADAPTERS_HH
#define GBEXTERNALADAPTERS_HH

#include "GbDefines.hh"

// external polymorphism behaviour template methods

// common adapters
template <class T, class D>
void 
setNormal(T *t, const GbVec3<D>& n)
{
  t->setNormal(n);
}

template <class T, class D>
GbVec3<D> 
getNormal(const T *t)
{
  return t->getNormal();
}

template <class T, class D>
void 
setId(T *t, int i)
{
  t->setId(i);
}

template <class T, class D>
int 
getId(const T *t)
{
  return t->getId();
}

template <class T, class D>
void
shrink(T *t)
{
  t->shrink();
}

template <class T, class D>
void
statistic(const T *t)
{
  t->poolStatistic();
}


// adapters for elements
template <class T, class D>
void 
setVertex(T *t, int i, GoVertex<D> *v)
{
  t->setVertex(i,v);
}

template <class T, class D>
GoVertex<D> *
getVertex(const T *t, int i)
{
  return t->getVertex(i);
}

template <class T, class D>
void 
setEdge(T *t, int i, GoEdge<D> *e)
{
  t->setEdge(i,e);
}

template <class T, class D>
GoEdge<D> *
getEdge(const T *t, int i)
{
  return t->getEdge(i);
}

template <class T, class D>
void 
setFlag(T *t, GbGeoStatusFlag f)
{
  t->setFlag(f);
}

template <class T, class D>
void
delFlag(T *t, GbGeoStatusFlag f)
{
  t->delFlag(f);
}

template <class T, class D>
GbBool 
testFlag(const T *t, GbGeoStatusFlag f)
{
  return t->testFlag(f);
}

template <class T, class D>
GbBool
subdivide(T *t)
{
  return t->subdivide();
}

template <class T, class D>
void 
setPartition(T *t, int i)
{
  t->setPartition(i);
}

template <class T, class D>
int 
getPartition(const T *t)
{
  return t->getPartition();
}

template <class T, class D>
void 
setNeighbour(T *t, int i, GoGeometryElement<D> *face)
{
  t->setNeighbour(i,face);
}

template <class T, class D>
GoGeometryElement<D> *
getNeighbour(const T *t, int i)
{
  return t->getNeighbour(i);
}

template <class T, class D>
GoGeometryElement<D> *
getChild(const T *t, int i)
{
  return t->getChild(i);
}

template <class T, class D>
void
setChild(T *t, int i, GoGeometryElement<D> *face)
{
  t->setChild(i, face);
}

template <class T, class D>
GoGeometryElement<D> *
getParent(const T *t)
{
  return t->getParent();
}

template <class T, class D>
void
setParent(T *t, GoGeometryElement<D> *face)
{
  t->setParent(face);
}


// adapters for vertices
template <class T, class D>
void 
setPosition(T *t, D x, D y, D z)
{
  t->setPosition(x,y,z);
}

template <class T, class D>
void 
setPosition(T *t, const GbVec3<D>& v)
{
  t->setPosition(v);
}

template <class T, class D>
void 
getPosition(const T *t, D& xx, D& yy, D& zz) 
{
  t->getPosition(xx,yy,zz);
}

template <class T, class D>
GbVec3<D> 
getPosition(const T *t) 
{
  return t->getPosition();
}

template <class T, class D>
void 
computeNormal(T *t, std::vector<GoGeometryElement<D> *>& S)
{
  t->computeNormal(S);
}

template <class T, class D>
void 
setNormal(T *t, D x, D y, D z)
{
  t->setNormal(x,y,z);
}

template <class T, class D>
void 
getNormal(const T *t, D& xx, D& yy, D& zz) 
{
  t->getNormal(xx,yy,zz);
}

template <class T, class D>
void 
setFlag(T *t, GbVertexFlag f)
{
  t->setFlag(f);
}

template <class T, class D>
void 
delFlag(T *t, GbVertexFlag f)
{
  t->delFlag(f);
}

template <class T, class D>
GbBool 
testFlag(const T *t, GbVertexFlag f) 
{
  return t->testFlag(f);
}

template <class T, class D>
void 
setFlag(T *t, GbEdgeStatusFlag f)
{
  t->setFlag(f);
}

template <class T, class D>
void
delFlag(T *t, GbEdgeStatusFlag f)
{
  t->delFlag(f);
}

template <class T, class D>
GbBool 
testFlag(const T *t, GbEdgeStatusFlag f)
{
  return t->testFlag(f);
}

template <class T, class D>
void 
setElement(T *t, GoGeometryElement<D> *f)
{
  t->setElement(f);
}

template <class T, class D>
GoGeometryElement<D> *
getElement(const T *t) 
{
  return t->getElement();
}

template <class T, class D>
void 
setValence(T *t, int v)
{
  t->setValence(v);
}

template <class T, class D>
int 
getValence(const T *t) 
{
  return t->getValence();
}

template <class T, class D>
void 
setValue(T *t, D v)
{
  t->setValue(v);
}

template <class T, class D>
D
getValue(const T *t) 
{
  return t->getValue();
}

template <class T, class D>
void 
setDensity(T *t, D v)
{
  t->setDensity(v);
}

template <class T, class D>
D
getDensity(const T *t) 
{
  return t->getDensity();
}

template <class T, class D>
void 
setEnergy(T *t, D v)
{
  t->setEnergy(v);
}

template <class T, class D>
D
getEnergy(const T *t) 
{
  return t->getEnergy();
}

template <class T, class D>
void 
setMomentum(T *t, const GbVec3<D>& n)
{
  t->setMomentum(n);
}

template <class T, class D>
GbVec3<D> 
getMomentum(const T *t)
{
  return t->getMomentum();
}

template <class T, class D>
void 
setVertex(T *t, GoVertex<D> *v)
{
  t->setVertex(v);
}

template <class T, class D>
GoVertex<D> *
getVertex(const T *t)
{
  return t->getVertex();
}

template <class T, class D>
GoVertex<D> *
getMidPoint(const T *t)
{
  return t->getMidPoint();
}

template <class T, class D>
void
setLevel(T *t, int l)
{
  t->setLevel(l);
}

template <class T, class D>
int
getLevel(const T *t)
{
  return t->getLevel();
}

template <class T, class D>
void
addElement(T *t)
{
  t->addElement();
}

template <class T, class D>
void
removeElement(T *t)
{
  t->removeElement();
}

template <class T, class D>
int
getNumElements(const T *t)
{
  return t->getNumElements();
}

template <class T, class D>
void
addRefinement(T *t)
{
  t->addRefinement();
}

template <class T, class D>
void
removeRefinement(T *t)
{
  t->removeRefinement();
}

template <class T, class D>
void
clearRefinement(T *t)
{
  t->clearRefinement();
}

template <class T, class D>
GbBool
isMarkedForRefinement(T *t)
{
  return t->isMarkedForRefinement();
}

template <class T, class D>
void
setVertices(T *t, GoVertex<D> *one, GoVertex<D> *two)
{
  t->setVertices(one, two);
}

template <class T, class D>
void
setVertices(T *t, GoVertex<D> *one, GoVertex<D> *two, GoVertex<D> *midPoint)
{
  t->setVertices(one, two, midPoint);
}

#endif // GBEXTERNALADAPTERS_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:56  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.5  2002/03/18 09:58:54  prkipfer
| refactored element structure
|
| Revision 1.4  2001/09/12 09:28:40  prkipfer
| introduced adaptive tet subdivision
|
| Revision 1.3  2001/02/13 11:03:57  prkipfer
| introduced boundary vertices
|
| Revision 1.2  2001/01/02 15:21:33  prkipfer
| introduced cloning and support for new base classes
|
| Revision 1.1  2000/11/07 17:27:25  prkipfer
| introduced external polymorphism for vertices
|
|
+---------------------------------------------------------------------*/
