/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOTETRAHEDRONELEMENT_HH
#define GOTETRAHEDRONELEMENT_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GbPlane.hh"
#include "GbVec3.hh"
#include "GoDefaultElement.hh"
#include "GoVertex.hh"
#include "GoEdge.hh"
#include "GoTetrahedronElementBase.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

typedef enum {
  NO_REFINEMENT =       0x00000001,
  REGULAR_REFINEMENT =  0x00000002,
  COARSING =            0x00000004,
  
  /* needed to store irregular refinement type as well as for test if marked according to refinement */
  
  IRR_REF_TYP_1_013 =   0x00000010,
  IRR_REF_TYP_1_025 =   0x00000020,
  IRR_REF_TYP_1_345 =   0x00000040,
  IRR_REF_TYP_1_124 =   0x00000080,
  
  IRR_REF_TYP_2_0 =     0x00000100,
  IRR_REF_TYP_2_1 =     0x00000200,
  IRR_REF_TYP_2_2 =     0x00000400,
  IRR_REF_TYP_2_3 =     0x00000800,
  IRR_REF_TYP_2_4 =     0x00001000,
  IRR_REF_TYP_2_5 =     0x00002000,
  
  IRR_REF_TYP_3_01 =    0x00010000,
  IRR_REF_TYP_3_02 =    0x00020000,
  IRR_REF_TYP_3_03 =    0x00040000,
  IRR_REF_TYP_3_05 =    0x00080000,
  IRR_REF_TYP_3_10 =    0x00100000,
  IRR_REF_TYP_3_12 =    0x00200000,
  IRR_REF_TYP_3_13 =    0x00400000,
  IRR_REF_TYP_3_14 =    0x00800000,
  IRR_REF_TYP_3_20 =    0x01000000,
  IRR_REF_TYP_3_21 =    0x02000000,
  IRR_REF_TYP_3_24 =    0x04000000,
  IRR_REF_TYP_3_25 =    0x08000000,
  IRR_REF_TYP_3_30 =    0x10000000,
  IRR_REF_TYP_3_31 =    0x20000000,
  IRR_REF_TYP_3_34 =    0x40000000,
  IRR_REF_TYP_3_35 =    0x80000000,

  /* to be able to test for irregular refinement (without considering type) */
  F1_IRREGULAR_REFINEMENT = 0xffff3ff0,
  /* to be able to test for type */
  F1_IRR_REF_TYPE_1  =     0x000000f0,
  F1_IRR_REF_TYPE_2  =     0x00003f00,
  F1_IRR_REF_TYPE_3  =     0xffff0000,
  F1_IRR_REF_TYPE_4  =     0x00000000,
  /* to be able to test for edges */
  F1_REFINE_EDGE_0   =     0x111f0132,
  F1_REFINE_EDGE_1   =     0x02f10252,
  F1_REFINE_EDGE_2   =     0x2f220492,
  F1_REFINE_EDGE_3   =     0xf40408a2,
  F1_REFINE_EDGE_4   =     0x40481062,
  F1_REFINE_EDGE_5   =     0x888020c2,

  /* to be able to use the formulation as in the paper for checking actual refinement of tetrahedra */
  LEAF_ELEMENT        = NO_REFINEMENT,
  REFINED_REGULARLY   = REGULAR_REFINEMENT

} GbTetraRefinementTypesFlag1;

typedef enum {
  IRR_REF_TYP_3_41 =    0x001,
  IRR_REF_TYP_3_42 =    0x002,
  IRR_REF_TYP_3_43 =    0x004,
  IRR_REF_TYP_3_45 =    0x008,
  IRR_REF_TYP_3_50 =    0x010,
  IRR_REF_TYP_3_52 =    0x020,
  IRR_REF_TYP_3_53 =    0x040,
  IRR_REF_TYP_3_54 =    0x080,

  IRR_REF_TYP_4_04 =    0x100,
  IRR_REF_TYP_4_15 =    0x200,
  IRR_REF_TYP_4_32 =    0x400,

  /* to be able to test for irregular refinement (without considering type) */
  F2_IRREGULAR_REFINEMENT = 0x7ff,
  /* to be able to test for type */
  F2_IRR_REF_TYPE_1  =     0x000,
  F2_IRR_REF_TYPE_2  =     0x000,
  F2_IRR_REF_TYPE_3  =     0x0ff,
  F2_IRR_REF_TYPE_4  =     0x700,
  /* to be able to test for edges */
  F2_REFINE_EDGE_0   =     0x101,
  F2_REFINE_EDGE_1   =     0x212,
  F2_REFINE_EDGE_2   =     0x420,
  F2_REFINE_EDGE_3   =     0x244,
  F2_REFINE_EDGE_4   =     0x48f,
  F2_REFINE_EDGE_5   =     0x1f8

} GbTetraRefinementTypesFlag2;

typedef enum {
  IRREGULAR_REFINEMENT,
  /* to be able to test for type */
  IRR_REF_TYPE_1,
  IRR_REF_TYPE_2,
  IRR_REF_TYPE_3,
  IRR_REF_TYPE_4,
  /* to be able to test for edges */
  REFINE_EDGE_0,
  REFINE_EDGE_1,
  REFINE_EDGE_2,
  REFINE_EDGE_3,
  REFINE_EDGE_4,
  REFINE_EDGE_5,

  /* to be able to use the formulation as in the paper for checking actual refinement of tetrahedra */
  REFINED_IRREGULARLY = IRREGULAR_REFINEMENT

} GbTetraRefinementTypesCombination;


/*!
  class GoTetrahedronElement
 */
class GoTetrahedronElement 
  : public GoDefaultElement
{
public: 

  GoTetrahedronElement(int i=-1, int p=0, int l=0, const GbVec3<float>* n=NULL)
    : GoDefaultElement(i,p,l,n)
    { 
      tet_ = new GoTetrahedronElementBase<float>; 
      init();
    }
  GoTetrahedronElement(GoDefaultElementBase<float> *t)
    : GoDefaultElement(t)
    {
      tet_ = new GoTetrahedronElementBase<float>; 
      init();
    }

  virtual ~GoTetrahedronElement() { delete tet_; }


  virtual void init() {
    markedFor1_ = NO_REFINEMENT;
    markedFor2_ = (GbTetraRefinementTypesFlag2)0;
    refinement1_ = LEAF_ELEMENT;
    refinement2_ = (GbTetraRefinementTypesFlag2)0;
  }


  // Operations on vertices of the geometry object
  virtual void setVertex(int i, GoVertex<float> *v) {
    ::setVertex<GoTetrahedronElementBase<float>,float> (tet_,i,v);
  }
  virtual GoVertex<float> *getVertex(int i) const {
    return ::getVertex<GoTetrahedronElementBase<float>,float> (tet_,i);
  }
  virtual int getNumVertices() const { return 4; }

  // The "normal"
  virtual void computeNormal() {
  }

  // The layout of the object
  virtual int getEdgeList(int **l) const {
    static int table_tetra[12] = {0,1,0,2,0,3,1,2,2,3,3,1};
    *l = table_tetra;
    return 12;
  }
  virtual void setEdge(int i, GoEdge<float> *e) {
    ::setEdge<GoTetrahedronElementBase<float>,float> (tet_, i, e);
  }
  virtual GoEdge<float> *getEdge(int i) const {
    return ::getEdge<GoTetrahedronElementBase<float>,float> (tet_,i);
  }
  virtual int getNumEdges() const { return 6; }

  // The neighboring face objects
  virtual void setNeighbour(int i, GoGeometryElement<float> *face) {
    ::setNeighbour<GoTetrahedronElementBase<float>,float> (tet_,i,face);
  }
  virtual GoGeometryElement<float> *getNeighbour(int i) const {
    return ::getNeighbour<GoTetrahedronElementBase<float>,float> (tet_,i);
  }
  virtual int getNumNeighbours() const { return 4; }
  
  virtual void getPlanes (GbPlane<float> akPlane[4]) const 
    {
      GbVec3<float> kEdge01 = getEdgeVec(0);
      GbVec3<float> kEdge02 = getEdgeVec(1);
      GbVec3<float> kEdge03 = getEdgeVec(2);
      GbVec3<float> kEdge12 = getEdgeVec(3);
      GbVec3<float> kEdge13 = -getEdgeVec(5);

      GbVec3<float> c = kEdge13.cross(kEdge12);
      c.normalize();
      akPlane[0].setNormal(c);  // <v1,v3,v2>
      c = kEdge02.cross(kEdge03);
      c.normalize();
      akPlane[1].setNormal(c);  // <v0,v3,v2>
      c = kEdge03.cross(kEdge01);
      c.normalize();
      akPlane[2].setNormal(c);  // <v0,v1,v3>
      c = kEdge01.cross(kEdge02);
      c.normalize();
      akPlane[3].setNormal(c);  // <v0,v2,v1>

      float fDet = akPlane[3].getNormal().dot(kEdge01);
      int i;
      if ( fDet < -std::numeric_limits<float>::epsilon() ) {
        // normals are inner pointing, reverse their directions
	warningmsg("inner pointing normals ? eps="<<-std::numeric_limits<float>::epsilon()<<" det="<<fDet<<" n="<<akPlane[3].getNormal());
        for (i = 0; i < 4; ++i) {
	  c = -akPlane[i].getNormal();
	  akPlane[i].setNormal(c);
	}
      }

      for (i = 0; i < 4; ++i)
        akPlane[i].setConstant( getVertex(i)->getPosition().dot(akPlane[i].getNormal()) );
    }

  // Parents and children objects
  virtual GoGeometryElement<float> *getChild(int i) const {
    return ::getChild<GoTetrahedronElementBase<float>,float> (tet_,i);
  }
  virtual void setChild(int i, GoGeometryElement<float> *face) {
    ::setChild<GoTetrahedronElementBase<float>,float> (tet_,i,face);
  }
  virtual int getNumChildren() const { return 8; }

  // Get memory statistics and reduce memory space needed
  virtual void statistic() const {
    ::statistic<GoDefaultElementBase<float>,float> (this_);
    ::statistic<GoTetrahedronElementBase<float>,float> (tet_);
  }
  virtual void shrink() {
    ::shrink<GoDefaultElementBase<float>,float> (this_);
    ::shrink<GoTetrahedronElementBase<float>,float> (tet_);
  }

//  friend std::ostream& operator<<(std::ostream&, const GoTetrahedronElement&);

  virtual void setMarkedFor(GbTetraRefinementTypesFlag1 m) {
    markedFor1_ = m;
  }

  virtual void setMarkedFor(GbTetraRefinementTypesFlag2 m) {
    markedFor2_ = m;
  }

  virtual void setMarkedFor(GbTetraRefinementTypesCombination m) {
    switch (m) {
    case IRREGULAR_REFINEMENT:
      setMarkedFor(F1_IRREGULAR_REFINEMENT);
      setMarkedFor(F2_IRREGULAR_REFINEMENT);
      break;
    case IRR_REF_TYPE_1:
      setMarkedFor(F1_IRR_REF_TYPE_1);
      setMarkedFor(F2_IRR_REF_TYPE_1);
      break;
    case IRR_REF_TYPE_2:
      setMarkedFor(F1_IRR_REF_TYPE_2);
      setMarkedFor(F2_IRR_REF_TYPE_2);
      break;
    case IRR_REF_TYPE_3:
      setMarkedFor(F1_IRR_REF_TYPE_3);
      setMarkedFor(F2_IRR_REF_TYPE_3);
      break;
    case IRR_REF_TYPE_4:
      setMarkedFor(F1_IRR_REF_TYPE_4);
      setMarkedFor(F2_IRR_REF_TYPE_4);
      break;
    case REFINE_EDGE_0:
      setMarkedFor(F1_REFINE_EDGE_0);
      setMarkedFor(F2_REFINE_EDGE_0);
      break;
    case REFINE_EDGE_1:
      setMarkedFor(F1_REFINE_EDGE_1);
      setMarkedFor(F2_REFINE_EDGE_1);
      break;
    case REFINE_EDGE_2:
      setMarkedFor(F1_REFINE_EDGE_2);
      setMarkedFor(F2_REFINE_EDGE_2);
      break;
    case REFINE_EDGE_3:
      setMarkedFor(F1_REFINE_EDGE_3);
      setMarkedFor(F2_REFINE_EDGE_3);
      break;
    case REFINE_EDGE_4:
      setMarkedFor(F1_REFINE_EDGE_4);
      setMarkedFor(F2_REFINE_EDGE_4);
      break;
    case REFINE_EDGE_5:
      setMarkedFor(F1_REFINE_EDGE_5);
      setMarkedFor(F2_REFINE_EDGE_5);
      break;
    default:
      break;
    }
  }

  virtual GbBool isMarkedFor(GbTetraRefinementTypesFlag1 m) const {
    return ((markedFor1_ & m) != 0);
  }

  virtual GbBool isMarkedFor(GbTetraRefinementTypesFlag2 m) const {
    return ((markedFor2_ & m) != 0);
  }

  virtual GbBool isMarkedFor(GbTetraRefinementTypesCombination m) const {
    switch (m) {
    case IRREGULAR_REFINEMENT:
      return isMarkedFor(F1_IRREGULAR_REFINEMENT) || isMarkedFor(F2_IRREGULAR_REFINEMENT);
    case IRR_REF_TYPE_1:
      return isMarkedFor(F1_IRR_REF_TYPE_1) || isMarkedFor(F2_IRR_REF_TYPE_1);
    case IRR_REF_TYPE_2:
      return isMarkedFor(F1_IRR_REF_TYPE_2) || isMarkedFor(F2_IRR_REF_TYPE_2);
    case IRR_REF_TYPE_3:
      return isMarkedFor(F1_IRR_REF_TYPE_3) || isMarkedFor(F2_IRR_REF_TYPE_3);
    case IRR_REF_TYPE_4:
      return isMarkedFor(F1_IRR_REF_TYPE_4) || isMarkedFor(F2_IRR_REF_TYPE_4);
    case REFINE_EDGE_0:
      return isMarkedFor(F1_REFINE_EDGE_0) || isMarkedFor(F2_REFINE_EDGE_0);
    case REFINE_EDGE_1:
      return isMarkedFor(F1_REFINE_EDGE_1) || isMarkedFor(F2_REFINE_EDGE_1);
    case REFINE_EDGE_2:
      return isMarkedFor(F1_REFINE_EDGE_2) || isMarkedFor(F2_REFINE_EDGE_2);
    case REFINE_EDGE_3:
      return isMarkedFor(F1_REFINE_EDGE_3) || isMarkedFor(F2_REFINE_EDGE_3);
    case REFINE_EDGE_4:
      return isMarkedFor(F1_REFINE_EDGE_4) || isMarkedFor(F2_REFINE_EDGE_4);
    case REFINE_EDGE_5:
      return isMarkedFor(F1_REFINE_EDGE_5) || isMarkedFor(F2_REFINE_EDGE_5);
    default:
      break;
    }
    return false;
  }

  virtual void setRefinement(GbTetraRefinementTypesFlag1 r) {
    refinement1_ = r;
  }

  virtual void setRefinement(GbTetraRefinementTypesFlag2 r) {
    refinement2_ = r;
  }

  virtual void setRefinement(GbTetraRefinementTypesCombination r) {
    switch (r) {
    case IRREGULAR_REFINEMENT:
      setRefinement(F1_IRREGULAR_REFINEMENT);
      setRefinement(F2_IRREGULAR_REFINEMENT);
      break;
    case IRR_REF_TYPE_1:
      setRefinement(F1_IRR_REF_TYPE_1);
      setRefinement(F2_IRR_REF_TYPE_1);
      break;
    case IRR_REF_TYPE_2:
      setRefinement(F1_IRR_REF_TYPE_2);
      setRefinement(F2_IRR_REF_TYPE_2);
      break;
    case IRR_REF_TYPE_3:
      setRefinement(F1_IRR_REF_TYPE_3);
      setRefinement(F2_IRR_REF_TYPE_3);
      break;
    case IRR_REF_TYPE_4:
      setRefinement(F1_IRR_REF_TYPE_4);
      setRefinement(F2_IRR_REF_TYPE_4);
      break;
    case REFINE_EDGE_0:
      setRefinement(F1_REFINE_EDGE_0);
      setRefinement(F2_REFINE_EDGE_0);
      break;
    case REFINE_EDGE_1:
      setRefinement(F1_REFINE_EDGE_1);
      setRefinement(F2_REFINE_EDGE_1);
      break;
    case REFINE_EDGE_2:
      setRefinement(F1_REFINE_EDGE_2);
      setRefinement(F2_REFINE_EDGE_2);
      break;
    case REFINE_EDGE_3:
      setRefinement(F1_REFINE_EDGE_3);
      setRefinement(F2_REFINE_EDGE_3);
      break;
    case REFINE_EDGE_4:
      setRefinement(F1_REFINE_EDGE_4);
      setRefinement(F2_REFINE_EDGE_4);
      break;
    case REFINE_EDGE_5:
      setRefinement(F1_REFINE_EDGE_5);
      setRefinement(F2_REFINE_EDGE_5);
      break;
    default:
      break;
    }
  }

  virtual GbBool testRefinement(GbTetraRefinementTypesFlag1 r) const {
    return ((refinement1_ & r) != 0);
  }

  virtual GbBool testRefinement(GbTetraRefinementTypesFlag2 r) const {
    return ((refinement2_ & r) != 0);
  }

  virtual GbBool testRefinement(GbTetraRefinementTypesCombination r) const {
    switch (r) {
    case IRREGULAR_REFINEMENT:
      return testRefinement(F1_IRREGULAR_REFINEMENT) || testRefinement(F2_IRREGULAR_REFINEMENT);
    case IRR_REF_TYPE_1:
      return testRefinement(F1_IRR_REF_TYPE_1) || testRefinement(F2_IRR_REF_TYPE_1);
    case IRR_REF_TYPE_2:
      return testRefinement(F1_IRR_REF_TYPE_2) || testRefinement(F2_IRR_REF_TYPE_2);
    case IRR_REF_TYPE_3:
      return testRefinement(F1_IRR_REF_TYPE_3) || testRefinement(F2_IRR_REF_TYPE_3);
    case IRR_REF_TYPE_4:
      return testRefinement(F1_IRR_REF_TYPE_4) || testRefinement(F2_IRR_REF_TYPE_4);
    case REFINE_EDGE_0:
      return testRefinement(F1_REFINE_EDGE_0) || testRefinement(F2_REFINE_EDGE_0);
    case REFINE_EDGE_1:
      return testRefinement(F1_REFINE_EDGE_1) || testRefinement(F2_REFINE_EDGE_1);
    case REFINE_EDGE_2:
      return testRefinement(F1_REFINE_EDGE_2) || testRefinement(F2_REFINE_EDGE_2);
    case REFINE_EDGE_3:
      return testRefinement(F1_REFINE_EDGE_3) || testRefinement(F2_REFINE_EDGE_3);
    case REFINE_EDGE_4:
      return testRefinement(F1_REFINE_EDGE_4) || testRefinement(F2_REFINE_EDGE_4);
    case REFINE_EDGE_5:
      return testRefinement(F1_REFINE_EDGE_5) || testRefinement(F2_REFINE_EDGE_5);
    default:
      break;
    }
    return false;
  }

  virtual GbBool isMarkedAccordingToRefinement() const {
    return ( (markedFor1_ == refinement1_) && (markedFor2_ == refinement2_) );
  }

protected:
  GoTetrahedronElementBase<float> *tet_;

  GbTetraRefinementTypesFlag1 markedFor1_, refinement1_;
  GbTetraRefinementTypesFlag2 markedFor2_, refinement2_;
};

#endif // GOTETRAHEDRONELEMENT_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.12  2002/03/18 09:58:56  prkipfer
| refactored element structure
|
| Revision 1.11  2001/12/18 13:23:51  prkipfer
| changed edge enumeration to be consistent throughout the lib
|
| Revision 1.10  2001/09/12 09:28:43  prkipfer
| introduced adaptive tet subdivision
|
| Revision 1.9  2001/01/02 15:21:35  prkipfer
| introduced cloning and support for new base classes
|
| Revision 1.8  2000/11/07 17:27:27  prkipfer
| introduced external polymorphism for vertices
|
| Revision 1.7  2000/11/07 09:58:53  uflabsik
| docu added
|
| Revision 1.6  2000/10/17 14:33:03  uflabsik
| added new primitives
|
| Revision 1.5  2000/09/07 16:57:17  prkipfer
| added subdivision
|
| Revision 1.4  2000/09/03 14:06:14  uflabsik
| functionality added
|
| Revision 1.3  2000/08/28 09:25:21  prkipfer
| added localized subdivision mechanism
|
| Revision 1.2  2000/07/21 10:40:43  prkipfer
| dropped g++ support
|
| Revision 1.1  2000/07/03 16:15:01  prkipfer
| extended IO subsystem, created class registry, added tetrahedron
|
|
+---------------------------------------------------------------------*/


