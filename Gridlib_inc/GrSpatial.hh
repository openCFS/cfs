/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRSPATIAL_HH
#define GRSPATIAL_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GbVec3.hh"
#include "GbMatrix3.hh"
#include "GrRenderState.hh"
#include "GrBound.hh"

class GrNode;

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GrSpatial
{
public:
  // local transform access
  INLINE GbMatrix3<float>& rotate ();
  INLINE const GbMatrix3<float>& rotate () const;
  INLINE GbVec3<float>& translate ();
  INLINE const GbVec3<float>& translate () const;
  INLINE float& scale ();
  INLINE const float& scale () const;

  // world transform access ('set' should be used only by controllers)
  INLINE GbMatrix3<float>& worldRotate ();
  INLINE const GbMatrix3<float>& worldRotate () const;
  INLINE GbVec3<float>& worldTranslate ();
  INLINE const GbVec3<float>& worldTranslate () const;
  INLINE float& worldScale ();
  INLINE const float& worldScale () const;
  INLINE void setWorldTransformToIdentity ();

  // world bound access
  INLINE GrBound<float>& worldBound ();
  INLINE const GrBound<float>& worldBound () const;

  // culling
  INLINE GbBool& forceCull ();
  INLINE const GbBool& forceCull () const;

  // render state
  GrRenderState* setRenderState (GrRenderState* pkState);
  GrRenderState* getRenderState (GrRenderState::Type eType);
  GrRenderState* removeRenderState (GrRenderState::Type eType);
  void removeAllStates ();

  // updates (GS = geometric state, RS = renderer state)
  void updateGS (float fAppTime, GbBool bInitiator = true);
  void updateRS (GbBool bInitiator = true);

  // parent access
  INLINE GrNode* getParent ();

  // support for searching by name
  INLINE void setName(GbString& s);
//  virtual Object* getObjectByName (const char* acName);

protected:
  // allow internal calls to updates and drawing
  friend class GrBSPNode;
  friend class GrConvexRegionManager;
  friend class GrNode;
  friend class GrRenderer;
//  friend class SwitchNode;

  // construction (abstract base class)
  GrSpatial ();
  virtual ~GrSpatial ();

  // parent access
  INLINE void setParent (GrNode* pkParent);

  // geometric updates
  virtual void updateWorldData (float fAppTime);
  virtual void updateWorldBound () {}; //!! = 0;
  void propagateBoundToRoot ();

  // render state updates
  virtual void updateRenderState () {}; //!! = 0;
  void propagateStateFromRoot ();
  void restoreStateToRoot ();

  // drawing
  void onDraw (GrRenderer& rkRenderer);
  virtual void draw (GrRenderer& rkRenderer) = 0;

  // parents
  GrNode* parent_;
  GbString name_;

  // local transforms
  GbMatrix3<float> rotate_;
  GbVec3<float> translate_;
  float scale_;

  // world transforms
  GbMatrix3<float> worldRotate_;
  GbVec3<float> worldTranslate_;
  float worldScale_;

  // world bound
  GrBound<float> worldBound_;
  
  // culling
  GbBool forceCull_;

  // render state
  GrRenderState::List* stateList_;
};

#ifndef OUTLINE
#include "GrSpatial.in"
#endif

#endif // GRSPATIAL_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:58  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.3  2001/06/26 12:25:26  prkipfer
| minor fixes for clean Linux compile
|
| Revision 1.2  2001/06/19 16:30:22  prkipfer
| fixed typos and Linux compiler target
|
| Revision 1.1  2001/06/18 11:23:16  prkipfer
| introduced new classes for spatial partitioning and hierarchical culling
|
|
+---------------------------------------------------------------------*/
