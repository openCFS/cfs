/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GBMESHFUNCTIONS_HH
#define GBMESHFUNCTIONS_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

struct Delete {
  void operator()(GoGeometryElement<float> *t) {
    if (t->getParent())
      t->getParent()->clearChildren();
    delete t;
  }
  void operator()(GoEdge<float> *e) {
    delete e;
  }
  void operator()(GoVertex<float> *t) {
    delete t;
  }
};

struct DeleteMarked {
  GbBool operator()(GoGeometryElement<float> *t) {
    if(t->testFlag(DELETED_G)) {
      if (t->getParent())
	t->getParent()->clearChildren();
      delete t;
      return true;
    }
    return false;
  }

  GbBool operator()(GoEdge<float> *e) {
    if(e->testFlag(DELETED_E)) {
      delete e;
      return true;
    }
    return false;
  }

  GbBool operator()(GoVertex<float> *t) {
    if(t->testFlag(DELETED_V)) {
      delete t;
      return true;
    }
    return false;
  }
};

struct ComputeFaceNormal {
  void operator()(GoGeometryElement<float> *t) {
    t->computeNormal();
  }
};

struct DeleteEmptyVectors {
  GbBool operator()(std::vector<GoGeometryElement<float>*> *vec) {
    if (vec->empty()) {
      delete vec;
      return true;
    }
   
    return false;
  }

  GbBool operator()(std::vector<GoEdge<float>*> *vec) {
    if (vec->empty()) {
      delete vec;
      return true;
    }
   
    return false;
  }

  GbBool operator()(std::vector<GoVertex<float>*> *vec) {
    if (vec->empty()) {
      delete vec;
      return true;
    }
   
    return false;
  }
};

struct ComputeVertexNormal {
  void operator()(GoVertex<float> *t) {
    std::vector<GoGeometryElement<float>*> S;

    getNeighbourFaces(S,t);
    t->computeNormal(S);
    
    S.clear();

  }
  // get a list of neighboring face objects
  void getNeighbourFaces(std::vector<GoGeometryElement<float>*>& S, 
			 GoVertex<float> *v) const
    {
      GoGeometryElement<float> *f2, *f1, *f;
      GoVertex<float> *p0, *p1, *p2;
      int i, j=0;

//!! this should be disabled, because normals can be calculated even
//!! if there is no closed fan of faces around the vertex
//        if (v->testFlag(BOUNDARY_V))
//  	return;
  
      f = v->getElement();

      if (f != NULL) {
	i = f->findVertex(v);
	p0 = f->getVertex( (i+2)%3);
	p1 = f->getVertex( (i+1)%3);
    
	while (j < v->getValence()) {
	  j++;
	  f1 = f->otherFace(p0);

	  if (f1==NULL) 
	    return;

	  f2 = f1;
      
	  
	  p2 = f1->otherVertex(f);
	  
	  S.push_back(f2);   
      
	  p0 = p1;
	  p1 = p2;
      
	  f = f1;
      
	}
      }
      else
      {
	warningmsg("no face yet associated with this vertex");
      }
      
    }
};



struct ComputeUmbrella {
  void operator()(GoVertex<float> *t) {
    std::vector<GoVertex<float>*> S;
    GoVertex<float> *pi;
    GbVec3<float> p(0,0,0);
    float alpha = 0.5;

    if (t->testFlag(BOUNDARY_V))
      return;


    getNeighbourVertices(S,t);

    // here the computation of the umbrella
    int valence = S.size();
    for (int i=0; i<valence; i++) {
      pi = S[i];
      
      p += pi->getPosition();

    }

    p = t->getPosition() * alpha + p /(float)valence * (1-alpha);
    t->setPosition(p);

    S.clear();

  }
  // get a list of neighboring vertices
  void getNeighbourVertices(std::vector<GoVertex<float>*>& S, 
			    GoVertex<float> *v) const
    {
      GoGeometryElement<float> *f1, *f;
      GoVertex<float> *p0, *p1, *p2, *pstart;
      int i, j=0;

      if (v->testFlag(BOUNDARY_V))
	return;
  
      f = v->getElement();

      if (f != NULL) {
	i = f->findVertex(v);
	p0 = f->getVertex( (i+2)%3);
	p1 = f->getVertex( (i+1)%3);
    
	pstart = p0;

	do {
	  j++;
	  f1 = f->otherFace(p0);

	  if (f1==NULL) 
	    return;

	  p2 = f1->otherVertex(f);

	  S.push_back(p2);   
      
	  p0 = p1;
	  p1 = p2;
      
	  f = f1;
      
	} while (p0 != pstart);
      }
      else
      {
	warningmsg("no face yet associated with this vertex");
      }
      
    }

};


struct ComputeCurvatureFlow {
  void operator()(GoVertex<float> *t) {
    std::vector<GoVertex<float>*> S;
    GbVec3<float> p(0,0,0), pi, pj, pj_m1, pj_p1;
    float area=0;

    if (t->testFlag(BOUNDARY_V))
      return;

    pi = t->getPosition();
    
    getNeighbourVertices(S,t);

    // here the computation of the umbrella
    int valence = S.size();
    for (int i=0; i<valence; i++) {
      pj_m1 = S[ (i-1+valence) % valence]->getPosition();
      pj = S[i]->getPosition();
      pj_p1 = S[ (i+1) % valence]->getPosition();
      
      // area += 0.5 * ((pj_p1-pi).cross(pj-pi)).getNorm();

      
      float cota = ((pj_p1-pi)|(pj_p1-pj)) / ((pj_p1-pi).cross(pj_p1-pj)).getNorm();
      float cotb = ((pj_m1-pj)|(pj_m1-pi)) / ((pj_m1-pj).cross(pj_m1-pi)).getNorm();

      // std::cerr << cota << " " << cotb << std::endl;


      area += cota+cotb;
      p += (pi-pj)*(cota+cotb);

      //std::cerr << cota+cotb << std::endl;
      
    }

    p = pi - p/area * (float)0.2;
    t->setPosition(p);

    S.clear();
  }
  
  // get a list of neighboring vertices
  void getNeighbourVertices(std::vector<GoVertex<float>*>& S, 
			    GoVertex<float> *v) const
    {
      GoGeometryElement<float> *f1, *f;
      GoVertex<float> *p0, *p1, *p2, *pstart;
      int i, j=0;

      if (v->testFlag(BOUNDARY_V))
	return;
  
      f = v->getElement();

      if (f != NULL) {
	i = f->findVertex(v);
	p0 = f->getVertex( (i+2)%3);
	p1 = f->getVertex( (i+1)%3);
    
	pstart = p0;

	do {
	  j++;
	  f1 = f->otherFace(p0);

	  if (f1==NULL) 
	    return;

	  p2 = f1->otherVertex(f);

	  S.push_back(p2);   
      
	  p0 = p1;
	  p1 = p2;
      
	  f = f1;
      
	} while (p0 != pstart);
      }
      else
      {
	warningmsg("no face yet associated with this vertex");
      }
      
    }

};

#endif // GBMESHFUNCTIONS_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:07  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.15  2002/03/18 09:57:33  prkipfer
| refactored element structure
|
| Revision 1.14  2001/09/12 11:53:01  prkipfer
| introduced adaptive tet subdivision
|
| Revision 1.13  2001/08/16 16:56:11  prkipfer
| improved type safety for template parameter
|
| Revision 1.12  2001/02/13 11:12:21  prkipfer
| fixed vertex normal bug
|
| Revision 1.11  2000/11/07 17:27:54  prkipfer
| introduced external polymorphism for vertices
|
| Revision 1.10  2000/10/31 16:29:06  uflabsik
| new version
|
| Revision 1.9  2000/09/07 16:55:49  prkipfer
| added subdivision
|
| Revision 1.8  2000/08/28 09:25:24  prkipfer
| added localized subdivision mechanism
|
| Revision 1.7  2000/08/22 17:15:23  prkipfer
| added possibility to query mesh contents
|
| Revision 1.6  2000/08/18 11:16:45  uflabsik
| Bug fixed
|
| Revision 1.5  2000/08/03 13:07:36  uflabsik
| curvature flow operator added
|
| Revision 1.4  2000/07/28 15:25:28  uflabsik
| functionality added
|
| Revision 1.3  2000/07/18 09:38:30  uflabsik
| some names changed
|
| Revision 1.2  2000/06/16 13:32:20  prkipfer
| started mesh class hierarchy
|
| Revision 1.1  2000/06/14 15:39:18  prkipfer
| improved base classes and added funcstruct processing for mesh
|
|
+---------------------------------------------------------------------*/
