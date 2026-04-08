#ifndef FILE_CFS_EDGE_HH
#define FILE_CFS_EDGE_HH

#include <bitset>

#include "General/Environment.hh"
#include "Utils/StdVector.hh"

namespace CoupledField {

  // Forward class declaration
  struct Elem;


  //! Defines a 1-dimensional grid entity, defined by two vertices
  
  //! This class models a 1-dimensional edge, defined by two vertices.
  //! Each edge has a unique orientation, where the first vertic always
  //! the node with the smaller number.
  //! Elements contain edges and often store their number signed integer,
  //! where the sign implies the orientation of their edge. 
  struct Edge 
  {
  //! Vector containing the neighbor elements
    
    //! This vector stores all the neighboring elements. Thus we can 
    //! easily find all connected elements for an edge.
    StdVector<Elem*> neighbors;
  };

  //! Defines a 2-dimensional grid entity, defined by 3 or 4 vertices
  struct Face {
    //! Normalize face orientation to match global one and return flagset
    
    //! This method re-orientates a given face to match the global orientation
    //! of the face. It returns the three-entry bitset which contains the
    //! orientation flags as defined in 
    //! Solin,Segeth: Higher-Order Finite Elements, p.165
    //! See also the faceFlags bitset of the class Elem.
    //! nodes is input and output variable - not it is not always sorted (BucklingPlate3D test elem 146)
    //! flags and nodes are output variables and are modified in-place.
    void Normalize( std::bitset<5>& flags, StdVector<UInt>& nodes) const; 

    
    //! Return permutation vector for face nodes according to bit-field
    static void GetSortedIndices( StdVector<UInt>& sorted, 
                                  const StdVector<UInt>& unsorted,
                                  UInt numVertices,
                                  const std::bitset<5>& flags );
    
    //! List with neighboring elements
    StdVector<Elem*> neighbors;
    
  private:

    //! static map for bits of quadrilateral faces
    static int quadBits[4][4][4];
    
    //! static map for bits of triangular face
    static int triaBits[3][3];
    
  };
}



#endif
