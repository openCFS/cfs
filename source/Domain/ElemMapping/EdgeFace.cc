#include "EdgeFace.hh"
#include "MatVec/Vector.hh"

namespace CoupledField {

// static variable initialization

char Face::quadBits[4][4] = { { -1, 7, -1, 3 }, { 6, -1, 2, -1 },
                              { -1, 0, -1, 4 }, { 1, -1, 5, -1 } };

char Face::triaBits[3][3] = { { -1, 4, 0 }, { 1, -1, 6 }, { 5, 2, -1 } };

bool operator<( const Face& a, const Face& b ) {

  UInt i = 0;
  while( i < a.nodes.GetSize() - 1 && ( a.nodes[i] == b.nodes[i] ) ) {
    i++;
  };

  return a.nodes[i] < b.nodes[i];
}

bool operator<( const Edge& a, const Edge& b ) {
  if( a.nodes[0] < b.nodes[0] ) {
    return true;
  }
  if( a.nodes[0] > b.nodes[0] ) {
    return false;
  }

  if( a.nodes[1] < b.nodes[1] ) {
    return true;
  }
  else {
    return false;
  }
}

//  permutation matrix for quadrilateral faces
//  char quadPerm[2][2][2][4] =
//  {
//   {
//    {
//     { 2,1,0,3 }, // 0
//     { 3,0,1,2 }  // 1
//    },
//    {
//     { 1,2,3,0 }, // 2
//     { 0,3,2,1 }  // 3
//    },
//   },
//   {
//    {
//     { 2,3,0,1 }, // 4
//     { 3,2,1,0 }  // 5
//    },
//    {
//     { 1,0,3,2 }, // 6
//     { 0,1,2,3 }  // 7
//    }
//   }
//  };

char quadPerm[8][4] = 
{ 
 { 2, 1, 0, 3 }, // 0
 { 3, 0, 1, 2 }, // 1
 { 1, 2, 3, 0 }, // 2
 { 0, 3, 2, 1 }, // 3
 { 2, 3, 0, 1 }, // 4
 { 3, 2, 1, 0 }, // 5
 { 1, 0, 3, 2 }, // 6
 { 0, 1, 2, 3 }, // 7
};

//  permutation matrix for triangular faces
char triaPerm[8][3] = 
{ { 0, 2, 1 }, // 0
  { 1, 0, 2 }, // 1
  { 2, 1, 0 }, // 2
  { -1, -1, -1 }, // 3 -> not defined
  { 0, 1, 2 }, // 4
  { 2, 0, 1 }, // 5
  { 1, 2, 0 }, // 6
  { -1, -1, -1 } // 7 -> not defined
};

void Face::GetSortedIndices( StdVector<UInt>& sorted, 
                             const StdVector<UInt>& unsorted,
                             UInt numVertices,
                             const std::bitset<3>& flags ) {

  Integer num = flags.to_ulong();
  // indices contains the
  if( numVertices == 4 ) {
    sorted.Resize( 4 );
//    std::cerr << "permutation array is " <<
//        int(quadPerm[num][0]) << ", " << int(quadPerm[num][1]) << ", " 
//        << int(quadPerm[num][2]) << ", " << int(quadPerm[num][3]) << std::endl;
    sorted[0] = unsorted[quadPerm[num][0]]-1;
    sorted[1] = unsorted[quadPerm[num][1]]-1;
    sorted[2] = unsorted[quadPerm[num][2]]-1;
    sorted[3] = unsorted[quadPerm[num][3]]-1;
  }
  else {
    sorted.Resize( 3 );
    sorted[0] = unsorted[triaPerm[num][0]]-1;
    sorted[1] = unsorted[triaPerm[num][1]]-1;
    sorted[2] = unsorted[triaPerm[num][2]]-1;
  }
}

void Face::Normalize( std::bitset<3>& flags ) {

  /* 
   In the case of higher order elements, we have to guarantee a global
   orientation of faces, described by the direction I and II. 
   By convention, we assume that the global orientation of a face is defined 
   by an ascending ordering of the nodal connectivity in such a way:
   
   C +----+ C   ^ II    Global orientation:  index(A) < index(B) < index(C)
     |    |     |       
     |    |     +--> I         I-direction: A -> B
   A +----+ B                 II-direction: A -> C                    
   
   Based on: 
   Solin, Segeth: "Higher-Order Finite-Element Methods", 2004, p. 164-170
   
   
   For quadrilateral faces, we thus have 8 possible orientations of the face,
   i.e. 8 different possibilities how the local element directions (xi/eta)
   relate to the global directions (I/II).  
   
   To describe them, we introduce for every face a 3-bit array, where the 
   flags have the following meaning:
   
   2   1   0     indices( 2 = most significant bit)
   +-----------+
   | 1 | 1 | 1 |  faceFlags(can be also set in integer representation, 
   +-----------+            i.e. 111=7)
   |   |   |
   |   |   \-- true  (=1), if  xi = +dirI/dirII
   |   |       false (=0), if  xi = -dirI/dirII
   |   \------ true  (=1), if eta = +dirI/dirII
   |           false (=0), if eta = -dirI/dirII
   \---------- true  (=1), if local and global directions match in same  
   order (|xi| = |I|,  |eta| = |II|)
   false (=0), if local and global directions match in reverse 
   order (|xi| = |II|, |eta| = |I|)
   
   Here are two examples for the orientation (Note: numbers in the interior
   denote the orientation of the reference elements, numbers on the outside 
   denote the nodal connectivity):
   
   Case 1
   ======
   
   connect:   1,2,3,4
   faceFlags: [111] 
   
   4+---------+3             
   |4  ^eta 3|     ^ I      xi  =   I
   |   |     |     |        eta =  II
   |   +->xi |     +-> II   
   |         |              faceFlags[2]=1/true // order: |xi|=|I|, |eta|=|II|
   |1       2|              faceFlags[1]=1/true // eta = eta
   1+---------+2             faceFlags[0]=1/true //  xi = xi
   

   indices                       indices
   0 1 2 3   permutate to match  0 1 2 3    take first two indices 0 1 2
   +-------+  global orientation +-------+   as key to quadBits    +-----+
   |1|2|3|4|       =>            |1|2|3|4|     =>   0: row     =>  |1|1|1| = 7
   +-------+                     +-------+          1: column      +-----+
   connect                       connect                         faceFlags
   (reordered)
   
   Case 2
   ======
   
   connect:   4,1,2,3
   faceFlags: [0,1,0] 
   
   3+---------+2
   |4  ^eta 3|     ^ I      xi  = -II
   |   |     |  II |        eta =   I
   |   +->xi |   <-+   
   |         |              faceFlags[2]=0/false // order: |xi|=|II|, |eta|=|I|
   |1       2|              faceFlags[1]=1/true  // xi = xi 
   4+---------+1             faceFlags[0]=0/false // eta = -eta
   
   
   indices                       indices
   0 1 2 3   permutate to match  1 2 3 0    take first two indices 2 1 0
   +-------+  global orientation +-------+   as key to quadBits    +-----+
   |4|1|2|3|       =>            |1|2|3|4|     =>   1: row     =>  |0|1|0| = 2
   +-------+                     +-------+          2: column      +-----+
   connect                       connect                         faceFlags
   (reordered)
   */

//  std::cerr << "\n=====================================\n"
//      << " Face Orientation\n" << "=====================================\n";
//  std::cerr << "\tconnect:" << nodes.ToString() << std::endl;

  StdVector<UInt> indices( nodes.GetSize() );
  UInt size = nodes.GetSize();

  // copy unsorted node Vector
  StdVector<UInt> unsorted;
  unsorted = nodes;
  
  // initialize indices array
  for( UInt i = 0; i < size; i++ ) {
    indices[i] = i;
  }

  // Insertion sort algorithm
  // ------------------------
  UInt j, comp;

  for( UInt i = 1; i < size; i++ ) {
    comp = nodes[i];
    j = i;
    while( ( j > 0 ) && ( nodes[j - 1] > comp ) ) {
      nodes[j] = nodes[j - 1];
      indices[j] = indices[j - 1];
      j = j - 1;
    }
    nodes[j] = comp;
    indices[j] = i;
  }
  // -----------------------

  // fetch orientation flags
  if( size == 4 ) {
    flags = std::bitset<3>( quadBits[indices[0]][indices[1]] );

    // switch 3rd and 4th node, so that we obtain a conforming face 
    // afterwards
    UInt val = nodes[3];
    nodes[3] = nodes[2];
    nodes[2] = val;

  }
  else if( size == 3 ) {
    flags = std::bitset<3>( triaBits[indices[0]][indices[1]] );
  }
  
//  // Check flags for orientation
//  std::cerr << "\tsorted: " << nodes.ToString() << std::endl;
//  std::cerr << "\tindices: " << indices.ToString() << std::endl;
//  std::cerr << "\tflag: " << flags << " (" << flags.to_ulong() << ")"
//      << std::endl;
//  std::cerr << "\tflag: [0] = " << flags.test( 0 ) << std::endl
//      << "\t      [1] = " << flags.test( 1 ) << std::endl << "\t      [2] = "
//      << flags.test( 2 ) << std::endl;

//  // obtain re-ordering array
//  StdVector<UInt> ind( size );
//  GetPermutation( ind,unsorted, flags );
//  std::cerr << "\tpermutation: " << ind.ToString() << std::endl;
//  GetPermutation( ind,nodes, flags );
//  std::cerr << "\t2ndpermutation: " << ind.ToString() << std::endl;
//
//  // print debug information
//  std::string xiString, etaString;
//  if( flags.test( 2 ) == true ) {
//    xiString = flags.test( 0 ) ? "I" : "-I";
//    etaString = flags.test( 1 ) ? "II" : "-II";
//  }
//  else {
//    xiString = flags.test( 0 ) ? "II" : "-II";
//    etaString = flags.test( 1 ) ? "I" : "-I";
//  }
//  std::cerr << "\txi =  " << xiString << std::endl;
//  std::cerr << "\teta =  " << etaString << std::endl;

}

}
