#include "EdgeFace.hh"
#include "MatVec/Vector.hh"
#include "Domain/ElemMapping/Elem.hh"

namespace CoupledField {


// static variable initialization


// the following array holds the face flags (see Elem::faceFlags)
// in numeric form. As key (3 indices, ranging from 0..3), we
// take the sorting indices of the face connectivity, i.e. we sort
// the initial connectivity array in ascending order and use the 
// first three indices as key into this array. This array is sparse
// (only 24 entries = 4-noded surface, where every position can
// be permutated = 4*3*2*1 = 24 entries). Invalid entries contain
// a -1, which can be used for consistency checks.
// -----------------------
//  0   1     2    3
// -----------------------
int Face::quadBits[4][4][4] =
{           
 {                 
  {-1,  -1,  -1,  -1}, // (0,0,*)
  {-1,  -1,   7,  23}, // (0,1,*)
  {-1,  15,  -1,  11}, // (0,2,*)
  {-1,  19,   3,  -1}  // (0,3,*)
 },
 {
  {-1,  -1,  22,   6}, // (1,0,*)
  {-1,  -1,  -1,  -1}, // (1,1,*)
  {18,  -1,  -1,   2}, // (1,2,*)
  {14,  -1,  10,  -1}  // (1,3,*)
 },
 {
  {-1,   8,  -1,  12}, // (2,0,*)
  { 0,  -1,  -1,  16}, // (2,1,*)
  {-1,  -1,  -1,  -1}, // (2,2,*)
  { 4,  20,  -1,  -1}  // (2,3,*)
 },
 {
  {-1,   1,  17,  -1}, // (3,0,*)
  { 9,  -1,  13,  -1}, // (3,1,*)
  {21,   5,  -1,  -1}, // (3,2,*)
  {-1,  -1,  -1,  -1}  // (3,3,*)
 }
};

int Face::triaBits[3][3] = 
{ 
  { -1,  4,  0 }, 
  {  1, -1,  6 }, 
  {  5,  2, -1 } 
};


// This array holds for each faceFlag value (5-bit array, but only
// 24 values are used) the permutation index set for
// quadrilateral shaped faces. If this 
// permutation is applied to the nodes of a quad-face, 
// its orientation will match the global one.
const char quadPerm[24][4] = 
{ 
 { 2, 1, 0, 3 }, //  0
 { 3, 0, 1, 2 }, //  1
 { 1, 2, 3, 0 }, //  2
 { 0, 3, 2, 1 }, //  3
 { 2, 3, 0, 1 }, //  4
 { 3, 2, 1, 0 }, //  5
 { 1, 0, 3, 2 }, //  6
 { 0, 1, 2, 3 }, //  7
 // ------------------
 { 2, 1, 0, 3 }, //  8
 { 3, 0, 1, 2 }, //  9
 { 1, 2, 3, 0 }, // 10
 { 0, 3, 2, 1 }, // 11
 { 2, 3, 0, 1 }, // 12
 { 3, 2, 1, 0 }, // 13
 { 1, 0, 3, 2 }, // 14
 { 0, 1, 2, 3 }, // 15
 // ------------------
 { 2, 1, 0, 3 }, // 16
 { 3, 0, 1, 2 }, // 17
 { 1, 2, 3, 0 }, // 18
 { 0, 3, 2, 1 }, // 19
 { 2, 3, 0, 1 }, // 20
 { 3, 2, 1, 0 }, // 21
 { 1, 0, 3, 2 }, // 22
 { 0, 1, 2, 3 }  // 23
};

//  permutation matrix for triangular faces
const char triaPerm[8][3] = 
{ 
  {  0,  2,  1 }, // 0
  {  1,  0,  2 }, // 1
  {  2,  1,  0 }, // 2
  { -1, -1, -1 }, // 3 -> not defined
  {  0,  1,  2 }, // 4
  {  2,  0,  1 }, // 5
  {  1,  2,  0 }, // 6
  { -1, -1, -1 } // 7 -> not defined
};

void Face::GetSortedIndices( StdVector<UInt>& sorted, 
                             const StdVector<UInt>& unsorted,
                             UInt numVertices,
                             const std::bitset<5>& flags ) {

  const Integer num = flags.to_ulong();
  // indices contains the
  if( numVertices == 4 ) {
    sorted.Resize( 4 );
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

void Face::Normalize(std::bitset<5>& flags, StdVector<unsigned int>& nodes) const
{
  assert(nodes.GetSize() <= 4);
  unsigned int numNodes = nodes.GetSize();
  std::array<unsigned int, 4> indices_array; // we need to sort the face connectivity to get the correct flag
  StdVector<unsigned int> indices(indices_array.data(), numNodes); // might be 3 or 4

  // copy unsorted node Vector but as we call normalize quite often, avoid dynamic allocation
  std::array<unsigned int, 4> unsorted_array; // we need to sort the face connectivity to get the correct flag
  StdVector<UInt> unsorted(unsorted_array);
  unsorted.Replace(nodes.GetPointer(), numNodes, COPY); // copy given nodes as unsorted array
    
  // initialize indices array
  for( UInt i = 0; i < numNodes; i++ ) {
    indices[i] = i;
  }

  // -------------------------
  // Insertion sort algorithm
  // ------------------------
  UInt j, comp;

  for( UInt i = 1; i < numNodes; i++ ) {
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
  if( nodes.GetSize() == 4 ) {
    // security check: we must always get a positive number (= valid bitset)
    assert( quadBits[indices[0]][indices[1]][indices[2]] >= 0 ); 
    
    // obtain permutation bit from the unique two first permutation indices
    flags = std::bitset<5>( quadBits[indices[0]][indices[1]][indices[2]] );
  }
  else if( nodes.GetSize() == 3 ) {
    flags = std::bitset<5>( triaBits[indices[0]][indices[1]] );
  }
  
  
  // re-sort in the end the facenodes
  std::array<unsigned int, 4> sorted_array; // we need to sort the face connectivity to get the correct flag
  StdVector<UInt> sorted(sorted_array);
  this->GetSortedIndices( sorted, unsorted, numNodes, flags ); // sets also flags
  for( UInt i = 0; i < numNodes; ++i ) {
    nodes [i] = sorted[i] + 1;
  }
  
  // ==============================================
  // UNCOMMENT THE FOLLOWING SECTION FOR DETAILED 
  // DEBUGGING INFORMATION
  // ==============================================
//  std::cerr 
//  << "\n=====================================\n"
//  << " Face Orientation\n" 
//  << "=====================================\n";
//  std::cerr << "\tconnect:" << unsorted.ToString() << std::endl;
//  
//  // Check flags for orientation
//  std::cerr << "\tsorted: " << nodes.ToString() << std::endl;
//  std::cerr << "\tindices: " << indices.ToString() << std::endl;
//  std::cerr << "\tflag: " << flags << " (" << flags.to_ulong() << ")"
//      << std::endl;
//  std::cerr << "\tflag: [0] = " << flags.test( 0 ) << std::endl
//      << "\t      [1] = " << flags.test( 1 ) 
//        << std::endl << "\t      [2] = "
//      << flags.test( 2 ) << std::endl;
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

} // end of namespace
