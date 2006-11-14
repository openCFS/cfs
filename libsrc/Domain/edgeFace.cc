#include "edgeFace.hh"
#include "Utils/vector.hh"

namespace CoupledField {


  // static variable initialization
  
  Char Face::quadBits[4][4] = 
    { { 0, 7, 0, 6,},
      { 3, 0, 2, 0,},
      { 0, 0, 0, 1,},
      { 4, 0, 5, 0,} };
  
  Char Face::triaBits[3][3] = 
    { { 0, 0, 0, },
      { 0, 0, 0, },
      { 0, 0, 0, } };
    
  bool operator<( const Face& a, const Face& b) {

    UInt i = 0; 
    while( i < a.nodes.GetSize()-1 
           && (a.nodes[i] == b.nodes[i]) ) {i++;};
    
    return a.nodes[i] < b.nodes[i];
  }
  
  bool operator<( const Edge& a, const Edge& b ) {
    if( a.nodes[0] < b.nodes[0] ) {
      return true;
    }
    if (a.nodes[0] > b.nodes[0] ) {
      return false;
    }

    if( a.nodes[1] < b.nodes[1] ) {
      return true;
    } else {
      return false;
    }
  }


  void Face::Normalize( std::bitset<3>& flags) {
    ENTER_FCN( "GridCFS:::NormalizeFace" );
    
    Vector<UInt> indices( nodes.GetSize() );
    UInt size = nodes.GetSize();
    
    // initialize indices array
    for( UInt i = 0; i < size; i++ ) {
      indices[i] = i+1;
    }

 
    // Insertion sort algorithm
    // ------------------------
    UInt j, comp;
    
    for (UInt i = 1; i < size; i++) {
      comp = nodes[i];
      j = i;
      while ((j > 0) && (nodes[j-1] > comp)) {
        nodes[j] = nodes[j-1];
        indices[j] = indices[j-1];
        j = j - 1;
      }
      nodes[j] = comp;
      indices[j] = i+1;
    }
    // -----------------------

  
    
    // fetch orientation flags
    if( size == 4 ) {

      flags = std::bitset<3>( quadBits[indices[0]-1][indices[1]-1]);

      //std::cerr << "flag for  [" << indices[0] << "][" 
      //          << indices[1] << "] is " 
      //          <<  (short)quadBits[indices[0]-1][indices[1]-1] << std::endl;

      // switch 3rd and 4th node
      UInt val = nodes[3];
      nodes[3] = nodes[2];
      nodes[3] = val;
      
    } else if( size == 3 ) {
      flags = std::bitset<3>( triaBits[indices[0]-1][indices[1]-1]);
    }
    // Check flags for orientation

  }

}
