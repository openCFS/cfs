// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "edgeFace.hh"
#include "MatVec/vector.hh"

namespace CoupledField {


  // static variable initialization
  
  char Face::quadBits[4][4] = 
    { { -1,  7, -1,  3,},
      {  6, -1,  2, -1,},
      { -1,  0, -1,  4,},
      {  1, -1,  5, -1,} };
  
  char Face::triaBits[3][3] = 
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
    
  
//    
//    std::cerr << "\n=====================================\n"
//        << " Face Orientation\n"
//        << "=====================================\n";
//    std::cerr << "\tconnect:" << nodes.ToString() << std::endl;
    
    StdVector<UInt> indices( nodes.GetSize() );
    UInt size = nodes.GetSize();
    
    // initialize indices array
    for( UInt i = 0; i < size; i++ ) {
      indices[i] = i;
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
      indices[j] = i;
    }
    // -----------------------

  
    
    // fetch orientation flags
    if( size == 4 ) {
      flags = std::bitset<3>( quadBits[indices[0]][indices[1]]);
//      std::cerr << "\tsorted: " << nodes.ToString() << std::endl;
//      std::cerr << "\tindices: " << indices.ToString() << std::endl;
//     
//      
//      
//      
//      
//      std::cerr << "\tflag: " << flags << " (" << flags.to_ulong() << ")" << std::endl;
//      std::cerr << "\tflag: [0] = " << flags.test(0)
//                           << "[1] = " << flags.test(1)
//                           << "[2] = " << flags.test(2) << std::endl;
      
      // print debug information
      std::string xiString, etaString;
      if( flags.test(2) == true) {
        xiString = flags.test(0) ? "I" :"-I";
        etaString = flags.test(1) ? "II" :"-II";
      }else {
        xiString = flags.test(0) ? "II" :"-II";
        etaString = flags.test(1) ? "I" :"-I";
      }
//      std::cerr << "xi =  " << xiString << std::endl;
//      std::cerr << "eta =  " << etaString << std::endl;

      // switch 3rd and 4th node
      UInt val = nodes[3];
      nodes[3] = nodes[2];
      nodes[3] = val;
      
    } else if( size == 3 ) {
      flags = std::bitset<3>( triaBits[indices[0]][indices[1]]);
    }
    // Check flags for orientation

  }

}
