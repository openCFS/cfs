// - C++ -
/***************************************************************************
    File        : TriSkeleton.h
    Description : 

 ---------------------------------------------------------------------------
    Begin       : Fri Feb 22 2002
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#ifndef TRISKELETON_H
#define TRISKELETON_H

namespace grd {

struct TriSkeleton {
  public:
  // Constructors
  TriSkeleton() {
    int ch;
    for (ch = 0; ch < 4; ch++)
      child[ch] = 0;

    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        pos[i][j] = 0.0;
  }

  // Destructor
  virtual ~TriSkeleton() { };

  // Method
  bool isRefined() {
    if (child[0] != 0)
      return true;
    else
      return false;
  }

  double pos[3][3];
  TriSkeleton* child[4];
};

}

#endif // TRISKELETON_H