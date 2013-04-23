#ifndef PARAMTOOLS_HH_
#define PARAMTOOLS_HH_

// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/exception.hh"
#include "Utils/tools.hh"
#include "MatVec/matrix.hh"
#include "boost/lexical_cast.hpp"

namespace CoupledField
{
  /** This class provides an extension tt the ParamNode class. It is kept
   * seperate from the ParamNode as the inline implementation is convenient
   * for the template function but the ParamNode is better lightwight. Actually
   * it is also a necessary to use the ParamNode both in OLAS and CFS++ */
  class ParamTools
  {
    public:
    
    
    /** gets the only direct child which has the name and stores its value in ret as Matrix<br>
     * Is only valid, if the corresponding GetList() would return one value. Exception if 0 or greater 1
     * ,if the found value can not be converted into type of return value or if the dimensions
     * of of the matrix provided do not match.<br>
     * If no matching element is found an exception will be thrown if throwException is set to true. 
     * Otherwise the function will return silently and the original value in ret will be kept.
     * example: "optimization" is a complex element which is a direct child of the root: param.Get("optimization") 
     * @throws exception if there is not such a direct child, e.g. if this is a leaf node OR if there more than only
     * one of such elements (e.g. simple xml elements). */     

    /** Interpretes the string content of this ParamNode as Matrix of dimension dim1 x dim2. Example:
     * <pre>
      <elasticity>
        <tensor dim1="6">
          <real>
            1.65682E+11 1.84091E+10 1.84091E+10 0.00000E+00 0.00000E+00 0.00000E+00
            1.84091E+10 1.65682E+11 1.84091E+10 0.00000E+00 0.00000E+00 0.00000E+00
            1.84091E+10 1.84091E+10 1.65682E+11 0.00000E+00 0.00000E+00 0.00000E+00
            0.00000E+00 0.00000E+00 0.00000E+00 7.36364E+10 0.00000E+00 0.00000E+00
            0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 7.36364E+10 0.00000E+00
            0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 7.36364E+10
          </real>
        </tensor>
      </elasticity></pre> 
     * Call this like the following:<pre>
     * PtrParamNode elast = pn->Get("elasticity");
     * PtrParamNode tensor = elast->Get("tensor", "dim1", "6");
     * Matrix<double> mat(6,6);
     * ParamTools::AsTensor<double>(tensor->Get("real"), 6, 6, mat);</pre> */
    template <class TYPE>
    static void AsTensor(PtrParamNode node, unsigned int dim1, unsigned int dim2, Matrix<TYPE>& ret)
    {
      StdVector<std::string> strVec;

      // do not use cfs-implemented SplitStringList() as it cannot handle newline separated SAX parsed data
      SplitStringListWhitespace(node->As<std::string>(), strVec);

      ret.Resize( dim1, dim2 );
      ret.Init();
      
      if (strVec.GetSize() != dim1*dim2) 
      {
         EXCEPTION("Wrong size of matrix '" << node->GetName() << "'. It contains "
                   << strVec.GetSize() << " entries and should be " << dim1 
                   << " x " << dim2);
      }

      for ( UInt i = 0; i < dim1; i++ ) {
        for ( UInt j = 0; j < dim2; j++ ) {
          ret[i][j]= boost::lexical_cast<TYPE>(strVec[i*dim2+j]);
        }
      }
    }

  }; 

} // end of namespace

#endif /*PARAMTOOLS_HH_*/
