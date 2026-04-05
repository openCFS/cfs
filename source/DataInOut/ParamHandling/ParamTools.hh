#ifndef PARAMTOOLS_HH_
#define PARAMTOOLS_HH_

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Exception.hh"
#include "Utils/ToolsFull.hh"
#include "MatVec/Matrix.hh"
#include "boost/algorithm/string.hpp"
#include "Domain/CoefFunction/CoefFunction.hh"

namespace CoupledField
{
  /** This class provides an extension to the ParamNode class. It is kept
   * separate from the ParamNode as the inline implementation is convenient
   * for the template function but the ParamNode is better lightweight. */
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
      SplitStringList(node->As<std::string>(), strVec, ' ' );

      ret.Resize( dim1, dim2 );
      ret.Init();
      
      if (strVec.GetSize() != dim1*dim2) 
      {
         EXCEPTION("Wrong size of matrix '" << node->GetName() << "'. It contains of " 
                   << strVec.GetSize() << " entries and should be " << dim1 
                   << " x " << dim2);
      }
      for ( UInt i = 0; i < dim1; i++ ) {
        for ( UInt j = 0; j < dim2; j++ ) {
          boost::algorithm::trim(strVec[i*dim2+j]);
          ret[i][j]= boost::lexical_cast<TYPE>(strVec[i*dim2+j]);

        }
      }
    }

    template <class TYPE>
    static void AsMatrix(PtrParamNode node, Matrix<TYPE>& ret)
    {
      StdVector<std::string> strVec;
      unsigned int dim1 = node->Get("dim1")->As<unsigned int>();
      unsigned int dim2 = node->Get("dim2")->As<unsigned int>();

      if ( dim1 == 1 ) {
        ret.Resize( dim2, dim1 );
      } else {
        ret.Resize( dim1, dim2 );
      }
      ret.Init();

      ParamNodeList rows, cols;
      unsigned int row, col;

      rows = node->GetList("row");
      for ( unsigned int i=0; i<rows.GetSize(); i++ ) {
        row = rows[i]->Get("id")->As<unsigned int>() - 1;
        cols = rows[i]->GetList("col");
        for ( unsigned int j=0; j<cols.GetSize(); j++ ) {
          col = cols[j]->Get("id")->As<unsigned int>() - 1;
          if ( dim1 == 1 ) {
            ret[col][row] = cols[j]->Get("data")->As<TYPE>();
          } else {
            ret[row][col] = cols[j]->Get("data")->As<TYPE>();
          }
        }
      }
    }
    
    //! Read the immediate child of a node and return it as a string vector with given dimension
    static void AsStringTensor(PtrParamNode node, UInt dim, StdVector<std::string>& ret )
    {
      SplitStringList(node->As<std::string>(), ret, ' ' );
      
      // be careful and wrap every entry in braces
      for( UInt i = 0; i < ret.GetSize(); ++i ) {
        ret[i] = Bracket( ret[i] );
      }
      if (ret.GetSize() != dim)  {
        EXCEPTION("Wrong size of tensor '" << node->GetName() << "'. It contains of " 
                  << ret.GetSize() << " entries and should be " << dim ); 
      }
    }

    //! Read the immediate child of a node and return it as a string vector usable to create a
    //! tensor valued coefficient function via CoefFunction::Generate
    static void AsStringVoigtTensor(PtrParamNode node, StdVector<std::string>& ret )
    {
        StdVector<std::string> elems;
        SplitStringList(node->As<std::string>(), elems, ' ' );
        UInt nelems = elems.GetSize();
        if ( nelems == 6 ) {
            ret.Resize(9); // 3x3 components
            ret[0] = Bracket( elems[0] ); // 11
            ret[1] = Bracket( elems[5] ); // 12
            ret[2] = Bracket( elems[4] ); // 13
            ret[3] = Bracket( elems[5] ); // 21 = 12
            ret[4] = Bracket( elems[1] ); // 22
            ret[5] = Bracket( elems[3] ); // 23
            ret[6] = Bracket( elems[4] ); // 31 = 13
            ret[7] = Bracket( elems[3] ); // 32 = 23
            ret[8] = Bracket( elems[2] ); // 33
        }
        else {
            EXCEPTION("Your specified '" << nelems << "' but symmetric tensors from Voigt notation only implemented for 3D (=6 components)!");
        }
    }

    template <class TYPE>
    static void AsVector(PtrParamNode node, Vector<TYPE>& ret)
    {
      StdVector<std::string> strVec;
      SplitStringList(node->As<std::string>(), strVec, ' ' );
      UInt length = strVec.GetSize();
      ret.Resize( length );
      ret.Init();


      for ( UInt i = 0; i < length; i++ ) {
        ret[i]= boost::lexical_cast<TYPE>(strVec[i]);
      }
    }

    static PtrCoefFct AsScalarCoefFct(MathParser* mp, PtrParamNode node) {
        PtrCoefFct ret;
        std::string sR,sI;
        bool complex = false;
        if ( node->Has("real") ) {
            sR = node->Get("real")->As<std::string>();
            //tRef = CoefFunction::Generate( mp_,
           // material->SetScalar(refT->Get("real")->As<std::string>(),MECH_TE_REFTEMPERATURE, Global::REAL );
        }
        if ( node->Has("imag")) {
            sI = node->Get("imag")->As<std::string>();
            complex = true;
        }
        if (complex) {
            ret = CoefFunction::Generate( mp, Global::REAL, sR, sI);
        }
        else {
            ret = CoefFunction::Generate( mp, Global::REAL, sR);
        }
        //Global::COMPLEX
        return ret;
    }

  }; 

} // end of namespace

#endif /*PARAMTOOLS_HH_*/
