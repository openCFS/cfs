// =====================================================================================
// 
//       Filename:  SimOutputParsed.hh
// 
//    Description:  This class defines the general output writer for the pos file format 
//                  GMSH. It offers the possibility to store higher order element results.
//                  These can later be viewed by gmsh using an adaptive visualization
//                  unfortunetly the format is not very well suited for a high number of 
//                  elements.
//                  The main task is to find the appropriate interpolation matrices. These
//                  have to be computed in the reference elements. and describe basically 
//                  the monomial representation of the shape function.
//                  Here is an example pos file for two second order lagrange quads with 
//                  are displayed adaptively
//
//                  x1 = 1.01;
//                  x2 = 2.50;  
//
//                  View "Quad" {
//                  SQ(
//                  0.0,0.0,0.0,
//                  1.0,0.0,0.0,
//                  1.0,1.0,0.0,
//                  0.0,1.0,0.0
//                  ){
//                  0.0, 0.0, x1, 01.0, 10.0, 1.0, 0.0, 6.0, 7.1,
//                  0.0, 0.0, x1, 01.0, 1.0, 5.0, 0.0, 1.0, 2.1
//                  };
//                  SQ(
//                  2.0,0.0,0.0,
//                  3.0,0.0,0.0,
//                  3.0,1.0,0.0,
//                  2.0,1.0,0.0
//                  ){
//                  0.0, 0.0, x2, 0.0, 0.0, 0.0, 0.0, 0.0, 10.5,
//                  0.0, 0.0, x1, 01.0, 1.0, 5.0, 0.0, 1.0, 2.1
//                  };
//                  TIME { 1, 25};
//                  INTERPOLATION_SCHEME
//                  {
//                   {0.000, 0.000, 0.000, 0.000, 0.250,-0.250, 0.000,-0.250, 0.250},
//                   {0.000, 0.000, 0.000, 0.000,-0.250, 0.250, 0.000,-0.250, 0.250},
//                   {0.000, 0.000, 0.000, 0.000, 0.250, 0.250, 0.000, 0.250, 0.250},
//                   {0.000, 0.000, 0.000, 0.000,-0.250,-0.250, 0.000, 0.250, 0.250},
//                   {0.000,-0.500, 0.500, 0.000, 0.000, 0.000, 0.000, 0.500,-0.500},
//                   {0.000, 0.000, 0.000, 0.500, 0.000,-0.500, 0.500, 0.000,-0.500},
//                   {0.000, 0.500, 0.500, 0.000, 0.000, 0.000, 0.000,-0.500,-0.500},
//                   {0.000, 0.000, 0.000,-0.500, 0.000, 0.500, 0.500, 0.000,-0.500},
//                   {1.000, 0.000,-1.000, 0.000, 0.000, 0.000,-1.000, 0.000, 1.000}
//                  }
//                  { {0,0,0},
//                    {0,1,0},
//                    {0,2,0},
//                    {1,0,0},
//                    {1,1,0},
//                    {1,2,0},
//                    {2,0,0},
//                    {2,1,0},
//                    {2,2,0}
//                  };
//                  };
//                  View[0].AdaptVisualizationGrid = 1;
//                  View[0].MaxRecursionLevel = 8;
//                  View[0].TargetError = 0.0005;
// 
//        Version:  1.0
//        Created:  11/19/2011 04:48:27 PM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
// 
// =====================================================================================

#ifndef FILE_POSWRITER
#define FILE_POSWRITER

#include "DataInOut/SimOutput.hh"
#include "General/Enum.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "MatVec/Matrix.hh"
#include "MatVec/Vector.hh"
#include "FeBasis/BaseFE.hh"
#include "FeBasis/FeFunctions.hh"
#include "FeBasis/FeSpace.hh"
#include <list>

#include <map>

namespace CoupledField {



  typedef enum{
    //type  typeLong                   #list-of-coords     #list-of-values
    // --------------------------------------------------------------------
     UNDEF = 0,
     SP=1,   //Scalar point                     3            1  * nb-time-steps
     VP=2,   //Vector point                     3            3  * nb-time-steps
     TP=3,   //Tensor point                     3            9  * nb-time-steps
     SL=4,   //Scalar line                      6            2  * nb-time-steps
     VL=5,   //Vector line                      6            6  * nb-time-steps
     TL=6,   //Tensor line                      6            18 * nb-time-steps
     ST=7,   //Scalar triangle                  9            3  * nb-time-steps
     VT=8,   //Vector triangle                  9            9  * nb-time-steps
     TT=9,   //Tensor triangle                  9            27 * nb-time-steps
     SQ=10,   //Scalar quadrangle                12           4  * nb-time-steps
     VQ=11,   //Vector quadrangle                12           12 * nb-time-steps
     TQ=12,   //Tensor quadrangle                12           36 * nb-time-steps
     SS=13,   //Scalar tetrahedron               12           4  * nb-time-steps
     VS=14,   //Vector tetrahedron               12           12 * nb-time-steps
     TS=15,   //Tensor tetrahedron               12           36 * nb-time-steps
     SH=16,   //Scalar hexahedron                24           8  * nb-time-steps
     VH=17,   //Vector hexahedron                24           24 * nb-time-steps
     TH=18,   //Tensor hexahedron                24           72 * nb-time-steps
     SI=19,   //Scalar prism                     18           6  * nb-time-steps
     VI=20,   //Vector prism                     18           18 * nb-time-steps
     TI=21,   //Tensor prism                     18           54 * nb-time-steps
     SY=22,   //Scalar pyramid                   15           5  * nb-time-steps
     VY=23,   //Vector pyramid                   15           15 * nb-time-steps
     TY=24,   //Tensor pyramid                   15           45 * nb-time-steps
     T2=25,   //2D text                          3            arbitrary
     T3=26    //3D text                          4            arbitrary
  }GmeshParsedElemTypes;

  static Enum<GmeshParsedElemTypes> GmeshParsedElemTypesEnum;

  struct ElemInterpolation{
      std::list<UInt> eIndices_;
      GmeshParsedElemTypes curType;
      UInt saveBegin_;
      UInt saveInc_;
      UInt saveEnd_;
      Matrix<Double> coefs_;
      Matrix<Integer> exponents_;
      bool matrixComputed_;

      ElemInterpolation(){
        matrixComputed_ = false;
      }
  };

  class SimOutputParsed : public SimOutput {

    public:

      SimOutputParsed(std::string fileName, PtrParamNode inputNode,
                      PtrParamNode infoNode, bool isRestart );
      virtual ~SimOutputParsed();

      virtual void Init(Grid* ptGrid, bool printGridOnly );

      //! Begin multisequence step
      virtual void BeginMultiSequenceStep( UInt step,
                                           BasePDE::AnalysisType type,
                                           UInt numSteps);

      //! Register result (within one multisequence step)
      virtual void RegisterResult( shared_ptr<BaseResult> sol,
                                   UInt saveBegin, UInt saveInc,
                                   UInt saveEnd,
                                   bool isHistory );

      //! Begin single analysis step
      virtual void BeginStep( UInt stepNum, Double stepVal );

      //! Add result to current step
      virtual void AddResult( shared_ptr<BaseResult> sol );

      //! End single analysis step
      virtual void FinishStep( );

      //! End multisequence step
      virtual void FinishMultiSequenceStep( );

      //! Finalize the output
      virtual void Finalize();

    protected:



    private:

      void FillElemMapping();

      void PrepareResultFile(shared_ptr<BaseResult> sol);

      void WriteDummyResults(const Elem* elem,ElemInterpolation& eInterpol,
                                shared_ptr<BaseFeFunction> feFnc, std::fstream* out);

      void GetInterpolationString(std::string & interp, ElemInterpolation eInterpol);

      void PutVarsToResultFile(std::fstream* outfile, std::fstream* infile, std::string vars,long& destination);

      std::map< Elem::ShapeType, GmeshParsedElemTypes > eTypeMap_;

      typedef std::map<RegionIdType, std::list<ElemInterpolation> > RegionInterpMap;
      typedef std::map<SolutionType, RegionInterpMap > SolutionRegionMap;

      SolutionRegionMap mySolutions_;

      Grid* myGrid_;

      UInt numSteps_;

      UInt curStep_;

      bool firstStep_;

      std::map<SolutionType, std::fstream*> outfiles_;
      std::map<SolutionType, std::fstream*> infiles_;

      Double curTime_;

      std::map<SolutionType,long> resFilePosition_;

  };
}
#endif
