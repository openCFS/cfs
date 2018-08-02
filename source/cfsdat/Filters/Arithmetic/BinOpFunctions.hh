// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     BinOpFilter.hh
 *       \brief    <Description>
 *
 *       \date     Apr 13, 2018
 *       \author   mtautz
 */
//================================================================================================

#ifndef SOURCE_CFSDAT_FILTERS_ARITHMETIC_BINOPFUNCTIONS_HH_
#define SOURCE_CFSDAT_FILTERS_ARITHMETIC_BINOPFUNCTIONS_HH_

#include "cfsdat/Filters/BaseFilter.hh"
#include <boost/bimap.hpp>

namespace CFSDat {
  
  namespace BinOpFunctions {
    
    template<typename DATA_TYPE>
    struct Operand {
      
      ResultManager* resultManager;
      
      uuids::uuid id;
      
      Vector<DATA_TYPE>& constant;
      
      Operand(ResultManager* resultmMan, uuids::uuid redid, Vector<DATA_TYPE>& consta):
        resultManager(resultmMan),
        id(redid),
        constant(consta) {
      }

      bool IsConstant() {
        return constant.GetSize() > 0;
      }
      
      bool IsStatic() {
        return IsConstant() || resultManager->IsStatic(id);
      }

      bool IsScalar() {
        if (IsConstant()) {
          return constant.GetSize() == 1;
        }
        return resultManager->GetEntryType(id) == ExtendedResultInfo::SCALAR;
      }
      
      void SetScalar() {
        if (IsConstant()) {
          EXCEPTION("Can not set Entrytype for constant");
        }
        StdVector<std::string> newDofNames;
        newDofNames.Push_back("-");
        resultManager->SetDofNames(id, newDofNames);
        resultManager->SetEntryType(id, ExtendedResultInfo::SCALAR);
      }

      bool IsVector() {
        if (IsConstant()) {
          return constant.GetSize() == 2 || constant.GetSize() == 3;
        }
        return resultManager->GetEntryType(id) == ExtendedResultInfo::VECTOR;
      }
      
      void SetVector(UInt numVecDofs) {
        if (IsConstant()) {
          EXCEPTION("Can not set Entrytype for constant");
        }
        StdVector<std::string> newDofNames;
        newDofNames.Resize(numVecDofs);
        newDofNames[0] = "x";
        newDofNames[1] = "y";
        if (numVecDofs == 3) {
          newDofNames[2] = "z";
        }
        resultManager->SetDofNames(id, newDofNames);
        resultManager->SetEntryType(id, ExtendedResultInfo::VECTOR);
      }

      bool IsTensor() {
        if (IsConstant()) {
          return constant.GetSize() == 4 || constant.GetSize() == 9;
        }
        return resultManager->GetEntryType(id) == ExtendedResultInfo::TENSOR;
      }
      
      void SetTensor(UInt numVecDofs) {
        if (IsConstant()) {
          EXCEPTION("Can not set Entrytype for constant");
        }
        StdVector<std::string> newDofNames;
        if (numVecDofs == 3) {
          newDofNames.Resize(9);
          newDofNames[0] = "xx";
          newDofNames[1] = "xy";
          newDofNames[2] = "xz";
          newDofNames[3] = "yx";
          newDofNames[4] = "yy";
          newDofNames[5] = "yz";
          newDofNames[6] = "zx";
          newDofNames[7] = "zy";
          newDofNames[8] = "zz";
        } else if (numVecDofs == 2) {
          newDofNames.Resize(2);
          newDofNames[0] = "xx";
          newDofNames[1] = "xy";
          newDofNames[2] = "yx";
          newDofNames[3] = "yy";
        }
        resultManager->SetDofNames(id, newDofNames);
        resultManager->SetEntryType(id, ExtendedResultInfo::TENSOR);
      }

      UInt GetNumDofs() {
        if (IsConstant()) {
          return constant.GetSize();
        }
        StdVector<std::string> dofnames = resultManager->GetDofNames(id);
        return dofnames.GetSize();
      }
      
      UInt GetVecNumDofs() {
        if (IsTensor()) {
          return GetNumDofs() == 9 ? 3 : 2;
        } else if (IsVector()) {
          return GetNumDofs();
        }
        return 0;
      }
      
    };
  
    //! elementary addition operator
    template<typename DATA_TYPE>
    struct Plus {
      static inline DATA_TYPE Apply(DATA_TYPE& inA, DATA_TYPE& inB) {
        return inA + inB;
      }
    };

    //! elementary subtraction operator
    template<typename DATA_TYPE>
    struct Minus {
      static inline DATA_TYPE Apply(DATA_TYPE& inA, DATA_TYPE& inB) {
        return inA - inB;
      }
    };

    //! elementary product operator
    template<typename DATA_TYPE>
    struct Mult {
      static inline DATA_TYPE Apply(DATA_TYPE& inA, DATA_TYPE& inB) {
        return inA * inB;
      }
    };

    //! elementary divison operator
    template<typename DATA_TYPE>
    struct Div {
      static inline DATA_TYPE Apply(DATA_TYPE& inA, DATA_TYPE& inB) {
        return inA / inB;
      }
    };

    //! entrywise operations
    template<typename DATA_TYPE, class OPERATOR, UInt NUMDOFS>
    struct EntrywiseOperator {
      static inline void Apply(DATA_TYPE* out, DATA_TYPE* inA, DATA_TYPE* inB) {
        for (UInt i = 0; i < NUMDOFS; i++) {
          out[i] = OPERATOR::Apply(inA[i], inB[i]);
        }
      }
    };
    
    //! entry - scalar operations
    template<typename DATA_TYPE, class OPERATOR, UInt NUMDOFS>
    struct EntryScalarOperator {
      static inline void Apply(DATA_TYPE* out, DATA_TYPE* inA, DATA_TYPE* inB) {
        for (UInt i = 0; i < NUMDOFS; i++) {
          out[i] = OPERATOR::Apply(inA[i], inB[0]);
        }
      }
    };

    //! scalar - entry operations
    template<typename DATA_TYPE, class OPERATOR, UInt NUMDOFS>
    struct ScalarEntryOperator {
      static inline void Apply(DATA_TYPE* out, DATA_TYPE* inA, DATA_TYPE* inB) {
        for (UInt i = 0; i < NUMDOFS; i++) {
          out[i] = OPERATOR::Apply(inA[0], inB[i]);
        }
      }
    };
    
    //! scalar operations
    template<typename DATA_TYPE, class OPERATOR>
    struct ScalarOperator {
      static inline void Apply(DATA_TYPE* out, DATA_TYPE* inA, DATA_TYPE* inB) {
        out[0] = OPERATOR::Apply(inA[0], inB[0]);
      }
    };

    //! scalar product
    template<typename DATA_TYPE, UInt NUMDOFS>
    struct ScalarProduct {
      static inline void Apply(DATA_TYPE* out, DATA_TYPE* inA, DATA_TYPE* inB) {
        out[0] =  inA[0] * inB[0];
        for (UInt i = 1; i < NUMDOFS; i++) {
          out[0] += inA[i] * inB[i];
        }
      }
    };
    
    //! scalar product
    template<typename DATA_TYPE, UInt VECNUMDOFS>
    struct CrossProduct {
      //! three dimensional
      static inline void Apply(DATA_TYPE* out, DATA_TYPE* inA, DATA_TYPE* inB) {
        if (VECNUMDOFS == 3) {
          out[0] = inA[1] * inB[2] - inA[2] * inB[1];
          out[1] = inA[2] * inB[0] - inA[0] * inB[2];
          out[2] = inA[0] * inB[1] - inA[1] * inB[0];
        } else {
          out[0] = inA[0] * inB[1] - inA[1] * inB[0]; // two dimensional
        }
      }
    };
    
    //! dyadic product
    template<typename DATA_TYPE, UInt VECNUMDOFS>
    struct DyadicProduct {
      static inline void Apply(DATA_TYPE* out, DATA_TYPE* inA, DATA_TYPE* inB) {
        UInt idx = 0;
        for (UInt i = 0; i < VECNUMDOFS; i++) {
          for (UInt j = 0; j < VECNUMDOFS; j++) {
            out[idx] = inA[i] * inB[j];
            idx++;
          }
        }
      }
    };
    

    //! matrix * vector product
    template<typename DATA_TYPE, UInt VECNUMDOFS>
    struct MatrixVectorProduct {
      static inline void Apply(DATA_TYPE* out, DATA_TYPE* inA, DATA_TYPE* inB) {
        for (UInt i = 0; i < VECNUMDOFS; i++) {
          out[i] = inA[i * 3] * inB[0];
          for (UInt j = 1; j < VECNUMDOFS; j++) {
            out[i] += inA[i * 3 + j] * inB[j];
          }
        }
      }
    };

    //! matrix * matrix product
    template<typename DATA_TYPE, UInt VECNUMDOFS>
    struct MatrixProduct {
      static inline void Apply(DATA_TYPE* out, DATA_TYPE* inA, DATA_TYPE* inB) {
        UInt idx = 0;
        for (UInt i = 0; i < VECNUMDOFS; i++) {
          for (UInt j = 0; j < VECNUMDOFS; j++) {
            out[idx] = inA[i * 3] * inB[j];
            for (UInt k = 1; k < VECNUMDOFS; k++) {
              out[idx] += inA[i * 3 + k] * inB[j + k * 3];
            }
            idx++;
          }
        }
      }
    };

    //! struct for a pointer to a function processing result vectors using binary operators
    template<typename DATA_TYPE>
    struct BinOpFctStruct {
      typedef void (*BinOpFctPtr)(Vector<DATA_TYPE>&,Vector<DATA_TYPE>&,Vector<DATA_TYPE>&,UInt);
    };

    //! template function used for the pointers to the functions
    //template<typename DATA_TYPE>
    //struct BinOpVectFuncStruct {
    //  template<class OPERATOR, class IDXOUT, class IDXA, class IDXB>
      template<typename DATA_TYPE, class OPERATOR, UInt IDXOUT, UInt IDXA, UInt IDXB>
      static inline void BinOpVectFunc(Vector<DATA_TYPE>& out, Vector<DATA_TYPE>& inA, Vector<DATA_TYPE>& inB, UInt size) {
        #pragma omp parallel for num_threads(CFS_NUM_THREADS)
        for (UInt i = 0; i < size; i++) {
          OPERATOR::Apply(&out[IDXOUT * i], &inA[IDXA * i], &inB[IDXB * i]);
        }
      };
    //};

    template<typename DATA_TYPE, class OPERATOR, UInt IDXOUT, UInt IDXA, UInt IDXB>
    typename BinOpFctStruct<DATA_TYPE>::BinOpFctPtr GetFunction() {
      //return &BinOpVectFuncStruct<DATA_TYPE>::BinOpVectFunc<OPERATOR,IDXOUT,IDXA,IDXB>;
      return &BinOpVectFunc<DATA_TYPE,OPERATOR,IDXOUT,IDXA,IDXB>;
    }
    
    template<typename DATA_TYPE, UInt VECNUMDOFS, class OPERATOR, UInt IDXOUT, UInt IDXA>
    typename BinOpFctStruct<DATA_TYPE>::BinOpFctPtr GetFunction(Operand<DATA_TYPE>& opB) {
      if(opB.IsConstant()) {
        return GetFunction<DATA_TYPE, OPERATOR, IDXOUT, IDXA, 0 >();
      }
      if (opB.IsTensor()) {
        return GetFunction<DATA_TYPE, OPERATOR, IDXOUT, IDXA, VECNUMDOFS * VECNUMDOFS >();
      } else if (opB.IsVector()) {
        return GetFunction<DATA_TYPE, OPERATOR, IDXOUT, IDXA, VECNUMDOFS >();
      }
      return GetFunction<DATA_TYPE, OPERATOR, IDXOUT, IDXA, 1 >();
    };
    
    template<typename DATA_TYPE, UInt VECNUMDOFS, class OPERATOR, UInt IDXOUT>
    typename BinOpFctStruct<DATA_TYPE>::BinOpFctPtr GetFunction(Operand<DATA_TYPE>& opA, Operand<DATA_TYPE>& opB) {
      if(opA.IsConstant()) {
        return GetFunction<DATA_TYPE, VECNUMDOFS, OPERATOR, IDXOUT, 0 >(opB);
      }
      if (opA.IsTensor()) {
        return GetFunction<DATA_TYPE, VECNUMDOFS, OPERATOR, IDXOUT, VECNUMDOFS * VECNUMDOFS >(opB);
      } else if (opA.IsVector()) {
        return GetFunction<DATA_TYPE, VECNUMDOFS, OPERATOR, IDXOUT, VECNUMDOFS >(opB);
      }
      return GetFunction<DATA_TYPE, VECNUMDOFS, OPERATOR, IDXOUT, 1 >(opB);
    };
    
    template<typename DATA_TYPE, UInt VECNUMDOFS, class OPERATOR>
    typename BinOpFctStruct<DATA_TYPE>::BinOpFctPtr GetFunction(Operand<DATA_TYPE>& opOut, Operand<DATA_TYPE>& opA, Operand<DATA_TYPE>& opB) {
      if (opOut.IsTensor()) {
        return GetFunction<DATA_TYPE, VECNUMDOFS, OPERATOR, VECNUMDOFS * VECNUMDOFS >(opA, opB);
      } else if (opOut.IsVector()) {
        return GetFunction<DATA_TYPE, VECNUMDOFS, OPERATOR, VECNUMDOFS >(opA, opB);
      }
      return GetFunction<DATA_TYPE, VECNUMDOFS, OPERATOR, 1 >(opA, opB);
    };
    
    template<typename DATA_TYPE, UInt VECNUMDOFS, class ELEMOPERATOR>
    typename BinOpFctStruct<DATA_TYPE>::BinOpFctPtr GetEntrywiseFunction(Operand<DATA_TYPE>& opOut, Operand<DATA_TYPE>& opA, Operand<DATA_TYPE>& opB) {
      if (opA.IsScalar()) {
        if (opB.IsScalar()) {
          opOut.SetScalar();
          return GetFunction<DATA_TYPE, VECNUMDOFS, ScalarOperator<DATA_TYPE, ELEMOPERATOR> >(opOut, opA, opB);
        } else if (opB.IsVector()) {
          opOut.SetVector(VECNUMDOFS);
          return GetFunction<DATA_TYPE, VECNUMDOFS, ScalarEntryOperator<DATA_TYPE, ELEMOPERATOR, VECNUMDOFS> >(opOut, opA, opB);
        } else if (opB.IsTensor()) {
          opOut.SetTensor(VECNUMDOFS);
          return GetFunction<DATA_TYPE, VECNUMDOFS, ScalarEntryOperator<DATA_TYPE, ELEMOPERATOR, VECNUMDOFS * VECNUMDOFS> >(opOut, opA, opB);
        }
      } else if (opA.IsVector()){
        if (opB.IsScalar()) {
          opOut.SetVector(VECNUMDOFS);
          return GetFunction<DATA_TYPE, VECNUMDOFS, EntryScalarOperator<DATA_TYPE, ELEMOPERATOR, VECNUMDOFS> >(opOut, opA, opB);
        } else if (opB.IsVector()) {
          opOut.SetVector(VECNUMDOFS);
          return GetFunction<DATA_TYPE, VECNUMDOFS, EntrywiseOperator<DATA_TYPE, ELEMOPERATOR, VECNUMDOFS> >(opOut, opA, opB);
        } else if (opB.IsTensor()) {
          EXCEPTION("Can not perform binary operation for oparands vector and tensor");
        }
      }
      if (opB.IsScalar()) {
        opOut.SetTensor(VECNUMDOFS);
        return GetFunction<DATA_TYPE, VECNUMDOFS, EntryScalarOperator<DATA_TYPE, ELEMOPERATOR, VECNUMDOFS * VECNUMDOFS> >(opOut, opA, opB);
      } else if (opB.IsVector()) {
        EXCEPTION("Can not perform binary operation for oparands vector and tensor (except matrix * vector)");
      }
      opOut.SetTensor(VECNUMDOFS);
      return GetFunction<DATA_TYPE, VECNUMDOFS, EntrywiseOperator<DATA_TYPE, ELEMOPERATOR, VECNUMDOFS * VECNUMDOFS> >(opOut, opA, opB);
    }
    
    template<typename DATA_TYPE, UInt VECNUMDOFS>
    typename BinOpFctStruct<DATA_TYPE>::BinOpFctPtr GetFunction(std::string opType, std::string multType, 
                                   Operand<DATA_TYPE>& opOut, Operand<DATA_TYPE>& opA, Operand<DATA_TYPE>& opB) {
      if (opType == "mult") {
        if (opA.IsTensor()) { // first operand is tensor
          if (opB.IsTensor()) { // tensor * tensor product
            if (multType == "default" || multType == "matrix") {
              opOut.SetTensor(VECNUMDOFS);
              return GetFunction<DATA_TYPE, VECNUMDOFS, MatrixProduct<DATA_TYPE, VECNUMDOFS> >(opOut, opA, opB);
            } else if (multType == "inner" || multType == "dot" || multType == "scalar") {
              opOut.SetScalar();
              return GetFunction<DATA_TYPE, VECNUMDOFS, ScalarProduct<DATA_TYPE, VECNUMDOFS * VECNUMDOFS> >(opOut, opA, opB);
            } else if (multType == "matrixvector" || multType == "cross" || multType == "outer" || multType == "dyadic") {
              EXCEPTION("multType " << multType << " not supported/implemented for matrix * matrix product");
            }
          } else if (opB.IsVector()) { // matrix * vector product
            if (multType != "default" && multType != "matrix" && multType != "matrixvector") {
              EXCEPTION("Invalid multType given for matrix * vector product, use either default, matrix of matrixvector");
            }
            opOut.SetVector(VECNUMDOFS);
            return GetFunction<DATA_TYPE, VECNUMDOFS, MatrixVectorProduct<DATA_TYPE, VECNUMDOFS> >(opOut, opA, opB);
          }            
        } else if (opA.IsVector() && opB.IsVector()) { // vector * vector product 
          if (multType == "default" || multType == "inner" || multType == "dot" || multType == "scalar") {
            opOut.SetScalar();
            return GetFunction<DATA_TYPE, VECNUMDOFS, ScalarProduct<DATA_TYPE, VECNUMDOFS> >(opOut, opA, opB);
          } else if (multType == "dyadic" || multType == "outer") {
            opOut.SetTensor(VECNUMDOFS);
            return GetFunction<DATA_TYPE, VECNUMDOFS, DyadicProduct<DATA_TYPE, VECNUMDOFS> >(opOut, opA, opB);
          } else if (multType == "cross") {
            if (VECNUMDOFS == 3) {
              opOut.SetVector(VECNUMDOFS);
            } else if (VECNUMDOFS == 2) {
              opOut.SetScalar();
            } else {
              EXCEPTION("Cross product only defined in R3 and R2");
            }
            return GetFunction<DATA_TYPE, VECNUMDOFS, CrossProduct<DATA_TYPE, VECNUMDOFS> >(opOut, opA, opB);
          } else if (multType == "matrixvector" || multType == "matrix") {
            EXCEPTION("multType " << multType << " not supported/implemented for matrix * matrix product");
          }
        }
      } else if (opType == "plus") {
        bool scalarTensorReplace = false;
        bool replaceB = false;
        if (opA.IsTensor() && opB.IsScalar() && opB.IsConstant()) {
          scalarTensorReplace = true;
          replaceB = true;
        } else if (opB.IsTensor() && opA.IsScalar() && opA.IsConstant()) {
          scalarTensorReplace = true;
          replaceB = false;
        }
        if (scalarTensorReplace) {
          Vector<DATA_TYPE>& constant = replaceB ? opB.constant : opA.constant;
          DATA_TYPE value = constant[0];
          if (VECNUMDOFS == 2) {
            constant.Resize(4);
            constant.Init(0.0);
            constant[0] = value;
            constant[2] = value;
          } else if (VECNUMDOFS == 3) {
            constant.Resize(9);
            constant.Init(0.0);
            constant[0] = value;
            constant[4] = value;
            constant[8] = value;
          }
        }
      }
      if (opType == "plus") {
        return GetEntrywiseFunction<DATA_TYPE, VECNUMDOFS, Plus<DATA_TYPE> >(opOut, opA, opB);
      } else if (opType == "minus") {
        return GetEntrywiseFunction<DATA_TYPE, VECNUMDOFS, Minus<DATA_TYPE> >(opOut, opA, opB);
      } else if (opType == "div") {
        return GetEntrywiseFunction<DATA_TYPE, VECNUMDOFS, Div<DATA_TYPE> >(opOut, opA, opB);
      }
      if (opType != "mult") {
        EXCEPTION("Invalid operator type given");
      }
      return GetEntrywiseFunction<DATA_TYPE, VECNUMDOFS, Mult<DATA_TYPE> >(opOut, opA, opB);
    };

    template<typename DATA_TYPE>
    typename BinOpFctStruct<DATA_TYPE>::BinOpFctPtr GetBinOpFunction(std::string opType, std::string multType, 
                                   Operand<DATA_TYPE>& opOut, Operand<DATA_TYPE>& opA, Operand<DATA_TYPE>& opB) {
      UInt vecNumDofA = opA.GetVecNumDofs();
      UInt vecNumDofB = opB.GetVecNumDofs();
      if (vecNumDofA > 0 && vecNumDofB > 0 && vecNumDofA != vecNumDofB) {
        EXCEPTION("vector/matrix sizes for binary operator not matching");
      } else if (vecNumDofA == 2 || vecNumDofB == 2) {
        return GetFunction<DATA_TYPE, 2>(opType, multType, opOut, opA, opB);
      }
      return GetFunction<DATA_TYPE, 3>(opType, multType, opOut, opA, opB);
    }

  }
  
}

#endif /* SOURCE_CFSDAT_FILTERS_ARITHMETIC_BINOPFUNCTIONS_HH_ */
