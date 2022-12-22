/*
 * CoefFunctionConvolution.cc
 *
 *  Created on: 29 Jan 2016
 *      Author: INTERN\ftoth
 */

#ifndef COEFFUNCTIONCONVOLUTION_HH
#define COEFFUNCTIONCONVOLUTION_HH

#include "CoefFunction.hh"
#include "FeBasis/FeFunctions.hh"
#include "FeBasis/FeSpace.hh"
#include "Forms/Operators/BaseBOperator.hh"
#include "Domain/ElemMapping/EntityLists.hh"


namespace CoupledField {

  // Forward declatarions
  //class BaseFE;
  //class EntityList;

  //! Provide a coefficient function for convolution
  class CoefFunctionConvolution : public CoefFunction {

  public:

    //! Constructor
    //CoefFunctionConvolution(shared_ptr<NodeList> actList, shared_ptr<BaseFeFunction> feFct);
    CoefFunctionConvolution(Double tau, shared_ptr<ElemList> actList, shared_ptr<FeSpace> ptFeSpace, BaseBOperator * bOp );

    //! Destructor
    virtual ~CoefFunctionConvolution();

    virtual string GetName() const { return "CoefFunctionConvolution"; }


    //! evaluate the recursive convolution
    void RecursiveConvolution(Double dt);

    //! Return real-valued vector at integration point
    void GetVector(Vector<Double>& vec, const LocPointMapped& lpm );

    //! return vector size
    UInt GetVecSize() const;

  protected:

    void Init();

    // return the id in the u_ or histvals_ vector
    UInt GetId(UInt nodeId) ;

    void ComputeRecursion();

    Vector<Double> GetSavedVals(const StdVector<UInt> nids, const Vector<Double> * vals);
    void ComputeCurrentVals(Vector<Double>& res, const Vector<Double> * h0, const Vector<Double> * u0, const Vector<Double> * u);

    //! vector size
    UInt vecSize_;

    shared_ptr<ElemList> entityList_;
    //NodeList nodeList_;
    shared_ptr<FeSpace> ptFeSpace_;
    shared_ptr<IntScheme> intScheme_;

    shared_ptr<BaseFeFunction> ptFeFct_;// nedded more than once (ptFeSpace_->GetFeFunction());

    //! history variables map[element][integration point]
    typedef std::map<UInt, StdVector<UInt> > IpMap;
    typedef std::map<UInt, IpMap > ElMap;
    ElMap ipHist_;
    std::map<UInt, UInt> indexMap_;
    Vector<Double> histvals_; // vector of memory variables (of last increment)
    Vector<Double> u_; // displacement vector (of last increment)
    BaseBOperator * bOp_;

    //! Pointer to math parser instance
    MathParser* mp_;
    //! Handle for expression
    unsigned int mHandleStep_;
    unsigned int mHandleTime_;
    unsigned int mHandleDt_;
    //Double t_; // current time
    Double tau_; // relaxation time
    Double lam_; // ratio of dt/tau
    Double f1_; // factor for first term: exp(-dt/tau)
    Double f2_; // factor for second term: tau*(1-f1)/dt
    //shared_ptr<BaseFeFunction> feFct_;(ptFeSpace_->GetFeFunction());
  };

}
#endif
