// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_SPACEERROR_2002
#define FILE_SPACEERROR_2002

#include "Utils/tools.hh"

#include "PDE/basePDE.hh"

#ifdef ADAPTGRID
#include "Vertex.h"
#include "Edge.h"
#include "Element.h"
#endif

namespace CoupledField
{

  //! class for space error estimation
  class SpaceErrorEstimator
  {
  public:
    //! constructor 
    SpaceErrorEstimator();

    //! Deconstructor
    virtual ~SpaceErrorEstimator();

    //! Init pointer to PDE
    /*!
      \param aptPDE pointer to the class of the PDE
    */
    virtual void Init(BasePDE * aptPDE){ ptPDE_=aptPDE;}

    //! calculation of error map
    /*!
      \param sol solution 
      \param subdoms vector with names of subdomains, on which we do calculation
      \param ptgrid pointer to the Grid
      \param errorMap out: vector with calculated relative error for each element.
      \f$\|\mbox{grad}_{\mbox{SPR}} - \mbox{grad}_{\mbox{FEM}}\|_{L2}\f$
      \param  atotalErr out: total error.\f$\|\mbox{grad}_{\mbox{SPR}}\|_{L2}\f$
    */
    void CalcErrorMap(const Vector<Double> & sol, std::vector<std::string> & subdoms,
                      Grid * ptgrid, Vector<Double> & relErrorMap, Double & atotalErr, const Integer level);

    //! calculation of error map for harmonic analysis of acoustic equation
    /*!
      \param sol solution 
      \param subdoms vector with names of subdomains, on which we do calculation
      \param ptgrid pointer to the Grid
      \param errorMap out: vector with calculated relative error for each element.
      \f$\|\mbox{grad}_{\mbox{SPR}} - \mbox{grad}_{\mbox{FEM}}\|_{L2}\f$
      \param  atotalErr out: total error.\f$\|\mbox{grad}_{\mbox{SPR}}\|_{L2}\f$
    */
    void CalcErrorMapHarmonic(const Vector<Double> & solRe, const Vector<Double> & solIm, 
                              std::vector<std::string> & subdoms, Grid * ptgrid,
                              Vector<Double> & relErrorMap, Double & atotalErr,
                              const Integer level);

    //! recovery procedure for the elements-patch
    /*!
      \param Elems vector with information about elements
      \param ptgrid pointer to the Grid
      \param sol solution 
      \param aComponent number of component of gradient 
      \param result (output) result
      \param locations (output) global numbers of nodes in patch. needful for the vector with result
    */
    void RecoveryProcedure4ElemsPatch(const std::vector<Elem*> &Elems,
                                      Grid * ptgrid, const Vector<Double> & sol,
                                      const Integer aComponent,
                                      Vector<Double> &result,
                                      std::vector<Integer>&locations,
                                      const Integer level=0);

    //! kelly error estimation
    /*!
      \param sbdoms vector with names of subdomains, on which results are calculated
      \param level level of the Grid
    */
    //  void KellyError(Grid * ptgrid, std::vector<std::string> &sbdoms,  Vector<Double> & errorMap, const Integer level);

  protected:
    //!
    BasePDE * ptPDE_;  

    //!
    Grid * ptGrid_;  

    //! calculation error for one element
    /*!
      \param elem pointer to the element
      \param SPRgrad (output) returned calculated SPR gradient 
      \param error (output) retuned value of the error
      \param normGradSPR (output) L2-norm of SPR gradient
      \param sol (input) solution
      \param ptgrid pointer to Grid
      \param level level in grid hierarchy
    */
    void CalcErrorForElem(const Elem* elem, const Vector<Double>* SPRgrad,
                          Double & error, Double & normGradSPR,
                          const Vector<Double> & sol, Grid * ptgrid,
                          const Integer level);

    //! calculation error for one element
    /*!
      \param elem pointer to the element
      \param SPRgrad (output) returned calculated SPR gradient 
      \param error (output) retuned value of the error
      \param normGradSPR (output) L2-norm of SPR gradient
      \param sol (input) solution
      \param ptgrid pointer to Grid
      \param level level in grid hierarchy
    */
    void CalcErrorForElemHarmonic(const Elem* elem, const Vector<Double>* SPRgrad,
                                  const Vector<Double>* SPRgradIm,
                                  Double & error, Double & normGradSPR,
                                  const Vector<Double> & sol,
                                  const Vector<Double> & solIm,
                                  Grid * ptgrid,
                                  const Integer level);

  private:  

    //! auxiliary function for calculation of Kelly error
    //  Double IntegralOverRegularFace_KellyError(const Elem * ptElem, Elem ** ptNeighborsOfFace);
  
    //! auxialary procedure for calculating number of nodes in patch
    /* direct calculating is implemented */
    Integer CalcNumberOfNodesInPatch(const std::vector<Elem*> & patch);

    //   //! form list of faces for the elem
    //   void FormFacesForElem(const Elem * ptelem, std::vector<Elem*> & afaces);

    //   //! find neighbours elems for elem
    //   void FormNeighForFace(Elem ** neighFcs, Grid * grid, Elem * face);
  };

  inline SpaceErrorEstimator::~SpaceErrorEstimator()
  {
    ENTER_FCN( "SpaceErrorEstimator::~SpaceErrorEstimator" );
    ;
  }

} // end of namespace

#endif
