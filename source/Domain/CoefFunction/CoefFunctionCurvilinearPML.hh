// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     CoefFunctionCurvilinearPML.hh
 *       \brief    Implements a custom CoefFunction for a curvilinear PML that is supposed to 
 *                 function in arbitrarily-shaped, convex domains. 
 *                 The CoefFunctionCurvilinearPML requires geometric properties (e.g, normal 
 *                 vectors, principal-direction vectors and principal curvatures) that must be 
 *                 provided or can be computed via the automaticLayerGeneration property.
 *       \date     Mar 2023
 *       \author   pheidegger
 */
//================================================================================================

#ifndef COEFFUNCTIONCURVILINEARPML_HH_
#define COEFFUNCTIONCURVILINEARPML_HH_

#include "CoefFunctionPML.hh"
#include "Domain/Mesh/Grid.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "FeBasis/H1/H1Elems.hh"

namespace CoupledField{
  //================================================================================================
  // Curvilinear PML
  //================================================================================================
  /* A curvilinear PML that is supposed to function in arbitrarily-shaped, convex 
  * domains. 
  * The CoefFunctionCurvilinearPML currently works in combination with an automatic layer generation
  * that triggers the computation of the bounary geometry (normal vectors, principal-direction 
  * vectors and principal curvatures) in the grid class. 
  * The computed nodal data is interpolated to the integration points and the scalars/tensors
  * are assembled directly there.
  */

  //! the type of output the coefFunction is computing (currently only for curvilinearPML):
  //! TENSOR (tensor output) ... the whole factorized PML matrix J^-1 = A^T (I-D)^-1 A
  //! DETERMINANT (scalar output) ... the determinant of the PML matrix
  //! DAMP_FACTOR (vector output) ... the eigenvalues of the PML matrix (entries of the 'inner' diagonal matrix outputted as vector)
  //! DISTANCE (scalar output) ... the distance of a point to its closest point on the PML interface
  enum OutputType {TENSOR, DETERMINANT, DAMP_FACTOR, DISTANCE};

  template<typename T>
  class CoefFunctionCurvilinearPML : public CoefFunctionPMLBase<T> {
  public:
    //! constructor
    //! \param pmlDef (in) param node that holds the PML parameters defined in the XML
    //! \param matCoef (in) coef function with the speed of wave (gets assigned to member matCoef_ in the base class)
    //! \param EntList (in) entity list containing all entities of the assigned PML region
    //! \param pdeDomains (in) StdVector holding all domains of the current PDE
    //! \param outputType (in) enum to set which output should be computed
    CoefFunctionCurvilinearPML(PtrParamNode pmlDef, PtrCoefFct matCoef, shared_ptr<EntityList> EntList,
                          StdVector<RegionIdType> pdeDomains, OutputType outputType, CoefFunction::CoefDimType dimType);

    //! destructor
    virtual ~CoefFunctionCurvilinearPML();

    //! Return the complex-valued tensor of the PML at integration point. The Tensor will be used in 
    //! TensorScaledGradientOperator() to compute the scaled derivative (for AcousticPDE).
    //! This function gets only called when the outputType_ is TENSOR
    //! \param tensor (out) the tensor holding the complex-valued inverse Jacobi matrix of
    //!                     the curvilinear coordinate stretch
    //! \param lpm (in) the local point for which the tensor is returned
    virtual void GetTensor(Matrix<Complex>& tensor, const LocPointMapped& lpm ) override; 

    //! Return a real-valued tensor at integration point. This is not implemented as the 
    //! coordinate stretching operation is generally complex-valued.
    virtual void GetTensor(Matrix<Double>& tensor, const LocPointMapped& lpm ) override {
      EXCEPTION("CoefFunctionCurvilinearPML::GetTensor() is only implemented for Complex values.");
    }

    //! Return complex-valued vector at integration point. Returns one of the four vector-valued outputTypes, 
    //! depending on which type is set.
    //! This function gets only called when the outputType_ is DAMP_FACTOR
    virtual void GetVector(Vector<Complex>& vec, const LocPointMapped& lpm ) override;

    //! Return scalar at integration point (lpm). Returns the determinant of the PML Jacobi matrix,
    //! This function gets only called when the outputType_ is DETERMINANT
    //! \param val (out) the scalar holding the complex-valued Jacobi determinant of the 
    //!                  curvilinear coordinate stretch
    //! \param lpm (in) the local point for which the scalar is returned
    virtual void GetScalar(Complex& val, const LocPointMapped& lpm ) override;

    //! Return scalar at integration point (lpm). Returns the distance of the lpm to the 
    //! PML interface.
    //! This function gets only called when the outputType_ is DISTANCE
    //! \param val (out) the distance of the lpm to the interface
    //! \param lpm (in) the local point for which the scalar is returned
    virtual void GetScalar(Double& val, const LocPointMapped& lpm ) override;

    //! use the base functions
    using CoefFunctionPMLBase<T>::UpdateOmega;
    using CoefFunctionPMLBase<T>::CreateDampFunction;

  private:
    //! read data from the paramNode and store
    void ReadDataPML(PtrParamNode pmlDef);

    //! check if autoLayerGeneration is specified and if it fits to the current PML region.
    //! If true, sets the layerGenNode_, numLayers_, elemHeight_, numSurfNodes_, and layerThickness_.
    void CheckForLayerGenerationNode(PtrParamNode pmlDef);

    //! check if invalid parameters are set and warn the user
    void CheckForInvalidParams(PtrParamNode pmlDef);

    //! compute distances of every node in the layer to its closest point on the PML interface.
    //! Stores the values in thicknessOnNodes_.
    void GetThicknessOnNodes();

    //! returns the idx position of a given nodeId. Exploits the knowledge of the auto layer 
    //! generation params to speed up the search.
    //! \param nodeId (in) the node id for which the index needs to be searched
    //! \param guess (in) a guess of the index, for the find algorithm when searching unstructured nodes
    //! \return the index of the desired node id within the PML layer
    UInt GetIdxByNodeId(const UInt& nodeId, unsigned int guess) const;

    //! Create identity operators for mapping the nodal values (e.g. geometry data) to local points.
    //! Defines one operator for vectors: vectorMappingOperator_; and one for scalars: scalarMappingOperator_;
    void CreateMappingOperators();

    //! Set the dimType_ of the coefFunction according to the assigned outputType_
    void SetDimType();

    //! member to set which quantity should be outputted/computed by the coefFunction
    OutputType outputType_;

    //! pointer to the current grid
    Grid* grid_;

    //! pointer to the layer generation node that keeps auto-mesh-generation parameters 
    //! (elemHeight and numLayers). Is empty if not specified.
    PtrParamNode layerGenNode_;

    //! the volume region to act on
    RegionIdType volRegion_;

    //! all iso surface regions in the volume (the 0th entry should be the PML interface)
    StdVector<RegionIdType> surfRegions_;

    //! the PML interface
    RegionIdType ifSurfRegion_;

    //! the lowest and highest node Id generated by the automatic layer generation
    StdVector<UInt> minIds_;
    StdVector<UInt> maxIds_;

    //! pointer to the stored geometry of every node in the PML layer
    shared_ptr<Grid::NodeGeometry> ptrNodeGeom_;

    //! vector containing the relative distance of every node to its closest 
    //! point on the PML interface
    StdVector<Double> thicknessOnNodes_;

    //! BOperator to map solutions to arbitrary points. Used is an identity operator,
    //! defined for interpolating Double-valued vectors and scalars. Tensors must be 
    //! stacked into one dimension 
    shared_ptr<BaseBOperator > tensorMappingOperator_;
    shared_ptr<BaseBOperator > vectorMappingOperator_;
    shared_ptr<BaseBOperator > scalarMappingOperator_;
  
    //! Layer generation parameters
    UInt numLayers_;      //number of iteratively generated layers
    Double elemHeight_;     //height of the elements in the layer
    UInt numSurfNodes_;   //number of nodes on the interface
    Double layerThickness_; //total thickness of the PML layer
  };
}
#endif /* COEFFUNCTIONCURVILINEARPML_HH_ */