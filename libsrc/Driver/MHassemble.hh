#ifndef FILE_MHASSEMBLE
#define FILE_MHASSEMBLE

#include "assemble.hh"

#include "Utils/StdVector.hh"
#include "Forms/baseForm.hh" 
#include "DataInOut/timefunc.hh"
#include "Utils/nodestoresol.hh"
#include "Utils/elemstoresol.hh"
#include "PDE/nodeEQN.hh"
#include "PDE/SinglePDE.hh"

#include "olas.hh"

namespace CoupledField
{

  // forward class declaration
  class MHMaterialData;

 class MHassemble : public Assemble
  {
  public:
    MHassemble(BaseSystem * algsys, Grid * agrid, 
               const std::string bcSequenceTag );
    
    virtual ~MHassemble();


    //!
    virtual void SetFrequency(Double actFreq);
   

    //! transform element vector to account for harmonic analysis
    virtual void TransformVector2Harmonic(Vector<Double>& harmMat, 
                                          Vector<Double> origVec,
                                          const Double valPhase);

    // void MHassembleMatrices(StdVector<RegionIdType> subdoms_,
//                             StdVector< StdVector<IntegratorDescriptor *>* > integrators_,
//                             NodeEQN * ptEQN1_,
//                             NodeEQN * ptEQN2_);

     virtual void AssembleMatrices();
    
    StdVector<RegionIdType> MHsubdoms_;  //!< subdomain-levels belongig to PDE
    //    SinglePDE * ptMyPDE_;

    //! transform element matrix to account for harmonic analysis
    virtual void TransformMatrix2Harmonic(Vector<Double>& harmMat,
                                          Matrix<Double> origMat,
                                          const FEMatrixType matrixType,
                                          const piezoMaterialType 
                                          piezoMatType);

    /// adds integrators to the pde
    virtual void AddIntegrator(IntegratorDescriptor * intDescr, 
                               const RegionIdType subdomain);

    /// adds surface integrators to the pde
    virtual void AddSurfIntegrator(IntegratorDescriptor * intDescr, 
				   const RegionIdType subdomain);

    //! define discrete PDE
    virtual void MatrixSettings(){};

    //! set information for algebraic system about PDE. set matrix factors
    virtual void SetMatrixFactors(){};

    // void SetPtr2EQNData(NodeEQN * aPtNodeEQN1,
    //                  NodeEQN * aPtNodeEQN1 = NULL );



  private: 
    BaseSystem * algsys_;                //!< pointer to algebraic system  
    //Grid * ptgrid_;                      //!< pointer to Grid

   //  NodeEQN * ptEQN1_;                    //!< pointer to equation data of pde1
//     NodeEQN * ptEQN2_;                    //!< pointer to equation data of pde2

    Integer nrMultHarms_;

//     Boolean systemMatrix_;               //!< need system matrix (TRUE/FALSE)
//     Boolean stiffnessMatrix_;            //!< need stiffness matrix (TRUE/FALSE)
//     Boolean massMatrix_;                 //!< need mass matrix (TRUE/FALSE)
//     Boolean convectionMatrix_;           //!< need convective matrix (TRUE/FALSE)

    //  StdVector< StdVector<IntegratorDescriptor *>* > integrators_;

    /// vector of all needed integrators (every subdomain needs one "list of integrators")
    //    StdVector< StdVector<BaseIntDescriptor *>* > rhsIntegrators_;


    /// vector of all needed RHS src-intergators (not every subdomain needs a "list of rhs_source_integrators")
    //    StdVector< StdVector<BaseIntDescriptor *>* > rhsSrcIntegrators_;

    //! needed to identify the PDE uniquely in the algebraic system
        PdeIdType MHpdeId1_;

    //! identifier for the second PDE

    //! needed to identify the PDE uniquely in the algebraic system
        PdeIdType MHpdeId2_;

    MHMaterialData *ptMHMat_;





 };

}
#endif

