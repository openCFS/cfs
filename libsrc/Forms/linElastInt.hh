#ifndef FILE_LINELASTINT
#define FILE_LINELASTINT

#include <Elements/basefe.hh>
#include <Forms/bdbInt.hh>
#include <DataInOut/MaterialData.hh>
#include <General/environment.hh>

namespace CoupledField {
  
  
  //! base class for calculation of linear elasticity
  class linElastInt : public BDBInt {

  public:

    //! Constructor
    linElastInt( BaseFE *aptelem, MaterialData &matData );

    //! Constructor
    linElastInt( MaterialData &matData );
  
    //! Destructor
    virtual ~linElastInt();

  protected:    

    //! calculates the material data for the axisymmetric case
    void CalcAxiMaterialMat( Matrix<Double> &dMat,
                             enum orientation2D actOrientation );

    //! calculates the material data for the axisymmetric case
    void CalcPlaneStrainMaterialMat( Matrix<Double> &dMat,
                                     enum orientation2D actOrientation );

    //! calculates the material data for the axisymmetric case
    void Calc3DMaterialMat( Matrix<Double> &dMat );

    //! returns B - matrix for BDB
    virtual void calcBMat( Matrix<Double> &bMat, UInt ip,
                           Matrix<Double> &ptCoord );

    //! Orientation for 2D simulations like axi or plane strain

    //! This represents the orientation for 2D simulations like in the
    //! axi-symmetric or plane strain case. Orientation actually is a
    //! misnomer since what we need is to define the 2D plane of computation
    //! and (for the axi-symmetric case) also the axis of symmetry.

    //! \note There are no setter methods for this property and it is
    //!       hardcoded to yz for technical reasons until further notice.
    orientation2D actOrientation;

  };
  

  //! class for calculation of mechanical plain strain state
  class mechPlainStrainInt : public linElastInt {  

  public:

    //! Constructor
    mechPlainStrainInt( BaseFE *aptelem, MaterialData &matDat );

    //! Constructor
    mechPlainStrainInt( MaterialData &matDat );
  
    //! Destructor
    virtual ~mechPlainStrainInt();
  
  protected:
  
    //! calculate the data-matrix for 2D plain-strain
    virtual void calcDMat(Matrix<Double> & dMat);

    //! returns dimension of D matrix
    virtual UInt getDimD(){
      return 3;
    };

    //! returns nr. of degrees of freedom
    virtual UInt getNrDofs(){
      return 2;
    };
  };


  //! class for calculation of mechanical axisymmetric state
  class mechAxiInt : public linElastInt {  

  public:

    //! Constructor
    mechAxiInt(BaseFE * aptelem, MaterialData & matDat);

    //! Constructor
    mechAxiInt(MaterialData & matDat);

    //! Destructor
    virtual ~mechAxiInt();
  
  protected:
  
    //! calculate the data-matrix for 2D axi
    virtual void calcDMat(Matrix<Double> & dMat);

    //! returns dimension of D matrix
    virtual UInt getDimD(){
      return 4;
    };
  
    //! returns nr. of degrees of freedom
    virtual UInt getNrDofs(){
      return 2;
    };
  };


  //! class for calculation of mechanical plain strain state
  class mech3DInt : public linElastInt {  

  public:

    //! Constructor
    mech3DInt(BaseFE * aptelem, MaterialData & matDat);

    //! Constructor
    mech3DInt(MaterialData & matDat);
  
    //! Destructor
    virtual ~mech3DInt();
  
  protected:
  
    //! returns D - matrix for BDB
    virtual void calcDMat(Matrix<Double> & dMat);

    //! returns dimension of D matrix
    virtual UInt getDimD(){
      return 6;
    };

    //! returns nr. of degrees of freedom
    virtual UInt getNrDofs(){
      return 3;
    };
  };


} //end namespace

#endif // FILE_LINELASTINT
