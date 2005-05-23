#ifndef FILE_LINPIEZOINT
#define FILE_LINPIEZOINT


#include "Elements/basefe.hh"
#include "Forms/bdbInt.hh"
#include "DataInOut/MaterialData.hh"

namespace CoupledField
{

  //! base class for linear piezoelectric simulations

  //! The main objective of this class is to implement the pure virtual
  //! method calcBMat from BDBInt for the case of a linear piezoelectric
  //! simulation. In this case the BDB operator looks as follows
  //! \f[
  //! [BDB] = \left(\begin{array}{cc} B_a & 0 \\ 0 & \tilde{B}_a \end{array}
  //! \right)^T
  //! \left(\begin{array}{cc} c^E & e^T \\ e & -\varepsilon^s \end{array}
  //! \right)
  //! \left(\begin{array}{cc} B_a & 0 \\ 0 & \tilde{B}_a \end{array} \right)
  //! \enspace.
  //! \f]
  //! Here \f$B_a\f$ is derived from the mechanical part, while
  //! \f$\tilde{B}_a\f$ comes from the electrical part,
  //! \f$c^E\f$ is the tensor of mechanical modulus,
  //! \f$\varepsilon^s\f$ is the tensor of electric permittivity and
  //! \f$e\f$ is the tensor of piezoelectric coupling.
  class linPiezoInt : public BDBInt
  {
  public:

    //!
    linPiezoInt(): BDBInt(), actOrientation(yz),isDamping_(FALSE), factorDamp_(1.0) {;}
    

    //! Constructor
    linPiezoInt(BaseFE * aptelem, MaterialData & matData) 
      : BDBInt(aptelem, matData), actOrientation(yz)
    {
      ENTER_FCN( "linPiezoInt::linPiezoInt" );
      isDamping_  = FALSE;
      factorDamp_ = 1.0;
    }
  
    //! Constructor
    linPiezoInt(MaterialData & matData)
      : BDBInt(matData), actOrientation(yz)
    {
      ENTER_FCN( "linPiezoInt::linPiezoInt" );
      isDamping_  = FALSE;
      factorDamp_ = 1.0;
    }
  
    //! Destructor
    ~linPiezoInt()
    {
      ENTER_FCN( "linPiezoInt::~linPiezoInt" );
    };

    /// calculates mechanical stresses (vector notation), overloaded for real and complexvalued stresses
    virtual void CalcStressVec(Vector<Double>& StressVec, Integer ip, Matrix<Double> & ptCoord);  
    virtual void CalcStressVec(Vector<Double>& StressVec, Integer ip, Matrix<Double> & ptCoord,
                               Integer comp, Double Dval); 
    virtual void CalcStressVec(Vector<Complex>& StressVec, Integer ip, Matrix<Double> & ptCoord);  

    //     virtual void SetActElemSol(Matrix<Double>& disp) {
    //       ENTER_FCN( "linPiezoInt::SetActElemSol 1" );
    //       elemSol_ = disp;};

    //     virtual void SetActElemSol(CFSMatrix& disp) {
    //         ENTER_FCN( "linPiezoInt::SetActElemSol 2" );
    //         Matrix<Double> & helpSol = dynamic_cast <Matrix<Double>&> (disp);
    //         elemSol_ = helpSol;};

    virtual void SetActElemSol(CFSMatrix& disp){
      ENTER_FCN("linPiezoInt::SetActElemSol 3");
      
      if(disp.IsComplex())
        //        Matrix<Complex> & elemSolComplex_ = dynamic_cast <Matrix<Complex>&> (disp);
        elemSol_ = new Matrix<Complex> (dynamic_cast<Matrix<Complex>&> (disp));
      else
        elemSol_ = new Matrix<Double> (dynamic_cast<Matrix<Double>&>(disp));
      //     Matrix<Double> & elemSol_ = dynamic_cast <Matrix<Double>&> (disp);


    };

    //! set multiplicative factor for matrix
    virtual void SetFactor(Double factor) 
    { if (factor <= 0) {
      Error("Additional damping factor cannot be zero");
    }
    factorDamp_ = factor;
    };

  protected:

    //! returns B - matrix for BDB

    //! The method computes the matrix B defined as
    //! \f$
    //! B = \left(\begin{array}{cc} B_a & 0 \\ 0 & \tilde{B}_a \end{array}
    //! \right)
    //! \f$.
    //! For a 3D simulation e.g. we have
    //! \f[ B_a = \left( \begin{array}{*{6}{c}}
    //! \frac{\partial N_a}{\partial x} & 0 & 0 & 0 &
    //! \frac{\partial N_a}{\partial z} & \frac{\partial N_a}{\partial y} \\
    //! 0 & \frac{\partial N_a}{\partial y} & 0 &
    //! \frac{\partial N_a}{\partial z} & 0 & \frac{\partial N_a}
    //! {\partial x}\\
    //! 0 & 0 & \frac{\partial N_a}{\partial z} &
    //! \frac{\partial N_a}{\partial y} & \frac{\partial N_a}{\partial x}
    //! & 0\\
    //! \end{array}\right)^T
    //! \enspace,\quad
    //! \tilde{B }_a = \left(\begin{array}{c}
    //! \frac{\partial N_a}{\partial x} \\[1ex]
    //! \frac{\partial N_a}{\partial y} \\[1ex]
    //! \frac{\partial N_a}{\partial z}
    //! \end{array}\right) \f]
    //! where \f$N_a\f$ are the Finite Element ansatz functions.
    //! To be more precise the above matrix B is computed for the given
    //! integration point ip and every FE ansatz function belonging to a node
    //! of the element. These partial matrices are appended one after another
    //! in a row-wise fashion to form the return matrix bMat.
    void calcBMat(Matrix<Double> & bMat, Integer ip,
                  Matrix<Double> & ptCoord);


    /// calculates the material data for the axisymmetric case
    void CalcAxiMaterialMat(Matrix<Double> & dMat);
  
    /// calculates the material data for the axisymmetric case
    void CalcPlaneStrainMaterialMat(Matrix<Double> & dMat);

    /// calculates the material data for the 3d case
    void Calc3DMaterialMat(Matrix<Double> & dMat);

   

    //! If set to true, stiffnessmatrix is computed for damping
    //! (just mechanical part)
    Boolean isDamping_;

    //!
    orientation2D actOrientation;

    //! displacement of all nodes of actual element
    CFSMatrix * elemSol_;
    //    Matrix<Double> elemSol_;
    //    Matrix<Complex> elemSolComplex_;

    Double factorDamp_;
  };


  //! class for piezoelectric simulations in 3D
  class linPiezo3DInt : public linPiezoInt
  {
  public:

    //! Constructor
    linPiezo3DInt(BaseFE * aptelem, MaterialData & matDat);
  
    //! Constructor
    linPiezo3DInt( MaterialData & matData, Boolean isdamping = FALSE ) :
      linPiezoInt(matData)
    {
      ENTER_FCN( "linPiezo3DInt::linPiezo3DInt" );
      isDamping_ = isdamping;
    }

    //! Destructor
    ~linPiezo3DInt();

   
  protected:

    //! calculate the data-matrix D.

    //! The method computes the matrix D of material properties.
    //! The latter is given by
    //! \f[
    //! D = \left(\begin{array}{cc} c^E & e^T \\ e & \varepsilon^s
    //! \end{array}\right)
    //! \f]
    //! where \f$c^E\f$ is the tensor of mechanical modulus,
    //! \f$\varepsilon^s\f$ is the tensor of electric permittivity and
    //! \f$e\f$ is the tensor of piezoelectric coupling.
    void calcDMat(Matrix<Double> & dMat);

    //! Returns dimension of data-matrix D.

    //! The method returns the dimension of the data-matrix D. In a 3D
    //! piezoelectric simulation D is a square 9x9 matrix, so the return
    //! value is 9.
    virtual Integer getDimD(){return 9;};
  
    //! Returns number of degrees of freedom per node.

    //! Returns the number of degrees of freedom per individual FE node.
    //! The current return value is fixed to 4, since in 3D simulations
    //! we have three mechanical and one potential component.
    virtual Integer getNrDofs(){ return 4; };

    //! Damps DMat with the factor \f$(1+j*\omega*\beta)\f$
    void calcDMaterialMatWithComplexDamping(Matrix<Complex> &dMat, Double &beta, Double &omega);

  };


  //! base class for linear piezoelectric simulations

  //! The main objective of this class is to implement the pure virtual
  //! method calcBMat from BDBInt for the case of a linear piezoelectric
  //! simulation. In this case the BDB operator looks as follows
  //! \f[
  //! [BDB] = \left(\begin{array}{cc} B_a & 0 \\ 0 & \tilde{B}_a \end{array}
  //! \right)^T
  //! \left(\begin{array}{cc} c^E & e^t \\ e & \varepsilon^s \end{array}
  //! \right)
  //! \left(\begin{array}{cc} B_a & 0 \\ 0 & \tilde{B}_a \end{array} \right)
  //! \enspace.
  //! \f]
  //! Here \f$B_a\f$ is derived from the mechanical part, while
  //! \f$\tilde{B}_a\f$ comes from the electrical part,
  //! \f$c^E\f$ is the tensor of mechanical modulus,
  //! \f$\varepsilon^s\f$ is the tensor of electric permittivity and
  //! \f$e\f$ is the tensor of piezoelectric coupling.
  
  class piezoAxiInt : public linPiezoInt
  {
  public:

    //! Constructor
    piezoAxiInt(BaseFE * aptelem, MaterialData & matData)
      : linPiezoInt(aptelem, matData)
    {
      ENTER_FCN( "piezoAxiInt::piezoAxiInt" );
      isaxi_ = TRUE;
    }

    //! Constructor
    piezoAxiInt(MaterialData & matData, Boolean isdamping=FALSE)
      : linPiezoInt(matData)
    {
      ENTER_FCN( "piezoAxiInt::piezoAxiInt" );
      isDamping_ = isdamping;
      isaxi_ = TRUE;
    }

    //! Destructor
    ~piezoAxiInt()
    {
      ENTER_FCN( "piezoAxiInt::~piezoAxiInt" );
    };

       
  protected:

    //if set to true, stiffnessmatrix is computed for damping (just mechanical part)
    //! Returns B - matrix for BDB
    //! The method computes the matrix B defined as
    //! \f$
    //! B = \left(\begin{array}{*{2}{c}} B_a & 0 \\ 0 & \tilde{B}_a \end{array}
    //! \right)
    //! \f$. 
    //! For an axisymmetric simulation we have
    //! \f[ B_a = \left( \begin{array}{*{2}{c}}
    //! \frac{\partial N_a}{\partial r} & 0 \\
    //! 0 & \frac{\partial N_a}{\partial z} \\
    //! \frac{\partial N_a}{\partial z} & \frac{\partial N_a}{\partial r}\\
    //! \frac{1}{r} & 0 \\
    //! \end{array}\right)
    //! \enspace,\quad
    //! \tilde{Bš}_a = \left(\begin{array}{c}
    //! \frac{\partial N_a}{\partial r} \\
    //! \frac{\partial N_a}{\partial z} \\
    //! \end{array}\right) \f]
    //! where \f$N_a\f$ are the Finite Element ansatz functions.
    //! To be more precise the above matrix B is computed for the given
    //! integration point ip and every FE ansatz function belonging to a node
    //! of the element. These partial matrices are appended one after another
    //! in a row-wise fashion to form the return matrix bMat.
    
    //!
    virtual void calcBMat(Matrix<Double> & bMat, Integer ip,
                          Matrix<Double> & ptCoord);
        
    //!   
    void  calcDMat(Matrix<Double> & dMat);

    //!
    virtual Integer getDimD(){return 6;};

    //!
    virtual Integer getNrDofs(){ return 3; };
    
    //! Damps DMat with the factor \f$(1+j*\omega*\beta)\f$
    void calcDMaterialMatWithComplexDamping(Matrix<Complex> & dMat, Double & beta, Double & omega);

  };


  //! class for  plane strain case
  class piezoPlainStrainInt : public linPiezoInt
  {
  public:

    //! Constructor
    piezoPlainStrainInt(BaseFE * aptelem, MaterialData & matData)
      : linPiezoInt(aptelem, matData)
    {
      ENTER_FCN( "piezoPlainStrainInt::piezoPlainStrainInt" );
    }
  
    //! Constructor
    piezoPlainStrainInt(MaterialData & matData, Boolean isdamping=FALSE)
      : linPiezoInt(matData)
    {
      ENTER_FCN( "piezoPlainStrainInt::piezoPlainStrainInt" );
      isDamping_=isdamping;
    }
  
    //! Destructor
    ~piezoPlainStrainInt()
    {
      ENTER_FCN( "piezoPlainStrainInt::~piezoPlainStrainInt" );
    };

  protected:
  
    //! returns B - matrix for BDB
    //! The method computes the matrix B defined as
    //! \f$
    //! B = \left(\begin{array}{cc} B_a & 0 \\ 0 & \tilde{B}_a \end{array}
    //! \right)
    //! \f$.
    //! For a 2D simulation e.g. we have
    //! \f[
    //! B_a = \left( \begin{array}{*{2}{c}}
    //! \displaystyle\frac{\partial N_a}{\partial x} & \displaystyle 0 \\[2ex]
    //! \displaystyle 0 & \displaystyle\frac{\partial N_a}{\partial y} \\[2ex]
    //! \displaystyle\frac{\partial N_a}{\partial y} &
    //! \displaystyle\frac{\partial N_a}{\partial x}
    //! \end{array}\right)  
    //! \enspace,\quad
    //! \tilde{B }_a = \left(\begin{array}{c}
    //! \displaystyle\frac{\partial N_a}{\partial x} \\[2ex]
    //! \displaystyle\frac{\partial N_a}{\partial y} 
    //! \end{array}\right) \f]  
    //! where \f$N_a\f$ are the Finite Element ansatz functions.
    //! To be more precise the above matrix B is computed for the given
    //! integration point ip and every FE ansatz function belonging to a node
    //! of the element. These partial matrices are appended one after another
    //! in a row-wise fashion to form the return matrix bMat.
    
    //!
    void calcBMat(Matrix<Double> & bMat, Integer ip,
                  Matrix<Double> & ptCoord);
 
    //!
    void calcDMat(Matrix<Double> & dMat);

    //! Returns dimension of data-matrix D.
    //! The method returns the dimension of the data-matrix D. In a 2D
    //! piezoelectric simulation D is a square 5x5 matrix, so the return
    //! value is 5.
    virtual Integer getDimD(){return 5;};

    
    //! Returns number of degrees of freedom per node.
    //! Returns the number of degrees of freedom per individual FE node.
    //! The current return value is fixed to 3, since in 2D simulations
    //! we have two mechanical and one potential component.
    virtual Integer getNrDofs(){ return 3; };  

    //! Damps DMat with the factor \f$(1+j*\omega*\beta)\f$
    void calcDMaterialMatWithComplexDamping(Matrix<Complex> &dMat, Double &beta, Double &omega);

  };


}

#endif
