#ifndef FILE_FLUIDMECHSTIFFINT
#define FILE_FLUIDMECHSTIFFINT

#include "fluidMechInt.hh"

namespace CoupledField
{


//**************************************************************************************************************************
//***************** Stabilized FEM ****************************************************************************************
//**************************************************************************************************************************
class FluidMechPlaneStiffInt_UV : public FluidMechInt
{
public:	

  FluidMechPlaneStiffInt_UV(Double density, Double kinematicViscosity,Matrix<Double> & stabilParams);

  virtual ~FluidMechPlaneStiffInt_UV();
  
  void CalcElementMatrix( Matrix<Double>& elemMat, EntityIterator& ent1, EntityIterator& ent2 );

protected:    
  Matrix<Double> & stabilParams_;
};


//**************************************************************************************************************************
class FluidMechPlaneStiffInt_PQ : public FluidMechInt
{
public:

  FluidMechPlaneStiffInt_PQ(Double density, Double kinematicViscosity, Matrix<Double> & stabilParams);

  virtual ~FluidMechPlaneStiffInt_PQ();

  void CalcElementMatrix( Matrix<Double>& elemMat, EntityIterator& ent1, EntityIterator& ent2 );

protected:    
  Matrix<Double> & stabilParams_;
};


//**************************************************************************************************************************
class FluidMechPlaneStiffInt_UQ : public FluidMechInt
{
public:

  FluidMechPlaneStiffInt_UQ(Double density, Double kinematicViscosity, Matrix<Double> & stabilParams);
  
  virtual ~FluidMechPlaneStiffInt_UQ();

  void CalcElementMatrix( Matrix<Double>& elemMat, EntityIterator& ent1, EntityIterator& ent2 );

protected:    
  Matrix<Double> & stabilParams_;
};

//**************************************************************************************************************************
class FluidMechPlaneStiffInt_PV : public FluidMechInt
{
public:
  FluidMechPlaneStiffInt_PV(Double density, Double kinematicViscosity, Matrix<Double> & stabilParams);
  virtual ~FluidMechPlaneStiffInt_PV();
  void CalcElementMatrix( Matrix<Double>& elemMat, EntityIterator& ent1, EntityIterator& ent2 );
protected:    
  Matrix<Double> & stabilParams_;
};


//**************************************************************************************************************************
//***************** mixed FEM **********************************************************************************************
//**************************************************************************************************************************
class FluidMechPlaneMixedStiffInt_UV : public FluidMechInt
{
public:	

  FluidMechPlaneMixedStiffInt_UV(Double density, Double kinematicViscosity,Matrix<Double> & stabilParams);

  virtual ~FluidMechPlaneMixedStiffInt_UV();
  
  void CalcElementMatrix( Matrix<Double>& elemMat, EntityIterator& ent1, EntityIterator& ent2 );

protected:    
  Matrix<Double> & stabilParams_;
};


//**************************************************************************************************************************
class FluidMechPlaneMixedStiffInt_UQ : public FluidMechInt
{
public:

  FluidMechPlaneMixedStiffInt_UQ(Double density, Double kinematicViscosity, Matrix<Double> & stabilParams);
  
  virtual ~FluidMechPlaneMixedStiffInt_UQ();

  void CalcElementMatrix( Matrix<Double>& elemMat, EntityIterator& ent1, EntityIterator& ent2 );

protected:    
  Matrix<Double> & stabilParams_;
};

//**************************************************************************************************************************
class FluidMechPlaneMixedStiffInt_PV : public FluidMechInt
{
public:
  FluidMechPlaneMixedStiffInt_PV(Double density, Double kinematicViscosity, Matrix<Double> & stabilParams);
  virtual ~FluidMechPlaneMixedStiffInt_PV();
  void CalcElementMatrix( Matrix<Double>& elemMat, EntityIterator& ent1, EntityIterator& ent2 );
protected:    
  Matrix<Double> & stabilParams_;
};



//   /// Class for calculation the contribution to the element stiffness matrix in 3D
// class FluidMech3DLaplNonConsStiffInt : public FluidMechInt
// {
// public:
//   /// Constructor
//   FluidMech3DLaplNonConsStiffInt(Double density, Double dynamicViscosity, Double kinematicViscosity);

//   /// 
//   virtual ~FluidMech3DLaplNonConsStiffInt();

//   /// Calculation of stiffmess matrix
//   void CalcElementMatrix( Matrix<Double>& elemMat,
//                           EntityIterator& ent1, 
//                           EntityIterator& ent2 );
// protected:    
// };


}

#endif // FILE_FLUIDMECHSTIFFINT
