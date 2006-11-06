#ifndef FILE_CFSVECTOR_2004
#define FILE_CFSVECTOR_2004

#include <iostream>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/is_abstract.hpp>

#include "tools.hh"
#include <vector>

namespace CoupledField
{


  //! Abstract interface class for a general dense vector
  /*! Abstract interface class for a general dense vector
    \note Altough this class provides a general interface
    to the vector, one should always prefer a cast into
    the current type of vector, e.g. 
    dynamic_cast<Vector<Double>*>(ptBaseVector)
  */
  class CFSVector {

  public:

    //! Default constructor
    CFSVector(){};

    //! Constructor with size
    explicit CFSVector(int size){};

    //! Default destructor
    virtual ~CFSVector(){};

    //! Get entry type of vector
    virtual EntryType::ScalarType GetEntryType() const = 0;

    //! Return the length of the vector
    virtual UInt GetSize() const = 0;

    //! Set the lenght of the vector
    /*!
      \param size (input) Lengh of vector
    */
    //! \note the entries are set to zero afterwards!
    virtual void Resize(const UInt size) = 0;

    //! Initalizes the vector with a given entry
    /*!
      \param entry (input) Entry vector gets inialized with
    */
    //! \note this method does not change the size of the vector!
    virtual void Init(const Double entry = 0.0)
    {Error("CFSVector::Init(): Not implemented here",__FILE__,__LINE__);}
  
    //! Set the entry i of the vector to the given value (x[i] = s)
    /*!
      \param i (input) Index of entry s
      \param s (input) Entry to be set on position i
    */
    virtual void SetEntry(const UInt i, const Double &s)
    {Error("CFSVector::SetEntry(): Not implemented here",__FILE__,__LINE__);}

    //! Get the entry i of the vector on the given value (ret = x[i])
    /*!
      \param i (input) Index of entry s
      \param ret (output) Entry on position i
    */
    virtual void GetEntry(const UInt i, Double &ret) const
    {Error("CFSVector::GetEntry(): Not implemented here",__FILE__,__LINE__);}
  
    //! Add s to i-th vector entry (x[i] += s)
    /*!
      \param i (input) Index of entry s
      \param s (input) Value to be added to x[i]
    */
    virtual void AddEntry(const UInt i, const Double &s)
    {Error("CFSVector::AddEntry(): Not implemented here",__FILE__,__LINE__);}
  
    //! Multiply the i-th vector entry with s (x[i] *= s)
    /*!
      \param i (input) Index of entry s
      \param s (input) Factor, which i-the entry gets multiplied with
    */
    virtual void MultEntry(const UInt i, const Double &s)
    {Error("CFSVector::MultEntry(): Not implemented here",__FILE__,__LINE__);} 
  
    //! Multiply the i-th vector entry with a and add s on it (x[i] = a*x[i]+s)
    /*!
      \param i (input) Index of entry s
      \param a (input) Factor the i-the entry gets multiplied with
      \param s (input) Value to be added to a*x[i]
    */
    virtual void MultAddEntry(const UInt i, const Double &a, const Double &s)
    {Error("CFSVector::MultAddEntry(): Not implemented here",__FILE__,__LINE__);}   

    //! Adds another basevector to itself (x = x+y)
    /*! 
      \param y (input) Addend to the vector
    */
    virtual void Add(const CFSVector& y) = 0;

    //@{
    //! Adds the multiple of another vector to itself (x = x +a*y)
    /*!
      \param a (input) Factor for scaling vector y
      \param y (input) Addend to the vector
    */
    virtual void Add(const Double a, const CFSVector &y)
    {Error("CFSVector::Add(): Not implemented here",__FILE__,__LINE__);}

    virtual void Add(const Complex a, const CFSVector &y)
    {Error("CFSVector::Add(): Not implemented here",__FILE__,__LINE__);}
    //@}

    //@{
    //! Replaces the vector by the sum of two scaled vectors (x = a*y+b*z)
    /*!
      \param a (input) Factor for vector y
      \param y (input) Vector scaled by factor a
      \param b (input) Factor for vector z
      \param z (input) Vector scaled by factor b
    */
    virtual void Add( const Double a, const CFSVector& y,
                      const Double b, const CFSVector& z ) 
    {Error("CFSVector::Add(): Not implemented here",__FILE__,__LINE__);}
    virtual void Add( const Complex a, const CFSVector& y,
                      const Complex b, const CFSVector& z )
    {Error("CFSVector::Add(): Not implemented here",__FILE__,__LINE__);}
    //@}
  
    //@{
    //! Scales the vector itself and adds another one y (x = a*x + y)
    /*!
      \param a (input) Factor for scaling the vector itself
      \param y (input) Addend for the vector (gets not scaled)
    */
    virtual void Axpy(const Double a, const CFSVector &y)
    {Error("CFSVector::Axpy(): Not implemented here",__FILE__,__LINE__);}
    virtual void Axpy(const Complex a, const CFSVector &y)
    {Error("CFSVector::Axpy(): Not implemented here",__FILE__,__LINE__);}
    //@}
  
    //! Performes the dot/inner product of the vector itself with y (=x^T*y)
    /*! 
      \param y (input) Vector to perform innter product with
      \param result (output) Result of x^T * y
    */
    virtual void Inner(const CFSVector &y, Double &result) const
    {Error("CFSVector::Add(): Not implemented here",__FILE__,__LINE__);}
  
  
    //! Compute Euclidean L2 norm of this vector object
    virtual Double NormL2() const  = 0;
      
#define DECL_BASEVECTOR_FCN(TYPE)                                                               \
 virtual void Init(const TYPE entry = TYPE())                                           \
  {Error("CFSVector::Init(): Not implemented here",__FILE__,__LINE__);}         \
  virtual void SetEntry(const UInt i, const TYPE &s)                                 \
  {Error("CFSVector::SetEntry(): Not implemented here",__FILE__,__LINE__);}             \
  virtual void GetEntry(const UInt i,  TYPE &ret) const                              \
  {Error("CFSVector::GetEntry(): Not implemented here",__FILE__,__LINE__);}             \
  virtual void AddEntry(const UInt i, const TYPE &s)                                 \
  {Error("CFSVector::AddEntry(): Not implemented here",__FILE__,__LINE__);}             \
  virtual void MultEntry(const UInt i, const TYPE &s)                                \
  {Error("CFSVector::MultEntry(): Not implemented here",__FILE__,__LINE__);}            \
  virtual void MultAddEntry(const UInt i, const TYPE &a, const TYPE &s)              \
  {Error("CFSVector::MultAddEntry(): Not implemented here",__FILE__,__LINE__);}        \
  virtual void Inner(const CFSVector &y, TYPE &result)                                 \
  {Error("CFSVector::Inner(): Not implemented here",__FILE__,__LINE__);}

    DECL_BASEVECTOR_FCN(Integer)
    DECL_BASEVECTOR_FCN(Complex)

  private:    
    // =======================================================================
    // SERIALIZATION FUNCTIONS
    // =======================================================================
    // These functions allow us to write a vector directly
    // into an boost::archive, for saving on a disk or in a 
    // iostream object
    
    //! allow serialization class to access vector entries
    friend class boost::serialization::access;
    
    //! Saving internal state into a boost::archive
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const {}
    
    //! Reading internal state from a boost::archive
    template<class Archive>
    void load(Archive & ar, const unsigned int version) {}
    
    //! The following statement is needed for boost
    BOOST_SERIALIZATION_SPLIT_MEMBER()

  };
  
} // end of namespace

//! Indicate that this is a pure virtual class
BOOST_IS_ABSTRACT( CoupledField::CFSVector )
#endif  


