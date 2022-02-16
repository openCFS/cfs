#ifndef CONTAINER_HH_
#define CONTAINER_HH_

#include "MatVec/Vector.hh"
#include "Utils/StdVector.hh"

namespace CoupledField
{
/** common wrapper for base features of Vector and StdVector (extend for std::vector, ...)
 * Note that possibly a template as in std::string ToStringCont(const Cont& cont, const std::string& sep = " ")
 * better suits your usage.
 *
 * The purpose of this class is to prevent copy and paste code for a Vector and StdVector implementation. */
template<class TYPE>
class Container
{
public:
  Container(Vector<TYPE>* in)    { vec = in;}
  Container(Vector<TYPE>& in)    { vec = &in;}
  Container(StdVector<TYPE>* in) { stv = in;}
  Container(StdVector<TYPE>& in) { stv = &in;}

  Container(const Vector<TYPE>* in) { vec = in; }
  Container(const Vector<TYPE>& in) { vec = &in; }
  Container(const StdVector<TYPE>* in) { stv = in; }
  Container(const StdVector<TYPE>& in) { stv = &in; }

  unsigned int GetSize() const {
    return vec != NULL ? vec->GetSize() : stv->GetSize();
  }

  void Resize(unsigned int size) {
    if(vec != NULL)
      const_cast<Vector<TYPE>* >(vec)->Resize(size);
    else
      const_cast<StdVector<TYPE>* >(stv)->Resize(size);
  }

  inline const TYPE& operator[](unsigned int i) const {
    return vec != NULL ? (*vec)[i] : (*stv)[i];
  }

  inline TYPE& operator[](unsigned int i) {
    if(vec != NULL)
      return const_cast<Vector<TYPE>& >(*vec)[i];
    else
      return const_cast<StdVector<TYPE>& >(*stv)[i];
  }


private:
  const Vector<TYPE>*    vec = NULL;
  const StdVector<TYPE>* stv = NULL;
};

} // end of namespace
#endif /* CONTAINER_HH_ */
