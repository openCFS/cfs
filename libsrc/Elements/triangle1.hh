#ifndef FILE_TRIANGLE1_2001
#define FILE_TRIANGLE1_2001
 
namespace CoupledField
{
  class BaseElem;
/// Description of triangle element with 1 order interpolation 
  class Triangle1 : public BaseElem
{
public:
  ///
  Triangle1(ShortInt aintegtype);
 
  virtual ~Triangle1();
 
  ShortInt numnodes;
  ShortInt numedges;
  ShortInt numfaces;
  ShortInt numintpoints;
 
  virtual void Init();
  virtual void SetIntPoints();
  virtual void SetShapeFnc();
  virtual void SetDShapeFnc();
 
protected:
 
private:
  ShortInt IntegType;
 
};
 
}
 
#endif // FILE_TRIANGLE1_2001
