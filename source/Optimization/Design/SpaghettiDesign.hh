#ifndef OPTIMIZATION_DESIGN_SPAGHETTIDESIGN_HH_
#define OPTIMIZATION_DESIGN_SPAGHETTIDESIGN_HH_

#include "Optimization/Design/FeaturedDesign.hh"

namespace CoupledField
{
/** Feature mapping variant spaghetti optimization. See also spaghetti.py */
class SpaghettiDesign : public FeaturedDesign
{
public:
  /** @method either SPAGHETTI or SPAGHETTI_PARAM_MAT */
  SpaghettiDesign(StdVector<RegionIdType>& regionIds, PtrParamNode pn, ErsatzMaterial::Method method);

  virtual ~SpaghettiDesign();

    /** read the design and do a conditional mapping */
  int ReadDesignFromExtern(const double* space_in, bool doNotSetAndWriteCurrent = false) override;

  void ToInfo(ErsatzMaterial* em) override;

  /** Called from DensityFile::ReadErsatzMaterial() with load ersatz material (-x)
   * See SetupDesign for expected ordering
   * @param set the set from the density.xml
   * @param lower_violation the maximal violation */
  void ReadDensityXml(PtrParamNode set, double& lower_violation, double& upper_violation) override;

  /** sets radius for spaghetti.py to visualize the stuff */
  void AddToDensityHeader(PtrParamNode pn) override;

  /** the spaghetti python script */
  PyObject* GetPythonModule() override;

  struct Noodle : public Feature 
  {
  public:

    /** NORMAL means we have vector a - which is for 2D, POINTS are 3D as a normal cannot be defined in 3D by a scalar */
    typedef enum { NO_TYPE = -1, NORMAL = 0, POINTS = 1 } Type;

    Noodle();

    const char* GetName() override { return "noodle"; } 

    void ToInfo(PtrParamNode info) override;

    /** helper for a and r */
    void ToInfo(PtrParamNode in, FeatureVariable::Type type, const std::string& label, int dim = -1) const  override;

    /** for debug purpose */
    std::string ToString() const override;

    void Parse(PtrParamNode noodle, int idx) override;

    int GetTotalVariables() const override;

    bool IsExtended() const override { return a.GetSize() > 0; }

    /** Helper for FeatureVariable::Parser() in case of inner */
    double Interpolate(double start, double end, unsigned int idx, unsigned int n_segments) const
    { return start + (end-start)/(n_segments-1) * idx;}

    /** to make it easier generic also for 2d */
    bool HasAlpha() const { return type == POINTS; }

    /** do we have normals or inner points */
    bool HasInner() const { return dim_ == 3; }

    /** the scalar scaling in the geometry projection style */
    FeatureVariable alpha;

    /** the segments-1 normal variables only for NORMAL mode */
    StdVector<FeatureVariable> a;

    /** the segments-1 radius variables only for POINTS, could be extended also for 2D */
    StdVector<FeatureVariable> r;

    /** up to now obtained by dim_ */
    Type type = NO_TYPE;
  };

  /** structured access to the full design space. Note that spaghetti is a plural of noodle.
   * Within the noodles are all Variables referenced in FeaturedDesign::shape_param_ and FeaturedDesign::opt_shape_param_
   * @See SetupDesign(), also for ordering within (opt_)shape_param_ */
  StdVector<Noodle> spaghetti; // instances for FeaturedDesign::features_;

private:

  /** Set up the design, the instances are within Noodle objects stored in spaghetti,
   * FeaturedDesgin::shape_param_ and opt_shape_param_ get references.
   * Ordering via type: First node, then profile, then normal */
  void SetupDesign(PtrParamNode pn) override;


  void SetupParsedFeatures(PtrParamNode base) override;

  /** our detailed implementation for SetupVirtualShapeElementMap() */
  void SetupVirtualShapeElementMapBending(Feature* f, StdVector<Function::Local::Identifier>& vem, StdVector<BaseDesignElement*>& nodes, bool two_signs, int sign_1, int sign_2) override;

  /** Map structure to rho (DesignSpace::data). Sets DesignSpace::data.
   *  Shall be called by ReadDesignFromExtern(). */
  void MapFeatureToDensity() override;

  /** Takes the density gradients and sums it up on the shape variables.
   *  To be called within WriteGradientToExtern().
   *  @param f the function we add the stuff to the gradient. */
  void MapFeatureGradient(Function* f) override;


  /** get generic special results from python.
   * MapFeatureGradient() is called for every function and the generic design might be for the next call only. */
  void PrepareSpecialResults() override;

  /** call cfs_init in python file */
  void PythonInit(PtrParamNode pn);

  /** create/update all spaghetti data in python */
  void PythonUpdateSpaghetti();

  /** gives the result of the python function cfs_info_field_keys */
  StdVector<std::string> PythonGetInfoFieldKeys();

  /** spaghetti radius - needs to confirm profile and curvature of normals */
  double radius = -1;

  /** for boundary functions linear and poly this is the full transition zone 2*h -> move to FeaturedDesign */
  double transition = -1;

  // penalty for smoothmax
  double pen = -1;
  // penalty for RAMP type inverse penalization of python_volume
  double q = 0;

  /** the rhomin we use, extracted from the first density variable. */
  double rhomin = -1;
  double rhomax = -1;

  /** in the 'python' element, the 'option' elements. See PythonKernel::CreatePythonDict */
  StdVector<std::pair<std::string, std::string> > pyopts;

  /** this contains the (later optional) Python script */
  PyObject* module_ = NULL;

  PtrParamNode sp_info_; // our own info

  /** created in ToInfo() */
  boost::shared_ptr<Timer> py_timer;

};

} // end of name space


#endif /* OPTIMIZATION_DESIGN_SPAGHETTIDESIGN_HH_ */
