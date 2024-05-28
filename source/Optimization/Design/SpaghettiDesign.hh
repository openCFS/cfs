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

  /** Set's up distance neighborhood (all nodes of a shape) */
  void SetupVirtualShapeElementMap(Function* f, StdVector<Function::Local::Identifier>& virtual_element_map, Function::Local::Locality locality) override;

  /** sets radius for spaghetti.py to visualize the stuff */
  void AddToDensityHeader(PtrParamNode pn) override;

  /** the spaghetti python script */
  PyObject* GetPythonModule() override;

  typedef enum { NO_TIP = -1, START, END, INNER } Tip;

  static Enum<Tip> tip;

  /** This represents a single design variable (which can be fixed)
   * set dof and tip if applicable
   * Public as used in DensityFile  */
  class Variable : public ShapeParamElement
  {
  public:

    static std::string ToString(const StdVector<Variable>& vec, bool show_fixed = true);

    static Vector<double> AsVector(const StdVector<Variable>& vec);

    /** optional late constructor */
    void InitInnerNode(int noodle_idx, Dof dof, double val, double lower, double upper);

    void Parse(PtrParamNode pn, int noodle, double interpolate_value = -12.34);

    /** all other BaseDesignElement children do this in the constructor.
     * @see SetOptIndex() */
    void SetIndex(int idx) { index_ = idx; }

    void ToInfo(PtrParamNode in) const;

    /** for debug info */
    std::string ToString() const override;


    std::string GetLabel() const override;

    /** does not apply for all */
    Tip tip = NO_TIP;

    /** subject to optimization or fixed */
    bool fixed = false;

    /** noodle index */
    int noodle = -1;
  };

  struct Noodle
  {
  public:

    /** NORMAL means we have vector a - which is for 2D, POINTS are 3D */
    typedef enum { NO_TYPE = -1, NORMAL = 0, POINTS = 1 } Type;

    Noodle();

    void Parse(PtrParamNode pn, int idx);

    void ToInfo(PtrParamNode info);

    /** helper for a and r */
    void ToInfo(PtrParamNode in, Variable::Type type, const std::string& label, int dim = -1) const;
    void CompareToInfoHelper(const Variable& v0, const Variable& test, std::string& fixed, std::string& lower, std::string& upper) const;

    /** for debug purpose */
    std::string ToString();

    int GetTotalVariables() const;

    int GetOptVariables() const { return opt_variables_; }

    /** check all components of point */
    bool IsFixed(const StdVector<Variable>& point) const;

    int CountNonFixed(const StdVector<Variable>& point) const;

    /** Helper for Variable::Parser() in case of inner */
    double Interpolate(double start, double end, unsigned int idx, unsigned int n);

    /** to make it easier generic also for 2d */
    bool HasAlpha() const { return type == POINTS; }

    /** do we have normals or inner points */
    bool HasInner() const { return dim_ == 3; }

    /** Variable objects are the read design variables in a technical sense.
      points contains all nodal coordinates with first point (2d or 3d) the tip start P and the last the tip end Q.
      In the POINTS case (3d) the points in between are the alternative variables to normals a.
      The in between points are either automatically created with interpolated values or a single inner point
      gives initial/upper/lower - all also by mathparser. inner is only necessary if the bounds for P and Q are not clear */
    StdVector<StdVector<Variable> > points;

    /** the scalar profile variable */
    Variable p;

    /** the scalar scaling in the geometry projection style */
    Variable alpha;

    /** the segments-1 normal variables only for NORMAL mode */
    StdVector<Variable> a;

    /** the segments-1 radius variables only for POINTS, could be extended also for 2D */
    StdVector<Variable> r;

    /** shape index */
    int idx = -1;

    /** up to now obtained by dim_ */
    Type type = NO_TYPE;

  private:
    int opt_variables_ = -1;
  };

  /** structured access to the full design space. Note that spaghetti is a plural of noodle.
   * Within the noodles are all Variables referenced in FeaturedDesign::shape_param_ and FeaturedDesign::opt_shape_param_
   * @See SetupDesign(), also for ordering within (opt_)shape_param_ */
  StdVector<Noodle> spaghetti;
  StdVector<int> opt_indices;

private:

  /** Set up the design, the instances are within Noodle objects stored in spaghetti,
   * FeaturedDesgin::shape_param_ and opt_shape_param_ get references.
   * Ordering via type: First node, then profile, then normal */
  void SetupDesign(PtrParamNode pn) override;

  /** Map structure to rho (DesignSpace::data). Sets DesignSpace::data.
   *  Shall be called by ReadDesignFromExtern(). */
  void MapFeatureToDensity() override;

  /** Takes the density gradients and sums it up on the shape variables.
   *  To be called within WriteGradientToExtern().
   *  @param f the function we add the stuff to the gradient. */
  void MapFeatureGradient(const Function* f) override;

  /** get generic special results from python.
   * MapFeatureGradient() is called for every function and the generic design might be for the next call only. */
  void PrepareSpecialResults() override;

  /** Helper for SetupDesign(). Stores the variable to shape_param_ and opt_shape_param_ */
  void AddVariable(Variable* var);

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
