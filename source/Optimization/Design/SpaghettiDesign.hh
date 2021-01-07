#ifndef OPTIMIZATION_DESIGN_SPAGHETTIDESIGN_HH_
#define OPTIMIZATION_DESIGN_SPAGHETTIDESIGN_HH_

#include "Optimization/Design/FeaturedDesign.hh"

#include <def_use_embedded_python.hh>
#ifdef USE_EMBEDDED_PYTHON
  #define PY_SSIZE_T_CLEAN // https://docs.python.org/3/c-api/intro.html
  #include <Python.h>
#endif

namespace CoupledField
{
/** Feature mapping variant spaghetti optimization. See also spaghetti.py */
class SpaghettiDesign : public FeaturedDesign
{
public:
  SpaghettiDesign(StdVector<RegionIdType>& regionIds, PtrParamNode pn, ErsatzMaterial::Method method = ErsatzMaterial::SPAGHETTI);

  virtual ~SpaghettiDesign();

  /** read the design and do a conditional mapping */
  int ReadDesignFromExtern(const double* space_in) override;

  void ToInfo(ErsatzMaterial* em) override;

  /** Called from DensityFile::ReadErsatzMaterial() with load ersatz material (-x)
   * See SetupDesign for expected ordering
   * @param set the set from the density.xml
   * @param lower_violation the maximal violation */
  void ReadDensityXml(PtrParamNode set, double& lower_violation, double& upper_violation) override;

  /** This is the variant of Function::Local::SetupVirtualElementMap() for slope constraints on ShapeParamElements.
   * This function is called within Function::Local() constructor, therefore Function::GetLocal() cannot work yet!
   * @param locality just given to assert() it is PREV_AND_NEXT */
  void SetupVirtualShapeElementMap(Function* f, StdVector<Function::Local::Identifier>& virtual_element_map, Function::Local::Locality locality) override {}

  /** sets radius for spaghetti.py to visualize the stuff */
  void AddToDensityHeader(PtrParamNode pn) override;

  typedef enum { NO_TIP = -1, START, END } Tip;

  static Enum<Tip> tip;

  /** This represents a single design variable (which can be fixed)
   * set dof and tip if applicable
   * Public as used in DensityFile  */
  class Variable : public ShapeParamElement
  {
  public:
    void Parse(PtrParamNode pn, int noodle);

    /** all other BaseDesignElement childs do this in the constructor.
     * @see SetOptIndex() */
    void SetIndex(int idx) { index_ = idx; }

    void ToInfo(PtrParamNode in) const;

    /** does not apply for all */
    Tip tip = NO_TIP;

    /** subject to optimization or fixed */
    bool fixed = false;

    /** noodle index */
    int noodle = -1;
  };

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

  struct Noodle
  {
  public:
    void Parse(PtrParamNode pn, int idx);

    void ToInfo(PtrParamNode info);

    /** priliminarily for debug purpose */
    std::string ToString();

    int GetTotalVariables() const { return 5 + a_var.GetSize(); }
    int GetOptVariables() const { return opt_variables_; }

    /** copy from Variable objects to convenience variables P,Q,p,a */
    void Update();

    // Variable objects are the read design variables in a technical sense

    Variable px; // make this vectors for 2D/3D
    Variable py;
    Variable qx;
    Variable qy;
    Variable p_var;
    StdVector<Variable> a_var;

    // the following data are convenience objects filled by Update() from the Variable objects.
    // we use capital letters for vectors and points following math and spaghetti.py
    Vector<double> P; // start point coordinate
    Vector<double> Q; // end point coordinate
    double         p = -1; // profile is doubled width
    StdVector<double> a; // list normals



    /** shape index */
    int idx = -1;
  private:
    int opt_variables_ = -1;
  };

  /** Helper for SetupDesign(). Stores the variable to shape_param_ and opt_shape_param_ */
  void AddVariable(Variable* var);

  /** call cfs_init in python file */
  void PythonInit(PtrParamNode pn);

  /** python part of the desctructor */
  void PythonDestructor();

  /** create/update all spaghetti data in python */
  void PythonUpdateSpaghetti();

  /** structured access to the full design space. Note that spaghetti is a plural of noodle.
   * Within the noodles are all Variables referenced in FeaturedDesign::shape_param_ and FeaturedDesign::opt_shape_param_
   * @See SetupDesign(), also for ordering within (opt_)shape_param_ */
  StdVector<Noodle> spaghetti;

  /** spaghetti radius - needs to confirm profile and curvature of normals */
  double radius = -1;

  /** for boundary functions linear and poly this is the full transition zone 2*h -> move to FeaturedDesign */
  double transition = -1;

  /** the rhomin we use, extracted from the first density variable. */
  double rhomin = -1;

  /** in the 'python' element, the 'file' attribute */
  std::string pythonfile;

  /** in the 'python' element, the 'option' elements. See PythonTools::CreatePythonDict */
  StdVector<std::pair<std::string, std::string> > pyopts;

  /** this contains the (later optional) Python module */
  PyObject* python = NULL;



  /** created in ToInfo() */
  boost::shared_ptr<Timer> py_timer;

};

} // end of name space


#endif /* OPTIMIZATION_DESIGN_SPAGHETTIDESIGN_HH_ */
