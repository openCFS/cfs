#ifndef DENSITYFILE_HH_
#define DENSITYFILE_HH_

#include <stddef.h>
#include <string>

#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField {
class DesignSpace;
}  // namespace CoupledField

namespace CoupledField
{

/** The density file encapsulates the stuff correlated with the density.xml
 * file. Feel free to make it complete :) */
class DensityFile
{
public:

  /** holds service methods for writing the density file */
  DensityFile(DesignSpace* ersatzMaterial,
               PtrParamNode export_pn,
               ParamNodeList& des,
               ParamNodeList& tfs,
               PtrParamNode regularize);

  /** write the file when finally_only */
  ~DensityFile();

  /** are we in the read density file case? Either because
   * "loadErsatzMaterial" is defined in the problem xml file or the
   * command line option -x rsp. --ersatz is used. */
  static bool NeedLoadErsatzMaterial();

  /** Reads an ersatz material file.
   * @param ersatzMaterial if given the data is overwritten otherwise it is created
   * @return the parameter if given or a new one where one needs remove it! */
  static DesignSpace* ReadErsatzMaterial(DesignSpace* space = NULL);
  
  /** set the current iteration and eventually write it. Save to call multiple times with the same parameter
   * @param current_iteration this info is used as set id */
  void SetAndWriteCurrent(int current_iteration);

  /** this actually stores the data which is exported as pseudo density file */
  PtrParamNode data;

  private:

  /** Creates the pseudo density node and stores the header */
  PtrParamNode Create(ParamNodeList& des, ParamNodeList& tfs, PtrParamNode regularize, bool non_desig_vicinity);

  /** create a DesignSpace in case we load ersatz material for a pure simulation. Does not work for shape mapping, ...*/
  static DesignSpace* CreateDesignSpace(bool force_region, const PtrParamNode& pn, const ParamNodeList& elems, const PtrParamNode& xml);
  static bool ReadDensity(PtrParamNode pn, const ParamNodeList& elems, bool force_region, DesignSpace* space,
      double& lower_violation, double& upper_violation);

  /** shall we write the densities for all iterations or overwrite? */
  bool all_iterations_;

  /** shall we write the density file each iteration or only in the destructor.
   * The difference is superfluous file writing .*/
  bool finally_only_;
  
  /** we don't own this data, therefore we also don't delete it! */
  DesignSpace* space_;

  /** our filename */
  std::string name_;

  /** the last set iteration */
  int last_set_iter;
};

}


#endif /* DENSITYFILE_HH_ */
