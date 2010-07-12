#ifndef DENSITYFILE_HH_
#define DENSITYFILE_HH_

#include "General/Enum.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Condition.hh"

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
  static DesignSpace* ReadErsatzMaterial(DesignSpace* ersatzMaterial = NULL);

  /** set the current iteration and eventually write it
   * @param current_iteration this info is used as set id */
  void SetCurrent(int current_iteration);

  /** this actually stores the data which is exported as pseudo density file */
  PtrParamNode data;

  private:

  /** Creates the pseudo density node and stores the header */
  PtrParamNode Create(ParamNodeList& des, ParamNodeList& tfs, PtrParamNode regularize);

  /** shall we write the densities for all iterations or overwrite? */
  bool all_iterations_;

  /** shall we write the density file each iteration or only in the destructor.
   * The difference is superfluous file writing .*/
  bool finally_only_;

  /** we don't own this data, therefore we also don't delete it! */
  DesignSpace* ersatzMaterial_;

  /** our filename */
  std::string name_;
};

}


#endif /* DENSITYFILE_HH_ */
