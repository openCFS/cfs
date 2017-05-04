#ifndef OPTIMIZATION_EIGENINFO_HH_
#define OPTIMIZATION_EIGENINFO_HH_


namespace CoupledField
{

/** For bloch optimization we usually search for minimal and maximal ev within wave vectors.
 * This struct collects some data to allow more detailed output */
struct EigenInfo
{
  EigenInfo()
  {
    col = 0;
    min = 0.0;
    max = 0.0;
  }
  /** the col_idx we found the extremal value */
  unsigned int col;
  /** the minimal value within this ev. Even for bloch=full. It depends on the the bound if this corresponds to col */
  double min;
  /** @see min */
  double max;
};

}



#endif /* OPTIMIZATION_EIGENINFO_HH_ */
