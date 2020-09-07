#ifndef DATAINOUT_DEPENCIES_HH_
#define DATAINOUT_DEPENCIES_HH_

#include <string>
#include "Utils/StdVector.hh"

using std::string;

namespace CoupledField
{

/* the purpose of this class is to organize the dependencies, cfs was built with */
class Dependencies
{
public:
  Dependencies();

  /** write a file in CMake format with all set(USE_* ON/OFF) to be included via CMake
   * @return success */
  bool WriteCMakeUSE(const string& filename);

  typedef enum { NOT_SET = 0, NOT_KNOWN = 1, EASY = 2, CFS = 3, COMMERCIAL = 4, CLOSED = 5, BSD = 6, LGPL = 7,
                 GPL = 8, ISSL = 9, BOOST = 10, APACHE2 = 11, MIT = 12, ECLIPSE = 13 } License;

  struct Dependency
  {
    // for StdVector only
    Dependency() {};

    Dependency(const string& name, const string& cmake, License lic);

    // this constructor sets automatically to active
    Dependency(const string& name, const string& cmake, const string& version, License lic, const string comment = "", const string date = "");

    /** is this a user switchable option? */
    bool IsSwitchable() const;

    // set the version and automatically sets active
    void SetVersion(const string& val) {
      this->version = val;
      this->active = true;
    }
    void SetVersion(const string& version, const string& sub, const string& minor);
    void SetVersion(int version, int sub = -1, int minor = -1);

    string name;     // like boost
    string cmake;    // like USE_CGAL or empty when implicit
    string version;  // like 1.73.0
    string date;     // optional date information
    string comment;  // any additional information
    bool   active = false;

    License lic = NOT_SET;
  };

  // Allows the binary a free distribution? Not for GPL, COMMERICAL and CLOSED
  static bool IsDistributable(License lic);

  // are all deps distributable?
  bool IsDistributable() const;

  // do we have a specific licence type?
  bool HasLicense(License test) const;

  // for debugging
  void Dump() const;

  StdVector<Dependency> data;

private:
  /** fill data by cecking all the build dependencies. Called by constructor */
  void ReadSetting();

};


} // end of namespace

#endif /* DATAINOUT_DEPENCIES_HH_ */
