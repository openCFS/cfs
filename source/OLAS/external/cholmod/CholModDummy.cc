// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

namespace CoupledField {
  // This Dummy function is needed to build a static library for Cholmod
  // support even if explicit template instantiation is turned off and
  // the information about external dependencies of Cholmod could not be
  // communicated to the linker otherwise. 

  void CholModDummyFunc() {}
}
