namespace CoupledField {
  // This Dummy function is needed to build a static library for Cholmod
  // support even if explicit template instantiation is turned off and
  // the information about external dependencies of Cholmod could not be
  // communicated to the linker otherwise. 

  void CholModDummyFunc() {}
}
