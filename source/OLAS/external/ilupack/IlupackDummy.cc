namespace CoupledField {
  // This Dummy function is needed to build a static library for Ilupack
  // support even if explicit template instantiation is turned off and
  // the information about external dependencies of Ilupack could not be
  // communicated to the linker otherwise. 

  void IlupackDummyFunc() {}
}
