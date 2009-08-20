# Sample tcl scripting file demonstrating how to print specific
# element matrices / vectors at runtime
# author: Andreas Hauck
# Last modified: 20.8.2009

proc CFS_PdeInit { pdeNname } {
  cfs assemble printElemMatVec {1}
}

proc CFS_AssembleMat { form elemNum1 elemNum2 elemMat } {
  puts -nonewline "Elementmatrix of bilinearform '$form'"
  puts "of elements #$elemNum1 and #$elemNum2 is"
  puts $elemMat
}

proc CFS_AssembleRhs { form elemNum elemVec } {
  puts -nonewline "Elementmatrix of linearform '$form'"
  puts "of element #$elemNum is"
  puts $elemVec
}
