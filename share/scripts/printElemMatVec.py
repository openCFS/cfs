# Sample tcl scripting file demonstrating how to print specific
# element matrices / vectors at runtime
# author: Andreas Hauck
# Last modified: 20.8.2009

def CFS_PdeInit( pdeNname ):
  cfs ("assemble", "printElemMatVec", "1 2 3 4" )


def CFS_AssembleMat ( form, elemNum1, elemNum2, elemMat ):
  print "Elementmatrix of bilinearform",form
  print " of elements",elemNum1,"and",elemNum2," is"
  print elemMat


def CFS_AssembleRhs ( form, elemNum, elemVec ):
  print "Elementmatrix of linearform",form
  print " of element",elemNum," is"
  print elemVec
  