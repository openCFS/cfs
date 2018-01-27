#include "BaseMatrix.hh"
#include "General/Enum.hh"

namespace CoupledField {

// Definition of matrix structure type mappings
static EnumTuple structureTypeTuples[] = 
{
 EnumTuple( BaseMatrix::NOSTRUCTURETYPE, "noStructureType" ),
 EnumTuple( BaseMatrix::SBM_MATRIX, "SBM" ),
 EnumTuple( BaseMatrix::SPARSE_MATRIX, "sparse"),
 EnumTuple( BaseMatrix::DENSE_MATRIX, "dense")
};

Enum<BaseMatrix::StructureType> BaseMatrix::structureType = \
Enum<BaseMatrix::StructureType>("Matrix Structure Types",
    sizeof(structureTypeTuples) / sizeof(EnumTuple),
    structureTypeTuples); 


std::string BaseMatrix::ToInfoString() const
{
  return "st=" + structureType.ToString(GetStructureType()) + " et=" + entryType.ToString(GetEntryType()) + " mu=" + boost::lexical_cast<std::string>(GetMemoryUsage());
}


// Definition of matrix entry type mappings
static EnumTuple entryTypeTuples[] = 
{
 EnumTuple( BaseMatrix::NOENTRYTYPE, "noEntryType" ),
 EnumTuple( BaseMatrix::INTEGER, "integer" ),
 EnumTuple( BaseMatrix::DOUBLE, "double" ),
 EnumTuple( BaseMatrix::COMPLEX, "complex" ),
 EnumTuple( BaseMatrix::F77REAL8, "F77REAL8" ),
 EnumTuple( BaseMatrix::F77COMPLEX16, "F77COMPLEX16")
};

Enum<BaseMatrix::EntryType> BaseMatrix::entryType = \
Enum<BaseMatrix::EntryType>("Matrix Entry Types",
    sizeof(entryTypeTuples) / sizeof(EnumTuple),
    entryTypeTuples); 

// Definition of matrix storage type mappings
static EnumTuple storageTypeTuples[] = 
{
 EnumTuple( BaseMatrix::NOSTORAGETYPE, "noStorageType" ),
 EnumTuple( BaseMatrix::SPARSE_SYM, "sparseSym" ),
 EnumTuple( BaseMatrix::SPARSE_NONSYM, "sparseNonSym"),
 EnumTuple( BaseMatrix::SKYLINE_SYM, "skylineSym" ),
 EnumTuple( BaseMatrix::SKYLINE_NONSYM, "skylineNonSym"),
 EnumTuple( BaseMatrix::DIAG, "diag" ),
 EnumTuple( BaseMatrix::VAR_BLOCK_ROW, "variableBlockRow" ),
 EnumTuple( BaseMatrix::LAPACK_GBMATRIX, "lapackGBMatrix" )
};

Enum<BaseMatrix::StorageType> BaseMatrix::storageType = \
Enum<BaseMatrix::StorageType>("Matrix Storage Types",
    sizeof(storageTypeTuples) / sizeof(EnumTuple),
    storageTypeTuples); 

// Definition of matrix storage type mappings
static EnumTuple outputFormatTuples[] = 
{
 EnumTuple( BaseMatrix::MATRIX_MARKET, "matrix-market" ),
 EnumTuple( BaseMatrix::HARWELL_BOEING, "harwell-boeing" ),
 EnumTuple( BaseMatrix::PLAIN, "plain" )
};


Enum<BaseMatrix::OutputFormat> BaseMatrix::outputFormat = \
Enum<BaseMatrix::OutputFormat>("Matrix Output Formats",
    sizeof(outputFormatTuples) / sizeof(EnumTuple),
    outputFormatTuples); 

} // end of namespace
