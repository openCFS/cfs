/** This HDF5 wrapper for cfs is adapted from ParaView's CFSReader, which itself is based an a legacy
    HDF5 wrapper from cfs. It has beed upgraded in CFSReader when it was merged upstream (as H5CPP was not
    available any more). 
    
    We use almost the plain hdf5 C-api and add a thin own C++ wrapper around it. The original H5CPP layer
    is not used here any more, as it is too outdated. */ 

#ifndef HDF5COMMON_H
#define HDF5COMMON_H

#include <string>
#include <H5public.h>
#include <H5Ipublic.h>
#include <H5Gpublic.h>

#include "Domain/Results/ResultInfo.hh"
#include "General/Enum.hh"
#include "Utils/StdVector.hh"

namespace CoupledField
{

/** for auto compression: below this size we enforce contiguous data */
const hsize_t HDF5_COMPRESSION_THRESHOLD = 64*1024; 

/** we have only time series for growing datasets where data is actually very small.
 * Too large chunk sizes waste space when we have many results which becomes relevant for small files (e.g. test cases).
 * Note that normal data has 1 exact chunk, as a chunk is the unit for compression */
const hsize_t HDF5_GROWING_CHUNK_SIZE = 1*1024;

/** Wrapper for H5 c-api H5Lget_name_by_idx */
std::string GetObjNameByIdx(hid_t loc, hsize_t idx);

/** Opens group. If it does not exist, an error is throw. 
 * In case check before with H5Lexists
 * @see CreateGroup(..., exists = true) for lazy opening */
hid_t OpenGroup(hid_t loc, const std::string& name);

hid_t OpenDataSet(hid_t loc, const std::string& name);

/** Wrapper for H5 c-api H5Gget_info */
H5G_info_t GetInfo(hid_t group_id);

/** Read h5 attribute via H5LTget_attribute_uint, and other types.
 *
 * Implemented via template specific implementations.
 * @param loc_id group or file
 * @param obj_name element of loc_id object
 * @param data output */
template <typename TYPE>
void ReadAttribute(hid_t loc_id, const std::string& obj_name, const std::string& attrName, TYPE& data);

/** Convenience ReadAttribute which provides own variable space */
template <typename TYPE>
TYPE ReadAttribute(hid_t loc_id, const std::string& obj_name, const std::string& attrName)
{
  TYPE val;
  ReadAttribute(loc_id, obj_name, attrName, val);
  return val;
}

template <typename TYPE>
TYPE ReadAttribute(hid_t loc_id, const std::string& attrName)
{
  TYPE val;
  ReadAttribute(loc_id, ".", attrName, val); // object name is this name
  return val;
}

template <typename TYPE>
void ReadDataSet(hid_t loc, const std::string& name, TYPE* out);

template <typename TYPE>
TYPE ReadDataSet(hid_t loc, const std::string& name);

/** Retrieve rank and dimensionality of a dataset */
StdVector<unsigned int> GetArrayDims(hid_t loc, const std::string& name);

/** Return number of entries of a dataset / array  */
unsigned int GetNumberOfEntries(hid_t id, const std::string& name);

/** Retrieve array data from a dataset */
template <typename TYPE>
void ReadArray(hid_t loc, const std::string& name, StdVector<TYPE>& data);

template <typename TYPE>
StdVector<TYPE> ReadArray(hid_t loc, const std::string& name) {
  StdVector<TYPE> data;
  ReadArray(loc, name, data);
  return data;
}

/** Obtain grid result group for specified multisequence step.
 *
 * @param msStep 1-based multisequence step. There is always at least the first one.
 * @param isHistory indicated region result instead of element/node result */
hid_t GetMultiStepGroup(hid_t root, unsigned int msStep, bool isHistory);

/** Obtain grid result group for specified step in a given multistep element/node result
 * @param stepNum the 0-based step for nodal and element results
 * @see GetMultiStepGroup() */
hid_t GetStepGroup(hid_t root, unsigned int msStep, unsigned int stepNum);

// =============================================================================
// WRITE helpers
// =============================================================================

/** Creates a group with the full path.
  * @param use_existing if false, an exception is thrown if already exists. */
hid_t CreateGroup(hid_t loc, const std::string& name, bool use_existing = false);

/** Write scalar data type. */
template <typename TYPE>
void WriteAttribute(hid_t loc_id, const std::string& obj_name, const std::string& attrName, const TYPE& data);

/** Write scalar data type. Assumes current object ('.') */
template <typename TYPE>
void WriteAttribute(hid_t loc_id, const std::string& attrName, const TYPE& data) {
  WriteAttribute(loc_id, ".", attrName, data);
}

/** this is the (not meant for public use) generalization of WriteDataSet1D/2D */
template <typename T>
void WriteDataSetND(hid_t loc, const std::string& name,
                           int rank, hsize_t* dims, hsize_t totalElems,
                           const T* data, int compressionLevel);


/** Write 1D data by pointer and number of entried. Optionally with compression.
 * This is not meant for growing datasets.
 * @param compressionLevel 0 = off, 1 = fast, anything larger is much slower with usually only marginally less size effect
                           when data is below HDF5_COMPRESSION_THRESHOLD we silently enforce 0. Check with h5info.py */
template <typename TYPE>
void WriteDataSet1D(hid_t loc, const std::string& name, const TYPE* data, hsize_t count, int compressionLevel = 1) {
  hsize_t dims[1] = {count};
  WriteDataSetND(loc, name, 1, dims, count, data, compressionLevel);
}

/** convenience function for unsigned int count */
template <typename TYPE>
void WriteDataSet1D(hid_t loc, const std::string& name, const TYPE* data, unsigned int count, int compressionLevel = 1) {
  WriteDataSet1D(loc, name, data, (hsize_t) count, compressionLevel);
}

/** write a single value dataset - e.g. where attribute would be more appropriate 
  * but it is dataset in the legacy code. */
template <typename TYPE>
void WriteSingleDataSet(hid_t loc, const std::string& name, TYPE data) {
  WriteDataSet1D(loc, name, &data, (hsize_t) 1);
}

void WriteSingleDataSet(hid_t loc, const std::string& name, const std::string& data);

/** Convenience function to write a string array dataset */
void WriteStringArray(hid_t loc, const std::string& name, const StdVector<std::string>& data);

/** Write a 2-D dataset (rows x cols). 
 * We assume that we exactly know how much to write with no later addition. 
 * @see WriteDataSet for compressionLevel parameter. */ 
template <typename TYPE>
void WriteDataSet2D(hid_t loc, const std::string& name, hsize_t rows, hsize_t cols, const TYPE* data, int compressionLevel = 1) {
  hsize_t dims[2] = {rows, cols};
  WriteDataSetND(loc, name, 2, dims, rows * cols, data, compressionLevel);
}

template <typename TYPE>
void WriteDataSet2D(hid_t loc, const std::string& name, unsigned int rows, unsigned int cols, const TYPE* data, int compressionLevel = 1) {
  WriteDataSet2D(loc, name, (hsize_t) rows, (hsize_t) cols, data, compressionLevel);
}

/** Create or append to a 1-D extensible dataset.
 * If the dataset does not yet exist it is created (empty). If data != nullptr
 * and count > 0, the elements are appended to the (possibly new) dataset.*/
template <typename TYPE>
void WriteGrowingDataSet1D(hid_t loc, const std::string& name,
                         const TYPE* data = nullptr, hsize_t count = 0);

/** Create or append rows to a 2-D extensible dataset.
 * For first calls, newRows can be 0 and data nullptr nut numCols must always be the same
 * @see WriteGrowingDataSet1D
 * @param numCols number of columns (fixed, must match existing dataset) */
template <typename TYPE>
void WriteGrowingDataSet2D(hid_t loc, const std::string& name, hsize_t newRows, hsize_t numCols, const TYPE* data);

/** Maps EntityUnknownType to the HDF5 group name string (e.g. NODE -> "Nodes").
 *  SURF_ELEM and ELEMENT both map to "Elements"; REGION_AVERAGE and REGION both map to "Regions". */
extern Enum<ResultInfo::EntityUnknownType> entityGroupNameEnum;

/** Maps EntityUnknownType to the legacy integer code stored in the HDF5 file. */
int MapEntityTypeToInt(ResultInfo::EntityUnknownType t);

/** Inverse of MapEntityTypeToInt. */
ResultInfo::EntityUnknownType MapEntityTypeFromInt(int t);

/** Maps EntryType to the legacy integer code stored in the HDF5 file. */
int MapEntryTypeToInt(ResultInfo::EntryType t);

/** Inverse of MapEntryTypeToInt. */
ResultInfo::EntryType MapEntryTypeFromInt(int t);

} // end of namespace CoupledField

#endif // HDF5COMMON_H
