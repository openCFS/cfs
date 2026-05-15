#include <algorithm>
#include <cassert>

#include <hdf5.h>
#include <hdf5_hl.h>

#include "CommonHDF5.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Domain/Results/ResultInfo.hh"
#include "General/Enum.hh"
#include "General/Exception.hh"
#include "Utils/StdVector.hh"

namespace CoupledField
{

DEFINE_LOG(h5c, "h5common")

// Map C++ type -> HDF5 native type constant - string usually needs special code
template <typename T> hid_t NativeType();
template <> hid_t NativeType<double>()       { return H5T_NATIVE_DOUBLE; }
template <> hid_t NativeType<int>()          { return H5T_NATIVE_INT; }
template <> hid_t NativeType<unsigned int>() { return H5T_NATIVE_UINT; }

std::string GetObjNameByIdx(hid_t loc, hsize_t idx)
{
  // We first obtain the length of the name, then by the same function write the name to a buffer
  // and third copy it to a std::string we return. There is no shortcut possible

  // nullptr for name gives the length
  ssize_t nameLen = H5Lget_name_by_idx(loc, ".", H5_INDEX_NAME, H5_ITER_NATIVE, idx, nullptr, 0, 0);
  if (nameLen < 0)
    EXCEPTION("Was not able to determine name");
  
  // now, allocate buffer to get the name, we may not write to std::string data directly.
  // std::make_unique is not C++11
  // https://stackoverflow.com/questions/1042940/writing-directly-to-stdstring-internal-buffers
  StdVector<char> buf(nameLen + 1); // the trailing NULL is written but not counted

  if (H5Lget_name_by_idx(loc, ".", H5_INDEX_NAME, H5_ITER_NATIVE, idx, buf.GetPointer(), nameLen + 1, 0) < 0)
  {
    EXCEPTION("error obtaining obj name with index " << idx << " of size " << nameLen);
  }

  std::string name(buf.GetPointer());
  return name;
}

hid_t OpenGroup(hid_t loc, const std::string& name)
{
  hid_t ret = H5Gopen(loc, name.c_str(), H5P_DEFAULT);
  if (ret < 0)
    throw Exception("cannot open hdf5 group '" + name + "'");
  return ret;
}


hid_t OpenDataSet(hid_t loc, const std::string& name)
{
  hid_t ret = H5Dopen(loc, name.c_str(), H5P_DEFAULT);
  if (ret < 0)
    throw Exception("cannot open hdf5 data set '" + name + "'");
  return ret;
}

//-----------------------------------------------------------------------------
H5G_info_t GetInfo(hid_t group_id)
{
  H5G_info_t info;
  if (H5Gget_info(group_id, &info) < 0)
  {
    throw Exception("cannot get hdf5 group info");
  }

  return info;
}

//-----------------------------------------------------------------------------
// FileInfo() is a callback function for H5Giterate().
// By comment we prevent warnings for unused loc_id
herr_t FileInfo(hid_t /* loc_id */, const char* name, void* opdata)
{
  StdVector<std::string>* names = reinterpret_cast<StdVector<std::string>*>(opdata);
  assert(names != nullptr);
  names->Push_back(name);

  return 0;
}

//-----------------------------------------------------------------------------
StdVector<std::string> ParseGroup(hid_t loc_id, const std::string& name)
{
  StdVector<std::string> result;

  H5Giterate(loc_id, name.c_str(), nullptr, FileInfo, &result);

  return result;
}

//-----------------------------------------------------------------------------
bool TestGroupChild(hid_t loc_id, const std::string& group, const std::string& child)
{
  StdVector<std::string> list = ParseGroup(loc_id, group);

  return std::find(list.begin(), list.end(), child) != list.end();
}

//-----------------------------------------------------------------------------
int GetDataSetRank(hid_t loc_id, const std::string& dsetName)
{
  int val = 0;
  if (H5LTget_dataset_ndims(loc_id, dsetName.c_str(), &val) < 0)
    throw Exception("cannot get rank of data set '" + dsetName + "'");
  
  return val;
}


template <>
void ReadAttribute<int>(hid_t loc_id, const std::string& objName, const std::string& attrName, int& data)
{
  if (H5LTget_attribute_int(loc_id, objName.c_str(), attrName.c_str(), &data) < 0)
    throw Exception("cannot read int attribute " + objName + "/" + attrName);
}

//-----------------------------------------------------------------------------
template <>
void ReadAttribute<unsigned int>(hid_t loc_id, const std::string& objName, const std::string& attrName, unsigned int& data)
{
  if (H5LTget_attribute_uint(loc_id, objName.c_str(), attrName.c_str(), &data) < 0)
    throw Exception("cannot read uint attribute " + objName + "/" + attrName);
}

//-----------------------------------------------------------------------------
template <>
void ReadAttribute<double>(hid_t loc_id, const std::string& objName, const std::string& attrName, double& data)
{
  if (H5LTget_attribute_double(loc_id, objName.c_str(), attrName.c_str(), &data) < 0)
    throw Exception("cannot read double attribute " + objName + "/" + attrName);
}

//-----------------------------------------------------------------------------
template <>
void ReadAttribute<std::string>(hid_t loc_id, const std::string& objName, const std::string& attrName, std::string& data)
{
  hid_t obj_id = H5Oopen(loc_id, objName.c_str(), H5P_DEFAULT);
  if (obj_id < 0)
    throw Exception("ReadAttribute: cannot open object " + objName);

  hid_t attr_id = H5Aopen(obj_id, attrName.c_str(), H5P_DEFAULT);
  if (attr_id < 0) {
    H5Oclose(obj_id);
    throw Exception("ReadAttribute: cannot open attribute " + objName + "/" + attrName);
  }

  hid_t file_type = H5Aget_type(attr_id);
  htri_t is_vl = H5Tis_variable_str(file_type);
  H5Tclose(file_type);

  herr_t err = 0;
  if (is_vl > 0) {
    // Variable-length string (written e.g. by H5CPP)
    char* vl_buf = nullptr;
    hid_t file_type2 = H5Aget_type(attr_id);
    hid_t mem_type = H5Tcopy(file_type2);
    H5Tclose(file_type2);
    err = H5Aread(attr_id, mem_type, &vl_buf);
    H5Tclose(mem_type);
    if (err >= 0 && vl_buf != nullptr)
      data = std::string(vl_buf);
    H5free_memory(vl_buf);
  } else {
    // Fixed-length string (written by H5LTset_attribute_string)
    hid_t file_type2 = H5Aget_type(attr_id);
    size_t size = H5Tget_size(file_type2);
    H5Tclose(file_type2);

    hid_t mem_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(mem_type, size);
    H5Tset_strpad(mem_type, H5T_STR_NULLTERM);

    std::vector<char> buf(size + 1, '\0');
    err = H5Aread(attr_id, mem_type, buf.data());
    H5Tclose(mem_type);
    if (err >= 0)
      data = std::string(buf.data());
  }

  H5Aclose(attr_id);
  H5Oclose(obj_id);

  if (err < 0)
    throw Exception("cannot obtain string attribute value for " + objName + "/" + attrName);
}

//-----------------------------------------------------------------------------
StdVector<unsigned int> GetArrayDims(hid_t loc, const std::string& name)
{
  hid_t set = OpenDataSet(loc, name);
  hid_t space_id = H5Dget_space(set);

  if (H5Sis_simple(space_id) <= 0)
    throw Exception("no simple data space " + name);

  const int ndims = H5Sget_simple_extent_ndims(space_id);

  StdVector<hsize_t> dims(ndims);
  // nullptr means no max dims
  if (H5Sget_simple_extent_dims(space_id, dims.GetPointer(), nullptr) != ndims)
    throw Exception("read dimensions not as expected for " + name);

  H5Sclose(space_id);
  H5Dclose(set);

  StdVector<unsigned int> ret(ndims);
  for (int i = 0; i < ndims; i++)
    ret[i] = (unsigned int) dims[i]; // dims is long long, so we cannot directly write to reg.GetPointer()
  
  return ret;
}

//-----------------------------------------------------------------------------
unsigned int GetNumberOfEntries(hid_t id, const std::string& name)
{
  hid_t set = OpenDataSet(id, name);
  hid_t dspace = H5Dget_space(set);

  if (H5Sis_simple(dspace) <= 0)
    throw Exception("no simple data space " + name);
  
  hssize_t npoints = H5Sget_simple_extent_npoints(dspace);

  H5Sclose(dspace);
  H5Dclose(set);

  if (npoints < 0)
    throw Exception("cannot get number of elements for " + name);
  
  return (unsigned int) npoints;
}

//-----------------------------------------------------------------------------
template <typename T>
void ReadDataSet(hid_t loc, const std::string& name, T* out)
{
  assert(out != nullptr);
  hid_t dset = H5Dopen2(loc, name.c_str(), H5P_DEFAULT);
  if (dset < 0)
    EXCEPTION("cannot open dataset " << name << " for reading: " << dset);
  herr_t err = H5Dread(dset, NativeType<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, out);
  H5Dclose(dset);
  if (err < 0)
    EXCEPTION("error reading dataset " << name << ": " << err);
}

template void ReadDataSet<double>(hid_t, const std::string&, double*);
template void ReadDataSet<int>(hid_t, const std::string&, int*);
template void ReadDataSet<unsigned int>(hid_t, const std::string&, unsigned int*);

/** the string pointer is used by TYPE ReadDataSet<std::string>() */
template <>
void ReadDataSet<std::string>(hid_t loc, const std::string& name, std::string* out)
{
  assert(out != nullptr);
  char* buffer = nullptr;
  if (H5LTread_dataset_string(loc, name.c_str(), reinterpret_cast<char*>(&buffer)) < 0)
    throw Exception("cannot read string dataset " + name);
  
  *out = std::string(buffer); // convert to string and copy
  free(buffer);
}

template <typename TYPE>
TYPE ReadDataSet(hid_t loc, const std::string& name)
{
  assert(GetNumberOfEntries(loc, name) == 1);
  TYPE val;
  ReadDataSet(loc, name, &val);
  return val;
}

//-----------------------------------------------------------------------------
hid_t GetMultiStepGroup(hid_t root, unsigned int msStep, bool isHistory)
{
  std::string path = (isHistory ? "/Results/History/MultiStep_" : "/Results/Mesh/MultiStep_") + std::to_string(msStep);
  return OpenGroup(root, path);
}

//-----------------------------------------------------------------------------
hid_t GetStepGroup(hid_t root, unsigned int msStep, unsigned int stepNum)
{
  std::string path = "/Results/Mesh/MultiStep_" + std::to_string(msStep) + "/Step_" + std::to_string(stepNum);
  return OpenGroup(root, path);
}


template <typename TYPE>
void ReadArray(hid_t loc, const std::string& name, StdVector<TYPE>& data)
{
  data.Resize(GetNumberOfEntries(loc, name));
  ReadDataSet(loc, name, data.GetPointer());
}

template void ReadArray<int>(hid_t loc, const std::string& name, StdVector<int>& data);
template void ReadArray<unsigned int>(hid_t loc, const std::string& name, StdVector<unsigned int>& data);
template void ReadArray<double>(hid_t loc, const std::string& name, StdVector<double>& data);

template <>
void ReadArray<std::string>(hid_t loc, const std::string& name, StdVector<std::string>& data)
{
  const unsigned int n = GetNumberOfEntries(loc, name);
  // We use RAII instead of char** buffer = new char*[n] as suggested in h5 documentation
  StdVector<char*> buffer(n);
  if (H5LTread_dataset_string(loc, name.c_str(), reinterpret_cast<char*>(buffer.GetPointer())) < 0)
    throw Exception("cannot read string array dataset " + name);

  data.Resize(n);
  for (unsigned int i = 0; i < n; i++)
  {
    data[i] = std::string(buffer[i]);
    free(buffer[i]);
  }
}





template unsigned int ReadDataSet<unsigned int>(hid_t loc, const std::string& name);
template std::string ReadDataSet<std::string>(hid_t loc, const std::string& name);

// =============================================================================
// WRITE helpers
// =============================================================================

hid_t CreateGroup(hid_t loc, const std::string& name, bool use_existing)
{
  if(H5Lexists(loc, name.c_str(), H5P_DEFAULT) > 0)
  {
    if(!use_existing)
      EXCEPTION("Group '" << name << "' already exists");
    return OpenGroup(loc, name);
  }
  else
  {
    // we use the intermediate group creation property to create all groups in the path if necessary
    // H5P_LINK_CREATE with intermediate group creation so callers can pass deep paths
    hid_t lcpl = H5Pcreate(H5P_LINK_CREATE);
    H5Pset_create_intermediate_group(lcpl, 1);
    hid_t ret = H5Gcreate2(loc, name.c_str(), lcpl, H5P_DEFAULT, H5P_DEFAULT);
    H5Pclose(lcpl);
    if (ret < 0)
      EXCEPTION("cannot create group '" << name << "': " << ret);
    return ret;
  }
}

//-----------------------------------------------------------------------------
template <>
void WriteAttribute<int>(hid_t loc_id, const std::string& objName, const std::string& attrName, const int& data)
{
  if (H5LTset_attribute_int(loc_id, objName.c_str(), attrName.c_str(), &data, 1) < 0)
    EXCEPTION("cannot write int attribute " << objName << "/" << attrName);
}

//-----------------------------------------------------------------------------
template <>
void WriteAttribute<unsigned int>(hid_t loc_id, const std::string& objName, const std::string& attrName, const unsigned int& data)
{
  if (H5LTset_attribute_uint(loc_id, objName.c_str(), attrName.c_str(), &data, 1) < 0)
    EXCEPTION("cannot write uint attribute " << objName << "/" << attrName);
}

template <>
void WriteAttribute<bool>(hid_t loc_id, const std::string& objName, const std::string& attrName, const bool& data)
{
  WriteAttribute<int>(loc_id, objName, attrName, (int) data);
}


template <>
void WriteAttribute<double>(hid_t loc_id, const std::string& objName, const std::string& attrName, const double& data)
{
  if (H5LTset_attribute_double(loc_id, objName.c_str(), attrName.c_str(), &data, 1) < 0)
    EXCEPTION("cannot write double attribute " << objName << "/" << attrName);
}

//-----------------------------------------------------------------------------
template <>
void WriteAttribute<std::string>(
  hid_t loc_id, const std::string& objName, const std::string& attrName, const std::string& data)
{
  // Write as VL string (same encoding as hdf5io.cc / H5CPP) so CFSReader can read it
  // H5LTset_attribute_string would crash CFSReader but cfs's ReadAttribute() reads both
  hid_t strtype = H5Tcopy(H5T_C_S1);
  H5Tset_size(strtype, H5T_VARIABLE);

  hid_t obj_id = H5Oopen(loc_id, objName.c_str(), H5P_DEFAULT);
  if (obj_id < 0) {
    H5Tclose(strtype);
    EXCEPTION("cannot write string attribute: cannot open object " << objName);
  }

  if (H5Aexists(obj_id, attrName.c_str()) > 0)
    H5Adelete(obj_id, attrName.c_str());

  hid_t space = H5Screate(H5S_SCALAR);
  hid_t attr_id = H5Acreate2(obj_id, attrName.c_str(), strtype, space, H5P_DEFAULT, H5P_DEFAULT);
  H5Sclose(space);
  if (attr_id < 0) {
    H5Oclose(obj_id);
    H5Tclose(strtype);
    EXCEPTION("cannot create string attribute " << objName << "/" << attrName);
  }

  const char* ptr = data.c_str();
  herr_t err = H5Awrite(attr_id, strtype, &ptr);
  H5Aclose(attr_id);
  H5Oclose(obj_id);
  H5Tclose(strtype);

  if (err < 0)
    EXCEPTION("cannot write string attribute " << objName << "/" << attrName);
}

//-----------------------------------------------------------------------------

// Common core for WriteDataSet1D / WriteDataSet2D
template <typename T>
void WriteDataSetND(hid_t loc, const std::string& name,
                           int rank, hsize_t* dims, hsize_t totalElems,
                           const T* data, int compressionLevel)
{
  assert(data != nullptr);
  assert(compressionLevel >=0 && compressionLevel <= 9);
  compressionLevel = totalElems * sizeof(T) < HDF5_COMPRESSION_THRESHOLD ? 0 : compressionLevel; // disable compression if data is small, otherwise the overhead may be larger than the gain
  hid_t space = H5Screate_simple(rank, dims, nullptr);
  hid_t dcpl = H5P_DEFAULT;
  if (compressionLevel > 0)
  {
    dcpl = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl, rank, dims);
    H5Pset_shuffle(dcpl);
    H5Pset_deflate(dcpl, compressionLevel); // gzip level 1-9 - anything > 1 usually too slow
  }
  hid_t dset = H5Dcreate2(loc, name.c_str(), NativeType<T>(), space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
  H5Sclose(space);
  if (dcpl != H5P_DEFAULT)
    H5Pclose(dcpl);
  if (dset < 0)
    throw Exception("cannot create dataset '" + name + "'");
  herr_t err = H5Dwrite(dset, NativeType<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  H5Dclose(dset);
  if (err < 0)
    throw Exception("cannot write dataset '" + name + "'");
}

template void WriteDataSetND<double>(hid_t, const std::string&, int, hsize_t*, hsize_t, const double*, int);
template void WriteDataSetND<int>(hid_t, const std::string&, int, hsize_t*, hsize_t, const int*, int);
template void WriteDataSetND<unsigned int>(hid_t, const std::string&, int, hsize_t*, hsize_t, const unsigned int*, int);

// String specialization remains HL - we ignore compression!
void WriteSingleDataSet(hid_t loc, const std::string& name, const std::string& data)
{
  StdVector<std::string> array{data};
  WriteStringArray(loc, name, array);
}

//-----------------------------------------------------------------------------
void WriteStringArray(hid_t loc, const std::string& name, const StdVector<std::string>& data)
{
  hid_t strtype = H5Tcopy(H5T_C_S1);
  H5Tset_size(strtype, H5T_VARIABLE);

  hsize_t dims[1] = {data.size()};
  hid_t space = H5Screate_simple(1, dims, nullptr);
  hid_t dset =  H5Dcreate2(loc, name.c_str(), strtype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dset < 0)
  {
    H5Sclose(space);
    H5Tclose(strtype);
    throw Exception("cannot create string array dataset " + name);
  }

  StdVector<const char*> ptrs(data.size());
  for (size_t i = 0; i < data.size(); i++)
    ptrs[i] = data[i].c_str();

  if (H5Dwrite(dset, strtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, ptrs.GetPointer()) < 0)
  {
    H5Dclose(dset);
    H5Sclose(space);
    H5Tclose(strtype);
    throw Exception("cannot write string array dataset " + name);
  }

  H5Dclose(dset);
  H5Sclose(space);
  H5Tclose(strtype);
}

// =============================================================================
// EXTENSIBLE dataset helpers (chunked, H5S_UNLIMITED — for growing arrays)
// =============================================================================

template <typename T>
hid_t CreateExtendable(hid_t loc, const std::string& name,
                       const StdVector<hsize_t>& initDims,
                       const StdVector<hsize_t>& maxDims,
                       const StdVector<hsize_t>& chunkDims)
{
  const int rank = static_cast<int>(initDims.size());
  assert(!chunkDims.Contains(0));

  hid_t space = H5Screate_simple(rank, initDims.GetPointer(), maxDims.GetPointer());
  hid_t dcpl  = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(dcpl, rank, chunkDims.GetPointer());
  hid_t dset  = H5Dcreate2(loc, name.c_str(), NativeType<T>(),
                            space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
  H5Pclose(dcpl);
  H5Sclose(space);
  return dset;
}

/** extend 1D and 2D arrays. 
 * @param currRows the old size of the dataset. Elements in 1D and (rows, numCols) for 2D
 * @param count number of entried to be added. For 2D #rows, numCols (which must always be the same) */
template <typename T>
void ExtendAndWriteSlab(hid_t dset, StdVector<hsize_t> currRows, StdVector<hsize_t> count, const T* data)
{
  const int rank = static_cast<int>(currRows.size());
  StdVector<hsize_t> newDims = currRows;   // fixed dims stay unchanged
  newDims[0] += count[0];                    // only dim-0 (rows) grows
  StdVector<hsize_t> offset(rank);
  offset.Init(0);
  offset[0] = currRows[0];

  H5Dextend(dset, newDims.GetPointer());
  hid_t fileSpace = H5Dget_space(dset);
  H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, offset.GetPointer(), nullptr, count.GetPointer(), nullptr);
  hid_t memSpace = H5Screate_simple(rank, count.GetPointer(), nullptr);
  H5Dwrite(dset, NativeType<T>(), memSpace, fileSpace, H5P_DEFAULT, data);
  H5Sclose(memSpace);
  H5Sclose(fileSpace);
}


/** Create or append to an extendable 1D dataset. 
 * The dataset is always searched by name - kept open datasets and handles are not implemented 
 * We use auto-chunks of HDF5_GROWING_CHUNK_SIZE KB */
template <typename T>
void WriteGrowingDataSet1D(hid_t loc, const std::string& name, const T* data, hsize_t count)
{
  bool exists = H5Lexists(loc, name.c_str(), H5P_DEFAULT) > 0;
  hid_t dset = -1;
  if (exists)
  {
    assert(count > 0 && data != nullptr); // if dataset already exists, we need to add something
    dset = H5Dopen2(loc, name.c_str(), H5P_DEFAULT);
  }
  else
  {
    dset = CreateExtendable<T>(loc, name, {0}, {H5S_UNLIMITED}, {HDF5_GROWING_CHUNK_SIZE/sizeof(T)});
  }
  if (dset < 0)
    throw Exception("cannot open/create 1D dataset '" + name + "'");

  if (data != nullptr && count > 0) {
    hid_t sp = H5Dget_space(dset);
    hsize_t oldSize = 0; // dummy initial value, will be overwritten
    H5Sget_simple_extent_dims(sp, &oldSize, nullptr);
    H5Sclose(sp);
    ExtendAndWriteSlab<T>(dset, {oldSize}, {count}, data);
  }
  H5Dclose(dset);
}

template void WriteGrowingDataSet1D<double>(hid_t, const std::string&, const double*, hsize_t);
template void WriteGrowingDataSet1D<int>(hid_t, const std::string&, const int*, hsize_t);
template void WriteGrowingDataSet1D<unsigned int>(hid_t, const std::string&, const unsigned int*, hsize_t);

/** @see WriteGrowingDataSet1D 
 * @param numCols must always be the same, usually 2 or 3 */
template <typename T>
void WriteGrowingDataSet2D(hid_t loc, const std::string& name, hsize_t newRows, hsize_t numCols, const T* data)
{
  assert(data != nullptr && newRows > 0 && numCols > 0);
  hid_t dset = -1;
  if (H5Lexists(loc, name.c_str(), H5P_DEFAULT) > 0)
    dset = H5Dopen2(loc, name.c_str(), H5P_DEFAULT);
  else
    dset = CreateExtendable<T>(loc, name, {0, numCols}, {H5S_UNLIMITED, numCols}, {HDF5_GROWING_CHUNK_SIZE/(numCols*sizeof(T)), numCols});
  if (dset < 0)
    throw Exception("cannot open/create growing 2D dataset '" + name + "'");

  hid_t sp = H5Dget_space(dset);
  hsize_t tmp[2];
  H5Sget_simple_extent_dims(sp, tmp, nullptr);
  H5Sclose(sp);
  hsize_t currRows = tmp[0];
  ExtendAndWriteSlab<T>(dset, {currRows, numCols}, {newRows, numCols}, data);
  H5Dclose(dset);
}

template void WriteGrowingDataSet2D<double>(hid_t, const std::string&, hsize_t, hsize_t, const double*);
template void WriteGrowingDataSet2D<int>(hid_t, const std::string&, hsize_t, hsize_t, const int*);
template void WriteGrowingDataSet2D<unsigned int>(hid_t, const std::string&, hsize_t, hsize_t, const unsigned int*);

//-----------------------------------------------------------------------------
// Entity group name mapping: EntityUnknownType -> HDF5 group name string

Enum<ResultInfo::EntityUnknownType> entityGroupNameEnum;

static void SetupEntityGroupNameEnum()
{
  entityGroupNameEnum.SetName("HDF5 entity group names");
  entityGroupNameEnum.Add(ResultInfo::NODE,           "Nodes");
  entityGroupNameEnum.Add(ResultInfo::EDGE,           "Edges");
  entityGroupNameEnum.Add(ResultInfo::FACE,           "Faces");
  entityGroupNameEnum.Add(ResultInfo::ELEMENT,        "Elements");
  entityGroupNameEnum.Add(ResultInfo::SURF_ELEM,      "Elements",      false); // same string as ELEMENT
  entityGroupNameEnum.Add(ResultInfo::REGION,         "Regions");
  entityGroupNameEnum.Add(ResultInfo::REGION_AVERAGE, "Regions",       false); // same string as REGION
  entityGroupNameEnum.Add(ResultInfo::SURF_REGION,    "ElementGroup");
  entityGroupNameEnum.Add(ResultInfo::NODELIST,       "NodeGroup");
  entityGroupNameEnum.Add(ResultInfo::COIL,           "Coils");
  entityGroupNameEnum.Add(ResultInfo::FREE,           "Unknowns");
}

// trick to enforce SetupEntityGroupNameEnum() is called - we have no constructor
static const bool entityGroupNameEnum_initialized_ = (SetupEntityGroupNameEnum(), true);

//-----------------------------------------------------------------------------
// Legacy integer code mappings (enum->int, no Enum<T> possible)

int MapEntityTypeToInt(ResultInfo::EntityUnknownType t)
{
  switch (t) {
    case ResultInfo::NODE:           return 1;
    case ResultInfo::EDGE:           return 2;
    case ResultInfo::FACE:           return 3;
    case ResultInfo::ELEMENT:        return 4;
    case ResultInfo::SURF_ELEM:      return 4;  // plug-in requires same code as ELEMENT
    case ResultInfo::REGION:         return 7;
    case ResultInfo::REGION_AVERAGE: return 7;  // same code as REGION
    case ResultInfo::SURF_REGION:    return 8;
    case ResultInfo::NODELIST:       return 9;
    case ResultInfo::COIL:           return 10;
    case ResultInfo::FREE:           return 11;
  }
  return 0;
}

ResultInfo::EntityUnknownType MapEntityTypeFromInt(int t)
{
  switch (t) {
    case 1:  return ResultInfo::NODE;
    case 2:  return ResultInfo::EDGE;
    case 3:  return ResultInfo::FACE;
    case 4:  return ResultInfo::ELEMENT;
    case 5:  return ResultInfo::SURF_ELEM;
    case 7:  return ResultInfo::REGION;
    case 8:  return ResultInfo::SURF_REGION;
    case 9:  return ResultInfo::NODELIST;
    case 10: return ResultInfo::COIL;
    case 11: return ResultInfo::FREE;
  }
  return ResultInfo::NODE;
}

int MapEntryTypeToInt(ResultInfo::EntryType t)
{
  switch (t) {
    case ResultInfo::UNKNOWN: return 0;
    case ResultInfo::SCALAR:  return 1;
    case ResultInfo::VECTOR:  return 3;
    case ResultInfo::TENSOR:  return 6;
    case ResultInfo::STRING:  return 32;
  }
  return 0;
}

ResultInfo::EntryType MapEntryTypeFromInt(int t)
{
  switch (t) {
    case 0:  return ResultInfo::UNKNOWN;
    case 1:  return ResultInfo::SCALAR;
    case 3:  return ResultInfo::VECTOR;
    case 6:  return ResultInfo::TENSOR;
    case 32: return ResultInfo::STRING;
  }
  return ResultInfo::UNKNOWN;
}

} // end of namespace CoupledField
