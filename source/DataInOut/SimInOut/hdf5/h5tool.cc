#include <stddef.h>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "boost/algorithm/string/trim.hpp"
#include "boost/tokenizer.hpp"

//#include <General/environment.hh>

#include "General/defs.hh"
#include "General/exception.hh"
#include "Utils/StdVector.hh"
#include "hdf5io.hh"

using namespace CoupledField;

template <typename TYPE>
CoupledField::Integer CreateAttribute(H5::H5File& file,
                        std::string path,
                        std::string attr_name,
                        TYPE attr_value,
                        UInt mode)
{
  H5::H5Object *object = NULL;
  H5::DataSet dataset;
  H5::Group group;
  H5::Attribute attrib;
  H5::DataSpace space;

  try
  {
    group = file.openGroup(path);
    object = &group;
  }
  catch (H5::Exception&) { }

  if(!object)
  {
    try
    {
      dataset = file.openDataSet(path);
      object = &dataset;
    }
    catch (H5::Exception&)
    {
      std::cerr << "CreateAttribute: Failed to open dataset or group!" << std::endl;
      return 1;
    }
  }
  

  UInt numAttrs = object->getNumAttrs();
  for(UInt i=0; i<numAttrs; i++)
  {
    H5::Attribute attr = object->openAttribute(i);
    
    //std::cout << "PATH: " << path << "/" << attr.getName() << std::endl;
    if(attr.getName() == attr_name)
    {
      switch(mode) {
      case 0: // mode fail
        return 1;
      case 1: // mode ignore
        return 0;
      case 2: // mode overwrite
        try {
          attr.close();
          object->removeAttr(attr_name);
        }
        catch (H5::Exception&)
        {
          std::cerr << "CreateAttribute: Removal of old attribute failed!"
                    << std::endl;
          return 1;
        }
        break;
        
      default:
        return 1;
      }
      
      break;
    }
    
    attr.close();
  }
  

  try
  {
    H5IO::WriteAttribute( *object, attr_name, attr_value );
    
    object->close();
    
    return 0;
  }
  catch (H5::Exception& h5ex)
  {
    std::cerr << "CreateAttribute: " << h5ex.getCDetailMsg() << std::endl;
    return 1;
  }
}


template <typename TYPE>
CoupledField::Integer ReadAttribute(H5::H5File& file, std::string path,
    std::string attr_name, TYPE& value)
{
  H5::H5Object *object = NULL;
  H5::DataSet dataset;
  H5::Group group;

  try
  {
    group = file.openGroup(path);
    object = &group;
  }
  catch (H5::Exception& h5ex) { }

  if(!object)
  {
    try
    {
      dataset = file.openDataSet(path);
      object = &dataset;
    }
    catch (H5::Exception& h5ex)
    {
      std::cerr << "ReadAttribute: Failed to open dataset or group!\n";
      return 1;
    }
  }

  try {
    H5IO::ReadAttribute(*object, attr_name, value);
  } catch (Exception& ex) {
    std::cerr << ex.GetMsg() << std::endl;
    return 1;
  }

  return 0;
}


CoupledField::Integer RenameAttribute(H5::H5File& file,
                        std::string path,
                        std::string old_name,
                        std::string new_name)
{
  H5::H5Object *object = NULL;
  H5::DataSet dataset;
  H5::Group group;

  try
  {
    group = file.openGroup(path);
    object = &group;
  }
  catch (H5::Exception& h5ex) { }

  if(!object)
  {
    try
    {
      dataset = file.openDataSet(path);
      object = &dataset;
    }
    catch (H5::Exception& h5ex)
    {
      std::cerr << "RenameAttribute: Failed to open dataset or group!\n";
      return 1;
    }
  }

  try {
    object->renameAttr(old_name, new_name);
  } catch (H5::AttributeIException& h5aex) {
    std::cerr << "RenameAttribute: Unable to rename attribute '" << old_name
              << "'\n";
    return 1;
  }
  
  return 0;
}


CoupledField::Integer DeleteAttribute(H5::H5File& file,
                        std::string path,
                        std::string attr_name)
{
  H5::H5Object *object = NULL;
  H5::DataSet dataset;
  H5::Group group;

  try
  {
    group = file.openGroup(path);
    object = &group;
  }
  catch (H5::Exception& h5ex) { }

  if(!object)
  {
    try
    {
      dataset = file.openDataSet(path);
      object = &dataset;
    }
    catch (H5::Exception& h5ex)
    {
      std::cerr << "DeleteAttribute: Failed to open dataset or group!\n";
      return 1;
    }
  }

  try {
    object->removeAttr(attr_name);
  } catch (H5::AttributeIException& h5aex) {
    std::cerr << "DeleteAttribute: Unable to delete attribute '" << attr_name
              << "'\n";
    return 1;
  }
  
  return 0;
}


CoupledField::Integer CreateGroup(H5::H5File& file,
                    std::string path,
                    std::string group_name,
                    UInt mode)
{
  H5::Group base_group, group;

  try
  {
    base_group = file.openGroup(path);
  }
  catch (H5::Exception&)
  {
    std::cerr << "CreateGroup: Failed to open base group!" << std::endl;
  }

  UInt numObjs = static_cast<UInt>( base_group.getNumObjs() );
  
  for(UInt i=0; i<numObjs; i++)
  {
    std::string objName = H5IO::GetObjNameByIdx( base_group, i );
    //std::cout << "PATH: " << path << "/" << objName << std::endl;
    
    if(objName == group_name)
    {
      switch(mode) {
      case 0: // mode fail
        return 1;
      case 1: // mode ignore
        return 0;
      case 2: // mode overwrite
        try {
          base_group.unlink(group_name);
        }
        catch (H5::Exception&)
        {
          std::cerr << "CreateGroup: Unlinking of old group failed!"
                    << std::endl;
          return 1;
        }
        break;
        
      default:
        return 1;
      }
      
      break;
    }
  }

  try
  {
    group = base_group.createGroup(group_name);
  }
  catch (H5::Exception&)
  {
    std::cerr << "CreateGroup: Failed to create group!" << std::endl;
    return 1;
  }

  return 0;
}


CoupledField::Integer ListGroup(H5::H5File& file, std::string path) {
  H5::Group base_group;
  
  try {
    base_group = file.openGroup(path);
  } catch (H5::Exception& h5ex) {
    std::cerr << "ListGroup: failed to open group '" << path << "'\n";
    return 1;
  }
  
  UInt numObjs = base_group.getNumObjs();
  
  for (UInt i=0; i < numObjs; ++i) {
    
    if (base_group.getObjTypeByIdx(i) == H5G_GROUP) {
      std::cout << H5IO::GetObjNameByIdx(base_group, i) << std::endl;
    }
  
  }
  
  return 0;
}


CoupledField::Integer CreateStringArrayDataset(H5::H5File& file,
                                 std::string path,
                                 std::string ds_name,
                                 std::string ds_value,
                                 int mode = 0)
{
  H5::Group base_group;

  try
  {
    base_group = file.openGroup(path);
  }
  catch (H5::Exception&)
  {
    std::cerr << "CreateStringArrayDataset: Failed to open base group!" << std::endl;
  }


  UInt numObjs = static_cast<UInt>( base_group.getNumObjs() );
  
  StdVector<std::string> dataset_vec;

  for(UInt i=0; i<numObjs; i++)
  {
    std::string objName;

    std::cout << "Trying to access object nr " << i << std::endl;
    
    try {
      objName = H5IO::GetObjNameByIdx( base_group, i );
    }
    catch (CoupledField::Exception &ex) {
      std::cerr << ex.what() << std::endl;
      return 1;
    }
    
    //std::cout << "PATH: " << path << "/" << objName << std::endl;
    
    if(objName == ds_name)
    {
      switch(mode) {
      case 0: // mode fail
        return 1;
      case 1: // mode ignore
        return 0;
      case 2: // mode overwrite
        try {
          base_group.unlink(ds_name);
        }
        catch (H5::Exception&)
        {
          std::cerr << "CreateStringArrayDataset: Failed unlink old dataset"
                    << " while trying to overwrite it!" << std::endl;
          return 1;
        }
        break;
      case 3: // mode append
      case 4: // mode insert
        try {
          H5IO::ReadArray( base_group, ds_name, dataset_vec );
          base_group.unlink(ds_name);
        }
        catch (H5::Exception&)
        {
          std::cerr << "CreateStringArrayDataset: Error while trying to read"
                    << " old dataset for appending!" << std::endl;
          return 1;
        }
        
        break;
        
      default:
        return 1;
      }
      
      break;
    }
  }

  std::string dummy;
  
  typedef boost::tokenizer<char_separator<char> > Tok;
  boost::char_separator<char> sep(";| ");
  Tok t(ds_value, sep);
  Tok::iterator it, end;
  it = t.begin();
  end = t.end();
  
  for( ; it != end; it++)
  {
    dummy = *it;
    boost::algorithm::trim(dummy);
    dataset_vec.Push_back(dummy);
  }

  if(mode == 4) {
    std::set< std::string > strset;
    std::set< std::string >::const_iterator it;

    strset.insert(dataset_vec.Begin(), dataset_vec.End());

    dataset_vec.Resize(strset.size());

    UInt i = 0;
    UInt n = strset.size();
    it = strset.begin();

    for( ; i<n; i++, it++)
      dataset_vec[i] = *it;
  }
  
  std::cout << dataset_vec << std::endl;
  
  H5IO::Write1DArray( base_group, ds_name, dataset_vec.GetSize(), &dataset_vec[0] );

  return 0;
}

CoupledField::Integer CreateStringDataset(H5::H5File& file,
                            std::string path,
                            std::string ds_name,
                            std::string ds_value,
			    UInt mode)
{
  H5::Group base_group;

  try
  {
    base_group = file.openGroup(path);
  }
  catch (H5::Exception&)
  {
    std::cerr << "CreateStringDataset: Failed to open base group!" << std::endl;
  }


  UInt numObjs = static_cast<UInt>( base_group.getNumObjs() );
  
  for(UInt i=0; i<numObjs; i++)
  {
    std::string objName = H5IO::GetObjNameByIdx( base_group, i );
    //    std::cout << "PATH: " << path << "/" << objName << std::endl;
    if(objName == ds_name) {

      switch(mode) {
      case 0: // mode fail
        return 1;
      case 1: // mode ignore
        return 0;
      case 2: // mode overwrite
        try {
          base_group.unlink(ds_name);
        }
        catch (H5::Exception&)
        {
          std::cerr << "CreateStringDataset: Failed unlink old dataset"
                    << " while trying to overwrite it!" << std::endl;
          return 1;
        }
        break;
      default:
        return 1;
      }
      
      break;
    }
  }

#if 0
  try
  {
    H5::DataSpace space;
    
    H5::DataSet dataset;
    dataset = base_group.createDataSet(ds_name,
                                       H5::StrType(H5::PredType::C_S1,
                                                   ds_value.length()),
                                       space);
    
    dataset.write(ds_value.c_str(), H5::StrType(H5::PredType::C_S1, ds_value.length()));
    
    base_group.close();
  }
  catch (H5::Exception&)
  {
    return 1;
  }
#endif

  std::vector<std::string> dataset_vec;
  
  boost::algorithm::trim(ds_value);

  dataset_vec.push_back(ds_value);
  
  H5IO::Write1DArray( base_group, ds_name, dataset_vec.size(), &dataset_vec[0] );

  return 0;
}

CoupledField::Integer ExistsDataset(H5::H5File& file,
                      std::string path,
                      std::string ds_name)
{
  H5::Group base_group;

  try
  {
    base_group = file.openGroup(path);
  }
  catch (H5::Exception&)
  {
    // std::cerr << "ExistsDataset: Failed to open base group!" << std::endl;
    return 1;
  }


  UInt numObjs = static_cast<UInt>( base_group.getNumObjs() );
  
  for(UInt i=0; i<numObjs; i++)
  {
    std::string objName = H5IO::GetObjNameByIdx( base_group, i );
    //    std::cout << "PATH: " << path << "/" << objName << std::endl;
    if(objName == ds_name)
    {
      H5::DataSet ds;

      try
      {
        ds = base_group.openDataSet(ds_name);
      }
      catch (H5::Exception&)
      {
        return 1;
      }
      
      ds.close();
      
      return 0;
    }
  }

  return 1;
}

CoupledField::Integer DeleteObj(H5::H5File& file,
                      std::string path,
                      std::string ds_name)
{
  H5::Group base_group;

  try
  {
    base_group = file.openGroup(path);
  }
  catch (H5::Exception&)
  {
    // std::cerr << "ExistsDataset: Failed to open base group!" << std::endl;
    return 0;
  }


  UInt numObjs = static_cast<UInt>( base_group.getNumObjs() );
  
  for(UInt i=0; i<numObjs; i++)
  {
    std::string objName = H5IO::GetObjNameByIdx( base_group, i );
    //    std::cout << "PATH: " << path << "/" << objName << std::endl;
    if(objName == ds_name) {
      try {
        base_group.unlink(ds_name);
      }
      catch (H5::Exception&)
      {
        std::cerr << "DeleteObj: Unlinking of group/dataset '" << ds_name
                  << "' failed" << std::endl;
        return 1;
      }
      break;
    }
  }

  return 0;
}

CoupledField::Integer MoveObj(H5::H5File& file,
                  std::string path,
                  std::string name_old,
                  std::string name_new)
{
  H5::Group base_group;
  
  try {
    base_group = file.openGroup(path);
  } catch (H5::Exception& h5ex) {
    std::cerr << "MoveObj: failed to open base group!\n";
    return 1;
  }
  
  try {
    base_group.move(name_old, name_new);
  } catch (H5::Exception& h5ex) {
    std::cerr << "MoveObj: failed to move '" << name_old << "'\n";
    return 1;
  }
  
  return 0;
}

int main(int argc, char** argv)
{
  std::string op_target;
  std::string op_name;
  std::string file_name;
  std::string path;
  int ret = 0;
  std::map<std::string, UInt> modeMap;

  modeMap["fail"]      = 0;
  modeMap["ignore"]    = 1;
  modeMap["overwrite"] = 2;
  modeMap["append"]    = 3;
  modeMap["insert"]    = 4;


  std::cout << "Operation target:" << std::endl;
  getline(std::cin, op_target);
  std::cout << "Operation name: " << std::endl;
  getline(std::cin, op_name);
  std::cout << "File name: " << std::endl;
  getline(std::cin, file_name);

  // Do not print HDF5 exceptions by default
  H5::Exception::dontPrint();

  H5::H5File *file = NULL;

  if(op_target == "file")
  {
    if(op_name == "create")
    {
      try {
        file = new H5::H5File( file_name, H5F_ACC_TRUNC );  
      }
      catch (H5::Exception&)
      {
        std::cerr << "New HDF5 File '" << file_name
                  << "' could not be openend!" << std::endl;
        return 1;
      }
      file->close();
      return ret;
    }
  }

  if(!file)
  {
    try {
    file = new H5::H5File( file_name, H5F_ACC_RDWR );
    }
    catch (H5::Exception&)
    {
      std::cerr << "HDF5 File '" << file_name
                << "' could not be openend!" << std::endl;
      return 1;
    }

  }
  
  std::cout << "Path: " << std::endl;
  getline(std::cin, path);

  if(op_target == "attribute")
  {
    if(op_name == "create")
    {
      std::string datatype;
      std::string attr_name;
      std::string attr_value;
      std::string attr_opmode;

      std::cout << "Data type: " << std::endl;
      getline(std::cin, datatype);
      std::cout << "Attribute name: " << std::endl;
      getline(std::cin, attr_name);
      std::cout << "Operation mode: " << std::endl;
      getline(std::cin, attr_opmode);
      std::cout << "Attribute value: " << std::endl;
      getline(std::cin, attr_value);
      
      if(modeMap.find(attr_opmode) == modeMap.end()) {
        std::cerr << "Operation mode '" << attr_opmode
                  << "' not permitted!" << std::endl;
        file->close();
        return 1;
      }
      
      if(datatype == "uint")
      {
        UInt val;
        std::stringstream str;
        str << attr_value;
        str >> val;

        ret = CreateAttribute(*file,
                              path,
                              attr_name,
                              val,
                              modeMap[attr_opmode]);
      }
      
      if(datatype == "double")
      {
        Double val;
        std::stringstream str;
        str << attr_value;
        str >> val;

        ret = CreateAttribute(*file,
                              path,
                              attr_name,
                              val,
                              modeMap[attr_opmode]);
      }
      
      if(datatype == "string")
      {
        ret = CreateAttribute(*file,
                              path,
                              attr_name,
                              attr_value,
                              modeMap[attr_opmode]);
      }
    }
    else if (op_name == "read") {
      std::string attr_name, datatype;
      
      std::cout << "Attribute name:\n";
      getline(std::cin, attr_name);
      std::cout << "Data type:\n";
      getline(std::cin, datatype);
      
      if (datatype == "double") {
        Double value;
        ret = ReadAttribute(*file, path, attr_name, value);
        if (ret == 0)
          std::cout << value;
      }
      else if (datatype == "string") {
        std::string value;
        ret = ReadAttribute(*file, path, attr_name, value);
        if (ret == 0)
          std::cout << value;
      }
      else if (datatype == "uint") {
        UInt value;
        ret = ReadAttribute(*file, path, attr_name, value);
        if (ret == 0)
          std::cout << value;
      }
    }
    else if (op_name == "rename") {
      std::string old_name, new_name;
      
      std::cout << "Current attribute name:\n";
      getline(std::cin, old_name);
      std::cout << "New attribute name:\n";
      getline(std::cin, new_name);
      
      ret = RenameAttribute(*file, path, old_name, new_name);
    }
    else if (op_name == "delete") {
      std::string attr_name;
      
      std::cout << "Attribute name:\n";
      getline(std::cin, attr_name);
      
      ret = DeleteAttribute(*file, path, attr_name);
    }
  }

  if(op_target == "group")
  {
    if(op_name == "create")
    {
      std::string group_name;
      std::string group_opmode;

      std::cout << "Group name: " << std::endl;
      getline(std::cin, group_name);
      std::cout << "Operation mode: " << std::endl;
      getline(std::cin, group_opmode);

      if(modeMap.find(group_opmode) == modeMap.end()) {
        std::cerr << "Operation mode '" << group_opmode
                  << "' not permitted!" << std::endl;
        file->close();
        return 1;
      }
      
      ret = CreateGroup(*file, path, group_name, modeMap[group_opmode]);
    }
    
    if (op_name == "move") {
      std::string group_old;
      std::string group_new;
      
      std::cout << "Current group name: \n";
      getline(std::cin, group_old);
      std::cout << "New group name: \n";
      getline(std::cin, group_new);
      
      ret = MoveObj(*file, path, group_old, group_new);
    }
    
    if (op_name == "delete") {
      std::string group_name;
      
      std::cout << "Group name: \n";
      getline(std::cin, group_name);
      
      ret = DeleteObj(*file, path, group_name);
    }
    
    if (op_name == "list") {
      ret = ListGroup(*file, path);
    }
  }

  if(op_target == "dataset")
  {
    if(op_name == "create")
    {
      std::string datatype;
      std::string ds_name;
      std::string ds_value;
      std::string ds_opmode;

      std::cout << "Data type: " << std::endl;
      getline(std::cin, datatype);
      std::cout << "Dataset name: " << std::endl;
      getline(std::cin, ds_name);
      std::cout << "Operation mode: " << std::endl;
      getline(std::cin, ds_opmode);

      if(modeMap.find(ds_opmode) == modeMap.end()) {
        std::cerr << "Operation mode '" << ds_opmode
                  << "' not permitted!" << std::endl;
        file->close();
        return 1;
      }
      
      //      std::cout << "Dataset value: " << std::endl;
      //      getline(std::cin, ds_value);

      std::cout << "Dataset Value: " << std::endl;
      
      std::string line;
      std::stringstream sstr;
      
      while(getline(std::cin, line))
      {
        //        std::cout << "getline: " << line << std::endl;
        sstr << line << std::endl;
      }
      ds_value = sstr.str();

      std::cout << "Dataset Value read: " << ds_value << std::endl;
      

      if(datatype == "string_array")
      {
                
        ret = CreateStringArrayDataset(*file,
                                       path,
                                       ds_name,
                                       ds_value,
                                       modeMap[ds_opmode]);
      }

      if(datatype == "string")
      {
        ret = CreateStringDataset(*file,
				  path,
				  ds_name,
				  ds_value,
				  modeMap[ds_opmode]);
      }
    }

    if(op_name == "exists")
    {
      std::string ds_name;

      std::cout << "Dataset name: " << std::endl;
      getline(std::cin, ds_name);
      
      ret = ExistsDataset(*file, path, ds_name);
    }

    if(op_name == "delete")
    {
      std::string ds_name;

      std::cout << "Dataset name: " << std::endl;
      getline(std::cin, ds_name);

      ret = DeleteObj(*file, path, ds_name);
    }
    
    if (op_name == "move") {
      std::string ds_old;
      std::string ds_new;
      
      std::cout << "Current dataset name: \n";
      getline(std::cin, ds_old);
      std::cout << "New dataset name: \n";
      getline(std::cin, ds_new);
      
      ret = MoveObj(*file, path, ds_old, ds_new);
    }
  }


  file->close();

  return ret;
}

