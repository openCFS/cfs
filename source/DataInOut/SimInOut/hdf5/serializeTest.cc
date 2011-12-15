/**
 * play around with boost serialization engine
 *
 */

#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// #include <tree.hpp>
#include "boost/serialization_hdf5/hdf5_oarchive.hpp"
#include "boost/serialization_hdf5/hdf5_iarchive.hpp"
#include "boost/serialization_hdf5/hdf5_attribute.hpp"
// #include <serialization/tree.hpp>

#include "H5Cpp.h"

#include <cstdio> // remove
#include "boost/config.hpp"
#if defined(BOOST_NO_STDC_NAMESPACE)
namespace std{ 
    using ::remove;
}
#endif

#include "boost/variant.hpp"
#include "boost/shared_ptr.hpp"

#include "boost/archive/tmpdir.hpp"
#include "boost/archive/xml_iarchive.hpp"
#include "boost/archive/xml_oarchive.hpp"
//#include <boost/archive/text_iarchive.hpp>
//#include <boost/archive/text_oarchive.hpp>
//#include <boost/archive/binary_iarchive.hpp>
//#include <boost/archive/binary_oarchive.hpp>

#include "boost/serialization/vector.hpp"
#include "boost/serialization/variant.hpp"

int main(int argc, char *argv[]) {   

  // serialize very simple objects
  //std::vector<std::size_t> vUints(1024);
  //vUints[0] = 1; vUints[1] = 10; vUints[2] = 100;
  //std::vector<std::size_t> vUints(5);
  std::vector<double> vUints(5);
  for(std::size_t i=0; i < vUints.size(); ++i)
    vUints[i] = double(i)/double(vUints.size());

  //std::vector<std::size_t> *pvUints = &vUints;
  std::vector<double> *pvUints = &vUints;

  std::string filename(boost::archive::tmpdir());
  filename += "/sTest.xml";
  //filename += "/sTest.txt";
  //filename += "/sTest3.bin";

  std::string empty = "";

  std::ofstream ofs(filename.c_str());
  assert(ofs.good());
  boost::archive::xml_oarchive oa(ofs);
  //boost::archive::text_oarchive oa(ofs);
  //boost::archive::binary_oarchive oa(ofs);

  oa << BOOST_SERIALIZATION_NVP(vUints);
  //oa << BOOST_SERIALIZATION_NVP(vUints);
  oa << BOOST_SERIALIZATION_NVP(pvUints);
  oa << BOOST_SERIALIZATION_NVP(filename);
  
  // close the file
  //oa.close();
  ofs.close();

  std::ifstream ifs(filename.c_str());
  assert(ifs.good());
  boost::archive::xml_iarchive ia(ifs);
  //boost::archive::text_iarchive ia(ifs);
  //boost::archive::binary_iarchive ia(ifs);

 // restore data from the archive
  vUints.clear();
  ia >> BOOST_SERIALIZATION_NVP(vUints);
  ia >> BOOST_SERIALIZATION_NVP(pvUints);

  std::cout << "prime-time: HDF5 serialization code ..." << std::endl;

  // create HDF5 file
  H5::H5File* h5file = new H5::H5File( std::string("sTest.hdf5").c_str(), H5F_ACC_TRUNC );

  boost::archive::hdf5_oarchive ho(h5file);

  // test compression
  ho.set_compression(6);

  double avar = 3.14;
  std::size_t count = 5000;

  std::cout << "count = " << count << "; avar = " << avar << std::endl;

  ho <<  BOOST_SERIALIZATION_NVP(vUints);
  ho << boost::serialization::make_h5attr("vUints-0", "id", &count);
  ho << BOOST_SERIALIZATION_NVP(pvUints);
  ho << BOOST_SERIALIZATION_NVP(filename);
  ho << BOOST_SERIALIZATION_NVP(empty);

  // make a sub-archive
  boost::archive::hdf5_oarchive hoSub(&ho, "subgroup");

  hoSub << BOOST_SERIALIZATION_NVP(avar);
  
  delete h5file;

  //return 0;

  // open the hdf5 file and see what is there
  h5file = new H5::H5File( std::string("sTest.hdf5").c_str(), H5F_ACC_RDWR);

  H5::Group *rootGroup = new H5::Group( h5file->openGroup("/") );

  std::size_t numElems = rootGroup->getNumObjs();
  std::cout << "wrote " << numElems << " entries to hdf5 archive" << std::endl;

  for(std::size_t i = 0; i < numElems; ++i)
    std::cout << rootGroup->getObjnameByIdx(i) << std::endl;
  

  std::cout << "restoring objects from hdf5-file" << std::endl;
  boost::archive::hdf5_iarchive hi(h5file);

  avar = count = 0;
  pvUints = NULL;
  filename = "";
  std::cout << "before load: count = " << count << "; avar = " << avar << "; filename = " << filename << std::endl;

  hi >>  BOOST_SERIALIZATION_NVP(vUints);
  hi >> boost::serialization::make_h5attr("vUints-0", "id", &count);
  hi >>  BOOST_SERIALIZATION_NVP(pvUints);
  hi >> BOOST_SERIALIZATION_NVP(filename);
  hi >> BOOST_SERIALIZATION_NVP(empty);

  std::cout << "pvUints[0] = " << (*pvUints)[0] << std::endl;

  boost::archive::hdf5_iarchive hiSub(&hi, "subgroup");

  hiSub >> BOOST_SERIALIZATION_NVP(avar);

  std::cout << "after load: count = " << count << "; avar = " << avar << "; filename = " << filename << std::endl;

  delete rootGroup;
  delete h5file;

  return 0;
}
