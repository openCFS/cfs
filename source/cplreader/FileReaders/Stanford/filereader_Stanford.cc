
#include <string>
#include <iostream>
#include <fstream>
#include "stdio.h"
#include <iomanip>

#ifdef MpCCI
#include "mpcci.h"
#endif

#include "filereader_Stanford.hh"


namespace CoupledField
{

    FileReader_Stanford::FileReader_Stanford(const std::string& name,
                                             const int dim,
                                             const int numFiles) :
      FileReader(name, dim, numFiles, 1) {

//       std::vector<Realtype> NodalCoords;
//       ReadNodalCoords(NodalCoords, 0);

//       std::vector<int> Topology;
//       std::vector<int> numNodesPerElem;
//       std::vector<int> elemTypes;
//       ReadTopology(Topology, numNodesPerElem,
//                    elemTypes,0);

//       std::vector<double> flowdata;
//       ReadNodalValues(flowdata,0,0);
    }

    FileReader_Stanford::~FileReader_Stanford()
    {
    }

    void FileReader_Stanford::Init()
    {
        std::string filename;
        char buf[128];
        int dummy;
        
        filename = basename_;
        sprintf(buf, "%02i", 1);
        filename+= buf;
        filename+= ".coord";

        std::cout << "filename:" << filename << std::endl;

        infile.clear();
        infile.open(filename.c_str());  
        if (!infile) {
            std::cerr << "ERROR(" << __FILE__ << " " << __LINE__
                      << ") Can't open " << filename << std::endl;
            exit(1);
        }

        infile >> dummy;
        infile >> numPartitions_;
        infile >> dummy;
        numFiles_ = dummy < numFiles_ ? dummy : numFiles_;
        if (settings.GetInt("numsteps"))
        {
          const UInt tmp = (UInt) settings.GetInt("numsteps");
          /* only take argument if tmo does not exceed the maximal number of timesteps possible */
          if (tmp < numFiles_)
          {
            numFiles_ = tmp;
          }
        }

        infile.close();
        infile.clear();

        filename = basename_;
        sprintf(buf, "%02i", 1);
        filename+= buf;
        filename+= "_";
        sprintf(buf, "%04i", 1);
        filename+= buf;
        filename+= ".dat";

        infile.clear();
        infile.open(filename.c_str());  
        if (!infile) {
            std::cerr << "ERROR(" << __FILE__ << " " << __LINE__
                      << ") Can't open " << filename << std::endl;
            exit(1);
        }

        infile >> numResults_;

        infile.close();
        infile.clear();

        std::cout << "Name: " << name_ << " Dim: " << dim_
                  << " numfiles: " << numFiles_
                  << " numPartitions: " << numPartitions_
                  << " numResults: " << numResults_
                  <<std::endl;

        elsize_.resize(numPartitions_);
        MpCCInodes_.resize(numPartitions_);
        MpCCIelems_.resize(numPartitions_);

        for(int i=0; i<numPartitions_; i++)
        {
            filename = basename_;
            sprintf(buf, "%02i", i+1);
            filename+= buf;
            filename+= ".coord";

            infile.clear();
            infile.open(filename.c_str());  
            if (!infile) {
                std::cerr << "ERROR(" << __FILE__ << " " << __LINE__
                          << ") Can't open " << filename << std::endl;
                exit(1);
            }

            infile >> dummy;
            infile >> dummy;
            infile >> dummy;
            infile >> MpCCInodes_[i];

            infile.close();
            infile.clear();


            filename = basename_;
            sprintf(buf, "%02i", i+1);
            filename+= buf;
            filename+= ".topo";

            infile.open(filename.c_str());  
            if (!infile) {
                std::cerr << "ERROR(" << __FILE__ << " " << __LINE__
                          << ") Can't open " << filename << std::endl;
                exit(1);
            }

//             infile >> dummy;
//             infile >> dummy;
//             infile >> dummy;
            infile >> MpCCIelems_[i];
            infile >> elsize_[i];

            infile.close();
            infile.clear();

            std::cout << "Partition " << (i+1)
                      << " nodes: " << MpCCInodes_[i]
                      << " elems: " << MpCCIelems_[i]
                      << " elsize: " << elsize_[i]
                      <<std::endl;

        }
    }


    void FileReader_Stanford::ReadNodalCoords(std::vector<Realtype> & NODECOORD,
                                             const int partitionIdx)
    {
     #ifdef TRACE
        (*trace) << "entering FileReader_Stanford::ReadNodalCoords" << std::endl;
     #endif

        std::string filename;
        char buf[128];
        int dummy;
        
        filename = basename_;
        sprintf(buf, "%02i", partitionIdx+1);
        filename+= buf;
        filename+= ".coord";

        infile.clear();
        infile.open(filename.c_str());  
        if (!infile) {
            std::cerr << "ERROR(" << __FILE__ << " " << __LINE__
                      << ") Can't open " << filename << std::endl;
            exit(1);
        }

        /* Set pointer to beginning of file: */
        std::string::size_type pos=0;
        infile.seekg(pos, std::ios::beg);

        infile >> dummy; 
        infile >> dummy; 
        infile >> dummy; 
        infile >> dummy; 

        NODECOORD.resize(3*MpCCInodes_[partitionIdx]);

        std::vector<double> max, min;

        if (dim_ == 3) {
          max.resize(3);
          min.resize(3);
        }
        else {
          max.resize(2);
          min.resize(2);
        }

        if ( dim_ == 3 ) {
          infile >> dummy >> NODECOORD[0] >> NODECOORD[1] >> NODECOORD[2];
          min[0] = max[0] = NODECOORD[0];
          min[1] = max[1] = NODECOORD[1];  
          min[2] = max[2] = NODECOORD[2];
        }
        else {
          infile >> dummy >> NODECOORD[0] >> NODECOORD[1];
          NODECOORD[2] = 0.0;
          min[0] = max[0] = NODECOORD[0];
          min[1] = max[1] = NODECOORD[1];  
        }

        for (int i=1; i < MpCCInodes_[partitionIdx]; i++) {
          if ( dim_ == 3 ) {
            for (int j=0; j < 3; j++) {
              infile >> dummy >> NODECOORD[3*i+j];
              
              if(NODECOORD[3*i+j] > max[j])
                max[j] = NODECOORD[3*i+j];
              
              if(NODECOORD[3*i+j] < min[j])
                min[j] = NODECOORD[3*i+j];
            }
          }
          else {
            infile >> dummy;
            for (int j=0; j < 2; j++) {
          
              infile >> NODECOORD[3*i+j];
                            
              if(NODECOORD[3*i+j] > max[j])
                max[j] = NODECOORD[3*i+j];
              
              if(NODECOORD[3*i+j] < min[j])
                min[j] = NODECOORD[3*i+j];
            }
            NODECOORD[3*i+2] = 0.0;
          }
        }

        if ( dim_ == 3) {
          std::cout << "MAX COORD: ( "<< max[0] << ", " << max[1] << ", "<< max[2] << ")"<< std::endl;
          std::cout << "MIN COORD: ( "<< min[0] << ", " << min[1] << ", "<< min[2] << ")"<< std::endl;
        }
        else {
          std::cout << "MAX COORD: ( "<< max[0] << ", " << max[1] << ")"<< std::endl;
          std::cout << "MIN COORD: ( "<< min[0] << ", " << min[1] << ")"<< std::endl;
        }

        infile.close();
    }

    void FileReader_Stanford::ReadTopology(std::vector<int> & TOPOLOGYDATA,
                                           std::vector<int> & numNodesPerElem,
                                           std::vector<int> & elemTypes,
                                           const int partitionIdx)
    {
     #ifdef TRACE
        (*trace) << "entering FileReader_Stanford::ReadTopology" << std::endl;
     #endif

        std::string filename;
        char buf[128];
        int dummy, numNodes, elemType;
        
        filename = basename_;
        sprintf(buf, "%02i", partitionIdx+1);
        filename+= buf;
        filename+= ".topo";

        infile.clear();
        infile.open(filename.c_str());  
        if (!infile) {
            std::cerr << "ERROR(" << __FILE__ << " " << __LINE__
                      << ") Can't open " << filename << std::endl;
            exit(1);
        }

        /* Set pointer to beginning of file: */
        std::string::size_type pos=0;
        infile.seekg(pos, std::ios::beg);

        infile >> dummy >> dummy;

        numNodes = elsize_[partitionIdx];
        switch(numNodes)
        {
	case 2:
	  {
	    elemType = CCI_ELEM_LINE;
	    break;
	  }
        case 3:
          {
            elemType = CCI_ELEM_TRIANGLE;
            break;
          }
        case 4:
          {
            if(dim_ == 3)
              elemType = CCI_ELEM_TETRAHEDRON;
            else
              elemType = CCI_ELEM_QUAD;
            break;
          }
	case 8:
	  {
	    if (dim_ == 3)
	      elemType = CCI_ELEM_HEXAHEDRON;
	    else
	      elemType = CCI_ELEM_QUAD8;
	    break;
          }
        case 5:
	  {
	    elemType = CCI_ELEM_PYRAMID;
	    break;
	  }
	case 6:
	  {
	    if (dim_ == 3)
	      elemType = CCI_ELEM_PRISM;
	    else
	      elemType = CCI_ELEM_TRIANGLE6;
	    break;
	  }
        }

        TOPOLOGYDATA.resize(elsize_[partitionIdx]*MpCCIelems_[partitionIdx]);
        numNodesPerElem.resize(MpCCIelems_[partitionIdx]);
        elemTypes.resize(MpCCIelems_[partitionIdx]);

        for (int i=0; i < MpCCIelems_[partitionIdx]; i++) {
          //element number
          infile >> dummy;
          for (int j=0; j < elsize_[partitionIdx]; j++) {
            infile >> TOPOLOGYDATA[elsize_[partitionIdx]*i+j];
          }
          numNodesPerElem[i] = numNodes;
          elemTypes[i] = elemType;
        }

        infile.close();
    }

    void FileReader_Stanford::ReadNodalValues(std::vector<double> & flowdata,
                                             const int partitionIdx,
                                             const int timeStepIdx)
    {
     #ifdef TRACE
        (*trace) << "entering FileReader_Stanford::ReadNodalValues" << std::endl;
     #endif


        std::string filename;
        char buf[128];
        int dummy;

        
        filename = basename_;
        sprintf(buf, "%02i", partitionIdx+1);
        filename+= buf;
        filename+= "_";
        sprintf(buf, "%04i", timeStepIdx+1);
        filename+= buf;
        filename+= ".dat";

        infile.clear();
        infile.open(filename.c_str());
        if (!infile) {
            std::cerr << "ERROR(" << __FILE__ << " " << __LINE__
                      << ") Can't open " << filename << std::endl;
            exit(1);
        }

        /* Set pointer to beginning of file: */
        std::string::size_type pos=0;
        infile.seekg(pos, std::ios::beg);

        infile >> dummy;

        flowdata.resize(7*MpCCInodes_[partitionIdx]);
  
        double soundSpeedNormalized, src1, src2, src3, totalSrc;

        for (int i=0; i < MpCCInodes_[partitionIdx]; i++) {
          if ( numResults_ == 4 ) {
            infile >> soundSpeedNormalized >> src1 >> src2  >> src3;
            totalSrc = src1 + src2 + src3;
            flowdata[i*7+0] = totalSrc;
            flowdata[i*7+1] = soundSpeedNormalized;
            flowdata[i*7+2] = 0.0;
            flowdata[i*7+3] = 0.0;
            flowdata[i*7+4] = 0.0;
            flowdata[i*7+5] = 0.0;
            flowdata[i*7+6] = 0.0;
            //            std::cout << "Q=" << flowdata[i*7] << "  cn=" << flowdata[i*7+1] << std::endl; 
          }
          else {
            std::cerr << "Number of Results not supported" << std::endl;
            exit(1);
          }
        }

        infile.close();

        
    }

} // end of namespace
