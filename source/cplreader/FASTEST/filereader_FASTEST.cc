
#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <sstream> 

#include "../params.hh"
#include "../settings.hh"
#include "filereader_FASTEST.hh"


namespace CoupledField
{

    FileReader_FASTEST::FileReader_FASTEST(const std::string& name,
                           const UInt dim,
                           const UInt numFiles) :
        FileReader(name, dim, numFiles)
    {
    }

    FileReader_FASTEST::~FileReader_FASTEST()
    {
    }

    void FileReader_FASTEST::Init()
    {
      Settings& settings = Settings::Instance();
        std::string filename;
        char buf[128];
        UInt dummy;
        UInt i;
        
        
        if(settings.GetDouble("timeStep") < 0)
          EXCEPTION("No proper time step has been specified! Use --timestep X.");
        
        // Let's first determine the format if filenames        
        std::stringstream sstr;
        std::string regionFormatStr;
        std::string timeStepFormatStr;
        
        std::cout << "Trying to determine file name format for "
                  << ".coord/.node files..." << std::endl;
        for(i=1; i<10; i++)
        {
          filename = basename_;
          sstr.str("");
          sstr << "%0" << i << "i";
          regionFormatStr = sstr.str();
          sprintf(buf, regionFormatStr.c_str(), 1);
          filename+= buf;
          filename+= ".coord";

          infile.clear();
          infile.open(filename.c_str());  
          if (!infile) {
            std::cerr << "Can't open " << filename << std::endl;
          }
          else
            break;
        }
        
        if( i == 10 )
          EXCEPTION("Could not determine file name format for .coord/.node files");
        
        infile >> dummy;
        infile >> numPartitions_;
        infile >> dummy;
        numFiles_ = dummy < numFiles_ ? dummy : numFiles_;

        infile.close();
        infile.clear();

        std::cout << "Trying to determine file name format for .dat files..."
                  << std::endl;
        for(i=1; i<10; i++)
        {
          filename = basename_;
          sprintf(buf, regionFormatStr.c_str(), 1);
          filename+= buf;
          filename+= "_";
          sstr.str("");
          sstr << "%0" << i << "i";
          timeStepFormatStr = sstr.str();
          sprintf(buf, timeStepFormatStr.c_str(), 1);
          filename+= buf;
          filename+= ".dat";

          infile.clear();
          infile.open(filename.c_str());  
          if (!infile) {
            std::cerr << "Can't open " << filename << std::endl;
          }
          else
            break;
        }
        
        if( i == 10 )
          EXCEPTION("Could not determine file name format for .dat files");

        infile >> numResults_;

        infile.close();
        infile.clear();

        // Initialize mapping between columns in FASTEST result files and
        // internal variables.
        Integer dataCol;
        std::vector<std::string> fastestDOFs;
        fastestDOFs.push_back("lhsrc");
        fastestDOFs.push_back("vx");
        fastestDOFs.push_back("vy");
        fastestDOFs.push_back("vz");
        fastestDOFs.push_back("pres");
        
        for(UInt i=0; i<fastestDOFs.size(); i++)
        {
          std::string colStr = settings.GetString(fastestDOFs[i]);
          
          if(colStr != "")
          {
            sstr.str("");
            sstr.clear();
            colStr = colStr.substr(3);
            sstr << colStr;
            sstr >> dataCol;
            
            if(sstr.fail())
            {
              std::cerr << "Error while trying to init FASTEST column mapping"
              << std::endl;
              exit(1);
            }
            
            dataColumns_.push_back(dataCol-1);
          }
          else
          {
            dataColumns_.push_back(-1);            
          }
        }
        
        std::cout << " Name: " << name_ << std::endl
                  << " Dim: " << dim_ << std::endl
                  << " numfiles: " << numFiles_ << std::endl
                  << " numPartitions: " << numPartitions_ << std::endl
                  << " numResults: " << numResults_ << std::endl
                  << " timeStep: " << settings.GetDouble("timeStep") << std::endl
                  << " Lighthill source term column: " << dataColumns_[0]+1 << std::endl
                  << " vx column: " << dataColumns_[1]+1 << std::endl
                  << " vy column: " << dataColumns_[2]+1 << std::endl
                  << " vz column: " << dataColumns_[3]+1 << std::endl
                  << " pressure column: " << dataColumns_[4]+1 << std::endl;

        elsize_.resize(numPartitions_);
        MpCCInodes_.resize(numPartitions_);
        MpCCIelems_.resize(numPartitions_);

        for(UInt i=0; i<numPartitions_; i++)
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

            infile >> dummy;
            infile >> dummy;
            infile >> dummy;
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


    void FileReader_FASTEST::ReadNodalCoords(std::vector<Double> & NODECOORD,
                                             const UInt partitionIdx)
    {
     #ifdef TRACE
        (*trace) << "entering FileReader_FASTEST::ReadNodalCoords" << std::endl;
     #endif

        std::string filename;
        char buf[128];
        UInt dummy;
        
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
        max.resize(3);
        min.resize(3);

        infile >> NODECOORD[0] >> NODECOORD[1] >> NODECOORD[2];
        min[0] = max[0] = NODECOORD[0];
        min[1] = max[1] = NODECOORD[1];  
        min[2] = max[2] = NODECOORD[2];

        for (UInt i=1; i < MpCCInodes_[partitionIdx]; i++)
        {
            for (UInt j=0; j < 3; j++)
            {
                infile >> NODECOORD[3*i+j];
                if(NODECOORD[3*i+j] > max[j])
                    max[j] = NODECOORD[3*i+j];
                if(NODECOORD[3*i+j] < min[j])
                    min[j] = NODECOORD[3*i+j];
                //                if (i<=10) std::cout<<"NODECOORD["<<3*i+j<<"]: "<<NODECOORD[3*i+j]<<std::endl;
            }
        }

        std::cout << "MAX COORD: ( "<< max[0] << ", " << max[1] << ", "<< max[2] << ")"<< std::endl;
        std::cout << "MIN COORD: ( "<< min[0] << ", " << min[1] << ", "<< min[2] << ")"<< std::endl;
  

        infile.close();
    }

    void FileReader_FASTEST::ReadTopology(std::vector<UInt> & TOPOLOGYDATA,
                                          std::vector<UInt> & numNodesPerElem,
                                          std::vector<UInt> & elemTypes,
                                          const UInt partitionIdx)
    {
     #ifdef TRACE
        (*trace) << "entering FileReader_FASTEST::ReadTopology" << std::endl;
     #endif
        //        std::cout<<"in FileReader_FASTEST::ReadTopology"<<std::endl;  

        std::string filename;
        char buf[128];
        UInt dummy, numNodes, elemType;
        
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

        infile >> dummy >> dummy >> dummy >> dummy >> dummy;

        numNodes = elsize_[partitionIdx];
        switch(numNodes)
        {
	case 2:
	  {
	    elemType = ET_LINE2;
	    break;
	  }
        case 3:
          {
            elemType = ET_TRIA3;
            break;
          }
        case 4:
          {
            if(dim_ == 3)
              elemType = ET_TET4;
            else
              elemType = ET_QUAD4;
            break;
          }
	case 8:
	  {
	    if (dim_ == 3)
	      elemType = ET_HEXA8;
	    else
	      elemType = ET_QUAD8;
	    break;
          }
        case 5:
	  {
	    elemType = ET_PYRA5;
	    break;
	  }
	case 6:
	  {
	    if (dim_ == 3)
	      elemType = ET_WEDGE6;
	    else
	      elemType = ET_TRIA6;
	    break;
	  }
        }

        TOPOLOGYDATA.resize(elsize_[partitionIdx]*MpCCIelems_[partitionIdx]);
        numNodesPerElem.resize(MpCCIelems_[partitionIdx]);
        elemTypes.resize(MpCCIelems_[partitionIdx]);

        for (UInt i=0; i < MpCCIelems_[partitionIdx]; i++)
        {
            for (UInt j=0; j < elsize_[partitionIdx]; j++)
            {
                infile >> TOPOLOGYDATA[elsize_[partitionIdx]*i+j];
            }

            numNodesPerElem[i] = numNodes;
            elemTypes[i] = elemType;
            
            //            if (i<=10) std::cout<<"topology["<<i<<"] :"<<TOPOLOGYDATA[elsize_[partitionIdx]*i+0]
            //                                <<" "<<TOPOLOGYDATA[elsize_[partitionIdx]*i+1]<<" "
            //                                <<TOPOLOGYDATA[elsize_[partitionIdx]*i+2]<<std::endl;

        }

        infile.close();
    }

    void FileReader_FASTEST::ReadNodalValues(std::vector<double> & flowdata,
                                             const UInt partitionIdx,
                                             const UInt timeStepIdx)
    {
     #ifdef TRACE
        (*trace) << "entering FileReader_FASTEST::ReadNodalValues" << std::endl;
     #endif


        std::string filename;
        char buf[128];
        UInt dummy;

        
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

        std::vector<Double> tempVec;
        tempVec.resize(numResults_);
        
        for (UInt i=0; i < MpCCInodes_[partitionIdx]; i++)
        {
          std::fill(&flowdata[i*7+0], &flowdata[i*7+7], 0.0);
          
          for(UInt j=0; j<numResults_; j++)
            infile >> tempVec[j];

          for(UInt j=0; j<dataColumns_.size(); j++)
          {
            if(dataColumns_[j] > -1)
              flowdata[i*7+j] = tempVec[dataColumns_[j]];
          }
        }

        infile.close();

        
    }

} // end of namespace
