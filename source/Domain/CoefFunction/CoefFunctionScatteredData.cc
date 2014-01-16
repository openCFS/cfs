#include <boost/filesystem/fstream.hpp>
#include "boost/filesystem/operations.hpp"
#include <boost/tokenizer.hpp>

#ifdef __MINGW64__
#include <intrin.h>
#endif

#include "CoefFunctionScatteredData.hh"
#include "FeBasis/FeSpace.hh"

#include "def_use_ccmio.hh"
#include "DataInOut/ScatteredDataInOut/ScatteredDataReaderCSV.hh"
#ifdef USE_CCMIO
#include "DataInOut/ScatteredDataInOut/ScatteredDataReaderCCM.hh"
#endif

namespace CoupledField{

  template<typename T, UInt DOFS>
  CoefFunctionScatteredData<T,DOFS>::CoefFunctionScatteredData(PtrParamNode& scatteredDataNode)
    : CoefFunction(),
      fileName_(""),
      factor_(1.0),
      interpolAlgo_(SHEPARD),
      numNeighbors_(20),
      p_(2)
  {
    dimType_ = VECTOR;
    dependType_ = CoefFunction::GENERAL;

    fileName_ = scatteredDataNode->Get("fileName")->As<std::string>();

    format_ = scatteredDataNode->Get("format")->As<std::string>();

    if(scatteredDataNode->Has("interpolAlgo")) 
    {
      std::string interpolAlgo;
      interpolAlgo = scatteredDataNode->Get("interpolAlgo")->As<std::string>();
      if(interpolAlgo == "nearest-neighbor") 
      {
        interpolAlgo_ = NEAREST_NEIGHBOR;
      }
    }    

    if(scatteredDataNode->Has("numNeighbors")) 
    {
      numNeighbors_ = scatteredDataNode->Get("numNeighbors")->As<UInt>();
      
      if(numNeighbors_ < 2) 
      {
        interpolAlgo_ = NEAREST_NEIGHBOR;
        numNeighbors_ = 1;
      }
    }

    if(scatteredDataNode->Has("p")) 
    {
      p_ = scatteredDataNode->Get("p")->As<Double>();
    }
    
    if(format_ == "csv") 
    {
      ParamNodeList coordList;
      coordList = scatteredDataNode->Get("coordinates")->GetList("comp");
      for(UInt i=0, n=coordList.GetSize(); i<n; i++) {
        std::string dof = coordList[i]->Get("dof")->As<std::string>();
        
        if( dof == "x" ) {
          dof2CoordColumn_[0] = coordList[i]->Get("col")->As<UInt>();
        }
        if( dof == "y" ) {
          dof2CoordColumn_[1] = coordList[i]->Get("col")->As<UInt>();
        }
        if( dof == "z" ) {
          dof2CoordColumn_[2] = coordList[i]->Get("col")->As<UInt>();
        }
      }
      
      ParamNodeList valueList;
      valueList = scatteredDataNode->Get("values")->GetList("comp");
      for(UInt i=0, n=valueList.GetSize(); i<n; i++) {
        std::string dof = valueList[i]->Get("dof")->As<std::string>();
        
        if( dof == "x" ) {
          dof2ValueColumn_[0] = valueList[i]->Get("col")->As<UInt>();
        }
        if( dof == "y" ) {
          dof2ValueColumn_[1] = valueList[i]->Get("col")->As<UInt>();
        }
        if( dof == "z" ) {
          dof2ValueColumn_[2] = valueList[i]->Get("col")->As<UInt>();
        }      
      }
    } else 
    {
      dof2CoordColumn_[0] = 0;
      dof2CoordColumn_[1] = 1;
      dof2CoordColumn_[2] = 2;
      
      dof2ValueColumn_[0] = 3;
      dof2ValueColumn_[1] = 4;
      dof2ValueColumn_[2] = 5;
    }    
    
    factor_ = scatteredDataNode->Get("factor")->As<Double>();
  }
  
  template<typename T, UInt DOFS>
  void CoefFunctionScatteredData<T,DOFS>::Read() 
  {
    if(!boost::filesystem::exists(fileName_))
    {
      EXCEPTION("Scattered data file '" << fileName_ << "' does not exist!")
    }

#ifdef USE_CGAL
    ScatteredDataReaderPtr sdataReader;

    if(format_ == "csv") 
    {
      ScatteredDataReaderCSV* SCRCSV = new ScatteredDataReaderCSV(fileName_);
      SCRCSV->SetNumSkipLines(1);
      sdataReader.reset(SCRCSV);
    } else if (format_ == "ccm") {
#ifdef USE_CCMIO
      ScatteredDataReaderCCM* SCRCCM = new ScatteredDataReaderCCM(fileName_);
      std::vector<std::string> componentShortNames(3);
      componentShortNames[0] = "SU";
      componentShortNames[1] = "SV";
      componentShortNames[2] = "SW";
      SCRCCM->SetComponentShortNames(componentShortNames);
      sdataReader.reset(SCRCCM);
#else
      EXCEPTION("STAR-CCM+ files not supported! Compile with USE_CCMIO=ON.");
#endif
    } else {
      EXCEPTION("No format for scattered data file specified!");
    }      
    sdataReader->Read(scatteredData_);

    std::list<Point> points;

    UInt n = scatteredData_.size();
    for(UInt i=0; i<n; i++)
    {
      switch(dof2CoordColumn_.size())
      {
      case 2:
        points.push_back(Point(scatteredData_[i][dof2CoordColumn_[0]],
        scatteredData_[i][dof2CoordColumn_[1]],
        0.0,
        scatteredData_[i][dof2ValueColumn_[0]] * factor_,
        scatteredData_[i][dof2ValueColumn_[1]] * factor_,
        0.0));
        break;        
      case 3:
        points.push_back(Point(scatteredData_[i][dof2CoordColumn_[0]],
        scatteredData_[i][dof2CoordColumn_[1]],
        scatteredData_[i][dof2CoordColumn_[2]],
        scatteredData_[i][dof2ValueColumn_[0]] * factor_,
        scatteredData_[i][dof2ValueColumn_[1]] * factor_,
        scatteredData_[i][dof2ValueColumn_[2]] * factor_));
        break;        
      }    
    }
    if(n < numNeighbors_) 
    {
      numNeighbors_ = n;
    }    

    searchTree_.reset(new Tree(points.begin(), points.end()));
#else
    EXCEPTION("Switch on USE_CGAL for nearest-neighbor mapping!")
#endif
  }
  
  template<typename T, UInt DOFS>
  void CoefFunctionScatteredData<T,DOFS>::GetVector( Vector<T>& vec, 
                                                     const LocPointMapped& lpm ) {
    if(!scatteredData_.size())
    {
      Read();
    }    

#ifdef USE_CGAL
    vec.Resize(DOFS);
    //    UInt size = scatteredData_.size();
    
    Vector<double> globPoint(DOFS);
    Vector<T> diff(DOFS);
    lpm.shapeMap->Local2Global(globPoint, lpm.lp);

    Point query(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    if(DOFS == 2) 
    {
      Point query2(globPoint[0], globPoint[1], 0.0, 0.0, 0.0, 0.0);
      query = query2;
    }
    else 
    {
      Point query3(globPoint[0], globPoint[1], globPoint[2], 0.0, 0.0, 0.0);
      query = query3;
    }

    K_neighbor_search search(*searchTree_.get(), query, numNeighbors_);

    K_neighbor_search::iterator it = search.begin();

    if(it == search.end())
    {
      EXCEPTION("Could not find a nearest neighbor for " << globPoint << "!");
    }

    Double dmin = it->second;
    Double dmax = dmin;
    for( ; it != search.end(); ++it) {
      dmin = dmin < it->second ? dmin : it->second;
      dmax = dmax > it->second ? dmax : it->second;
    }
    dmin = std::sqrt(dmin);
    dmax = std::sqrt(dmax);

    if(numNeighbors_ == 1 ||
        interpolAlgo_ == NEAREST_NEIGHBOR ||
        std::distance(search.begin(), search.end()) == 1 ||
        dmin/dmax < 1e-6)
    {
      it = search.begin();

      // Apply nearest neigbor interpolation.
      if(DOFS == 2) 
      {
        vec[0] = it->first.vx();
        vec[1] = it->first.vy();
      }
      else 
      {
        vec[0] = it->first.vx();
        vec[1] = it->first.vy();
        vec[2] = it->first.vz();
      }
    }
    else 
    {    
      // Apply Shepard interpolation cf. Numerical Recipes 3rd ed. p. 143ff.
      // or http://www.ems-i.com/gmshelp/Interpolation/Interpolation_Schemes \
      // /Inverse_Distance_Weighted/Shepards_Method.htm
      Vector<T> sum(3);
      Double weights = 0.0;
      // The point which is farthest away, should at least have a non-zero
      // weight of 0.01. If we would choose R = dmax, it would not contribute
      // at all.
      Double R = 1.01 * dmax;

      // report the N nearest neighbors and their distance
      // This should sort all N points by increasing distance from origin
      for(it = search.begin(); it != search.end(); ++it) {
        Double d = std::sqrt(it->second);
        Double w = std::pow((R-d)/(R*d), p_);
        weights += w;
        sum[0] += it->first.vx() * w;
        sum[1] += it->first.vy() * w;
        sum[2] += it->first.vz() * w;
      }

      if(DOFS == 2) 
      {
        vec[0] = sum[0] / weights;
        vec[1] = sum[1] / weights;
      }
      else 
      {
        vec[0] = sum[0] / weights;
        vec[1] = sum[1] / weights;
        vec[2] = sum[2] / weights;
      }
    }
    
#else
    EXCEPTION("CoefFunctionScatteredData needs to be compiled with USE_CGAL=ON!");
#endif
  }

  template<typename T, UInt DOFS>
  std::string CoefFunctionScatteredData<T,DOFS>::ToString() const {
  std::stringstream out;
  out << "CoefFunctionScatteredData<" << typeid(T).name() << ", " << DOFS << ">";
  return out.str();
}


#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class CoefFunctionScatteredData<Double,2>;
  template class CoefFunctionScatteredData<Complex,2>;

  template class CoefFunctionScatteredData<Double,3>;
  template class CoefFunctionScatteredData<Complex,3>;
#endif
}// end of namespace

