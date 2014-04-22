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
      p_(2),
      knnLib_(CGAL)
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

    if(scatteredDataNode->Has("knnLib")) 
    {
      std::string knnLib;
      knnLib = scatteredDataNode->Get("knnLib")->As<std::string>();
      if(knnLib == "flann") 
      {
        knnLib_ = FLANN;
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
    

    if(scatteredDataNode->Has("factor")) 
    {
      factor_ = scatteredDataNode->Get("factor")->As<Double>();
    }
  }
  
  template<typename T, UInt DOFS>
  void CoefFunctionScatteredData<T,DOFS>::Read() 
  {
    if(!boost::filesystem::exists(fileName_))
    {
      EXCEPTION("Scattered data file '" << fileName_ << "' does not exist!")
    }

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
    UInt n = scatteredData_.size();

    switch(knnLib_)
    {
    case CGAL:
#ifdef USE_CGAL
      {
      std::list<Point> points;

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

      searchTree_.reset(new Tree(points.begin(), points.end()));
    }
#else
      EXCEPTION("CGAL not supported! Compile with USE_CGAL=ON.");
#endif
      break;

    case FLANN:
#ifdef USE_FLANN
      {
      dataset_.reset(new flann::Matrix<Double>(new Double[n*3], n, 3));
      Double *dPtr = dataset_->ptr();
      for(UInt i=0; i<n; i++)
      {
        UInt idx = i*3;
        
        switch(dof2CoordColumn_.size())
        {
        case 2:
          dPtr[idx+0] = scatteredData_[i][dof2CoordColumn_[0]];
          dPtr[idx+1] = scatteredData_[i][dof2CoordColumn_[1]];
          dPtr[idx+2] = 0.0;        
          break;        
        case 3:
          dPtr[idx+0] = scatteredData_[i][dof2CoordColumn_[0]];
          dPtr[idx+1] = scatteredData_[i][dof2CoordColumn_[1]];
          dPtr[idx+2] = scatteredData_[i][dof2CoordColumn_[2]];
          break;        
        }    
      }
      
      // construct an randomized kd-tree index using a single kd-tree
      index_.reset(new flann::Index<flann::L2<Double> >(*dataset_.get(), flann::KDTreeSingleIndexParams(12)));
      index_->buildIndex();
      }
#else
      EXCEPTION("FLANN not supported! Compile with USE_FLANN=ON.")
#endif
      break;
    default:
      break;
    }

    if(n < numNeighbors_) 
    {
      numNeighbors_ = n;
    }    
  }
  
  template<typename T, UInt DOFS>
  void CoefFunctionScatteredData<T,DOFS>::GetVector( Vector<T>& vec, 
                                                     const LocPointMapped& lpm ) {
    Vector<double> globPoint(DOFS);
    lpm.shapeMap->Local2Global(globPoint, lpm.lp);
    this->InterpolateVector(globPoint,vec);
  }

  template<typename T, UInt DOFS>
  void CoefFunctionScatteredData<T,DOFS>::GetScalar( T & value,
                                                       const LocPointMapped& lpm ){
    //if we have no interpolation, we just throw an exception as everything is
    //only implemented for vector valued data....
    if(this->derivType_ == NONE){
      EXCEPTION("CoefFunctionScatteredData supports only vector valued input data")
    }

    if(this->derivType_ == VECTOR_DIVERGENCE){
      //here we do a very dirty thing. we just compute derivatives by
      //application of a four point stencil...
      //this is very dirty hack implementation to ensure compatibility
      // we will need the derivative w.r.t. the interpolation procedure used here...


      //obtain the global coordinate
      Vector<double> globPoint(DOFS);

      T divergence = 0.0;
      lpm.shapeMap->Local2Global(globPoint, lpm.lp);
      LocPointMapped tmpLocPoint = lpm;

      //temorarily change the dimType to Vector again
      this->dimType_ = VECTOR;

      // determine eps according to element jacobian
      // perhaps we should take the volume of the element instead?
      Vector<Double> dia;
      lpm.shapeMap->CalcDiameter(dia);
      Double eps = 1e-4;

      //First, the x-value

      for(UInt d = 0;d<DOFS;d++){
        T val1,val2,val3,val4;
        Double buffer = globPoint[d];
        Vector<T> curValue;

        globPoint[d] = buffer + 2*eps * dia[d];
        this->InterpolateVector(globPoint,curValue);
        val1 = curValue[d];

        globPoint[d] = buffer + 1*eps * dia[d];
        this->InterpolateVector(globPoint,curValue);
        val2 = curValue[d];

        globPoint[d] = buffer - 1*eps * dia[d];
        this->InterpolateVector(globPoint,curValue);
        val3 = curValue[d];

        globPoint[d] = buffer - 2*eps * dia[d];
        this->InterpolateVector(globPoint,curValue);
        val4 = curValue[d];

        divergence += (-val1 + 8.0 * val2 - 8.0 * val3 + val4 ) / (12.0 * eps * dia[d]);

        globPoint[d] = buffer;
      }

      this->dimType_ = SCALAR;
      value = divergence;
    }
  }


  template<typename T, UInt DOFS>
  void CoefFunctionScatteredData<T,DOFS>::InterpolateVector(Vector<Double> globPoint, Vector<T> & vec){
    if(!scatteredData_.size())
    {
      Read();
    }

    StdVector< Vector<Double> > neighbors;
    StdVector< Double > l2dists;
    StdVector< Vector<T> > vectors;

    switch(knnLib_) 
    {
    case CGAL:
#ifdef USE_CGAL
      KNNSearch_CGAL(globPoint, neighbors, l2dists, vectors);
#else
      EXCEPTION("CoefFunctionScatteredData needs to be compiled with USE_CGAL=ON!");
#endif
      break;
      
    case FLANN:
#ifdef USE_FLANN
      KNNSearch_FLANN(globPoint, neighbors, l2dists, vectors);
#else
      EXCEPTION("CoefFunctionScatteredData needs to be compiled with USE_FLANN=ON!");
#endif
      break;
    default:
      break;
    }

    vec.Resize(DOFS);

    Double dmin = l2dists[0];
    Double dmax = dmin;
    StdVector< Double >::iterator it, end;
    it = l2dists.Begin();
    end = l2dists.End();    

    for( ; it != end; ++it) {
      Double dist = (*it);
      dmin = dmin < dist ? dmin : dist;
      dmax = dmax > dist ? dmax : dist;
    }

    if(numNeighbors_ == 1 ||
       interpolAlgo_ == NEAREST_NEIGHBOR ||
       l2dists.GetSize() == 1 ||
       dmin/dmax < 1e-6)
    {
      // Apply nearest neigbor interpolation.
      for(UInt dof=0; dof < DOFS; dof++) 
      {
        vec[dof] = vectors[0][dof];
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
      it = l2dists.Begin();
      for(UInt i=0; it != end; ++it, i++) {
        Double d = *it;
        Double w = std::pow((R-d)/(R*d), p_);
        weights += w;
        
        for(UInt dof=0; dof < DOFS; dof++) 
        {
          sum[dof] += vectors[i][dof] * w;
        }
      }

      for(UInt dof=0; dof < DOFS; dof++) 
      {
        vec[dof] = sum[dof] / weights;
      }
    }
  }


  template<typename T, UInt DOFS>
    std::string CoefFunctionScatteredData<T,DOFS>::ToString() const {
    std::stringstream out;
    out << "CoefFunctionScatteredData<" << typeid(T).name() << ", " << DOFS << ">";
    return out.str();
  }
  
#ifdef USE_CGAL
  template<typename T, UInt DOFS>
    void CoefFunctionScatteredData<T,DOFS>::KNNSearch_CGAL(const Vector<Double> globPoint,
      StdVector< Vector<Double> >& neighbors,
      StdVector< Double >& l2Distances,
      StdVector< Vector<T> >& vectors)
  {
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

    UInt nn = std::distance(search.begin(), search.end());    
    neighbors.Resize(nn);
    l2Distances.Resize(nn);
    vectors.Resize(nn);
    
    for(UInt i=0 ; it != search.end(); ++it, i++) {
      l2Distances[i] = std::sqrt(it->second);
      neighbors[i].Resize(DOFS);
      vectors[i].Resize(DOFS);

      if(DOFS == 2)
      {
        vectors[i][0] = it->first.vx();
        vectors[i][1] = it->first.vy();

        neighbors[i][0] = it->first.x();
        neighbors[i][1] = it->first.y();
      }
      else
      {
        vectors[i][0] = it->first.vx();
        vectors[i][1] = it->first.vy();
        vectors[i][2] = it->first.vz();

        neighbors[i][0] = it->first.x();
        neighbors[i][1] = it->first.y();
        neighbors[i][2] = it->first.z();
      }
    }
  }
#endif  

#ifdef USE_FLANN
  template<typename T, UInt DOFS>
    void CoefFunctionScatteredData<T,DOFS>::KNNSearch_FLANN(const Vector<Double> globPoint,
      StdVector< Vector<Double> >& neighbors,
      StdVector< Double >& l2Distances,
      StdVector< Vector<T> >& vectors) 
  {
    Double q[3];

    if(DOFS == 2)
    {
      q[0] = globPoint[0];
      q[1] = globPoint[1];
      q[2] = 0.0;
    }
    else
    {
      q[0] = globPoint[0];
      q[1] = globPoint[1];
      q[2] = globPoint[2];
    }

    flann::Matrix<Double> query(q, 1,3);

    flann::Matrix<int> indices(new int[query.rows*numNeighbors_], query.rows, numNeighbors_);
    flann::Matrix<Double> dists(new Double[query.rows*numNeighbors_], query.rows, numNeighbors_);
    // do a knn search, using 3 checks
    index_->knnSearch(query, indices, dists, numNeighbors_, flann::SearchParams(flann::FLANN_CHECKS_UNLIMITED));

    neighbors.Resize(numNeighbors_);
    l2Distances.Resize(numNeighbors_);    
    vectors.Resize(numNeighbors_);    
    for(UInt i=0; i<indices.rows; i++) 
    {
      for(UInt j=0; j<numNeighbors_; j++) 
      {
        l2Distances[j] = std::sqrt(dists[i][j]);

        UInt idx = indices[i][j];
        
        neighbors[j].Resize(DOFS);
        vectors[j].Resize(DOFS);

        for(UInt d=0; d<DOFS; d++) 
        {
          vectors[j][d] = scatteredData_[idx][dof2ValueColumn_[d]] * factor_;
          neighbors[j][d] = scatteredData_[idx][dof2CoordColumn_[d]];
        }
      }
    }

    delete[] indices.ptr();
    delete[] dists.ptr();
  }
#endif  


#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class CoefFunctionScatteredData<Double,2>;
  template class CoefFunctionScatteredData<Complex,2>;

  template class CoefFunctionScatteredData<Double,3>;
  template class CoefFunctionScatteredData<Complex,3>;
#endif
}// end of namespace

