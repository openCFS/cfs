#include <boost/filesystem/fstream.hpp>
#include "boost/filesystem/operations.hpp"
#include <boost/tokenizer.hpp>  //TODO

#include "CoefFunctionScatteredData.hh"
#include "FeBasis/FeSpace.hh"

#include "DataInOut/ScatteredDataInOut/ScatteredDataReader.hh"

namespace CoupledField{

  template<typename T, UInt DOFS>
  CoefFunctionScatteredData<T,DOFS>::CoefFunctionScatteredData(PtrParamNode& scatteredDataNode)
    : CoefFunction(),
      qid_("default"),
      factor_(1.0),
      interpolAlgo_(SHEPARD),
      numNeighbors_(20),
      p_(2),
      knnLib_(CGAL)
  {
    dimType_ = VECTOR;
    dependType_ = CoefFunction::GENERAL;

    isComplex_ =  std::is_same<T,Complex>::value;

    // Obtain id of quantity this CoefFunctionScatteredData should handle.
    qid_ = scatteredDataNode->Get("quantityId")->As<std::string>();

    // Obtain base node for all scattered data inputs.
    PtrParamNode scatteredInputNode = scatteredDataNode->GetRoot()
      ->Get("fileFormats")->Get("scatteredData");

    // Create all listed scattered data input readers.
    ScatteredDataReader::CreateReaders(scatteredInputNode);

    // Now let's register our desired quantity with the corresponding reader.
    ScatteredDataReader::RegisterQuantity(qid_);

    // Get hold of quantity node for obtaining our parameters.
    ParamNodeList scatteredNodes;
    scatteredNodes = scatteredInputNode->GetChildren();
    for(UInt i=0, n=scatteredNodes.GetSize(); i<n; i++) {
      std::string id = scatteredNodes[i]->Get("id")->As<std::string>();

      ParamNodeList quantityNodes;
      quantityNodes = scatteredNodes[i]->GetList("quantity");
      // Let's first check if we have unique ids.
      for(UInt j=0, m=quantityNodes.GetSize(); j<m; j++) {
        std::string qid = quantityNodes[j]->Get("id")->As<std::string>();

        if(qid == qid_)
        {
          quantityNode_ = quantityNodes[j];
        }
      }
    }

    if(quantityNode_->Has("interpolAlgo")) 
    {
      std::string interpolAlgo;
      interpolAlgo = quantityNode_->Get("interpolAlgo")->As<std::string>();
      if(interpolAlgo == "nearest-neighbor") 
      {
        interpolAlgo_ = NEAREST_NEIGHBOR;
      }
    }    

    if(quantityNode_->Has("knnLib")) 
    {
      std::string knnLib;
      knnLib = quantityNode_->Get("knnLib")->As<std::string>();
      if(knnLib == "flann") 
      {
        knnLib_ = FLANN;
      } 
    }    

    if(quantityNode_->Has("numNeighbors")) 
    {
      numNeighbors_ = quantityNode_->Get("numNeighbors")->As<UInt>();
      
      if(numNeighbors_ < 2) 
      {
        interpolAlgo_ = NEAREST_NEIGHBOR;
        numNeighbors_ = 1;
      }
    }

    if(quantityNode_->Has("p")) 
    {
      p_ = quantityNode_->Get("p")->As<Double>();
    }
    
    if(quantityNode_->Has("factor")) 
    {
      factor_ = quantityNode_->Get("factor")->As<Double>();
    }

    if(quantityNode_->Has("searchRadius"))
    {
      searchRadius_ = quantityNode_->Get("searchRadius")->As<Double>();
    }
  }
  
  template<typename T, UInt DOFS>
  void CoefFunctionScatteredData<T,DOFS>::GetQuantityData(bool updateMode)
  {
    ScatteredDataReader::Read(updateMode);

    if(quantityNode_->Has("bbox")) 
    {
      PtrParamNode bboxNode = quantityNode_->Get("bbox");

      boost::array<Double,6> bbox;
      bbox[0] = bboxNode->Get("xmin")->As<Double>();
      bbox[1] = bboxNode->Get("ymin")->As<Double>();
      bbox[2] = bboxNode->Get("zmin")->As<Double>();
      bbox[3] = bboxNode->Get("xmax")->As<Double>();
      bbox[4] = bboxNode->Get("ymax")->As<Double>();
      bbox[5] = bboxNode->Get("zmax")->As<Double>();

      std::vector< std::vector<double> > coords;
      std::vector< std::vector<T> > data;  // CHANGED
      ScatteredDataReader::GetQuantity(qid_, coords, data);

      UInt n = data.size();

      for(UInt i=0; i<n; i++)
      {
        Double& x = coords[i][0];
        Double& y = coords[i][1];
        Double& z = coords[i][2];
        
        if( x >= bbox[0] && x <= bbox[3] &&
            y >= bbox[1] && x <= bbox[4] &&
            z >= bbox[2] && z <= bbox[5]) 
        {
          coordinates_.push_back(coords[i]);
          scatteredData_.push_back(data[i]);
        }      
      }
    }
    else
    {
      ScatteredDataReader::GetQuantity(qid_, coordinates_, scatteredData_);
    }
  }


  template<typename T, UInt DOFS>
  void CoefFunctionScatteredData<T,DOFS>::DumpData() 
  {
    if(!quantityNode_->Has("dump")) 
    {
      return;
    }
    
    PtrParamNode dumpNode = quantityNode_->Get("dump");
    std::string fileName = dumpNode->Get("fileName")->As<std::string>();
    std::string format = "csv";
    if(dumpNode->Has("format")) 
    {
      format = dumpNode->Get("format")->As<std::string>();
    }
    
    std::ofstream csv(fileName.c_str(), std::ios_base::binary);
    if(!csv) 
    {
      EXCEPTION("Could not open scattered data dump file '" << fileName << "'.");
    }
    
    // Write title line to CSV
    UInt m=scatteredData_[0].size();
    csv << "x,y,z";
    for(UInt i=0; i<m; i++)
    {
      csv << "," << qid_ << ":" << i;
    }    
    csv << std::endl;

    // Write actual data to CSV
    UInt n = scatteredData_.size();
    for(UInt i=0; i<n; i++)
    {
      Double& x = coordinates_[i][0];
      Double& y = coordinates_[i][1];
      Double& z = coordinates_[i][2];
      
      csv << x << ","
          << y << ","
          << z;
      for(UInt j=0; j<m; j++)
      {
        csv << "," << scatteredData_[i][j] * factor_;
      }    
      csv << std::endl;
    }
    csv.close();
  }

  template<typename T, UInt DOFS>
  void CoefFunctionScatteredData<T,DOFS>::Read(bool updateMode)
  {
    GetQuantityData(updateMode);
    DumpData();



    UInt n = scatteredData_.size();

    switch(knnLib_)
    {
    case CGAL:
#ifdef USE_CGAL
      {
      std::vector<CGAL::Point> points;
      points.resize(n);
      for(UInt i=0; i<n; i++)
      {
        points[i] = (CGAL::Point(coordinates_[i][0],
                               coordinates_[i][1],
                               coordinates_[i][2],
                               scatteredData_[i][0] * factor_,
                               scatteredData_[i][1] * factor_,
                               scatteredData_[i][2] * factor_));
      }

      searchTree_.reset(new Tree(points.begin(), points.end()));
      if(updateMode)
        return;
    }
#else
      EXCEPTION("CGAL not supported! Compile with USE_CGAL=ON.");
#endif
      break;

    case FLANN:
#ifdef USE_FLANN
      {
        if(updateMode)
          return;
      dataset_.reset(new flann::Matrix<Double>(new Double[n*3], n, 3));
      Double *dPtr = dataset_->ptr();
      for(UInt i=0; i<n; i++)
      {
        UInt idx = i*3;
        
        dPtr[idx+0] = coordinates_[i][0];
        dPtr[idx+1] = coordinates_[i][1];
        dPtr[idx+2] = coordinates_[i][2];
      }
      
      // construct an randomized kd-tree index using a single kd-tree
      index_.reset(new flann::Index<flann::L2<Double> >(*dataset_.get(),
                                                        flann::KDTreeSingleIndexParams(12)));
      index_->buildIndex();
      }
#else
      EXCEPTION("FLANN not supported! Compile with USE_FLANN=ON.")
#endif
      break;
    default:
      EXCEPTION("Unknown k-nearest neighbors library '" << knnLib_ 
                << "' specified for quantity id '" << qid_ << "'.")
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
    lpm.GetGlobal(globPoint);
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
      lpm.GetGlobal(globPoint);
      T divergence = 0.0;
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

    // Thread-safe Read: protect concurrent access during parallel assembly
    // - Initial read (updateMode=false) is done exactly once via call_once
    // - Update reads (updateMode=true) for time-varying data are serialized via mutex
    //   The CSVT reader checks if time step changed internally, so most update calls are no-ops
    {
      std::lock_guard<std::mutex> lock(readMutex_);
      bool updateMode = !scatteredData_.empty();
      if (!updateMode) {
        // First read - use call_once to ensure exactly one thread does it
        std::call_once(readOnce_, [this]() {
          Read(false);
        });
      } else {
        // Update mode - needed for time-varying data (CSVT)
        // The reader internally checks if data needs updating for current time step
        Read(true);
      }
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
    CGAL::Point query(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    if(globPoint.GetSize() == 2)
    {
      CGAL::Point query2(globPoint[0], globPoint[1], 0.0, 0.0, 0.0, 0.0);
      query = query2;
    }
    else
    {
      CGAL::Point query3(globPoint[0], globPoint[1], globPoint[2], 0.0, 0.0, 0.0);
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
      neighbors[i].Resize(globPoint.GetSize());
      vectors[i].Resize(globPoint.GetSize());

      if(quantityNode_->Has("searchRadius"))
      {
        if(l2Distances[i]<searchRadius_){
          if(globPoint.GetSize() == 2)
          {
            it->first.vx(vectors[i][0]);
            it->first.vy(vectors[i][1]);
            neighbors[i][0] = it->first.x();
            neighbors[i][1] = it->first.y();
          }
          else
          {
            it->first.vx(vectors[i][0]);
            it->first.vy(vectors[i][1]);
            it->first.vz(vectors[i][2]);

            neighbors[i][0] = it->first.x();
            neighbors[i][1] = it->first.y();
            neighbors[i][2] = it->first.z();
          }
        }else{
          if(globPoint.GetSize() == 2)
          {
            vectors[i][0] = 0.0;
            vectors[i][1] = 0.0;
            neighbors[i][0] = it->first.x();
            neighbors[i][1] = it->first.y();
          }
          else
          {
            vectors[i][0] = 0.0;
            vectors[i][1] = 0.0;
            vectors[i][2] = 0.0;

            neighbors[i][0] = it->first.x();
            neighbors[i][1] = it->first.y();
            neighbors[i][2] = it->first.z();
          }
        }
      }
      else{
        if(globPoint.GetSize() == 2)
        {
          it->first.vx(vectors[i][0]);
          it->first.vy(vectors[i][1]);
          neighbors[i][0] = it->first.x();
          neighbors[i][1] = it->first.y();
        }
        else
        {
          it->first.vx(vectors[i][0]);
          it->first.vy(vectors[i][1]);
          it->first.vz(vectors[i][2]);

          neighbors[i][0] = it->first.x();
          neighbors[i][1] = it->first.y();
          neighbors[i][2] = it->first.z();
        }

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

    if(globPoint.GetSize()==2)
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

        if(quantityNode_->Has("searchRadius"))
        {
          if(l2Distances[j]<searchRadius_){
            UInt idx = indices[i][j];

            neighbors[j].Resize(DOFS);
            vectors[j].Resize(DOFS);

            for(UInt d=0; d<DOFS; d++)
            {
              vectors[j][d] = scatteredData_[idx][d] * factor_;
              neighbors[j][d] = coordinates_[idx][d];
            }
          }else{
            UInt idx = indices[i][j];

            neighbors[j].Resize(DOFS);
            vectors[j].Resize(DOFS);

            for(UInt d=0; d<DOFS; d++)
            {
              vectors[j][d] = scatteredData_[idx][d] * 0.0;
              neighbors[j][d] = coordinates_[idx][d];
            }
          }
        }else{
          UInt idx = indices[i][j];

          neighbors[j].Resize(DOFS);
          vectors[j].Resize(DOFS);

          for(UInt d=0; d<DOFS; d++)
          {
            vectors[j][d] = scatteredData_[idx][d] * factor_;
            neighbors[j][d] = coordinates_[idx][d];
          }
        }

      }
    }

    delete[] indices.ptr();
    delete[] dists.ptr();
  }
#endif  


  template class CoefFunctionScatteredData<Double,1>;
  template class CoefFunctionScatteredData<Complex,1>;

  template class CoefFunctionScatteredData<Double,2>;
  template class CoefFunctionScatteredData<Complex,2>;

  template class CoefFunctionScatteredData<Double,3>;
  template class CoefFunctionScatteredData<Complex,3>;
}// end of namespace

