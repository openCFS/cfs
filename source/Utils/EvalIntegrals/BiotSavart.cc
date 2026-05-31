#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>

#include "BiotSavart.hh"
#include "Utils/mathParser/mathParser.hh"
#include "Domain/Domain.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

namespace CoupledField {

// Define / declare logging stream
DEFINE_LOG(bisa, "biotSavart")

  BiotSavart::BiotSavart( ) {
    formulation_ = BiotSavart::UNDEF;
    fieldIsComputed_ = false;
    isAxi_ = false;
    dim_ = 0;
    numVecComponents_ =  0;
  }
  
  BiotSavart::~BiotSavart() {
  }

  void BiotSavart::Init( PtrParamNode currentNode, 
                         Grid *ptGrid ,
                         shared_ptr<EqnMap> map ) {
    
    LOG_TRACE(bisa) << "BiotSavart: Started initialization" ;
    myParam_ = currentNode;
    ptGrid_ = ptGrid;
    eqnMap_ = map;
    
    dim_ = ptGrid_->GetDim();
    isAxi_ = ptGrid_->IsAxi();
    
    // Read number of coils
    ParamNodeList coils = myParam_->GetList("coil");
    coils_.Reserve(coils.GetSize());
    
    ParamNodeList::iterator it = coils.Begin();
    std::string name, file, value;
    
    for( ; it != coils.End(); it++ ) {

      (*it)->GetValue( "name", name );
      (*it)->GetValue( "file", file );
      (*it)->GetValue( "value", value );

      LOG_DBG(bisa) << "Defining coil '" << name << "'";

      // initialize coil struct
      BsCoil actCoil;
      actCoil.name = name;
      actCoil.mHandle_ = domain->GetMathParser()->GetNewHandle();
      domain->GetMathParser()->SetExpr( actCoil.mHandle_, value );

      // read in real coil coordinates
      ReadFile( file, actCoil );

      // store coil struct
      coils_.Push_back(actCoil);
    }
    
    // Print definition of coils to grid, if activated by the 
    // user
    if(myParam_->Get("showInGrid")->As<bool>()) {
      GenGridRepresentation();
    }
    
    
    // If formulation was already set, re-set it
    if( formulation_ != UNDEF ) {
      SetFormulation(formulation_);
    }
    LOG_TRACE(bisa) << "BiotSavart: Finished initialization";
  }

  void BiotSavart::GenGridRepresentation() {

    LOG_TRACE(bisa) << "Starting grid representation of coils";
    
    // loop over all coils
    StdVector<BsCoil>::iterator coilIt = coils_.Begin();
    for( ; coilIt != coils_.End(); coilIt++ ) {
      const BsCoil & actCoil = *coilIt;

      LOG_DBG(bisa) << "adding coil '" << actCoil.name << "'";
      UInt actElem = ptGrid_->GetNumElems()+1;
      UInt actNodeNum = ptGrid_->GetNumNodes()+1;
      StdVector<UInt> elemNums;
      
      // loop over all points
      UInt numPts = actCoil.points_.GetSize();
      for( UInt i = 1; i < numPts; i++ ) {
        
        // in 2D grids, we can only have x and y component
        Vector<Double> actCoord  = actCoil.points_[i]; 
        if( dim_ == 2) 
          actCoord[2] = 0.0;

        // Add nodes to grid
        if( i == 1 ) {
          ptGrid_->AddNodes(2);
          Vector<Double> prevCoord = actCoil.points_[i];
          // in 2D grids, we can only have x and y component
          if( dim_ == 2) 
            prevCoord[2] = 0.0;
          
          ptGrid_->SetNodeCoordinate( actNodeNum++, prevCoord );
        } else {
          ptGrid_->AddNodes(1);
        }
        ptGrid_->SetNodeCoordinate( actNodeNum, actCoord );
        
        // Add line element
        StdVector<UInt> connect(2);
        connect[0] = actNodeNum-1;
        connect[1] = actNodeNum;
        ptGrid_->AddElems(1);
        ptGrid_->SetElemData( actElem, Elem::ET_LINE2, NO_REGION_ID, &connect[0] );
        LOG_DBG3(bisa) << "added line element " << actElem 
            << " with connectivity " << connect.ToString();
        elemNums.Push_back(actElem++);
        actNodeNum++;
        
      } // loop over points
      
      // Add named elements
      ptGrid_->AddNamedElems( actCoil.name, elemNums );
      
    } // loop over coils
    
    LOG_TRACE(bisa) << "Finished grid representation of coils";
  }
  
  void BiotSavart::SetFormulation( FormulationType formulation ) {

    // re-set flag for field computation, if formulation has changed
    if( formulation != formulation_) {
      formulation_ = formulation;
      fieldIsComputed_ = false;
    }

    // Set function pointer for evaluating a single segment

    switch( formulation_ ) {
      case VEC_POT:
        if( isAxi_ ) {
          ptSegmentFunc_ = &BiotSavart::ArcIntegralVecPot;
        } else {
          ptSegmentFunc_ = &BiotSavart::LineIntegralVecPot;
        }
        // determine number of vector components for A
        if (dim_ == 0) // class is not yet initialized
          break;
        if( dim_ == 2 ) {
          numVecComponents_ = 1;
        } else if( dim_ == 3 ) {
          numVecComponents_ = 2;
        } else {
          EXCEPTION("Biot-Savart formulation only possible for 2D and 3D");
        }
        break;
      case MAG_FIELD_STRENGTH:
        if( isAxi_ ) {
          WARN("Arc representation of Biot-Savart coil for H-field"
              "formulation is not implemented yet.")
        } else {
          ptSegmentFunc_ = &BiotSavart::LineIntegralHField;
        }
        // determine number of vector components for H
        if (dim_ == 0) // class is not yet initialized
          return;
        if( dim_ == 2 ) {
          numVecComponents_ = 2;
        } else if( dim_ == 3 ) {
          numVecComponents_ = 3;
        } else {
          EXCEPTION("Biot-Savart formulation only possible for 2D and 3D");
        }
        break;
      default:
        EXCEPTION("No valid formulation set for Biot-Savart coil class");
    }
  }

  Vector<Double>& BiotSavart::CalcFieldAllEqns( bool normalized ) {

    LOG_TRACE(bisa) << "Starting computation of initial field";
    REFACTOR;

    return field_;
  }

  Vector<Double>& BiotSavart::CalcFieldDeriv1AllEqns(FormulationType type) {
   EXCEPTION("Not implemented yet"); 
  }

  Double BiotSavart::CalcFieldSingleEqn( UInt eqn ) {
    return CalcFieldAllEqns(false)[eqn];
  }

  void BiotSavart::CalcFieldAtPoint( Vector<Double>& field, 
                                     const Vector<Double>& observer ) {
    
    field.Resize( numVecComponents_ );
    field.Init();
    Vector<Double> tempField(numVecComponents_), tempField2;
    
    // obtain math parser
    MathParser * mp = domain->GetMathParser();
    Double factor = 0.0;
    
    // loop over all coils
    StdVector<BsCoil>::iterator coilIt = coils_.Begin();
    for( ; coilIt != coils_.End(); coilIt++ ) {
      BsCoil & actCoil = *coilIt;

      // Trigger for every coil the calculation of the contribution
      // to the observer point
      this->ComputeFieldNormalized( tempField, observer, actCoil );
      
      // scale it by expression and add it to own solution
      factor = mp->Eval(actCoil.mHandle_);
      
      

      // sum up contribution
      for(UInt i = 0; i < numVecComponents_; i++ ) {
        field[i] += tempField[i] * factor;
      }
    }
    // Loop over all coils
  }

  void BiotSavart::ReadFile(std::string fileName, BsCoil& coil ) {

    LOG_DBG(bisa) << "Reading coil definition from file '" << fileName << "'";
    // open file
    std::ifstream file;
    file.open(fileName.c_str());
    if ( !file ) {
      EXCEPTION( "Failed to open file '" << fileName
                 << "' containing Biot-Savart coil definition!" );
    }
    file.clear(); // clear flags

    // we don't trust .eof()
    file.seekg(0,std::ios::end);
    std::string::size_type pos_end = file.tellg();

    // start from the beginning
    file.seekg(0,std::ios::beg);
    char* buf;

    // read whole file content in one string
    buf = new char[pos_end+1];
    file.read(buf, pos_end);
    file.close();
    buf[pos_end]=0;
    std::string data(buf);
    delete[] buf;

    // split lines and omit comments
    typedef boost::tokenizer<boost::char_separator<char> > Tok;
    boost::char_separator<char> sep("\n\r");
    Tok t(data, sep);
    Tok::iterator it, end;
    it = t.begin();
    end = t.end();

    std::string line;
    std::stringstream sstr;
    for(; it != end; it++) {
      line = *it;

      // strip whitespaces
      boost::trim(line);

      // big choice of signs for comment's
      if (line.length() == 0 || line[0] == '#' ||
          line[0] == '%' || line[0] == '!')
        continue;

      // put line into a string stream
      sstr.clear();
      sstr.str(line);     

      Vector<Double> p(3);
      // read x value from string stream
      sstr >> p[0];
      if(!sstr)
        EXCEPTION("A problem occurred while reading from '" << fileName
                  << "'.\nInvalid entry: '" << line <<"'");

      // read y value from string stream
      sstr >> p[1];
      if(!sstr)
        EXCEPTION("A problem occurred while reading from '" << fileName
                  << "'.\nInvalid entry: '" << line <<"'");

      // read z value from string stream (only in 3D)
      // ToDo: This is not true in general, as the coil is always
      // defined in 3D. Only the resulting field could be in a 2D domain.
      //if( dim_ == 3) {
        sstr >> p[2];
        if(!sstr)
          EXCEPTION("A problem occurred while reading from '" << fileName
                    << "'.\nInvalid entry: '" << line <<"'");
      //}

      // store vector in current coil struct
      LOG_DBG2(bisa) << "Adding point " << p.ToString();
      coil.points_.Push_back(p);
    }
  }

  
  
  void BiotSavart::ComputeFieldNormalized( Vector<Double>& field,
                                           const Vector<Double>& observer,
                                           const BsCoil& coil) {

    // calculate field (A / H)  due to Biot-Savart
    field.Resize(3);
    field.Init();

    Vector<Double> start(dim_), end(dim_);
    Vector<Double> partField(3);
    UInt numPoints = coil.points_.GetSize();
    for( UInt iSeg = 1; iSeg < numPoints; iSeg++ ) {
      start = coil.points_[iSeg-1];
      end   = coil.points_[iSeg];
      
      // Call function pointer for evaluating the field along a
      // a given segment.
      LOG_DBG3(bisa) << "Calculating field for stick "
                            << start.ToString() << " -> " << end.ToString();
      (this->*ptSegmentFunc_)( partField, observer, start, end );
      field += partField;
    }
  }

  void BiotSavart::ArcIntegralVecPot( Vector<Double>& partMagVec, 
                                      const Vector<Double>& observer,
                                      const Vector<Double>& start, 
                                      const Vector<Double>& end ) {
    
    EXCEPTION( "Not yet implemented" );
  }

  

  void BiotSavart::LineIntegralVecPot( Vector<Double>& partMagVec, 
                                       const Vector<Double>& observer,
                                       const Vector<Double>& start, 
                                       const Vector<Double>& end ) {
    
    // calculate vector potential due to Biot-Savart for
    // a straight line element ("current stick")
    Double const piTimes4=16.0*atan(1.0), mue0=piTimes4*1e-7;

    partMagVec.Resize(3);
    partMagVec.Init(3);
    Vector<Double> vec1(3), vec2(3), vec3(3), ob;

    vec1 = end - start;
    vec2 = start - observer;
    vec3 = end - observer;

    Double norm1,norm2,norm3,dotProd13,dotProd12,fac;

    norm1 = vec1.NormL2();
    norm2 = vec2.NormL2();
    norm3 = vec3.NormL2();
    dotProd13 = vec1 * vec3;
    dotProd12 = vec1 * vec2;

    Double denom = (dotProd12 + norm1*norm2);
    
    // check, if observer point is located directly on coil 
    if( denom < EPS) {
      EXCEPTION("Distance from coil segment ("
          << start.ToString() << ") <-> (" << end.ToString() 
          << ") to observer (" << observer.ToString() << ") is 0");
    }
    fac = log((dotProd13 + norm1*norm3)/ denom);
    fac *= mue0/piTimes4;
    partMagVec = vec1 / norm1 * fac;
  }

  void BiotSavart::KernelVecPot( const Vector<Double>& observer, 
                                   Vector<Double>& p,
                                   Vector<Double>& dir, 
                                   Vector<Double>& kernelVP) {
      EXCEPTION( "Not yet implemented" );
  }
  
  void BiotSavart::LineIntegralHField( Vector<Double>& partMagVec, 
                                       const Vector<Double>& observer,
                                       const Vector<Double>& start, 
                                       const Vector<Double>& end ) {
    
    // calculate vector potential due to Biot-Savart for
    // a straight line element ("current stick")
    Double const piTimes4=16.0*atan(1.0);
    partMagVec.Resize(3);
    partMagVec.Init(3);
    Vector<Double> vec1(3), vec2(3), vec3(3), cross, ob;
    
    // ensure that observer has 3 entries
    ob = observer;
    if( ob.GetSize() == 2 ) {
      ob.Push_back(0.0);
    }
    
    vec1 = end - start;
    vec2 = start - ob;
    vec3 = end -ob;
    
    // Double norm1;
    Double norm2,norm3,dotProd13,dotProd12,fac;
    // norm1 = vec1.NormL2();
    norm2 = vec2.NormL2();
    norm3 = vec3.NormL2();
    dotProd13 = vec1 * vec3;
    dotProd12 = vec1 * vec2;
    
    // calculate cross product
    vec3.CrossProduct(vec1,cross);
    fac = dotProd13 / norm3 - dotProd12/norm2;
    fac *= 1.0/piTimes4;
    partMagVec = cross  * fac / (cross*cross);
    
  }

  void BiotSavart::SetTimeFncValue() {
    EXCEPTION( "Not yet implemented" );
  }

} // end of namespace
