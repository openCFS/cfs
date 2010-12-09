// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>

#include "biotSavart.hh"
#include "Domain/domain.hh"
#include "DataInOut/Logging/cfslog.hh"

namespace CoupledField {

// Define / declare logging stream
DECLARE_LOG(bisa)
DEFINE_LOG(bisa, "biotSavart")

  BiotSavart::BiotSavart( FormulationType type ) {
    formulation_ = type;
    fieldIsComputed_ = false;
    isAxi_ = false;
    dim_ = 0;
  }
  
  BiotSavart::~BiotSavart() {
  }

  void BiotSavart::Init( PtrParamNode currentNode, 
                         Grid *ptGrid ,
                         shared_ptr<EqnMap> map ) {
    
    // Currently, we just have the vector potential formulation
    if( formulation_ == SCALAR_POT ) {
      EXCEPTION( "Scalar potential formulation for Biot-Savart formulation "
          << "is currently not implemented!");
    }
    
    LOG_TRACE(bisa) << "BiotSavart: Started initialization" ;
    myParam_ = currentNode;
    ptGrid_ = ptGrid;
    eqnMap_ = map;
    
    dim_ = ptGrid_->GetDim();
    isAxi_ = param->Get("domain")->
        Get("geometryType")->As<std::string>() == "axi";
    
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
    
    // Print definition of coils to grid
    GenGridRepresentation();
    
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
        const Vector<Double> actCoord = actCoil.points_[i];  
        
        
        // Add nodes to grid
        if( i == 1 ) {
          ptGrid_->AddNodes(2);
          const Vector<Double> prevCoord = actCoil.points_[0];
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
        ptGrid_->SetElemData( actElem, Elem::LINE2, NO_REGION_ID, &connect[0] );
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
 
  Double BiotSavart::GetMagVec( UInt eqn, BasePDE::AnalysisType analysis ) {
    return GetMagVec(false)[eqn];
  }

  Vector<Double>& BiotSavart::GetMagVec(bool normalized) {
    LOG_TRACE(bisa) << "Starting GetMagVec";
    
    // Perform calculation of normalized field
    // Note: This is just done once, but we have to do it as late as possible
    // as we require already a reordered equation map object.
    ComputeBiotSavartField();
    
    // obtain math parser
    MathParser * mp = domain->GetMathParser();
    Double factor = 0.0;
    
    // initialize solution
    magVec_.Init();
    
    // loop over all coils
    StdVector<BsCoil>::iterator coilIt = coils_.Begin();
    for( ; coilIt != coils_.End(); coilIt++ ) {
      const BsCoil & actCoil = *coilIt;

      if (normalized) {
        magVec_ += actCoil.magVecNormalized_;
      }else {
        // scale it by expression and add it to own solution
        factor = mp->Eval(actCoil.mHandle_);
        magVec_ += actCoil.magVecNormalized_ * factor;
      }
    }
    LOG_TRACE(bisa) << "Starting GetMagVec";
    
    return magVec_;
  }

  Vector<Double>& BiotSavart::GetMagVecDeriv1(){
    EXCEPTION( "Not yet implemented" );
  }

  void BiotSavart::CalcH( Vector<Double>& HField,Vector<Double>& intPoint) {
    EXCEPTION( "Not yet implemented" );
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
    typedef boost::tokenizer<char_separator<char> > Tok;
    boost::char_separator<char> sep("\n\r");
    Tok t(data, sep);
    Tok::iterator it, end;
    it = t.begin();
    end = t.end();

    std::string line;
    std::stringstream sstr;
    for(UInt lineNum=1; it != end; it++, lineNum++) {
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
        EXCEPTION("A problem occured while reading from '" << fileName
                  << "'.\nInvalid entry: '" << line <<"'");

      // read y value from string stream
      sstr >> p[1];
      if(!sstr)
        EXCEPTION("A problem occured while reading from '" << fileName
                  << "'.\nInvalid entry: '" << line <<"'");

      // read z value from string stream (only in 3D)
      if( dim_ == 3) {
        sstr >> p[2];
        if(!sstr)
          EXCEPTION("A problem occured while reading from '" << fileName
                    << "'.\nInvalid entry: '" << line <<"'");
      }

      // store vector in current coil struct
      LOG_DBG2(bisa) << "Adding point " << p.ToString();
      coil.points_.Push_back(p);
    }
  }

  void BiotSavart::ComputeBiotSavartField() {
    
    LOG_TRACE(bisa) << "Starting computation of initial field";
    
    // Leave, if normalized field was already calculated
    if( fieldIsComputed_) 
      return;
    
    if( formulation_ == SCALAR_POT ) {
      EXCEPTION( "Scalar potential formulation for Biot-Savart formulation "
          << "is currently not implemented!");
    }
    
    
    UInt numNodes = eqnMap_->GetNumLocalNodes();
    UInt numEqns  = eqnMap_->GetNumEqns();

  

    // loop over all coils
    StdVector<BsCoil>::iterator coilIt = coils_.Begin();
    for( ; coilIt != coils_.End(); coilIt++ ) {
      BsCoil & actCoil = *coilIt;

      //allocate memory
      actCoil.magVecNormalized_.Resize(numEqns);
      actCoil.magVecNormalized_.Init();

      // run over all "nodes of calculation region"
      Vector<Double> pCoordCFS, pCoord(3), magVec(3);
      for ( UInt inode = 1; inode <= numNodes; inode++) {


        //get the global node number
        UInt globNodeNr =   eqnMap_->Pde2MeshNode(inode);

        // get the corresponding coordinates (observer point)
        ptGrid_->GetNodeCoordinate( pCoordCFS, globNodeNr );
        if ( dim_ == 2 ) {
          pCoord[0] = pCoordCFS[0];
          pCoord[1] = pCoordCFS[1];
          pCoord[2] = 0.0;
        }
        else
          pCoord = pCoordCFS;

        try {
        // get contribution to actual point from one coil
        ComputeMagVecNormalized( magVec, pCoord, actCoil );
        } catch ( Exception e ) {
          //RETHROW_EXCEPTION(e, "Could not calculate Biot-Savart solution of coil '"
          //                  << actCoil.name << "' as node #" << globNodeNr
          //                  << " (coords: " << pCoord.ToString() 
          //                  << ") is located within the coil!");
          //
          WARN( "Could not calculate Biot-Savart solution of coil '"
              << actCoil.name << "' as node #" << globNodeNr
              << " (coords: " << pCoord.ToString() 
              << ") is located within the coil!: "
              << e.what() );
                    
        }
        //store normalized magnetic vector potential
        if ( dim_ == 3 ) {
          for ( UInt dof=1; dof<=3; dof++ ) {
            UInt eqn = eqnMap_->GetNodeEqn( globNodeNr, dof);
            if ( eqn > 0 )
              actCoil.magVecNormalized_[eqn-1] += magVec[dof-1];
          }
        }
        else {
          UInt eqn = eqnMap_->GetNodeEqn( globNodeNr, 1);
          if ( eqn > 0)  
            actCoil.magVecNormalized_[eqn-1] += magVec[2];
        }      
      } // loop over all nodes
    } // loop over all coils

    // Resize also "the one" vector potential array
    magVec_.Resize( numEqns );
    
    // set flag
    fieldIsComputed_ = true;
    
    LOG_TRACE(bisa) << "Finished computation of initial field";
  }
  
  void BiotSavart::ComputeMagVecNormalized( Vector<Double>& magVec, 
                                            const Vector<Double>& observer,
                                            const BsCoil& coil) {

    // calculate vectorpotential due to Biot-Savart
    magVec.Resize(3);
    magVec.Init();

    Vector<Double> start(dim_), end(dim_);
    Vector<Double> partMagVec(3);
    UInt numPoints = coil.points_.GetSize();
    for( UInt iSeg = 1; iSeg < numPoints; iSeg++ ) {
      start = coil.points_[iSeg-1];
      end   = coil.points_[iSeg];
      
      if (isAxi_) {
        // calculate vector potential of a curved current stick (circle segment)
        ArcIntegralVecPot( partMagVec,observer,start,end );
      } else {
        // calculate vector potential of a straight current stick
        LOG_DBG3(bisa) << "Calculating field for stick "
                       << start.ToString() << " -> " << end.ToString();  
        LineIntegralVecPot( partMagVec,observer,start,end );
      }
      magVec += partMagVec;
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
    Vector<Double> vec1(3), vec2(3), vec3(3);

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
  
  void BiotSavart::SetTimeFncValue() {
    EXCEPTION( "Not yet implemented" );
  }

} // end of namespace









// ===========================================================================
//   OLD IMPLEMENTATION
// ===========================================================================


//#include <fstream>
//#include <iostream>
//#include <map>
//#include <boost/algorithm/string/trim.hpp>
//#include <boost/algorithm/string.hpp>
//#include <boost/tokenizer.hpp>
////#include <string>
//#include <math.h>
//#include "DataInOut/ParamHandling/ParamNode.hh"
//#include "Domain/grid.hh"// just for the visualization
//
//#include "DataInOut/Logging/cfslog.hh"
////#include "Driver/stdSolveStep.hh"
////#include "Driver/solveStepMagHyst.hh"
//#include "Utils/Coil.hh"
////#include "Utils/SmoothSpline.hh"
////#include "Utils/LinInterpolate.hh"
////#include "Forms/curlfieldop.hh"
//#include "Forms/curlCurlNodeInt.hh"
////#include "Forms/nonConformingInt.hh"
////#include "Forms/nLincurlCurlNodeInt.hh"
////#include "Forms/laplaceInt.hh"
////#include "Forms/linearForm.hh"
////#include "Forms/massInt.hh"
////#include "trapezoidal.hh"
////#include "CoupledPDE/pdecoupling.hh"
////#include "Domain/ansatzFct.hh"
////#include "Driver/assemble.hh"
//#include "Utils/coordSystem.hh"
////#include "Utils/mathParser/mathParser.hh"
//#include "MatVec/vector.hh"
//#include "Domain/elem.hh"
//
//
//#ifdef USE_SCRIPTING
//#include "DataInOut/Scripting/cfsmessenger.hh"
//#endif
//
//
//
//#include "Utils/Bisa.hh" //the head file for Biot-Savart Law
//
//namespace CoupledField {
//DECLARE_LOG(bisa)
//DEFINE_LOG(bisa, "biotSavart")
// // ------------------------------------
//  //  Constructor
//  // ------------------------------------
//  BiotSavart::Bisa(){
//    
//  }
//
//  // ----------------------
//  //   Default Destructor
//  // ----------------------
//  BiotSavart::~Bisa() {
//
//  }
//
//
//  // *********************************
//  //  Read special boundary conitions
//  // *********************************
//
////  void BiotSavart::ReadSpecialBCs(StdVector<std::string>& CoilName,StdVector<std::string>& FileName,StdVector<Double>& Current) {
//
//  void BiotSavart::ReadSpecialBCs(ParamNode * CurrentNode,Grid* currgrid){
//    // --------------------------------------------------------------------
//    //   Get information about coils and open files for measurement coils
//    // --------------------------------------------------------------------
//
//     // get the coordinates 
//		  std::string coilname, tmpfilename;
//		  Double tmpvalue;//reading the current value from the xml file
//		  
//            Bisap_=CurrentNode;
//            gidgrid_=currgrid;
//
//          
//		        // try to get coils from the xml coil part
//	          ParamNode * coilNodes = Bisap_->Get("BiotSavart");
//
//	          ParamNode * model1 = Bisap_->Get("BiotSavart")->Get("BsCoil",false);
//
//	          if(!model1)
//	            return;
//	          else{
//	                     
//	                     StdVector<ParamNode*> Sticks = coilNodes->GetList("BsCoil","name","filename");
//	                     //Sticks.GetSize() get the size of the total elements number of "BsCoil"
//	                     LOG_DBG2(bisa)<<"List total xml Elements Number :"<<Sticks.GetSize();
//
//	                     UInt count=0;
//	                     count=coilNodes->Count("CurrentStick","name","filename");
//	                     LOG_DBG2(bisa)<<"Coils number  :"<<count;
//	                     UInt it=0;
//	                     UInt leap= Sticks.GetSize()/count;
//	                     LOG_DBG2(bisa)<<"leap:(how many elements a coil has)   "<<leap;
//
//	                      for(UInt i = 0; i < Sticks.GetSize(); i+=leap ) {
//										                  Sticks[i]->Get( "name", coilname );
//	                                    Sticks[i]->Get( "filename", tmpfilename);
//	                                    Sticks[i]->Get( "value", tmpvalue);
////	                                    std::cout<<"-------------DEBUG----"<<it+1<<"-start---------"<<std::endl;
//	                                    LOG_DBG2(bisa)<<"coilname:   "<<coilname;
//	                                    LOG_DBG2(bisa)<<"filename:   "<<tmpfilename;
//	                                    LOG_DBG2(bisa)<<"current:   "<<tmpvalue;
//
//	                                    Bisacoil coil_;
//                                                                        
//                                      coil_.cname_ = coilname;
////                                      std::cout<<"-------------DEBUG----"<<coil_.cname_<<"-end 1---------"<<std::endl;
////                                      std::cout<<"-------------DEBUG----"<<it+1<<"-end 1---------"<<std::endl;
//                                      coil_.fname_ = tmpfilename;
////                                      std::cout<<"-------------DEBUG----"<<coil_.fname_<<"-end 2---------"<<std::endl;
////                                      std::cout<<"-------------DEBUG----"<<it+1<<"-end 2---------"<<std::endl;
//                                      coil_.curr_ = tmpvalue;
////                                      std::cout<<"-------------DEBUG----"<<coil_.curr_<<"-end---------"<<std::endl;
////                                      std::cout<<"-------------DEBUG----"<<it+1<<"-end---------"<<std::endl;
////	                                    
//	                                    coils_.Push_back(coil_);  
//	                                    
////	                                    std::cout<<"-------------DEBUG----"<<it+1<<"-end-------"<<std::endl;
////	                                    std::cout<<coils_[it].fname_<<std::endl;
//
//                                      it++;
//
//	                      }             
//	              }
//	         
//	          
//	          for(UInt i=0;i<coils_.GetSize();i++){// read the text file      
//                       coilp_= i;
//	                     ReadFile(coils_[i].fname_);
//	                     if(dim_==3) {
//	                       std::cout<<"start to view"<<std::endl;
//	                       View(coils_[i].fname_);
//	                     }
//	          }
//	          
//	       
//	          
//
//}
//
//
//// template<class TYPE>//just let the debug version could work better
// void BiotSavart::CalcH( Vector<Double>& HField,Vector<Double>& intPoint) {
//
//	 Double inner1,inner2,vnorm2;
//	 Vector<Double> a,b,c;
//	 Vector<Double> v,sum,H;
//   sum.Resize(3,true);//init to zero
////   std::cout<<"how many points do we have in the current coil now  "<<coils_[coilp_].coilnodesno_<<std::endl;
//   for(UInt index=0;index<coils_.GetSize();index++){// iterate 
//               coilp_=index;
//           for(UInt j=0;j<(coils_[coilp_].coilnodesno_-1);j++){
//             
//                      Vector<Double> p1,p2;
//                      p1.Resize(0);
//                      p2.Resize(0);
//                      
//                      p1 = coils_[coilp_].coords[j];
//                      p2 = coils_[coilp_].coords[j+1];
//                               
//                      LOG_DBG2(bisa)<<"Here we are calculate the test case of "<<dim_<<" dimension";
//                      LOG_DBG2(bisa)<<"intPoint before:  "<<intPoint;
//                      
//                      if(dim_==3){
//                        a = p2-p1;
//                        b = p1-intPoint;
//                        c = p2-intPoint;
//                      }
//                      else{                    
//                        if(index==0)intPoint.Push_back(0.0);
//                        a = p2-p1;
//                        b = p1-intPoint;
//                        c = p2-intPoint;
//                      }
//                      LOG_DBG2(bisa)<<"intPoint after:  "<<intPoint;
//                    //  std::cout<<"a:  "<<a<<std::endl;
//                    //  std::cout<<"b:  "<<b<<std::endl;
//                    //  std::cout<<"c:  "<<c<<std::endl;
//                      
//                      c.CrossProduct(a,v);
//                      
//                      a.Inner(c,inner1);
//                      a.Inner(b,inner2);
//                      
//                     //  std::cout<<"coil current point at coil:  "<<(coil_.coilp_+1)<<std::endl;
//                      
//                      vnorm2 = v.NormL2()*v.NormL2();
//                     
//                     // the formular from the book 
//                     // vtmp = v.ScalarDiv(vnorm2);
//                     // H =(coil_.curr_[coil_.coilp_]/(4*M_PI))*(vtmp)*(inner1/c.NormL2()-inner2/b.NormL2());
//                    
//                                        
//                      
//        //            v.ScalarMult(coil_.curr_[coil_.coilp_]);
//                     
//                     v.ScalarMult(coils_[coilp_].curr_);
//                     
//                     H = v /(4*M_PI)/vnorm2*(inner1/c.NormL2()-inner2/b.NormL2());
//                     
//                    //  std::cout<<"H:  "<<H<<std::endl;
//                    //  std::cout<<"sum1:  "<<sum<<std::endl;
//                      sum+=H;
//                               
//                    //  std::cout<<"sum2:  "<<sum<<std::endl;
////                          HField = sum;
//
//                      LOG_DBG2(bisa)<<"HField for coil("<<coilp_+1<<") is :"<<sum;
//                      
//                      if(j==(coils_[coilp_].coilnodesno_-2)){ LOG_DBG2(bisa)<<"-----------the coil end-------------- ";}
//            
//           }
//           
//           HField = sum;
//
//               
//
//     }
//
//
// }
//
//
// void BiotSavart::ReadFile( std::string fileName){
////     std::cout<<"-------------DEBUG--ReadFile---start- 1------"<<std::endl;
////     std::cout<<fileName<<std::endl;
////     std::cout<<fileName.c_str()<<std::endl;
//     // open file
//     std::ifstream file;
//     sampleData.open(fileName.c_str());
//  
////     std::cout<<"-------------DEBUG--ReadFile---start- 2------"<<std::endl;
//     
//     // check if file 'fileName' exists
//     if ( !sampleData ) {
//       EXCEPTION( "Failed to open file '" << fileName
//                  << "' containing data of sample function!" );
//     }
//     sampleData.clear(); // clear flags
//
//     // we don't trust .eof()
//     sampleData.seekg(0,std::ios::end);
//     std::string::size_type pos_end = sampleData.tellg();
//
//     // start from the beginning
//     sampleData.seekg(0,std::ios::beg);
//     char* buf;
//     Double  x, y, z;
//
//     buf = new char[pos_end+1];
//     sampleData.read(buf, pos_end);
//     sampleData.close();
//     buf[pos_end]=0;
//     std::string data(buf);
//     delete[] buf;
//
//     typedef boost::tokenizer<char_separator<char> > Tok;
//     boost::char_separator<char> sep("\n\r");
//     Tok t(data, sep);
//     Tok::iterator it, end;
//     it = t.begin();
//     end = t.end();
//
//     std::string line;
//     std::stringstream sstr;
//
////     std::cout<<"-------------DEBUG--ReadFile---start- 3------"<<std::endl;
//     coils_[coilp_].coilnodesno_= 0;//initilization
//      
//     for(UInt lineNum=1; it != end; it++, lineNum++) {
//       line = *it;
//
//       // strip whitespaces
//       boost::trim(line);
//
//       // big choice of signs for comment's
//       if (line.length() == 0 || line[0] == '#' ||
//           line[0] == '%' || line[0] == '!')
//         continue;
//
//       // put line into a string stream
//       sstr.clear();
//       sstr.str(line);     
//       
//
//       // read x value from string stream
//       sstr >> x;
//       if(!sstr)
//         EXCEPTION("A problem occured while reading from '" << fileName
//                   << "'.\nInvalid entry: '" << line <<"'");
//
//       // read y value from string stream
//       sstr >> y;
//       if(!sstr)
//         EXCEPTION("A problem occured while reading from '" << fileName
//                  << "'.\nInvalid entry: '" << line <<"'");
//
//       sstr >> z;
//       if(!sstr)
//         EXCEPTION("A problem occured while reading from '" << fileName
//             << "'.\nInvalid entry: '" << line <<"'");
//    
//     
//       // store values 
//       
//       //std::cout<<"point"<<coil_.coilnodesno_+1<<": "<<std::endl;
//       //std::cout<<"x :"<<x<<std::endl;
//       //std::cout<<"y :"<<y<<std::endl;
//       //std::cout<<"z :"<<z<<std::endl;
//       
//              
//       Vector<Double> p;
//       p.Resize(0);
//       p.Push_back(x);
//       p.Push_back(y);
//       p.Push_back(z);
//       
////       std::cout<<"-------------DEBUG--ReadFile"<<p<<"-end 1-------"<<std::endl;
//       coils_[coilp_].coords.Push_back(p);
////       std::cout<<"-------------DEBUG--ReadFile"<<p<<"-end 2-------"<<std::endl;
//       coils_[coilp_].coilnodesno_++;
//     
//     }
//   }
// 
// 
// void BiotSavart::View(std::string fileName){
//   // loop over coils
//    
//    // // get number of nodes in the grid
//    
//    // // loop over nodes of coil
//    // // // add new node 
//    
//    // // get coil name
//    // // obtain regionId 
//    
//    // // query number of elements in grid (GetNumElems)
//    // // add nodes-1 elements (AddElems)
//    
//    // // loop over elements ( = nodes of coil -1)
//    // // // set element definition (SetElemData)
//
////   std::cout<<"start to view--but this version will view for each coil file----"<<std::endl;
//
//       for(UInt j=0;j<(coils_[coilp_].coilnodesno_-1);j++){
//       
//       StdVector<UInt> nodeNums;//for the connectivity of gid biot savart coil visualization
//       StdVector<UInt> elemNums;
//  //     std::cout<<"Add the nodes start"<<std::endl;
//       
//        nodeNums.Resize(0);
//        elemNums.Resize(0);
//        Vector<Double> p1,p2;
//        p1.Resize(0);
//        p2.Resize(0);
//        p1 = coils_[coilp_].coords[j];
//        p2 = coils_[coilp_].coords[j+1];     
//  
//        UInt NumNodes=gidgrid_->GetNumNodes();
//  //      std::cout<<"NumNodes before: "<<NumNodes<<std::endl;  
//  //      std::cout<<"add the nodes"<<std::endl;
//        gidgrid_->AddNode(p1,NumNodes);
//        nodeNums.Push_back(gidgrid_->GetNumNodes());
//  //      std::cout<<"NumNodes after: "<<gidgrid_->GetNumNodes()<<std::endl; 
//  //      std::cout<<"add the nodes"<<std::endl;
//        gidgrid_->AddNode(p2,NumNodes);
//        nodeNums.Push_back(gidgrid_->GetNumNodes());
//  //      std::cout<<"NumNodes after: "<<gidgrid_->GetNumNodes()<<std::endl;
//  //      std::cout<<"Add the nodes end"<<std::endl;
//   
//       
//  //     std::cout<<"Add the region start"<<std::endl;
//       UInt NumRegions=gidgrid_->GetNumRegions();
//       fileName.append(j,'1');
//       LOG_DBG2(bisa)<<"NumRegions before:  "<<NumRegions;
//       RegionIdType idg=NumRegions;  
//       gidgrid_->AddRegion(fileName,idg);
//  //     std::cout<<"NumRegions after:  "<<gidgrid_->GetNumRegions()<<std::endl;
//  //     std::cout<<"Add the region end"<<std::endl;     
//      
//  //     std::cout<<"Add the element start"<<std::endl;
//  //     std::cout<<"NumElems before:  "<<NumElems<<std::endl;
//       UInt counter= 1;
//       gidgrid_->AddElems(counter);
//  //     std::cout<<"NumElems after:  "<<gidgrid_->GetNumElems()<<std::endl;
//       UInt * connect= & nodeNums[0];
//                     
//       gidgrid_->SetElemData(gidgrid_->GetNumElems(),Elem::LINE2,idg,connect);
//       fileName.append(j+1,'2');
//       elemNums.Push_back(gidgrid_->GetNumElems());
//      
//       gidgrid_->AddNamedElems(fileName,elemNums);   
//       
//     
//   }
//   
//        
//  
////   std::cout<<"Add the element end"<<std::endl;
//   
////   std::cout<<"end of view"<<std::endl;
//
//  
//   
//
// }
