#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>

#include "grid_struct.hh"

#include <Domain/grid.hh>
#include <Elements/elements_header.hh>
#include <DataInOut/ParamHandling/BaseParamHandler.hh>
#include "DataInOut/GMV/outGMV.hh"
#include "DataInOut/WriteInfo.hh"

namespace CoupledField
{

  template<UInt DIM>
  GridStruct<DIM> :: GridStruct(UInt adim)
  {
    ENTER_FCN( "GridStruct::GridStruct" );

    dim_     = adim;
    pdename_ = "acoustic";

    isInitialized_ = FALSE;

    /*
      DEFINED IN BASEPARAMHANDLER.hh && XMLPARAMHANDLER.hh
      The method will try to find the specified keyword in the parameter tree and will
      return the values of the respective elements or attributes in a vector of Doubles.
      The search can be restricted to certain subtrees by specifying keywords for section
      and subsection. The methjod will return an empty vector, if there is no match at all.
      It will issue an error message, if there are matches for both, elements and attributes.
      "name" = keyword
      sd_    = Vector of Doubles (output)
      "domain" = Name of a section in which to look for keyword (optional)
      "region" = Name of a subsection in which to look for keyword (optional)
    */
    params->GetList( "name", regionNames_, "domain", "region" );
    if (regionNames_.GetSize() != 1) {
        Info->Error( "Slicing techique just supports one region!!", 
                     __FILE__, __LINE__ );
    }

    // check for absorbing boundary conditions
    UInt numRegions;
    StdVector<std::string> abcBCs;
    params->GetList( "name", abcBCs, pdename_, "absorbingBCs" );

    if ( abcBCs.GetSize() == 0 ) {
      numRegions = 1;
    }
    else if ( abcBCs.GetSize() == 1 ) {
      numRegions = 2;
      regionNames_.Push_back(abcBCs[0]);
    }
    else {
      Error("Just one ABC-region allowed",__FILE__,__LINE__);
    }

    //std::cout << "reg1=" << regionNames_[0] << "  reg2=" << regionNames_[1] << std::endl;

    //we just allow one volume and one surface region
    //therefore, we can hardcode the Ids level
    volRegionIds_.Resize(1);
    volRegionIds_[0] = 0;

    surfRegionIds_.Resize(1);
    surfRegionIds_[0] = 1;

    //Get slice data
    StdVector<std::string> keyVec;

    //The information about the Dimension itself, should hopefully be
    //provided by the grid-sturct class itself

    keyVec  = pdename_, "sliceData", "dimension", "dimX";
    params->Get( keyVec, dimX_);
    keyVec  = pdename_, "sliceData", "dimension", "dimY";
    params->Get( keyVec, dimY_);
    keyVec  = pdename_, "sliceData", "dimension", "dimZ";
    params->Get( keyVec, dimZ_);

    keyVec  = pdename_, "sliceData", "waveLength";
    params->Get( keyVec, waveLength_);

    keyVec  = pdename_, "sliceData", "elementsPerWavelength";
    params->Get( keyVec, elementsPerWavelength_);

    keyVec  = pdename_, "sliceData", "pulseTime";
    params->Get( keyVec, pulseTime_);

    keyVec  = pdename_, "sliceData", "soundSpeed";
    params->Get( keyVec, soundSpeed_);

    keyVec  = pdename_, "sliceData", "safetyRegion";
    params->Get( keyVec, safetyRegion_);

    keyVec  = pdename_, "sliceData", "addNrOfWaveLength";
    params->Get( keyVec, tol_);

//     keyVec  = pdename_, "region", "name";
//     params->Get( keyVec, subname_);

  }


  template<UInt DIM>
  void GridStruct<DIM> :: Read()
  {
    ENTER_FCN( "GridStruct::Read" );

    // now we construct our own grid
    GenGridStruct(1,3,20);
    SetNamedNodes();

    isInitialized_ = TRUE;
  }


  template<UInt DIM>
  void GridStruct<DIM> :: GenGridStruct(const UInt elemx, const UInt elemy,
					const UInt elemz)
  {
    ENTER_FCN( "GridStruct::GenGridStruct");

    maxnumelemz_ = elemz;	//Number of elements in z-direction (moving- direction)
    maxnumelemy_ = elemy;	//Number of elements in y-direction
    maxnumelemx_ = elemx;	//Number of elements in x-direction

    Integer innodes_;		//Number of nodes per element
    Integer i,j,k,l,m,temp,temp2,temp3;		//index- variables
    Double meshsizez_ = 1.5E-4;	//
    Double meshsizey_ = 1.5E-4;	//
    Double meshsizex_ = 0.0;	//
    Double x = 0.0;			//Coordinates
    Double y = 0.0;			//Coordinates
    Double z = 0.0;			//Coordinates
    Double start_=0.0;
    Double totalLength;

    Integer absbcs=0;


    switch(dim_) {

    case 2:
      //calculate number of elements for x- and y-direction
      meshsizez_   = waveLength_/elementsPerWavelength_;
      maxnumelemy_ = (Integer) (dimY_/meshsizez_);

      totalLength  = 2*pulseTime_*soundSpeed_ + (tol_ * waveLength_);
      //	           + (safetyRegion_ * waveLength_);
      maxnumelemz_ = (Integer) ( totalLength / meshsizez_ + 10*meshsizez_ ) ;

      std::cout << std::endl;
      std::cout << "Size of slice in  z- direction:  " << maxnumelemz_* meshsizez_ << std:: endl;
      std::cout << "Size of slice in  y- direction:  " << maxnumelemy_* meshsizez_ << std:: endl;
      std::cout << "elements in  z- direction:  " << maxnumelemz_ << std:: endl;
      std::cout << "elements in  y- direction:  " << maxnumelemy_ << std:: endl;
      std::cout << "elements in  x- direction:  " << maxnumelemx_ << std:: endl;
      break;

    case 3:

      meshsizez_   = waveLength_/elementsPerWavelength_;
      meshsizex_   = meshsizez_;
      meshsizey_   = meshsizez_;
      maxnumelemy_ = (Integer) (dimY_/meshsizez_);
      maxnumelemx_ = (Integer) (dimX_/meshsizez_);

      totalLength  = 2*pulseTime_*soundSpeed_ + (tol_ * waveLength_);
	//	           + (safetyRegion_ * waveLength_);
      maxnumelemz_ = (Integer) ( totalLength / meshsizez_ + 10*meshsizez_ );
      break;
    }

    //Calculate maximum number of nodes and elements
    std::cout << dim_ << "Dimensional Problem" << std:: endl;
    switch(dim_) {

    case 2:
      innodes_ = 4;
      numNodes_    = (maxnumelemz_ + 1) * (maxnumelemy_ + 1);
      numVolElems_ = maxnumelemz_ * maxnumelemy_;
      break;

    case 3:
      innodes_ = 8;
      numNodes_    = (maxnumelemz_ + 1) * (maxnumelemy_+1) * (maxnumelemx_+1);
      numVolElems_ = maxnumelemz_ * maxnumelemy_ * maxnumelemx_;
      std::cout << "Maxnumnodes:" << numNodes_ << std:: endl;
      std::cout << "Maxnumelements:" << numVolElems_ << std:: endl;
      std::cout << "Maxnumelements in z-direction:" << maxnumelemz_ << std:: endl;
      std::cout << "Maxnumelements in x_direction:" << maxnumelemx_ << std:: endl;
		
      break;
    }

    //allocate for coordinates
    ptCoordinate_ = new Point<DIM> [numNodes_];

    //allocate for one region 
    volElems_.Resize(1);  
    surfElems_.Resize(1); 

    //absorbing boundary conditions
    Boolean ABC = TRUE;

    //Init
    switch(dim_) {

    case 2:

      k=0;

      for(i=0;i<maxnumelemz_+1;i++) {
        for(j=0;j<maxnumelemy_+1;j++) {

          ptCoordinate_[k][0] = start_ + z;
          ptCoordinate_[k][1] = start_ + y;
          y+= meshsizey_;
          k++;
        }
        z+=meshsizez_;
        y=0.0;
      }

      k = 0;
      m=1;

      for(i=0;i<maxnumelemz_;i++) {
        for(j=0;j<maxnumelemy_;j++) {

          Elem * el= new Elem();
          el->elemNum  = m;
          el->ptElem   = ptQ1;
          el->regionId = volRegionIds_[0] ;
          el->connect.Resize(innodes_);

          el->connect[0] = 1+ k;		 	     //Connect elements
          el->connect[1] = 1+ k +     (maxnumelemy_+1);//Connect elements
          el->connect[2] = 1+ k + 1 + (maxnumelemy_+1);//Connect elements
          el->connect[3] = 1+ k + 1;		     //Connect elements

          volElems_[0].Push_back(el);

          k++;
          m++;
        }
        k++;
      }

      //absorbing BCs
      //lower
      temp = m;

      if ( ABC ) {
	Elem * volEl = volElems_[0][0];
        for(i=0; i< maxnumelemz_;i++) {
          SurfElem* el= new SurfElem();
          el->elemNum = temp;
          el->ptElem = ptL1;
          el->regionId = surfRegionIds_[0];
          el->connect.Resize(2);
          
          el->connect[0] = 1+i*maxnumelemy_;
          el->connect[1] = 1+(i+1)*maxnumelemy_;

	  //set arbitrary volume element
          el->ptVolElem1 = volEl;

          surfElems_[0].Push_back(el);
          
          temp++;
        }
        
        //upper
        temp2 = temp;
        for(i=0; i< maxnumelemz_;i++) {
          SurfElem* el= new SurfElem();
          el->elemNum = temp2;
          el->ptElem = ptL1;
          el->regionId = surfRegionIds_[0];
          el->connect.Resize(2);
          
          el->connect[0] = maxnumelemy_ +i*maxnumelemy_;
          el->connect[1] = maxnumelemy_ +(i+1)*maxnumelemy_;
          
	  //set arbitrary volume element
          el->ptVolElem1 = volEl;

          surfElems_[0].Push_back(el);
          
          temp2++;
        }
        
        absbcs = maxnumelemy_-((maxnumelemy_+1)/2);
        //left
        temp3 = temp2;
        for(i=0; i< absbcs;i++) {
          SurfElem* el= new SurfElem();
          el->elemNum = temp3;
          el->ptElem = ptL1;
          el->regionId = surfRegionIds_[0];
          el->connect.Resize(2);
          
          el->connect[0] = (maxnumelemy_+1)/2+i;
          el->connect[1] = (maxnumelemy_+1)/2+1+i;
          
	  //set arbitrary volume element
          el->ptVolElem1 = volEl;

          surfElems_[0].Push_back(el);
          
          temp3++;
        }
      }

      break;
	
	
      // 3-dimensional grid generation
    case 3: 
		
      k=0;
      //calculate the coordinates of the points
      for(i=0;i<maxnumelemz_+1;i++) {
        for(j=0;j<maxnumelemx_+1;j++) {
          for(l=0;l<maxnumelemy_+1;l++) {

            ptCoordinate_[k][0] = start_ + z;
            ptCoordinate_[k][1] = start_ + y;
            ptCoordinate_[k][2] = start_ + x;
            y+= meshsizey_;
            k++;
          }
          x+=meshsizex_;
          y=0.0;
        }
        z+=meshsizez_;
        x=0.0;
      }
      k=0;
      m=1;
      std::cout << maxnumelemz_ << " " << maxnumelemx_ << " " <<  maxnumelemy_ << std::endl;
      //defining, what point belong to what element
      for(l=0;l<maxnumelemz_;l++) {
        for(i=0; i<maxnumelemx_; i++){
          for(j=0; j<maxnumelemy_; j++){
            Elem* el = new Elem();
            el -> elemNum  = m;
            el -> ptElem   = ptHexa1;
            el->  regionId = volRegionIds_[0] ;
            el -> connect.Resize(innodes_);

            //Connect elements
            //std::cerr << "Element Nr.: " << k << ": ";
            el->connect[0] = j+i*(maxnumelemy_+1)+l*(maxnumelemy_+1)*(maxnumelemx_+1)+1;
            //std::cerr <<  j+i*(maxnumelemy_+1)+l*(maxnumelemy_+1)*(maxnumelemx_+1)<< ", ";
            el->connect[3] = j+i*(maxnumelemy_+1)+(l+1)*(maxnumelemy_+1)*(maxnumelemx_+1)+1;
            //std::cerr <<   j+i*(maxnumelemy_+1)+(l+1)*(maxnumelemy_+1)*(maxnumelemx_+1)<< ", ";
            el->connect[2] =j+(i+1)*(maxnumelemy_+1)+(l+1)*(maxnumelemy_+1)*(maxnumelemx_+1)+1;
            //std::cerr <<  j+(i+1)*(maxnumelemy_+1)+(l+1)*(maxnumelemy_+1)*(maxnumelemx_+1) << ", ";
            el->connect[1] = j+(i+1)*(maxnumelemy_+1)+l*(maxnumelemy_+1)*(maxnumelemx_+1)+1;
            //std::cerr <<  j+(i+1)*(maxnumelemy_+1)+l*(maxnumelemy_+1)*(maxnumelemx_+1) << ", ";
            el->connect[4] = (j+1)+i*(maxnumelemy_+1)+l*(maxnumelemy_+1)*(maxnumelemx_+1)+1;
            //std::cerr << (j+1)+i*(maxnumelemy_+1)+l*(maxnumelemy_+1)*(maxnumelemx_+1) << ", ";
            el->connect[7] = (j+1)+i*(maxnumelemy_+1)+(l+1)*(maxnumelemy_+1)*(maxnumelemx_+1)+1;
            //std::cerr <<   (j+1)+i*(maxnumelemy_+1)+(l+1)*(maxnumelemy_+1)*(maxnumelemx_+1)<< ", ";
            el->connect[6] = (j+1)+(i+1)*(maxnumelemy_+1)+(l+1)*(maxnumelemy_+1)*(maxnumelemx_+1)+1;
            //std::cerr <<   (j+1)+(i+1)*(maxnumelemy_+1)+(l+1)*(maxnumelemy_+1)*(maxnumelemx_+1)<< ", ";
            el->connect[5] = (j+1)+(i+1)*(maxnumelemy_+1)+l*(maxnumelemy_+1)*(maxnumelemx_+1)+1;
            //std::cerr <<   (j+1)+(i+1)*(maxnumelemy_+1)+l*(maxnumelemy_+1)*(maxnumelemx_+1)<< std::endl;
            //add the Element el to the Element-Vector

            volElems_[0].Push_back(el);
            
            k++;
            m++;
          }
        }
      }

      /*******************************************************************
        ABSORBING BOUNDARY CONDITIONS HANDLING
        adding the surface planes to socalled boundary elements
      ********************************************************************/
      // referring to the TOP-BOUNDARY

      if ( ABC ) {
	Elem * volEl = volElems_[0][0];
        temp = 0;
        for(i=0;i<maxnumelemz_;i++) {
          for(j=0;j<maxnumelemx_;j++) {
            SurfElem* el = new SurfElem();
            el->elemNum = m;
            el->ptElem = ptQ1;
            el->regionId = surfRegionIds_[0];
            el->connect.Resize(4);
            
            el->connect[0] = maxnumelemy_+(maxnumelemy_+1)*j+(maxnumelemy_+1)*(maxnumelemx_+1)*i+1;
            //std::cout << maxnumelemy_+(maxnumelemy_+1)*j+(maxnumelemy_+1)*(maxnumelemx_+1)*i+1 << ", ";
            el->connect[1] =  maxnumelemy_+(maxnumelemy_+1)*j+(maxnumelemy_+1)*(maxnumelemx_+1)*(i+1)+1; 
            //std::cout <<   maxnumelemy_+(maxnumelemy_+1)*j+(maxnumelemy_+1)*(maxnumelemx_+1)*(i+1)+1 <<", ";
            el->connect[2] =  maxnumelemy_+(maxnumelemy_+1)*(j+1)+(maxnumelemy_+1)*(maxnumelemx_+1)*(i+1)+1; 
            //std::cout <<  maxnumelemy_+(maxnumelemy_+1)*(j+1)+(maxnumelemy_+1)*(maxnumelemx_+1)*(i+1)+1<< ", ";
            el->connect[3] =  maxnumelemy_+(maxnumelemy_+1)*(j+1)+(maxnumelemy_+1)*(maxnumelemx_+1)*i+1;
            //std::cout << maxnumelemy_+(maxnumelemy_+1)*(j+1)+(maxnumelemy_+1)*(maxnumelemx_+1)*i+1<< std::endl;
            
	    //set arbitrary volume element
	    el->ptVolElem1 = volEl;

            surfElems_[0].Push_back(el);
            
            m++;
            temp++;
          }
        }
        //std::cerr << "soviele top BE: " << temp << std::endl;
        // referring to the BOTTOM-BOUNDARY
        
       //  k=0;
//         temp = 0;
//         for(i=0;i<maxnumelemz_;i++) {
//           for(j=0;j<maxnumelemx_;j++) {
//             SurfElem* el = new SurfElem();
//             el->elemNum = m;
//             el->ptElem = ptQ1;
//             el->regionId = surfRegionIds_[0];
//             el->connect.Resize(4);
            
//             el->connect[0] = (maxnumelemy_+1)*j + (maxnumelemy_+1)*(maxnumelemx_+1)*i+1;
//             //std::cout << (maxnumelemy_+1)*j + (maxnumelemy_+1)*(maxnumelemx_+1)*i<< ", ";
//             el->connect[1] = (maxnumelemy_+1)*j + (maxnumelemy_+1)*(maxnumelemx_+1)*(i+1)+1;
//             //std::cout << (maxnumelemy_+1)*j + (maxnumelemy_+1)*(maxnumelemx_+1)*(i+1)<< ", ";
//             el->connect[2] = (maxnumelemy_+1)*(j+1) + (maxnumelemy_+1)*(maxnumelemx_+1)*(i+1)+1;
//             //std::cout << (maxnumelemy_+1)*(j+1) + (maxnumelemy_+1)*(maxnumelemx_+1)*(i+1)<< ", ";
//             el->connect[3] = (maxnumelemy_+1)*(j+1) + (maxnumelemy_+1)*(maxnumelemx_+1)*i+1;
//             //std::cout << (maxnumelemy_+1)*(j+1) + (maxnumelemy_+1)*(maxnumelemx_+1)*i<< std::endl;
            
// 	    //set arbitrary volume element
// 	    el->ptVolElem1 = volEl;

// 	    surfElems_[0].Push_back(el);

//             k++;
//             m++;
//             temp++;
//           }
//         }
        //std::cerr << "soviele bottom  BE: " << temp << std::endl;
        // referring to the RIGHT-BOUNDARY
       //  k=0;
//         temp = 0;
//         for(i=0;i<maxnumelemz_;i++) {
//           for(j=0;j<maxnumelemy_;j++) {
//             SurfElem* el = new SurfElem();
//             el->elemNum = m;
//             el->ptElem = ptQ1;
//             el->regionId = surfRegionIds_[0];
//             el->connect.Resize(4);
            
//             el->connect[0] = j + (maxnumelemx_+1)*(maxnumelemy_+1)*i+1;
//             //std::cout << j + (maxnumelemx_+1)*(maxnumelemy_+1)*i << ", ";
//             el->connect[1] = j + (maxnumelemx_+1)*(maxnumelemy_+1)*(i+1)+1;
//             //std::cout <<  j + (maxnumelemx_+1)*(maxnumelemy_+1)*(i+1)<< ", ";
//             el->connect[2] = (j+1) + (maxnumelemx_+1)*(maxnumelemy_+1)*(i+1)+1;
//             //std::cout <<  (j+1) + (maxnumelemx_+1)*(maxnumelemy_+1)*(i+1)<< ", ";
//             el->connect[3] = (j+1) + (maxnumelemx_+1)*(maxnumelemy_+1)*i+1;
//             //std::cout <<  (j+1) + (maxnumelemx_+1)*(maxnumelemy_+1)*i<< std::endl;

// 	    //set arbitrary volume element
// 	    el->ptVolElem1 = volEl;  
          
//             surfElems_[0].Push_back(el);

//             k++;
//             m++;
//             temp++;
//           }
//         }
        //std::cerr << "soviele right BE: " << temp << std::endl;
        // referring to the LEFT_BOUNDARY
        
        k=1;
        temp = 0;
        for(i=0;i<maxnumelemz_;i++) {
          for(j=0;j<maxnumelemy_;j++) {
            SurfElem* el = new SurfElem();
            el->elemNum = m;
            el->ptElem = ptQ1;
            el->regionId = surfRegionIds_[0];
            el->connect.Resize(4);
            
            el->connect[0] = (maxnumelemy_+1)*maxnumelemx_ + j + (maxnumelemy_+1)*(maxnumelemx_+1)*i+1;
            //std::cout << (maxnumelemy_+1)*maxnumelemx_ + j + (maxnumelemy_+1)*(maxnumelemx_+1)*i+1<< ", ";
            el->connect[1] = (maxnumelemy_+1)*maxnumelemx_ + j + (maxnumelemy_+1)*(maxnumelemx_+1)*(i+1)+1;
            //std::cout << (maxnumelemy_+1)*maxnumelemx_ + j + (maxnumelemy_+1)*(maxnumelemx_+1)*(i+1)+1<< ", ";
            el->connect[2] = (maxnumelemy_+1)*maxnumelemx_ + (j+1) + (maxnumelemy_+1)*(maxnumelemx_+1)*(i+1)+1;
            //std::cout << (maxnumelemy_+1)*maxnumelemx_ + (j+1) + (maxnumelemy_+1)*(maxnumelemx_+1)*(i+1)+1<< ", ";
            el->connect[3] = (maxnumelemy_+1)*maxnumelemx_ + (j+1) + (maxnumelemy_+1)*(maxnumelemx_+1)*i+1;
            //std::cout <<  (maxnumelemy_+1)*maxnumelemx_ + (j+1) + (maxnumelemy_+1)*(maxnumelemx_+1)*i+1<< std::endl;

	    //set arbitrary volume element
	    el->ptVolElem1 = volEl;  
          
            surfElems_[0].Push_back(el);

            k++;
            m++;
            temp++;
          }
        }
	//referring to the BOUNDARY AT THE END OF THE SLICE
	temp = 0;
	double shift = (maxnumelemy_+1)*(maxnumelemx_+1);
        for(i=0;i<maxnumelemy_;i++) {
          for(j=0;j<maxnumelemx_;j++) {
            SurfElem* el = new SurfElem();
            el->elemNum = m;
            el->ptElem = ptQ1;
            el->regionId = surfRegionIds_[0];
            el->connect.Resize(4);
	    el->connect[0] =  maxnumelemz_*shift + (maxnumelemy_+1)*j + i + 1;
	    //std::cout << maxnumelemz_*shift + (maxnumelemy_+1)*j + i + 1 << ", ";
            el->connect[1] =  maxnumelemz_*shift + (maxnumelemy_+1)*j + i + 2;
	    //std::cout << maxnumelemz_*shift + (maxnumelemy_+1)*j + i + 2<< ", ";
            el->connect[2] =  maxnumelemz_*shift + (maxnumelemy_+1)*(j+1) + i + 2;
	    //std::cout <<maxnumelemz_*shift + (maxnumelemy_+1)*(j+1) + i + 2 << ", ";
            el->connect[3] =  maxnumelemz_*shift + (maxnumelemy_+1)*(j+1) + i + 1; 
            //std::cout <<maxnumelemz_*shift + (maxnumelemy_+1)*(j+1) + i + 1 << std::endl;
            
	    //set arbitrary volume element
	    el->ptVolElem1 = volEl;

            surfElems_[0].Push_back(el);
            
            m++;
            temp++;
          }
        }	
      }

      break;
    }
  }


  template<UInt DIM>
  void GridStruct<DIM>::SetNamedNodes() {

   ENTER_FCN( "GridStruct::SetNamedNodes");

   //check xml-file
    // vectors for parameter handling
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    StdVector<std::string> bcs_id;

    // Get names of node sets, values and filenames for inhomogenous
    // Dirichlet boundary conditions
    keyVec  = pdename_, "bcsAndLoads", "dirichletInhom", "name";
    attrVec = "", "tag", "";
    valVec  = "", "", "";
    params->GetList(keyVec, attrVec, valVec, bcs_id);

    if ( bcs_id.GetSize() == 0 ) {
      Error("GridStruct: You have to specify excat one inh. Dirichlet BCs", 
	    __FILE__,__LINE__);
    }
    else if ( bcs_id.GetSize() > 1 ) {
      Error("GridStruct: You can just specify one inh. Dirichlet BC", 
	    __FILE__,__LINE__);
    }

 
    //inhomogeneous Dirichlet boundary conditions
    namedNodes_.Resize(1);
    namedNodeNames_.Resize(1);   
    namedNodeNames_[0] = bcs_id[0];

    UInt numbc;
    
    switch(dim_) {
    case 2:
      numbc = (UInt)(maxnumelemy_+1);
      break;
      
    case 3:
      numbc = (UInt) ( (maxnumelemy_+1) * (maxnumelemx_+1) );
      break;
    }
    
    namedNodes_[0].Resize(numbc);
    for ( UInt i=0; i<numbc; i++) {
      namedNodes_[0][i] = i+1;
    }
    
  }



  /****************************************************************************
  Transforms the grid, like described in SSOUTSCHEK Thesis:
	- shift the last part of the Grid, to the front of the 
	  in built grid
	- write the knew coordinates of the Grid in ptCoordinate_ (Point-
	  class)
  ****************************************************************************/
  template<UInt DIM>
 //  void GridStruct<DIM> :: TransformGridStruct(UInt& nodeShift, UInt & shiftFactor, 
// 					      const UInt flag)
  void GridStruct<DIM> :: TransformGridStruct(UInt & shiftFactor, UInt & nodeShift,
		UInt & elemgrid, Double &  meshsize, const UInt flag)
  {
    ENTER_FCN( "GridStruct::TransformGridStruct");
    //Variablen
    //Numshift*2 + 2*sik = SizeOfGrid
    Integer innodes_;
    Integer numelemshift_, maxnumelem;	
    Integer shiftfactor_;	//
    Double a,b,x,y,z,shiftLength, sliceLength;
    Double meshsizez_;	//
    Double meshsizey_;	//
    Double meshsize_;	// it should be enough to have one variable meshsize,
                        // because the grid we are using is structured.
    Integer i,j,k,m;		//Index
    Double start_;

    k=0;
    x=0.0;
    y=0.0;
    z=0.0;

    shiftLength   = pulseTime_*soundSpeed_ + (tol_ * waveLength_);
    shiftLength  -= safetyRegion_ * waveLength_;
    meshsizez_    = waveLength_/elementsPerWavelength_;
    numelemshift_ = (Integer) ( shiftLength / meshsizez_ );

   
    //calculate shiftfaktor
    switch(dim_) {
    case 2 :
      shiftfactor_ = maxnumelemy_+1;
      innodes_ = 4;
      //Calculate meshsize
      a = ptCoordinate_[0][1];
      b = ptCoordinate_[1][1];
      //meshsizey_ = b-a;
      break;
    case 3 :
      //stimmt?
      // shiftfactor is "one slice" of the grid, so i*shiftfactor
      // shifts i slices in z-direction
      shiftfactor_ = (maxnumelemy_+1)*(maxnumelemx_+1);
      innodes_ = 8;
      //Calculate meshsizes for x,y,z direction
      break;
    }
   
 //If flag = 1 the function TransformGridStruct is writing
 //information for the saveNodes-Routine
    if(flag){
      sliceLength = 2*pulseTime_*soundSpeed_ + (tol_ * waveLength_);
      maxnumelem = (Integer) (sliceLength / meshsizez_);
      elemgrid    = maxnumelem;
      shiftFactor = shiftfactor_;
      meshsize    = meshsizez_;
      return;
    }
    Info->StartProgress(" \n Shift Executed \n");
    std::cout << "Number of elements shifted in z-direction  - gridstruct = " 
	      << numelemshift_ << std::endl;
    


    // switch(dim_) {

//     case 2 :
//       //Shift of known values
//       for(i=0; i<(((numelemshift_) + 1 + (maxnumelemz_%2))*(maxnumelemy_+1));i++) {
//         ptCoordinate_[i][0] = ptCoordinate_[i + (numelemshift_*shiftfactor_) ][0];
//       }

//       //calc meshsize
//       a = ptCoordinate_[0][0];
//       b = ptCoordinate_[1+shiftfactor_][0];
//       meshsizez_ = b-a;


//       start_ = ptCoordinate_[(((numelemshift_) + 1 
// 			       + (maxnumelemz_%2))*(maxnumelemy_+1))-1][0] + meshsizez_;

//       k = ((numelemshift_)+(maxnumelemz_%2)+1)*(maxnumelemy_+1);

//       for(i=0;i<numelemshift_;i++) {
//         for(j=0;j<maxnumelemy_+1;j++) {

//           ptCoordinate_[k][0] = start_ + z;
//           ptCoordinate_[k][1] = y;
//           y+= meshsizey_;
//           k++;
//         }
//         z+=meshsizez_;
//         y=0.0;
//       }

//       break;

//     case 3:
//       std::cerr << "Transformation of Slice 3D" << std::endl;
//       //Shift of known values
//       //the only values that change, are the  one´s referring to the 
//       //z-direction, the rest remains!
//       for(i=0; i<(((numelemshift_)+(maxnumelemz_%2))*shiftfactor_);i++) {
//         ptCoordinate_[i][0] = ptCoordinate_[i + (numelemshift_*shiftfactor_)][0];
//       }
//       start_ = ptCoordinate_[(((numelemshift_)+(maxnumelemz_%2))*shiftfactor_)-1][0] + meshsize_;
//       k = ((numelemshift_)+(maxnumelemz_%2))*shiftfactor_;
//       std::cout << "k= " << k << std::endl;
//       for(i=0;i<numelemshift_;i++) {
//         for(j=0;j<maxnumelemx_+1;j++) {
//           for(m=0;m<maxnumelemy_+1;m++) {

//             ptCoordinate_[k][0] = start_ + z;
//             ptCoordinate_[k][1] =  y;
//             ptCoordinate_[k][2] =  x;
//             y+= meshsize_;
//             k++;
//           }
//           x+=meshsize_;
//           y=0.0;
//         }
//         z+=meshsize_;
//         x=0.0;
//       }
//     }
    //std::cout << "k= " << k << std::endl;
    //std::cout << "z at point: " << ptCoordinate_[2500][0] << std::endl;
    //std::cout << "z at point: " << ptCoordinate_[3375][0] << std::endl;
    //std::cout << "z at point: " << ptCoordinate_[3380][0] << std::endl;
    //std::cout << "z at point: " << ptCoordinate_[3355][0] << std::endl;
    //std::cout << "z at point: " << ptCoordinate_[0][0] << std::endl;
	
    nodeShift = numelemshift_*shiftfactor_;
    std::cout << "nodeshift = " << nodeShift << std::endl;
  }



  //===================================================================================================================================

  template<UInt DIM>
  Integer GridStruct<DIM> :: GetMaxElem(std::string dir)
  {
    ENTER_FCN( "GridStruct::GetMaxElem");

    Integer var;
    var = 0;
    if(dir=="x") {
      var= maxnumelemx_;
      //std::cout << "x= " << var << std::endl;
      //      return var;
    } else
      if(dir=="y") {
        var= maxnumelemy_;
        //std::cout << "y= " << var << std::endl;
	//        return var;
      } else if(dir=="z") {
        var= maxnumelemz_;
        //std::cout << "z= " << var << std::endl;
	//        return var;
      }

    return var;
  }


  template<UInt DIM>
  void GridStruct<DIM> :: GetAllRegionNames( StdVector<std::string> & regionNames ) {

    ENTER_FCN( "GridStruct::GetMaxElem");

    regionNames = regionNames_;

//     params->GetList( "name", regionNames, "domain", "region" );

//     if ( regionNames.GetSize() > 1) {
//       Info->Error( "Slicing techique just supports one region!!", 
//                      __FILE__, __LINE__ );
//     }


  }
  

  template<UInt DIM>
  void GridStruct<DIM> :: GetRegionIds( StdVector<RegionIdType> & regions ) {
    ENTER_FCN( "GridStruct::GetRegionIds");
    regions.Clear();
    regions.Reserve(volRegionIds_.GetSize() + surfRegionIds_.GetSize());
    
    for( UInt i = 0; i < volRegionIds_.GetSize(); i++) {
      regions.Push_back(volRegionIds_[i]);
    }

    for( UInt i = 0; i < surfRegionIds_.GetSize(); i++) {
      regions.Push_back(surfRegionIds_[i]);
    }
  }

  // ======================================================
  // NODE ACCESS FUNCTIONS
  // ======================================================
  template<UInt DIM>
  void GridStruct<DIM>::GetNodesByName( StdVector<UInt> & nodeList,
					const std::string & name ) {
    ENTER_FCN( "GridStruct::GetNodesByName" );

    Integer index = namedNodeNames_.Find(name);
    if ( index != -1 ) {
      nodeList = namedNodes_[index];
    } else {
      (*error) << "GridStruct: There are no nodes with name'" << name
               << "' in the grid!";
      Error( __FILE__, __LINE__ );
    }
    
  }
  
  template<UInt DIM>
  UInt GridStruct<DIM>::GetNumNodes( const std::string & nodesName ) {
    ENTER_FCN( "GridCFS::GetNumNodes" );

    UInt numNodes = 0;
    Integer index = namedNodeNames_.Find(nodesName);

    if ( index != -1 ) {
      numNodes = namedNodes_[index].GetSize();
    } else {
      (*error) << "GridCFS: The Nodes with name '" << nodesName
               << "' were not found in the grid!";
      Error( __FILE__, __LINE__ );
    }

    return numNodes;
  }


  template<UInt DIM>
  void GridStruct<DIM>::GetNodesByRegion( StdVector<UInt> & nodeList,
                                       const RegionIdType regionId ) {
    ENTER_FCN( "GridStruct::GetNodesByRegion" );

    Integer index = 0;
    
    // look in volume regions
    index = volRegionIds_.Find(regionId);
    if ( index != -1 ) {
      nodeList = volElemNodes_[index];
    } else {
      
      // look in surface regions
      index = surfRegionIds_.Find(regionId);
      if ( index != -1 ) {
        nodeList  = surfElemNodes_[index];
      } else {
        (*error) << "GridStruct: The region with id '" << regionId
                 << "' was not found in the grid!";
        Error( __FILE__, __LINE__ );
      }
    }
    
  }
  
  template<UInt DIM>
  void GridStruct<DIM>::GetNodeCoordinate( Point<DIM> & rfPoint,
                                        const UInt inode ) {
    ENTER_FCN( "GridStruct::GetNodeCoordinate" );
    
    if ( inode > numNodes_ || inode < 0 ) {
      (*error) << "GridStruct: There are only " << numNodes_
               << " nodes in the grid. You requested coordinates for "
               << "node number " << inode <<". Go check your mesh file!";
      Error( __FILE__, __LINE__ );
    }
    
    rfPoint = ptCoordinate_[inode-1];
  }
  
  template<UInt DIM>
  void GridStruct<DIM>::GetNodeCoordinate( Vector<Double> & rfPoint,
                                        const UInt inode ) {
    ENTER_FCN( "GridStruct::GetNodeCoordinate" );

    if ( inode > numNodes_ || inode < 0 ) {
      (*error) << "GridStruct: There are only " << numNodes_
               << " nodes in the grid. You requested coordinates for "
               << "node number " << inode <<". Go check your mesh file!";
      Error( __FILE__, __LINE__ );
    }
    rfPoint = ptCoordinate_[inode-1];
  }
    
  // ======================================================
  // ELEMENT ACCESS FUNCTIONS
  // ======================================================
  template<UInt DIM>
  const Elem * GridStruct<DIM>::GetElem( UInt elemNr ) {
    ENTER_FCN( "GridStruct::GetElem" );
    
    //!!!!!!!!!!!!!!!! currently just volume elements!!!!
    
    if ( elemNr > numVolElems_ || elemNr < 0) {  
      (*error) << "GridStruct: There are only " << numVolElems_ 
               << " elements in the grid! You requested element number " 
               << elemNr << ". Go check your mesh file!";
      Error( __FILE__, __LINE__ );
    }
    
    return volElems_[0][elemNr-1];

  }

  template<UInt DIM>
  void GridStruct<DIM>::GetElems( StdVector<Elem*> & elems, 
                                  const RegionIdType regionId ) {
    ENTER_FCN( "GridStruct::GetElems" );
    elems.Clear();

    // check if region Id is ALL_REGIONS
    if ( regionId == ALL_REGIONS ) {
      elems.Reserve( GetNumVolElems() + GetNumSurfElems() );

      // iterate over all volume elements
      for ( UInt i = 0; i < volElems_.GetSize(); i++) {
        for (UInt iElem = 0; iElem < volElems_[i].GetSize(); iElem++ ) {
          elems.Push_back(volElems_[i][iElem]);
        }
      }
      // iterate over all surface elements
      for ( UInt i = 0; i < surfElems_.GetSize(); i++) {
        for (UInt iElem = 0; iElem < surfElems_[i].GetSize(); iElem++ ) {
          elems.Push_back(surfElems_[i][iElem]);
        }
      }
      
    } else {
      // look in volume regions
      Integer index = volRegionIds_.Find(regionId);
      if ( index != -1 ) {
        elems = volElems_[index];
      } else {    
        // look in surface regions
        index = surfRegionIds_.Find(regionId);
        if ( index != -1 ) {
          elems.Reserve( surfElems_[index].GetSize());
          for (UInt iElem=0; iElem<surfElems_[index].GetSize(); iElem++ ) {
            elems.Push_back(surfElems_[index][iElem]);
          }
        } else {    
          (*error) << "GridStruct: The region with id '" << regionId
                   << "' was not found in the grid!";
          Error( __FILE__, __LINE__ );
        }
      }
    }
  }


  template<UInt DIM>
  void GridStruct<DIM>::GetVolElems( StdVector<Elem*> & elems, 
                                  const RegionIdType regionId ) {
    ENTER_FCN( "GridStruct::GetVolElems" );
    
    // check if region Id is ALL_REGIONS
    if ( regionId == ALL_REGIONS ) {
      elems.Reserve( GetNumVolElems() );
      for ( UInt i = 0; i < volElems_.GetSize(); i++) {
        for (UInt iElem = 0; iElem < volElems_[i].GetSize(); iElem++ ) {
          elems.Push_back(volElems_[i][iElem]);
        }
      }
    } else {
      Integer index = volRegionIds_.Find(regionId);
      if ( index != -1 ) {
        elems = volElems_[index];
      } else {    
        (*error) << "GridStruct: The volume region with id '" << regionId
                 << "' was not found in the grid!";
        Error( __FILE__, __LINE__ );
      }
    }
  }
  
  template<UInt DIM>
  void GridStruct<DIM>::GetSurfElems( StdVector<SurfElem*> & elems, 
                                   const RegionIdType regionId ) {
    ENTER_FCN( "GridStruct::GetSurfElems" );
    
    Integer index = surfRegionIds_.Find(regionId);
    if ( index != -1 ) {
      elems = surfElems_[index];
    } else {    
      (*error) << "GridStruct: The surface region with id '" << regionId
               << "' was not found in the grid!";
      Error( __FILE__, __LINE__ );
    }
  }

  template<UInt DIM>
  void GridStruct<DIM>::GetElemsByName( StdVector<Elem*> & elems,
                                     const std::string & elemsName ) {
    ENTER_FCN( "GridStruct::GetElemsByName" );
    Error( "Not implemented", __FILE__, __LINE__ );
  }

  template<UInt DIM>
  void GridStruct<DIM>::GetElemNodes( StdVector<UInt> & connect, 
                                   const UInt iElem ) {
    ENTER_FCN( "GridStruct::GetElemNodes" );
    

    //!!!!!!!!!!!!!!!! currently just volume elements!!!!

    if ( iElem > numVolElems_ || iElem < 0) {  
      (*error) << "GridStruct: There are only " << numVolElems_ 
               << " elements in the grid! You requested element number " 
               << iElem << ". Go check your mesh file!";
      Error( __FILE__, __LINE__ );
    }
    
    connect = volElems_[0][iElem-1]->connect;
    
  }

  template<UInt DIM>
  void GridStruct<DIM>::GetElemNodesCoord( Matrix<Double> & coordMat,  
                                        const StdVector<UInt> & connect ) {
    ENTER_FCN( "GridStruct::GetElemNodesCoord" );

    coordMat.Resize(dim_, connect.GetSize());
    for (UInt k=0; k < connect.GetSize(); k++)    
      for (UInt actDim=0; actDim < dim_; actDim++)
        coordMat[actDim][k] = ptCoordinate_[connect[k]-1][actDim];
  }

  template<UInt DIM>
  void GridStruct<DIM>::CalcSurfNormal( Vector<Double> & n, 
                                     const Elem & surfElem ) {
    ENTER_FCN( "GridStruct::CalcSurfNormal" );
    
    //compute normal vector
    if (dim_ ==2 ) {
      n.Resize(2);
      n[0] = 0.0;
      n[1] = 1.0;
    }
    else {
      n.Resize(3);
      n[0] = 0.0;
      n[1] = 1.0;
      n[2] = 0.0;
    }
  }

#ifdef __sgi
#pragma instantiate GridStruct<3>
#pragma instantiate GridStruct<2>
#endif
} // end namespace


