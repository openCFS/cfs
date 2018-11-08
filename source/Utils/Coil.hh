// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
//=================================
/*
 * \file   Coil.hh
 * \brief  This class is basically a container for all parameters required for
 *         describing the different types of coils used in our coupled field
 *         computations.
 *         The file contains some code from NACS. The code for automatic
 *         calculation of the direction of the current density is obsolete due
 *         to the implementation of the possibility to use current densities
 *         calculated by the electric current PDE. It can be chosen to use the
 *         current density or only its direction, which makes the difference
 *         between ideal (also called stranded) coils and solid conductors.
 *
 * \date   unknown
 * \author ahauck, dperchto
 */
//=================================

#ifndef COIL_FILE_HH
#define COIL_FILE_HH

#include "General/Environment.hh"
#include "MatVec/Vector.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Domain/ElemMapping/EntityLists.hh"
#include "Utils/mathParser/mathParser.hh"

namespace CoupledField {

  // Forward class definition
  class Grid;
  
  class CoordSystem;
  class MathParser;
  class EntityList;

   //! Class for describing coils

   //! This class is basically a container for all parameters required for
   //! describing the different types of coils used in our coupled field
   //! computations.
   //! Besides its container functionality it also provides via its constructor
   //! a method to interface with the XML parameter handling object from which
   //! it obtains the parameter values.
   class Coil {

   public:
     
     //! Type for identifying coils
     typedef std::string IdType;

     //! Enumeration type for distinguishing the different source types of coils
     typedef enum {NO_SOURCE_TYPE, CURRENT, CURRENT_MULTHARM, VOLTAGE, EXTERNAL } SourceType;
     
     //! Constructor for coils

     //! This is the only allowed constructor for the Coil class. It expects
     //! as first input argument the region identifier of the coil. 
     //! As second input argument
     //! it expects the name of the PDE in whose definition the coil appears.
     //! The constructor will use these two symbolic names to query the XML
     //! parameter handling object for the coil's parameters.
     //! If the coil is a measurement coil it will also, if desired, open the
     //! files for storing the current/voltages and the inductivity.
     Coil( PtrParamNode coilNode, 
           PtrParamNode infoNode,
           Grid * ptGrid,
           MathParser * mp,
           Global::ComplexPart type );

     
     /*//! Simplified constructor for coils, just defined by string
     Coil( const std::string& id, Grid * ptGrid );*/
     
     //! Default destructor

     //! The default destructor is responsible for closing the output files
     //! in the case that the object describes a measurement coil.
     ~Coil();
     
     //! Returns all elements of the coil, i.e. collects elements from
     //! all regions of all parts.
     shared_ptr<EntityList> GetElems();

     //! Id of coil
     IdType coilId_;
     
     //! Type of coil excitation
     SourceType sourceType_;
     
     //! Scalar value of excitation (voltage, current)
     PtrCoefFct srcVal_;
     
     //! Values of excitation in different harmonics.
     //! Key is the harmonic (only prescribe positive harmonics,
     //! because the mirroring to the negative spectrum is performed
     //! internally).
     std::map<UInt, PtrCoefFct> srcValMH_;

     //! Flag if we have multiharmonic excitation
     bool isMultHarm_;

     //! Define part of a coil
     struct Part {
       
       //! Constructor
       Part();
       
       //! Destructor
       ~Part();

       //! Regions for part of coil
       StdVector<RegionIdType> regions;

       //! Coefficient function with unit vector for current density direction
       PtrCoefFct jUnitVec;

       /*//! Identifier of input surface region
       StdVector<std::string> inputSurfRegions;
       
       //! Identifier of input surface region
       StdVector<std::string> outputSurfRegions;*/
       
       //! Orientation flag (+/- 1)
       Integer orientFlag;

       //! Coordinate system for current density
       shared_ptr<CoordSystem> coordSys;
       
       //! Contains the source value multiplied by the orientation flag
       
       //! This string contains the source value (current / voltage)
       //! multiplied by the orientation flag (which can be different per 
       //! coil).
       PtrCoefFct sourceVal;
       
       //! Resistance of coil part
       std::string resistance;
       
       //! Number of coil turns
       UInt numTurns;
       
       /*! Cross section of coil
       Double coilCrossSect;*/
       
       //! Cross section of wire
       Double wireCrossSect;
       
       /*! Fill factor kappa
       Double fillFactor;*/
       
       /*! Assume uniform current density
       bool uniformCurrentDens_;
       
       //! Flag if current density direction is analytical
       bool hasAnalyticalDir_;*/

     private:
       //! Prevent usage of copy constructor
       Part(const Part&){ 
         EXCEPTION("Not allowed");
       }
       
     };

     //! Part for each regionId
     std::map<RegionIdType, shared_ptr<Part> > parts_;
     
     //! Contains parts having direction from external simulation
     std::map<shared_ptr<Part>, PtrParamNode > partsExtJDir_;

   private:

     //! The default constructor is not allowed
     Coil();

     //! The copy constructor is not allowed
     Coil( const Coil &c );

     /* Helper method for setting up a free coordinate system
     shared_ptr<CoordSystem>
     SetupCoosy( UInt iPart, 
                 const Part& actPart,
                 StdVector<RegionIdType>& regions,
                 const StdVector<std::string>& inSurfaces,
                 const StdVector<std::string>& outSurfaces );*/
                 
     //! Parameter node
     PtrParamNode myParam_;
     
     //! Info node
     PtrParamNode myInfo_;
     
     //! Pointer to grid
     Grid * ptGrid_;
     
     //! Math parser instance
     MathParser * mParser_;
   };
}

#endif
