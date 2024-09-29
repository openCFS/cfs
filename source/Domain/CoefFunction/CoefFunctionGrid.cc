// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*- 
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     CoefFunctionGrid.cc 
 *       \brief    Implementation file for the Grid interpolation CoefFunction
 *
 *       \date     Jan. 23, 2012
 *       \author   Andreas Hueppe
 */
//================================================================================================


/*
 * Folgende Konzepte sollen verwendet werden:
 *  1. Normale Interpolation:
 *     Es gibt einen zusaetzlichen tag im XML ueber den der user eine Ordnung einstellen kann.
 *     Die CoefFunction wird immer LagrangeElemente verwenden. Dei Ordnung richtet sich dann dabei
 *     nach der Ansatzordnung des Zielgitters.
 *     Der Vorteil liegt darin, dass man nur noch die DOF des Zielgitters speichert und dennoch schnell an
 *     das Element bzw dessen loesung kommt bei GetScalar.
 *     Ist das Zielgitter mit GRID order gegeben, werden die Knotenfreiheitsgrade als zielknoten verwendet
 *     Ist das Zielgitter mit einer (anistropen) ordnung gegeben wird versucht auf Knoten einen entsprechenden
 *     Lagrange elements der gegebenen Ordnung zu interpolieren.
 *     Folgender Ablauf:
 *     1. Beschaffe die Koordinaten des Zielgitters
 *        a. Wenn Knoteninterpolation gewuenscht, nimm einfach die Gitterknoten
 *        b. Bei hoeheren ordungen versuchen die DOF Koordinaten zu beschaffen. (Es werden nur nodal Spaces erstellt!)
 *     2. Hole elementliste des Quellgitters
 *        Wir sind hierbei auf Knotenergebnisse fixiert, benutzen also die Referenzelemente und Ansatzfunktionen des H1NodalExpl
 *     3. Vielleicht kann man sogar noch die Moeglichkeit schaffen direkt Ableitungsoperatoren an die Coeffunction zu geben
 *     4. Starten der Vorberechnung, welche Knoten des Zeilgitters liegen in welchem Element des Quellgitters
 *     5. Interpolation der Quellergebnisse auf die Knotenfreiheitsgrade des Zielgitters
 *     6. In GetScalar/Vector/Tensor
 *        a. Hol dir das Elementergebnis zum lpm
 *        b. Bereche Operator am LPM zu den gegebenen Ansatzfunktionen
 *        c. Extrahiere Elementergebnis und mutipliziere mit Operator.
 *
 *     Offene Fragen:
 *      - Wie sind die Referenzelemente zu erstellen. Brauch ich nen Space?
 *        -> Eigentlich braucht man keinen space, man kann einfach Knotenelemente nehmen
 *           Nochmal Map aus elemType auf referenzelement
 *      - Wie erhalte ich die Koordinaten der Freiheitsgrade bei hoeheren Ordnungen
 *        -> Speizalisierte Funktion fuer die H1LagExpl und H1LagVar Elemente
 *
 *  2. Conservativ:
 *     Vorberechnungen
 *     1. Hole Knotenkoordinaten des Quellgitters
 *     2. Hole Elementlisten des Zielgitters
 *     3. Fuehre Punktsuche aus
 *     4. Fuer jedes Element des Zielgitters:
 *        a. Hole die Gleichungsnummern
 *        b. Fuer jeden Punkt innerhalb des Zielelements
 *           I.  Bereche Ansatzfunktion an lokalem Punkt
 *           II. Fuege entsprechende Daten in CoordMatrix ein
 *        c. Convert CoordMatrix->CRS
 *
 *     InterpolateConservative:
 *     1. Matrix Vektor Multiplikation...
 *
 *     GetScalar/Vektor/Tensor
 *     -> Exception macht keinen Sinn hierbei. Man muss MapToFeFunction benutzen
 *
 *  3. NoInterpolation:
 *     GetScalar/Vector/Tensor:
 *      1. Hole elementnummer
 *      2. Extrahiere Elementloesung
 *      3. Inteproliere mit Operator
 *
 *  4. MapToCoefFunction mit Coordinaten.
 *     Das kann vom Space oder von FeFunction aufgerufen werden
 *     Konzeptionell besser waere von Space aber mal sehen
 *
 */

#include "CoefFunctionGrid.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "FeBasis/FeSpace.hh"
#include "DataInOut/ResultHandler.hh"
#include <boost/lexical_cast.hpp>

//include headers of subclasses for factory method
#include "CoefFunctionGridNodalInterp.hh"
#include "CoefFunctionGridNodalDefault.hh"
#include "CoefFunctionGridNodalSource.hh"
//element based version
#include "CoefFunctionGridElemDefault.hh"
//#include "CoefFunctionGridHigherDefault.hh"
//#include "CoefFunctionGridHigherInterp.hh"
#include "Driver/BaseDriver.hh"


namespace CoupledField{
    
  
PtrCoefFct CoefFunctionGrid::Generate( Domain* ptDomain,
                                       Global::ComplexPart format, 
                                       PtrParamNode infoNode, 
                                       PtrParamNode configNode,
									   shared_ptr<RegionList> regions,
									   ResultInfo::EntryType type
									   ){
  shared_ptr<CoefFunctionGrid> ret;
  PtrParamNode tmpNode  =  infoNode->Get("externalData");

  // verify the solution type (nodal/element data)

  // parse the quantity and check if it is known
  std::string definedOnString;
  configNode->Get("defaultGrid")->GetValue("definedOn",definedOnString);
  
  if (definedOnString=="node") {
    if(configNode->Has("defaultGrid")) {
      std::string dependString;
      PtrParamNode tmpNodeType  =  configNode->Get("defaultGrid");
      tmpNodeType->GetValue("dependtype",dependString);

      if ( dependString == "INVSOURCE" || dependString == "INVMEASURE") {
        if(format == Global::COMPLEX){
          ret.reset(new CoefFunctionGridNodalSource<Complex>(ptDomain,
              configNode->Get("defaultGrid"), tmpNode, regions));
        }
        else{
          EXCEPTION("CoefFunctionGridNodalSource can just be complex");
        }
      }
      else {
        if (format == Global::COMPLEX) {
          ret.reset(new CoefFunctionGridNodalDefault<Complex>(ptDomain,
                configNode->Get("defaultGrid"), tmpNode, regions,type));
        }
        else {
          ret.reset(new CoefFunctionGridNodalDefault<Double>(ptDomain,
                      configNode->Get("defaultGrid"), tmpNode, regions,type));
        }
      }

    }
    else if(configNode->Has("externalGrid")){
      if(format == Global::COMPLEX){
        ret.reset(new CoefFunctionGridNodalInterp<Complex>(ptDomain,
                      configNode->Get("externalGrid", ParamNode::INSERT), tmpNode, regions));
      }
      else {
        ret.reset(new CoefFunctionGridNodalInterp<Double>(ptDomain,
                      configNode->Get("externalGrid", ParamNode::INSERT), tmpNode, regions));
        ret->PrepareInterpolation();
      }
    } else {
      EXCEPTION("CoefFunctionGrid generator called with invalid xml tag. This is a serious bug, please report!");
    }
  } else if (definedOnString=="element") {
    if(configNode->Has("defaultGrid")) {
      std::string dependString;
      PtrParamNode tmpNodeType  =  configNode->Get("defaultGrid");
      tmpNodeType->GetValue("dependtype",dependString);

      if ( dependString == "INVSOURCE" || dependString == "INVMEASURE") {
        EXCEPTION("CoefFunctionGridElemSource is not implemented yet!");
      }
      else {
        if (format == Global::COMPLEX) {
          ret.reset(new CoefFunctionGridElemDefault<Complex>(ptDomain,
                configNode->Get("defaultGrid"), tmpNode, regions,type));
        }
        else {
          ret.reset(new CoefFunctionGridElemDefault<Double>(ptDomain,
                      configNode->Get("defaultGrid"), tmpNode, regions,type));
        }
      }

    }
    else if(configNode->Has("externalGrid")){
      EXCEPTION("CoefFunctionGridElemInterp is not implemented yet!");
    } else {
      EXCEPTION("CoefFunctionGrid generator called with invalid xml tag. This is a serious bug, please report!");
    }
  }                

  return ret;
}

CoefFunctionGrid::CoefFunctionGrid(Domain* ptDomain, PtrParamNode configNode, shared_ptr<RegionList> regions){
  dimDof_ = 0;
  inputId_ = "";
  gridId_ = "";
  srcGrid_ = NULL;
  domain_ = ptDomain;
  solType_ = NO_SOLUTION_TYPE;
  aSeqStep_ = 0;
  curInterpType_ = NO_INTERPOLATION;
  myConfigNode_ = configNode;
  verbose_ = false;
  //obtain the sequence step for result
  //plz note we take here the same value as for the current simulation
  //TODO: Parametrize by XML
  //check if the usseqenceSteper specifies it
  if(configNode->Get("sequenceStep")){
    this->aSeqStep_  = configNode->Get("sequenceStep")->As<UInt>();
  }else{
    this->aSeqStep_ = domain_->GetDriver()->GetActSequenceStep();
    WARN("external data did not specify its sequence step. assuming current step...")
  }
  SetEntitiesByRegions(regions);
}

CoefFunctionGrid::~CoefFunctionGrid(){

}

void CoefFunctionGrid::GetTensor(Matrix<Double>& CoefMat,
                                 const LocPointMapped& lpm ){
  EXCEPTION("CoefFunctionGrid: GetTensor is not implemented here"); 

}
void CoefFunctionGrid::GetTensor(Matrix<Complex>& CoefMat,
                                 const LocPointMapped& lpm ){
  EXCEPTION("CoefFunctionGrid: GetTensor is not implemented here"); 

}

void CoefFunctionGrid::GetVector(Vector<Complex>& CoefMat,
                                            const LocPointMapped& lpm ){
  EXCEPTION("GetVector is not implemented here"); 
}

void CoefFunctionGrid::GetVector(Vector<Double>& CoefMat,
                                            const LocPointMapped& lpm ){
  EXCEPTION("GetVector is not implemented here"); 
//=======
//  //1. call to grid to find the element
//  //2. compute solution of this element
//  //3. apply operator
//  assert(dimType_ != TENSOR);
//  if(dimType_ == SCALAR)
//    CoefMat.Resize(1);
//  else if (dimType_ == VECTOR)
//    CoefMat.Resize(domain_->GetGrid(gridId_)->GetDim());
//
//  CoefMat.Init();
//  //const Elem* targetElem = lpm.ptEl;
//  Vector<Double> globCoord;
//  Vector<Double> localCoord;
//  Vector<DATA_TYPE> elemSol;
//
//  lpm.shapeMap->Local2Global(globCoord,lpm.lp);
//  LocPoint lp;
//  const Elem* sourceElem = domain_->GetGrid(gridId_)->
//      GetElemAtGlobalCoord(globCoord,lp,srcRegions_);
//
//  if(!sourceElem){
//    WARN("Could not find element for coord " << globCoord);
//    return;
//  }
//
//  shared_ptr<ElemShapeMap> esm = domain->GetGrid(gridId_)->GetElemShapeMap( sourceElem, true );
//  LocPoint myLp = localCoord;
//  LocPointMapped lpmS;
//  lpmS.Set(myLp,esm,lpm.weight);
//
//  feFunct_->GetElemSolution( elemSol, sourceElem);
//  BaseFE * ptFe = feFunct_->GetFeSpace()->GetFe(sourceElem->elemNum);
//  myOperator_->ApplyOp(CoefMat,lpmS,ptFe,elemSol);
//>>>>>>> .r12248
}

void CoefFunctionGrid::GetScalar(Complex& CoefMat,
                                           const LocPointMapped& lpm ){
  EXCEPTION("GetScalar is not implemented here"); 
//
//  assert(dimType_ == SCALAR );
//
//  //const Elem* targetElem = lpm.ptEl;
//  Vector<Double> globCoord;
//  Vector<Double> localCoord;
//  Vector<DATA_TYPE> elemSol;
//  Vector<DATA_TYPE> ptSol(1);
//  ptSol.Init();
//
//  lpm.shapeMap->Local2Global(globCoord,lpm.lp);
//  StdVector< Vector<Double> > locCoords;
//  LocPoint lp;
//  const Elem* sourceElem = domain_->
//      GetGrid(gridId_)->GetElemAtGlobalCoord(globCoord,lp,srcRegions_);
//  if(!sourceElem){
//    return;
//  }
//
//  shared_ptr<ElemShapeMap> esm = domain_->GetGrid(gridId_)->GetElemShapeMap( sourceElem, true );
//  LocPoint myLp = localCoord;
//  LocPointMapped lpmS;
//  lpmS.Set(myLp,esm,lpm.weight);
//  feFunct_->GetElemSolution( elemSol, sourceElem);
//  BaseFE * ptFe = feFunct_->GetFeSpace()->GetFe(sourceElem->elemNum);
//  myOperator_->ApplyOp(ptSol,lpmS,ptFe,elemSol);
//  CoefMat = ptSol[0];
//>>>>>>> .r12248
}
void CoefFunctionGrid::GetScalar(Double& CoefMat,
                                           const LocPointMapped& lpm ){
  EXCEPTION("GetScalar is not implemented here"); 
}



UInt CoefFunctionGrid::GetVecSize() const {
  return dimDof_;
}

void CoefFunctionGrid::GetTensorSize( UInt& numRows, UInt& numCols ) const {
  numRows = numRowTensor_;
  numCols = numColTensor_;
  //EXCEPTION("No tensor valued data available yet");
}


void CoefFunctionGrid::DetermineResult(std::string inputID,UInt seqStep){
  //obtain availResults and search for the requested one
  StdVector<shared_ptr<ResultInfo> > results;
  domain_->GetResultHandler()->GetResultTypes(inputID,seqStep,results,false);
  //now we search for the appropriate result
  for(UInt i = 0;i<results.GetSize();i++){
	  if( results[i]->resultType == solType_ ) {
      resultInfo_ = results[i];
	    break;
	  }
  }
  std::string solString = SolutionTypeEnum.ToString(solType_);
  if(!resultInfo_.get())
	  EXCEPTION("Can not find result " << solString << " in given file!");
}


}

// Definition of coefficient function dependency type
static EnumTuple coefInterpTuples[] =
{
 EnumTuple(CoefFunctionGrid::STANDARD, "standard"),
 EnumTuple(CoefFunctionGrid::CONSERVATIVE, "conservative"),
 EnumTuple(CoefFunctionGrid::NO_INTERPOLATION,  "none")
};

Enum<CoefFunctionGrid::InterpType>CoefFunctionGrid::InterpType_ = \
    Enum<CoefFunctionGrid::InterpType>("Interpolation type to be used",
                                       sizeof(coefInterpTuples) / sizeof(EnumTuple),
                                       coefInterpTuples);


