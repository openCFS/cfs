#include <assert.h>
#include <cmath>
#include <map>
#include <ostream>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Domain/Domain.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "Domain/ElemMapping/EntityLists.hh"
#include "Domain/Mesh/Grid.hh"
#include "Domain/ElemMapping/SurfElem.hh"
#include "FeBasis/FeSpace.hh"
#include "Driver/Assemble.hh"
#include "Driver/BaseDriver.hh"
#include "Driver/FormsContexts.hh"
#include "Driver/TransientDriver.hh"
#include "FeBasis/BaseFE.hh"
#include "Forms/LinForms/LinearForm.hh"
#include "General/Environment.hh"
#include "MatVec/SingleVector.hh"
#include "MatVec/exprt/xpr1.hh"
#include "MatVec/exprt/xpr2.hh"
#include "MatVec/Vector.hh"
#include "Optimization/Condition.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/ShapeDesign.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Excitation.hh"
#include "Optimization/Function.hh"
#include "Optimization/Objective.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/ShapeOpt.hh"
#include "PDE/SinglePDE.hh"
#include "Utils/StdVector.hh"
#include "Utils/mathParser/mathParser.hh"

namespace CoupledField {

DEFINE_LOG(ShOpt, "shapeOpt")

ShapeOpt::ShapeOpt() : ParamMat()
{
  shapedesign = dynamic_cast<ShapeDesign*>(design);

  PtrParamNode sopn = pn->Get("shapeOpt");
  shapedesign->Configure(sopn, objectives.data.GetSize(), constraints.view->GetNumberOfActiveConstraints());
  exoprt_fe_design_ = shapedesign->AlsoMatOpt();

  Assemble* assemble = context->pde->GetAssemble();

  // all (bi)linear forms need to use updated coordinates
  for (const auto& [_, ctxts] : assemble->GetBiLinForms())
    for(BiLinFormContext* ctxt : ctxts)
      ctxt->GetIntegrator()->SetCoordUpdate(true);

  // set the linearForms used in multiple excitations, note that this does contain all linear forms (some even several times)
  for(unsigned int i = 0; i < me->excitations.GetSize(); i++)
  {
    StdVector<LinearFormContext*>& linForms = me->excitations[i].forms;
    for(unsigned int j = 0; j < linForms.GetSize(); ++j)
      linForms[j]->GetIntegrator()->SetCoordUpdate(true);
  }
}

double ShapeOpt::CalcVolume(Objective* c, Condition* g, bool derivative, bool normalized){
  // the exponent is used in Ersatzmaterial for the volume cost function
  // if an exponent != 1.0 at this point makes any sense is unknown
/*  Function* f = Function::GetFunction(c, g);

  if(derivative)
  {
    StdVector<double> der; // solution
    Matrix<double> CornerCoords;
    Matrix<double> J;
    Matrix<double> iJ;
    Matrix<double> dCornerCoords;
    Matrix<double> dJ;
    Matrix<double> diJ;

    int np = shapedesign->GetNumberOfAuxParameters();
    der.Resize(np, 0);
    if(exoprt_fe_design_){
      // this needs to be done before, we do use fraction
      ErsatzMaterial::CalcVolume(c, g, derivative, normalized);
    }
    if(!exoprt_fe_design_ || (g && g->GetDesignType() == DesignElement::UNITY))
    {
      if(!normalized)
      {
        // this is independent of material optimization, simply the derivative of the real volume
        Grid* grd = domain->GetGrid();
        StdVector<RegionIdType> regs;
        grd->GetVolRegionIds(regs);
        for(unsigned int ri = 0; ri < regs.GetSize(); ri++) // loop regions
        {
          RegionIdType rid = regs[ri];
          if(!g || g->IsForRegion(rid))
          {
            StdVector<Elem*> elems;
            grd->GetElems(elems,rid);
            for( UInt i = 0; i < elems.GetSize(); i++ ) // loop elements
            {
              const Elem* elem = elems[i];
              if(shapedesign->IsElemDependentAtAll(elem->connect))
              {
                // if this element does not depend on any parameters, we can simply skip all the calculations
                shared_ptr<ElemShapeMap> esm = grd->GetElemShapeMap(elem, true); // updated ?!
                shared_ptr<BaseFeFunction> fe = context->pde->GetFeFunction(context->pde->GetNativeSolutionType());
                shared_ptr<IntScheme> intScheme = fe->GetFeSpace()->GetIntScheme();
                // Obtain FE element from feSpace and integration scheme
                ElemList elemList(grd);
                elemList.SetElement(elem);
                IntegOrder order;
                IntScheme::IntegMethod method;
                BaseFE* ptFe = fe->GetFeSpace()->GetFe(elemList.GetIterator(), method, order );
                // Get integration points
                StdVector<LocPoint> intPoints;
                StdVector<Double> weights;
                intScheme->GetIntPoints(Elem::GetShapeType(elem->type), IntScheme::GAUSS, order, intPoints, weights);
                // Loop over all integration points

                LocPointMapped lpm;
                for(unsigned int i = 0; i < intPoints.GetSize(); i++) // loop integration points
                {
                  lpm.Set(intPoints[i], esm, weights[i]);
                  lpm.jac;
                  lpm.jacInv;

                  double w = lpm.jacDet * lpm.weight;

                  // this is for the derived corner coords

                  shared_ptr<ElemShapeMap> der_esm = grd->GetElemShapeMap(elem, true); // updated ?!
                  LocPointMapped der_lpm;

                  for(int p = 0; p < np; p++) // loop over all parameters
                  {
                    // sets dCornerCoords
                    if(shapedesign->GetElemNodesCoordDerivative(dCornerCoords, elem->connect, p))
                    { // returns false if dCornerCoords == 0
                      der_esm->SetElem(dCornerCoords, elem);
                      der_lpm.Set(intPoints[i], esm, weights[i]);

                      ptelem->CalcJacobianAtIp(dJ, ip, dCornerCoords, elem);
                      diJ.Resize(dimJ, dimJ);
                      iJ.Mult(dJ, diJ); // diJ = iJ * dJ;
                      // double tr = diJ.Trace(); // tr = trace(iJ*dJ) = trace(dJ*iJ)
                      der[p] += w * diJ.Trace();
                    } // if dCornerCoords
                  }
                }



                BaseFE* ptelem = elem->ptElem;
                grd->GetElemNodesCoord(CornerCoords, elem->connect, true );
                const int nip = ptelem->GetNumIntPoints();
                const Vector<Double> & intWeights = ptelem->GetIntWeights();

                for(int ip = 1; ip <= nip; ip++){ // loop over all integration points
                  ptelem->CalcJacobianAtIp(J, ip, CornerCoords, elem);
                  const int dimJ(J.GetNumRows());
                  J.Invert(iJ);
                  double det;
                  J.Determinant(det);
                  double w = det * intWeights[ip-1];
                  for(int p = 0; p < np; p++) { // loop over all parameters
                    if(shapedesign->GetElemNodesCoordDerivative(dCornerCoords, elem->connect, p)){ // returns false if dCornerCoords == 0 
                      ptelem->CalcJacobianAtIp(dJ, ip, dCornerCoords, elem);
                      diJ.Resize(dimJ, dimJ);
                      iJ.Mult(dJ, diJ); // diJ = iJ * dJ;
                      // double tr = diJ.Trace(); // tr = trace(iJ*dJ) = trace(dJ*iJ)
                      der[p] += w * diJ.Trace();
                    } // if dCornerCoords
                  } // params
                } // int.points
              }
            } // elems
          } // if region
        } // region
      } // normalized
    }else{ // volume is average or sum of design
      // this is similar to ErsatzMaterial::CalcVolume but calculates derivatives w.r.t. shape
      Grid* grd = domain->GetGrid();
      bool isObjective = g == NULL;
      double fraction = isObjective ? volumeFraction_ : g->volume_fraction; // this already considers everything
      double volume = 0.0;
      if(!normalized){  // needed for derivative in normalized versions
        volume = CalcVolume(c, g, false, normalized);
      }
      bool allDesignsRelevant = g == NULL || g->GetDesignType() == DesignElement::TENSOR_TRACE || g->GetDesignType() == DesignElement::DEFAULT;
      bool ersatzMaterialTensor = domain->HasErsatzMaterialTensor() && allDesignsRelevant;
      unsigned int upper = ersatzMaterialTensor ? design->GetNumberOfElements() : design->data.GetSize();
      Matrix<double> material;
      for(unsigned int i = 0; i < upper; i++) {
        DesignElement* de = &design->data[i];
        bool relevant = (allDesignsRelevant || g->GetDesignType() == de->GetType()) && (isObjective || g->IsForRegion(de->elem->regionId));
        const Elem* elem = de->elem;
        if(relevant && shapedesign->IsElemDependentAtAll(elem->connect)){
          double des;
          if(ersatzMaterialTensor){ // use the trace of the stiffness Tensor as "volume"
            design->GetErsatzMaterialTensor(material, pde->GetSubTensorType(), de->elem, DesignElement::NO_DERIVATIVE, f->GetNotation());
            des = material.Trace();
          }else{
            des = de->GetDesign(DesignElement::PLAIN);
          }
          BaseFE* ptelem = elem->ptElem;
          grd->GetElemNodesCoord(CornerCoords, elem->connect, true );
          const int nip = ptelem->GetNumIntPoints();
          const Vector<Double> & intWeights = ptelem->GetIntWeights();
          for(int ip = 1; ip <= nip; ip++){ // loop over all integration points
            double intWeight = intWeights[ip-1];
            ptelem->CalcJacobianAtIp(J, ip, CornerCoords, elem);
            const int dimJ(J.GetNumRows());
            J.Invert(iJ);
            double det;
            J.Determinant(det);
            for(int p = 0; p < np; p++) { // loop over all parameters
              if(shapedesign->GetElemNodesCoordDerivative(dCornerCoords, elem->connect, p)){ // returns false if dCornerCoords == 0
                ptelem->CalcJacobianAtIp(dJ, ip, dCornerCoords, elem);
                diJ.Resize(dimJ, dimJ);
                iJ.Mult(dJ, diJ); // diJ = iJ * dJ;
                double tr = diJ.Trace();
                der[p] += fraction * intWeight * tr * det * (des - volume); // fraction * intWeight * dArea/dalpha * (d - volume)    ( fraction = 1 und volume = 0 in non-normalized version)
              }
            } // params
          } // int.points
        }
      }
    }
    // derivative in direction of our parameters is always needed and calculated here
    // derivative in direction of design-element parameters only if on that regions
    // these derivatives are independent of our parameters and can be calculated as before
    shapedesign->AddAuxDerivatives(c, g, der, 1.0);
  }else{ // derivative
    if(!exoprt_fe_design_ || (g && g->GetDesignType() == DesignElement::UNITY)){ // this is the real volume, not multiplied by design, we also use this if no design available
      // if design is unity, we use the grid instead of designspace
      if(normalized) return(1.0);
      Grid* grd = domain->GetGrid();
      StdVector<RegionIdType> regs;
      grd->GetVolRegionIds(regs);
      double s = 0.0;
      for(unsigned int i = 0; i < regs.GetSize(); i++){
        int rid = regs[i];
        if(!g || g->IsForRegion(rid)){
          s += grd->CalcVolumeOfRegion(rid, true);
        }
      }
      return(s);
    }else{ // working on a design, exoprt_fe_design_ must be true
      return(ErsatzMaterial::CalcVolume(c, g, derivative, normalized));
    }
  } */
  return 0.0;
}

void ShapeOpt::CalcMinusU1dKU2(StateContainer& forward, StateContainer& adjoint, Objective* f, Condition* constraint, const Matrix<double>* tensor_diff)
{
  /*
  StdVector<double> der; // solution
  int np = shapedesign->GetNumberOfAuxParameters();
  der.Resize(np, 0.0);
  const bool homogenization = tensor_diff != NULL;
  const unsigned int ex_size(me->excitations.GetSize());
  double rcubevol(1.0);
  const bool transient = IsTransient();
  if(homogenization){
    rcubevol = 1.0 / grid->CalcVolumeSpannedByNamedNodes();    
  }
  UInt timesteps(domain->GetDriver()->GetNumSteps());
  double dt = 1.0, gamma = 1.0, beta = 1.0;
  if(transient){
    dt = dynamic_cast<TransientDriver*>(domain->GetDriver())->GetDeltaT();
    gamma = pde->getTimeStepping()->GetNewmarkGamma();
    beta = pde->getTimeStepping()->GetNewmarkBeta();
  }
  
  MathParser* parser = domain->GetMathParser();
  unsigned int mathParserHandle = parser->GetNewHandle();
  
  Grid* grd = domain->GetGrid();

  Matrix<double> CornerCoords;
  Matrix<double> D;
  Matrix<double> J;
  Matrix<double> iJ;
  Matrix<double> dPhi;
  Matrix<double> B;
  Matrix<double> BD;
  Matrix<double> BDB;
  Matrix<double> dCornerCoords;
  Matrix<double> dJ;
  Matrix<double> A1;
  Matrix<double> A2;
  Matrix<double> A3;
  Matrix<double> A4;
  Matrix<double> dK;
  Matrix<double> M;
  Vector<double> dKp;
  Vector<double> dMp;
  Vector<double> dDp;
  Vector<double> shapeFncAtIp;
  Matrix<double> partM;
  
  Matrix<double> tmp_strain(dim, dim); // homogenization
  Matrix<double> tmp_displacement;
  Vector<double> u1diff;
  Vector<double> u2diff;
  
  LOG_DBG2(ShOpt) << "np=" << np << ", ex_size=" << ex_size << ", timesteps=" << timesteps 
      << ", dt=" << dt << ", gamma=" << gamma << ", beta=" << beta;
  
  //caching (see CalcNewmarkDerivative for why)
  StdVector<StdVector<StdVector<SingleVector*>* > > forwards;
  StdVector<StdVector<StdVector<SingleVector*>* > > forwarddt;
  StdVector<StdVector<StdVector<SingleVector*>* > > forwarddtt;
  StdVector<StdVector<StdVector<SingleVector*>* > > adjoints;
  forwards.Resize(ex_size);
  adjoints.Resize(ex_size);
  for(unsigned int e = 0; e < ex_size; ++e){
    forwards[e].Resize(timesteps);
    adjoints[e].Resize(timesteps);
    for(unsigned int t = 0; t < timesteps; ++t){
      forwards[e][t] = &forward.Get(e, NULL, t)->gridelem[App::MECH];
      adjoints[e][t] = &adjoint.Get(e, f, t)->gridelem[App::MECH];
    }
  }
  if(transient){
    forwarddt.Resize(ex_size);
    forwarddtt.Resize(ex_size);
    for(unsigned int e = 0; e < ex_size; ++e){
      forwarddt[e].Resize(timesteps);
      forwarddtt[e].Resize(timesteps);
      for(unsigned int t = 0; t < timesteps; ++t){
        forwarddt[e][t] = &forward.Get(e, NULL, t, FIRST_DERIV)->gridelem[App::MECH];
        forwarddtt[e][t] = &forward.Get(e, NULL, t, SECOND_DERIV)->gridelem[App::MECH];
      }
    }
  }
  

  StdVector<BiLinFormContext*>& biLinForms = assemble_->GetBiLinForms();
  for(unsigned int i = 0; i < biLinForms.GetSize(); ++i){ // loop over all linElastInt bilinear forms (as assemble does)
    BiLinFormContext* biLinForm = biLinForms[i];
    if(biLinForm->GetFirstPde()->GetName() != pde->GetName()) continue;
    if(biLinForm->GetSecondPde()->GetName() != pde->GetName()) continue;
    if(biLinForm->GetIntegrator()->GetName() != "LinElastInt") continue;
    linElastInt* form = (linElastInt*)(biLinForm->GetIntegrator());
    MassInt* massInt = NULL;
    EntityIterator it = biLinForm->GetFirstEntities()->GetIterator();
    for(it.Begin(); !it.IsEnd(); it++){ // loop over all corresponding elements
      const Elem* elem = it.GetElem();
      if(shapedesign->IsElemDependentAtAll(elem->connect)){ // if this element does not depend on any parameters, we can simply skip all the calculations
        int e = elem->elemNum - 1; // index for u and z which are 0-based
        BaseFE* ptelem = elem->ptElem;
        grd->GetElemNodesCoord(CornerCoords, elem->connect, true);
        form->ExtractElemInfo(it); // this is needed before calcBMat

        form->calcDMat(D, elem, DesignElement::NO_DERIVATIVE);
        LOG_DBG2(ShOpt) << "calcDMat returned D=" << D;
        double dampingAlpha = 0.0; double dampingBeta = 0.0;
        if(transient){
          BiLinFormContext* massIntCtxt = GetFormContext(elem->regionId, pde, pde, "MassInt");
          massInt = (MassInt*)massIntCtxt->GetIntegrator();
          if(!design->GetErsatzMaterialDamping(dampingAlpha, dampingBeta, elem)){ // check whether damping is also design and if get it from there
            if(biLinForm->GetSecDestMat() != NOTYPE){
              parser->SetExpr(mathParserHandle, biLinForm->GetSecMatFac());
              dampingBeta = parser->Eval(mathParserHandle);
            }
            if(massIntCtxt->GetSecDestMat() != NOTYPE){
              parser->SetExpr(mathParserHandle, massIntCtxt->GetSecMatFac());
              dampingAlpha = parser->Eval(mathParserHandle);
            }
          }
        }
        LOG_DBG2(ShOpt) << "damping: alpha=" << dampingAlpha << ", beta= " << dampingBeta;
        const bool damping = dampingAlpha > 0.0 || dampingBeta > 0.0;
        
        const unsigned int dimD(D.GetNumRows());

        const int nip = ptelem->GetNumIntPoints();
        const Vector<Double> & intWeights = ptelem->GetIntWeights();
        for(int ip = 1; ip <= nip; ip++){ // loop over all integration points

          double intWeight(intWeights[ip-1]), jacdet;

          ptelem->CalcJacobianAtIp(J, ip, CornerCoords, elem);
          const unsigned int dimJ(J.GetNumCols());

          J.Invert(iJ);
          ptelem->GetGlobDerivShFncAtIp(dPhi, ip, CornerCoords, jacdet, elem); // really is already dPhi * J~
          const unsigned int rowPhi(dPhi.GetNumRows());
          
          if(transient){
            ptelem->GetShFncAtIp(shapeFncAtIp, ip, elem );
            partM.DyadicMult(shapeFncAtIp);
            partM *= massInt->GetErsatzMaterialMass(elem, DesignElement::NO_DERIVATIVE);
            massInt->MassMultiDof(M, partM); // this is M without jacdet & intWeights
          }

          form->CalcBMatOnly(B, ip, ptelem, CornerCoords);

          const unsigned int colB(B.GetNumCols()); 

          BD.Resize(colB, dimD);
          B.MultT(D, BD); // BD = B^T D, this is needed later

          BDB.Resize(colB, colB);
          BD.Mult(B, BDB); // BDB = B^T D B

          for(int p = 0; p < np; p++) { // loop over all parameters
            if(shapedesign->GetElemNodesCoordDerivative(dCornerCoords, elem->connect, p)){ // if dCornerCoords is zero for all nodes, this returns false

              ptelem->CalcJacobianAtIp(dJ, ip, dCornerCoords, elem);

              A1.Resize(dimJ, dimJ);
              dJ.Mult(iJ, A1); // A1 = dJ J~ 

              double trA1 = A1.Trace();

              A2.Resize(rowPhi, dimJ);
              dPhi.Mult(A1, A2); // A2 = (dPhi J~) dJ J~
              form->ReorderBLikeMatrix(A2, A3, ip, ptelem, CornerCoords); // A3 = reorder(dPhi J~ dJ J~) 

              A4.Resize(colB, colB);
              BD.Mult(A3, A4); // A4 = B'D reorder(dPhi J~ dJ J~)

              const unsigned int r = A4.GetNumRows();
              assert(r == A4.GetNumCols());
              // dK = ( -tr(J~ dJ) B'DB + reorder(dPhi J~ dJ J~)'DB + B'D reorder(dPhi J~ dJ J~) ), we compute here with wrong sign
              dK.Resize(r, r);
              for(unsigned int i = 0; i < r; ++i){
                for(unsigned int j = 0; j < r; ++j){
                  dK[i][j] = trA1 * BDB[i][j] - A4[i][j] - A4[j][i];
                }
              }
              LOG_DBG2(ShOpt) << "-tr(iJ dJ) B'DB + reorder(dPhi iJ dJ iJ)'DB + B'D reorder(dPhi iJ dJ iJ)";
              
              if(!homogenization){
                for(unsigned int ex = 0; ex < ex_size; ++ex){
                  StdVector<StdVector<SingleVector*>* >& forward_ex = forwards[ex];
                  StdVector<StdVector<SingleVector*>* >& adjoint_ex = adjoints[ex];
                  StdVector<StdVector<SingleVector*>* >* forwarddt_ex = NULL;
                  StdVector<StdVector<SingleVector*>* >* forwarddtt_ex = NULL;
                  if(transient){
                    forwarddt_ex = &forwarddt[ex];
                    forwarddtt_ex = &forwarddtt[ex];
                  }
                  double v = 0.0;
                  for(unsigned int t = 0; t < timesteps; ++t){ // loop over all timesteps, static analysis has 1 timestep
                    Vector<double>& p_vec = dynamic_cast<Vector<double>&>(*(*adjoint_ex[t])[e]);
                    Vector<double>& u_vec = dynamic_cast<Vector<double>&>(*(*forward_ex[t])[e]);
                    dKp = dK * p_vec;
                    v -= u_vec * dKp;
                    if(transient && (t > 0 || !IsFirstTransientStepStatic())){
                      Vector<double>& utt_vec = dynamic_cast<Vector<double>&>(*(*(*forwarddtt_ex)[t])[e]);
                      dMp = M * p_vec; // note that dM = rho * trA1 * jacdet * NN', M here is rho * NN' (no jacdet), and no trA1
                      v -= utt_vec * dMp * trA1; // keep in mind, we calculate - u1 dA u2
                      if(damping){
                        Vector<double>& ut_vec = dynamic_cast<Vector<double>&>(*(*(*forwarddt_ex)[t])[e]);
                        dDp = dampingAlpha * dMp + dampingBeta * dKp;
                        v -= ut_vec * dDp;
                      }
                    } // if transient
                  } // loop over timesteps
                  der[p] += intWeight * jacdet * me->excitations[ex].weight * v;
                  LOG_DBG2(ShOpt) << "der[" << p << "] += " << intWeight << " * " << jacdet << " * " << me->excitations[ex].weight << " * " << v << " = " << intWeight * jacdet * me->excitations[ex].weight * v << " -> " << der[p];
                } // loop over excitations
              }else{ // we calculate homogenization 
                for(unsigned int ij = 0; ij < ex_size; ++ij){
                  u1diff = *dynamic_cast<Vector<double>* >((*forwards[ij][0])[e]); // assign is needed here
                  SubtractTestDisplacement(ij, CornerCoords, u1diff, tmp_strain, tmp_displacement);
                  for(unsigned int kl = ij; kl < ex_size; ++kl){
                    u2diff = *dynamic_cast<Vector<double>* >((*adjoints[kl][0])[e]);
                    SubtractTestDisplacement(kl, CornerCoords, u2diff, tmp_strain, tmp_displacement);
                    // description see above, is needed twice for SPEED
                    double v1 = 0.0;
                    for(unsigned int i = 0; i < r; ++i){
                      double v2 = 0.0;
                      for(unsigned int j = 0; j < r; ++j){
                        v2 += (A4[i][j] + A4[j][i] -trA1 * BDB[i][j]) * u2diff[j];
                      }
                      v1 += u1diff[i] * v2;
                    }
                    der[p] -= intWeight * jacdet * v1 * (*tensor_diff)[ij][kl] * (ij == kl ? 1.0 : 2.0) * rcubevol; // diagonal is doubled, note homogenization needs PlusUdKu 
                  }
                }                
              }

            } // if GetElemNodesCoordDerivative

          } // parameter loop
        } // integration point loop
      } // if element is dependent on parameters at all
    } // element loop
  } // biLinForm loop
  
  shapedesign->AddAuxDerivatives(f, constraint, der, 1.0);
  
  parser->ReleaseHandle(mathParserHandle);
  */
}

void ShapeOpt::CalcUdF(StateContainer& adjoint, Objective* f, Condition* constraint, double w)
{
  /*
  int np(shapedesign->GetNumberOfAuxParameters());
  const unsigned int ex_size(me->excitations.GetSize());
  UInt timesteps(domain->GetDriver()->GetNumSteps());

  StdVector<double> der; // solution
  der.Resize(np, 0);

  MathParser* parser = domain->GetMathParser();
  double dt = 0.0;
  if(IsTransient()){
    dt = dynamic_cast<TransientDriver*>(domain->GetDriver())->GetDeltaT();
  }
  parser->SetValue(MathParser::GLOB_HANDLER, "dt", dt);

  Grid* grd = domain->GetGrid();
  
  Matrix<double> CornerCoords;
  Vector<double> normal;
  Vector<double> elemVec;
  Vector<double> delemVec;
  Matrix<double> J;
  Vector<double> shapeFnc;
  Matrix<double> dCornerCoords;
  Matrix<double> dJ;
  Vector<double> dnormal;
  Vector<Double> vec1(3), vec2(3), dvec1(3), dvec2(3);
  
  unsigned int dim = grd->GetDim();

  for(unsigned int ex = 0; ex < ex_size; ++ex){
    for(unsigned int t = 0; t < timesteps; ++t){
      parser->SetValue(MathParser::GLOB_HANDLER, "t", dt*(t+1)); // GetPressureFactor uses this
      parser->SetValue(MathParser::GLOB_HANDLER, "step", t+1);
      StateSolution* u = adjoint.Get(ex, f, t);
      StdVector<SingleVector*>* u_vec = &u->gridelem[App::MECH];
      
      Excitation& excite = me->excitations[ex];
      
      StdVector<LinearFormContext*>& linForms = excite.GetLinForms();
      for(unsigned int i = 0; i < linForms.GetSize(); ++i){ // loop over all pressure linear forms (as assemble does)
        LinearFormContext* linForm = linForms[i];
        if(linForm->GetPde()->GetName() != pde->GetName()) continue;
        if(linForm->GetIntegrator()->GetName() != "PressureLinForm") continue;
        PressureLinForm* form = (PressureLinForm*)(linForm->GetIntegrator());
        EntityIterator it = linForm->GetEntities()->GetIterator();
        for(it.Begin(); !it.IsEnd(); it++){ // loop over all corresponding elements
          const SurfElem* elem = it.GetSurfElem();
          int e = elem->elemNum - 1; // index for u and z which are 0-based
          BaseFE* ptelem = elem->ptElem;
          grd->GetElemNodesCoord(CornerCoords, elem->connect, true);

          Vector<double>& uelem = dynamic_cast<Vector<double>& >(*(*u_vec)[e]);
          double pres = form->GetPressureFactor(it);
          grd->CalcSurfNormal(normal, *elem, true);

          elemVec.Resize(uelem.GetSize());
          delemVec.Resize(uelem.GetSize());

          const int nip = ptelem->GetNumIntPoints();
          const Vector<Double> & intWeights = ptelem->GetIntWeights();
          for(int ip = 1; ip <= nip; ip++){ // loop over all integration points
            double intWeight = intWeights[ip-1];

            ptelem->CalcJacobianAtIp(J, ip, CornerCoords, elem);
            // note that the "Jacobian" is a 2x1 matrix, for calculation of |J| see linefe.cc LineFE::CalcJacobianAtIp
            double detJ = sqrt(J[0][0]*J[0][0] + J[1][0]*J[1][0]);
            ptelem->GetShFncAtIp(shapeFnc, ip, elem);
            for (UInt i=0; i<dim; i++) { // see linPressureInt
              double n = -normal[i];
              for (UInt j=0; j<shapeFnc.GetSize(); j++) {
                elemVec[j*(dim) +i] = shapeFnc[j] * n;
              }
            }
            double f = intWeight * pres * excite.weight;
            double v = uelem * elemVec;

            for(int p = 0; p < np; p++) { // loop over all parameters
              if(shapedesign->GetElemNodesCoordDerivative(dCornerCoords, elem->connect, p)){ // this returns false if that nodes do not depend on p at all
                ptelem->CalcJacobianAtIp(dJ, ip, dCornerCoords, elem);
                double ddetJ = (dJ[0][0]*J[0][0] + dJ[1][0]*J[1][0]) / detJ;

                // calculate derivative of normal ( see grid_cfs.cc CalcSurfNormal )
                if (ptelem->GetDim() == 1) {
                  double dx  = CornerCoords[0][1] - CornerCoords[0][0]; double ddx = dCornerCoords[0][1] - dCornerCoords[0][0];
                  double dy  = CornerCoords[1][1] - CornerCoords[1][0]; double ddy = dCornerCoords[1][1] - dCornerCoords[1][0];
                  double len = sqrt(dx*dx+dy*dy);
                  double dilen = 2.0 * (dx*ddx + dy*ddy) / (len*len*len); // derivative of 1/len
                  dnormal.Resize(2);
                  dnormal[0] = ddy/len + dy*dilen;
                  dnormal[1] = -ddx/len - dx*dilen;
                } else {
                  // note that we have to assume a tetrahedral element
                  // what is the normal of a quadrilateral el.??? as changing one point would result in no quadrilateral el. any more
                  // however if meshdeformations ensure the quadr. el. stay quadr. el., this derivative is correct
                  UInt surfCorners = ptelem->GetNumCorners();

                  //compute the two vectors in the plane
                  for (UInt i=0; i<3; i++) {
                    vec1[i] = CornerCoords[i][1]             - CornerCoords[i][0];
                    vec2[i] = CornerCoords[i][surfCorners-1] - CornerCoords[i][0];
                    dvec1[i] = dCornerCoords[i][1]             - dCornerCoords[i][0];
                    dvec2[i] = dCornerCoords[i][surfCorners-1] - dCornerCoords[i][0];
                  }
                  //compute cross product
                  Vector<Double> n(3);
                  n[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
                  n[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
                  n[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
                  dnormal[0] = dvec1[1]*vec2[2] - dvec1[2]*vec2[1] + vec1[1]*dvec2[2] - vec1[2]*dvec2[1];
                  dnormal[1] = dvec1[2]*vec2[0] - dvec1[0]*vec2[2] + vec1[2]*dvec2[0] - vec1[0]*dvec2[2];
                  dnormal[2] = dvec1[0]*vec2[1] - dvec1[1]*vec2[0] + vec1[0]*dvec2[1] - vec1[1]*dvec2[0];
                  //normalize the length to 1
                  Double length = n.NormL2();
                  Double rlength = 1/length;
                  Double sumEntries = 0;
                  for(UInt i=0; i<3; i++){
                    sumEntries += n[i] + dnormal[i];
                  }
                  n *= -sumEntries*rlength*rlength;
                  dnormal += n;
                  dnormal *= rlength;
                }
                for (UInt i=0; i<dim; i++) { // see linPressureInt
                  double n = -dnormal[i];
                  for (UInt j=0; j<shapeFnc.GetSize(); j++) {
                    delemVec[j*(dim) +i] = shapeFnc[j] * n;
                  }
                }

                const double vadd = f * (ddetJ * v + detJ * (uelem * delemVec)); 
                der[p] += vadd;
              }
            } // parameter loop
          } // integration point loop
        } // element loop
      } // linForm loop
    } // timesteps
  } // excitations  

  shapedesign->AddAuxDerivatives(f, constraint, der, w);
  */
}

double ShapeOpt::CalcCompliance(Excitation& excite, Objective* f, Condition* constraint, bool derivative)
{
  if(derivative){
    // the derivative of tracking w.r.t. shape is: - u' dA/dShape u + 2 u dF/dShape
    if(excite.index == (int) me->excitations.GetSize() - 1){
      CalcMinusU1dKU2(forward, forward, f, constraint); // todo: Is this correct for transient?
      CalcUdF(IsTransient() ? adjoint : forward, f, constraint, 2.0); // in Transient case the system is not self-adjoint any more
    }
    // this however is done per excite
    if(exoprt_fe_design_){
      ErsatzMaterial::CalcCompliance(excite, f, constraint, true);
    }
  }else{
    return(ErsatzMaterial::CalcCompliance(excite, f, constraint, derivative));
  }
  return 0.0;
}

double ShapeOpt::CalcTracking(Excitation& excite, Objective* f, Condition* constraint, bool derivative){
  if(derivative){
    // the derivative of tracking w.r.t. shape is: - z' dA/dShape u + z dF/dShape
    if(excite.index == (int) me->excitations.GetSize() - 1){ // CalcMinusU1dKU2 runs over the excitations, this speeds up things a little
      CalcMinusU1dKU2(forward, adjoint, f, constraint);
      CalcUdF(adjoint, f, constraint);
    }
    // this however is done per excite
    if(exoprt_fe_design_){
      ErsatzMaterial::CalcTracking(excite, f, constraint, true);
    }
  }else{
    return(ErsatzMaterial::CalcTracking(excite, f, constraint, derivative));
  }
  return 0.0;
}

void ShapeOpt::CalcHomogenizedTrackingGradient(const Matrix<double>& target, const Matrix<double>& hom, Function* f){
  Matrix<double> tensor_diff;
  tensor_diff = hom - target;
  CalcMinusU1dKU2(forward, forward, dynamic_cast<Objective*>(f), dynamic_cast<Condition*>(f), &tensor_diff);
}

Matrix<double> ShapeOpt::CalcHomogenizedTensor(Function* f)
{
  const unsigned int ex_size(me->excitations.GetSize());
  assert((dim == 2 && ex_size == 3) || (dim == 3 && ex_size == 6));
  
  
  Matrix<double> result(ex_size, ex_size);
  result.Init();
/*
  assert(false);
  double rcubevol(1.0 / grid->CalcVolumeSpannedByNamedNodes());

  Matrix<double> elemMat;
  Matrix<double> tmp_strain(dim, dim);
  Matrix<double> tmp_displacement;
  Matrix<double> CornerCoords;
  Vector<double> u1diff;
  Vector<double> u2diff;
  Vector<double> Ku;
  
  StdVector<BiLinFormContext*>& biLinForms = assemble_->GetBiLinForms();
  for(unsigned int i = 0; i < biLinForms.GetSize(); ++i){ // loop over all linElastInt bilinear forms (as assemble does)
    BiLinFormContext* biLinForm = biLinForms[i];
    if(biLinForm->GetFirstPde()->GetName() != pde->GetName()) continue;
    if(biLinForm->GetSecondPde()->GetName() != pde->GetName()) continue;
    if(biLinForm->GetIntegrator()->GetName() != "LinElastInt") continue;
    linElastInt* form = (linElastInt*)(biLinForm->GetIntegrator());
    EntityIterator it = biLinForm->GetFirstEntities()->GetIterator();
    for(it.Begin(); !it.IsEnd(); it++){ // loop over all corresponding elements
      const Elem* elem = it.GetElem();
      int e = elem->elemNum - 1;
      grid->GetElemNodesCoord(CornerCoords, elem->connect, true);
      form->CalcElementMatrix(elemMat, it, it);
      for(unsigned int ij = 0; ij < ex_size; ++ij){
        u1diff = *dynamic_cast<Vector<double>* >(forward.Get(ij)->gridelem[App::MECH][e]);
        SubtractTestDisplacement(ij, CornerCoords, u1diff, tmp_strain, tmp_displacement);
        Ku = elemMat * u1diff;
        for(unsigned int kl = ij; kl < ex_size; ++kl){ // only upper triangle
          u2diff = *dynamic_cast<Vector<double>* >(forward.Get(kl)->gridelem[App::MECH][e]);
          SubtractTestDisplacement(kl, CornerCoords, u2diff, tmp_strain, tmp_displacement);
          result[ij][kl] += Ku * u2diff;
        }
      }
    } // elem loop
  }
  // copy the rest of the tensor and divide the upper triangle (and diag)
  for(unsigned int ij = 0; ij < ex_size; ++ij){
    for(unsigned int kl = 0; kl < ex_size; ++kl){
      if(ij <= kl){
        result[ij][kl] *= rcubevol;
      }else{
        result[ij][kl] = result[kl][ij];
      }
    }
  }

  homogenizedTensor.Assign(result, 1.0);
*/
  return result;
}

void ShapeOpt::StorePDESolution(StateContainer& solutions, Excitation &excite, Function* f, UInt timestep, bool read_sol, bool read_rhs, bool save_sol, TimeDeriv derivative, const std::string& comment){
  ParamMat::StorePDESolution(solutions, excite, f, timestep, read_sol, read_rhs, save_sol, derivative, comment);
  if(read_sol){
    solutions.Get(excite, f, timestep, derivative)->Read(StateSolution::GRIDELEM_VECTORS, Optimization::context->pde, App::MECH, false, derivative);
  }
}

void ShapeOpt::SubtractTestDisplacement(unsigned int idx, Matrix<double>& CornerCoords, Vector<double>& result, Matrix<double>& tmp_strain, Matrix<double>& tmp_displacement){
  SetTestStrainMatrix(App::MECH,tmp_strain, me->excitations[idx].test_strain);
  unsigned int cols = CornerCoords.GetNumCols();
  tmp_displacement.Resize(dim, cols);
  tmp_strain.Mult(CornerCoords, tmp_displacement);
  for(unsigned int out = 0; out < cols; ++out){
    for(unsigned int in = 0; in < dim; ++in){
      result[out*dim + in] -= tmp_displacement[in][out];
    }
  }
}

} // namespace
