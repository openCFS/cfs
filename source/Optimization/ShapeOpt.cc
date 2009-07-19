#include "Optimization/ShapeOpt.hh"
#include "Optimization/ShapeDesign.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "Forms/linElastInt.hh"
#include "Forms/linPressureInt.hh"
#include "PDE/mechPDE.hh"
#include "Optimization/DesignElement.hh"
#include "Driver/assemble.hh"

using namespace CoupledField;

ShapeOpt::ShapeOpt() : ParamMat() {
  shapedesign = dynamic_cast<ShapeDesign*>(design);

  ParamNode* sopn = pn->Get("shapeOpt");
  shapedesign->Configure(sopn, constraints.GetSize());
  alsomatopt_ = shapedesign->AlsoMatOpt();

  // all (bi)linear forms need to use updated coordinates
  std::set<BiLinFormContext*>* biLinForms = assemble_->GetBiLinForms();
  for(std::set<BiLinFormContext*>::iterator iBiLinForm = biLinForms->begin(); iBiLinForm != biLinForms->end(); iBiLinForm++){
    (*iBiLinForm)->GetIntegrator()->SetUseCoordUpdate(true);
  }
  // set the linearForms used in multiple excitations, note that this does contain all linearforms (some even several times)
  for(unsigned int i = 0; i < excitations.GetSize(); i++){
    std::set<LinearFormContext*>* linForms = excitations[i].GetLinForms();
    for(std::set<LinearFormContext*>::iterator iForm = linForms->begin(); iForm != linForms->end(); iForm++){
      (*iForm)->GetIntegrator()->SetUseCoordUpdate(true);
    }
  }
}

double ShapeOpt::CalcVolume(bool derivative, Condition* constraint, bool normalized){
  if(derivative){
    StdVector<double> der; // solution
    int np = shapedesign->GetNumberOfShapeParameters();
    der.Resize(np, 0);
    if(alsomatopt_){
      // this needs to be done before, we do use fraction
      ErsatzMaterial::CalcVolume(derivative, constraint, normalized);
    }
    if(!alsomatopt_ || (constraint && constraint->design == DesignElement::UNITY)){
      if(!normalized){
        // this is independent of material optimization, simply the derivative of the real volume
        Grid* grd = domain->GetGrid();
        StdVector<RegionIdType> regs;
        grd->GetVolRegionIds(regs);
        for(unsigned int ri = 0; ri < regs.GetSize(); ri++){
          RegionIdType rid = regs[ri];
          if(!constraint || constraint->IsForRegion(rid)){
            StdVector<Elem*> elems;
            grd->GetElems(elems,rid);
            for( UInt i = 0; i < elems.GetSize(); i++ ) {
              const Elem* elem = elems[i];
              Matrix<double> CornerCoords;
              BaseFE* ptelem = elem->ptElem;
              grd->GetElemNodesCoord(CornerCoords, elem->connect, true );
              const int nip = ptelem->GetNumIntPoints();
              const Vector<Double> & intWeights = ptelem->GetIntWeights();

              for(int ip = 1; ip <= nip; ip++){ // loop over all integration points
                Matrix<double> J;
                ptelem->CalcJacobianAtIp(J, ip, CornerCoords, elem);
                Matrix<double> iJ;
                J.Invert(iJ);
                double det;
                J.Determinant(det);
                double w = det * intWeights[ip-1];
                for(int p = 0; p < np; p++) { // loop over all parameters
                  Matrix<double> dCornerCoords;
                  shapedesign->GetElemNodesCoordDerivative(dCornerCoords, elem->connect, p);
                  Matrix<double> dJ;
                  ptelem->CalcJacobianAtIp(dJ, ip, dCornerCoords, elem);
                  Matrix<double> diJ;
                  diJ = iJ * dJ;
                  double tr;
                  diJ.Trace(tr); // tr = trace(iJ*dJ) = trace(dJ*iJ)
                  der[p] += w * tr;
                } // params
              } // int.points
            } // elems
          } // if region
        } // region
      } // normalized
    }else{ // volume is average or sum of design
      // this is similar to ErsatzMaterial::CalcVolume but calculates derivatives w.r.t. shape
      Grid* grd = domain->GetGrid();
      bool isObjective = constraint == NULL;
      double fraction = isObjective ? volume_fraction_ : constraint->volume_fraction_; // this already considers everything
      double volume = 0.0;
      if(!normalized){  // needed for derivative in normalized versions
        volume = isObjective ? cost->GetValue() : constraint->value;
        if(volume == 0.0){ // if function was never evaluated before
          CalcVolume(false, constraint, normalized);
          volume = isObjective ? cost->GetValue() : constraint->value;
        }
      }
      bool allDesignsRelevant = constraint == NULL || constraint->design == DesignElement::TENSOR_TRACE || constraint->design == DesignElement::DEFAULT;
      bool ersatzMaterialTensor = domain->HasErsatzMaterialTensor() && allDesignsRelevant;
      unsigned int upper = ersatzMaterialTensor ? design->GetNumberOfElements() : design->data.GetSize();
      for(unsigned int i = 0; i < upper; i++) {
        DesignElement* de = &design->data[i];
        bool relevant = (allDesignsRelevant || constraint->design == de->GetType()) && (isObjective || constraint->IsForRegion(de->elem->regionId));
        if(relevant){
          double des;
          if(ersatzMaterialTensor){ // use the trace of the stiffness Tensor as "volume"
            Matrix<double> material;
            GetErsatzMaterialTensor(material, de->elem);
            material.Trace(des);
          }else{
            des = de->GetDesign(DesignElement::PLAIN);
          }
          const Elem* elem = de->elem;
          Matrix<double> CornerCoords;
          BaseFE* ptelem = elem->ptElem;
          grd->GetElemNodesCoord(CornerCoords, elem->connect, true );
          const int nip = ptelem->GetNumIntPoints();
          const Vector<Double> & intWeights = ptelem->GetIntWeights();

          for(int ip = 1; ip <= nip; ip++){ // loop over all integration points
            double intWeight = intWeights[ip-1];
            Matrix<double> J;
            ptelem->CalcJacobianAtIp(J, ip, CornerCoords, elem);
            Matrix<double> iJ;
            J.Invert(iJ);
            for(int p = 0; p < np; p++) { // loop over all parameters
              Matrix<double> dCornerCoords;
              shapedesign->GetElemNodesCoordDerivative(dCornerCoords, elem->connect, p);
              Matrix<double> dJ;
              ptelem->CalcJacobianAtIp(dJ, ip, dCornerCoords, elem);
              Matrix<double> diJ;
              diJ = iJ * dJ;
              double tr, det;
              diJ.Trace(tr);
              J.Determinant(det);
              der[p] += fraction * intWeight * tr * det * (des - volume); // fraction * intWeight * dArea/dalpha * (d - volume)    ( fraction = 1 und volume = 0 in non-normalized version)
            } // params
          } // int.points
        }
      }
    }
    // derivative in direction of our parameters is always needed and calculated here
    // derivative in direction of design-element parameters only if on that regions
    // these derivatives are independent of our parameters and can be calculated as before
    shapedesign->SetShapeDerivatives(constraint, der);
  }else{ // derivative
    if(!alsomatopt_ || (constraint && constraint->design == DesignElement::UNITY)){ // this is the real volume, not multiplied by design, we also use this if no design available
      // if design is unity, we use the grid instead of designspace
      if(normalized) return(1.0);
      Grid* grd = domain->GetGrid();
      StdVector<RegionIdType> regs;
      grd->GetVolRegionIds(regs);
      double s = 0.0;
      for(unsigned int i = 0; i < regs.GetSize(); i++){
        int rid = regs[i];
        if(!constraint || constraint->IsForRegion(rid)){
          s += grd->CalcVolumeOfRegion(rid, false, true);
        }
      }
      return(s);
    }else{ // working on a design, alsomatopt_ must be true
      return(ErsatzMaterial::CalcVolume(derivative, constraint, normalized));
    }
  }
  return 0.0;
}

void ShapeOpt::CalcMinusU1dKU2(StdVector<SingleVector*>& u1, StdVector<SingleVector*>& u2, Condition* constraint, double w){
  StdVector<double> der; // solution
  int np = shapedesign->GetNumberOfShapeParameters();
  der.Resize(np, 0);

  Grid* grd = domain->GetGrid();

  std::set<BiLinFormContext*>* biLinForms = assemble_->GetBiLinForms();
  for(std::set<BiLinFormContext*>::iterator iBiLinForm = biLinForms->begin(); iBiLinForm != biLinForms->end(); iBiLinForm++){ // loop over all linElastInt bilinear forms (as assemble does)
    BiLinFormContext* biLinForm = *iBiLinForm;
    if(biLinForm->GetFirstPde()->GetName() != mech->GetName()) continue;
    if(biLinForm->GetSecondPde()->GetName() != mech->GetName()) continue;
    if(biLinForm->GetIntegrator()->GetName() != "linElastInt") continue;
    linElastInt* form = (linElastInt*)(biLinForm->GetIntegrator());
    EntityIterator it = biLinForm->GetFirstEntities()->GetIterator();
    for(it.Begin(); !it.IsEnd(); it++){ // loop over all corresponding elements
      const Elem* elem = it.GetElem();
      int e = elem->elemNum - 1; // index for u and z which are 0-based
      Matrix<double> CornerCoords;
      BaseFE* ptelem = elem->ptElem;
      grd->GetElemNodesCoord(CornerCoords, elem->connect, true);

      form->ExtractElemInfo(it); // this is needed before calcBMat
      
      Matrix<double> D;
      form->calcDMat(D, elem, DesignElement::NO_DERIVATIVE);

      Vector<double>& u1elem = dynamic_cast<Vector<double>& >(*u1[e]);
      Vector<double>& u2elem = dynamic_cast<Vector<double>& >(*u2[e]);
      
      const int nip = ptelem->GetNumIntPoints();
      const Vector<Double> & intWeights = ptelem->GetIntWeights();
      for(int ip = 1; ip <= nip; ip++){ // loop over all integration points

        double intWeight = intWeights[ip-1], jacdet;

        Matrix<double> J;
        ptelem->CalcJacobianAtIp(J, ip, CornerCoords, elem);
        Matrix<double> iJ;
        J.Invert(iJ);
        Matrix<double> dPhi;
        ptelem->GetGlobDerivShFncAtIp(dPhi, ip, CornerCoords, jacdet, elem); // really is already dPhi * J~

        Matrix<double> B;
        form->calcBMatOnly(B, ip, ptelem, CornerCoords);
        Matrix<double> BT;
        B.Transpose(BT);

        Matrix<double> BD;
        BD = BT * D; // B^T D

        Matrix<double> BDB;
        BDB = BD * B; // B^T D B
        
        for(int p = 0; p < np; p++) { // loop over all parameters
          Matrix<double> dCornerCoords;
          shapedesign->GetElemNodesCoordDerivative(dCornerCoords, elem->connect, p);

          Matrix<double> dJ;
          ptelem->CalcJacobianAtIp(dJ, ip, dCornerCoords, elem);
          
          Matrix<double> A1;
          A1 = dJ * iJ;

          double trA1;
          A1.Trace(trA1);

          Matrix<double> A2;
          A2 = dPhi * A1;
          form->ReorderBLikeMatrix(A2, A1, ip, ptelem, CornerCoords);

          A2 = BD * A1;
          A2.Transpose(A1);
          A2.Add(1.0, A1);
          A2.Add(-trA1, BDB);

          // A2 = ( -tr(J~ dJ) B'DB + reorder(dPhi J~ dJ J~)'DB + B'D reorder(dPhi J~ dJ J~) )
          Vector<double> A2u;
          A2u.Resize(u2elem.GetSize());
          A2.Mult(u2elem, A2u);
          const double vadd = intWeight * jacdet * (u1elem * A2u); 
          der[p] += vadd;

        } // parameter loop
      } // integration point loop
    } // element loop
  } // biLinForm loop
  
  shapedesign->AddShapeDerivatives(constraint, der, w);
}

void ShapeOpt::CalcUdF(Excitation& excite, StdVector<SingleVector*>& u, Condition* constraint, double w){
  StdVector<double> der; // solution
  int np = shapedesign->GetNumberOfShapeParameters();
  der.Resize(np, 0);

  Grid* grd = domain->GetGrid();
  
  unsigned int dim = grd->GetDim();
  std::set<LinearFormContext*>* linForms = excite.GetLinForms();
  for(std::set<LinearFormContext*>::iterator iLinForm = linForms->begin(); iLinForm != linForms->end(); iLinForm++){ // loop over all pressure linear forms (as assemble does)
    LinearFormContext* linForm = *iLinForm;
    if(linForm->GetPde()->GetName() != mech->GetName()) continue;
    if(linForm->GetIntegrator()->GetName() != "PressureLinForm") continue;
    PressureLinForm* form = (PressureLinForm*)(linForm->GetIntegrator());
    EntityIterator it = linForm->GetEntities()->GetIterator();
    for(it.Begin(); !it.IsEnd(); it++){ // loop over all corresponding elements
      const SurfElem* elem = it.GetSurfElem();
      int e = elem->elemNum - 1; // index for u and z which are 0-based
      Matrix<double> CornerCoords;
      BaseFE* ptelem = elem->ptElem;
      grd->GetElemNodesCoord(CornerCoords, elem->connect, true);

      Vector<double>& uelem = dynamic_cast<Vector<double>& >(*u[e]);
      double pres = form->GetPressureFactor(elem);
      Vector<double> normal;
      grd->CalcSurfNormal(normal, *elem, true);

      Vector<double> elemVec;
      elemVec.Resize(uelem.GetSize());
      Vector<double> delemVec;
      delemVec.Resize(uelem.GetSize());

      const int nip = ptelem->GetNumIntPoints();
      const Vector<Double> & intWeights = ptelem->GetIntWeights();
      for(int ip = 1; ip <= nip; ip++){ // loop over all integration points
        double intWeight = intWeights[ip-1];

        Matrix<double> J;
        ptelem->CalcJacobianAtIp(J, ip, CornerCoords, elem);
        // note that the "Jacobian" is a 2x1 matrix, for calculation of |J| see linefe.cc LineFE::CalcJacobianAtIp
        double detJ = sqrt(J[0][0]*J[0][0] + J[1][0]*J[1][0]);
        Vector<double> shapeFnc;
        ptelem->GetShFncAtIp(shapeFnc, ip, elem);
        for (UInt i=0; i<dim; i++) { // see linPressureInt
          double n = -normal[i];
          for (UInt j=0; j<shapeFnc.GetSize(); j++) {
            elemVec[j*(dim) +i] = shapeFnc[j] * n;
          }
        }
        double f = intWeight*pres;
        double v = uelem * elemVec;

        for(int p = 0; p < np; p++) { // loop over all parameters
          Matrix<double> dCornerCoords;
          shapedesign->GetElemNodesCoordDerivative(dCornerCoords, elem->connect, p);
          Matrix<double> dJ;
          ptelem->CalcJacobianAtIp(dJ, ip, dCornerCoords, elem);
          double ddetJ = (dJ[0][0]*J[0][0] + dJ[1][0]*J[1][0]) / detJ;

          // calculate derivative of normal ( see grid_cfs.cc CalcSurfNormal )
          Vector<double> dnormal;
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
            Vector<Double> vec1(3), vec2(3), dvec1(3), dvec2(3);
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

        } // parameter loop
      } // integration point loop
    } // element loop
  } // linForm loop

  shapedesign->AddShapeDerivatives(constraint, der, w);
}

double ShapeOpt::CalcCompliance(Excitation& excite, bool derivative, Condition* constraint){
  if(derivative){
    forward.data[excite.index]->Read(Solution::GRIDELEM_VECTORS, mech, MECH);
    StdVector<SingleVector*>& u = forward.data[excite.index]->gridelem[MECH];
    
    // the derivative of tracking w.r.t. shape is: - u' dA/dShape u + 2 u dF/dShape 
    CalcMinusU1dKU2(u, u, constraint, excite.weight);
    CalcUdF(excite, u, constraint, 2*excite.weight);
    if(alsomatopt_){
      ErsatzMaterial::CalcCompliance(excite, true, constraint);
    }
  }else{
    return(ErsatzMaterial::CalcCompliance(excite, derivative, constraint));
  }
  return 0.0;
}

double ShapeOpt::CalcTracking(Excitation& excite, bool derivative, Condition* constraint, bool solveproblem){
  if(derivative){
    if(solveproblem){
      SolveTrackingProblem(excite, alsomatopt_, true);
    }

    StdVector<SingleVector*>& z = tracking_->gridelem[MECH];
    forward.data[excite.index]->Read(Solution::GRIDELEM_VECTORS, mech, MECH);
    StdVector<SingleVector*>& u = forward.data[excite.index]->gridelem[MECH];
    
    // the derivative of tracking w.r.t. shape is: - z' dA/dShape u + z dF/dShape 
    CalcMinusU1dKU2(z, u, constraint, excite.weight);
    CalcUdF(excite, u, constraint, excite.weight);
    
    if(alsomatopt_){
      ErsatzMaterial::CalcTracking(excite, true, constraint, false);
    }
  }else{
    return(ErsatzMaterial::CalcTracking(excite, derivative, constraint));
  }
  return 0.0;
}

