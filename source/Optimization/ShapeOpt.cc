#include "Optimization/ShapeOpt.hh"
#include "Optimization/Design/ShapeDesign.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "Forms/linElastInt.hh"
#include "Forms/linPressureInt.hh"
#include "PDE/mechPDE.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Driver/assemble.hh"

using namespace CoupledField;

ShapeOpt::ShapeOpt() : ParamMat() {
  shapedesign = dynamic_cast<ShapeDesign*>(design);

  PtrParamNode sopn = pn->Get("shapeOpt");
  shapedesign->Configure(sopn, objectives.data.GetSize(), constraints.view->GetNumberOfActiveConstraints());
  alsomatopt_ = shapedesign->AlsoMatOpt();

  // all (bi)linear forms need to use updated coordinates
  StdVector<BiLinFormContext*>& biLinForms = assemble_->GetBiLinForms();
  for(unsigned int i = 0; i < biLinForms.GetSize(); ++i){
    biLinForms[i]->GetIntegrator()->SetUseCoordUpdate(true);
  }
  // set the linearForms used in multiple excitations, note that this does contain all linearforms (some even several times)
  for(unsigned int i = 0; i < me->excitations.GetSize(); i++){
    StdVector<LinearFormContext*>& linForms = me->excitations[i].GetLinForms();
    for(unsigned int j = 0; j < linForms.GetSize(); ++j){
      linForms[j]->GetIntegrator()->SetUseCoordUpdate(true);
    }
  }
}

double ShapeOpt::CalcVolume(Objective* f, Condition* constraint, bool derivative, bool normalized){
  // the exponent is used in Ersatzmaterial for the volume cost function
  // if an exponent != 1.0 at this point makes any sense is unknown
  
  if(derivative){
    StdVector<double> der; // solution
    Matrix<double> CornerCoords;
    Matrix<double> J;
    Matrix<double> iJ;
    Matrix<double> dCornerCoords;
    Matrix<double> dJ;
    Matrix<double> diJ;

    int np = shapedesign->GetNumberOfShapeParameters();
    der.Resize(np, 0);
    if(alsomatopt_){
      // this needs to be done before, we do use fraction
      ErsatzMaterial::CalcVolume(f, constraint, derivative, normalized);
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
              if(shapedesign->IsElemDependentAtAll(elem->connect)){ // if this element does not depent on any parameters, we can simply skip all the calculations
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
                      double tr;
                      diJ.Trace(tr); // tr = trace(iJ*dJ) = trace(dJ*iJ)
                      der[p] += w * tr;
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
      bool isObjective = constraint == NULL;
      double fraction = isObjective ? volume_fraction_ : constraint->volume_fraction; // this already considers everything
      double volume = 0.0;
      if(!normalized){  // needed for derivative in normalized versions
        volume = CalcVolume(f, constraint, derivative, normalized);
      }
      bool allDesignsRelevant = constraint == NULL || constraint->design == DesignElement::TENSOR_TRACE || constraint->design == DesignElement::DEFAULT;
      bool ersatzMaterialTensor = domain->HasErsatzMaterialTensor() && allDesignsRelevant;
      unsigned int upper = ersatzMaterialTensor ? design->GetNumberOfElements() : design->data.GetSize();
      Matrix<double> material;
      for(unsigned int i = 0; i < upper; i++) {
        DesignElement* de = &design->data[i];
        bool relevant = (allDesignsRelevant || constraint->design == de->GetType()) && (isObjective || constraint->IsForRegion(de->elem->regionId));
        const Elem* elem = de->elem;
        if(relevant && shapedesign->IsElemDependentAtAll(elem->connect)){
          double des;
          if(ersatzMaterialTensor){ // use the trace of the stiffness Tensor as "volume"
            GetErsatzMaterialTensor(material, de->elem);
            material.Trace(des);
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
                double tr;
                diJ.Trace(tr);
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
    shapedesign->AddShapeDerivatives(f, constraint, der, 1.0);
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
      return(ErsatzMaterial::CalcVolume(f, constraint, derivative, normalized));
    }
  }
  return 0.0;
}

void ShapeOpt::CalcMinusU1dKU2(Solutions& u1, Solutions& u2, Objective* f, Condition* constraint, const Matrix<double>* tensor_diff){
  StdVector<double> der; // solution
  int np = shapedesign->GetNumberOfShapeParameters();
  der.Resize(np, 0);
  const bool homogenization = tensor_diff != NULL;
  const unsigned int ex_size(me->excitations.GetSize());
  double rcubevol(1.0);
  if(homogenization){
    rcubevol = 1.0 / grid->CalcVolumeSpannedByNamedNodes();    
  }

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
  
  Matrix<double> tmp_strain(dim, dim); // homogenization
  Matrix<double> tmp_displacement;
  Vector<double> u1diff;
  Vector<double> u2diff;

  StdVector<BiLinFormContext*>& biLinForms = assemble_->GetBiLinForms();
  for(unsigned int i = 0; i < biLinForms.GetSize(); ++i){ // loop over all linElastInt bilinear forms (as assemble does)
    BiLinFormContext* biLinForm = biLinForms[i];
    if(biLinForm->GetFirstPde()->GetName() != pde->GetName()) continue;
    if(biLinForm->GetSecondPde()->GetName() != pde->GetName()) continue;
    if(biLinForm->GetIntegrator()->GetName() != "linElastInt") continue;
    linElastInt* form = (linElastInt*)(biLinForm->GetIntegrator());
    EntityIterator it = biLinForm->GetFirstEntities()->GetIterator();
    for(it.Begin(); !it.IsEnd(); it++){ // loop over all corresponding elements
      const Elem* elem = it.GetElem();
      if(shapedesign->IsElemDependentAtAll(elem->connect)){ // if this element does not depend on any parameters, we can simply skip all the calculations
        int e = elem->elemNum - 1; // index for u and z which are 0-based
        BaseFE* ptelem = elem->ptElem;
        grd->GetElemNodesCoord(CornerCoords, elem->connect, true);
        form->ExtractElemInfo(it); // this is needed before calcBMat

        form->calcDMat(D, elem, DesignElement::NO_DERIVATIVE);
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

          form->calcBMatOnly(B, ip, ptelem, CornerCoords);

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

              double trA1;
              A1.Trace(trA1);

              A2.Resize(rowPhi, dimJ);
              dPhi.Mult(A1, A2); // A2 = (dPhi J~) dJ J~
              form->ReorderBLikeMatrix(A2, A3, ip, ptelem, CornerCoords); // A3 = reorder(dPhi J~ dJ J~) 

              A4.Resize(colB, colB);
              BD.Mult(A3, A4); // A4 = B'D reorder(dPhi J~ dJ J~)

              const unsigned int r = A4.GetNumRows();
              const unsigned int c = A4.GetNumCols();
              // A5 = ( -tr(J~ dJ) B'DB + reorder(dPhi J~ dJ J~)'DB + B'D reorder(dPhi J~ dJ J~) ), this is not needed completely
              
              if(!homogenization){
                for(unsigned int ex = 0; ex < ex_size; ++ex){
                  // Get() requires f exclusively for adjoint solutions.
                  Vector<double>& u1elem = dynamic_cast<Vector<double>& >(*u1.Get(ex, &u1 == &u2 ? NULL : f)->gridelem[MECH][e]);
                  Vector<double>& u2elem = dynamic_cast<Vector<double>& >(*u2.Get(ex)->gridelem[MECH][e]);
                  double v1 = 0.0; // v1 = u1elem' * (A5 * u2elem)
                  for(unsigned int i = 0; i < r; ++i){
                    double v2 = 0.0; // v2 = component i of (A5 * u2elem)
                    for(unsigned int j = 0; j < c; ++j){
                      v2 += (A4[i][j] + A4[j][i] -trA1 * BDB[i][j]) * u2elem[j];
                    }
                    v1 += u1elem[i] * v2;
                  }
                  der[p] += intWeight * jacdet * v1 * me->excitations[ex].weight;
                }
              }else{ // we calculate homogenization 
                for(unsigned int ij = 0; ij < ex_size; ++ij){
                  u1diff = *dynamic_cast<Vector<double>* >(u1.Get(ij)->gridelem[MECH][e]); // assign is needed here
                  SubtractTestDisplacement(ij, CornerCoords, u1diff, tmp_strain, tmp_displacement);
                  for(unsigned int kl = ij; kl < ex_size; ++kl){
                    u2diff = *dynamic_cast<Vector<double>* >(u2.Get(kl)->gridelem[MECH][e]);
                    SubtractTestDisplacement(kl, CornerCoords, u2diff, tmp_strain, tmp_displacement);
                    // description see above, is needed twice for SPEED
                    double v1 = 0.0;
                    for(unsigned int i = 0; i < r; ++i){
                      double v2 = 0.0;
                      for(unsigned int j = 0; j < c; ++j){
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
  
  shapedesign->AddShapeDerivatives(f, constraint, der, 1.0);
}

void ShapeOpt::CalcUdF(Excitation& excite, StdVector<SingleVector*>& u, Objective* f, Condition* constraint, double w){
  StdVector<double> der; // solution
  int np = shapedesign->GetNumberOfShapeParameters();
  der.Resize(np, 0);

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

      Vector<double>& uelem = dynamic_cast<Vector<double>& >(*u[e]);
      double pres = form->GetPressureFactor(elem);
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
        double f = intWeight*pres;
        double v = uelem * elemVec;

        for(int p = 0; p < np; p++) { // loop over all parameters
          shapedesign->GetElemNodesCoordDerivative(dCornerCoords, elem->connect, p);
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

        } // parameter loop
      } // integration point loop
    } // element loop
  } // linForm loop

  shapedesign->AddShapeDerivatives(f, constraint, der, w);
}

double ShapeOpt::CalcCompliance(Excitation& excite, Objective* f, Condition* constraint, bool derivative){
  if(derivative){
    // the derivative of tracking w.r.t. shape is: - u' dA/dShape u + 2 u dF/dShape
    if(excite.index == 0){      
      CalcMinusU1dKU2(forward, forward, f, constraint);
    }
    // this however is done per excite
    StdVector<SingleVector*>& u = forward.Get(excite)->gridelem[MECH];
    CalcUdF(excite, u, f, constraint, 2*excite.weight);
    if(alsomatopt_){
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
    if(excite.index == 0){ // CalcMinusU1dKU2 runs over the excitations, this speeds up things a little
      CalcMinusU1dKU2(adjoint, forward, f, constraint);
    }
    // this however is done per excite
    StdVector<SingleVector*>& u = forward.Get(excite)->gridelem[MECH];
    CalcUdF(excite, u, f, constraint, excite.weight);
    
    if(alsomatopt_){
      ErsatzMaterial::CalcTracking(excite, f, constraint, true);
    }
  }else{
    return(ErsatzMaterial::CalcTracking(excite, f, constraint, derivative));
  }
  return 0.0;
}

void ShapeOpt::CalcHomogenizedTrackingGradient(const Matrix<double>& target, const Matrix<double>& hom, Objective* f, Condition* g){
  Matrix<double> tensor_diff;
  tensor_diff = hom - target;
  CalcMinusU1dKU2(forward, forward, f, g, &tensor_diff);
}

Matrix<double> ShapeOpt::CalcHomogenizedTensor(){
  const unsigned int ex_size(me->excitations.GetSize());
  assert((dim == 2 && ex_size == 3) || (dim == 3 && ex_size == 6));
  
  double rcubevol(1.0 / grid->CalcVolumeSpannedByNamedNodes());
  
  Matrix<double> result(ex_size, ex_size);
  result.Init();

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
    if(biLinForm->GetIntegrator()->GetName() != "linElastInt") continue;
    linElastInt* form = (linElastInt*)(biLinForm->GetIntegrator());
    EntityIterator it = biLinForm->GetFirstEntities()->GetIterator();
    for(it.Begin(); !it.IsEnd(); it++){ // loop over all corresponding elements
      const Elem* elem = it.GetElem();
      int e = elem->elemNum - 1;
      grid->GetElemNodesCoord(CornerCoords, elem->connect, true);
      form->CalcElementMatrix(elemMat, it, it);
      for(unsigned int ij = 0; ij < ex_size; ++ij){
        u1diff = *dynamic_cast<Vector<double>* >(forward.Get(ij)->gridelem[MECH][e]);
        SubtractTestDisplacement(ij, CornerCoords, u1diff, tmp_strain, tmp_displacement);
        Ku = elemMat * u1diff;
        for(unsigned int kl = ij; kl < ex_size; ++kl){ // only upper triangle
          u2diff = *dynamic_cast<Vector<double>* >(forward.Get(kl)->gridelem[MECH][e]);
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

  return result;
}

void ShapeOpt::StorePDESolution(Solutions& solutions, Excitation &excite, Function* f, UInt timestep, bool read_sol, bool read_rhs, bool save_sol, const std::string& comment)
{
  ParamMat::StorePDESolution(solutions, excite, f, timestep, read_sol, read_rhs, save_sol, comment);
  if(read_sol){
    solutions.Get(excite, f, timestep)->Read(Solution::GRIDELEM_VECTORS, pde, MECH);
  }
}

void ShapeOpt::SubtractTestDisplacement(unsigned int idx, Matrix<double>& CornerCoords, Vector<double>& result, Matrix<double>& tmp_strain, Matrix<double>& tmp_displacement){
  SetTestStrainMatrix(tmp_strain, me->excitations[idx].test_strain);
  unsigned int cols = CornerCoords.GetNumCols();
  tmp_displacement.Resize(dim, cols);
  tmp_strain.Mult(CornerCoords, tmp_displacement);
  for(unsigned int out = 0; out < cols; ++out){
    for(unsigned int in = 0; in < dim; ++in){
      result[out*dim + in] -= tmp_displacement[in][out];
    }
  }
}

