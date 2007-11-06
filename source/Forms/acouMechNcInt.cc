#include "acouMechNcInt.hh"
#include "Domain/ncElem.hh"
#include "Domain/grid.hh"
#include "Domain/domain.hh"

namespace CoupledField {

  AcouMechNcInt::AcouMechNcInt(UInt dofsPerNode, bool isAxi) {
    name_ = "AcouMechNcInt";
    dofs_ = dofsPerNode;
    isaxi_ = isAxi;
    isSymmetric_ = false;
    formulation_ = NO_SOLUTION_TYPE;
  }

  AcouMechNcInt::~AcouMechNcInt() {
  }

  void AcouMechNcInt::SetFormulation( SolutionType aformulation) {
    if ((aformulation != ACOU_POTENTIAL) && (aformulation != ACOU_PRESSURE)) {
      EXCEPTION(name_ << ": formulation may only be set to ACOU_POTENTIAL or ACOU_PRESSURE");
      return;
    }

    formulation_ = aformulation;
  }
    
  void AcouMechNcInt::CalcElementMatrix(Matrix<Double> &elemMat,
      EntityIterator &ent1,
      EntityIterator &ent2)
  {
    // make sure that formulation is set
    if (formulation_ == NO_SOLUTION_TYPE)
      EXCEPTION(name_
          << ": formulation must be set before calcution of element matrix");

    // cast to NCElem
    actElem_ = ent1.GetSurfElem();
    const NCElem *actNCElem = dynamic_cast< const NCElem* >(actElem_);
    if (actNCElem == NULL)
      EXCEPTION(name_ << " can only handle NCElems.");

    bool fluidOnSlaveSide, isCoplanar;
    Double fluidDensity, jacDet, coordAtIp;
    Grid *ptGrid = domain->GetGrid();
    Matrix<Double> helpMat, partHelpMat, globalIP, localMasterIP,
      localSlaveIP, matNodeCoord;
    std::map<RegionIdType, BaseMaterial*> *acouMaterials;
    std::map<RegionIdType, BaseMaterial*>::iterator itMat;
    SurfElem *elemMaster = actNCElem->ptSurfParent,
             *elemSlave = actNCElem->ptLagrangeParent;
    UInt dim = ptGrid->GetDim(), i, j,
         numMasterShpFncs, numSlaveShpFncs, noIntPoints;
    Vector<Double> point, *intPoints, intWeights,
      shpFncMaster, shpFncSlave, shpFncNCElem, tmp;

    // extract node coordinates
    ptGrid->GetElemNodesCoord(ptCoord_, actElem_->connect);
    // calculate surface normal
    ptGrid->CalcSurfNormal(normal_, *elemSlave);
    normal_ *= (Double) elemSlave->normalSign;

    // is our interface coplanar?
    isCoplanar = ptGrid->IsNcInterfaceCoplanar(actElem_->regionId);

    // get pointer to materials of acoustic pde
    if (firstPDEName_ == "acoustic")
      acouMaterials = &firstMaterials_;
    else
      acouMaterials = &secondMaterials_;
    // does slave element belong to acoustic domain?
    itMat = acouMaterials->find(elemSlave->ptVolElem1->regionId);
    if (itMat == acouMaterials->end()) { // no
      // get material of master element (in acoustic domain)
      itMat = acouMaterials->find(elemMaster->ptVolElem1->regionId);
      fluidOnSlaveSide = false;
    } else { // reverse surface normal (must point into acoustic domain)
      normal_ *= -1.0;
      fluidOnSlaveSide = true;
    }
    // get the density of the fluid
    itMat->second->GetScalar(fluidDensity, DENSITY, REAL);
    if (formulation_ == ACOU_PRESSURE) {
      if (firstPDEName_ == "mechanic")
        fluidDensity = 1.0;
      else
        fluidDensity *= -fluidDensity;
    }

    // obtain no. of shape functions
    elemMaster->ptElem->SetAnsatzFct(ansatzFct1_);
    numMasterShpFncs = elemMaster->ptElem->GetNumFncs(ansatzFct1_);

    elemSlave->ptElem->SetAnsatzFct(ansatzFct2_);
    numSlaveShpFncs = elemSlave->ptElem->GetNumFncs(ansatzFct2_);

    // obtain integration points and weights
    noIntPoints = actNCElem->ptElem->GetNumIntPoints();
    intPoints = actNCElem->ptElem->GetIntPoints();
    intWeights = actNCElem->ptElem->GetIntWeights();

    // resize coordinate matrices
    globalIP.Resize(dim, noIntPoints);
    localMasterIP.Resize(dim - 1, noIntPoints);
    localSlaveIP.Resize(dim - 1, noIntPoints);
    point.Resize(dim);

    // transform integration points to global coordinates
    for (i = 0; i < noIntPoints; ++i) {
      actNCElem->ptElem->Local2GlobalCoord(point, intPoints[i],
          ptCoord_, NULL);
      for (j = 0; j < dim; ++j)
        globalIP[j][i] = point[j];
    }

    // transform global to local coordinates of slave element
    ptGrid->GetElemNodesCoord(matNodeCoord, elemSlave->connect);
    elemSlave->ptElem->Global2LocalCoords(localSlaveIP, globalIP,
        matNodeCoord);

    // transform global to local coordinates of master element
    ptGrid->GetElemNodesCoord(matNodeCoord, elemMaster->connect);
    if (!isCoplanar) { // projection is necessary
      Double scale, dist, sign;
      Vector<Double> normal, p0, line;

      // calculate surface normal of master element
      ptGrid->CalcSurfNormal(normal, *elemMaster);

      // get first node of master element
      p0.Resize(dim);
      for (i = 0; i < dim; ++i)
        p0[i] = matNodeCoord[i][0];

      // project integration points onto master element
      for (i = 0; i < noIntPoints; ++i) {
        // get global coordinates of i-th integration point
        for (j = 0; j < dim; ++j)
          point[j] = globalIP[j][i];

        line = point - p0;
        
        // calculate distance of int. point to master element
        normal.Inner(line, dist);

        // determine orientation of normal on slave element
        normal_.Inner(line, sign);
        sign /= fabs(sign);
        
        // scale the distance for the normal projection
        normal_.Inner(normal, scale);
        
        // do the projection
        for (j = 0; j < dim; ++j) {
          globalIP[j][i] -= sign * normal_[j] * fabs(dist) * fabs(scale);
        }
      }
    }
    elemMaster->ptElem->Global2LocalCoords(localMasterIP, globalIP,
        matNodeCoord);

    // resize matrix according to master/slave config
    if (fluidOnSlaveSide)
      helpMat.Resize(numMasterShpFncs, numSlaveShpFncs);
    else
      helpMat.Resize(numSlaveShpFncs, numMasterShpFncs);
    helpMat.Init();

    // first calculate single-DoF helper matrix
    for (i = 0; i < noIntPoints; ++i) {
      // get values of shape functions
      for (j = 0; j < dim - 1; ++j)
        point[j] = localMasterIP[j][i];
      try {
        elemMaster->ptElem->GetShFnc(shpFncMaster, point, elemMaster);
      } catch (Exception& ex) {
        (*warning) << ex.GetMsg() << "\n NCElem: " << actNCElem->elemNum
          << "\tmaster: " << elemMaster->elemNum
          << "\tslave: " << elemSlave->elemNum;
        Warning(__FILE__, __LINE__);
      }

      for (j = 0; j < dim - 1; ++j)
        point[j] = localSlaveIP[j][i];
      try {
        elemSlave->ptElem->GetShFnc(shpFncSlave, point, elemSlave);
      } catch (Exception& ex) {
        (*warning) << ex.GetMsg() << "\n NCElem: " << actNCElem->elemNum
          << "\tmaster: " << elemMaster->elemNum
          << "\tslave: " << elemSlave->elemNum;
        Warning(__FILE__, __LINE__);
      }

      // calculate normal mass matrix
      if (fluidOnSlaveSide)
        partHelpMat.DyadicMult(shpFncMaster, shpFncSlave);
      else
        partHelpMat.DyadicMult(shpFncSlave, shpFncMaster);

      // calculate Jacobian determinant
      jacDet = actNCElem->ptElem->CalcJacobianDet(intPoints[i],
          ptCoord_, actNCElem);

      if (isaxi_) {
        for (j = 0; j < dim - 1; ++j)
          point[j] = intPoints[i][j];
        actNCElem->ptElem->GetShFnc(shpFncNCElem, point, actNCElem);
        coordAtIp = 0.0;
        for (j = 0; j < noIntPoints; ++j)
          coordAtIp += ptCoord_[0][j] * shpFncNCElem[j];
        partHelpMat *= 2 * PI * intWeights[i] * fluidDensity *
            std::fabs(jacDet) * coordAtIp;
      }
      else 
      { // multiply matrix by density of fluid
        partHelpMat *= intWeights[i] * fluidDensity * std::fabs(jacDet);
      }
      
      helpMat += partHelpMat;
    }

    // generate multi-dof matrix (scalar product with surface normal)
    if (firstPDEName_ == "acoustic") {
      elemMat.Resize( numMasterShpFncs, numSlaveShpFncs * dofs_ );
      for ( UInt iRow = 0; iRow < numMasterShpFncs; iRow++ ) {
        for ( UInt iCol = 0; iCol < numSlaveShpFncs; iCol++ ) {
          for ( UInt iDof = 0; iDof < dofs_; iDof++ ) {
            elemMat[iRow][iCol * dofs_ + iDof] =
              normal_[iDof] * helpMat[iCol][iRow];
          }
        }
      }
    } else {
      elemMat.Resize( numMasterShpFncs * dofs_, numSlaveShpFncs );
      for ( UInt iRow = 0; iRow < numMasterShpFncs; iRow++ ) {
        for ( UInt iCol = 0; iCol < numSlaveShpFncs; iCol++ ) {
          for ( UInt iDof = 0; iDof < dofs_; iDof++ ) {
            elemMat[iRow * dofs_ + iDof][iCol] =
              normal_[iDof] * helpMat[iRow][iCol];
          }
        }
      }
    }
  }

} // namespace CoupledField
