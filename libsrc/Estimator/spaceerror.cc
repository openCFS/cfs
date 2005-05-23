#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "Domain/grid.hh"
#include "spaceerror.hh"
 
#include "AlgebraicSystem/interface_piles.hh"
#include "Forms/forms_header.hh"
#include "AlgebraicSystem/LinAlg/linsystem.hh"

namespace CoupledField
{

  SpaceErrorEstimator::SpaceErrorEstimator()
  {
    ;
  }

  void SpaceErrorEstimator::CalcErrorMap(const Vector<Double> & sol,
                                         std::vector<std::string> &subdoms,
                                         Grid * ptgrid,
                                         Vector<Double> & errorMap,
                                         Double & atotalErr, const Integer level)
  {
#ifdef TRACE
    *trace << "entering SpaceErrorEstimation::CalcErrorMap" << std::endl;
#endif
  
    // loop over elements
    Vector<Double>          result;
    std::vector<Elem*>      elemssd;
    std::vector<Elem*>      * neighbours;
    //vector, where we store global number of nodes for which result is obtained in fnc RecoveryProcedure4ElemsPatch
    std::vector<Integer>    locationsInResult;
    Integer                 dim=ptgrid->GetDim();
    Vector<Double>          areaElms;
    //vector, where we store SPR gradient
    Vector<Double>          * SPRgrad=new Vector<Double>[dim];
    //vector, where we store number of summered value of SPR at node. this vector is needed for averaging procedure of SPRgrad
    Vector<Integer>         * numberAvergVals=new Vector<Integer>[dim];
    Integer                 nrNodes = ptgrid->GetMaxnumnodes(level);
    Integer                 maxelem = ptgrid->GetMaxnumElem(level,subdoms);
    Integer                 idm,icmp,isd, ind, jel, ires, ic, igr, iel;
    Integer                 nrSubdoms = subdoms.size();
    Integer                 size=sol.size();
    std::string             type_patch;
    Integer                 counterElems = 0;
    Vector<Double>          gradSPRElemL2norm;
    Double                  sumErrMap = 0;
    Double                  glNrmErr = 0;
    Double                  areaElm, sumArea=0;
    Double                  errEl;
    Double                  normGradSPR;

    errorMap.Resize(maxelem);
    gradSPRElemL2norm.Resize(maxelem);
    areaElms.Resize(maxelem);

    for (idm=0;idm<dim;idm++)
      {
        SPRgrad[idm].Resize(size);
        numberAvergVals[idm].Resize(size);
      }
  
    // computation of SPR gradient
    if (conf->ifget("type_of_patch",type_patch) && type_patch=="nodebased")
      {
#ifdef DEBUG
        (*debug) << " step: error map : node based patches are used \n";
#endif
        Warning("Node based patch is working only for 1 subdomain",__FILE__,__LINE__);

        for (icmp=0; icmp<dim; icmp++)
          { // loop over components of gradient
            for (ind=0; ind< nrNodes; ind++)
              {
                // get patch
                ptgrid->GetNeighboursOfNode(ind, neighbours);
                // 
                RecoveryProcedure4ElemsPatch(*neighbours,ptgrid,
                                             sol,icmp,
                                             result,locationsInResult, level);
              
                // arrays for averaging values
                // form arrays for averaging values of SPRgrad over nodes
                for (ires=0;ires<result.size();ires++)
                  {
                    SPRgrad[icmp][locationsInResult[ires]-1]+=result[ires];
                    numberAvergVals[icmp][locationsInResult[ires]-1]++;
                  }              
              }
          } // end: loop over components of gradient
      } // end of if: nodebased patch
    // elementbased
    else
      {
#ifdef DEBUG
        (*debug) << " step: error map : element based patches are used \n";
#endif  
        for (icmp=0; icmp<dim; icmp++)
          { // loop over components of gradient     
            for (isd=0; isd< nrSubdoms; isd++) // loop over all subdomains
              { 
                ptgrid->GetElemSD(elemssd,subdoms[isd],level);
              
                for (jel=0; jel < elemssd.size(); jel++) // loop over elements of subdomain
                  {
                    //  std::vector<Elem*> * neighbours;
                    // in: j - noOfElem, suddoms_[isd] - name of subdomain;
                    // out: pt to vector with neughbors-element
                    neighbours=ptgrid->GetNeighboursOfElem(jel,subdoms[isd]);
                  
                    // add element in list of neighbors to get a full patch of elements
                    neighbours->push_back(elemssd[jel]);

                    // in: elemssd[j] - list with elements of this subdomain
                    RecoveryProcedure4ElemsPatch(*neighbours,ptgrid,
                                                 sol,icmp,result,
                                                 locationsInResult,level);

                    // form arrays for averaging values of SPRgrad over nodes
                    Integer nrVals = result.size();
                    for (ires=0; ires< nrVals;ires++)
                      {
                        SPRgrad[icmp][locationsInResult[ires]-1]+=result[ires];
                        numberAvergVals[icmp][locationsInResult[ires]-1]++;
                      }
         
                  } // loop over elements of subdomains
              } // loop over subdomains
          } // end: loop over component of gradient
      } // end of else: element patch
  
    //average values of SPRgrad over nodes
    for (ic=0; ic<dim; ic++)
      { // loop over components of gradient         
        for (igr=0; igr<nrNodes;igr++) 
          SPRgrad[ic][igr]/=numberAvergVals[ic][igr];     
      } // end of loop over components of gradient
  
    // delete 
    delete [] numberAvergVals;                             
  
    //   std::cout << " SPR grad: " << std::endl;
    //    for (igr=0; igr<sol.size();igr++)
    //     std::cout << SPRgrad[0][igr] << " " <<  SPRgrad[1][igr] << std::endl;

    // calculation of error for each element
    for (isd=0; isd< nrSubdoms; isd++) // loop over all subdomains
      {      
        // get vector of Elem of subdomain with color: subdoms[isd]
        ptgrid->GetElemSD(elemssd,subdoms[isd],level);
       
        // loop over elements of subdomain
        Integer nrElemsSd = elemssd.size();
        for (iel=0; iel< nrElemsSd; iel++, counterElems++)
          {
            CalcErrorForElem(elemssd[iel], SPRgrad, errEl, normGradSPR,
                             sol, ptgrid, level);
           
            // calculation of area of element
            areaElm=CalcArea(elemssd[iel], ptgrid,level);
            sumArea+=areaElm;
 
            gradSPRElemL2norm[counterElems]=normGradSPR;
            glNrmErr+=normGradSPR;
            sumErrMap+=errEl;
            errorMap[counterElems]=errEl/areaElm;
          } // end of loop over elems in subdomain
      } // end of loop over subdomains

    delete [] SPRgrad;

    atotalErr=sumErrMap/glNrmErr;
    errorMap/=glNrmErr;
    errorMap*=sumArea;   
  
    //   std::cout << errorMap << std::endl;
    //   exit(1);

#ifdef TRACE
    *trace << "leaving SpaceErrorEstimator::CalcErrorMap" << std::endl;
#endif    
  }

  void SpaceErrorEstimator::CalcErrorForElem(const Elem* elem,
                                             const Vector<Double>* SPRgrad,  Double & error, Double & normGradSPR,
                                             const Vector<Double> & sol, Grid * ptgrid, const Integer level)
  {
#ifdef TRACE
    (*trace) << "entering SpaceErrorEstimator::CalcErrorForElem" << std::endl;
#endif
#ifdef DEBUG
    (*debug) << " step: error map \n";
#endif

    BaseFE          * ptelem = elem->ptElem;
    Vector<Integer> connectVec = elem->connect;
    Integer         nrNodes = ptelem->GetNumNodes();
    const Integer   nrIntPnts = ptelem->GetNumIntPoints();
    const Integer   dim = ptelem->GetDim();
    Matrix<Double>       dvShFnc;
    const std::vector<Double> & intWeights = ptelem->GetIntWeights(); 
    std::vector<double> shFnc;
    Matrix<Double>      ptCoord;
    Double              valsol;
    Double              jacDet;
    Double              diffGradL2norm;
    Double              SPRGradL2norm; 


    error = 0;
    normGradSPR = 0;

    ptgrid->GetCoordNodesElemMat(connectVec,ptCoord,level); 

    //do a loop over int. points in element
    Integer iIntPnts;
    for (iIntPnts=0; iIntPnts<nrIntPnts; iIntPnts++)
      {   
     
        jacDet = ptelem->CalcJacobianDetAtIp(iIntPnts+1, ptCoord);
     
        // calculation of FEM gradient and SPR at intgr. points in local coord.
        Double     gradFEM[dim];
        Double     gradSPR[dim];
        Integer    idm,iShFnc;
      
        for (idm=0; idm<dim; idm++)
          {
            gradFEM[idm]=0;
            gradSPR[idm]=0;
          }
     
        ptelem->GetGlobDerivShFncAtIp(dvShFnc, iIntPnts + 1, ptCoord, jacDet); 
        ptelem->GetShFncAtIp(shFnc,iIntPnts+1); 
          
        for (iShFnc = 0; iShFnc < nrNodes; iShFnc++)
          for (idm = 0; idm < dim; idm++)
            {
              gradFEM[idm]+=sol[connectVec[iShFnc]-1]*dvShFnc[iShFnc][idm];
              gradSPR[idm]+=SPRgrad[idm][connectVec[iShFnc]-1]*shFnc[iShFnc];
            }
          
#ifdef DEBUG
#ifdef DEBUGFULL     
        (*debug) << gradSPR[0] << " " << gradFEM[0] << " y: " << gradSPR[1] << " " << gradFEM[1] << std::endl;
        if (dim==3) (*debug) << " z: " << gradSPR[2] << " " << gradFEM[2] << std::endl;
#endif
#endif

        // calculation of L2 norm for difference of SPR and FEM gradient,
        // and L2 norm of SPR gradient
        diffGradL2norm = 0;
        SPRGradL2norm  = 0; 

        for (idm = 0; idm<dim; idm++)
          {
            diffGradL2norm+=sqr(gradSPR[idm]-gradFEM[idm]);
            SPRGradL2norm+=sqr(gradSPR[idm]);
          }
     
        diffGradL2norm=sqrt(diffGradL2norm);
        SPRGradL2norm=sqrt(SPRGradL2norm);

        error+=diffGradL2norm*jacDet*intWeights[iIntPnts];
        normGradSPR+=SPRGradL2norm*jacDet*intWeights[iIntPnts];

      } // end of loop over integration points

#ifdef TRACE
    (*trace) << "leaving SpaceErrorEstimator::CalcErrorForElem" << std::endl;
#endif
  }

  void SpaceErrorEstimator::RecoveryProcedure4ElemsPatch(const std::vector<Elem*> &Elems,
                                                         Grid * ptgrid,
                                                         const Vector<Double> & sol, 
                                                         const Integer aComponent,
                                                         Vector<Double> &result,
                                                         std::vector<Integer>&locations,
                                                         const Integer level)
  {
#ifdef TRACE
    (*trace) << "entering SpaceErrorEstimator::RecoveryProcedure4ElemsPatch " << std::endl;
#endif
    
    Vector<Integer>  vec_connect; // auxilialary vector, stores connection of element
    Matrix<Double>   elemmat; // sys matrix for element
    Vector<Double>   vecRHS; // rhs vector for element
    Integer          counter=0; // counter nodes in element-patch
    Matrix<Double>   SysMatrix; // system matrix for elements-patch 
    Vector<Double>   RHSVector; // rhs vector for elements-patch
   
    // define dimension of system matrix and rhs 
    Integer          nnodes=CalcNumberOfNodesInPatch(Elems);
    //    Warning("For mass element matrices material data is not defined",__FILE__,__LINE__);
    SysMatrix.Resize(nnodes,nnodes);
    RHSVector.Resize(nnodes);
    // in this vector we store global numbers of nodes that occure in the patch.
    // number of position is ,accordingly, local number in patch
    std::vector<Integer>  mapGlbNumNodes2LocInPatch;
    mapGlbNumNodes2LocInPatch.resize(nnodes);

    Integer          iel;
    Integer          nrElems = Elems.size();
    for(iel=0; iel< nrElems; iel++)
      {   // loop over elements in patch 

        // element system matrix
        BaseFE         * ptBaseFE = Elems[iel]->ptElem;

        BaseForm       * bilinear_mass  = new MassInt(ptBaseFE,1);
        Matrix<Double> ptCoord;
        vec_connect = Elems[iel]->connect;

        Integer        nrNodesElem = ptBaseFE->GetNumNodes();

        ptgrid->GetCoordNodesElemMat(vec_connect,ptCoord,level); 
        bilinear_mass->CalcElementMatrix(ptCoord, elemmat);

        // element rhs
        RHSForRecoveryProcedure * rhs=new RHSForRecoveryProcedure(ptBaseFE);
        Vector<Double> valFnc(nrNodesElem);  // get val of sol at nodes of the element
    
        Integer in;
        for (in=0; in<nrNodesElem; in++) 
          valFnc[in]=sol[vec_connect[in]-1];

        rhs->CalcElemVectorRHSForSPR(ptCoord,valFnc,aComponent,vecRHS);
        
        // do renumbering of nodes for this element for assembling procedure
        Boolean UseOldNumber;
   
        std::vector<Integer> assembleConnect(vec_connect.size());
        Integer ivc,imp;
        // do a loop over connection vector 
        for (ivc=0; ivc<vec_connect.size(); ivc++)
          {    
            UseOldNumber=FALSE;

            // loop over vector with global nodes for previous elements
            for (imp=0; imp<mapGlbNumNodes2LocInPatch.size(); imp++)
              {
                // check that this node is not new
                if (mapGlbNumNodes2LocInPatch[imp] == vec_connect[ivc])
                  {
                    // if node is new, then we use previous local number in assembling
                    assembleConnect[ivc] = imp;      
                    UseOldNumber=TRUE;
                  }
              }
            
            if (!UseOldNumber)
              {
                assembleConnect[ivc]=counter;
                mapGlbNumNodes2LocInPatch[counter]=vec_connect[ivc];
                counter++;
              }
          }
        
        // assemble sysmatrix
        Integer irow,icln;
        Integer ii, iii;
        for (ii=0; ii< nrNodesElem; ii++)
          for (iii=0; iii< nrNodesElem; iii++)
            {
              irow=assembleConnect[ii];
              icln=assembleConnect[iii];
              SysMatrix.Add(irow,icln,elemmat[ii][iii]);
            } 
   
        // assemble rhs
        Integer sizeRHS = vecRHS.size();
        for (ii=0; ii< sizeRHS;ii++) {
          irow=assembleConnect[ii];
          RHSVector[irow]+=vecRHS[ii];
        }    

        delete bilinear_mass;
        delete rhs;
    
      } // end of loop over elements in patch

    // solve lin system
    Double tol=1e-8;
    LinSystem<Double,Matrix<Double> > LinSys(tol);
    LinSys.SetSysMatrix(SysMatrix);
    LinSys.SetRHS(RHSVector);

    // LinSys.BiCGSTAB(100,Jacobi);
    // LinSys.CG(100,Jacobi);
    LinSys.CholerskyMethod();
    LinSys.GetSolution(result);

    locations=mapGlbNumNodes2LocInPatch;

#ifdef TRACE
    (*trace) << "leaving SpaceErrorEstimator::RecoveryProcedure4ElemsPatch " << std::endl;
#endif

  }  // end of fnc: RecoveryProcedure4ElemsPatch

  Integer SpaceErrorEstimator::CalcNumberOfNodesInPatch(const std::vector<Elem*> & patch)
  {
#ifdef TRACE
    (*trace) << "entering Elecst2dPDE::CalcNumberOfNodesInPatch" << std::endl;
#endif

    Integer              iels,ivc,imp;
    Vector<Integer>      vec_connect;
    Boolean              NewNode;
    std::vector<Integer> globalNumbers; // store global numbers of nodes of element
    // that occurs in patch

    for (iels=0; iels<patch.size(); iels++) // loop over elements in patch
      {     
        vec_connect=patch[iels]->connect;
      
        for (ivc=0; ivc<vec_connect.size(); ivc++)
          { 
            NewNode=TRUE;
            // loop over vector with global nodes for previous elements
            for (imp=0; imp<globalNumbers.size(); imp++)
              {
                // check that this node is not new
                if (globalNumbers[imp] == vec_connect[ivc])  
                  NewNode=FALSE;
          
              }
          
            if (NewNode) 
              globalNumbers.push_back(vec_connect[ivc]);
          
          } // end of loop over nodes in element

      } // end of loop over elements in patch

    return globalNumbers.size();
     
  }

  void SpaceErrorEstimator::CalcErrorMapHarmonic(const Vector<Double> & sol,
                                                 const Vector<Double> & solIm,
                                                 std::vector<std::string> &subdoms,
                                                 Grid * ptgrid,
                                                 Vector<Double> & errorMap,
                                                 Double & atotalErr, const Integer level)
  {
#ifdef TRACE
    *trace << "entering SpaceErrorEstimation::CalcErrorMapHarmonic" << std::endl;
#endif
  
    // loop over elements
    Vector<Double>          result, resultIm;
    std::vector<Elem*>      elemssd;
    std::vector<Elem*>      * neighbours;
    //vector, where we store global number of nodes for which result is obtained in fnc RecoveryProcedure4ElemsPatch
    std::vector<Integer>    locationsInResult;
    Integer                 dim=ptgrid->GetDim();
    Vector<Double>          areaElms;
    //vector's, where we store SPR gradient
    Vector<Double>          * SPRgrad=new Vector<Double>[dim];
    Vector<Double>          * SPRgradIm=new Vector<Double>[dim];
    //vector, where we store number of summered value of SPR at node. this vector is needed for averaging procedure of SPRgrad
    Vector<Integer>         * numberAvergVals=new Vector<Integer>[dim];
    Integer                 nrNodes = ptgrid->GetMaxnumnodes(level);
    Integer                 maxelem = ptgrid->GetMaxnumElem(level,subdoms);
    Integer                 idm,icmp,isd, ind, jel, ires, ic, igr, iel;
    Integer                 nrSubdoms = subdoms.size();
    Integer                 size=sol.size();
    std::string             type_patch;
    Integer                 counterElems = 0;
    Vector<Double>          gradSPRElemL2norm;
    Double                  sumErrMap = 0;
    Double                  glNrmErr = 0;
    Double                  areaElm, sumArea=0;
    Double                  errEl;
    Double                  normGradSPR;

    errorMap.Resize(maxelem);
    gradSPRElemL2norm.Resize(maxelem);
    areaElms.Resize(maxelem);

    for (idm=0;idm<dim;idm++)
      {
        SPRgrad[idm].Resize(size);
        SPRgradIm[idm].Resize(size);
        numberAvergVals[idm].Resize(size);
      }
  
    // computation of SPR gradient
    if (conf->ifget("type_of_patch",type_patch) && type_patch=="nodebased")
      {
#ifdef DEBUG
        (*debug) << " step: error map : node based patches are used \n";
#endif
        Warning("Node based patch is working only for 1 subdomain",__FILE__,__LINE__);

        for (icmp=0; icmp<dim; icmp++)
          { // loop over components of gradient
            for (ind=0; ind< nrNodes; ind++)
              {
                // get patch
                ptgrid->GetNeighboursOfNode(ind, neighbours);
                // 
                RecoveryProcedure4ElemsPatch(*neighbours,ptgrid,
                                             sol,icmp,
                                             result,locationsInResult, level);

                RecoveryProcedure4ElemsPatch(*neighbours,ptgrid,
                                             solIm,icmp,
                                             resultIm,locationsInResult, level);

                // arrays for averaging values
                // form arrays for averaging values of SPRgrad over nodes
                for (ires=0;ires<result.size();ires++)
                  {
                    SPRgrad[icmp][locationsInResult[ires]-1]+=result[ires];
                    SPRgradIm[icmp][locationsInResult[ires]-1]+=resultIm[ires];
                    numberAvergVals[icmp][locationsInResult[ires]-1]++;
                  } 

              }
          } // end: loop over components of gradient
      } // end of if: nodebased patch
    // elementbased
    else
      {
#ifdef DEBUG
        (*debug) << " step: error map : element based patches are used \n";
#endif  
        for (icmp=0; icmp<dim; icmp++)
          { // loop over components of gradient     
            for (isd=0; isd< nrSubdoms; isd++) // loop over all subdomains
              { 
                ptgrid->GetElemSD(elemssd,subdoms[isd],level);
              
                for (jel=0; jel < elemssd.size(); jel++) // loop over elements of subdomain
                  {
                    //  std::vector<Elem*> * neighbours;
                    // in: j - noOfElem, suddoms_[isd] - name of subdomain;
                    // out: pt to vector with neughbors-element
                    neighbours=ptgrid->GetNeighboursOfElem(jel,subdoms[isd]);
                  
                    // add element in list of neighbors to get a full patch of elements
                    neighbours->push_back(elemssd[jel]);

                    // in: elemssd[j] - list with elements of this subdomain
                    RecoveryProcedure4ElemsPatch(*neighbours,ptgrid,
                                                 sol,icmp,result,
                                                 locationsInResult,level);

                    RecoveryProcedure4ElemsPatch(*neighbours,ptgrid,
                                                 solIm,icmp,resultIm,
                                                 locationsInResult,level);

                    // form arrays for averaging values of SPRgrad over nodes
                    Integer nrVals = result.size();
                    for (ires=0; ires< nrVals;ires++)
                      {
                        SPRgrad[icmp][locationsInResult[ires]-1]+=result[ires];
                        SPRgradIm[icmp][locationsInResult[ires]-1]+=resultIm[ires];
                        numberAvergVals[icmp][locationsInResult[ires]-1]++;
                      }
         
                  } // loop over elements of subdomains
              } // loop over subdomains
          } // end: loop over component of gradient
      } // end of else: element patch
  
    //average values of SPRgrad over nodes
    for (ic=0; ic<dim; ic++)
      { // loop over components of gradient         
        for (igr=0; igr<nrNodes;igr++)
          {
            SPRgrad[ic][igr]/=numberAvergVals[ic][igr];     
            SPRgradIm[ic][igr]/=numberAvergVals[ic][igr]; 
          }
      } // end of loop over components of gradient
  
    // delete 
    delete [] numberAvergVals;                             
  
    //   std::cout << " SPR grad: " << std::endl;
    //    for (igr=0; igr<sol.size();igr++)
    //     std::cout << SPRgrad[0][igr] << " " <<  SPRgrad[1][igr] << std::endl;

    // calculation of error for each element
    for (isd=0; isd< nrSubdoms; isd++) // loop over all subdomains
      {      
        // get vector of Elem of subdomain with color: subdoms[isd]
        ptgrid->GetElemSD(elemssd,subdoms[isd],level);
       
        // loop over elements of subdomain
        Integer nrElemsSd = elemssd.size();
        for (iel=0; iel< nrElemsSd; iel++, counterElems++)
          {
            CalcErrorForElemHarmonic(elemssd[iel], SPRgrad, SPRgradIm,
                                     errEl, normGradSPR, sol, solIm, 
                                     ptgrid, level);
           
            // calculation of area of element
            areaElm=CalcArea(elemssd[iel], ptgrid,level);
            sumArea+=areaElm;
 
            gradSPRElemL2norm[counterElems]=normGradSPR;
            glNrmErr+=normGradSPR;
            sumErrMap+=errEl;
            errorMap[counterElems]=errEl/areaElm;
          } // end of loop over elems in subdomain
      } // end of loop over subdomains

    delete [] SPRgrad;
    delete [] SPRgradIm;

    atotalErr=sumErrMap/glNrmErr;
    errorMap/=glNrmErr;
    errorMap*=sumArea;   
  
#ifdef TRACE
    *trace << "leaving SpaceErrorEstimator::CalcErrorMapHarmonic" << std::endl;
#endif    
  }

  void SpaceErrorEstimator::CalcErrorForElemHarmonic(const Elem* elem,
                                                     const Vector<Double>* SPRgrad,  const Vector<Double>* SPRgradIm,
                                                     Double & error, Double & normGradSPR,
                                                     const Vector<Double> & sol, const Vector<Double> & solIm,                                              
                                                     Grid * ptgrid, const Integer level)
  {
#ifdef TRACE
    (*trace) << "entering SpaceErrorEstimator::CalcErrorForElem" << std::endl;
#endif
#ifdef DEBUG
    (*debug) << " step: error map \n";
#endif

    BaseFE          * ptelem = elem->ptElem;
    Vector<Integer> connectVec = elem->connect;
    Integer         nrNodes = ptelem->GetNumNodes();
    const Integer   nrIntPnts = ptelem->GetNumIntPoints();
    const Integer   dim = ptelem->GetDim();
    Matrix<Double>       dvShFnc;
    const std::vector<Double> & intWeights = ptelem->GetIntWeights(); 
    std::vector<double> shFnc;
    Matrix<Double>      ptCoord;
    Double              valsol;
    Double              jacDet;
    Double              diffGradL2norm,diffGradL2normIm;
    Double              SPRGradL2norm,SPRGradL2normIm; 


    error = 0;
    normGradSPR = 0;

    ptgrid->GetCoordNodesElemMat(connectVec,ptCoord,level); 

    //do a loop over int. points in element
    Integer iIntPnts;
    for (iIntPnts=0; iIntPnts<nrIntPnts; iIntPnts++)
      {   
     
        jacDet = ptelem->CalcJacobianDetAtIp(iIntPnts+1, ptCoord);
     
        // calculation of FEM gradient and SPR at intgr. points in local coord.
        Double     gradFEM[dim],gradFEMIm[dim];
        Double     gradSPR[dim],gradSPRIm[dim];
        Integer    idm,iShFnc;
      
        for (idm=0; idm<dim; idm++)
          {
            gradFEM[idm]=0;
            gradSPR[idm]=0;
            gradFEMIm[idm]=0;
            gradSPRIm[idm]=0;
          }
     
        ptelem->GetGlobDerivShFncAtIp(dvShFnc, iIntPnts + 1, ptCoord, jacDet); 
        ptelem->GetShFncAtIp(shFnc,iIntPnts+1); 
          
        for (iShFnc = 0; iShFnc < nrNodes; iShFnc++)
          for (idm = 0; idm < dim; idm++)
            {
              gradFEM[idm]+=sol[connectVec[iShFnc]-1]*dvShFnc[iShFnc][idm];
              gradFEMIm[idm]+=solIm[connectVec[iShFnc]-1]*dvShFnc[iShFnc][idm];
              gradSPR[idm]+=SPRgrad[idm][connectVec[iShFnc]-1]*shFnc[iShFnc];
              gradSPRIm[idm]+=SPRgradIm[idm][connectVec[iShFnc]-1]*shFnc[iShFnc];
            }
          
#ifdef DEBUG
#ifdef DEBUGFULL     
        (*debug) << gradSPR[0] << " " << gradFEM[0] << " y: " << gradSPR[1] << " " << gradFEM[1] << std::endl;
        if (dim==3) (*debug) << " z: " << gradSPR[2] << " " << gradFEM[2] << std::endl;
#endif
#endif

        // calculation of L2 norm for difference of SPR and FEM gradient,
        // and L2 norm of SPR gradient
        diffGradL2norm = 0;
        diffGradL2normIm = 0;
        SPRGradL2norm  = 0; 
        SPRGradL2normIm = 0;

        for (idm = 0; idm<dim; idm++)
          {
            diffGradL2norm+=sqr(gradSPR[idm]-gradFEM[idm]);
            diffGradL2normIm+=sqr(gradSPRIm[idm]-gradFEMIm[idm]);
            SPRGradL2norm+=sqr(gradSPR[idm]);
            SPRGradL2normIm+=sqr(gradSPRIm[idm]);   
          }
     
        diffGradL2norm=sqrt(diffGradL2norm+diffGradL2normIm);
        SPRGradL2norm=sqrt(SPRGradL2norm+SPRGradL2normIm);

        error+=diffGradL2norm*jacDet*intWeights[iIntPnts];
        normGradSPR+=SPRGradL2norm*jacDet*intWeights[iIntPnts];

      } // end of loop over integration points

#ifdef TRACE
    (*trace) << "leaving SpaceErrorEstimator::CalcErrorForElem" << std::endl;
#endif
  }

  // template <Integer Dim>
  // void SpaceErrorEstimator<Dim>::KellyError(Grid * ptgrid, std::vector<std::string> &sbdoms, Vector<Double> & errorMap,const Integer level) 
  // {
  // #ifdef TRACE
  //   (*trace) << "entering SpaceErrorEstimator::KellyError " << std::endl;
  // #endif
  // #ifdef DEBUG
  //   (*trace) << " step: kelly error map is calculated\n";
  // #endif
  //   Integer isd,jel;
  //   Elem * ptEl;

  //   // resize
  //   Integer maxelem=ptgrid->GetMaxnumElem(level,sbdoms);
  //   errorMap.Resize(maxelem);

  //   std::vector<Elem*> elemssd;
  //   Integer counter=0;
  //   for (isd=0; isd<sbdoms.size(); isd++){
 
  //     ptgrid->GetElemSD(elemssd,sbdoms[isd],level); 

  //     // loop over cells to calculate jump over faces
  //     // IntegralOverFace_KellyError();
  //     for (jel=0; jel < elemssd.size(); jel++, counter++)
  //       {
  //      std::vector<Elem*> faces;

  //      FormFacesForElem(elemssd[jel],faces);

  //      // loop over faces of Elem to calculate jump
  //      Integer noFcs=faces.size();
        
  //      //      std::cout << " size of faces: " << noFcs << std::endl;

  //      Double jump_sum=0;
  //      Double aux;
  //      Integer ifcs;
  //      for (ifcs=0; ifcs < noFcs; ifcs++) 
  //        {
  //          Elem * neighFcs[2];
  //          FormNeighForFace(neighFcs,ptgrid,faces[ifcs]);
  //          aux=IntegralOverRegularFace_KellyError(faces[ifcs],neighFcs);
  //          jump_sum+=aux;
  //        }

  //       errorMap[counter]=jump_sum;
  //       //      std::cout << " jump_ sum " << jump_sum << std::endl;
  //       } // loop over elems of subdomains
  //   } // loop over subdomains
  // } // end of KellyError

  // template <Integer Dim>
  // Double SpaceErrorEstimator<Dim>::IntegralOverRegularFace_KellyError(const Elem * ptface, Elem ** ptneigh)
  // {
  //   Double result=0;
  //   Integer ipoint,ifunc;
  //   Integer noQPoints, noQFunc; //number of integration points and shape functions 
  //   Integer noIntPnt=ptface->ptElem->GetNumIntPoints();
  //   Integer ip;
  //   for (ip=0; ip<noIntPnt; ip++) {
  //     ;
  //   }
  
  
  //   // 1. get gradients of the function at the first cell at the quadrature points. we know value of solution and shape functions in the quadrature points.
  //   // 2. the same for the second cell - neigbor
  //   // 3. calculate the difference. 
  //   // 4. take the normal
  //   // 5. now we have vector, each component of it is value at the quadrature point.
  //   // 6. take a square of each component of vector and sum up with integration weightes
  //   // 7. also here should be check is this face at the boundary or not
 
  //   // void calcNormal2Line(std::vector<Double> & normal,const Point2D a, const Point2D b) --> in tools.cc

  //   return result;
  // }

  // template<>
  // void SpaceErrorEstimator<2>::FormFacesForElem(const Elem * elem, std::vector<Elem*> & afaces)
  // {
  //   // 
  //   Elem * face;

  //   Integer type = (elem->ptElem)->feType();

  //   Integer k,l;
  //   switch(type) {
  //   case TRIA:
  //     // check that it's new
  //     // add to list
  //     //  afaces.resize(3);
  //     for (k=0; k<3; k++) {
  //       face=new Elem();
  //       face->ptElem=ptL1;
  //       face->connect.Resize(2);
  //       if (k<2)
  //       for (l=0; l<2; l++) 
  //      face->connect[l]=elem->connect[k+l];
  //       else {
  //      face->connect[0]=elem->connect[k];
  //      face->connect[1]=elem->connect[0];
  //       }
  //       face->namesd="face";      //// !!!!!
  //       //      afaces[k]=face;
  //       std::cout << face->connect[0] << " " << face->connect[1] << std::endl;
  //     }
  //     break;
  //   case QUAD:
  //     // check that it's new
  //     // add to list
  //     for (k=0; k<4; k++) {
  //       afaces.resize(4);
  //       face=new Elem();
  //       face->ptElem=ptL1;
  //       face->connect.Resize(2);
  //       for (l=0; l<2; l++) 
  //      face->connect[l]=elem->connect[k+l];
  //       afaces[k]=face;
  //     }
  //     break;
  //   default:
  //     Error("This type is not implemented",__FILE__,__LINE__);
  //     break;
  //   } // switch

  // } // end of FormFacesForElem

  // template<>
  // void SpaceErrorEstimator<3>::FormFacesForElem(const Elem * elem, std::vector<Elem*> & afaces)
  // {
  //   // 
  //   Elem * face;
  //   Integer type = elem->ptElem->feType();

  //       Integer k,l;
  //       switch(type) {
  //       case TET:
  //      // check that it's new
  //      // add to list
  //      for (k=0; k<3; k++) {
  //      face->ptElem=ptTr;
  //      face->connect.Resize(3);
  //      for (l=0; l<3; l++) 
  //        face->connect[l]=elem->connect[k+l];
  //      afaces.push_back(face);
  //      }
  //      break;
  //       case HEX:
  //      // check that it's new
  //      // add to list
  //      Error("Not implemented",__FILE__,__LINE__);
  //      for (k=0; k<4; k++) {
  //      face->ptElem=ptQ;
  //      face->connect.Resize(2);
  //      for (l=0; l<2; l++) 
  //        face->connect[l]=elem->connect[k+l];
  //      afaces.push_back(face);
  //      }
  //      break;
  //       default:
  //      Error("This type is not implemented",__FILE__,__LINE__);
  //      break;
  //       } // switch

  // } // end of FormFacesForElem

  // //
  // template<>
  // void SpaceErrorEstimator<2>::FormNeighForFace(Elem ** neighFcs, Grid * grid, Elem * face)
  // {
 
  //   Integer node=face->connect[0];
  //   std::vector<Elem*> * neighbours;
  //   neighbours=grid->GetNeighboursOfNode(node);

  //   Integer num=0;
  //   Integer i,j;
  //   for (i=0; i<neighbours->size(); i++) {
  //     Integer size=(*neighbours)[i]->connect.size();
  //     for (j=0; j < size; j++) {
  //       if (face->connect[1]==(*neighbours)[i]->connect[j]) {
  //      neighFcs[num]=(*neighbours)[i];
  //       num++;
  //     } // end if
      
  //     }
  //   } 
  
  // }

  // //
  // template<>
  // void SpaceErrorEstimator<3>::FormNeighForFace(Elem ** neighFcs, Grid * grid, Elem * face)
  // {
 
  //   Error("Not implemented",__FILE__,__LINE__);
  
  // }

} // end of namespace
