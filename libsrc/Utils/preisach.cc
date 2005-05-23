#include <iostream>
#include <fstream>

#include "preisach.hh"

namespace CoupledField
{ 

  Preisach :: Preisach(Integer numElem, Double xSat, Double ySat, Double xRem,
                       Boolean isVirgin) 
    : Hysteresis(numElem)
  {
    ENTER_FCN("Preisach::Preisach" );

    if (xSat > 0 ) {
      Xsaturated_  = xSat;
    }
    else {
      Xsaturated_  = 1.0;
    }

    YSaturated_  = ySat;
    YRemnant_    = xRem;
    isVirgin_    = isVirgin;

    lastVal_.Resize(numElem);
    preisachSum_.Resize(numElem);
    preisachSum_.Init(0);

    isMinMax_    = new std::vector<Integer>[numElem];
    extremaList_ = new std::vector<Double>[numElem];
  }

  Preisach :: ~Preisach()
  {
    delete [] isMinMax_;
    delete [] extremaList_;
  }

  Double Preisach :: computeValue(Double Xin, Integer idxElem) 
  {
    ENTER_FCN( "Preisach::computeValue" );

    Integer idx = idxElem - 1;
  
    //normalize input
    Double newX = normalizeInput(Xin);

    // get last entry of extremaList
    Integer posEL = extremaList_[idx].size()-1;
    Double lastX  = extremaList_[idx][posEL] ;

    Double val =  everett(lastX,newX);

    //  std::cerr << "lastX=" << lastX << " newX=" << newX << " val=" << val << std::endl;
    return ( (val + preisachSum_[idx])*YSaturated_ );
  }


  void Preisach :: updateMinMaxList(Double Xin, Integer nrEl)
  {
    ENTER_FCN( "Preisach::updateMinMaxList" );

    //==========================================================================//
    //
    //  modification of the extrema list                                         
    //
    //  empty list:                                                            
    //
    //      I) isvirgin = 0                                                    
    //         insert two extrema for polarization into the list and    
    //         insert the new timestepvalue as can be seen below.   
    //
    //     II) isvirgin = 1                                                    
    //         insert the new timestepvalue in the extremalist as       
    //         maximum if value is greater zero and as minimum if 
    //         value is less than zero                       
    //
    //  general algorithm:
    //
    //  if the new value is larger as the last value in the extrema list and
    //  this last value was a maximum, then overwrite the last value in the 
    //  extrema list with the new one
    //
    //  if the new value is smaller as the last value in the extrema list and
    //  this last value was a minimum, then overwrite the last value in the 
    //  extrema list with the new one
    //
    //  if the new value is larger as the last value in the extrema list and
    //  this last value was a minimum, then insert the new value as a maximum
    //
    //  if the new value is smaller as the last value in the extrema list and
    //  this last value was a maximum, then insert the new value as a minimum
    //
    //==========================================================================//

    //normalize input
    Double newX = normalizeInput(Xin);

    //array starts at 0!!
    nrEl -= 1;

    if ((isMinMax_[nrEl].size() == 0) && (isVirgin_ == TRUE) ) {
      //list of minimums and maximums is empty
      extremaList_[nrEl].push_back(newX);
      if (newX > 0) {
        isMinMax_[nrEl].push_back(1);
      }
      else if (newX < 0) {
        isMinMax_[nrEl].push_back(-1);
      }
      else {
        //exact zero: complicated to handle;
        isMinMax_[nrEl].push_back(0);
      }

      preisachSum_[nrEl] = 0.5 * everett(-newX, newX);
    }

    else {
      if ( isMinMax_[nrEl].size() == 0 ) {
        //insert polarisation state
        isMinMax_[nrEl].push_back(2);
        //   extremaList_[nrEl].push_back(Xsaturated_);
        extremaList_[nrEl].push_back(1);
        preisachSum_[nrEl] = 0.5*everett(-1,1);
        // preisachSum_[nrEl] = everett(-Xsaturated_,Xsaturated_);

        isMinMax_[nrEl].push_back(-2);
        extremaList_[nrEl].push_back(0);
        preisachSum_[nrEl] += everett(1,0);
        //      extremaList_[nrEl].push_back(XCoercitive_);
        //    preisachSum_[nrEl] += everett(XCoercitive_,Xsaturated_);
      }

      // get last entry of extremaList
      Integer posEL = extremaList_[nrEl].size()-1;
      Double lastX  = extremaList_[nrEl][posEL] ;

      // get last entry of minMaxList
      Integer posMM      = isMinMax_[nrEl].size()-1;
      Integer lastMinMax = isMinMax_[nrEl][posMM];

      Double prevLastX=0;
      if ( posMM > 0)
        prevLastX = extremaList_[nrEl][posMM-1];

      //if the first value was excat zero, lasMinMax was set to 0
      //now check if we should set it to 1 or to -1;
      if (lastMinMax == 0 ) {
        if (newX >  prevLastX) {
          lastMinMax = 1;
          isMinMax_[nrEl][posMM] = 1;
        }
        else {
          lastMinMax = -1;
          isMinMax_[nrEl][posMM] = -1;
        }
      }

      if ( lastMinMax > 0 ) {
        //last entry was a maximum
        if ( newX >= lastX ) {
          //overwrite last entry
          extremaList_[nrEl][posEL] = newX;
          if ( posMM == 0) {
            preisachSum_[nrEl] = 0.5 * everett(-newX, newX);
          }
          else {
            preisachSum_[nrEl] -= everett(prevLastX,lastX);
            preisachSum_[nrEl] += everett(prevLastX,newX);
          }

        }
        else {
          //add new value as a minimum
          isMinMax_[nrEl].push_back(-1);
          extremaList_[nrEl].push_back(newX);
          preisachSum_[nrEl] += everett(lastX,newX);
        }
      }
      else {
        //last entry was a minimum
        if ( newX <= lastX ) {
          //overwrite last entry
          extremaList_[nrEl][posEL] = newX;
          if ( posMM == 0) {
            preisachSum_[nrEl] = 0.5 * everett(-newX, newX);
          }
          else {
            preisachSum_[nrEl] -= everett(prevLastX,lastX);
            preisachSum_[nrEl] += everett(prevLastX,newX);
          }
        }
        else {
          //add new value as a maximum
          isMinMax_[nrEl].push_back(1);
          extremaList_[nrEl].push_back(newX);
          preisachSum_[nrEl] += everett(lastX,newX);
        }
      }
    }

    //  Integer size = extremaList_[nrEl].size();

    //   std::cout << "Size=" << extremaList_[nrEl].size() << std::endl;
    //   std::vector<Double>   &extremaList =  extremaList_[nrEl];
  
    //   std::cout << "List" << std::endl;
    //   for (Integer k=0; k < extremaList.size(); k++) {
    //     std::cout << " " << extremaList[k] << std::endl;
    //   }

    //now check, if we can cancle any extrema
    //   if (isMinMax_[nrEl].size() > 0) {
    //     wipout(nrEl);
    //   }

  }


  void Preisach :: wipout(Integer nrEl)
  {
    ENTER_FCN( "Preisach::wipout" );

    //==========================================================================//
    //  clear all non dominant maxima out of the list                            
    //     
    //  I) last (latest) value is a maxima: will be denoted as the actual 
    //                                      maximum
    //
    //         if the actual maximum is larger than the last maximum    
    //         search for a minimum that is smaller than the actual  
    //         minimum (= value before the actual maximum);
    //         test if all maxima in between are smaller than the actual 
    //         one. If such a minimum exists clear all extrema between 
    //         the actual maximum and this minimum.
    //
    //
    //  II) last (latest) value is a minimum: will be denotes as the actual 
    //                                        minimum
    //
    //         if the actual minimum is smaller than the last minimum    
    //         search for a maximum that is larger than the actual  
    //         maximum (= value before the actual minimum);
    //         test if all maxima in between are larger than the actual 
    //         one. If such a maximim exists clear all extrema between 
    //         the actual minimum and this maximum.
    //                                            
    //==========================================================================//
                                                                           

    std::vector<Double>   &extremaList =  extremaList_[nrEl];
    std::vector<Integer>  &isMinMax    = isMinMax_[nrEl];

    // get last entry of extremaList
    //  Double lastX  = extremaList.back();
  
    // get last entry of minMaxList
    Integer lastMinMax = isMinMax.back();
  
    Integer actPosMax, actPosMin, lastPos, actPos;
    std::vector<Double>::iterator iterExtremalList;
    std::vector<Integer>::iterator iterMinMax;

    std::vector<Double>::iterator endDeleteExtremList, startDeleteExtremList;
    std::vector<Integer>::iterator endDeleteMinMax, startDeleteMinMax ;

    std::cout << "ExtremaList before wipout: " << extremaList.size() << std::endl;

    Boolean kill = TRUE;
    if ( lastMinMax == 1 ) {
      //last value is a maximum
  
      Boolean minFound = FALSE;

      while ( kill ) {
        std::cout << "kill=" << kill << std::endl;

        //position of actual maximum
        actPosMax = extremaList.size()-1;

        //position of actual minmum
        actPosMin = actPosMax -1;

        //set the iterators
        iterExtremalList = extremaList.end();
        iterMinMax = isMinMax.end();

        if (actPosMin >= 0 && isMinMax[actPosMin] > -2 ) {
          lastPos = actPosMin;
          actPos = lastPos - 1;
          --iterExtremalList;
          --iterMinMax;
          std::cout << "actPosMin=" << actPosMin << std::endl;

          //now the iteration counter points to the last element
          //in Extrema- and MinmaxList, which is a maximum!!
          //this element is not allowed to be cleared!
          endDeleteExtremList = iterExtremalList;
          endDeleteMinMax     = iterMinMax;

          while ( kill == TRUE && minFound == FALSE ) {
            std::cout << "actPos=" << actPos << std::endl;
            if ( actPos >= 0) {
              if ( isMinMax[actPos] == 2 || isMinMax[actPos] == -2 ) {
                kill = FALSE;
              }
              else if ( isMinMax[actPos] == 1 && extremaList[actPos] > extremaList[actPosMax]) {
                kill = FALSE;
              }
              else if ( isMinMax[actPos] == -1 && extremaList[actPos] < extremaList[actPosMin]) {
                minFound = TRUE;
              }

              lastPos = actPos;
              actPos -= 1;
              --iterExtremalList;
              --iterMinMax;

              //if minFound = TRUE, then iterExtremalList (iterMinMax) points to
              //the first element, which will be wiped out
              startDeleteExtremList = iterExtremalList;
              startDeleteMinMax     = iterMinMax;
            }
            else {
              kill = FALSE;
            }
          }

          if ( (kill = TRUE) && (minFound == TRUE) ) {
            // kill all extrema in between
            extremaList.erase(startDeleteExtremList,endDeleteExtremList);
            isMinMax.erase(startDeleteMinMax,endDeleteMinMax);
            minFound = FALSE;
          }
        }

        else {
          kill = FALSE;
        }
   
      }

    }
  
    else if ( lastMinMax == -1 ) {
      //last value is a minimum
      
      Boolean maxFound = FALSE;
      while ( kill ) {
        std::cout << "kill=" << kill << std::endl;

        //position of actual minimum
        actPosMin = extremaList.size()-1;

        //position of actual maximum
        actPosMax = actPosMin -1;

        //set the iterators
        iterExtremalList = extremaList.end();
        iterMinMax = isMinMax.end();

        std::cout << "actPosMax=" << actPosMax << std::endl;

        if (actPosMax >= 0 && isMinMax[actPosMax] < 2 ) {
          lastPos = actPosMax;
          actPos  = lastPos - 1;
          --iterExtremalList;
          --iterMinMax;

          std::cout << "lastPos=" << lastPos << std::endl;

          //now the iteration counter point to the last element
          //in Extrema- and MinmaxList, which is a maximum!!
          //this element is not allowed to be cleared!
          endDeleteExtremList = iterExtremalList;
          endDeleteMinMax     = iterMinMax;

          while ( kill == TRUE && maxFound == FALSE ) {
            std::cout << "In while, actPos=" << actPos << std::endl;

            if ( actPos >= 0) {
              if ( isMinMax[actPos] == 2 || isMinMax[actPos] == -2 ) {
                kill = FALSE;
              }
              else if ( isMinMax[actPos] == -1 && extremaList[actPos] < extremaList[actPosMin]) {
                kill = FALSE;
              }
              else if ( isMinMax[actPos] == 1 && extremaList[actPos] > extremaList[actPosMax]) {
                maxFound = TRUE;
              }
              lastPos = actPos;
              actPos -= 1;
              --iterExtremalList;
              --iterMinMax;

              //if maxFound = TRUE, then iterExtremalList (iterMinMax) points to
              //the first element, which will be wiped out
              startDeleteExtremList = iterExtremalList;
              startDeleteMinMax     = iterMinMax;
            }
            else {
              kill = FALSE;
            }
          }

          if ( (kill = TRUE) && (maxFound == TRUE) ) {
            // kill all extrema in between
            extremaList.erase(startDeleteExtremList,endDeleteExtremList);
            isMinMax.erase(startDeleteMinMax,endDeleteMinMax);
            maxFound = FALSE;
          }
        }

        else {
          kill = FALSE;
        }
   
      }
    }

    std::cout << "ExtremaList after wipout: " << extremaList.size() << std::endl;

  }


  Double Preisach :: everett(Double X1, Double X2)
  {
    ENTER_FCN( "Preisach:: everett" );

    Double newY;
    Double diffX = X2 - X1;

    if ( diffX > 0) {
      newY = 0.5*diffX*diffX;
    }
    else {
      newY = -0.5*diffX*diffX;
    }

    return newY;
  }


  Double Preisach :: normalizeInput(Double Xin)
  {
    ENTER_FCN( "Preisach::normalizeInput" );

    Double Xout;

    if ( Xin > Xsaturated_ ) {
      //saturation achieved!!
      Xout = 1.0;
    }
    else if ( Xin < -Xsaturated_ ) {
      //saturation achieved!!
      Xout = -1.0;
    }
    else {
      //normalize input
      Xout = Xin / Xsaturated_;
    }

    return Xout;
  }

}
