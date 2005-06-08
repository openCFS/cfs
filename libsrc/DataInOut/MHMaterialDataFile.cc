#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <limits.h>
#include <string>

#include "MHMaterialDataFile.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"


namespace CoupledField
{

  MHMaterialDataFile::MHMaterialDataFile()
  {
    ENTER_FCN("MHMaterialDataFile::MHMaterialDataFile");
  }

  MHMaterialDataFile::~MHMaterialDataFile()
  {
    ENTER_FCN("MHMaterialDataFile::~MHMaterialDataFile");
    indexSetFile->close();
    indexSetFileIn->close();
    indexSetCountFile->close();
    indexSetCountFileIn->close();
  }

  void MHMaterialDataFile::createFiles()  {
    ENTER_FCN( "MHMaterialData::createFiles" );
  }
  
  void MHMaterialDataFile::computeIndexSet(Integer PP, UInt PPmax, UInt nrMultHarms_){
    ENTER_FCN("MultiHarmonicDriver::computeIndexSet");

    Integer countmax;
    Integer inddiff;
    Integer p;
    Integer delta;

    Matrix<Integer> count;

    count.Resize(PPmax+1, 4*nrMultHarms_+2);

    countmax = pow(2*nrMultHarms_+1,PPmax);

    std::string filename = "indexSetN";
    std::string filenameCount = "indexSetCountN";

   
    //    Integer indexset[PPmax+2][4*nrMultHarms_+2][countmax+1][2*nrMultHarms_+2];
    Integer indexset[PPmax+1][4*nrMultHarms_+2][countmax][2*nrMultHarms_+2];
    // Integer indexset[PPmax][4*nrMultHarms_][2*nrMultHarms_];
    for (UInt i=0;i<PPmax;i++)
      for (UInt j=0;j<4*nrMultHarms_+2;j++)
        for (UInt k=0;k<2*nrMultHarms_;k++)
          for (UInt l=0;l<PPmax;l++)
            indexset[i][j][k][l]=0;

    //Integer ****indexset;

    switch (nrMultHarms_){
    case 1:
      filename=filename+"1_";
      filenameCount=filenameCount+"1_";
      for (Integer pm1=0; pm1<=PP; pm1++)
        for (Integer p0=0; p0<=PP; p0++)
          for (Integer pp1=0; pp1<=PP; pp1++){
            p=pm1+p0+pp1;
            delta = -pm1+pp1;
            if(((p<=PP)&&(-2*(Integer)nrMultHarms_<=(Integer)delta))&&((Integer)delta<=2*(Integer)nrMultHarms_))
              {         
              inddiff=delta+2*nrMultHarms_+1;
              count[p][inddiff]=count[p][inddiff]+1;
              indexset[p][inddiff][count[p][inddiff]][0]=pm1;
              indexset[p][inddiff][count[p][inddiff]][1]=p0;
              indexset[p][inddiff][count[p][inddiff]][2]=pp1;
              indexset[p][inddiff][count[p][inddiff]][3]=delta+2*nrMultHarms_;
              }
          }
      break;
    case 2:
      filename=filename+"2_";
      filenameCount=filenameCount+"2_";
      for (Integer pm2=0; pm2<=PP; pm2++)
        for (Integer pm1=0; pm1<=PP; pm1++)
          for (Integer p0=0; p0<=PP; p0++)
            for (Integer pp1=0; pp1<=PP; pp1++)
              for (Integer pp2=0; pp2<=PP; pp2++){
                p=pm2+pm1+p0+pp1+pp2;
                delta = -2*pm2-pm1+pp1+2*pp2;
                if(((p<=PP)&&(-2*(Integer)nrMultHarms_<=(Integer)delta))&&((Integer)delta<=2*(Integer)nrMultHarms_))
                  {            
                  inddiff=delta+2*nrMultHarms_+1;
                  count[p][inddiff]=count[p][inddiff]+1;
                  indexset[p][inddiff][count[p][inddiff]][0]=pm2;
                  indexset[p][inddiff][count[p][inddiff]][1]=pm1;
                  indexset[p][inddiff][count[p][inddiff]][2]=p0;
                  indexset[p][inddiff][count[p][inddiff]][3]=pp1;
                  indexset[p][inddiff][count[p][inddiff]][4]=pp2;
                  indexset[p][inddiff][count[p][inddiff]][5]=delta+2*nrMultHarms_;
                }
              }
      break;
    case 3: 
      filename=filename+"3_";
      filenameCount=filenameCount+"3_";
      for (Integer pm3=0; pm3<=PP; pm3++)
        for (Integer pm2=0; pm2<=PP; pm2++)
          for (Integer pm1=0; pm1<=PP; pm1++)
            for (Integer p0=0; p0<=PP; p0++)
              for (Integer pp1=0; pp1<=PP; pp1++)
                for (Integer pp2=0; pp2<=PP; pp2++)
                  for (Integer pp3=0; pp3<=PP; pp3++){
                    p=pm3+pm2+pm1+p0+pp1+pp2+pp3;
                    delta = -3*pm3-2*pm2-pm1+pp1+2*pp2+3*pp3;
               
                    if(((p<=PP)&&(-2*(Integer)nrMultHarms_<=(Integer)delta))&&((Integer)delta<=2*(Integer)nrMultHarms_)){            
                      inddiff=delta+2*nrMultHarms_+1;
                      count[p][inddiff]=count[p][inddiff]+1;
                      indexset[p][inddiff][count[p][inddiff]][0]=pm3;
                      indexset[p][inddiff][count[p][inddiff]][1]=pm2;
                      indexset[p][inddiff][count[p][inddiff]][2]=pm1;
                      indexset[p][inddiff][count[p][inddiff]][3]=p0;
                      indexset[p][inddiff][count[p][inddiff]][4]=pp1;
                      indexset[p][inddiff][count[p][inddiff]][5]=pp2;
                      indexset[p][inddiff][count[p][inddiff]][6]=pp3;
                      indexset[p][inddiff][count[p][inddiff]][7]=delta+2*nrMultHarms_;
                    }
                  }           
      break;

    case 4: 
      filename=filename+"4_";
      filenameCount=filenameCount+"4_";
      for (Integer pm4=0; pm4<=PP; pm4++)
        for (Integer pm3=0; pm3<=PP; pm3++)
          for (Integer pm2=0; pm2<=PP; pm2++)
            for (Integer pm1=0; pm1<=PP; pm1++)
              for (Integer p0=0; p0<=PP; p0++)
                for (Integer pp1=0; pp1<=PP; pp1++)
                  for (Integer pp2=0; pp2<=PP; pp2++)
                    for (Integer pp3=0; pp3<=PP; pp3++)
                      for (Integer pp4=0; pp4<=PP; pp4++){
                        p=pm4+pm3+pm2+pm1+p0+pp1+pp2+pp3+pp4;
                        delta =-4*pp4 -3*pm3-2*pm2-pm1+pp1+2*pp2+3*pp3+4*pp4;                        
                        if(((p<=PP)&&(-2*(Integer)nrMultHarms_<=(Integer)delta))&&((Integer)delta<=2*(Integer)nrMultHarms_)){            
                          inddiff=delta+2*nrMultHarms_+1;                          
                          count[p][inddiff]=count[p][inddiff]+1;
                          indexset[p][inddiff][count[p][inddiff]][0]=pm4;
                          indexset[p][inddiff][count[p][inddiff]][1]=pm3;
                          indexset[p][inddiff][count[p][inddiff]][2]=pm2;
                          indexset[p][inddiff][count[p][inddiff]][3]=pm1;
                          indexset[p][inddiff][count[p][inddiff]][4]=p0;
                          indexset[p][inddiff][count[p][inddiff]][5]=pp1;
                          indexset[p][inddiff][count[p][inddiff]][6]=pp2;
                          indexset[p][inddiff][count[p][inddiff]][7]=pp3;
                          indexset[p][inddiff][count[p][inddiff]][8]=pp4;
                          indexset[p][inddiff][count[p][inddiff]][9]=delta+2*nrMultHarms_;
                        }
                      }           
      break;

    case 5: 
      filename=filename+"5_";
      filenameCount=filenameCount+"5_";
      for (Integer pm5=0; pm5<=PP; pm5++)
        for (Integer pm4=0; pm4<=PP; pm4++)
          for (Integer pm3=0; pm3<=PP; pm3++)
            for (Integer pm2=0; pm2<=PP; pm2++)
              for (Integer pm1=0; pm1<=PP; pm1++)
                for (Integer p0=0; p0<=PP; p0++)
                  for (Integer pp1=0; pp1<=PP; pp1++)
                    for (Integer pp2=0; pp2<=PP; pp2++)
                      for (Integer pp3=0; pp3<=PP; pp3++)
                        for (Integer pp4=0; pp4<=PP; pp4++)
                          for (Integer pp5=0; pp5<=PP; pp5++){
                            p=pm5+pm4+pm3+pm2+pm1+p0+pp1+pp2+pp3+pp4+pp5;
                            delta =-5*pp5-4*pp4 -3*pm3-2*pm2-pm1+pp1+2*pp2+3*pp3+4*pp4+pp5;                        
                            if(((p<=PP)&&(-2*(Integer)nrMultHarms_<=(Integer)delta))&&((Integer)delta<=2*(Integer)nrMultHarms_))
                              {            
                                inddiff=delta+2*nrMultHarms_+1;                          
                                count[p][inddiff]=count[p][inddiff]+1;
                                indexset[p][inddiff][count[p][inddiff]][0]=pm5;
                                indexset[p][inddiff][count[p][inddiff]][1]=pm4;
                                indexset[p][inddiff][count[p][inddiff]][2]=pm3;
                                indexset[p][inddiff][count[p][inddiff]][3]=pm2;
                                indexset[p][inddiff][count[p][inddiff]][4]=pm1;
                                indexset[p][inddiff][count[p][inddiff]][5]=p0;
                                indexset[p][inddiff][count[p][inddiff]][6]=pp1;
                                indexset[p][inddiff][count[p][inddiff]][7]=pp2;
                                indexset[p][inddiff][count[p][inddiff]][8]=pp3;
                                indexset[p][inddiff][count[p][inddiff]][9]=pp4;
                                indexset[p][inddiff][count[p][inddiff]][10]=pp5;
                                indexset[p][inddiff][count[p][inddiff]][11]=delta+2*nrMultHarms_;
                              }
                          }           
      break;
      
    default:
      std::cout<<"Your order N of harmonics considered is too high (for this implementation)"<<std::endl;
      break;
    }

    
    if(PP==1)
      {
        filename+="PP1.dat";
        filenameCount+="PP1.dat";
      }
    else if (PP==2)
      {
        filename+="PP2.dat";
        filenameCount+="PP2.dat";
      }
    else if (PP==3)
      {
        filename+="PP3.dat";
        filenameCount+="PP3.dat";
      }
    else if (PP==4)
      {
        filename+="PP4.dat";
        filenameCount+="PP4.dat";
      }
    else if (PP==5)
      {
        filename+="PP5.dat";
        filenameCount+="PP5.dat";
      }

      indexSetFile = new std::ofstream(filename.c_str(),std::basic_ios<char>::out);
      if (!indexSetFile){
        std::cerr << "\n indexSetFile could not be initialized" << std::endl;
      }
      indexSetCountFile = new std::ofstream(filenameCount.c_str(),std::basic_ios<char>::out);
      if (!indexSetCountFile){
        std::cerr << "\n indexSetCountFile could not be initialized" << std::endl;
      }

      for (Integer pInd=0;pInd<=PP;pInd++)
        {
          
          *indexSetFile<<"p="<<pInd<<std::endl;
          
              for(UInt jj=0;jj<4*nrMultHarms_+2;jj++){
                for(Integer countInd=1; countInd<count[pInd][jj]; countInd++){
                  *indexSetFile<<indexset[pInd][jj][count[pInd][jj]][2*nrMultHarms_+1]<<")"; // writes the delta+2N+1
                  for(UInt kk=0;kk<2*nrMultHarms_+1;kk++){
                    *indexSetFile<<indexset[pInd][jj][countInd][kk]<<" ";
                    if(kk==2*nrMultHarms_)
                      *indexSetFile<<","<<std::endl;
                  }
                }
              }
        }
      *indexSetFile<<"-"<<std::endl;
      
      for (UInt pInd=0;pInd<=(UInt)PP;pInd++)
        for(UInt jj=1;jj<4*nrMultHarms_+2;jj++){
          if (count[pInd][jj]>0)
            *indexSetCountFile<<count[pInd][jj]-1<<" ";
          else 
            *indexSetCountFile<<"0 ";
          if(jj==4*nrMultHarms_+1)
            *indexSetCountFile<<std::endl;
        }
      
    
  } // end compute IndexSet

  void MHMaterialDataFile::getExponentArray(Matrix<Integer> & exponent, Matrix<Integer> & count, UInt PP, UInt p, Integer delta)
  {
    ENTER_FCN("MultiHarmonicDriver::getExponentArray");
    
    std::string filename = "indexSetN";
    std::string filenameCount = "indexSetCountN";
    std::string fileExt;
    std::string fileExtCount;
    Integer nrMultHarms_;
    StdVector<std::string> keyVec, attrVec, valVec;
    attrVec = "tag";
    valVec = "MULTIHARMONIC";
    keyVec = "multiHarmonic", "nrMultiHarmonics";
    params->Get(keyVec, attrVec, valVec, nrMultHarms_);

    char mDataRow[256];
    
    if(PP==1)
      {
        fileExt="PP1.dat";
      }
    else if (PP==2)
      {
        fileExt="PP2.dat";
      }
    else if (PP==3)
      {
        fileExt="PP3.dat";
      }
    else if (PP==4)
      {
        fileExt="PP4.dat";
      }
    else if (PP==5)
      {
        fileExt="PP5.dat";
      }
    
    switch(nrMultHarms_)
      {
      case 1:

        filename=filename+"1_"+fileExt;
        filenameCount=filenameCount+"1_"+fileExt;
        indexSetFileIn = new std::ifstream(filename.c_str(),std::basic_ios<char>::out);
        if (*indexSetFileIn==NULL){
          std::cerr << "\n " <<filename<<  " does not exist! " << std::endl;
          getchar();
        }
        indexSetCountFileIn = new std::ifstream(filenameCount.c_str(),std::basic_ios<char>::out);
        if (*indexSetCountFileIn==NULL){
          std::cerr << "\n " <<filenameCount<<  " does not exist! " << std::endl;
          getchar();
        }
      break;

    
      case 2: 
        filename=filename+"2_"+fileExt;
        filenameCount=filenameCount+"2_"+fileExt;
        indexSetFileIn = new std::ifstream(filename.c_str(),std::basic_ios<char>::in);
        if (*indexSetFileIn==NULL){
          std::cerr << "\n " <<filename<<  " does not exist! " << std::endl;
          getchar();
        }
        indexSetCountFileIn = new std::ifstream(filenameCount.c_str(),std::basic_ios<char>::in);
        if (*indexSetCountFileIn==NULL){
          std::cerr << "\n " <<filenameCount<<  " does not exist! " << std::endl;
          getchar();
        }
        break;

      case 3:
        filename=filename+"3_"+fileExt;
        filenameCount=filenameCount+"3_"+fileExt;
        indexSetFileIn = new std::ifstream(filename.c_str(),std::basic_ios<char>::in);
        if (*indexSetFileIn==NULL){
          std::cerr << "\n " <<filename<<  " does not exist! " << std::endl;
          getchar();
        }
        indexSetCountFileIn = new std::ifstream(filenameCount.c_str(),std::basic_ios<char>::in);
        if (*indexSetCountFileIn==NULL){
          std::cerr << "\n " <<filenameCount<<  " does not exist! " << std::endl;
          getchar();
        }

        break;

      case 4:
        filename=filename+"4_"+fileExt;
        filenameCount=filenameCount+"4_"+fileExt;
        indexSetFileIn = new std::ifstream(filename.c_str(),std::basic_ios<char>::in);
        if (*indexSetFileIn==NULL){
          std::cerr << "\n " <<filename<<  " does not exist! " << std::endl;
          getchar();
        }
        indexSetCountFileIn = new std::ifstream(filenameCount.c_str(),std::basic_ios<char>::in);
        if (*indexSetCountFileIn==NULL){
          std::cerr << "\n " <<filenameCount<<  " does not exist! " << std::endl;
          getchar();
        }
        break;

      case 5:
        filename=filename+"5_"+fileExt;
        filenameCount=filenameCount+"5_"+fileExt;
        indexSetFileIn = new std::ifstream(filename.c_str(),std::basic_ios<char>::in);
        if (indexSetFileIn==NULL){
          std::cerr << "\n " <<filename<<  " does not exist! " << std::endl;
          getchar();
        }
        indexSetCountFileIn = new std::ifstream(filenameCount.c_str(),std::basic_ios<char>::in);
        if (*indexSetCountFileIn==NULL){
          std::cerr << "\n " <<filenameCount<<  " does not exist! " << std::endl;
          getchar();
        }
        break;

      default:
        std::cout<<"Your order N of harmonics considered is too high (for this implementation)"<<std::endl;
        break;
    } // end switch

    
    
    char cChar[2];
    count.Resize(PP+1,4*nrMultHarms_+2);
    for (UInt pInd=0;pInd<=PP;pInd++){
      indexSetCountFileIn->getline(mDataRow,265);
      for(Integer jj=0;jj<4*nrMultHarms_+2;jj++){
        cChar[0]=mDataRow[2*jj];
        cChar[1]=mDataRow[2*jj+1];
        count[pInd][jj]=atoi(cChar);
      }
    }

    while(indexSetFileIn->getline(mDataRow, 265)){
      if (mDataRow[0]=='p'){
        char pChar[1];
        pChar[0]=mDataRow[2];
        UInt pNum=atoi(pChar);
        if(pNum==p){
          Integer cCounter=0;
          exponent.Resize(count[p][delta+2*nrMultHarms_],2*nrMultHarms_+1);
          do{
            indexSetFileIn->getline(mDataRow,265);
            char pChar2[2];
            pChar2[0]=mDataRow[0];
            pChar2[1]=mDataRow[1];
            Integer deltaNum = atoi (pChar2);
            if (deltaNum==delta+2*nrMultHarms_){
              UInt counter=0;
              UInt expCounter=0;
              while(mDataRow[2+counter]!=','){
                if ((mDataRow[2+counter]!=',')&&(mDataRow[2+counter]!=' ')){
                  exponent[cCounter][expCounter] = (mDataRow[2+counter])-48;
                  expCounter++;
                }
                counter++;
              }
              cCounter++;
            }                
          }
          while(mDataRow[0]!='-');
        }
      }
    }
    
  }
  

} // end of namespace
