// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

extern "C" {

#define LE true

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/*********************************************************************
 *
 * endianSwap(addr,num_el,size_el): Changes between big- and 
 * little-endian by address-swapping.
 *
 * addr: (long int) address of variable to be swapped.
 * num_el: (int) number of elements in array pointed to by addr
 * size_el: (int) size of individual elements of *addr
 *
 * Written by Daniel J. Bodony (bodony@Stanford.EDU)
 * Copyright (c) 2001
 *
 * WARNING: This is only works on arrays that contiguous in memory.
 *
 *********************************************************************/

#define SWAP(a,b) temp=(a); (a)=(b); (b)=temp;

void endianSwap (long int addr, int num_el, int size_el) 
{
   int i, x;
   char *a;
   char temp;

   a = (char *)addr;

   for (i = 0; i < num_el; i++) {
      a = (char *)(addr + i*size_el);
      for (x = 0; x < size_el/2; x++) {
         SWAP(a[x],a[size_el-x-1]);
      }
   }
}


  void read_in(const char *filename, double** data, int nNodes, int nVars, int* nodesNr){

  int     ii, kk;
  int     nNodesIn;
  int     nVarsIn;
  double *ftmp;
  int    *ntmp;
  FILE   *in;
  int     record;

  if((in = fopen(filename, "rb")) == NULL){
    perror(filename);
  }

  // read header
	int return_value(0);
  return_value = fread(&record,sizeof(int),1,in);
  return_value = fread(&nNodesIn,sizeof(int),1,in);
  return_value = fread(&nVarsIn,sizeof(int),1,in);
  return_value = fread(&record,sizeof(int),1,in);

  // swap binary data from little endian to bigendian
#ifdef LE
  endianSwap((long int)&nNodesIn,1,sizeof(int));
  endianSwap((long int)&nVarsIn,1,sizeof(int));
#endif

  // check consistency of data
  if( nNodes != nNodesIn){
    printf("nNodes   = %d\n",nNodes);
    printf("nNodesIn = %d\n",nNodesIn);
    printf("===> nNodes != nNodesIn - check data file %s\n\n", filename);
    exit(1);
  }
  if( nVars != nVarsIn){
    printf("nVars   = %d\n",nVars);
    printf("nVarsIn = %d\n",nVarsIn);
    printf("===> nVars != nVarsIn - check data file %s\n\n", filename);
    exit(1);
  }
  
  ftmp = (double *)malloc(nNodes*sizeof(double));
  ntmp = (int *)malloc(nNodes*sizeof(int));
  
  // read node indec
  return_value = fread(&record,sizeof(int),1,in);
  return_value = fread(ntmp,sizeof(int),nNodes,in);
  return_value = fread(&record,sizeof(int),1,in);
#ifdef LE
  endianSwap((long int)ntmp,nNodes,sizeof(int));
#endif
  for(ii=0; ii<nNodes; ii++){
    nodesNr[ii] = ntmp[ii];
    //printf("node: %d\n",nodesNr[ii]);
  }

  //  printf("numVar: %d\n",nVars);

  // read data
  for(kk=0; kk<nVars; kk++){
    return_value = fread(&record,sizeof(int),1,in);
    return_value = fread(ftmp,sizeof(double),nNodes,in);
#ifdef LE
      endianSwap((long int)ftmp,nNodes,sizeof(double));
#endif
    return_value = fread(&record,sizeof(int),1,in);

    if(kk == 0){
       //speed of sound
      for(ii=0; ii<nNodes; ii++){
        data[kk][ii] = ftmp[ii];
   //     printf("sos: %f\n",data[kk][ii]);
      }
    }
    else if(kk == 1){
      // reynold stress
      for(ii=0; ii<nNodes; ii++){
        data[kk][ii] = ftmp[ii];
    //    printf("stress: %f\n",data[kk][ii]);
      }
    }
    else if(kk == 2){
      // momentum flux fluctuation
      for(ii=0; ii<nNodes; ii++){
        data[kk][ii] = ftmp[ii];
     //   printf("momentum: %f\n",data[kk][ii]);
      }
    }
    else if(kk == 3){
       // chemcial
      for(ii=0; ii<nNodes; ii++){
        data[kk][ii] = ftmp[ii];
      //  printf("chemical: %f\n",data[kk][ii]);
      }
    }
    else if(kk == 4){
      // second time derivative of density
      for(ii=0; ii<nNodes; ii++){
        data[kk][ii] = ftmp[ii];
       // printf("density: %f\n",data[kk][ii]);
      }
    }
    else if(kk == 5){
      // heat release
      for(ii=0; ii<nNodes; ii++){
        data[kk][ii] = ftmp[ii];
       // printf("heat: %f\n",data[kk][ii]);
      }
    }
    else if(kk == 6){
      // second time derivative of gas constant R
      for(ii=0; ii<nNodes; ii++){
        data[kk][ii] = ftmp[ii];
        //printf("gas: %f\n",data[kk][ii]);
      }
    }
    else if(kk == 7){
     // shear term 
      for(ii=0; ii<nNodes; ii++){
        data[kk][ii] = ftmp[ii];
       // printf("shear: %f\n",data[kk][ii]);
      }
    }
    else{
      printf("nVar larger than expected - check nVar in SR read_in\n");
      exit(1);
    }
  }
  free(ftmp);
  free(ntmp);
  fclose(in);
  return;
}

}
