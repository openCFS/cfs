
#ifndef FILE_SIMBEM_PILES
#define FILE_SIMBEM_PILES

namespace CoupledField
{
//   void GetConnectivity(Integer connect[], Integer size);
  
//   void GetCoordinate(Double coord[], Double meshsize, Integer size);

//   void MakeDirichlet(Double dir[], Double meshsize, Integer size);

//   void MakeGnuFile(Double u[], Integer size);

//   void Print(Double mat[], Integer size);

  /// OSTBEM

  void MakeElement(Integer numnode, Integer connect[], Double coord[], Double des[]);

  void MakeSingleLaplace(Integer numelem, Integer numnode, Integer connect[],Double coord[], 
			 Integer ga, Integer nu, Double mat[], Double epar[]);

  void MakeDoubleLaplace(Integer numelem1, Integer numnode1, Integer connect1[], Double coord1[],
			 Integer numelem2, Integer numnode2, Integer connect2[], Double coord2[],
			 Integer ga, Integer nu, Integer mu, Double mat[], Double epar[]);

  void MakeHyperLaplace(Integer numelem, Integer numnode, Integer type, Integer nu, 
			Integer connect[], Double mat1[], Double mat2[], Double des[], 
			Double epar[]);

  void MakeMasse(Integer numelem, Integer nu, Integer ga, Integer connect1[], Double coord1[],
		 Double des[], Integer connect2[], Double coord2[], Double mat[], 
		 Integer numnode);

  void CalcData(Integer numelem1, Integer numnode1, Integer connect1[], Double coord1[], 
		Integer numelem2, Integer numnode2, Integer connect2[], Double coord2[],
		Integer type, Integer nu, Integer mu, Double dir[], Double neu[], Double mat[], 
		Double point[], Double val[]);
}

#endif // FILE_SIMBEM_PILES
