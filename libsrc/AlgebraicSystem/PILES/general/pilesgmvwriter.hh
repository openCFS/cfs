#ifndef FILE_PILESGMVWRITER_PILES
#define FILE_PILESGMVWRITER_PILES

/*****************************************************************************/
/* File:   pebblesgmvwriter.hh                                               */
/* Author: Stefan Reitzinger                                                 */
/* Date:   October 2001                                                      */
/*****************************************************************************/


/*
 *   this class provides an interfaces for writing files for gmv
*/

namespace CoupledField
{

class PILES_GMV_Writer : public GMV_Writer
{
public:
  /// Constructor, Destructor
  PILES_GMV_Writer (BaseMesh * amesh, const char * name, 
		    int abinary = 1, const char * afromfile = NULL);

  ///
  virtual ~PILES_GMV_Writer ();

  /// Write gridfunction to file
  void Write (double * sol);

protected:
  /// write nodal data
  void WriteNodes ();

  /// write cell data
  void WriteCells ();

  /// calcualte number of nodes
  virtual long NrNodes () const {return mesh->GetNumNode();};

  /// calculate number of cells
  virtual long NrCells () const {return mesh->GetNumElem();};

private:
  ///
  BaseMesh * mesh;
};
}


#endif    // FILE_PILESGMVWRITER_PILES
