#ifndef FILE_GMVWRITER_PILES
#define FILE_GMVWRITER_PILES

/*****************************************************************************/
/* File:   gmvwriter.hh                                                      */
/* Author: Wolfram Muehlhuber                                                */
/* Date:   November 1999                                                     */
/*****************************************************************************/


/*
 *   this class provides an interfaces for writing files for gmv
*/


class GMV_Writer
{

public:
  /// Constructor, Destructor
  GMV_Writer (const char * name, int abinary, const char * afromfile);

  ///
  virtual ~GMV_Writer ();

protected:
  /// write file header
  void WriteFileHeader ();

  /// write nodal data
  virtual void WriteNodes () = 0;

  /// write cell data
  virtual void WriteCells () = 0;

  /// calcualte number of nodes
  virtual long NrNodes () const = 0;

  /// calculate number of cells
  virtual long NrCells () const = 0;

  ///
  void WriteLong (long val);
  ///
  void WriteLong (long len, const long * values);
  ///
  void WriteString (const char * str);
  ///
  void WriteFloat (float val);
  ///
  void WriteFloat (long len, const float * val);


  /// name of the outputfile
  char * filename;
  /// ASCII or binary output
  int binary;

  /// take nodal information from a different file
  char * fromfile;

  /// file handle
  ostream * gmvfile;

  /// writing variable data?
  int variablemode;
  /// number of nodes in the gmv file
  long nrnodes;
  /// number of cells in the gmv file
  long nrcells;


};



#endif    // FILE_GMVWRITER_PILES











