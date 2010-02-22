/*
 * inttable.hh
 *
 *  Created on: 23.07.2009
 *      Author: ahauck
 */

#ifndef FILE_INTEGRATION_SCHEME_HH_
#define FILE_INTEGRATION_SCHEME_HH_

#include "Domain/elem.hh"
#include "Elements/elemShapeMap.hh"

namespace CoupledField {

  //! Integration scheme

  //! The integration scheme returns the integration points / elements for a 
  //! given element shape. As additional parameters one cann pass the 
  //! order (maybe anisotropic) and the type of integration used
  class IntegrationScheme {
  
  public:
    typedef std::map<Elem::FEType,StdVector<LocPoint> > IntegrationPoints;
    typedef std::map<Elem::FEType,StdVector<Double> > IntegrationWeights;
    
    //! Constructor
    IntegrationScheme();

    //! Destructor
    ~IntegrationScheme();
  
    //! Set method and order
    void SetOrder( std::string method, UInt order);

    //! Set method and order
    void SetOrder( IntegrationMethod method, UInt order);

    //! Get integration points / weights
    void GetIntPoints( Elem::FEType elemType,
                         StdVector<LocPoint>& intPts, 
                         StdVector<Double>& weights );  
private:
    //! Adds the Gauss Lobatto Points up to the given order to the Integration maps
    void FillGaussLobattoIntegPoints(UInt order);

    //! Calculate the Gauss-Lobatto points an weights for the given order
    void CalcGaussLobattoPointsWeights(UInt order,StdVector<Double>& points, StdVector<Double>& weights);
  
    //! Map with integration points for each element type According to the Template Paramter Integration Scheme
    std::map<IntegrationMethod, std::map< UInt, IntegrationPoints > > intPoints_;
    
    //! Map with integration weights for each element type
    std::map<IntegrationMethod, std::map< UInt, IntegrationWeights > > intWeights_;

    //! the current method of integration
    IntegrationMethod integMethod_;

    UInt order_;
      
    // Old Integration point imeplementation from BaseFE
//    
//    // ========================================================================
//         //  I N T E G R A T I O N   P O I N T   H A N D L I N G
//         // =======================================================================
//
//         /** Expicitly set and load the integration type and order */
//         void SetIntPoints(IntegrationMethod method, int order);
//
//         /** The child classes add here their integration point data to the elements map. Called in the constructors.
//         * Too avoid errors, any key (type+order) must be unique, otherwise exit!
//         * @param type e.g. ECONOMICAL
//         * @param order the polynomial order
//         * @param numberOfRows the first dimension of data
//         * @param data actually [][2] for 1D (coord + weight) and [][3] for 2D and  [][4] for 3D */
//         void AddIntegrationPoints(IntegrationMethod type, int order, int numberOfRows, Double* data);
//
//         /** returns the shape name as used for integration type in the XML file based on
//               * the abstract feType() */
//         const char* GetShapeName();
//         
//         //! Helper function for printing a coordinate matrix in a string
//         std::string CoordMatrix2String(const Matrix<Double> & coordMat);
//
//
//         /** Creates the integration points by cartesian product of 1D for rectangle and cube.
//         * The creation is quite expensive, but the results are cached in THIS element!
//         * For rectangle, hexaeder (cube) and wedge practically arbitrary high integration orders
//         * can be constructed via the cartesian product rule. For the hexaeder the ECONOMICAL
//         * ones are much more efficient! For anisotropic elements, RectangeFE::SetCartesianIntegration()
//         * and HexaFE::SetCartesianIntegration() creates the elements and stores them in
//         * the map or simple loads it from the map. The internal method type is CARTESIAN and
//         * the order is zyx where z is the order of the coordinate x3 with z in [1;9], ...
//         * This can also be set for testing in the XML file, always set for three coordinates (order > 111).
//         *
//         * When debugging or altering code be carefull, the calling sequence of
//         * SetCartesianIntegration() is not really straight forward. :(
//         *
//         * @param order3 if smaller 1 this dimension is omitted (rectangle)
//         * @param create_only false or leave to default for any user, true only for internal purpose (recursive call) */
//         void SetCartesianInteg(int order1, int order2, int order3, bool create_only);
//
//         /** Helper method that generates the key for the map */
//         void MakeKey(IntegrationMethod type, int order, std::string &out);
//
//         /** generate a key for a cartesian integration
//         * @param line3 optional, use < 1 for 2D */
//         void MakeKey(int order1, int order2, int order3, std::string &out);
//
//         /** Helper function for SetIntPoints. Makes a comfortable search and searches for default when fallback is true
//         * @param fallback caller do true, false for internal recursive search
//         * @return when fallback is true not null but exit with message */
//         StdVector<Double*>* GetIntegrationPoints(IntegrationMethod type, int order, bool search_upwards, bool search_downwards, bool fallback);
//
//         /** Creates the integration points via cartesian product for 2D and 3D elements. The result is stored in the map.
//         *  @param element1 IntegMethod and IntegOrder shall be set
//         *  @param element2 IntegMethod and IntegOrder shall be set
//         *  @param element3 the only optional element */
//         void CreateCartesianIntegration(BaseFE* line1, BaseFE* line2, BaseFE* line3);
//
//
//         /** Private Helper Method that dumps the currently set integration points */
//         void DumpIntegrationPoints();
//         virtual void FillIntegrationPoints();
//
//         /** Sets the integration points for the data filled by the virtual FillIntegrationPoints.
//         * This mehtod reads the data from the map with the help of the dimension */
//         void SetIntPoints();
//
//         /** Every Element has to set the preferred default integration mode */
//         virtual void SetDefaultIntegration();
//
//         /** Every Element has also a default reduced integration set -> see the tex documentation for details */
//         virtual void SetDefaultReducedIntegration();
//
//
//
//         /** encodes the orders of the cartesian product
//         * @param order3 is ignored if < 1 or Dim_ != 3
//         * @return zyx where z = x3 in [1..9], ... */
//         int EncodeCartesianOrder(int order1, int order2, int order3);
//
//         /** private order decoder
//         * @param order3 0 is written here if Dim_ != 3 */
//         void DecodeCartesianOrder(int encoded_order, int* order1, int* order2, int* order3);
//
//
//         
//         UInt NumIntPoints_;
//        Vector<Double> * IntPoints_;      //!< integration points
//        Vector<Double> IntWeights_;       //!< integration weights
//
//        //! AddIntegrationPoints stores here all availabe integrations
//        std::map<const std::string, StdVector<Double*>*> IntegrationPointsMap_;
//
//        enum IntegrationMethod IntegMethod; //! set in XML and SetStandard/ReducedIntegration
//
//        UInt IntegOrder;                    //! accompanies IntegMethod
//        
//        
//        // ========================================================================
//        //  E N D   O F   I N T E G R A T I O N   P O I N T   H A N D L I N G
//        // =======================================================================

};



} // namespace CoupledField

#endif /* INTTABLE_HH_ */
