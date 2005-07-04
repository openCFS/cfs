// GridAdaption.cpp

#include "GridAdaption.hh"
#include <algorithm>

namespace CoupledField
{
	GridAdaption::GridAdaption()
	{
		InitDefaults();
	}

	GridAdaption::GridAdaption(std::string strDataFile,
				   Grid* pGrid)
	{
		InitDefaults();
		Init(strDataFile, pGrid);
	}

	void GridAdaption::InitDefaults()
	{
		m_pDataFile = NULL;

		// indicator for WriteFile()
		m_bFirstStep = true;

		// don't use percentual fraction of node no.
		m_bPercentual = false;

		// use shepard's method by default
		m_nInterpolationType = NNB;

		// number of nodes that are taken into accout when interpolating
		// -1 means all nodes are considered
		m_nInfluencingNodes = -1;

		// the root of this size is being used for the weighting functions
		m_nP = 2;

		// threshold for comparing Doubles
		m_dEpsDouble = 1e-3;

		// epsilon for when one matches exactly
		m_dEpsDistance = 1e-10;
	}

	void GridAdaption::Init(std::string strDataFile, Grid* pGrid)
	{
		// file, we want to write to
		m_pDataFile = new std::ofstream(strDataFile.c_str());

		// save grid pointer
		m_pGrid = pGrid;
	}

	GridAdaption::~GridAdaption()
	{
		if(m_pDataFile->is_open()) {
			m_pDataFile->close();
		}
	}

	// sets the type of interpolation that is to be used
	void GridAdaption::SetInterpolationType(INTERPOLATION_TYPE it)
	{
		m_nInterpolationType = it;
	}

	// here we can say if we want to take all nodes into account or only some,
	// that will have an effect on the weights of each node
	void GridAdaption::SetInfluencingNodes(Integer nInfluencingNodes,
	                                            bool bPercentual)
	{
		m_nInfluencingNodes = nInfluencingNodes;
		m_bPercentual       = bPercentual;
	}

	inline Integer GridAdaption::GetNodeCount()
	{
		return m_vecReadNodes.size();
	}

	Integer GridAdaption::GetDistancesCount(t_vecNodeDistance& vec)
	{
		// consider all nodes?
		if(m_nInfluencingNodes == -1 || m_nInfluencingNodes >= GetNodeCount()) {
			return vec.size();
		}

		return m_nInfluencingNodes;
	}

	// get f(x,y,z) at timestep dTime
	Double GridAdaption::GetAt(Double x, Double y, Double z, Double dTime,
	                                INTERPOLATION_TYPE it)
	{
		ENTER_FCN("GridAdaption::GetAt()");

		return InterpolateAt(x, y, z, dTime, it);
	}

	// get f(x,y,z) at timestep dTime
	Double GridAdaption::InterpolateAt(Double x, Double y, Double z,
	                                        Double dTime, INTERPOLATION_TYPE it)
	{
		if(it != UNDEFINED) {
			SetInterpolationType(it);
		}

		// decide how to interpolate
		switch(m_nInterpolationType) {
			case NNB: return InterpolateNNB(x, y, z, dTime);
			case IDW: return InterpolateIDW(x, y, z, dTime);
			case SHP: return InterpolateShepard(x, y, z, dTime);
			case RBF: return InterpolateRBF(x, y, z, dTime);
			default:  return 0;
		}
	}

	// Nearest Neighbor
	Double GridAdaption::InterpolateNNB(Double x, Double y, Double z,
	                                         Double dTime)
	{
		ENTER_FCN("GridAdaption::InterpolateNNB()");

		// setup one nodeinfo struct
		struct NodeInfo niGiven;
		niGiven.x = x;
		niGiven.y = y;
		niGiven.z = z;

		// calculate relevant distances
		t_vecNodeDistance vecDistances;
		CalcDistances(niGiven, vecDistances);

		struct TimeInfo *ti;
		GetTimeInfo(dTime, vecDistances[0].ni, &ti);

		return ti->dValue;
	}

	// Inverse Distance Weighted
	//
	// f(x,y) = \sum(f_i/d_i) / \sum(1/d_i)
	// d_i = \sqrt((x-x_i)^2 + (y-y_i)^2 + (z-z_i)^2)
	Double GridAdaption::InterpolateIDW(Double x, Double y, Double z,
	                                         Double dTime)
	{
		ENTER_FCN("GridAdaption::InterpolateIDW()");

		// setup one nodeinfo struct
		struct NodeInfo niGiven;
		niGiven.x = x;
		niGiven.y = y;
		niGiven.z = z;

		// calculate relevant distances
		t_vecNodeDistance vecDistances;
		CalcDistances(niGiven, vecDistances);

		Double dDistanceCurrent   = 0;
		Double dDistanceTotal     = 0;
		Double dDistancesWeighted = 0;
		Integer nDistances = GetDistancesCount(vecDistances);

		struct TimeInfo *ti;

		// now in vecDistances are either only the nearest nodes, or all
		// => we can perform our calculation
		for(Integer i=0; i<nDistances; i++) {
			dDistanceCurrent = vecDistances[i].dDistance;

			// get info about this timestep
			GetTimeInfo(dTime, vecDistances[i].ni, &ti);

			// check if wanted node is exactly the value of an already known one
			if(dDistanceCurrent < m_dEpsDistance) { // < 1.0e-10
				return ti->dValue;
			}

			// sum up all distances
			dDistanceTotal += 1/dDistanceCurrent;

			dDistancesWeighted += ti->dValue / dDistanceCurrent;
		}

		if(dDistancesWeighted == 0) {
			return 0;
		}

		return dDistancesWeighted / dDistanceTotal;
	}

	// Modified SHEPARD's method
	//
	// f(x,y) = \sum_i(w_i*f_i)
	//
	// w_i = ( ((R-d_i) / (R*d_i)) ^2) / \sum_j( ((R-d_j) / (R*d_j)) ^2)
	// d_i = \sqrt((x-x_i)^2 + (y-y_i)^2 + (z-z_i)^2)
	// R = distance from interpolation point to most distant scatter point
	Double GridAdaption::InterpolateShepard(Double x, Double y, Double z,
	                                             Double dTime)
	{
		ENTER_FCN("GridAdaption::InterpolateShepard()");

		// setup one nodeinfo struct
		struct NodeInfo niGiven;
		niGiven.x = x;
		niGiven.y = y;
		niGiven.z = z;

		// calculate relevant distances
		t_vecNodeDistance vecDistances;
		CalcDistances(niGiven, vecDistances);

		Integer nDistances = GetDistancesCount(vecDistances);
		Double  dDistanceTotal   = 0;
		Double  dDistanceCurrent = 0;

		// now in vecDistances are either only the nearest nodes, or all
		// => we can perform our calculation

		// first, sum up all distances
		for(Integer i=0; i<nDistances; i++) {
			dDistanceCurrent = vecDistances[i].dDistance;

			// check if wanted node is exactly the value of an already known one
			if(dDistanceCurrent < m_dEpsDistance) { // < 1.0e-10
				struct TimeInfo *ti;
				GetTimeInfo(dTime, vecDistances[i].ni, &ti);

				return ti->dValue;
			}

			// when regarding only one node R would be equal to dDistanceCurrent
			if(nDistances > 1) {
				// weight calculation
				Double R   = vecDistances[nDistances-1].dDistance;
				Double d_i = vecDistances[i].dDistance;
				dDistanceCurrent = pow((R - d_i) / (R * d_i), 2);

				// overwrite value in vecDistances, we don't need the actual
				// distance any more
				vecDistances[i].dDistance = dDistanceCurrent;
			}

			// add to total distance
			dDistanceTotal += dDistanceCurrent;
		}

		Double dResult = 0;
		Double dWeight = 0;
		struct TimeInfo *ti;

		// only now we can efficiently calculate the weight for each regarded node
		for(Integer i=0; i<nDistances; i++) {
			dWeight = vecDistances[i].dDistance / dDistanceTotal;

			GetTimeInfo(dTime, vecDistances[i].ni, &ti);

			// the contribution of the current node
			dResult += dWeight * ti->dValue;
		}

		return dResult;
	}

	// Radial Basis Functions
	//
	// NOT IMPLEMENTED YET
	//
	// f(x,y) = \sum c_i*h_i
	// h_i = h(r) = R^2 + r^2 | gaussian | ...
	Double GridAdaption::InterpolateRBF(Double x, Double y, Double z,
	                                         Double dTime)
	{
		Double dResult = 0;
		return dResult;
	}

	// adds a node to our intern data structure
	void GridAdaption::AddNode(struct NodeInfo ni, bool bValue,
	                                Double dValue, Double dTime)
	{
		// look for node with correct coordinates and add value of ni-struct
		if(bValue) {
			Integer nIndex = m_mapNodes[ni.nNode];
			// simply add new timestep info (we're only adding, not inserting!)
			struct TimeInfo ti;
			ti.dTime  = dTime;
			ti.dValue = dValue;
			m_vecReadNodes[nIndex].vecTimeInfo.push_back(ti);
		}

		// add coordinate
		else {
			Integer nNodes = m_vecReadNodes.size();

			// add to node vector
			m_vecReadNodes.push_back(ni);

			// save position in map
			m_mapNodes[ni.nNode] = nNodes;
		}
	}

	bool GridAdaption::ReadFile(std::string strFile)
	{
		ENTER_FCN("GridAdaption::ReadFile()");

		std::ifstream *pInFile = new std::ifstream(strFile.c_str());

		if(!pInFile->is_open()) {
			std::cerr << "ERROR: cannot read from file '" << strFile << "'!" << std::endl;
			exit(EXIT_FAILURE);
			return false;
		}

		bool bCoordinates = false, bTimestep = false;
		std::string strLine;
		struct NodeInfo ni;
		Double dValue, dTime;

		// read in all lines of file
		while(getline(*pInFile, strLine)) {
			// check for comment line
			if(strLine[0] == '#' || strLine == "") {
				continue;
			}

			// now we check for section separators, they begin with '['
			if(strLine[0] == '[') {
				std::string::size_type pos;

				if(strLine.find(FORMAT_BEGINSECTION_COORD) == 0) {
					bCoordinates = true;
				}
				else if(strLine.find(FORMAT_ENDSECTION_COORD) == 0) {
					bCoordinates = false;
				}
				else if(strLine.find(FORMAT_ENDSECTION_TIMESTEP) == 0) {
					bTimestep = false;
				}
				// here we need to extract the data from the attributes "step" and "dt"
				else if(  (pos = strLine.find("[Timestep dt=")) == 0
				        && pos != std::string::npos) {
					Integer iRet = sscanf(strLine.c_str(), FORMAT_BEGINSECTION_TIMESTEP, &dTime);

					if(iRet == EOF) {
						std::cerr << "ERROR: reading line '" << strLine << "'!" << std::endl;
						exit(EXIT_FAILURE);
					}

					bTimestep = true;
				}
				// if we're here, this line isn't valid!
				else {
					std::cerr << "ERROR: invalid input line '" << strLine << "'!" << std::endl;
					exit(EXIT_FAILURE);
				}

				// and continue with reading data till end of section
				continue;
			}

			// if we come this far, strLine (should) contain(s) actual data info
			Integer iRet = 0;
			if(bCoordinates) {
				iRet = sscanf(strLine.c_str(), FORMAT_COORD, &ni.nNode, &ni.x, &ni.y, &ni.z);

				// check if z- or/and y-coordinates are missing => set to 0
				if(iRet < 4) {
					ni.z = 0;

					// not even y-coordinate there?
					if(iRet < 3) {
						ni.y = 0;
					}
				}

				// false indicates coordinate node
				AddNode(ni, false);
			}
			else if(bTimestep) {
				iRet = sscanf(strLine.c_str(), FORMAT_TIMESTEP, &ni.nNode, &dValue);
				// the true-flag means to lookup node and add value
				AddNode(ni, true, dValue, dTime);
			}

			// s.th. bad happended while reading line
			// or line is malformed (e.g. x-coordinate missing => iRet < 2)
			if(iRet == EOF || iRet < 2) {
				std::cerr << "ERROR: invalid input line '" << strLine << "'!" << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		// close file again
		pInFile->close();

		// only now we can calculate percentual fraction of read nodes
		if(m_bPercentual) {
			m_nInfluencingNodes = (int) ((m_nInfluencingNodes/100.0) * GetNodeCount());
		}

		return true;
	}

	bool GridAdaption::Add2DataFile(NodeStoreSol<Double>& NodeData, Double dTime)
	{
		ENTER_FCN("GridAdaption::Add2DataFile()");

		Integer nDim = m_pGrid->GetDim();

		if(nDim != 2 && nDim != 3) {
			std::cerr << "ERROR: GridAdaption only possible for 2D problems!" << std::endl;
			std::cerr << "\tdimension given: " << m_pGrid->GetDim() << std::endl;
			exit(EXIT_FAILURE);
			return false;
		}

		if(m_pDataFile == NULL || !m_pDataFile->is_open()) {
			std::cerr << "ERROR: DataFile hasn't been opened correctly!" << std::endl;
			exit(EXIT_FAILURE);
			return false;
		}

		Integer  nNodeNo, nDofs;
		Double   val;
		Point<3> pt3D;

		// retrieve solution types
		StdVector<SolutionType> vecSolTypes;
		NodeData.GetSolutionTypes(vecSolTypes);

		char szBuffer[512];

		// which dof?
		nDofs = NodeData.GetDof(vecSolTypes[0]);

		// in first step, print coordinates
		if(m_bFirstStep) {
			m_bFirstStep = false;

			// also do some things that need to be done only once
			InitXMLInfo();

			// here begins the info about coordinates
			(*m_pDataFile) << FORMAT_BEGINSECTION_COORD << std::endl;
			// add "table header"
			(*m_pDataFile) << COMMENT_BEGINSECTION_COORD << std::endl;

			// iterate over all history nodes
			for(UInt iNode = 0; iNode<m_nXMLNodes; iNode++) {
				// current node no.
				nNodeNo = m_vecXMLNodes[iNode];

				// get coordinates, as it's 0-based, decrement no.
				GetCoordinatesFromGridNode(nNodeNo, pt3D);

				// setup line 2 print
				sprintf(szBuffer, FORMAT_COORD, nNodeNo, pt3D[0], pt3D[1], pt3D[2]);
				(*m_pDataFile) << szBuffer << std::endl;
			}

			// done with this section
			(*m_pDataFile) << FORMAT_ENDSECTION_COORD << std::endl << std::endl;
		}

		// we're at timestep dTime
		sprintf(szBuffer, FORMAT_BEGINSECTION_TIMESTEP, dTime);
		(*m_pDataFile) << szBuffer << std::endl;

		// add "table header"
		(*m_pDataFile) << COMMENT_BEGINSECTION_TIMESTEP << std::endl;

		// iterate over all special nodes
		for(Integer iNode = 0; iNode<m_nXMLNodes; iNode++) {
			// current node no.
			nNodeNo = m_vecXMLNodes[iNode];

			// iterate over all dofs
			for(Integer iDof = 0; iDof<nDofs; iDof++) {
				// get current value
				NodeData.Get(nNodeNo-1, iDof, val);

				// print output line
				sprintf(szBuffer, FORMAT_TIMESTEP, nNodeNo, val);
				(*m_pDataFile) << szBuffer << std::endl;
			}
		}

		// done with this timestep
		(*m_pDataFile) << FORMAT_ENDSECTION_TIMESTEP << std::endl << std::endl;

		return true;
	}

	bool GridAdaption::InitXMLInfo()
	{
		ENTER_FCN("GridAdaption::InitXMLInfo()");

		StdVector<std::string> vecKey, vecAttr, vecVals;
		StdVector<std::string> vecNodes;

		vecAttr = "", "";
		vecVals = "", "";

		// extract information from xml-file
		vecKey = "storeResults", "specialNodes", "saveNodes";
		params->GetList(vecKey, vecAttr, vecVals, vecNodes);

		/*
		vecKey = "storeResults", "nodeSpecial", "type";
		params->GetList(vecKey, vecAttr, vecVals, vecQuants);

		// we can only process acoustic potential
		if(vecQuants.GetSize() != 1) {
			std::cerr << "ERROR: only ACOUSTIC POTENTIAL can be stored!" << std::endl;
			return false;
		}
		*/

		// now read in our node numbers
		//		m_pXMLFile->ReadSaveNodes(m_vecXMLNodes, vecNodes[0]);
		m_pGrid->GetNodesByName(m_vecXMLNodes, vecNodes[0]);

		// how many nodes are there?
		m_nXMLNodes = m_vecXMLNodes.GetSize();

		return true;
	}

	void GridAdaption::InsertSorted(struct NodeDistance& nd,
	                                     t_vecNodeDistance& vec)
	{
		Integer nIndex = -1;
		Integer nSize = vec.size();

		for(Integer i=0; i<nSize; i++) {
			if(nd.dDistance < vec[i].dDistance) {
				nIndex = i;
				break;
			}
		}

		// found location to insert
		if(nIndex != -1) {
			// duplicate last node, if we consider all nodes (i.e.
			// m_nInfluencingNodes == -1),
			// otherwise, node with biggest distance falls away
			if(m_nInfluencingNodes == -1 || nSize < m_nInfluencingNodes) {
				vec.push_back(vec[nSize-1]);
			}

			for(Integer j=nSize-1; j>nIndex; j--) {
				vec[j] = vec[j-1];
			}

			// insert given node at wanted location
			vec[nIndex] = nd;
		}
		// if no insertion point found, append node
		else {
			vec.push_back(nd);
		}
	}

	// store info of nNode in ni
	// returns true on success, and false if nNode does'nt exist
	bool GridAdaption::GetNodeInfo(Integer nNode, struct NodeInfo** ni)
	{
		// fill struct
		(*ni) = &(m_vecReadNodes[nNode]);

		return true;
	}

	bool GridAdaption::GetTimeInfo(Double dTime, struct NodeInfo *ni,
	                                    struct TimeInfo **ti)
	{
		Integer nSteps = ni->vecTimeInfo.size() - 1;
		Integer iFirst = 0;
		Integer iLast  = nSteps;
		Integer iMiddle;
		Double dRelNorm, dCurrTime;

		// do binary search
		while(iFirst <= iLast) {
			iMiddle   = (iFirst + iLast) / 2;
			dCurrTime = ni->vecTimeInfo[iMiddle].dTime;
			dRelNorm  = fabs(dTime - dCurrTime) / dTime;

			// value has been found, bail out
			if(dRelNorm < m_dEpsDouble || dCurrTime == dTime) {
				(*ti) = &(ni->vecTimeInfo[iMiddle]);
				return true;
			}
			else if(dTime < dCurrTime) {
				iLast = iMiddle - 1;
			}
			else { //if(dTime > dCurrTime)
				iFirst = iMiddle + 1;
			}
		} // end while loop

		// timestep wasn't found
		// too few timesteps => pad last value
		if(iFirst >= nSteps) {
			(*ti) = &(ni->vecTimeInfo[nSteps-1]);
			//fprintf(stderr, "last step, f: %d, m: %d, l:%d (n: %d), dt: %e\n", iFirst, iMiddle, iLast, nSteps, dTime);
		}
		// timesteps is in between => interpolate
		else {
			// assign any timeinfo struct, doesn't matter, only dValue is used later on
			(*ti) = &(ni->vecTimeInfo[iMiddle]);

			/*
			// determine timesteps before and after dTime
			Integer iPrev, iNext;
			if(dCurrTime > dTime) {
				iPrev = iMiddle-1;
				iNext = iMiddle;
			}
			else {
				iPrev = iMiddle;
				iNext = iMiddle+1;
			}

			Double dPrevTime = ni->vecTimeInfo[iPrev].dTime;
			Double dPrevVal  = ni->vecTimeInfo[iPrev].dValue;
			Double dNextTime = ni->vecTimeInfo[iNext].dTime;
			Double dNextVal  = ni->vecTimeInfo[iNext].dValue;

			//fprintf(stderr, "interpolate (%d,%d,%d,%d), dt: %e, dt-1: %e, dt+1: %e\n", iFirst, iMiddle, iLast, nSteps, dTime, dPrevTime, dNextTime);

			// simply assign current struct, only the dValue member is used later on
			(*ti) = &(ni->vecTimeInfo[iMiddle]);

			// as timestep doesn't exist, we interpolate between previous and next
			Double u = (dTime - dPrevTime) / (dNextTime - dPrevTime);
			(*ti)->dValue = (1-u)*dPrevVal + u*dNextVal;
			//fprintf(stderr, "value = %e (-: %e, +: %e)\n", (*ti)->dValue, dPrevVal, dNextVal);
			*/
		}

		return false;
	}

	void GridAdaption::CalcDistances(struct NodeInfo& niGiven,
	                                      t_vecNodeDistance& vecDistances)
	{
		Integer nNodes = GetNodeCount();
		struct NodeInfo    *ni;
		struct NodeDistance nd;

		// generate key for map
		char szKey[64];
		sprintf(szKey, "%.10e,%.10e,%.10e", niGiven.x, niGiven.y, niGiven.z);

		// not calculated yet, add to given distance vector
		if(m_mapDistances.find(szKey) == m_mapDistances.end()) {
			for(Integer i=0; i<nNodes; i++) {
				GetNodeInfo(i, &ni);

				nd.ni        = ni;
				nd.dDistance = CalcDistance(niGiven, *ni);

				InsertSorted(nd, vecDistances);
			}

			// now save distances to map
			m_mapDistances[szKey] = vecDistances;
		}
		// already calculated, simply lookup
		else {
			vecDistances = m_mapDistances[szKey];
		}
	}

	// returns distance between specified two nodes
	Double GridAdaption::CalcDistance(struct NodeInfo niCenter,
	                                       struct NodeInfo niAway)
	{
		return CalcDistance(niCenter.x, niCenter.y, niCenter.z, niAway.x, niAway.y, niAway.z);
	}

	// calculate distance between (xCenter, yCenter, zCenter) and
	// (xAway, yCenter, zCenter)
	Double GridAdaption::CalcDistance(Double xCenter, Double yCenter, Double zCenter,
	                                       Double xAway,   Double yAway,   Double zAway)
	{
		Double dx = xCenter - xAway;
		Double dy = yCenter - yAway;
		Double dz = zCenter - zAway;

		return sqrt(dx*dx + dy*dy + dz*dz);
	}

	void GridAdaption::GetCoordinatesFromGridNode(UInt nNode, Point<3>& pt)
	{
		GetCoordinatesFromGridNode(m_pGrid, nNode, pt);
	}

	void GridAdaption::GetCoordinatesFromGridNode(Grid* pGrid, UInt nNode, Point<3>& pt)
	{
		ENTER_FCN("GridAdaption::GetCoordinatesFromGridNode()");

		std::cout << "node=" << nNode << std::endl;
		if(pGrid->GetDim() == 2) {
			Point<2> pt2D;
			pGrid->GetNodeCoordinate(pt2D, nNode);

			pt[0] = pt2D[0];
			pt[1] = pt2D[1];
			pt[2] = 0;
		}
		else {
			Point<3> pt3D;
			pGrid->GetNodeCoordinate(pt3D, nNode);

			pt[0] = pt3D[0];
			pt[1] = pt3D[1];
			pt[2] = pt3D[2];
		}
	}

} // namespace CoupledField
