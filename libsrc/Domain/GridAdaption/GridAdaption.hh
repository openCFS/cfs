// GridAdaption.hh

#ifndef __GRIDADAPTION_H__
#define __GRIDADAPTION_H__

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "Utils/nodestoresol.hh"
#include "General/environment.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"

//#include "utils.h"

#define FORMAT_BEGINSECTION_COORD     "[Coordinates]"
#define FORMAT_ENDSECTION_COORD       "[/Coordinates]"

#define FORMAT_BEGINSECTION_TIMESTEP  "[Timestep dt=%le]"
#define FORMAT_ENDSECTION_TIMESTEP    "[/Timestep]"

#define FORMAT_COORD                  "\t%d\t%le\t%le\t%le"
#define FORMAT_TIMESTEP               "\t%d\t%le"

#define COMMENT_BEGINSECTION_COORD     "#\tNo.\tx\t\ty\t\tz"
#define COMMENT_BEGINSECTION_TIMESTEP  "#\tNo.\tValue"

namespace CoupledField
{
	class GridAdaption
	{
		public:
			// con-/destructor(s)
			GridAdaption();
			GridAdaption(std::string strDataFile,
				     Grid* pGrid);
			~GridAdaption();

			// structs and enums
			struct TimeInfo
			{
				Double  dTime;
				Double  dValue;
			};
			typedef std::vector<struct TimeInfo> t_vecTimeInfo;

			struct NodeInfo
			{
				Integer nNode;
				Double  x, y, z;
				std::vector<struct TimeInfo> vecTimeInfo;
			};
			typedef std::vector<struct NodeInfo> t_vecNodeInfo;

			struct NodeDistance
			{
				Double dDistance;
				struct NodeInfo *ni;
			};
			typedef std::vector<struct NodeDistance> t_vecNodeDistance;

			enum INTERPOLATION_TYPE
			{
				UNDEFINED = 0,
				NNB,
				IDW,
				SHP,
				RBF
			};

			// methods
			void SetWriteFile(std::string strFile);
			void Init(std::string strDataFile, Grid* pGrid);
			bool ReadFile(std::string strFile);
			bool Add2DataFile(NodeStoreSol<Double>& NodeData, Double dTime);
			void AddNode(struct NodeInfo ni, bool bValue, Double dValue = 0,
			             Double dTime = 0);

			void SetInterpolationType(INTERPOLATION_TYPE it);
			void SetInfluencingNodes(Integer nInfluencingNodes, bool bPercent);
			Double GetAt(Double x, Double y, Double z, Double dTime,
			             INTERPOLATION_TYPE it = UNDEFINED);
			Double InterpolateAt(Double x, Double y, Double z, Double dTime,
			                     INTERPOLATION_TYPE it = UNDEFINED);
			void GetCoordinatesFromGridNode(Grid* pGrid, UInt nNode, Point<3>& pt);

		private:
			// attributes
			bool m_bFirstStep;
			bool m_bPercentual;
			Integer m_nInfluencingNodes;
			Integer m_nP;
			Integer m_nXMLNodes;
			Double m_dEpsDouble;
			Double m_dEpsDistance;
			std::ofstream* m_pDataFile;
			Grid*          m_pGrid;
			INTERPOLATION_TYPE m_nInterpolationType;
			t_vecNodeInfo m_vecReadNodes;
			StdVector<UInt> m_vecXMLNodes;
			std::map<int,int> m_mapNodes;
			std::map<std::string,t_vecNodeDistance> m_mapDistances;

			// methods
			void InitDefaults();
			bool InitXMLInfo();
			inline Integer GetNodeCount();
			Integer GetDistancesCount(t_vecNodeDistance& vec);
			Double InterpolateNNB(Double x, Double y, Double z, Double dTime);
			Double InterpolateIDW(Double x, Double y, Double z, Double dTime);
			Double InterpolateShepard(Double x, Double y, Double z, Double dTime);
			Double InterpolateRBF(Double x, Double y, Double z, Double dTime);

			Double CalcDistance(struct NodeInfo niCenter, struct NodeInfo niAway);
			Double CalcDistance(Double xCenter, Double yCenter, Double zCenter,
			                    Double xAway,   Double yAway,   Double zAway);
			void CalcDistances(struct NodeInfo& niGiven, t_vecNodeDistance& vecDistances);

			void GetCoordinatesFromGridNode(UInt nNode, Point<3>& pt);
			bool GetNodeInfo(Integer nNode, struct NodeInfo **ni);;
			bool GetTimeInfo(Double dTime, struct NodeInfo *ni, struct TimeInfo **ti);
			void InsertSorted(struct NodeDistance& nd, t_vecNodeDistance& vec);
	};

} // namespace CoupledField

#endif //#ifndef __GRIDADAPTION_H__
