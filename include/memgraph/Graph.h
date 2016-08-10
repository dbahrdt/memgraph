#ifndef MEMGRAPH_GRAPH_H
#define MEMGRAPH_GRAPH_H
#include <vector>
#include <limits>
#include <stdint.h>
#include <string>

namespace memgraph {

class Graph {
public:
	struct NodeInfo final {
		int64_t osmId;
		double lat;
		double lon;
		NodeInfo(int64_t osmId, double lat, double lon) : osmId(osmId), lat(lat), lon(lon) {}
	};
	
	struct Node final {
		//offset into edges array
		uint32_t begin;
		//offset into edges array, one-passed-the-end
		uint32_t end;
		Node() : begin(0xFFFFFFFF), end(begin) {}
		uint32_t edgeCount() const { return end-begin; }
	};
	
	struct Edge final {
		typedef enum {ET_INVALID=0,
			ET_BEGIN=1,
			ET_MOTORWAY=1, ET_TRUNK=2, ET_PRIMARY=3, ET_SECONDARY=4, ET_TERTIARY=5, ET_UNCLASSIFIED=6, ET_RESIDENTIAL=7, ET_SERVICE=8,
			ET_MOTORWAY_LINK=9, ET_TRUNK_LINK=10, ET_PRIMARY_LINK=11, ET_SECONDARY_LINK=12, ET_TERTIARY_LINK=13,
			ET_LIVING_STREET=14, ET_PEDESTRIAN=15, ET_TRACK=16, ET_BUS_GUIDE_WAY=17, ET_RACEWAY=18, ET_ROAD=19,
			ET_FOOTWAY=20, ET_BRIDLEWAY=21, ET_STEPS=22, ET_PATH=23, ET_CYCLEWAY=24,
			ET_BUS_STOP=25, ET_PLATFORM=26,
			ET_END=25
		} EdgeTypes;
		
		typedef enum {
			AT_NONE=0, AT_FOOT=0x1, AT_BIKE=0x2, AT_CAR=0x4, AT_ALL=0x7
		} AccessTypes;
		constexpr static const char * edge_type_2_osm_highway_value[] = {
			"invalid",
			"motorway", "trunk", "primary", "secondary", "tertiary", "unclassified", "residential", "service",
			"motorway_link", "trunk_link", "primary_link", "secondary_link", "tertiary_link",
			"living_street", "pedestrian", "track", "bus_guide_way", "raceway", "road",
			"footway", "bridleway", "steps", "path", "cycleway",
			"bus_stop", "platform"
		};
		//edge speeds in km/h
		constexpr static const int edge_type_2_speed[] = {
			0,
			130, 120, 100, 100, 80, 50, 30, 10,
			130, 120, 100, 100, 80,
			5, 5, 10, 5, 300, 10,
			5, 5, 5, 5, 15,
			5, 5
		};
		
		constexpr static const int edge_type_2_access_allowance[] = {
			AT_ALL,//allow all by default
			AT_CAR, AT_CAR, AT_CAR, AT_ALL, AT_ALL, AT_ALL, AT_ALL, AT_ALL,
			AT_CAR, AT_CAR, AT_CAR, AT_ALL, AT_ALL,
			AT_ALL, AT_FOOT, AT_ALL, AT_ALL, AT_CAR, AT_ALL,
			AT_FOOT|AT_BIKE, AT_ALL, AT_FOOT, AT_FOOT|AT_BIKE, AT_BIKE,
			AT_ALL, AT_FOOT
		};
		
		uint32_t source;//this is needed to draw the edge
		uint32_t target;
		uint32_t oneway:1;
		uint32_t type:5; //one of EdgeTypes
		uint32_t access:3; //one of AccessTypes
		uint32_t speed:20;//in mm/s
		uint32_t distance:30; //in mm
		Edge() : source(0xFFFFFFFF), target(0xFFFFFFFF), oneway(0), type(ET_INVALID), access(AT_NONE), speed(0), distance(0) {}
	};
	
	struct Route {
		std::vector<uint32_t> nodes;
		double time;
		double distance;
		int at;
		Route() : time(0), distance(0), at(Edge::AT_NONE) {}
		Route(const Route & other) : nodes(other.nodes), time(other.time), distance(other.distance), at(other.at) {}
		Route(Route && other) : nodes(std::forward< std::vector<uint32_t> >(other.nodes)), time(other.time), distance(other.distance), at(other.at) {}
		~Route() {}
	};
	
	typedef std::vector<NodeInfo> NodeInfoContainer;
	typedef std::vector<Node> NodeContainer;
	typedef std::vector<Edge> EdgeContainer;
	
	typedef EdgeContainer::const_iterator ConstEdgeIterator;
	
public:
	Graph();
	Graph(const Graph & other);
	Graph(Graph && other);
	virtual ~Graph();
	Graph & operator=(const Graph & other);
	Graph & operator=(Graph && other);
	
	inline uint32_t nodeCount() const { return  m_nodes.size(); }
	inline uint32_t edgeCount() const { return m_edges.size(); }
	
	const NodeInfo & nodeInfo(uint32_t pos) const { return m_nodeInfo.at(pos); }
	const Node & node(uint32_t pos) const { return m_nodes.at(pos); }
	ConstEdgeIterator edgesBegin(uint32_t nodeId) const { return m_edges.cbegin()+m_nodes.at(nodeId).begin; }
	ConstEdgeIterator edgesEnd(uint32_t nodeId) const { return m_edges.cbegin()+m_nodes.at(nodeId).end; }
	const Edge & edge(uint32_t pos) const { return m_edges.at(pos); }
	
	///@param vehicleMaxSpeed in km/h
	Route routeInfo(std::vector<uint32_t> && route, double vehicleMaxSpeed, int accessType);
	
	///complexity: O(|Nodes|)
	void bbox(double & minLat, double & maxLat, double & minLon, double & maxLon) const;
	
	void printStats(std::ostream & out) const;
	
	bool selfCheck();
	
public:
	static Graph fromPBF(const std::string & path);
private:
	NodeInfoContainer m_nodeInfo;
	NodeContainer m_nodes;
	EdgeContainer m_edges;
};

}//end namespace memgraph

#endif