

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Ray_2.h>
#include <iostream>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/Random.h>
#include <CGAL/algorithm.h>
#include <CGAL/Line_2.h>
#include <CGAL/Segment_3.h>
//Geomview 
#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Polyhedron_geomview_ostream.h>

#include <CGAL/IO/Triangulation_geomview_ostream_2.h>
#include <CGAL/IO/Triangulation_geomview_ostream_3.h>

//Random Point Set Creation
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

#include <cassert>
#include <CGAL/intersections.h>
#include <CGAL/Polyhedron_3.h>
#include <typeinfo>
#include <boost/optional/optional_io.hpp>
#include <CGAL/enum.h>
#include <fstream>
#include <CGAL/squared_distance_2.h> //for 2D functions
#include <CGAL/centroid.h>

#include <CGAL/convex_hull_3.h>
#include <CGAL/Tetrahedron_3.h>

#include <unordered_map>
#include <string>
#include <CGAL/centroid.h>

struct FaceInfo2
{

  int nesting_level;
  bool visited;
  bool in_domain(){ 
    return nesting_level%2 == 1;
  }
    FaceInfo2(){
  		visited = false;
 	 }
};
typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;
typedef CGAL::Triangulation_vertex_base_2<K>                      Vb;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2,K>    Fbb;
typedef CGAL::Constrained_triangulation_face_base_2<K,Fbb>        Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>               TDS;
typedef CGAL::Exact_predicates_tag                                Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>  CDT;
typedef CDT::Point                                                Point;
//typedef CGAL::Constrained_triangulation_plus_2<CDT> 			  CDTplus;	
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef K::Point_2 Point_2;
typedef K::Triangle_2 Triangle;
typedef std::vector<Point> Vector;
typedef K::Intersect_3 Intersect_3;
typedef K::Point_3 Point_3;
typedef K::Line_2 Line;
typedef K::Segment_3 Segment;
typedef K::Tetrahedron_3 Tetrahedron_3;


typedef std::list<Point>                                Container;
typedef CGAL::Polygon_2<K, Container>                     Polygon_2;
typedef K::Ray_2  Ray;
typedef CGAL::Creator_uniform_2<int, Point>             Creator;
typedef CGAL::Random_points_in_square_2<Point, Creator> Point_generator;

typedef CGAL::Delaunay_triangulation_2<K> Delaunay_triangulation_2;

using namespace std;
using std::cout; using std::endl;




std::vector<Point_3> getNodes(CGAL::Geomview_stream& gv) {
	ifstream input("example.1.node");
	double id,x, y,z;

	std::vector<Point_3> nodes;
	input >> id >> x >> y >> z ;
	while (input >> id >> x >> y >> z){
	    cout << id << ": ("<< x << "," << y <<  "," << z << ")"<<endl;
    	Point_3 p(x,y,z);
    	nodes.push_back(p);
	    
	}    
	return nodes;
}





std::vector<Point_3> PATH;

std::vector<Tetrahedron_3> getTetrahedra(std::vector<Point_3>& nodes,CGAL::Geomview_stream& gv,unordered_map<int,int>&  umap) {
	ifstream input("example.1.ele");
	double id, a,b,c,d;
	Point_3 n1,n2,n3,n4;
	std::vector<Tetrahedron_3> tetrahedra;

	input >> a >> b >> c;

	while (input >> id >> a >> b>> c >> d){
	    cout << id << ": ("<< a << "," << b <<  "," << c << ","<< d <<")"<<endl;
	    umap[a] = umap[b] = umap[c] = umap[d] = tetrahedra.size();
	    Point_3 n1(nodes[a-1]);
	    Point_3 n2(nodes[b-1]);
	    Point_3 n3(nodes[c-1]);
	    Point_3 n4(nodes[d-1]);
	    cout << id << ": [("<< n1 << "),(" << n2 <<  "),(" << n3 << "),(" << n4 <<")]"<<endl<<endl;
    	Tetrahedron_3 p(n1,n2,n3,n4);
    	tetrahedra.push_back(p);

	}    
	return tetrahedra;
}



std::vector< std::vector<int> > getEdges(CGAL::Geomview_stream& gv,int nodes_num ) {
	ifstream input("example.1.edge");
	double id,n1,n2;
	std::vector< std::vector<int> > edges;
	std::vector<int> v;
	cout << nodes_num << endl;
	for (int i = 0; i < nodes_num; ++i){
		edges.push_back(v);
	}
	input >> id;
	while (input >> id >> n1 >> n2){
	    cout << id << ": " << n1 << " - " << n2 << endl;
    	edges[n1-1].push_back(n2-1);
    	edges[n2-1].push_back(n1-1);
	}  
	return edges;
}



std::vector<Tetrahedron_3> getAdjacentTetrahedra(std::vector<Point_3> nodes, Tetrahedron_3& T, std::vector<Tetrahedron_3>& tetrahedra,
	std::vector< std::vector<int> >& edges,unordered_map<int,int>& nodeMap,std::vector<bool>& visited){
	
	std::vector<Tetrahedron_3> v;
	std::vector<Point_3> TP;
	Point_3 a = T[0];
	Point_3 b = T[1];
	Point_3 c = T[2];
	Point_3 d = T[3];
	TP.push_back(a);
	TP.push_back(b);
	TP.push_back(c);
	TP.push_back(d);
	for (int i = 0; i < 4; ++i){
		Point_3 p = T[i];
		int p_pos = 0;
		for (int j = 0; j < nodes.size(); ++j){
			if(nodes[j] == p){
				p_pos = j;
			}
		}
		vector<int> neighbors = edges[p_pos];
		for (int n = 0; n < neighbors.size() ; ++n)	{
			int tetid = nodeMap[neighbors[n]];
			if(visited[tetid] == true)continue;
			Tetrahedron_3 tetrahedron =  tetrahedra[tetid];
			v.push_back(tetrahedron);
			visited[tetid] = true;
		}


	}
	return v;

}



void DepthFirst(std::vector<Point_3> nodes, Tetrahedron_3& T, Point_3& G, std::vector<Tetrahedron_3>& tetrahedra,
std::vector< std::vector<int> >& edges,unordered_map<int,int>& nodeMap,std::vector<bool>& visited,CGAL::Geomview_stream& gv){
	
	gv << CGAL::RED;
	std::vector<Tetrahedron_3> neighbors = getAdjacentTetrahedra(nodes,T,tetrahedra,edges,nodeMap,visited);
	
	if(neighbors.empty() == true){
		if(PATH.empty() == false){
			PATH.pop_back();
		}
	}
	
	for (int i = 0; i < neighbors.size(); ++i){
		//sleep(1);
		//gv << neighbors[i];

		std::vector<Point_3> v;
		for(int p = 0 ; p < 4 ; p++){
			v.push_back(neighbors[i][p]);
			if(neighbors[i][p] == G){
				PATH.push_back(G);
				gv.clear();
				gv<<CGAL::GREEN;
				
				for(int x = 1 ; x < PATH.size() ; x++){
					//gv<< PATH[x];
					Segment connect(PATH[x-1],PATH[x]);
					gv << connect;
				}
				return;
			}

		}

		Point_3 center = CGAL::centroid(v.begin(),v.end());
		//gv << center;
		PATH.push_back(center);
		

		DepthFirst(nodes,neighbors[i],G,tetrahedra,edges,nodeMap,visited,gv);
	}
}



Polyhedron PDB2Poly(string typex1,string typex2,string typex3){

	std::vector<string> angletypes;
	angletypes.push_back(typex1);
	angletypes.push_back(typex2);
	angletypes.push_back(typex3);

	std::vector<double> minangles,maxangles;
	for(int i = 0 ; i < 3 ; i++){
		string type = angletypes[i];
		if (type == "g+"){
			minangles.push_back(0.0);
			maxangles.push_back(120.0);
		}
		else if (type == "t"){
			minangles.push_back(120.0);
			maxangles.push_back(240.0);			
		}
		else if (type == "g+"){
			minangles.push_back(240.0);
			maxangles.push_back(360.0);			
		}
	}

	std::vector<Point_3> pointset;

	for(int i = 0 ; i < 2 ; i ++){
		double c1,c2,c3;
		if(i == 0){
			c1 = minangles[0];
		}
		else{
			c1 = maxangles[0];
		}

		for(int j = 0 ; j < 2 ; j ++){

			if(j == 0){
				c2 = minangles[1];
			}
			else{
				c2 = maxangles[1];
			}

			for(int k = 0 ; k < 2 ; k ++){
				if (k == 0){
					c3 = minangles[2];
				}	
				else{
					c3 = maxangles[2];
				}
				Point_3 p(c1,c2,c3);
				cout<<"("<<c1<<","<<c2<<","<<c3<<")"<<endl;
				pointset.push_back(p);

			}
		}
	}

	Polyhedron poly;

  	CGAL::convex_hull_3(pointset.begin(), pointset.end(), poly);

 	std::cout << "The convex hull contains " << poly.size_of_vertices() << " vertices" << std::endl;
  	typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_facet_circulator;


 	for(Polyhedron::Facet_iterator pit = poly.facets_begin(); pit!= poly.facets_end() ; ++pit){

 		Halfedge_facet_circulator j = pit->facet_begin();
	    CGAL_assertion( CGAL::circulator_size(j) >= 3);
	    std::cout << CGAL::circulator_size(j) << ' ';
	    cout<<"facet"<<endl;
	    do{
	    		cout << "\tvertex:" << j->vertex()->point() << endl;

	      //std::cout << ' ' << std::distance(poly.vertices_begin(), j->vertex());
	    }while ( ++j != pit->facet_begin());

	    std::cout << std::endl;
	 }
	cout<<poly.is_valid()<<endl;

	return poly;


}

int main(void){
	
	CGAL::Geomview_stream gv(CGAL::Bbox_3(-20,-20,-20,20,20,20));
	gv.clear();
	gv.set_line_width(5);
	gv.set_bg_color(CGAL::Color(0, 200, 200));
	gv.set_wired(true);
	


	std::vector<Point_3> nodes = getNodes(gv);
	unordered_map<int,int> nodeMap ;
	std::vector<Tetrahedron_3> tetrahedra = getTetrahedra(nodes,gv,nodeMap);
	std::vector< std::vector<int> > edges = getEdges(gv,nodes.size());
	

	std::vector<bool> visited(tetrahedra.size());
	for (std::vector<bool>::iterator i = visited.begin(); i != visited.end(); ++i){
		*i = false;
	}
	gv << CGAL::RED;



	Point_3 start = nodes[53];
	Tetrahedron_3 start_T = tetrahedra[nodeMap[53]];
	cout << "mapped to " << nodeMap[53] << endl;

	Point_3 end = nodes[52];


	gv << start;

	DepthFirst(nodes,start_T,end,tetrahedra,edges,nodeMap,visited,gv);


	Polyhedron P = PDB2Poly("g+","t","g-");
	gv << CGAL::BLUE;
	gv << P;

	char ch;
	std::cout << "Enter a key to finish" << std::endl;
	std::cin >> ch;







	return 0;
}



