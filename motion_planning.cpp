//CGAL Initialization
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Point_2.h>


#include <CGAL/Triangulation_2.h>



//Random Point Set Creation
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */


//Geomview 
#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Polyhedron_geomview_ostream.h>
#include <CGAL/IO/Triangulation_geomview_ostream_2.h>
#include <CGAL/IO/Triangulation_geomview_ostream_3.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <cassert>

#include <CGAL/intersections.h>
#include <CGAL/Triangle_2.h>


typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

typedef CGAL::Projection_traits_xy_3<Kernel>  Gt;


typedef CGAL::Delaunay_triangulation_2<Kernel>  Delaunay;
typedef Delaunay::Point Point;

typedef Kernel::Point_3 Point3;

typedef Kernel::Point_2 Point_2;

typedef Kernel::Triangle_2 Triangle;
typedef std::vector<Point_2> Vector;

typedef Kernel::Intersect_2 Intersect_2;

Triangle FaceToTriangle(Delaunay::Face_handle& face){
	Point_2 a,b,c;
	a = face->vertex(0)->point();
	b = face->vertex(1)->point();
	c = face->vertex(2)->point();
	return Triangle(a,b,c);
}

void PrintFace(Delaunay::Face_handle& face){
	Point_2 a,b,c;
	a = face->vertex(0)->point();
	b = face->vertex(1)->point();
	c = face->vertex(2)->point();
	std::cout << "\t" <<  a  << std::endl;
	std::cout << "\t" <<  b  << std::endl;
	std::cout << "\t" <<  c << std::endl;

}

int main(int argc , char* argv[]){




	unsigned int box_size = (argc>2) && (strcmp(argv[1],"-size")==0) ? atoi(argv[2]) : 30;
	unsigned int num_obstacles = (argc>2) && (strcmp(argv[3],"-obstacles")==0) ? atoi(argv[4]) : box_size/6;


	std::cout << "Box :" << box_size << "x"<< box_size << std::endl;
	std::cout<<"Obstacles = " << num_obstacles << std::endl;
	

	Vector Terrain;
	Terrain.reserve(box_size*box_size);
	for(int i = 0 ; i < box_size ; i++){
		for(int j = 0 ; j < box_size ; j++){
			Point_2 p(i,j);
			std::cout<<"Push "<< p <<std::endl;
			Terrain.push_back(p);

		}
	}



	CGAL::Geomview_stream gv(CGAL::Bbox_3(0,0, 0,box_size,box_size,0));

	Delaunay T;
	T.insert( Terrain.begin(),Terrain.end() );


	
	gv.set_line_width(10);
	gv.set_bg_color(CGAL::Color(0, 200, 200));
	std::cout << "Drawing 2D Delaunay triangulation in wired mode.\n";
	gv.set_wired(true);
	gv << CGAL::BLACK;
	gv << T;

	std::cout << "OK" << std::endl;
	

	Point_2 StartPoint(0,0);
	Point_2 EndPoint(box_size-1,box_size-1);

	std::cout << "Locate starting point..." << std::endl;
	Delaunay::Face_handle start_face = T.locate(StartPoint,T.all_faces_begin());
	std::cout<<"Starting Triangle"<<std::endl;
	Triangle start_triangle = FaceToTriangle(start_face);
	PrintFace(start_face);

	std::cout << "Locate ending point..." << std::endl;
	Delaunay::Face_handle end_face = T.locate(EndPoint,T.all_faces_begin());
	std::cout<<"Ending Triangle"<<std::endl;
	Triangle end_triangle = FaceToTriangle(end_face);
	PrintFace(end_face);


	std::cout << "Place Obstacles .." << std::endl;
	Vector Obstacles;
	Obstacles.reserve(num_obstacles);


	srand (time(NULL));
	int radius =  box_size;

	bool obstacle_intersect = true;
	for (int i = 0; i < num_obstacles ; i++) {
		
		int trys = 0;

		Delaunay::Face_handle obstacle_face;
		bool obstacle_intersect = true;
		while(obstacle_intersect){
			
			if(trys == 3){
				std::cout<<"Grid size too small or too many obstacles:= unable to place obstacle #"<<i<<std::endl;
				return -1;
			}
			
			Point p(rand()%box_size,rand()%box_size);		
			Obstacles.push_back(p);
			std::cout <<  "\t " << p <<std::endl;
			obstacle_face = T.locate(p,T.all_faces_begin());
			std::cout<<"Obstacle #"<<i<<std::endl;
			Triangle obstacle = FaceToTriangle(obstacle_face);

			CGAL::cpp11::result_of<Intersect_2(Kernel::Triangle_2,Kernel::Triangle_2)>::type
		    result1 = CGAL::intersection(obstacle,start_triangle);

		    CGAL::cpp11::result_of<Intersect_2(Kernel::Triangle_2,Kernel::Triangle_2)>::type
			result2 = CGAL::intersection(obstacle,end_triangle);		

			if((!result1)&& (!result2)){
				obstacle_intersect = false;
				PrintFace(obstacle_face);
			}
			trys++;
		}

		

	} 



	assert(Obstacles.size() == num_obstacles);






	std::cout << "Enter a key to finish" << std::endl;
	char ch;
	std::cin >> ch;

	return 0;
}





