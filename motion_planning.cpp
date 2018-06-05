//CGAL Initialization
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Point_2.h>

//Random Point Set Creation
#include <CGAL/point_generators_2.h>
#include <CGAL/Random.h>


//Geomview 
#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Polyhedron_geomview_ostream.h>
#include <CGAL/IO/Triangulation_geomview_ostream_2.h>
#include <CGAL/IO/Triangulation_geomview_ostream_3.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <cassert>




typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

typedef CGAL::Projection_traits_xy_3<Kernel>  Gt;


typedef CGAL::Delaunay_triangulation_2<Kernel>  Delaunay;
typedef Delaunay::Point Point;

typedef Kernel::Point_3 Point3;

typedef Kernel::Point_2 Point_2;
typedef std::vector<Point_2> Vector;
typedef std::vector<Point_2>::iterator viterator;

int main(int argc , char* argv[]){



	unsigned int box_size = (argc>2) && (strcmp(argv[1],"-size")==0) ? atoi(argv[2]) : 30;
	unsigned int num_obstacles = (argc>2) && (strcmp(argv[1],"-obstacles")==0) ? atoi(argv[2]) : box_size/6;


	std::cout << "Box :" << box_size << "x"<< box_size << std::endl;
	std::cout<<"Obstacles = " << num_obstacles << std::endl;
	
/*
	Vector Obstacles;
	Obstacles.reserve(num_obstacles);


	double radius =  box_size;
	CGAL::Random_points_in_square_2<Point> gen (radius);
	for (int i = 0; i < num_obstacles ; i++) {
		Point p = *gen++;		
		Obstacles.push_back(p);
		std::cout <<  "\t " << p <<std::endl;
	} 


	assert(Obstacles.size() == num_obstacles);
*/

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

	gv << CGAL::GREEN;

  	Delaunay::Finite_vertices_iterator vit;
  	/*for (vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit){
  		;
    }
    */
	gv << CGAL::RED;
	gv.set_line_width(4);
	gv.set_bg_color(CGAL::Color(0, 200, 200));
	std::cout << "Drawing 2D Delaunay triangulation in wired mode.\n";
	gv.set_wired(true);
	gv << T;

	std::cout << "OK" << std::endl;

	//CGAL_assertion( T.number_of_vertices() == 6 );

	

	

	std::cout << "Enter a key to finish" << std::endl;
	char ch;
	std::cin >> ch;

	return 0;
}

