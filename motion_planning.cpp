
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Triangle_2.h>
#include <iostream>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/Random.h>
#include <CGAL/algorithm.h>

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
typedef CGAL::Exact_predicates_exact_constructions_kernel         K;
typedef CGAL::Triangulation_vertex_base_2<K>                      Vb;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2,K>    Fbb;
typedef CGAL::Constrained_triangulation_face_base_2<K,Fbb>        Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>               TDS;
typedef CGAL::Exact_predicates_tag                                Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>  CDT;
typedef CGAL::Constrained_triangulation_plus_2<CDT> 			  CDTplus;
typedef CDT::Point                                                Point;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef K::Point_2 Point_2;
typedef K::Triangle_2 Triangle;
typedef std::vector<Point_2> Vector;
typedef K::Intersect_2 Intersect_2;
typedef K::Point_3 Point_3;



typedef std::list<Point_2>                                Container;
typedef CGAL::Polygon_2<K, Container>                     Polygon_2;
typedef CGAL::Creator_uniform_2<int, Point_2>             Creator;
typedef CGAL::Random_points_in_square_2<Point_2, Creator> Point_generator;

Triangle FaceToTriangle(CDTplus::Face_handle& face){
	Point_2 a,b,c;
	a = face->vertex(0)->point();
	b = face->vertex(1)->point();
	c = face->vertex(2)->point();
	return Triangle(a,b,c);
}

void PrintFace(CDTplus::Face_handle& face){
	Point_2 a,b,c;
	a = face->vertex(0)->point();
	b = face->vertex(1)->point();
	c = face->vertex(2)->point();
	std::cout << "\t" <<  a  << std::endl;
	std::cout << "\t" <<  b  << std::endl;
	std::cout << "\t" <<  c << std::endl;
}




void  mark_domains(CDTplus& ct, 
             CDTplus::Face_handle start, 
             int index, 
             std::list<CDTplus::Edge>& border,
             CGAL::Geomview_stream& lvlgv )
{
  if(start->info().nesting_level != -1){
    return;
  }
  std::list<CDT::Face_handle> queue;
  queue.push_back(start);
  while(! queue.empty()){
    CDTplus::Face_handle fh = queue.front();
    queue.pop_front();
    if(fh->info().nesting_level == -1){
      fh->info().nesting_level = index;
      std::cout<<"INDEX "<<index<<std::endl;
      Triangle t = FaceToTriangle(fh);
      if(fh->info().in_domain()){
      	lvlgv << CGAL::BLUE;
      	lvlgv << t;
      }
      else{
      	lvlgv << CGAL::RED;
      	lvlgv << t;
      }      
      

      for(int i = 0; i < 3; i++){
        CDTplus::Edge e(fh,i);
        CDTplus::Face_handle n = fh->neighbor(i);
        if(n->info().nesting_level == -1){
          if(ct.is_constrained(e)) border.push_back(e);
          else queue.push_back(n);
        }
      }
    }
  }
}
//explore set of facets connected with non constrained edges,
//and attribute to each such set a nesting level.
//We start from facets incident to the infinite vertex, with a nesting
//level of 0. Then we recursively consider the non-explored facets incident 
//to constrained edges bounding the former set and increase the nesting level by 1.
//Facets in the domain are those with an odd nesting level.
void mark_domains(CDTplus& cdt,CGAL::Geomview_stream& lvlgv){
 
  for(CDTplus::All_faces_iterator it = cdt.all_faces_begin(); it != cdt.all_faces_end(); ++it){
    it->info().nesting_level = -1;
  }
 
  std::list<CDTplus::Edge> border;
  mark_domains(cdt, cdt.infinite_face(), 0, border,lvlgv);
 
  while(! border.empty()){
    CDTplus::Edge e = border.front();
    border.pop_front();
    CDTplus::Face_handle n = e.first->neighbor(e.second);
    if(n->info().nesting_level == -1){
      mark_domains(cdt, n, e.first->info().nesting_level+1, border,lvlgv);
    }
 
  }
}






void FindPath(CDTplus& ct, CDTplus::Face_handle start, Point_2 EndPoint,CGAL::Geomview_stream& gv,int sleep_time){

	std::list<CDTplus::Face_handle> queue;
	std::list<CDTplus::Face_handle> path;

	std::cout<<"START LEVEL : " << start->info().nesting_level<<std::endl;
	path.push_back(start);
	queue.push_back(start);
	std::cout<<"Push"<<std::endl;
	while(! queue.empty()){
		CDTplus::Face_handle fh = queue.front();
		queue.pop_front();
		std::cout<<"Pop" << std::endl;
		fh->info().visited = true;
		std::cout<<"\tVisit face"<<std::endl;
		bool dead_end = true;
		for(int i = 0; i < 3; i++){
			CDTplus::Face_handle n = fh->neighbor(i);
			
			//CDTplus::Edge e(fh,i);
			//std::cout<<"EDGE"<<std::endl;
			//std::cout << e->first() << std::endl;
			//PrintFace(e.first());
			if(n->info().visited){
				std::cout<<"\t\tneighbor visited"<<std::endl;
				continue;
			}
			std::cout<<"\tlevel "<<n->info().nesting_level<<std::endl;
			if( (n->info().nesting_level != -1)&&((n->info().nesting_level)%2) == 1){
				std::cout<<"\t\tneighbor available"<<std::endl;
				dead_end = false;	
				path.push_back(n);
				Triangle current = FaceToTriangle(n);				
				gv << CGAL::BLUE;
				gv << current;
				sleep(sleep_time);
				
				queue.push_back(n);
				std::cout<<"Push"<<std::endl;
				CGAL::cpp11::result_of<Intersect_2(Triangle,Point_2)>::type
			    result = CGAL::intersection(current,EndPoint);
				if(result){
					std::cout<<"PATH FOUND"<<std::endl;
					std::list<CDTplus::Face_handle>::iterator it;

					int step = 0;
					for( it = path.begin(); it != path.end(); it++){
						Triangle t = FaceToTriangle(*it);
						//gv << CGAL::BLUE;
						//gv << t;
						std::cout<<"Step "<<step++ <<std::endl;
						PrintFace(*it);										
						
						return;

					}
				}	
			}
			else if (n->info().nesting_level == 0 )
				std::cout<<"\t\tneighbor is the bounding box"<<std::endl;	
			else 
				std::cout<<"\t\tneighbor is an obstacle"<<std::endl;		
			Triangle m = FaceToTriangle(n);
			gv << CGAL:: PURPLE;
			gv << m;
			sleep(sleep_time);
		}
		if(dead_end){

			std::cout<<"\tDead end\n"<<std::endl;
			//path.pop_back();			

		}
	}
	if(queue.empty()){
		std::cout<<"NO VIABLE PATH FOUND" <<std::endl;
	}
	else{
		std::cout<<"Uknown Path"<<std::endl;
	}
	return;

}


double RADIUS;
unsigned int MAX_POLY_SIZE;

int main(int argc,char* argv[]){


	unsigned int box_size = (argc>2) && (strcmp(argv[1],"-size")==0) ? atoi(argv[2]) : 30;
	unsigned int num_obstacles = (argc>4) && (strcmp(argv[3],"-obstacles")==0) ? atoi(argv[4]) : box_size/6;
	int sleep_time = (argc>6) && (strcmp(argv[5],"-sleep")==0) ? atoi(argv[6]) : 0;
	RADIUS = box_size;
	MAX_POLY_SIZE = box_size/10;
	std::cout << "Box :" << box_size << "x"<< box_size << std::endl;
	std::cout<<"Obstacles = " << num_obstacles << std::endl;

	CGAL::Geomview_stream gv(CGAL::Bbox_3(0,0, 0,box_size-1,box_size-1,0));
	gv.clear();
	gv.set_line_width(10);
	gv.set_bg_color(CGAL::Color(0, 200, 200));
	gv.set_wired(true);



	CGAL::Geomview_stream lvlgv(CGAL::Bbox_3(0,0, 0,box_size-1,box_size-1,0));
	lvlgv.clear();
	lvlgv.set_line_width(10);
	lvlgv.set_bg_color(CGAL::Color(0, 200, 200));
	lvlgv.set_wired(false);

	Point_2 StartPoint(0,(box_size-1)/2);
	Point_2 EndPoint(box_size-1,(box_size-1)/2);
	gv << CGAL::RED;
	gv << StartPoint;
	gv << EndPoint;


	//Insert the polygons into a constrained triangulation
	


	Polygon_2 Bounding_box;
	Bounding_box.push_back(Point(0,0));
	Bounding_box.push_back(Point(0,box_size-1));
	Bounding_box.push_back(Point(box_size-1,0));
	Bounding_box.push_back(Point(box_size-1,box_size-1));


	CDTplus cdt;


	cdt.insert_constraint(Bounding_box.vertices_begin(), Bounding_box.vertices_end(), true);


	srand(time(NULL));
	std::vector<K::Iso_rectangle_2> obstacles_vector;

	//construct two non-intersecting nested polygons  

	bool obstacle_intersect = true;
	for (int i = 0; i < num_obstacles ; i++) {
		
		int tries = 0;
		obstacle_intersect = true;
		while(obstacle_intersect){
			
			if(tries == 10){
				std::cout<<"Grid size too small or too many obstacles:= unable to place obstacle #"<<i<<std::endl;
				return -1;
			}					
		   
		   Polygon_2 Obstacle;
		   std::list<Point_2>   initial_point_set;
		   std::list<Point_2>   positive_point_set;
		   CGAL::Random         rand;
		   //std::cerr << "Seed = " <<  rand.get_seed() << std::endl;
		   	int size = rand.get_int(4, MAX_POLY_SIZE);

	   		CGAL::copy_n_unique(Point_generator(RADIUS), size,
		                       std::back_inserter(initial_point_set));
			


			for(std::list<Point_2>::iterator it = initial_point_set.begin();it != initial_point_set.end() ; it++){
				if((it->hx() < 0)){
					if((it->hy() < 0)){
						positive_point_set.push_back(Point_2(-it->hx(),-it->hy()));
					}
					else{
						positive_point_set.push_back(Point_2(-it->hx(),it->hy()));
					}
				}
				else{
					if((it->hy() < 0)){
						positive_point_set.push_back(Point_2(it->hx(),-it->hy()));
					}
					else{
						positive_point_set.push_back(Point_2(it->hx(),it->hy()));
					}
				}
			}


		    CGAL::random_polygon_2(positive_point_set.size(), std::back_inserter(Obstacle),
		                          positive_point_set.begin());
						
	
			
			Point a(Obstacle.left_vertex()->hx(),Obstacle.left_vertex()->hy());		
			Point b(Obstacle.top_vertex()->hx(),Obstacle.top_vertex()->hy());
			Point c(Obstacle.right_vertex()->hx(),Obstacle.right_vertex()->hy());
			Point d(Obstacle.bottom_vertex()->hx(),Obstacle.bottom_vertex()->hy());
			

		   
			Point_3 p(Obstacle.left_vertex()->hx(),Obstacle.left_vertex()->hy(),0);
			Point_3 q(Obstacle.top_vertex()->hx(),Obstacle.top_vertex()->hy(),0);
			Point_3 r(Obstacle.right_vertex()->hx(),Obstacle.right_vertex()->hy(),0);
			Point_3 s(Obstacle.bottom_vertex()->hx(),Obstacle.bottom_vertex()->hy(),0);	

			Polyhedron P;
			P.make_tetrahedron( p, q, r, s);


			K::Iso_rectangle_2 rect_obstacle(a,b,c,d);

			//std::cout <<  "\t " << Obstacle <<std::endl;
			CGAL::cpp11::result_of<Intersect_2(Point_2,K::Iso_rectangle_2)>::type
		    result1 = CGAL::intersection(StartPoint,rect_obstacle);

			CGAL::cpp11::result_of<Intersect_2(Point_2,K::Iso_rectangle_2)>::type
		    result2 = CGAL::intersection(EndPoint,rect_obstacle);	


		    bool obstacles_overlap = false;
		   
		    for(int l = 0 ; l < obstacles_vector.size() ; l++){
		    	K::Iso_rectangle_2 prev_obstacle = obstacles_vector.at(l);

				CGAL::cpp11::result_of<Intersect_2(K::Iso_rectangle_2,K::Iso_rectangle_2)>::type
			    result = CGAL::intersection(prev_obstacle,rect_obstacle);	
			   // std::cout<<"check for obstacles overlapping fully"<<std::endl;
			    if(result){
			    	obstacles_overlap = true;
			    	break;
			    }
		    }
		    
			

			if((!result1)&& (!result2) && (!obstacles_overlap)){
				obstacle_intersect = false;
				cdt.insert_constraint(Obstacle.vertices_begin(), Obstacle.vertices_end(), true);
				obstacles_vector.push_back(rect_obstacle);
				std::cout<<"Placed obstacle after "<<tries<<" tries"<<std::endl;
				gv << CGAL::RED;
				gv << P;

			}
			tries++;
		}

	
	} 



	assert(Obstacles.size() == num_obstacles);


	//Mark facets that are inside the domain bounded by the polygon
	mark_domains(cdt,lvlgv);



/*	std::cout << "Enter a key to finish" << std::endl;
	char ch;
	std::cin >> ch;
*/


	
	int count=0;
	for (CDTplus::Finite_faces_iterator fit=cdt.finite_faces_begin(); fit!=cdt.finite_faces_end();++fit){
		if (fit->info().in_domain())
			count++;
	}

	std::cout << "There are " << count << " facets in the domain." << std::endl;


	for (CDTplus::Subconstraint_iterator scit = cdt.subconstraints_begin();
		scit != cdt.subconstraints_end();
		++scit)  ++count;
	std::cout << "The number of resulting constrained edges is  ";
	std::cout <<  count << std::endl;



	std::cout << "Drawing 2D Constrained Delaunay triangulation in wired mode.\n";
	
	gv << CGAL::RED;
	gv << cdt;

	std::cout << "OK" << std::endl;
	


	std::cout << "Locate starting point..." << std::endl;
	CDTplus::Face_handle start_face = cdt.locate(StartPoint,cdt.all_faces_begin());
	std::cout<<"Starting Triangle"<<std::endl;
	Triangle start_triangle = FaceToTriangle(start_face);
	lvlgv << CGAL::GREEN;
	lvlgv << start_triangle;
	PrintFace(start_face);

	std::cout << "Locate ending point..." << std::endl;
	CDTplus::Face_handle end_face = cdt.locate(EndPoint,cdt.all_faces_begin());
	std::cout<<"Ending Triangle"<<std::endl;
	Triangle end_triangle = FaceToTriangle(end_face);
	lvlgv << CGAL::ORANGE;
	lvlgv << end_triangle;
	PrintFace(end_face);

	FindPath(cdt,start_face,EndPoint,gv,sleep_time);

	std::cout << "Enter a key to finish" << std::endl;
	char ch;
	std::cin >> ch;


	return 0;

}