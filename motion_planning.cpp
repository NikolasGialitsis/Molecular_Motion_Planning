
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Triangle_2.h>
#include <iostream>
#include <CGAL/Iso_rectangle_2.h>


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
  FaceInfo2(){}
  int nesting_level;
  bool in_domain(){ 
    return nesting_level%2 == 1;
  }
};
typedef CGAL::Exact_predicates_exact_constructions_kernel         K;
typedef CGAL::Triangulation_vertex_base_2<K>                      Vb;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2,K>    Fbb;
typedef CGAL::Constrained_triangulation_face_base_2<K,Fbb>        Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>               TDS;
typedef CGAL::Exact_predicates_tag                                Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>  CDT;
typedef CDT::Point                                                Point;
typedef CGAL::Polygon_2<K>                                        Polygon_2;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef K::Point_2 Point_2;
typedef K::Triangle_2 Triangle;
typedef std::vector<Point_2> Vector;
typedef K::Intersect_2 Intersect_2;
typedef K::Point_3 Point_3;

void  mark_domains(CDT& ct, 
             CDT::Face_handle start, 
             int index, 
             std::list<CDT::Edge>& border )
{
  if(start->info().nesting_level != -1){
    return;
  }
  std::list<CDT::Face_handle> queue;
  queue.push_back(start);
  while(! queue.empty()){
    CDT::Face_handle fh = queue.front();
    queue.pop_front();
    if(fh->info().nesting_level == -1){
      fh->info().nesting_level = index;
      for(int i = 0; i < 3; i++){
        CDT::Edge e(fh,i);
        CDT::Face_handle n = fh->neighbor(i);
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
void mark_domains(CDT& cdt){
 
  for(CDT::All_faces_iterator it = cdt.all_faces_begin(); it != cdt.all_faces_end(); ++it){
    it->info().nesting_level = -1;
  }
 
  std::list<CDT::Edge> border;
  mark_domains(cdt, cdt.infinite_face(), 0, border);
 
  while(! border.empty()){
    CDT::Edge e = border.front();
    border.pop_front();
    CDT::Face_handle n = e.first->neighbor(e.second);
    if(n->info().nesting_level == -1){
      mark_domains(cdt, n, e.first->info().nesting_level+1, border);
    }
 
  }
}


Triangle FaceToTriangle(CDT::Face_handle& face){
	Point_2 a,b,c;
	a = face->vertex(0)->point();
	b = face->vertex(1)->point();
	c = face->vertex(2)->point();
	return Triangle(a,b,c);
}

void PrintFace(CDT::Face_handle& face){
	Point_2 a,b,c;
	a = face->vertex(0)->point();
	b = face->vertex(1)->point();
	c = face->vertex(2)->point();
	std::cout << "\t" <<  a  << std::endl;
	std::cout << "\t" <<  b  << std::endl;
	std::cout << "\t" <<  c << std::endl;
}





void FindPath(CDT& ct, CDT::Face_handle start, Point_2 EndPoint,CGAL::Geomview_stream& gv){

	std::list<CDT::Face_handle> queue;
	std::list<CDT::Face_handle> path;

	path.push_back(start);
	queue.push_back(start);

	while(! queue.empty()){
		CDT::Face_handle fh = queue.front();
		queue.pop_front();
		std::cout<<"New face"<<std::endl;
		for(int i = 0; i < 3; i++){
			CDT::Edge e(fh,i);
			CDT::Face_handle n = fh->neighbor(i);
			bool dead_end = true;

			if( ((n->info().nesting_level)%2) == 1){

				dead_end = false;	
				path.push_back(n);
				Triangle current = FaceToTriangle(n);				
				gv << CGAL::BLUE;
				gv << current;
				queue.push_back(n);

				CGAL::cpp11::result_of<Intersect_2(Triangle,Point_2)>::type
			    result = CGAL::intersection(current,EndPoint);
				if(result){
					std::cout<<"PATH FOUND"<<std::endl;
					std::list<CDT::Face_handle>::iterator it;

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
			if(dead_end){
				std::cout<<"Dead end"<<std::endl;
				path.pop_back();
			}
			
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



int main(int argc,char* argv[]){


	unsigned int box_size = (argc>2) && (strcmp(argv[1],"-size")==0) ? atoi(argv[2]) : 30;
	unsigned int num_obstacles = (argc>2) && (strcmp(argv[3],"-obstacles")==0) ? atoi(argv[4]) : box_size/6;


	std::cout << "Box :" << box_size << "x"<< box_size << std::endl;
	std::cout<<"Obstacles = " << num_obstacles << std::endl;

	CGAL::Geomview_stream gv(CGAL::Bbox_3(0,0, 0,box_size-1,box_size-1,0));
	gv.clear();
	gv.set_line_width(10);
	gv.set_bg_color(CGAL::Color(0, 200, 200));
	gv.set_wired(true);


	Point_2 StartPoint(0,0);
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


	CDT cdt;
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
			

			int obstacle_area = 0;
			bool obstacle_big = true;
			
			
			

			int ax,ay,bx,by,cx,cy,dx,dy;
			int times_resized = 0;
			while(obstacle_big){


				Polygon_2 Obstacle1;
				ax = rand()%box_size;				
				bx = rand()%box_size;
				cx = rand()%box_size;
				dx = rand()%box_size;

				ay = rand()%box_size;
				by = rand()%box_size;
				cy = rand()%box_size;
				dy = rand()%box_size;

				Point a(ax,ay);		
				Point b(bx,by);
				Point c(cx,cy);
				Point d(dx,dy);
				

				
				Obstacle1.push_back(a);
				Obstacle1.push_back(b);
				Obstacle1.push_back(c);
				Obstacle1.push_back(d);

				Polygon_2::FT pol_area = Obstacle1.area();

				if(pol_area < 0) pol_area = - pol_area;
			
				if(( pol_area <= ((box_size*box_size)/50.0) )&& pol_area > 20.0){
					std::cout << "Obstacle #"<< i<< "  area after resizing "<< times_resized<< " times is "<< pol_area << std::endl;
					obstacle_big = false;
				}
				else 
					times_resized++;
			}
			
	
			Point a(ax,ay);		
			Point b(bx,by);
			Point c(cx,cy);
			Point d(dx,dy);

			Polygon_2 Obstacle;
			Obstacle.push_back(a);			
			Obstacle.push_back(b);
			Obstacle.push_back(c);
			Obstacle.push_back(d);


			

			Point_3 p( a.hx(), a.hy(), 0.0);
			Point_3 q( b.hx(), b.hy(), 0.0);
			Point_3 r( c.hx(), c.hy(), 0.0);
			Point_3 s( d.hx(), d.hy(), 0.0);

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
	mark_domains(cdt);

	int count=0;
	for (CDT::Finite_faces_iterator fit=cdt.finite_faces_begin(); fit!=cdt.finite_faces_end();++fit){
		if (fit->info().in_domain())
			count++;
	}

	std::cout << "There are " << count << " facets in the domain." << std::endl;


	
	std::cout << "Drawing 2D Constrained Delaunay triangulation in wired mode.\n";
	
	gv << CGAL::RED;
	gv << cdt;

	//gv << CGAL::BLUE;
	//gv << polygon1;
	//gv << CGAL::PURPLE;
	//gv << Bounding_box;

	std::cout << "OK" << std::endl;
	


	std::cout << "Locate starting point..." << std::endl;
	CDT::Face_handle start_face = cdt.locate(StartPoint,cdt.all_faces_begin());
	std::cout<<"Starting Triangle"<<std::endl;
	Triangle start_triangle = FaceToTriangle(start_face);
	gv << CGAL::GREEN;
	gv << start_triangle;
	PrintFace(start_face);

	std::cout << "Locate ending point..." << std::endl;
	CDT::Face_handle end_face = cdt.locate(EndPoint,cdt.all_faces_begin());
	std::cout<<"Ending Triangle"<<std::endl;
	Triangle end_triangle = FaceToTriangle(end_face);
	gv << CGAL::ORANGE;
	gv << end_triangle;
	PrintFace(end_face);




	FindPath(cdt,start_face,EndPoint,gv);

	std::cout << "Enter a key to finish" << std::endl;
	char ch;
	std::cin >> ch;


	return 0;

}