/*NIKOLAS GIALITSIS 111520140027 UPOLOGISTIKI GEWMETRIA */
/*

Compilation:
cgal_create_CMakeLists -s robot
cmake .
make

Execution
./robot -size 100 -obstacles 2 -sleep 3 -freedom 0 
*/


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
typedef K::Intersect_2 Intersect_2;
typedef K::Point_3 Point_3;
typedef K::Line_2 Line;
typedef K::Segment_2 Segment;



typedef std::list<Point>                                Container;
typedef CGAL::Polygon_2<K, Container>                     Polygon_2;
typedef K::Ray_2  Ray;
typedef CGAL::Creator_uniform_2<int, Point>             Creator;
typedef CGAL::Random_points_in_square_2<Point, Creator> Point_generator;

typedef CGAL::Delaunay_triangulation_2<K> Delaunay_triangulation_2;

using std::cout; using std::endl;


Triangle FaceToTriangle(CDT::Face_handle&);

void PrintFace(CDT::Face_handle&);

void  mark_domains(CDT& ct, 
             CDT::Face_handle , 
             int , 
             std::list<CDT::Edge>& ,
             CGAL::Geomview_stream& );
void mark_domains(CDT& ,CGAL::Geomview_stream&) ;

std::list<CDT::Face_handle> FindPath(
	CDT&,CDT::Face_handle,Point,Point,std::vector<K::Iso_rectangle_2>,
	CGAL::Geomview_stream&,int
);


Triangle FaceToTriangle(CDT::Face_handle& face){
	Point a,b,c;
	a = face->vertex(0)->point();
	b = face->vertex(1)->point();
	c = face->vertex(2)->point();
	return Triangle(a,b,c);
}

void PrintFace(CDT::Face_handle& face){
	Point a,b,c;
	a = face->vertex(0)->point();
	b = face->vertex(1)->point();
	c = face->vertex(2)->point();
	std::cout << "\t" <<  a  << std::endl;
	std::cout << "\t" <<  b  << std::endl;
	std::cout << "\t" <<  c << std::endl;
}




void  mark_domains(CDT& ct, 
             CDT::Face_handle start, 
             int index, 
             std::list<CDT::Edge>& border,
             CGAL::Geomview_stream& indomain
 			)
{
  if(start->info().nesting_level != -1){
    return;
  }
  std::list<CDT::Face_handle> stack;
  stack.push_back(start);
  while(! stack.empty()){
    CDT::Face_handle fh = stack.front();
    stack.pop_front();
    if(fh->info().nesting_level == -1){
      fh->info().nesting_level = index;
      //std::cout<<"INDEX "<<index<<std::endl;
      Triangle t = FaceToTriangle(fh);
      if(fh->info().nesting_level == 1){
      	indomain << CGAL::BLUE;
      	indomain << t;
      }  
      else if(fh->info().nesting_level > 1){
      	indomain << CGAL::RED;
      	indomain << t;
      }  
      

      for(int i = 0; i < 3; i++){
        CDT::Edge e(fh,i);
        CDT::Face_handle n = fh->neighbor(i);
        if(n->info().nesting_level == -1){
          if(ct.is_constrained(e)) border.push_back(e);
          else stack.push_back(n);
        }
      }
    }
    else{
    	if(fh->info().nesting_level > 1){
    		Triangle t = FaceToTriangle(fh);
	      	indomain << CGAL::RED;
	      	indomain << t;
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
void mark_domains(CDT& cdt,CGAL::Geomview_stream& indomain){
 
  for(CDT::All_faces_iterator it = cdt.all_faces_begin(); it != cdt.all_faces_end(); ++it){
    it->info().nesting_level = -1;
  }
 
  std::list<CDT::Edge> border;
  mark_domains(cdt, cdt.infinite_face(), 0, border,indomain);
  
  while(! border.empty()){
    CDT::Edge e = border.front();
    border.pop_front();
    CDT::Face_handle n = e.first->neighbor(e.second);
    if(n->info().nesting_level == -1){
      mark_domains(cdt, n, e.first->info().nesting_level+1, border,indomain);
    }
 
  }
}



void PaintObstacles(std::vector<Point> obstacles_vector,CGAL::Geomview_stream& gv){
	gv << CGAL::RED;
	
	for(std::vector<Point>::iterator obit = obstacles_vector.begin();
		obit != obstacles_vector.end() ; obit ++)	
	{

		Polyhedron P;
		Point_3 a(obit->x(),obit->y(),0);
		obit++;
		Point_3 b(obit->x(),obit->y(),0);
		obit++;
		Point_3 c(obit->x(),obit->y(),0);
		obit++;
		Point_3 d(obit->x(),obit->y(),0);		
		P.make_tetrahedron(a,b,c,d);
		gv << P;		

	}	
	

}


std::list<CDT::Face_handle> FindPath(CDT& ct, CDT::Face_handle start, 
	Point StartPoint,
	Point  EndPoint,
	std::vector<Point> obstacles_vector,
	CGAL::Geomview_stream& gv,int sleep_time){

	std::list<CDT::Face_handle> stack;
	std::list<CDT::Face_handle> path;
	stack.push_back(start);
	std::cout<<"Push"<<std::endl;
	double current_dist = 100000000.0;
	double best_dist = 100000000.0;

	std::list<CDT::Face_handle> best_path;	

	std::list<std::list<CDT::Face_handle> > output;
	bool path_exists = false;

	sleep(sleep_time);
	while(! stack.empty()){
		CDT::Face_handle fh = stack.back();
		stack.pop_back();
		if(!output.empty())path = output.back();
		path.push_back(fh);
		output.push_back(path);
	

		std::list<CDT::Face_handle>::iterator pit;
		std::list<CDT::Face_handle> temp;					
		temp = output.back();					
		

		gv << CGAL::GREEN;
		gv << StartPoint;
		gv << CGAL::RED;
		gv << EndPoint;
		PaintObstacles(obstacles_vector,gv);
		for( pit = temp.begin(); pit != temp.end(); pit++){
			Triangle t = FaceToTriangle(*pit);
			gv << CGAL::BLUE;
			gv << t;
			
			//PrintFace(*it);
		}
	
		sleep(sleep_time);
		gv.clear();		
		gv << ct;
		PaintObstacles(obstacles_vector,gv);
		



		cout<<"OUTPUT SIZE "<<output.size()<<endl<<endl<<endl;
		std::cout<<"Pop" << std::endl;
		if(fh->info().visited == true){
			cout<<"\tBacktracking..."<<endl;
			output.pop_back();
			while(true){
				if(output.empty() || stack.empty())break;
				CDT::Face_handle  c = output.back().back();
				bool reset = true;
				for(int x = 0;  x < 3; x++){
					CDT::Face_handle n = c ->neighbor(x);			
					if (n == stack.back()){
						reset = false;
						break;
					}


					
				}			
				if(!reset)break;
				else {
					output.pop_back();			
				}							
			}		
			continue;
		}

		fh->info().visited = true;
		std::cout<<"\tVisit face"<<std::endl;
		bool dead_end = true;
		int free_neighbors = 0;
		for(int i = 0; i < 3; i++){
			CDT::Face_handle n = fh->neighbor(i);
			
			if(n->info().visited){
				std::cout<<"\t\tneighbor visited"<<std::endl;
				continue;
			}
			Triangle current = FaceToTriangle(n);	
			if(n->info().nesting_level == 1 ){
				free_neighbors++;
				dead_end = false;
				stack.push_back(n);
				std::cout<<"Push"<<std::endl;
				CGAL::cpp11::result_of<Intersect_2(Triangle,Point)>::type
			    result = CGAL::intersection(current,EndPoint);
				if(result){
					stack.pop_back();

					std::cout<<"PATH FOUND"<<std::endl;
					path_exists = true;
					std::list<CDT::Face_handle>::iterator it;
					std::list<CDT::Face_handle> final_path;					
					final_path = output.back();					
					final_path.push_back(n);
					int step = 0;

					gv.clear();
					gv << CGAL::GREEN;
					gv << StartPoint;
					gv << CGAL::RED;
					gv << EndPoint;

					PaintObstacles(obstacles_vector,gv);
					gv << ct;
					gv << CGAL::YELLOW;
					current_dist = 0;
					for( it = final_path.begin(); it != final_path.end(); it++){
						Triangle t = FaceToTriangle(*it);
						gv << t;
						for(std::vector<Point>::iterator obit = obstacles_vector.begin();
							obit != obstacles_vector.end() ; obit++)	
						{
							std::vector<Point> obstacle_vertices;
							for(int p = 0 ; p < 4 ; p ++){
								obstacle_vertices.push_back(*obit);
							}						

							Point obstacle_center = CGAL::centroid(obstacle_vertices.begin(),obstacle_vertices.end());
							current_dist += CGAL::squared_distance(t,obstacle_center);
						}			
					}
					if(current_dist <= best_dist){
						best_dist = current_dist;
						best_path.clear();
						std::copy(final_path.begin(),final_path.end(),std::back_inserter(best_path));
					}	

					dead_end = true;
					sleep(sleep_time);
					gv.clear();
					gv << CGAL::GREEN;
					gv << StartPoint;
					gv << CGAL::RED;
					gv << EndPoint;
					gv << ct;
								
					break;
				}
				else{

					std::cout<<"\t\tneighbor is an obstacle"<<std::endl;	
				}		
			}
			
		}
		if(dead_end){
			std::cout<<"\tDead end\n"<<std::endl;
			output.pop_back();	
			while(true){
				if(output.empty() || stack.empty())break;
				CDT::Face_handle  c = output.back().back();
				bool reset = true;
				for(int x = 0;  x < 3; x++){
					CDT::Face_handle n = c ->neighbor(x);			
					if (n == stack.back()){
						reset = false;
						break;
					}


					
				}			
				if(!reset)break;
				else {
					output.pop_back();			
				}							
			}		
		}
	}

	if(!path_exists){
		std::cout<<"NO VIABLE PATH FOUND" <<std::endl;		
		std::list<CDT::Face_handle> empty_list;
		return empty_list;
	}
	else{
		std::cout<<"Printing Best Path"<<endl;
		std::list<CDT::Face_handle>::iterator it;
		gv.clear();
		gv << ct;
		gv << CGAL::GREEN;
		gv << StartPoint;
		gv << CGAL::RED;
		gv << EndPoint;
		for( it = best_path.begin(); it != best_path.end(); it++){
			Triangle t = FaceToTriangle(*it);
			gv << CGAL::ORANGE;
			gv << t;
			
		}
		return best_path;
	}
	

}






double RADIUS;
unsigned int MAX_POLY_SIZE;

int main(int argc,char* argv[]){


	unsigned int box_size =  30;
	unsigned int num_obstacles = 2;
	int sleep_time =  0;
	unsigned int freedom = 0 ;

	for(int i = 0 ; i < argc ; i++){
		if(strcmp(argv[i],"-size")==0)box_size = atoi(argv[i+1]);
		else if (strcmp(argv[i],"-obstacles")==0)num_obstacles = atoi(argv[i+1]);
		else if (strcmp(argv[i],"-sleep")==0)sleep_time = atoi(argv[i+1]);
		else if (strcmp(argv[i],"-freedom")==0) freedom = atoi(argv[i+1]);
		else continue;
		i++;
	}

	RADIUS = box_size-1;
	MAX_POLY_SIZE = box_size/2.0;
	std::cout << "Box :" << box_size << "x"<< box_size << std::endl;
	std::cout<<"Obstacles = " << num_obstacles << std::endl;
	std::cout<<"Delay  = " << sleep_time << std::endl;
	std::cout<<"Freedom = " << freedom << std::endl;


	CGAL::Geomview_stream gv(CGAL::Bbox_3(0,0, 0,box_size-1,box_size-1,0));
	gv.clear();
	gv.set_line_width(5);
	gv.set_bg_color(CGAL::Color(0, 200, 200));
	gv.set_wired(true);



	CGAL::Geomview_stream indomain(CGAL::Bbox_3(0,0, 0,box_size-1,box_size-1,0));
	indomain.clear();
	indomain.set_line_width(5);
	indomain.set_bg_color(CGAL::Color(0, 200, 200));
	indomain.set_wired(true);



	Point StartPoint(0,(box_size-1)/2);
	Point EndPoint(box_size-1,(box_size-1)/2);

	gv << CGAL::GREEN;
	gv << StartPoint;

	gv << CGAL::RED;	
	gv << EndPoint;


	//Insert the polygons into a constrained triangulation
	

	Polygon_2 Bounding_box;
	Bounding_box.push_back(Point(0,0));
	Bounding_box.push_back(Point(box_size-1,0));
	Bounding_box.push_back(Point(box_size-1,box_size-1));
	Bounding_box.push_back(Point(0,box_size-1));

	CDT cdt;
	cdt.insert_constraint(Bounding_box.vertices_begin(),Bounding_box.vertices_end(),true);

	std::cout<<"Inserted bounding box"<<std::endl;
	std::vector<Point> obstacles_vector;



	//Dimiourgia Empodiwn
	for (int i = 0; i < num_obstacles ; i++) {
	
			Polygon_2 Obstacle;
			std::list<Point>   initial_point_set;
			std::list<Point>   positive_point_set;
			CGAL::Random         rand;
			std::cerr << "Seed = " <<  rand.get_seed() << std::endl;

			//Dimiourgise tessera monadika simeia se aktina isi me ton xwro
			CGAL::copy_n_unique(Point_generator(RADIUS), 4,
			               std::back_inserter(initial_point_set));

			//Metetrepse ta arnitika simeia se thetika,
			//An ta simeia apexoun apostasi < freedom apo to bounding box
			//sprwksta prosta auto,gia na diskoleutei to robot

			for(std::list<Point>::iterator it = initial_point_set.begin();it != initial_point_set.end() ; it++){
				if((it->x() < 0)){
					if((it->y() < 0)){					
						if(-it->x() < freedom)
							positive_point_set.push_back(Point(0,-it->y()));
						else if(-it->y() < freedom)
							positive_point_set.push_back(Point(-it->x(),0));
						else if(-it->x() > box_size-1-freedom)
							positive_point_set.push_back(Point(box_size-1,-it->y()));
						else if(-it->y() > box_size-1-freedom)
							positive_point_set.push_back(Point(-it->x(),box_size-1));
						else
							positive_point_set.push_back(Point(-it->x(),-it->y()));						
					}
					else{
						if(-it->x() < freedom)
							positive_point_set.push_back(Point(0,it->y()));
						else if(it->y() < freedom)
							positive_point_set.push_back(Point(-it->x(),0));
						else if(-it->x() > box_size-1-freedom)
							positive_point_set.push_back(Point(box_size-1,it->y()));
						else if(it->y() > box_size-1-freedom)
							positive_point_set.push_back(Point(-it->x(),box_size-1));
						else
							positive_point_set.push_back(Point(-it->x(),it->y()));						
					}
				}
				else{
					if((it->y() < 0)){
						if(it->x() < freedom)
							positive_point_set.push_back(Point(0,-it->y()));
						else if(-it->y() < freedom)
							positive_point_set.push_back(Point(it->x(),0));
						else if(it->x() > box_size-1-freedom)
							positive_point_set.push_back(Point(box_size-1,-it->y()));
						else if(-it->y() > box_size-1-freedom)
							positive_point_set.push_back(Point(it->x(),box_size-1));
						else
							positive_point_set.push_back(Point(it->x(),-it->y()));		
					}
					else{
						if(it->x() < freedom)
							positive_point_set.push_back(Point(0,it->y()));
						else if(it->y() < freedom)
							positive_point_set.push_back(Point(it->x(),0));
						else if(it->x() > box_size-1-freedom)
							positive_point_set.push_back(Point(box_size-1,it->y()));
						else if(it->y() > box_size-1-freedom)
							positive_point_set.push_back(Point(it->x(),box_size-1));
						else
							positive_point_set.push_back(Point(it->x(),it->y()));	
					}
				}
			}


	
			//Dimiourgise to polugono me ta tuxaia simeia
		    CGAL::random_polygon_2(4, std::back_inserter(Obstacle),
		                          positive_point_set.begin());


		    assert(Obstacle.vertices_num() == 4);

		    //Elegxos prosanatolismou polugonou
			if(Obstacle.orientation() == CGAL::CLOCKWISE)
				Obstacle.reverse_orientation();
		    
		    if(Obstacle.orientation() != CGAL::COUNTERCLOCKWISE){
		    	cout << "Invalid orientation"<<std::endl;
		    	return -1;
		    }
			//Elegxos kurtotitas empodiou
  			std::cout << "The polygon is " <<
    		(Obstacle.is_convex() ? "" : "not ") << "convex." << std::endl;
			

			//Apomonose ta stoixeia gia na topothetithoun ston obstacles_vector
			std::_List_iterator<Point> iter = Obstacle.vertices_begin();
			
			Point a(iter->x(),iter->y());	
			cout<<"A : "<< a<< endl;
			iter++;	
			Point b(iter->x(),iter->y());	
			cout<<"B : "<< b<< endl;
			iter++;
			Point c(iter->x(),iter->y());	
			cout<<"C : "<< c<< endl;
			iter++;	
			Point d(iter->x(),iter->y());	
			cout<<"D : "<< d<< endl;

			cout << std::endl;
		   	
				
			//Eisigage to empodio stin trigwnopoihsh
			cdt.insert_constraint(Obstacle.vertices_begin(),Obstacle.vertices_end(),true);

			//Elegxe tin orthotita tis trigwnopoihshs
			if(!cdt.is_valid()){
				std::cout<<"CDT not valid"<<std::endl;
				return -2;
			}

			obstacles_vector.push_back(a);		
			obstacles_vector.push_back(b);
			obstacles_vector.push_back(c);
			obstacles_vector.push_back(d);	
		
	} 



	assert(Obstacles.size() == num_obstacles);

	indomain << cdt;

	//Arithmise ta faces tis trigwnopoihshs 
	//ksekinontas apo to infinite vertex kai epipedo 0(diladi eksw ap to bounding box)
	//etsi wste ta faces epipedou 1 na apoteloun ton eleuthero xwro
	//se kathe epanalipsi to epipedo auksanetai

	mark_domains(cdt,indomain);

	int count=0;
	for (CDT::Finite_faces_iterator fit=cdt.finite_faces_begin(); fit!=cdt.finite_faces_end();++fit){
		if (fit->info().in_domain())
			count++;
	}

	std::cout << "There are " << count << " facets in the domain." << std::endl;
	std::cout << "Drawing 2D Constrained Delaunay triangulation in wired mode.\n";
	
	
	gv << cdt;
	

	std::cout << "OK" << std::endl;


	//Entopise pou vrisketai to simeio ekkinisis kai to simeio stoxou
	//stin trigwnopoihsh

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
	 

	std::list<CDT::Face_handle> output;
	
	//Vres an uparxei,to kalutero monopati 
	output = FindPath(cdt,start_face,StartPoint,EndPoint,obstacles_vector,gv,sleep_time);



	sleep(2*sleep_time);
	gv.clear();
	gv << cdt;
	gv << CGAL::GREEN;
	PaintObstacles(obstacles_vector,gv);



	//Dimiourgise to euthigrammo monopati tis diadromis tou robot ws eksis:
	//Ksekinontas apo to simeio ekkinisis upologizetai to epomeno simeio
	//tis kinisis , to opoio einai to meso tis pleuras pou enwnei to arxiko 
	//me to epomeno trigwno tou monopatiou.
	if(!output.empty()){
		bool first_point = true;
		Point last_point;
		for(std::list<CDT::Face_handle>::iterator it = output.begin();it!=output.end();it++){
		
			std::vector<Point> points;
			Triangle t = FaceToTriangle(*it);
			Point p;
			if(first_point){
				first_point = false;
				last_point = StartPoint;
				it--;
			}
			else if(t == end_triangle){
				Segment connect(last_point,EndPoint);
				gv << CGAL::GREEN;
				gv << connect;
				break;
			}
			else{
				it++;
				if(it == output.end()){
					it--;
					break;
				}

				Triangle next_t = FaceToTriangle(*it);
				it--;		
			

				CGAL::cpp11::result_of<Intersect_2(Triangle,Triangle)>::type result = CGAL::intersection(t,next_t);
				if(!result){
					break;
					cout<<"No intersection"<<endl;
				}

				std::vector<Point>* lv = boost::get<std::vector<Point> >(&*result);		
				Point* lp = boost::get<Point>(&*result);
				Segment* ls = boost::get<Segment>(&*result);
				Point midpoint;
				if(lv)
					midpoint = CGAL::midpoint(lv->front(),lv->back());
				else if(lp){
					midpoint = *lp;
				}
				else if(ls){
					midpoint = CGAL::midpoint(ls->source(),ls->target());
				}
				else
					cout<<"Unexpected return value"<<endl;

				Segment connect(last_point,midpoint);			
				gv << CGAL::GREEN;
				gv << connect;
				last_point = midpoint;
			}

			
		}
	}
	gv << CGAL::BLUE;
	gv << StartPoint;
	gv<< CGAL::RED;
	gv << EndPoint;



	
	std::cout << "Enter a key to finish" << std::endl;
	char ch;
	std::cin >> ch;


	return 0;

}