//#include "test.h"


//int mk_CGAL::test()
//{
	//cout << __FUNCTION__ << endl;
	//// construction from a list of points :
	//std::list<Point> L;
	//L.push_front(Point(0,0,0));
	//L.push_front(Point(1,0,0));
	//L.push_front(Point(0,1,0));
	//Triangulation T(L.begin(), L.end());
	//Triangulation::size_type n = T.number_of_vertices();
	//// insertion from a vector :
	//std::vector<Point> V(3);
	//V[0] = Point(0,0,1);
	//V[1] = Point(1,1,1);
	//V[2] = Point(2,2,2);
	//n = n + T.insert(V.begin(), V.end());
	//assert( n == 6 ); // 6 points have been inserted
	//assert( T.is_valid() ); // checking validity of T
	//Locate_type lt;
	//int li, lj;
	//Point p(0,0,0);
	//Cell_handle c = T.locate(p, lt, li, lj);
	//// p is the vertex of c of index li :
	//assert( lt == Triangulation::VERTEX );
	//assert( c->vertex(li)->point() == p );
	//Vertex_handle v = c->vertex( (li+1)&3 );
	//// v is another vertex of c
	//Cell_handle nc = c->neighbor(li);
	//// nc = neighbor of c opposite to the vertex associated with p
	//// nc must have vertex v :
	//int nli;
	//assert( nc->has_vertex( v, nli ) );
	//// nli is the index of v in nc
	//
	//std::ofstream oFileT("output",std::ios::out);
	//// writing file output;	
	//oFileT << T;

	//Polyhedron P;
	//std::cin >> P;
	//CGAL::VRML_1_ostream out(std::cout);
	//out << P;
	////return ( std::cin && std::cout) ? 0 : 1;

	//Triangulation T1;
	//std::ifstream iFileT("output",std::ios::in);
	//// reading file output;
	//iFileT >> T1;
	//assert( T1.is_valid() );
	//assert( T1.number_of_vertices() == T.number_of_vertices() );
	//assert( T1.number_of_cells() == T.number_of_cells() );

	//return 0;
//}