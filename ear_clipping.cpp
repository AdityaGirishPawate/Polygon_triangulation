/*
 * This is a simple polygon triangulation algorithm developed by me.
 * This is submitted as a part of assignment 2 of oosd lab iit kgp
 * My name is Aditya Girish Pawate and Roll Number is 18MA20054
 * The following code is a simple implementation of polygon triangulation in order(n^2) using ear clipping method
 * An ear of a polygon is a triangle formed by three consecutive vertices Vi0, Vi1, and Vi2 for which Vi1 is a
    convex vertex (the interior angle at the vertex is smaller than π radians), the line segment from Vi0 to Vi2
    lies completely inside the polygon, and no vertices of the polygon are contained in the triangle other than
    the three vertices of the triangle. In the computational geometry jargon, the line segment is_in_between Vi0 and Vi2 
    is a diagonal of the polygon. The vertex Vi1 is called the ear tip. A triangle consists of a single ear,
    although you can place the ear tip at any of the three vertices. A polygon of four or more sides always has at
    least two non overlapping ears. This suggests a recursive approach to the triangulation. If you can locate
    an ear in a polygon with n ≥ 4 vertices and remove it, you have a polygon of n − 1 vertices and can repeat
    the process. A straightforward implementation of this will lead to an O(n^3) algorithm.
    With some careful attention to details, the ear clipping can be done in O(n^2) time. The first step is to
    store the polygon as a doubly linked list so that you can quickly remove ear tips. Construction of this list
    is an O(n) process. The second step is to iterate over the vertices and find the ears. For each vertex Vi
    and corresponding triangle hVi−1, Vi, Vi+1i (indexing is modulo n, so Vn = V0 and V−1 = Vn−1), test
    all other vertices to see if any are inside the triangle. If none are inside, you have an ear. If at least one
    is inside, you do not have an ear. The actual implementation I provide tries to make this somewhat more
    efficient. It is sufficient to consider only reflex vertices in the triangle containment test. A reflex vertex is one
    for which the interior angle formed by the two edges sharing it is larger than 180 degrees. A convex vertex
    is one for which the interior angle is smaller than 180 degrees. The data structure for the polygon maintains
    four doubly linked lists simultaneously, using an array for storage rather than dynamically allocating and
    deallocating memory in a standard list data structure. The vertices of the polygon are stored in a cyclical
    list, the convex vertices are stored in a linear list, the reflex vertices are stored in a linear list, and the ear
    tips are stored in a cyclical list.
 *
*/
#include <cmath>
#include <iostream>
using namespace std;

/* Global variable definitions */
typedef	double Point_Double[2];   /* Type double point */
typedef	struct Vertex_Structure vertex_structure_t;
typedef	vertex_structure_t *vertex_t;
//this is a structure of a vertex
struct	Vertex_Structure {
    int		vertex_number;		/* Index */
    Point_Double	v;		/* Coordinates */
    bool 	ear;		/* true iff an ear */
    vertex_t 	next,prev;
};
vertex_t	vertices  = nullptr;	/* "Head" of circular list. */
int	number_of_vertices = 0;		/* Total number of polygon vertices. */
typedef struct {
    double x, y;
} point_coordinates;
point_coordinates ** output;
int size=0;


class polygon_triangulation{
public:
    static void	Triangulate( );
    static void	ear_clipping_initiate( );
    static bool	is_proper_diagonal( vertex_t a, vertex_t b );
    static bool	is_diagonal( vertex_t a, vertex_t b );
    static bool	is_internal_diagonal( vertex_t a, vertex_t b );
    static bool	Xor( bool x, bool y );
    static bool	Left( Point_Double a, Point_Double b, Point_Double c );
    static bool	LeftOn( Point_Double a, Point_Double b, Point_Double c );
    static bool	is_collinear( Point_Double a, Point_Double b, Point_Double c );
    static bool	is_in_between( Point_Double a, Point_Double b, Point_Double c );
    static bool	does_intersect( Point_Double a, Point_Double b, Point_Double c, Point_Double d );
    static bool	does_intersect_properly( Point_Double a, Point_Double b, Point_Double c, Point_Double d );
    static vertex_t create_new_vertex( );
    static int sign_of_area( const Point_Double a, const Point_Double b, const Point_Double c );
};

/*---------------------------------------------------------------------
Exclusive or: true iff exactly one argument is true.
---------------------------------------------------------------------*/
bool polygon_triangulation::Xor( bool x, bool y )
{
    /* The arguments are negated to ensure that they are 0/1 values. */
    return   !x ^ !y;
}

/*---------------------------------------------------------------------
Returns true iff ab properly intersects cd: they share
a point interior to both segments.  The properness of the
intersection is ensured by using strict leftness.
---------------------------------------------------------------------*/
bool polygon_triangulation::does_intersect_properly( Point_Double a, Point_Double b, Point_Double c, Point_Double d )
{
    /* Eliminate improper cases. */
    if (
            is_collinear(a,b,c) ||
            is_collinear(a,b,d) ||
            is_collinear(c,d,a) ||
            is_collinear(c,d,b)
            )
        return false;

    return
            Xor( Left(a,b,c), Left(a,b,d) )
            && Xor( Left(c,d,a), Left(c,d,b) );
}

/*---------------------------------------------------------------------
Returns true iff c is strictly to the left of the directed
line through a to b.
---------------------------------------------------------------------*/
bool	polygon_triangulation::Left( Point_Double a, Point_Double b, Point_Double c )
{
    return  sign_of_area( a, b, c ) > 0;
}

bool	polygon_triangulation::LeftOn( Point_Double a, Point_Double b, Point_Double c )
{
    return  sign_of_area( a, b, c ) >= 0;
}

bool	polygon_triangulation::is_collinear( Point_Double a, Point_Double b, Point_Double c )
{
    return  sign_of_area( a, b, c ) == 0;
}
int     polygon_triangulation::sign_of_area( const Point_Double a, const Point_Double b, const Point_Double c )
{
    double area2;

    area2 = ( b[0] - a[0] ) * (double)( c[1] - a[1] ) -
            ( c[0] - a[0] ) * (double)( b[1] - a[1] );

    /* The area should be an integer. */
    if      ( area2 >  0.5 ) return  1;
    else if ( area2 < -0.5 ) return -1;
    else                     return  0;
}
/*---------------------------------------------------------------------
Returns true iff point c lies on the closed segement ab.
First checks that c is is_collinear with a and b.
---------------------------------------------------------------------*/
bool	polygon_triangulation::is_in_between( Point_Double a, Point_Double b, Point_Double c )
{
    if ( ! is_collinear( a, b, c ) )
        return  false;

    /* If ab not vertical, check is_in_betweenness on x; else on y. */
    if ( a[0] != b[0] )
        return ((a[0] <= c[0]) && (c[0] <= b[0])) ||
               ((a[0] >= c[0]) && (c[0] >= b[0]));
    else
        return ((a[1] <= c[1]) && (c[1] <= b[1])) ||
               ((a[1] >= c[1]) && (c[1] >= b[1]));
}

/*---------------------------------------------------------------------
Returns true iff segments ab and cd intersect, properly or improperly.
---------------------------------------------------------------------*/
bool	polygon_triangulation::does_intersect( Point_Double a, Point_Double b, Point_Double c, Point_Double d ){
    if     ( does_intersect_properly( a, b, c, d )
             || is_in_between( a, b, c )
             || is_in_between( a, b, d )
             || is_in_between( c, d, a )
             || is_in_between( c, d, b )
            )
        return  true;
    else    return  false;
}

/*---------------------------------------------------------------------
Returns true iff (a,b) is a proper internal *or* external
diagonal of P, *ignoring edges incident to a and b*.
---------------------------------------------------------------------*/
bool   polygon_triangulation::is_diagonal( vertex_t a, vertex_t b ){
    vertex_t c, c1;

    /* For each edge (c,c1) of P */
    c = vertices;
    do {
        c1 = c->next;
        /* Skip edges incident to a or b */
        if (    ( c != a ) && ( c1 != a )
                && ( c != b ) && ( c1 != b )
                && does_intersect( a->v, b->v, c->v, c1->v )
                )
            return false;
        c = c->next;
    } while ( c != vertices );
    return true;
}
/*---------------------------------------------------------------------
This function initializes the data structures, and calls
Triangulate2 to clip off the ears one by one.
---------------------------------------------------------------------*/
void   polygon_triangulation::ear_clipping_initiate(){
    vertex_t v0, v1, v2;   /* three consecutive vertices */

    /* Initialize v1->ear for all vertices. */
    v1 = vertices;
    do {
        v2 = v1->next;
        v0 = v1->prev;
        v1->ear = is_proper_diagonal( v0, v2 );
        v1 = v1->next;
    } while ( v1 != vertices );

}
/*---------------------------------------------------------------------
Prints out n-3 diagonals (as pairs of integer indices)
which form a triangulation of P.
---------------------------------------------------------------------*/
void   polygon_triangulation::Triangulate(){
    vertex_t v0, v1, v2, v3, v4;	/* five consecutive vertices */
    int   n = number_of_vertices;		/* number of vertices; shrinks to 3. */
    bool earfound;		/* for debugging and error detection only. */

    ear_clipping_initiate();
    /* Each step of outer loop removes one ear. */
    while ( n > 3 ) {
        /* Inner loop searches for an ear. */
        v2 = vertices;
        earfound = false;
        do {
            if (v2->ear) {
                earfound = true;
                /* Ear found. Fill variables. */
                v3 = v2->next; v4 = v3->next;
                v1 = v2->prev; v0 = v1->prev;

                /* (v1,v3) is a diagonal */
                output[size]->x = v1->vertex_number;
                output[size++]->y = v3->vertex_number;

                /* Update earity of diagonal endpoints */
                v1->ear = is_proper_diagonal( v0, v3 );
                v3->ear = is_proper_diagonal( v1, v4 );

                /* Cut off the ear v2 */
                v1->next = v3;
                v3->prev = v1;
                vertices = v3;	/* In case the head was v2. */
                n--;
                break;   /* out of inner loop; resume outer loop */
            } /* end if ear found */
            v2 = v2->next;
        } while ( v2 != vertices );

        if ( !earfound ) {
            cout<<"%%Error in Triangulate:  No ear found."<<endl;
            exit(1);
        }
    } /* end outer while loop */
}
/*---------------------------------------------------------------------
Returns true iff the diagonal (a,b) is strictly internal to the
polygon in the neighborhood of the a endpoint.
---------------------------------------------------------------------*/
bool   polygon_triangulation::is_internal_diagonal( vertex_t a, vertex_t b ){
    vertex_t a0,a1;	/* a0,a,a1 are consecutive vertices. */

    a1 = a->next;
    a0 = a->prev;

    /* If a is a convex vertex ... */
    if( LeftOn( a->v, a1->v, a0->v ) )
        return    Left( a->v, b->v, a0->v )
                  && Left( b->v, a->v, a1->v );

    /* Else a is reflex: */
    return !(    LeftOn( a->v, b->v, a1->v )
                 && LeftOn( b->v, a->v, a0->v ) );
}
/*---------------------------------------------------------------------
Returns true iff (a,b) is a proper internal diagonal.
---------------------------------------------------------------------*/
bool	polygon_triangulation::is_proper_diagonal( vertex_t a, vertex_t b ){
    return is_internal_diagonal( a, b ) && is_internal_diagonal( b, a ) && is_diagonal( a, b );
}
/*---------------------------------------------------------------------
create_new_vertex: Makes a vertex.
---------------------------------------------------------------------*/
vertex_t   polygon_triangulation::create_new_vertex( ){
    vertex_t  v;
    if ((v=(vertex_structure_t *) malloc (sizeof(vertex_structure_t))) == nullptr) {
        cout<<"NEW: Out of Memory!"<<endl;
        exit(1);
    }
    if ( vertices )  {
        v->next = vertices;
        v->prev = vertices->prev;
        vertices->prev = v;
        v->prev->next = v;
    }
    else {
        vertices = v;
        vertices->next = vertices->prev = v;
    }
    return v;
}


//this is main function where all the routines are called
int main()
{
    //object of class declared
    polygon_triangulation obj;
    vertex_t	v;
    int	vertex_number = 0;
    //this is the number of vertices
    int n=9; // CHANGE THIS TO GIVE DIFFERENT INPUTS
    //this is the data taken in clockwise order as required by the assignment
    double data[9][2] = {{0,0},  //CHANGE THESE VALUES TO GIVE DIFFERENT INPUTS // GIVE INPUT IN CLOCKWISE ORDER
                         {0,1},
                         {1,1},
                         {1,2},
                         {2,2},
                         {2,3},
                         {3,3},
                         {3,1},
                         {2,0}};
    int count=n;
    number_of_vertices=n;
    output = (point_coordinates **) malloc((n-3)*sizeof(point_coordinates **));
    for (int j = 0; j < (n-3); ++j) {
        output[j] = new point_coordinates;
    }
    //make a null vertex
    v = obj.create_new_vertex();
    // put vertices into the null vertex
    v->v[0] = data[0][0];
    v->v[1] = data[0][1];
    // increase the vertex number
    v->vertex_number = vertex_number++;
    while (  --count ){
        v = obj.create_new_vertex();
        v->v[0] = data[count][0];
        v->v[1] = data[count][1];
        v->vertex_number = vertex_number++;
    }
    //call the triangulate method through the object
    obj.Triangulate();
    //format output as required by the assignment
    std::cout<<"{";
    for(int k=0;k<(n-3);k++){
        if(output[k]->x !=output[k]->y ) {
            std::cout << "{";
            if (output[k]->x == 0)
                std::cout << output[k]->x;
            else
                std::cout << n - output[k]->x;
            std::cout << ",";
            if (output[k]->y == 0)
                std::cout << output[k]->y;
            else
                std::cout << n - output[k]->y;
            std::cout << "},";
        }
    }
    std::cout<<"}";
    return 0;
}
