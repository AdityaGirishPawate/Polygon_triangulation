////////////////////////////////////////////////BEGIN////////////////////////////////////////////////////////////////////////
// This program is written as a submission to oosd lab assignment 2 on polygon trinagulation
// My ROLL number is 18MA20054, Name is Aditya Girish Pawate
// In this program I have used the Seidel's Algorithm to develop my program.
// I have developed the algorithm from scratch using c++
// Computing the triangulation of a polygon is a fundamental algorithm in computational geometry.
// I have first explained what the algorithm is and why it's order is O(Nlog(N))
/*  This is the complete description of the algorithm
 *  It is an incremental randomized algorithm whose expected complexity is O(n log*n). In practice, it is almost linear
 *  time for a simple polygon having n vertices.
 *  The triangulation does not introduce any additional vertices and decomposes the polygon into n-2 triangles.
 *  Seidel’s algorithm is a randomized way of computing a trapezoidation. The
    complexity analysis is based upon probabilistic expectancies about for instance the depth
    of a binary lookup tree. In order to use expectancies, we need to have some random
    variables. Nothing is assumed about the input, but it is randomized in the sense that every
    order of adding the edges is equally likely. At any stage it will then be equally likely for
    the next edge to be added to be any one of the remaining edges. In this way one can use
    expectancies to estimate the running time. This will not estimate the worst-case behavior,
    but it is unlikely that the actual running time will be significantly worse than the expected
    running time.
 *  Any trapezoidation routine starts with an empty plane, and adds in subset by subset
    of the polygon. In Seidel’s algorithm each subset consists of one edge only: The
    algorithm first checks to see if the endpoints of the edge have been added to the
    trapezoidation already. Those that have not been added yet are added one by one, and
    divide the regions they are added into in two by means of a horizontal ray extending in
    each direction from the point. Next the edge is added between the two points. This is
    done by starting from the top point, and moving down one region at a time until the
    bottom point is reached. For each region traversed in this manner, the region is divided in
    two by means of the new edge. Whenever this causes two regions on top of each other to
    have the same left and right hand boundary, these two regions are merged into one. Each
    of the regions created in this manner will be a trapezoid. Note that in doing this, we never
    need to explicitly calculate the intersections of edges and rays, it suffices for each
    trapezoid to know which edges and rays it is bounded by
 *  Since we treat points lexicographically, no two points will be considered to have the
    same height, i.e. y-coordinate, and therefore each trapezoid will have only one upper and
    one lower bounding vertex, and no more than two neighbors above or below. Further this
    means that no edge will be horizontal, and the meaning of “left of” an edge will be well
    defined. When the algorithm has split one region during the addition of an edge, it needs
    to find the correct region below to move into. If there is only one region, this is trivial, if
    there are two, then the correct one can be found by checking whether the endpoint of the
    edge being added is to the right or left of the edge that separates the two regions.
 *  The algorithm also needs to efficiently be able to tell where to place new points in
    the current trapezoidation. One way of doing this is by creating a tree lookup structure for
    the trapezoids. The tree starts with only a root representing the original empty plane, and
    each time a region is split in two, the corresponding leaf in the tree is changed into a node
    with two children. The node represents the edge or the vertex that split the trapezoid, the
    two children represent the two resulting trapezoids. In this way each leaf in the tree
    corresponds to a trapezoid of the trapezoidation and each node represents either a vertex
    or an edge. By querying whether the point we are adding is above or below a vertex in
    the tree, or to which side of an edge it is, we can move down the tree to find the correct
    trapezoid containing the point. Lemma 1 shows the expected query time of this tree
    structure.
 *  Lemma 1: If edges are added in a random order, then after i edges have been added the
    expected depth of the tree will be O(log i).
    Proof: The tree is built by splitting the leaves representing the
    trapezoids into which a vertex or an edge is added. By randomizing the order in which
    the edges are added, at any step, the vertex to be added next resides in any trapezoid with
    equal probability, and any leaf has an equal probability of being split. Any new trapezoid
    created by the introduction of a vertex thus has an equal probability of being put below
    any of the old leaves, so trapezoids will be distributed across_product the leaves of the tree with
    equal probability. An edge splits the trapezoids between its two end points, and since
    these are distributed across_product the tree, so will the new trapezoids be. In total this ensures us
    that all the trapezoids will be randomly distributed across_product the tree. When a tree is
    generated by randomly selecting a leaf and splitting it, after j splittings, its expected depth
    will be O(log j). Each point added will create one new trapezoid, and due to merging,
    each edge will only generate one new trapezoid as well. After the addition of i edges, the
    number of trapezoids will be O(3i), and the expected depth of the tree is then O(log 3i) =
    O(log i). Once both the points of the edge have been added to the trapezoidation, the edge
    should be added too. Lemma 2 shows that the expected number of trapezoidal boundaries
    an edge will have to traverse is O(1), which is at most a constant.
 *  Lemma 2: If edges are added in a random order, then after i edges have been added the
    expected number of horizontal rays intersecting an edge not added is at most 4.
    Proof: Each point added will make two horizontal rays. Thus after j edges have been
    added to the trapezoidation, there are at most 2j points, and at most 4j horizontal rays in
    the trapezoidation, each of which may abut upon one other edge. Because of the random
    ordering of the edges, any ray is equally likely to hit the k-th edge added as the l-th edge
    added, for some k and l j, and the expected number of rays hitting any one edge is at
    most 4j / j = 4, for any j. The number of horizontal rays that will intersect an unadded
    edge after i random edges have been added, is the same as the number of rays that would
    abut upon this edge after addition if this edge were added next, i.e. after j = i + 1 edges
    had been added. By the analysis above this number is expected to be at most 4.
    This means that the number of horizontal rays an edge has to cross_product when it is being
    added is O(1), and the number of trapezoids that have to be split is one more. In total this
    leads to an algorithm which has an expected running time of O(n log n); expected
    constant time for each of n edges, and expected logarithmic time for each of n point
    locations.
 *  By using the fact that in a polygon all the edges are connected, this time bound can
    be further improved upon. The algorithm can be divided into log* n phases, where each
    phase first adds in some edges, and afterwards the remaining points are located in the
    trapezoidation by traversing edges as described above, but this time without introducing
    new regions. Finding the correct trapezoids for new points will now be cheaper, as we no
    longer have to start from the root of the tree. We get the following algorithm:
    1. Generate a random ordering of the edges
    2. For each of log* n phases do
        2.1 For each edge in the current phase do
            2.1.1 If endpoints are not added, find them through the tree and add them
            2.1.2 Add the edge between the two endpoints
        2.2 If there are more points to add then
            Locate the remaining points in the tree by tracing edges through trapezoids
    Assume that step 1 can be done linearly (See [Seidel] for a discussion of this). Step
    2.2 will expectedly take O(n) time; the location of at least one point is known in the tree,
    and tracing an edge towards the next point is O(1) for each edge, by lemma 2. Step 2.1.2
    will take expected constant time per edge, for total expected time O(n), also by lemma 2.
    The total expected time for all tree queries during a phase becomes linear in n, and thus step 2.1.1 will also
    take expected time O(n). Thus the expected overall running time of step 2, and the
    algorithm in general is O(n log* n).
 *  The algorithm is summed up in three steps given below:
 *  1. Decompose the Polygon into Trapezoids.
    Let S be a set of non-horizontal, non-intersecting line segments of the polygon .
    The randomized algorithm is used to create the trapezoidal decomposition of the X-Y plane arising due the segments
    of set S.
    This is done by taking a random ordering s1 .. sN of the segments in S and adding one segment at a time to
    incrementally construct the trapezoids. This divides the polygon into trapezoids (which can degenerate into a
    triangle if any of the horizontal segments of the trapezoid is of zero euler_distance).
    The restriction that the segments be non-horizontal is necessary to limit the number of neighbors of any trapezoid.
    However, no generality is lost due to this assumption as it can be simulated using lexicographic ordering.
    That is, if two points have the same y-coordinate then the one with larger x-coordinate is considered higher.
    The number of trapezoids is linear in the number of segments. Seidel proves that if each permutation of s1 .. sN is
    equally likely then trapezoid formation takes O(n log*n) expected time.
 *  2. Decompose the Trapezoids into Monotone Polygons.
    A monotone polygon is a polygon whose boundary consists of two y-monotone chains. These polygons are computed from
    the trapezoidal decomposition by checking whether the two vertices of the original polygon lie on the same side.
    This is a linear time operation.
 *  3. Triangulate the Monotone Polygons.
    A monotone polygon can be triangulated in linear time by using a simple greedy algorithm which repeatedly cuts off
    the convex corners of the polygon [Fournier and Montuno 1984]. Hence, all the monotone polygons can be triangulated
    in O(n) time.
 *  So the total time complexity of the algorithm is O(n*Log(n))
 *
 */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// I have included header files
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
// I have defined some useful definations */
// Some trivial ones are not commented
#define segment_size 1000		/* max# of segments. Determines how many points can be specified as input. If the dataset is having large number of points, increase this value accordingly. */
#define max_table_size   8*segment_size	/* maximum table sizes */
#define max_num_trap  4*segment_size	/* max# trapezoids */
#define collinearity_tolerance 1.0e-7		/* tolerance value: Used for making all decisions about co-linearity or left/right of segment. Decrease this value if the input points are spaced very close together */
#define cross_product(v0, v1, v2) (((v1).x - (v0).x)*((v2).y - (v0).y) - ((v1).y - (v0).y)*((v2).x - (v0).x)) //implementation of cross_product product
#define dot_product(v0, v1) ((v0).x * (v1).x + (v0).y * (v1).y) // implementation of dot_product product
#define is_abs_less_than_col_tol(s, t) (fabs(s - t) <= collinearity_tolerance)
#define cross_product_SINE(v0, v1) ((v0).x * (v1).y - (v1).x * (v0).y)
#define euler_distance(v0) (sqrt((v0).x * (v0).x + (v0).y * (v0).y))

//Now I have defined some structures which will be helpful in the program

// this structure stores a vertex
typedef struct {
    double x, y;
} point_coordinates;

// This structure stores the trapezoid attributes
typedef struct {
    int lseg, rseg;		/* two adjoining segments */
    point_coordinates hi, lo;		/* max/min y-values */
    int u0, u1;
    int d0, d1;
    int sink;			/* pointer to corresponding in Q */
    int usave, uside;
    int state;
} trapezoid_struct;

// This structure stores node attributes for every node in the query structure
typedef struct {
    int nodetype;			/* Y-node or S-node */
    int segnum;
    point_coordinates yval;
    int trnum;
    int parent;			/* doubly linked DAG */
    int left, right;		/* children */
} node_struct;

// This structure stores Segment attributes */
typedef struct {
    point_coordinates v0, v1;		/* two endpoints */
    int is_inserted;		/* inserted in trapezoidation yet ? */
    int root0, root1;		/* root nodes in Q */
    int next;			/* Next logical segment */
    int prev;			/* Previous segment */
} segment_struct;

// This structure is used to store a monotone polygon
typedef struct {
    int vnum;
    int next;			/* Circularly linked list  */
    int prev;			/* describing the monotone */
    int marked;			/* polygon */
} monotone_chain_struct;

// This structure is used to store the vertex for the 4 chains
typedef struct {
    point_coordinates pt;
    int vnext[4];			/* next vertices for the 4 chains */
    int vpos[4];			/* position of v in the 4 chains */
    int nextfree;
} vertexchain_t;

//We define some global variables
node_struct qs[max_table_size];		/* Query structure */
trapezoid_struct tr[max_num_trap];		/* Trapezoid structure */
segment_struct seg[segment_size];		/* Segment table */

static int index_q;
static int index_of_trapezoid;
static int choose_index;
static int chain_index, op_index, mon_index;
static int permute[segment_size];
static monotone_chain_struct mchain[max_num_trap]; /* Table to hold all the monotone polygons . Each monotone polygon is a circularly linked list */
static vertexchain_t vert[segment_size]; /* chain init. information. This is used to decide which monotone polygon to split if there are several other polygons touching at the same vertex  */
static int monotone_chain[segment_size];	/* contains position of any vertex in the monotone chain for the polygon */
static int visited[max_num_trap];

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class polygon_triangulation{
public:
    /* Functions */
    int monotonate_trapezoids(int);
    static int triangulate_monotone_polygons(int, int, int (*)[3]);
    static int greater_than(point_coordinates *, point_coordinates *);
    static int equal_to(point_coordinates *, point_coordinates *);
    static int greater_than_equal_to(point_coordinates *, point_coordinates *);
    static int less_than(point_coordinates *, point_coordinates *);
    static int math_logstar_n(int);
    static int math_N(int, int);
    static int inside_polygon(trapezoid_struct *t);
    static double get_angle(point_coordinates *vp0, point_coordinates *vpnext, point_coordinates *vp1);
    static int get_vertex_positions(int v0, int v1, int *ip, int *iq);
    static int make_new_monotone_polygon(int mcur, int v0, int v1);
    int traverse_polygon(int mcur, int trnum, int from, int dir);
    static int triangulate_single_polygon(int nvert, int posmax, int side, int op[][3]);
    static int create_newnode();
    static int create_trapezoid();
    int triangulate_polygon(int n, int ncontours, const int cntr[], double (*vertices)[2], int (*triangles)[3]);
    static int _max(point_coordinates *yval, point_coordinates *v0, point_coordinates *v1);
    static int _min(point_coordinates *yval, point_coordinates *v0, point_coordinates *v1);
    static int init_query_structure(int segnum);
    static int i1_of(int segnum, point_coordinates *v);
    static int inserted(int segnum, int whichpt);
    int locate_endpoint(point_coordinates *v, point_coordinates *vo, int r);
    static int merge_trapezoids(int segnum, int tfirst, int tlast, int side);
    int add_segment(int segnum);
    int construct_trapezoids(int nseg);
    int find_new_roots(int segnum);
    static int initialise(int n);
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Get log*n for given n */
int polygon_triangulation::math_logstar_n(int n){
    int i;
    double v;
    for (i = 0, v = (double) n; v >= 1; i++)
        v = log2(v);
    return (i - 1);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int polygon_triangulation::math_N(int n, int h) {
    int i;
    double v;
    for (i = 0, v = (int) n; i < h; i++)
        v = log2(v);
    return (int) ceil((double) 1.0 * n / v);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Function returns 1 if the trapezoid lies inside the polygon */
int polygon_triangulation:: inside_polygon(trapezoid_struct *t) {
    int rseg = t->rseg;
    if (t->state == 2)
        return 0;

    if ((t->lseg <= 0) || (t->rseg <= 0))
        return 0;

    if (((t->u0 <= 0) && (t->u1 <= 0)) ||
        ((t->d0 <= 0) && (t->d1 <= 0))) /* triangle */
        return (greater_than(&seg[rseg].v1, &seg[rseg].v0));

    return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double polygon_triangulation::get_angle(point_coordinates *vp0, point_coordinates *vpnext, point_coordinates *vp1) {
    point_coordinates v0, v1;
    v0.x = vpnext->x - vp0->x;
    v0.y = vpnext->y - vp0->y;
    v1.x = vp1->x - vp0->x;
    v1.y = vp1->y - vp0->y;
    if (cross_product_SINE(v0, v1) >= 0)    /* sine is positive */
        return dot_product(v0, v1) / euler_distance(v0) / euler_distance(v1);
    else
        return (-1.0 * dot_product(v0, v1) / euler_distance(v0) / euler_distance(v1) - 2);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* (v0, v1) is the new diagonal to be added to the polygon. Find which chain to use and return the positions of v0 and v1 in p and q */
int polygon_triangulation::get_vertex_positions(int v0, int v1, int *ip, int *iq) {
    vertexchain_t *vp0, *vp1;
    int i;
    double angle, temp;
    int tp, tq;
    vp0 = &vert[v0];
    vp1 = &vert[v1];
/* p is identified as follows. Scan from (v0, v1) rightwards till you hit the first segment starting from v0. That chain is the chain of our interest */
    angle = -4.0;
    for (i = 0; i < 4; i++) {
        if (vp0->vnext[i] <= 0)
            continue;
        if ((temp = get_angle(&vp0->pt, &(vert[vp0->vnext[i]].pt),
                              &vp1->pt)) > angle) {
            angle = temp;
            tp = i;
        }
    }
    *ip = tp;
/* Do similar actions for q */
    angle = -4.0;
    for (i = 0; i < 4; i++) {
        if (vp1->vnext[i] <= 0)
            continue;
        if ((temp = get_angle(&vp1->pt, &(vert[vp1->vnext[i]].pt),
                              &vp0->pt)) > angle) {
            angle = temp;
            tq = i;
        }
    }
    *iq = tq;
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* v0 and v1 are specified in anti-clockwise order with respect to the current monotone polygon mcur. Split the current polygon into two polygons using the diagonal (v0, v1)*/
int polygon_triangulation::make_new_monotone_polygon(int mcur, int v0, int v1) {
    int p, q, ip, iq;
    int mnew = ++mon_index;
    int i, j, nf0, nf1;
    vertexchain_t *vp0, *vp1;
    vp0 = &vert[v0];
    vp1 = &vert[v1];
    get_vertex_positions(v0, v1, &ip, &iq);
    p = vp0->vpos[ip];
    q = vp1->vpos[iq];

/* At this stage, we have got the positions of v0 and v1 in the desired chain. Now modify the linked lists */
    i = ++chain_index;    /* for the new list */
    j = ++chain_index;

    mchain[i].vnum = v0;
    mchain[j].vnum = v1;

    mchain[i].next = mchain[p].next;
    mchain[mchain[p].next].prev = i;
    mchain[i].prev = j;
    mchain[j].next = i;
    mchain[j].prev = mchain[q].prev;
    mchain[mchain[q].prev].next = j;

    mchain[p].next = q;
    mchain[q].prev = p;

    nf0 = vp0->nextfree;
    nf1 = vp1->nextfree;

    vp0->vnext[ip] = v1;

    vp0->vpos[nf0] = i;
    vp0->vnext[nf0] = mchain[mchain[i].next].vnum;
    vp1->vpos[nf1] = j;
    vp1->vnext[nf1] = v0;

    vp0->nextfree++;
    vp1->nextfree++;

    monotone_chain[mcur] = p;
    monotone_chain[mnew] = i;
    return mnew;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/* Main routine to get monotone polygons from the trapezoidation of the polygon.*/
int polygon_triangulation::monotonate_trapezoids(int n){
    int i,tr_start;
    memset((void *)vert, 0, sizeof(vert));
    memset((void *)visited, 0, sizeof(visited));
    memset((void *)mchain, 0, sizeof(mchain));
    memset((void *)monotone_chain, 0, sizeof(monotone_chain));
    /* First locate a trapezoid which lies inside the polygon and which is triangular */
    for (i = 0; i < max_num_trap; i++)
        if (inside_polygon(&tr[i]))
            break;
    tr_start = i;

    /* Initialise the monotone_chain data-structure and start spanning all the */
    /* trapezoids within the polygon */

    for (i = 1; i <= n; i++) {
        mchain[i].prev = seg[i].prev;
        mchain[i].next = seg[i].next;
        mchain[i].vnum = i;
        vert[i].pt = seg[i].v0;
        vert[i].vnext[0] = seg[i].next; /* next vertex */
        vert[i].vpos[0] = i;    /* locn. of next vertex */
        vert[i].nextfree = 1;
    }

    chain_index = n;
    mon_index = 0;
    monotone_chain[0] = 1;			/* position of any vertex in the first */
    /* chain  */


    /* traverse the polygon */
    if (tr[tr_start].u0 > 0)
        traverse_polygon(0, tr_start, tr[tr_start].u0, 1);
    else if (tr[tr_start].d0 > 0)
        traverse_polygon(0, tr_start, tr[tr_start].d0, 2);

    /* return the number of polygons created */
    return ++mon_index;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* recursively visit all the trapezoids */
int polygon_triangulation::traverse_polygon(int mcur, int trnum, int from, int dir) {
    trapezoid_struct *t = &tr[trnum];
    int howsplit, mnew;
    int v0, v1, v0next, v1next;
    int retval, tmp;
    int do_switch = 0;

    if ((trnum <= 0) || visited[trnum])
        return 0;

    visited[trnum] = 1;

/* We have much more information available here. */
/* rseg: goes upwards   */
/* lseg: goes downwards */

/* Initially assume that dir = 2 (from the left) */
/* Switch v0 and v1 if necessary afterwards */


/* special cases for triangles with cusps at the opposite ends. */
/* take care of this first */
    if ((t->u0 <= 0) && (t->u1 <= 0)) {
        if ((t->d0 > 0) && (t->d1 > 0)) /* downward opening triangle */
        {
            v0 = tr[t->d1].lseg;
            v1 = t->lseg;
            if (from == t->d1) {
                do_switch = 1;
                mnew = make_new_monotone_polygon(mcur, v1, v0);
                traverse_polygon(mcur, t->d1, trnum, 1);
                traverse_polygon(mnew, t->d0, trnum, 1);
            } else {
                mnew = make_new_monotone_polygon(mcur, v0, v1);
                traverse_polygon(mcur, t->d0, trnum, 1);
                traverse_polygon(mnew, t->d1, trnum, 1);
            }
        } else {
            retval = -1;    /* Just traverse all neighbours */
            traverse_polygon(mcur, t->u0, trnum, 2);
            traverse_polygon(mcur, t->u1, trnum, 2);
            traverse_polygon(mcur, t->d0, trnum, 1);
            traverse_polygon(mcur, t->d1, trnum, 1);
        }
    } else if ((t->d0 <= 0) && (t->d1 <= 0)) {
        if ((t->u0 > 0) && (t->u1 > 0)) /* upward opening triangle */
        {
            v0 = t->rseg;
            v1 = tr[t->u0].rseg;
            if (from == t->u1) {
                do_switch = 1;
                mnew = make_new_monotone_polygon(mcur, v1, v0);
                traverse_polygon(mcur, t->u1, trnum, 2);
                traverse_polygon(mnew, t->u0, trnum, 2);
            } else {
                mnew = make_new_monotone_polygon(mcur, v0, v1);
                traverse_polygon(mcur, t->u0, trnum, 2);
                traverse_polygon(mnew, t->u1, trnum, 2);
            }
        } else {
            retval = -1;    /* Just traverse all neighbours */
            traverse_polygon(mcur, t->u0, trnum, 2);
            traverse_polygon(mcur, t->u1, trnum, 2);
            traverse_polygon(mcur, t->d0, trnum, 1);
            traverse_polygon(mcur, t->d1, trnum, 1);
        }
    } else if ((t->u0 > 0) && (t->u1 > 0)) {
        if ((t->d0 > 0) && (t->d1 > 0)) /* downward + upward cusps */
        {
            v0 = tr[t->d1].lseg;
            v1 = tr[t->u0].rseg;
            retval = 3;
            if (((dir == 2) && (t->d1 == from)) ||
                ((dir == 1) && (t->u1 == from))) {
                do_switch = 1;
                mnew = make_new_monotone_polygon(mcur, v1, v0);
                traverse_polygon(mcur, t->u1, trnum, 2);
                traverse_polygon(mcur, t->d1, trnum, 1);
                traverse_polygon(mnew, t->u0, trnum, 2);
                traverse_polygon(mnew, t->d0, trnum, 1);
            } else {
                mnew = make_new_monotone_polygon(mcur, v0, v1);
                traverse_polygon(mcur, t->u0, trnum, 2);
                traverse_polygon(mcur, t->d0, trnum, 1);
                traverse_polygon(mnew, t->u1, trnum, 2);
                traverse_polygon(mnew, t->d1, trnum, 1);
            }
        } else            /* only downward cusp */
        {
            if (equal_to(&t->lo, &seg[t->lseg].v1)) {
                v0 = tr[t->u0].rseg;
                v1 = seg[t->lseg].next;

                retval = 4;
                if ((dir == 1) && (t->u0 == from)) {
                    do_switch = 1;
                    mnew = make_new_monotone_polygon(mcur, v1, v0);
                    traverse_polygon(mcur, t->u0, trnum, 2);
                    traverse_polygon(mnew, t->d0, trnum, 1);
                    traverse_polygon(mnew, t->u1, trnum, 2);
                    traverse_polygon(mnew, t->d1, trnum, 1);
                } else {
                    mnew = make_new_monotone_polygon(mcur, v0, v1);
                    traverse_polygon(mcur, t->u1, trnum, 2);
                    traverse_polygon(mcur, t->d0, trnum, 1);
                    traverse_polygon(mcur, t->d1, trnum, 1);
                    traverse_polygon(mnew, t->u0, trnum, 2);
                }
            } else {
                v0 = t->rseg;
                v1 = tr[t->u0].rseg;
                retval = 5;
                if ((dir == 1) && (t->u1 == from)) {
                    do_switch = 1;
                    mnew = make_new_monotone_polygon(mcur, v1, v0);
                    traverse_polygon(mcur, t->u1, trnum, 2);
                    traverse_polygon(mnew, t->d1, trnum, 1);
                    traverse_polygon(mnew, t->d0, trnum, 1);
                    traverse_polygon(mnew, t->u0, trnum, 2);
                } else {
                    mnew = make_new_monotone_polygon(mcur, v0, v1);
                    traverse_polygon(mcur, t->u0, trnum, 2);
                    traverse_polygon(mcur, t->d0, trnum, 1);
                    traverse_polygon(mcur, t->d1, trnum, 1);
                    traverse_polygon(mnew, t->u1, trnum, 2);
                }
            }
        }
    } else if ((t->u0 > 0) || (t->u1 > 0)) /* no downward cusp */
    {
        if ((t->d0 > 0) && (t->d1 > 0)) /* only upward cusp */
        {
            if (equal_to(&t->hi, &seg[t->lseg].v0)) {
                v0 = tr[t->d1].lseg;
                v1 = t->lseg;
                retval = 6;
                if (!((dir == 2) && (t->d0 == from))) {
                    do_switch = 1;
                    mnew = make_new_monotone_polygon(mcur, v1, v0);
                    traverse_polygon(mcur, t->u1, trnum, 2);
                    traverse_polygon(mcur, t->d1, trnum, 1);
                    traverse_polygon(mcur, t->u0, trnum, 2);
                    traverse_polygon(mnew, t->d0, trnum, 1);
                } else {
                    mnew = make_new_monotone_polygon(mcur, v0, v1);
                    traverse_polygon(mcur, t->d0, trnum, 1);
                    traverse_polygon(mnew, t->u0, trnum, 2);
                    traverse_polygon(mnew, t->u1, trnum, 2);
                    traverse_polygon(mnew, t->d1, trnum, 1);
                }
            } else {
                v0 = tr[t->d1].lseg;
                v1 = seg[t->rseg].next;

                retval = 7;
                if ((dir == 2) && (t->d1 == from)) {
                    do_switch = 1;
                    mnew = make_new_monotone_polygon(mcur, v1, v0);
                    traverse_polygon(mcur, t->d1, trnum, 1);
                    traverse_polygon(mnew, t->u1, trnum, 2);
                    traverse_polygon(mnew, t->u0, trnum, 2);
                    traverse_polygon(mnew, t->d0, trnum, 1);
                } else {
                    mnew = make_new_monotone_polygon(mcur, v0, v1);
                    traverse_polygon(mcur, t->u0, trnum, 2);
                    traverse_polygon(mcur, t->d0, trnum, 1);
                    traverse_polygon(mcur, t->u1, trnum, 2);
                    traverse_polygon(mnew, t->d1, trnum, 1);
                }
            }
        } else            /* no cusp */
        {
            if (equal_to(&t->hi, &seg[t->lseg].v0) &&
                equal_to(&t->lo, &seg[t->rseg].v0)) {
                v0 = t->rseg;
                v1 = t->lseg;
                retval = 2;
                if (dir == 1) {
                    do_switch = 1;
                    mnew = make_new_monotone_polygon(mcur, v1, v0);
                    traverse_polygon(mcur, t->u0, trnum, 2);
                    traverse_polygon(mcur, t->u1, trnum, 2);
                    traverse_polygon(mnew, t->d1, trnum, 1);
                    traverse_polygon(mnew, t->d0, trnum, 1);
                } else {
                    mnew = make_new_monotone_polygon(mcur, v0, v1);
                    traverse_polygon(mcur, t->d1, trnum, 1);
                    traverse_polygon(mcur, t->d0, trnum, 1);
                    traverse_polygon(mnew, t->u0, trnum, 2);
                    traverse_polygon(mnew, t->u1, trnum, 2);
                }
            } else if (equal_to(&t->hi, &seg[t->rseg].v1) &&
                       equal_to(&t->lo, &seg[t->lseg].v1)) {
                v0 = seg[t->rseg].next;
                v1 = seg[t->lseg].next;

                retval = 1;
                if (dir == 1) {
                    do_switch = 1;
                    mnew = make_new_monotone_polygon(mcur, v1, v0);
                    traverse_polygon(mcur, t->u0, trnum, 2);
                    traverse_polygon(mcur, t->u1, trnum, 2);
                    traverse_polygon(mnew, t->d1, trnum, 1);
                    traverse_polygon(mnew, t->d0, trnum, 1);
                } else {
                    mnew = make_new_monotone_polygon(mcur, v0, v1);
                    traverse_polygon(mcur, t->d1, trnum, 1);
                    traverse_polygon(mcur, t->d0, trnum, 1);
                    traverse_polygon(mnew, t->u0, trnum, 2);
                    traverse_polygon(mnew, t->u1, trnum, 2);
                }
            } else            /* no split possible */
            {
                retval = -1;
                traverse_polygon(mcur, t->u0, trnum, 2);
                traverse_polygon(mcur, t->d0, trnum, 1);
                traverse_polygon(mcur, t->u1, trnum, 2);
                traverse_polygon(mcur, t->d1, trnum, 1);
            }
        }
    }

    return retval;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* For each monotone polygon, find the ymax and ymin (to determine the */
/* two y-monotone chains) and pass on this monotone polygon for greedy */
/* triangulation. */
/* Take care not to triangulate duplicate monotone polygons */

int polygon_triangulation:: triangulate_monotone_polygons(int nvert, int nmonpoly, int op[][3]) {
    int i;
    point_coordinates ymax, ymin;
    int p, vfirst, posmax, posmin, v;
    int vcount, processed;

    op_index = 0;
    for (i = 0; i < nmonpoly; i++) {
        vcount = 1;
        processed = 0;
        vfirst = mchain[monotone_chain[i]].vnum;
        ymax = ymin = vert[vfirst].pt;
        posmax = posmin = monotone_chain[i];
        mchain[monotone_chain[i]].marked = 1;
        p = mchain[monotone_chain[i]].next;
        while ((v = mchain[p].vnum) != vfirst) {
            if (mchain[p].marked) {
                processed = 1;
                break;        /* break from while */
            } else
                mchain[p].marked = 1;

            if (greater_than(&vert[v].pt, &ymax)) {
                ymax = vert[v].pt;
                posmax = p;
            }
            if (less_than(&vert[v].pt, &ymin)) {
                ymin = vert[v].pt;
                posmin = p;
            }
            p = mchain[p].next;
            vcount++;
        }

        if (processed)        /* Go to next polygon */
            continue;

        if (vcount == 3)        /* already a triangle */
        {
            op[op_index][0] = mchain[p].vnum;
            op[op_index][1] = mchain[mchain[p].next].vnum;
            op[op_index][2] = mchain[mchain[p].prev].vnum;
            op_index++;
        } else            /* triangulate the polygon */
        {
            v = mchain[mchain[posmax].next].vnum;
            if (equal_to(&vert[v].pt, &ymin)) {            /* LHS is a single line */
                triangulate_single_polygon(nvert, posmax, 1, op);
            } else
                triangulate_single_polygon(nvert, posmax, 2, op);
        }
    }
    return op_index;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* A greedy corner-cutting algorithm to triangulate a y-monotone polygon in O(n) time. I have taken the idea from
 * Joseph O-Rourke, Computational Geometry.*/
int polygon_triangulation::triangulate_single_polygon(int nvert, int posmax, int side, int op[][3]) {
    int v;
    int rc[segment_size], ri = 0;    /* reflex chain */
    int endv, tmp, vpos;

    if (side == 2)        /* RHS segment is a single segment */
    {
        rc[0] = mchain[posmax].vnum;
        tmp = mchain[posmax].next;
        rc[1] = mchain[tmp].vnum;
        ri = 1;

        vpos = mchain[tmp].next;
        v = mchain[vpos].vnum;

        if ((endv = mchain[mchain[posmax].prev].vnum) == 0)
            endv = nvert;
    } else                /* LHS is a single segment */
    {
        tmp = mchain[posmax].next;
        rc[0] = mchain[tmp].vnum;
        tmp = mchain[tmp].next;
        rc[1] = mchain[tmp].vnum;
        ri = 1;

        vpos = mchain[tmp].next;
        v = mchain[vpos].vnum;

        endv = mchain[posmax].vnum;
    }

    while ((v != endv) || (ri > 1)) {
        if (ri > 0)        /* reflex chain is non-empty */
        {
            if (cross_product(vert[v].pt, vert[rc[ri - 1]].pt,vert[rc[ri]].pt) > 0) { /* convex corner: cut if off */
                op[op_index][0] = rc[ri - 1];
                op[op_index][1] = rc[ri];
                op[op_index][2] = v;
                op_index++;
                ri--;
            } else        /* non-convex */
            {        /* add v to the chain */
                ri++;
                rc[ri] = v;
                vpos = mchain[vpos].next;
                v = mchain[vpos].vnum;
            }
        } else            /* reflex-chain empty: add v to the */
        {            /* reflex chain and advance it  */
            rc[++ri] = v;
            vpos = mchain[vpos].next;
            v = mchain[vpos].vnum;
        }
    } /* end-while */

/* reached the bottom vertex. Add in the triangle formed */
    op[op_index][0] = rc[ri - 1];
    op[op_index][1] = rc[ri];
    op[op_index][2] = v;
    op_index++;
    ri--;

    return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Return a new node to be added into the query tree */
int polygon_triangulation::create_newnode(){
    if (index_q < max_table_size)
        return index_q++;
    else
    {
        fprintf(stderr, "create_newnode: Query-table overflow\n");
        return -1;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Return a free trapezoid */
int polygon_triangulation::create_trapezoid(){
    if (index_of_trapezoid < max_num_trap)
    {
        tr[index_of_trapezoid].lseg = -1;
        tr[index_of_trapezoid].rseg = -1;
        tr[index_of_trapezoid].state = 1;
        return index_of_trapezoid++;
    }
    else
    {
        fprintf(stderr, "create_trapezoid: Trapezoid-table overflow\n");
        return -1;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Return the maximum of the two points into the yval structure */
int polygon_triangulation::_max(point_coordinates *yval, point_coordinates *v0, point_coordinates *v1) {
    if (v0->y > v1->y + collinearity_tolerance)
        *yval = *v0;
    else if (is_abs_less_than_col_tol(v0->y, v1->y)) {
        if (v0->x > v1->x + collinearity_tolerance)
            *yval = *v0;
        else
            *yval = *v1;
    } else
        *yval = *v1;

    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Return the minimum of the two points into the yval structure */
int polygon_triangulation::_min(point_coordinates *yval, point_coordinates *v0, point_coordinates *v1) {
    if (v0->y < v1->y - collinearity_tolerance)
        *yval = *v0;
    else if (is_abs_less_than_col_tol(v0->y, v1->y)) {
        if (v0->x < v1->x)
            *yval = *v0;
        else
            *yval = *v1;
    } else
        *yval = *v1;

    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int polygon_triangulation :: greater_than(point_coordinates *v0, point_coordinates *v1) {
    if (v0->y > v1->y + collinearity_tolerance)
        return 1;
    else if (v0->y < v1->y - collinearity_tolerance)
        return 0;
    else
        return (v0->x > v1->x);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int polygon_triangulation ::equal_to(point_coordinates *v0, point_coordinates *v1){
    return (is_abs_less_than_col_tol(v0->y, v1->y) && is_abs_less_than_col_tol(v0->x, v1->x));
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int polygon_triangulation ::greater_than_equal_to(point_coordinates *v0, point_coordinates *v1) {
    if (v0->y > v1->y + collinearity_tolerance)
        return 1;
    else if (v0->y < v1->y - collinearity_tolerance)
        return 0;
    else
        return (v0->x >= v1->x);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int polygon_triangulation ::less_than(point_coordinates *v0, point_coordinates *v1) {
    if (v0->y < v1->y - collinearity_tolerance)
        return 1;
    else if (v0->y > v1->y + collinearity_tolerance)
        return 0;
    else
        return (v0->x < v1->x);
}


/* Initilialise the query structure (Q) and the trapezoid table (T)
 * when the first segment is added to start the trapezoidation. The
 * query-tree starts out with 4 trapezoids, one S-node and 2 Y-nodes
 *
 *                4
 *   -----------------------------------
 *  		  \
 *  	1	   \        2
 *  		    \
 *   -----------------------------------
 *                3
 */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int polygon_triangulation::init_query_structure(int segnum) {
    int i1, i2, i3, i4, i5, i6, i7, root;
    int t1, t2, t3, t4;
    segment_struct *s = &seg[segnum];

    index_q = index_of_trapezoid = 1;
    memset((void *) tr, 0, sizeof(tr));
    memset((void *) qs, 0, sizeof(qs));

    i1 = create_newnode();
    qs[i1].nodetype = 2;
    _max(&qs[i1].yval, &s->v0, &s->v1); /* root */
    root = i1;

    qs[i1].right = i2 = create_newnode();
    qs[i2].nodetype = 3;
    qs[i2].parent = i1;

    qs[i1].left = i3 = create_newnode();
    qs[i3].nodetype = 2;
    _min(&qs[i3].yval, &s->v0, &s->v1); /* root */
    qs[i3].parent = i1;

    qs[i3].left = i4 = create_newnode();
    qs[i4].nodetype = 3;
    qs[i4].parent = i3;

    qs[i3].right = i5 = create_newnode();
    qs[i5].nodetype = 1;
    qs[i5].segnum = segnum;
    qs[i5].parent = i3;

    qs[i5].left = i6 = create_newnode();
    qs[i6].nodetype = 3;
    qs[i6].parent = i5;

    qs[i5].right = i7 = create_newnode();
    qs[i7].nodetype = 3;
    qs[i7].parent = i5;

    t1 = create_trapezoid();        /* middle left */
    t2 = create_trapezoid();        /* middle right */
    t3 = create_trapezoid();        /* bottom-most */
    t4 = create_trapezoid();        /* topmost */

    tr[t1].hi = tr[t2].hi = tr[t4].lo = qs[i1].yval;
    tr[t1].lo = tr[t2].lo = tr[t3].hi = qs[i3].yval;
    tr[t4].hi.y = (double) (INFINITY);
    tr[t4].hi.x = (double) (INFINITY);
    tr[t3].lo.y = (double) -1 * (double) (INFINITY);
    tr[t3].lo.x = (double) -1 * (double) (INFINITY);
    tr[t1].rseg = tr[t2].lseg = segnum;
    tr[t1].u0 = tr[t2].u0 = t4;
    tr[t1].d0 = tr[t2].d0 = t3;
    tr[t4].d0 = tr[t3].u0 = t1;
    tr[t4].d1 = tr[t3].u1 = t2;

    tr[t1].sink = i6;
    tr[t2].sink = i7;
    tr[t3].sink = i4;
    tr[t4].sink = i2;

    tr[t1].state = tr[t2].state = 1;
    tr[t3].state = tr[t4].state = 1;

    qs[i2].trnum = t4;
    qs[i4].trnum = t3;
    qs[i6].trnum = t1;
    qs[i7].trnum = t2;

    s->is_inserted = 1;
    return root;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Return 1 if the vertex v is to the left of line segment no. segnum. Takes care of the degenerate cases when both the vertices have the same y--cood, etc.*/
int polygon_triangulation::i1_of(int segnum, point_coordinates *v) {
    segment_struct *s = &seg[segnum];
    double area;

    if (greater_than(&s->v1, &s->v0)) /* seg. going upwards */
    {
        if (is_abs_less_than_col_tol(s->v1.y, v->y)) {
            if (v->x < s->v1.x)
                area = 1.0;
            else
                area = -1.0;
        } else if (is_abs_less_than_col_tol(s->v0.y, v->y)) {
            if (v->x < s->v0.x)
                area = 1.0;
            else
                area = -1.0;
        } else
            area = cross_product(s->v0, s->v1, (*v));
    } else                /* v0 > v1 */
    {
        if (is_abs_less_than_col_tol(s->v1.y, v->y)) {
            if (v->x < s->v1.x)
                area = 1.0;
            else
                area = -1.0;
        } else if (is_abs_less_than_col_tol(s->v0.y, v->y)) {
            if (v->x < s->v0.x)
                area = 1.0;
            else
                area = -1.0;
        } else
            area = cross_product(s->v1, s->v0, (*v));
    }

    if (area > 0.0)
        return 1;
    else
        return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Returns 1 if the corresponding endpoint of the given segment is already inserted into the segment tree. Use the simple test of whether the segment which shares this endpoint is already inserted */

int polygon_triangulation::inserted(int segnum, int whichpt) {
    if (whichpt == 1)
        return seg[seg[segnum].prev].is_inserted;
    else
        return seg[seg[segnum].next].is_inserted;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* This is query routine which determines which trapezoid does the point v lie in. The return value is the trapezoid number. */

int polygon_triangulation::locate_endpoint(point_coordinates *v, point_coordinates *vo, int r) {
    node_struct *rptr = &qs[r];
    switch (rptr->nodetype) {
        case 3:
            return rptr->trnum;

        case 2:
            if (greater_than(v, &rptr->yval)) /* above */
                return locate_endpoint(v, vo, rptr->right);
            else if (equal_to(v, &rptr->yval)) /* the point is already */
            {                      /* inserted. */
                if (greater_than(vo, &rptr->yval)) /* above */
                    return locate_endpoint(v, vo, rptr->right);
                else
                    return locate_endpoint(v, vo, rptr->left); /* below */
            } else
                return locate_endpoint(v, vo, rptr->left); /* below */

        case 1:
            if (equal_to(v, &seg[rptr->segnum].v0) ||
                equal_to(v, &seg[rptr->segnum].v1)) {
                if (is_abs_less_than_col_tol(v->y, vo->y)) /* horizontal segment */
                {
                    if (vo->x < v->x)
                        return locate_endpoint(v, vo, rptr->left); /* left */
                    else
                        return locate_endpoint(v, vo, rptr->right); /* right */
                } else if (i1_of(rptr->segnum, vo))
                    return locate_endpoint(v, vo, rptr->left); /* left */
                else
                    return locate_endpoint(v, vo, rptr->right); /* right */
            } else if (i1_of(rptr->segnum, v))
                return locate_endpoint(v, vo, rptr->left); /* left */
            else
                return locate_endpoint(v, vo, rptr->right); /* right */

        default:
            //fprintf(stderr, "Haggu !!!!!\n");
            break;
    }
    return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Thread in the segment into the existing trapezoidation. The limiting trapezoids are given by tfirst and tlast (which are the trapezoids containing the two endpoints of the segment. Merges all
 * possible trapezoids which flank this segment and have been recently divided because of its insertion */

int polygon_triangulation::merge_trapezoids(int segnum, int tfirst, int tlast, int side) {
    int t, tnext, cond;
    int ptnext;

/* First merge polys on the LHS */
    t = tfirst;
    while ((t > 0) && greater_than_equal_to(&tr[t].lo, &tr[tlast].lo)) {
        if (side == 1)
            cond = ((((tnext = tr[t].d0) > 0) && (tr[tnext].rseg == segnum)) ||
                    (((tnext = tr[t].d1) > 0) && (tr[tnext].rseg == segnum)));
        else
            cond = ((((tnext = tr[t].d0) > 0) && (tr[tnext].lseg == segnum)) ||
                    (((tnext = tr[t].d1) > 0) && (tr[tnext].lseg == segnum)));

        if (cond) {
            if ((tr[t].lseg == tr[tnext].lseg) &&
                (tr[t].rseg == tr[tnext].rseg)) /* good neighbours */
            {                          /* merge them */
/* Use the upper node as the new node i.e. t */

                ptnext = qs[tr[tnext].sink].parent;

                if (qs[ptnext].left == tr[tnext].sink)
                    qs[ptnext].left = tr[t].sink;
                else
                    qs[ptnext].right = tr[t].sink;    /* redirect parent */


/* Change the upper neighbours of the lower trapezoids */

                if ((tr[t].d0 = tr[tnext].d0) > 0) {
                    if (tr[tr[t].d0].u0 != tnext) {
                        if (tr[tr[t].d0].u1 == tnext)
                            tr[tr[t].d0].u1 = t;
                    } else tr[tr[t].d0].u0 = t;
                }

                if ((tr[t].d1 = tr[tnext].d1) > 0) {
                    if (tr[tr[t].d1].u0 != tnext) {
                        if (tr[tr[t].d1].u1 == tnext)
                            tr[tr[t].d1].u1 = t;
                    } else tr[tr[t].d1].u0 = t;
                }

                tr[t].lo = tr[tnext].lo;
                tr[tnext].state = 2; /* invalidate the lower trapezium */
            } else            /* not good neighbours */
                t = tnext;
        } else            /* do not satisfy the outer if */
            t = tnext;
    } /* end-while */
    return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/* Add in the new segment into the trapezoidation and update Q and T
 * structures. First locate the two endpoints of the segment in the
 * Q-structure. Then start from the topmost trapezoid and go down to
 * the  lower trapezoid dividing all the trapezoids in between .
 */

int polygon_triangulation::add_segment(int segnum){
    segment_struct s;
    segment_struct *so = &seg[segnum];
    int tu, tl, sk, tfirst, tlast, tnext;
    int tfirstr, tlastr, tfirstl, tlastl;
    int i1, i2, t, t1, t2, tn;
    point_coordinates tpt;
    int tritop = 0, tribot = 0, is_swapped = 0;
    int tmptriseg;

    s = seg[segnum];
    if (greater_than(&s.v1, &s.v0)) /* Get higher vertex in v0 */
    {
        int tmp;
        tpt = s.v0;
        s.v0 = s.v1;
        s.v1 = tpt;
        tmp = s.root0;
        s.root0 = s.root1;
        s.root1 = tmp;
        is_swapped = 1;
    }

    if ((is_swapped) ? !inserted(segnum, 2) :
        !inserted(segnum, 1))     /* insert v0 in the tree */
    {
        int tmp_d;

        tu = locate_endpoint(&s.v0, &s.v1, s.root0);
        tl = create_trapezoid();		/* tl is the new lower trapezoid */
        tr[tl].state = 1;
        tr[tl] = tr[tu];
        tr[tu].lo.y = tr[tl].hi.y = s.v0.y;
        tr[tu].lo.x = tr[tl].hi.x = s.v0.x;
        tr[tu].d0 = tl;
        tr[tu].d1 = 0;
        tr[tl].u0 = tu;
        tr[tl].u1 = 0;

        if (((tmp_d = tr[tl].d0) > 0) && (tr[tmp_d].u0 == tu))
            tr[tmp_d].u0 = tl;
        if (((tmp_d = tr[tl].d0) > 0) && (tr[tmp_d].u1 == tu))
            tr[tmp_d].u1 = tl;

        if (((tmp_d = tr[tl].d1) > 0) && (tr[tmp_d].u0 == tu))
            tr[tmp_d].u0 = tl;
        if (((tmp_d = tr[tl].d1) > 0) && (tr[tmp_d].u1 == tu))
            tr[tmp_d].u1 = tl;

        /* Now update the query structure and obtain the sinks for the */
        /* two trapezoids */

        i1 = create_newnode();		/* Upper trapezoid sink */
        i2 = create_newnode();		/* Lower trapezoid sink */
        sk = tr[tu].sink;

        qs[sk].nodetype = 2;
        qs[sk].yval = s.v0;
        qs[sk].segnum = segnum;	/* not really reqd ... maybe later */
        qs[sk].left = i2;
        qs[sk].right = i1;

        qs[i1].nodetype = 3;
        qs[i1].trnum = tu;
        qs[i1].parent = sk;

        qs[i2].nodetype = 3;
        qs[i2].trnum = tl;
        qs[i2].parent = sk;

        tr[tu].sink = i1;
        tr[tl].sink = i2;
        tfirst = tl;
    }
    else				/* v0 already present */
    {       /* Get the topmost intersecting trapezoid */
        tfirst = locate_endpoint(&s.v0, &s.v1, s.root0);
        tritop = 1;
    }


    if ((is_swapped) ? !inserted(segnum, 1) :
        !inserted(segnum, 2))     /* insert v1 in the tree */
    {
        int tmp_d;

        tu = locate_endpoint(&s.v1, &s.v0, s.root1);

        tl = create_trapezoid();		/* tl is the new lower trapezoid */
        tr[tl].state = 1;
        tr[tl] = tr[tu];
        tr[tu].lo.y = tr[tl].hi.y = s.v1.y;
        tr[tu].lo.x = tr[tl].hi.x = s.v1.x;
        tr[tu].d0 = tl;
        tr[tu].d1 = 0;
        tr[tl].u0 = tu;
        tr[tl].u1 = 0;

        if (((tmp_d = tr[tl].d0) > 0) && (tr[tmp_d].u0 == tu))
            tr[tmp_d].u0 = tl;
        if (((tmp_d = tr[tl].d0) > 0) && (tr[tmp_d].u1 == tu))
            tr[tmp_d].u1 = tl;

        if (((tmp_d = tr[tl].d1) > 0) && (tr[tmp_d].u0 == tu))
            tr[tmp_d].u0 = tl;
        if (((tmp_d = tr[tl].d1) > 0) && (tr[tmp_d].u1 == tu))
            tr[tmp_d].u1 = tl;

        /* Now update the query structure and obtain the sinks for the */
        /* two trapezoids */

        i1 = create_newnode();		/* Upper trapezoid sink */
        i2 = create_newnode();		/* Lower trapezoid sink */
        sk = tr[tu].sink;

        qs[sk].nodetype = 2;
        qs[sk].yval = s.v1;
        qs[sk].segnum = segnum;	/* not really reqd ... maybe later */
        qs[sk].left = i2;
        qs[sk].right = i1;

        qs[i1].nodetype = 3;
        qs[i1].trnum = tu;
        qs[i1].parent = sk;

        qs[i2].nodetype = 3;
        qs[i2].trnum = tl;
        qs[i2].parent = sk;

        tr[tu].sink = i1;
        tr[tl].sink = i2;
        tlast = tu;
    }
    else				/* v1 already present */
    {       /* Get the lowermost intersecting trapezoid */
        tlast = locate_endpoint(&s.v1, &s.v0, s.root1);
        tribot = 1;
    }

    /* Thread the segment into the query tree creating a new X-node */
    /* First, split all the trapezoids which are intersected by s into */
    /* two */

    t = tfirst;			/* topmost trapezoid */

    while ((t > 0) &&
           greater_than_equal_to(&tr[t].lo, &tr[tlast].lo))
        /* traverse from top to bot */
    {
        int t_sav, tn_sav;
        sk = tr[t].sink;
        i1 = create_newnode();		/* left trapezoid sink */
        i2 = create_newnode();		/* right trapezoid sink */

        qs[sk].nodetype = 1;
        qs[sk].segnum = segnum;
        qs[sk].left = i1;
        qs[sk].right = i2;

        qs[i1].nodetype = 3;	/* left trapezoid (use existing one) */
        qs[i1].trnum = t;
        qs[i1].parent = sk;

        qs[i2].nodetype = 3;	/* right trapezoid (allocate new) */
        qs[i2].trnum = tn = create_trapezoid();
        tr[tn].state = 1;
        qs[i2].parent = sk;

        if (t == tfirst)
            tfirstr = tn;
        if (equal_to(&tr[t].lo, &tr[tlast].lo))
            tlastr = tn;

        tr[tn] = tr[t];
        tr[t].sink = i1;
        tr[tn].sink = i2;
        t_sav = t;
        tn_sav = tn;

        /* error */

        if ((tr[t].d0 <= 0) && (tr[t].d1 <= 0)) /* case cannot arise */
        {
            fprintf(stderr, "add_segment: error\n");
            break;
        }

            /* only one trapezoid below. partition t into two and make the */
            /* two resulting trapezoids t and tn as the upper neighbours of */
            /* the sole lower trapezoid */

        else if ((tr[t].d0 > 0) && (tr[t].d1 <= 0))
        {			/* Only one trapezoid below */
            if ((tr[t].u0 > 0) && (tr[t].u1 > 0))
            {			/* continuation of a chain from abv. */
                if (tr[t].usave > 0) /* three upper neighbours */
                {
                    if (tr[t].uside == 1)
                    {
                        tr[tn].u0 = tr[t].u1;
                        tr[t].u1 = -1;
                        tr[tn].u1 = tr[t].usave;

                        tr[tr[t].u0].d0 = t;
                        tr[tr[tn].u0].d0 = tn;
                        tr[tr[tn].u1].d0 = tn;
                    }
                    else		/* intersects in the right */
                    {
                        tr[tn].u1 = -1;
                        tr[tn].u0 = tr[t].u1;
                        tr[t].u1 = tr[t].u0;
                        tr[t].u0 = tr[t].usave;

                        tr[tr[t].u0].d0 = t;
                        tr[tr[t].u1].d0 = t;
                        tr[tr[tn].u0].d0 = tn;
                    }

                    tr[t].usave = tr[tn].usave = 0;
                }
                else		/* No usave.... simple case */
                {
                    tr[tn].u0 = tr[t].u1;
                    tr[t].u1 = tr[tn].u1 = -1;
                    tr[tr[tn].u0].d0 = tn;
                }
            }
            else
            {			/* fresh seg. or upward cusp */
                int tmp_u = tr[t].u0;
                int td0, td1;
                if (((td0 = tr[tmp_u].d0) > 0) &&
                    ((td1 = tr[tmp_u].d1) > 0))
                {		/* upward cusp */
                    if ((tr[td0].rseg > 0) &&
                        !i1_of(tr[td0].rseg, &s.v1))
                    {
                        tr[t].u0 = tr[t].u1 = tr[tn].u1 = -1;
                        tr[tr[tn].u0].d1 = tn;
                    }
                    else		/* cusp going leftwards */
                    {
                        tr[tn].u0 = tr[tn].u1 = tr[t].u1 = -1;
                        tr[tr[t].u0].d0 = t;
                    }
                }
                else		/* fresh segment */
                {
                    tr[tr[t].u0].d0 = t;
                    tr[tr[t].u0].d1 = tn;
                }
            }

            if (is_abs_less_than_col_tol(tr[t].lo.y, tr[tlast].lo.y) &&
                is_abs_less_than_col_tol(tr[t].lo.x, tr[tlast].lo.x) && tribot)
            {		/* bottom forms a triangle */

                if (is_swapped)
                    tmptriseg = seg[segnum].prev;
                else
                    tmptriseg = seg[segnum].next;

                if ((tmptriseg > 0) && i1_of(tmptriseg, &s.v0))
                {
                    /* L-R downward cusp */
                    tr[tr[t].d0].u0 = t;
                    tr[tn].d0 = tr[tn].d1 = -1;
                }
                else
                {
                    /* R-L downward cusp */
                    tr[tr[tn].d0].u1 = tn;
                    tr[t].d0 = tr[t].d1 = -1;
                }
            }
            else
            {
                if ((tr[tr[t].d0].u0 > 0) && (tr[tr[t].d0].u1 > 0))
                {
                    if (tr[tr[t].d0].u0 == t) /* passes thru LHS */
                    {
                        tr[tr[t].d0].usave = tr[tr[t].d0].u1;
                        tr[tr[t].d0].uside = 1;
                    }
                    else
                    {
                        tr[tr[t].d0].usave = tr[tr[t].d0].u0;
                        tr[tr[t].d0].uside = 2;
                    }
                }
                tr[tr[t].d0].u0 = t;
                tr[tr[t].d0].u1 = tn;
            }

            t = tr[t].d0;
        }


        else if ((tr[t].d0 <= 0) && (tr[t].d1 > 0))
        {			/* Only one trapezoid below */
            if ((tr[t].u0 > 0) && (tr[t].u1 > 0))
            {			/* continuation of a chain from abv. */
                if (tr[t].usave > 0) /* three upper neighbours */
                {
                    if (tr[t].uside == 1)
                    {
                        tr[tn].u0 = tr[t].u1;
                        tr[t].u1 = -1;
                        tr[tn].u1 = tr[t].usave;

                        tr[tr[t].u0].d0 = t;
                        tr[tr[tn].u0].d0 = tn;
                        tr[tr[tn].u1].d0 = tn;
                    }
                    else		/* intersects in the right */
                    {
                        tr[tn].u1 = -1;
                        tr[tn].u0 = tr[t].u1;
                        tr[t].u1 = tr[t].u0;
                        tr[t].u0 = tr[t].usave;

                        tr[tr[t].u0].d0 = t;
                        tr[tr[t].u1].d0 = t;
                        tr[tr[tn].u0].d0 = tn;
                    }

                    tr[t].usave = tr[tn].usave = 0;
                }
                else		/* No usave.... simple case */
                {
                    tr[tn].u0 = tr[t].u1;
                    tr[t].u1 = tr[tn].u1 = -1;
                    tr[tr[tn].u0].d0 = tn;
                }
            }
            else
            {			/* fresh seg. or upward cusp */
                int tmp_u = tr[t].u0;
                int td0, td1;
                if (((td0 = tr[tmp_u].d0) > 0) &&
                    ((td1 = tr[tmp_u].d1) > 0))
                {		/* upward cusp */
                    if ((tr[td0].rseg > 0) &&
                        !i1_of(tr[td0].rseg, &s.v1))
                    {
                        tr[t].u0 = tr[t].u1 = tr[tn].u1 = -1;
                        tr[tr[tn].u0].d1 = tn;
                    }
                    else
                    {
                        tr[tn].u0 = tr[tn].u1 = tr[t].u1 = -1;
                        tr[tr[t].u0].d0 = t;
                    }
                }
                else		/* fresh segment */
                {
                    tr[tr[t].u0].d0 = t;
                    tr[tr[t].u0].d1 = tn;
                }
            }

            if (is_abs_less_than_col_tol(tr[t].lo.y, tr[tlast].lo.y) &&
                is_abs_less_than_col_tol(tr[t].lo.x, tr[tlast].lo.x) && tribot)
            {		/* bottom forms a triangle */
                int tmpseg;

                if (is_swapped)
                    tmptriseg = seg[segnum].prev;
                else
                    tmptriseg = seg[segnum].next;

                if ((tmpseg > 0) && i1_of(tmpseg, &s.v0))
                {
                    /* L-R downward cusp */
                    tr[tr[t].d1].u0 = t;
                    tr[tn].d0 = tr[tn].d1 = -1;
                }
                else
                {
                    /* R-L downward cusp */
                    tr[tr[tn].d1].u1 = tn;
                    tr[t].d0 = tr[t].d1 = -1;
                }
            }
            else
            {
                if ((tr[tr[t].d1].u0 > 0) && (tr[tr[t].d1].u1 > 0))
                {
                    if (tr[tr[t].d1].u0 == t) /* passes thru LHS */
                    {
                        tr[tr[t].d1].usave = tr[tr[t].d1].u1;
                        tr[tr[t].d1].uside = 1;
                    }
                    else
                    {
                        tr[tr[t].d1].usave = tr[tr[t].d1].u0;
                        tr[tr[t].d1].uside = 2;
                    }
                }
                tr[tr[t].d1].u0 = t;
                tr[tr[t].d1].u1 = tn;
            }

            t = tr[t].d1;
        }

            /* two trapezoids below. Find out which one is intersected by */
            /* this segment and proceed down that one */

        else
        {
            int tmpseg = tr[tr[t].d0].rseg;
            double y0, yt;
            point_coordinates tmppt;
            int tnext, i_d0, i_d1;

            i_d0 = i_d1 = 0;
            if (is_abs_less_than_col_tol(tr[t].lo.y, s.v0.y))
            {
                if (tr[t].lo.x > s.v0.x)
                    i_d0 = 1;
                else
                    i_d1 = 1;
            }
            else
            {
                tmppt.y = y0 = tr[t].lo.y;
                yt = (y0 - s.v0.y)/(s.v1.y - s.v0.y);
                tmppt.x = s.v0.x + yt * (s.v1.x - s.v0.x);

                if (less_than(&tmppt, &tr[t].lo))
                    i_d0 = 1;
                else
                    i_d1 = 1;
            }

            /* check continuity from the top so that the lower-neighbour */
            /* values are properly filled for the upper trapezoid */

            if ((tr[t].u0 > 0) && (tr[t].u1 > 0))
            {			/* continuation of a chain from abv. */
                if (tr[t].usave > 0) /* three upper neighbours */
                {
                    if (tr[t].uside == 1)
                    {
                        tr[tn].u0 = tr[t].u1;
                        tr[t].u1 = -1;
                        tr[tn].u1 = tr[t].usave;

                        tr[tr[t].u0].d0 = t;
                        tr[tr[tn].u0].d0 = tn;
                        tr[tr[tn].u1].d0 = tn;
                    }
                    else		/* intersects in the right */
                    {
                        tr[tn].u1 = -1;
                        tr[tn].u0 = tr[t].u1;
                        tr[t].u1 = tr[t].u0;
                        tr[t].u0 = tr[t].usave;

                        tr[tr[t].u0].d0 = t;
                        tr[tr[t].u1].d0 = t;
                        tr[tr[tn].u0].d0 = tn;
                    }

                    tr[t].usave = tr[tn].usave = 0;
                }
                else		/* No usave.... simple case */
                {
                    tr[tn].u0 = tr[t].u1;
                    tr[tn].u1 = -1;
                    tr[t].u1 = -1;
                    tr[tr[tn].u0].d0 = tn;
                }
            }
            else
            {			/* fresh seg. or upward cusp */
                int tmp_u = tr[t].u0;
                int td0, td1;
                if (((td0 = tr[tmp_u].d0) > 0) &&
                    ((td1 = tr[tmp_u].d1) > 0))
                {		/* upward cusp */
                    if ((tr[td0].rseg > 0) &&
                        !i1_of(tr[td0].rseg, &s.v1))
                    {
                        tr[t].u0 = tr[t].u1 = tr[tn].u1 = -1;
                        tr[tr[tn].u0].d1 = tn;
                    }
                    else
                    {
                        tr[tn].u0 = tr[tn].u1 = tr[t].u1 = -1;
                        tr[tr[t].u0].d0 = t;
                    }
                }
                else		/* fresh segment */
                {
                    tr[tr[t].u0].d0 = t;
                    tr[tr[t].u0].d1 = tn;
                }
            }

            if (is_abs_less_than_col_tol(tr[t].lo.y, tr[tlast].lo.y) &&
                is_abs_less_than_col_tol(tr[t].lo.x, tr[tlast].lo.x) && tribot)
            {
                /* this case arises only at the lowest trapezoid.. i.e.
               tlast, if the lower endpoint of the segment is
               already inserted in the structure */

                tr[tr[t].d0].u0 = t;
                tr[tr[t].d0].u1 = -1;
                tr[tr[t].d1].u0 = tn;
                tr[tr[t].d1].u1 = -1;

                tr[tn].d0 = tr[t].d1;
                tr[t].d1 = tr[tn].d1 = -1;

                tnext = tr[t].d1;
            }
            else if (i_d0)
                /* intersecting d0 */
            {
                tr[tr[t].d0].u0 = t;
                tr[tr[t].d0].u1 = tn;
                tr[tr[t].d1].u0 = tn;
                tr[tr[t].d1].u1 = -1;

                /* new code to determine the bottom neighbours of the */
                /* newly partitioned trapezoid */

                tr[t].d1 = -1;

                tnext = tr[t].d0;
            }
            else			/* intersecting d1 */
            {
                tr[tr[t].d0].u0 = t;
                tr[tr[t].d0].u1 = -1;
                tr[tr[t].d1].u0 = t;
                tr[tr[t].d1].u1 = tn;

                /* new code to determine the bottom neighbours of the */
                /* newly partitioned trapezoid */

                tr[tn].d0 = tr[t].d1;
                tr[tn].d1 = -1;

                tnext = tr[t].d1;
            }

            t = tnext;
        }

        tr[t_sav].rseg = tr[tn_sav].lseg  = segnum;
    } /* end-while */

    /* Now combine those trapezoids which share common segments. We can */
    /* use the pointers to the parent to connect these together. This */
    /* works only because all these new trapezoids have been formed */
    /* due to splitting by the segment, and hence have only one parent */

    tfirstl = tfirst;
    tlastl = tlast;
    merge_trapezoids(segnum, tfirstl, tlastl, 1);
    merge_trapezoids(segnum, tfirstr, tlastr, 2);

    seg[segnum].is_inserted = 1;
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Update the roots stored for each of the endpoints of the segment.
 * This is done to speed up the location-query for the endpoint when
 * the segment is inserted into the trapezoidation subsequently
 */
int polygon_triangulation::find_new_roots(int segnum) {
    segment_struct *s = &seg[segnum];

    if (s->is_inserted)
        return 0;

    s->root0 = locate_endpoint(&s->v0, &s->v1, s->root0);
    s->root0 = tr[s->root0].sink;

    s->root1 = locate_endpoint(&s->v1, &s->v0, s->root1);
    s->root1 = tr[s->root1].sink;
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Main routine to perform trapezoidation */
int polygon_triangulation::construct_trapezoids(int nseg){
    int i;
    int root, h;
/* Add the first segment and get the query structure and trapezoid */
/* list initialised */
/* the next segment in the generated random ordering of all the segments in S */
    root = init_query_structure(permute[choose_index++]);

    for (i = 1; i <= nseg; i++)
        seg[i].root0 = seg[i].root1 = root;

    for (h = 1; h <= math_logstar_n(nseg); h++) {
        for (i = math_N(nseg, h - 1) + 1; i <= math_N(nseg, h); i++)
            /*the next segment in the generated random ordering of all the segments in S */
            add_segment(permute[choose_index++]);

/* Find a new root for each of the segment endpoints */
        for (i = 1; i <= nseg; i++)
            find_new_roots(i);
    }

    for (i = math_N(nseg, math_logstar_n(nseg)) + 1; i <= nseg; i++)
        /*the next segment in the generated random ordering of all the segments in S */
        add_segment(permute[choose_index++]);

    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int polygon_triangulation::initialise(int n) {
    int i;

    for (i = 1; i <= n; i++)
        seg[i].is_inserted = 0;

    choose_index = 1;
    for (i = 1; i <= n; i++) {
        permute[i] = i;
    }
    return 0;
}

/* Input specified as contours. *
 * Every contour is specified by giving all its points in order. No
 * point shoud be repeated. i.e. if the outer contour is a square,
 * only the four distinct endpoints shopudl be specified in order.
 * ncontours: #contours
 * cntr: An array describing the number of points in each
 *	 contour. Thus, cntr[i] = #points in the i'th contour.
 * vertices: Input array of vertices. Vertices for each contour
 *           immediately follow those for previous one. Array location
 *           vertices[0] must NOT be used (i.e. i/p starts from
 *           vertices[1] instead. The output triangles are
 *	     specified  w.r.t. the indices of these vertices.
 * triangles: Output array to hold triangles.
 *
 * Enough space must be allocated for all the arrays before calling
 * this routine
 */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int polygon_triangulation::triangulate_polygon(int n, int ncontours, const int cntr[], double (*vertice)[2], int (*triangles)[3]) {
    int i;
    int nmonpoly, ccount, npoints, genus;
    auto ** vertices = (point_coordinates **) malloc((n+1)*sizeof(point_coordinates **));
    for (int j = 0; j < n+1; ++j) {
        vertices[j] = new point_coordinates;
    }
    vertices[0]={};
    vertices[1]->x =vertice[0][0];
    vertices[1]->y =vertice[0][1];
    for (int j = 1; j < n; ++j) {
        vertices[n-j+1]->x =vertice[j][0];
        vertices[n-j+1]->y =vertice[j][1];
    }
    memset((void *) seg, 0, sizeof(seg));
    ccount = 0;
    i = 1;

    while (ccount < ncontours) {
        int j;
        int first, last;

        npoints = cntr[ccount];
        first = i;
        last = first + npoints - 1;
        for (j = 0; j < npoints; j++, i++) {
            seg[i].v0.x = vertices[i]->x;
            seg[i].v0.y = vertices[i]->y;

            if (i == last) {
                seg[i].next = first;
                seg[i].prev = i - 1;
                seg[i - 1].v1 = seg[i].v0;
            } else if (i == first) {
                seg[i].next = i + 1;
                seg[i].prev = last;
                seg[last].v1 = seg[i].v0;
            } else {
                seg[i].prev = i - 1;
                seg[i].next = i + 1;
                seg[i - 1].v1 = seg[i].v0;
            }

            seg[i].is_inserted = 0;
        }

        ccount++;
    }

    genus = ncontours - 1;
    n = i - 1;

    initialise(n);
    construct_trapezoids(n);
    nmonpoly = monotonate_trapezoids(n);
    int number_of_triangles;
    number_of_triangles = triangulate_monotone_polygons(n, nmonpoly, triangles);

    return number_of_triangles;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//WARNING!!!!!!
// I DON'T KNOW DUE TO SOME MINOR ISSUE, THE FIRST TIME YOU RUN IT WITH NEW INPUT YOU MAY OR MAY NOT GET ERROR
// PLEASE RUN AGAIN TO GET CORRECT ANSWER

int main() {
    ////////////////////////////////////////////////////////////
    //change n and data in order to put different values
    //DON'T FORGET TO CHANGE THE VALUE FOR N
    int n=10, nmonpoly, npoints, first, last, op[segment_size][3], i, number_of_triangles;
    double data[10][2]={{0,0},//CHANGE THIS TO GIVE DIFFERENT INPUTS // GIVE INPUT IN CLOCKWISE ORDER
                       {0,5},
                       {1,2},
                       {2,5},
                       {3,2},
                       {4,5},
                       {4,0},
                       {3,1},
                       {2,0},
                       {1,1}};
    ///////////////////////////////////////////////////////////////////
    //this is the number of contours
    int cnt[1] = {n};
    //here I have declared the object of the class
    polygon_triangulation obj;
    number_of_triangles = obj.triangulate_polygon(10,1,cnt,data,op);
    //create a container for output
    auto ** output = (point_coordinates **) malloc(2*(n-3)*sizeof(point_coordinates **));
    for (int j = 0; j < 2*(n-3); ++j) {
        output[j] = new point_coordinates;
        output[j]->x =0;
        output[j]->y = 0;
    }
    int size=0;
    //add diagonals from triangles
    for (i = 0; i < number_of_triangles; i++) {
        if(i==number_of_triangles-1) {
            if (abs(op[i][0] - op[i][1]) % (n - 1) > 1 && abs(op[i][1] - op[i][2]) % (n - 1) > 1) {
                output[size]->x = op[i][0] - 1;
                output[size++]->y = op[i][1] - 1;
                output[size]->x = op[i][1] - 1;
                output[size++]->y = op[i][2] - 1;
            }
            else if (abs(op[i][0] - op[i][1]) % (n - 1) > 1){
                output[size]->x = op[i][0] - 1;
                output[size++]->y = op[i][1] - 1;
            }
            else if (abs(op[i][1] - op[i][2]) % (n - 1) > 1){
                output[size]->x = op[i][1] - 1;
                output[size++]->y = op[i][2] - 1;
            }
        }
        else{
            for (int k = 0; k < 2; ++k) {
                if (abs(op[i][k] - op[i][(k + 1)]) % (n - 1) > 1){
                    output[size]->x = op[i][k] - 1;
                    output[size++]->y = op[i][k+1] - 1;
                }
            }
        }
    }
    //remove duplicates
    for(i=0;i<size;i++){
        for(int j=i;j < size;j++){
            if((output[i]->x==output[j]->y && output[j]->x==output[i]->y)){
                output[i]->x = -1;
                output[i]->y = -1;
                size--;
            }
        }
    }
    //output the final values
    std::cout<<"{";
    for(int k=0;k<2*(n-3);k++){
        //std::cout << output[k]->x << " " << output[k]->y<< std::endl;//for debugging code
        if(output[k]->x != output[k]->y ) {
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

////////////////////////////////////////////////////END////////////////////////////////////////////////////////////////////





