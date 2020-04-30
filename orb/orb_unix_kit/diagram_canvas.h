#ifndef DIAGRAM_CANVAS_H
#define DIAGRAM_CANVAS_H

#include "diagram_to_trig.h"

#include <iostream>
//#include <qwidget.h>
//#include <qpainter.h>
//#include <qpixmap.h>
//#include <qpoint.h>
//#include <qevent.h>
#include <complex>
#include <string>
#include <vector>
#include <set>
//#include <qcolor.h>
//#include <qtextstream.h>
#include "graph_complement.h"
//#include "color.h"
//using namespace std;

#include <sstream>

class QPoint
{
public:
    int x() const { return _x; }
    int y() const { return _y; }
    
    QPoint(int x = 0, int y = 0) : _x(x), _y(y) { }

    QPoint operator-(const QPoint &other) const
    {
        return QPoint(_x - other._x, _y - other._y);
    }

    QPoint operator+(const QPoint &other) const
    {
        return QPoint(_x + other._x, _y + other._y);
    }

private:
    int _x, _y;
};

inline QPoint operator*(const double s, const QPoint &other)
{
    return QPoint(s * other.x(), s * other.y());
}

typedef std::stringstream QTextStream;

//class QTextStream { };
//public: 
//    QTextStream &operator<<(const char *) { return *this;}
//    QTextStream &operator<<(int) { return *this;}
//
//    QTextStream &operator>>(int &) { return *this;}
//    QTextStream &operator>>(double &) { return *this;}
//};

struct EndData;
struct Edge;
struct Vertex;
struct Crossing;


enum JoinState  { no_join, join };
enum EndType  {begin = 0, end = 1};
enum EdgeType {singular = 0, drilled = 1};

#define	THICK 3
#define	THIN 1 
#define ARROW 1 

struct EndData
{
    EndData( Edge *edge, EndType type );
    ~EndData() {;} 
    Edge *edge;
    EndType type;
    bool singular;
    double angle;
};

struct Vertex
{
    Vertex( QPoint point);
    ~Vertex();
    QPoint position;
    int connected_component;
    int vertex_id;
    int link_id;
    std::vector<EndData *> incidentEndData;
};

struct Crossing
{
    QPoint position;
    int crossing_id;
    int crossing_sign;
    Edge *over, *under;
    double position_on_overstrand, position_on_understrand;
    void switchCrossing();
};


struct Edge
{
    Edge( Vertex *vertex1, Vertex *vertex2 );
    ~Edge();
    Vertex *vertex[2];
    std::set<double> underpasses;
    std::vector<Crossing *> crossings; 
    int arc_id;
    int link_id;
    int cuff_id;
    int arc;
    int edge_id;
    int thickness;
    bool selected;
    int a, b, c;
    EdgeType edge_type;
    void compute_equation_coefficients();
    double distance_from_edge( QPoint point );
    double projection_to_edge( QPoint point );
};

class DiagramCanvas// : public QWidget 
{
//    Q_OBJECT

public:
    DiagramCanvas( /*QWidget *parent,*/ const char *name = 0, bool readOnly = TRUE );
    ~DiagramCanvas();
    void readDiagram( QTextStream &stream);
    void saveDiagram( QTextStream &stream);
    void highlightEdge( int e );
    void setReadOnly();
    void deselectEdges();
    bool isUntouched();
    void deleteArc( int a );
    void drillArc( int a );

//public slots:
	void drillToggle();
	void deleteSelectedEdges();
	void clearDiagram();
	void exportSlot();

private:
    int num_arcs;
    int num_links;
    int num_cuffs;
    bool edge_creation_in_progress;
    bool vertex_drag_in_progress;
    JoinState joinState;
    bool untouched;
    bool drill_on;
    //   QColor colorList[21];
    //   QPoint saved_vertex_position;
    std::vector<Crossing *> crossingList;
    std::vector<Edge *> edgeList;
    std::vector<Vertex *> vertexList;
    Crossing *currentCrossing;    
    Edge *currentEdge;
    Vertex *draggedVertex;
    bool		readOnly; 
//    QPixmap		buffer;

public:
    Triangulation *outputTriangulation();
private:
    bool proximity( QPoint point1, QPoint point2 );
    bool proximity( QPoint point, Edge *edge );
    bool new_point_general_position( QPoint point, int *edge_index );
    bool completed_edge_general_position( Edge *edge, int *edge_index );
    bool dragged_vertex_general_position();
    bool edge_intersect(Edge *edge1, Edge *edge2, double *t1, double *t2, QPoint *point);
    void assign_new_crossings(Edge *edge);
    void refresh_crossing( Crossing *c );
    void assign_crossings_to_edges();
    void assign_links();
    void reassign_links( bool arc_was_drilled );
    void assign_arcs();
    void reassign_arcs( int a);
    void assign_cuffs();
    void delete_dragged_crossings();
    double proximity_tolerance;
    void ed_angles();
    void getCrossingSigns();
    void invert_arc( Edge *e );
    void prepare_components_for_output();
//    void mousePressEvent( QMouseEvent *e );
//    void mouseReleaseEvent( QMouseEvent *e );
//    void mouseMoveEvent( QMouseEvent *e );
//    void resizeEvent( QResizeEvent * );
//    void paintEvent( QPaintEvent * );
//    void paint_all( QRect& r);
//    void paint( Edge *edge, bool paint_over );
//    void redrawMovedEdge( QMouseEvent *e );
//    void redrawDraggedVertex( QMouseEvent *e );
    void initialize_colors();
    bool isReadyForTriangulation();

//protected slots:

//signals:
    void drillToggled();
//    void statusBarChanged(const QString &message);
    void exporting(Triangulation *manifold ); 
    void selectionChanged( int arc_id );

};


#endif
