#ifndef _link_diagram_
#define _link_diagram_

/*
 ********************************
 *                              *
 *                              *
 *       link_diagram.h         *
 *                              *
 *                              *
 ********************************
 */

#include <vector>
#include <string>

struct Pass;
struct Crossing;

enum Level  { over, under };


class Link
{

public:
    Link( const std::string& s );
    ~Link();

    int millett[100][10], millett_component[100][2];
    int crossing_number, num_components;
    std::vector<int> component_length, component;
    
    bool make_geometry();
    void make_millett();


private:
    
    std::string link_record;
    
    std::vector<int> comp_start, comp_end;
    std::vector<int> dt;
    std::vector<int> dowker;
    std::vector<int> crossing_sign;           // left-right sign
    std::vector<int> alt_sign;                // over-under sign
    std::vector<int> next_index, prev_index;

    std::vector<Crossing *> crossingList;     // crossings of diagram
    std::vector<Pass *> passList;
    
    bool compute_crossing_signs();       // implement DT algorithm
    void clean_up_diagram();
};


struct Pass
{
    int index, new_index, crossing_sign, orientation, component, alt_sign;
    Crossing *crossing;
    Pass *next;
    Pass *prev;
    Pass *partner;
    Level level;
};


struct Crossing
{
    ~Crossing();
    int index;
    int sign, alt_sign, orientation;
    Pass *odd, *even, *over, *under;
};

#endif

