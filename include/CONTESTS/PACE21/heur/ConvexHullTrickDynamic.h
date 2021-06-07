//
// Created by sylwester on 4/1/21.
//

#ifndef ALGORITHMSPROJECT_CONVEXHULLTRICKDYNAMIC_H
#define ALGORITHMSPROJECT_CONVEXHULLTRICKDYNAMIC_H

#include "Makros.h"

//*********************************************************************************************************
namespace HullLine {
    using _type = long double; // it can be changed to long long if coefficients are always integers

    const _type is_query = 1.0 * -100'000'000 * 100'000'000;

    struct Line {
        _type m, b;
        int id;
        mutable function<const Line *()> succ;

        bool operator<(const Line &rhs) const {
            if (rhs.b != is_query) return m < rhs.m;
            const Line *s = succ();
            if (!s) return 0;
            _type x = rhs.m;
            return b - s->b < (s->m - m) * x;
        }
    };
}


/**
 * Class responsible for dynamically adding lines and checking for MAXIMUM in given point
 */
struct HullDynamic : public multiset<HullLine::Line> { // will maintain upper hull for maximum

    bool bad(iterator y) {
        auto z = next(y);
        if (y == begin()) {
            if (z == end()) return 0;
            return y->m == z->m && y->b <= z->b;
        }
        auto x = prev(y);
        if (z == end()) return y->m == x->m && y->b <= x->b;

        // **** May need long double typecasting here
        return (long double)(x->b - y->b)*(z->m - y->m) >= (long double)(y->b - z->b)*(y->m - x->m);
    }

    /**
     * Inserts line    y = mx + b
     */
    void insert_line(HullLine::_type m, HullLine::_type b, int id) { // version for doubles
        auto y = insert({ m, b, id });
        y->succ = [=] { return next(y) == end() ? 0 : &*next(y); };
        if (bad(y)) { erase(y); return; }
        while (next(y) != end() && bad(next(y))) erase(next(y));
        while (y != begin() && bad(prev(y))) erase(prev(y));
    }

    /**
     *
     * @param real
     * @return maximum value in point x over all added functions.
     */
    HullLine::_type eval(HullLine::_type x) {
        auto l = *lower_bound((HullLine::Line) { x, HullLine::is_query, -1 });
        return l.m * x + l.b;
    }

    /**
     * Finds and returns a line that maximizes  y = mx + b
     * @param x
     */
    HullLine::Line get_line(HullLine::_type x) {
        auto l = *lower_bound((HullLine::Line) { x, HullLine::is_query, -1 });
        return l;
    }
};

//************************************************************************************************************

class ConvexHullTrickDynamic{
public:

    /**
     * If [max] is true, then we will maximize added functions, otherwise they will be minimized. By default they are
     * minimized.
     * Coefficients can be real (double) value
     */
    ConvexHullTrickDynamic( const bool max = false ) : maximize(max){}

    /**
     * Adds line mx + b, with given id.
     * If [maximize] is false, then in fact line -mx-b is added to [hull].
     */
    void insert_line( HullLine::_type m, HullLine::_type b, int id ){
        if( maximize ) hull.insert_line(m,b,id);
        else hull.insert_line(-m,-b,id);
    }

    /**
     * Finds and returns id of the line that minimizes/maximizes value mx + b
     */
    int get_line_id(HullLine::_type x){
        return hull.get_line(x).id;
    }

    void clear(){ hull.clear(); }

    /**
     *
     * @return number of lines that are currently in the hull. It may be less then the number of lines added, if some
     * were removed.
     */
    int size(){ return hull.size(); }

private:

    HullDynamic hull;
    const bool maximize;
};

#endif //ALGORITHMSPROJECT_CONVEXHULLTRICKDYNAMIC_H
