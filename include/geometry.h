#ifndef GEOMETRY_H
#define GEOMETRY_H

#define _USE_MATH_DEFINES

#include <cstdio>
#include <cwchar>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>


#define EPS 1.e-5

#define sign(a) ((a == 0) ? 0.0 : (a < 0 ? -1.0 : 1.0))

#define sqr(a) ((a)*(a))

#define MinKM 1000 // meters in kilometers
#define KM2M(a) (a * MinKM)
#define M2KM(a) (a / MinKM)

struct point
{
	double x;
	double y;

	point() : x(0), y(0) {};
	point(double a, double b) : x(a), y(b) {};
	point(const point &p) : x(p.x), y(p.y) {};

	bool equal(const point &p);

	bool equal_eps(const point &p);

	double distance_to(point p) const;

	void to_related_cs(const point &origin);
	void to_global_cs(const point &origin);

	void to_dec_cs(const point &origin);
	void to_geo_cs(const point &origin);

	point at_geo_cs(const point &origin) const;
	point at_dec_cs(const point &origin) const;

	std::string toString(std::string name = std::string("")) const;

};

//todo если создать vec из двух point, а затем последние изменить, то изменится ли vec?
struct vec
{
	point start;
	point end;

	vec(point a, point b) : start(a), end(b) {};

	bool contains(point p);

	double length();

	void shorten(double len/*, double dl = 0.0*/);

	point middle();

	std::string toString(std::string name = std::string(""));

	vec at_geo_cs(const point &origin);

	std::string toGlanceFormat();
};

class Line // line segment
{
	double a, b, c;
	point p1, p2;

	double f(point p);
public:
	Line();
	Line(vec v);
	Line(point _p1, point _p2);
	Line(double _a, double _b, double _c);

	double get_a();
	double get_b();
	double get_c();

	point start();
	point end();

	vec interval();
	void set_interval(point _p1, point _p2);

	bool contains(point p);

	double distance_to(point p);
	double angle(vec v);

	point intersection(Line l);
	point projection_of(point p);

	vec perpendicular(point p);
	Line parallel(point p);

	std::string as_string(std::string name = std::string("")) const;
};

#endif // GEOMETRY_H