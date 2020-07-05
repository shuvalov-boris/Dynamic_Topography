#ifndef GEOMETRY_H
#define GEOMETRY_H

#define _USE_MATH_DEFINES

#include <cstdio>
#include <cwchar>
#include <fstream>
#include <string>
#include <math.h>

#define EPS 0.0000001

#define sign(a) (a == 0) ? 0.0 : (a < 0 ? -1.0 : 1.0)

using namespace std; //	delete (escape it)

struct point
{
	double x;
	double y;

	point() : x(0), y(0) {};
	point(double a, double b) : x(a), y(b) {};
	// point(const point &p) : x(p.x), y(p.y) {};

	bool equal(const point &p)
	{
		return (x == p.x && y == p.y);
	}

	double distance_to(point p)
	{
		return sqrt((x - p.x) * (x - p.x) + (y - p.y) * (y - p.y));
	}

	string toString(string name = string(""))
	{
		char str[30];
		sprintf(str, "%s(%2f, %2f)", name.c_str(), x, y);
		return string(str);
	}

};

struct vec
{
	point start;
	point end;

	vec(point a, point b) : start(a), end(b) {};

	bool contains(point p)
	{
		return sign(p.x - start.x) * sign(p.x - end.x) < 0;
	}

	double length()
	{
		return sqrt((start.x - end.x) * (start.x - end.x) + (start.y - end.y) * (start.y - end.y));
	}

/*	double sign_vv(vec ln, vec v)
	{
		double l1 = dist_vp(ln, v.start), l2 = dist_vp(ln, v.end);
		// if ((p2v(ln, v.start) < 0 && p2v(ln, v.end) > 0) || sign(p2v(ln, v.start)) * sign(l2 - l1) >= 0) return 1.0;
		if ((p2v(ln, v.start) < 0 && (p2v(ln, v.end) > 0 || (p2v(ln, v.end) < 0 && l2 - l1 > 0)))) return 1.0;
		if (p2v(ln, v.start) > 0 && l2 - l1 > 0) return 1.0;
		return -1.0;
	}*/

	void shorten(double len/*, double dl = 0.0*/)
	{
		double m = len, n;
		/*if (dl == 0.0) */n = start.distance_to(end) - len;
		/*else n = dl;*/
		double x1 = start.x, y1 = start.y, x2 = end.x, y2 = end.y;
		double l = m / n;
		double x = (x1 + l * x2) / (1 + l),
			   y = (y1 + l * y2) / (1 + l);
		end = point(x, y);
		return;
	}

	point middle()
	{
		return point((start.x + end.x) / 2, (start.y + end.y) / 2);
	}

	string toString(string name = string(""))
	{
		char str[100];
		sprintf(str, "%s[%s -> %s]", name.c_str(), start.toString().c_str(), end.toString().c_str());
		return string(str);
	}

	string toGlanceFormat()
	{
		char str[100];
		sprintf(str, "TYPE = VECTOR\tCOLOR = 14\tWIDTH = 1\tSCALE = 1.00\tGEO = (%2f,%2f %2f,%2f)\n",
					 start.x, start.y, end.x, end.y);
		return string(str);
	}
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
	point projection(point p);

	vec perpendicular(point p);
	Line parallel(point p);
};

Line::Line() : a(0), b(0), c(0), p1(point(0, 0)), p2(point(0, 0)) {}

Line::Line(vec v) : p1(v.start), p2(v.end)
{
	a = p2.y - p1.y;
	b = p1.x - p2.x;
	c = - p1.y * b - p1.x * a;
}

Line::Line(point _p1, point _p2) : p1(_p1), p2(_p2)
{
	a = p2.y - p1.y;
	b = p1.x - p2.x;
	c = - p1.y * b - p1.x * a;
}

Line::Line(double _a, double _b, double _c) : a(_a), b(_b), c(_c), p1(point(0, 0)), p2(point(0, 0))
{
	// p1 and p2 are NOT INITIALIZED
}

double Line::f(point p)
{
	return a * p.x + b * p.y + c;
}

double Line::get_a()
{
	return a;
}

double Line::get_b()
{
	return b;
}

double Line::get_c()
{
	return c;
}

point Line::start()
{
	return p1;
}

point Line::end()
{
	return p2;
}

vec Line::interval()
{
	return vec(p1, p2);
}

void Line::set_interval(point _p1, point _p2)
{
	// Line(_p1, _p2);
	p1 = _p1; p2 = _p2;
	a = p2.y - p1.y;
	b = p1.x - p2.x;
	c = - p1.y * b - p1.x * a;
}

bool is_t(double numerator, double denominator)
{
	return (denominator == 0 && numerator == 0) || 
		   (denominator != 0 && numerator / denominator >= 0 && numerator / denominator <= 1);
}

bool Line::contains(point p)
{
	return fabs(f(p)) < EPS && is_t(p.x - p1.x, p2.x - p1.x) && is_t(p.y - p1.y, p2.y - p1.y);
}


double Line::distance_to(point p)
{
	return f(p) / sqrt(a * a + b * b);
}

double Line::angle(vec v)
{
	Line l(v);
	return atan2(a * l.get_b() - l.get_a() * b, a * l.get_a() + b * l.get_b()) * 180 / M_PI;
}


point Line::intersection(Line l)
{
	double y = (a * l.get_c() - l.get_a() * c) / (l.get_a() * b - a * l.get_b());
	double x = (b * l.get_c() - l.get_b() * c) / (a * l.get_b() - l.get_a() * b);
	return point(x, y);
}

point Line::projection(point p)
{
	Line l = perpendicular(p);
	return intersection(l);
}

vec Line::perpendicular(point p)
{
	return vec(p, point(p.x + a, p.y + b));
}

Line Line::parallel(point p)
{
	return Line(get_a(), get_b(), -(get_a() * p.x + get_b() * p.y));
}
/*
double sign(double arg);
double sign_vec(double arg);
void set_kof(double &a, double &b, double &c, double x1, double y1, double x2, double y2);
void calc_kof(double &a, double &b, double &c, vec v);
double dist_lp(double a, double b, double c, point p);
double dist_pp(point a, point b);
double dist_vp(vec v, point p);
double dist_proj(point a, point b);
double p2v(vec v, point p);
double sign_vv(vec ln, vec v);
point get_cross_ll(double a1, double b1, double c1, double a2, double b2, double c2);
point get_projection(double a1, double b1, double c1, point p);
point get_perpend(point p, vec v);
point get_midpoint(point a, point b);
bool into(point p1, point p2, point p);
double angle(double a1, double b1, double a2, double b2);
double angle_p4(point p1, point p2, point p3, point p4);
string outp(point p, string name);
string outv(vec v, string name);
point short_vec(vec v, double len);
point short_v(vec v, double len, double dl);
string str_vec_format(vec v);



double sign(double arg)
{
	if (arg >  0) return  1.0;
	if (arg <  0) return -1.0;
	if (arg == 0) return  0.0;
}

double sign_vec(double arg)
{
	if (arg < 90.0 && arg > -90.0) return 1.0;
	else return -1.0;
}

void set_kof(double &a, double &b, double &c, double x1, double y1, double x2, double y2)
{
	a = y2 - y1;
	b = x1 - x2;
	c = y1 * (x2 - x1) - x1 * (y2 - y1);
}

void calc_kof(double &a, double &b, double &c, vec v)
{
	set_kof(a, b, c, v.start.x, v.start.y, v.end.x, v.end.y);
}

double dist_lp(double a, double b, double c, point p)
{
	return (a * p.x + b * p.y + c) / sqrt(a * a + b * b);
}

double dist_pp(point a, point b)
{
	return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

double dist_vp(vec v, point p)
{
	double a, b, c;
	calc_kof(a, b, c, v);
	return dist_lp(a, b, c, p);
}

double dist_proj(point a, point b)
{
	return dist_pp(a, b) * sign(b.x - a.x);
}

double p2v(vec v, point p)
{
	double a, b, c;
	calc_kof(a, b, c, v);
	return a * p.x + b * p.y + c;
}

double sign_vv(vec ln, vec v)
{
	double l1 = dist_vp(ln, v.start), l2 = dist_vp(ln, v.end);
	// if ((p2v(ln, v.start) < 0 && p2v(ln, v.end) > 0) || sign(p2v(ln, v.start)) * sign(l2 - l1) >= 0) return 1.0;
	if ((p2v(ln, v.start) < 0 && (p2v(ln, v.end) > 0 || (p2v(ln, v.end) < 0 && l2 - l1 > 0)))) return 1.0;
	if (p2v(ln, v.start) > 0 && l2 - l1 > 0) return 1.0;
	return -1.0;
}

point get_cross_ll(double a1, double b1, double c1, double a2, double b2, double c2)
{
	double y = (a1 * c2 - a2 * c1) / (a2 * b1 - a1 * b2);
	double x = (b1 * c2 - b2 * c1) / (a1 * b2 - a2 * b1);
	return point(x, y);
}

point get_projection(double a1, double b1, double c1, point p)
{
	point p2(p.x + a1, p.y + b1);
	double a2, b2, c2;
	calc_kof(a2, b2, c2, vec(p, p2));
	return get_cross_ll(a1, b1, c1, a2, b2, c2);
}

point get_perpend(point p, vec v)
{
	double a, b, c;
	calc_kof(a, b, c, v);
	return point(p.x + a, p.y + b);
}

point get_midpoint(point a, point b)
{
	return point((a.x + b.x) / 2, (a.y + b.y) / 2);
}

bool into(point p1, point p2, point p)
{
	return sign(p.x - p1.x) * sign(p.x - p2.x) < 0;
}

double angle(double a1, double b1, double a2, double b2)
{
	return atan((a1 * b2 - a2 * b1) / (a1 * a2 + b1 * b2)) * 180 / M_PI;
}

double angle_p4(point p1, point p2, point p3, point p4)
{
	double a1, b1, c1;
	calc_kof(a1, b1, c1, vec(p1, p2));

	double a2, b2, c2;
	calc_kof(a2, b2, c2, vec(p3, p4));

	return angle(a1, b1, a2, b2);
}

string outp(point p, string name = string(""))
{
	char str[30];
	sprintf(str, "%s(%2f, %2f)", name.c_str(), p.x, p.y);
	return string(str);
}

string outv(vec v, string name = string(""))
{
	char str[100];
	sprintf(str, "%s[%s -> %s]", name.c_str(), outp(v.start).c_str(), outp(v.end).c_str());
	return string(str);
}

point short_vec(vec v, double len)
{
	double fi = angle_p4(v.start, v.end, v.start, point(v.start.x + 10, v.start.y));
	return point(v.start.x + len * cos(fi), v.start.y + len * sin(fi));
}

point short_v(vec v, double len, double dl = 0.0)
{
	double m = len, n;
	if (dl == 0.0) n = dist_pp(v.start, v.end) - len;
	else n = dl;
	double x1 = v.start.x, y1 = v.start.y, x2 = v.end.x, y2 = v.end.y;
	// double x = (n * x1 + m * x2) / (n + m),
		   // y = (n * y1 + m * y2) / (n + m);
	double l = m / n;
	double x = (x1 + l * x2) / (1 + l),
		   y = (y1 + l * y2) / (1 + l);
	// fitg << "\tshort:  " << outp(v.start, "start") << " -> " << outp(v.end, "end") << "\tm is " << m << ", n is " << n 
			// << outp(point(x, y), "\t\tavr_p") << endl;
	return point(x, y);
}

string str_vec_format(vec v)
{
	char str[100];
	sprintf(str, "TYPE = VECTOR\tCOLOR = 14\tWIDTH = 1\tSCALE = 1.00\tGEO = (%2f,%2f %2f,%2f)\n",
				 v.start.x, v.start.y, v.end.x, v.end.y);
	return string(str);
}*/

#endif // GEOMETRY_H