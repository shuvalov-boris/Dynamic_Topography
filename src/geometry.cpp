#include "geometry.h"

////////////////////////////////////////////////////////////////////////////////
// ----------------------------- point struct --------------------------------//
////////////////////////////////////////////////////////////////////////////////

bool point::equal(const point &p)
{
	return (x == p.x && y == p.y);
}

bool point::equal_eps(const point &p)
{
	return (fabs(x - p.x) < EPS && fabs(y - p.y) < EPS);
}

double point::distance_to(point p) const
{
	return sqrt((x - p.x) * (x - p.x) + (y - p.y) * (y - p.y));
}

void point::to_related_cs(const point &origin)
{
	x -= origin.x, y-= origin.y;
}

void point::to_global_cs(const point &origin)
{
	x += origin.x, y+= origin.y;
}

point point::at_geo_cs(const point &origin) const
{
	point pt(this->x, this->y);
	pt.to_geo_cs(origin);
	return pt;
}

point point::at_dec_cs(const point &origin) const
{
	point pt(this->x, this->y);
	pt.to_dec_cs(origin);
	return pt;
}

std::string point::toString(std::string name) const
{
	char str[100];
	sprintf(str, "%s(%2f, %2f)", name.c_str(), x, y);
	return std::string(str);
}

////////////////////////////////////////////////////////////////////////////////
// ------------------------------ vec struct ---------------------------------//
////////////////////////////////////////////////////////////////////////////////

bool vec::contains(point p)
{
	return sign(p.x - start.x) * sign(p.x - end.x) < 0;
}

double vec::length()
{
	return sqrt((start.x - end.x) * (start.x - end.x) + (start.y - end.y) * (start.y - end.y));
}

void vec::shorten(double len/*, double dl = 0.0*/)
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

point vec::middle()
{
	return point((start.x + end.x) / 2, (start.y + end.y) / 2);
}

std::string vec::toString(std::string name)
{
	char str[100];
	sprintf(str, "%s[%s -> %s]", name.c_str(), start.toString().c_str(), end.toString().c_str());
	return std::string(str);
}

vec vec::at_geo_cs(const point &origin)
{
	return vec(this->start.at_geo_cs(origin), this->end.at_geo_cs(origin));
}

std::string vec::toGlanceFormat()
{
	char str[200];
	sprintf(str, "TYPE = VECTOR\tCOLOR = 14\tWIDTH = 1\tSCALE = 1.00\tGEO = (%2f,%2f %2f,%2f)\n",
				 start.x, start.y, end.x, end.y);
	return std::string(str);
}

////////////////////////////////////////////////////////////////////////////////
// ------------------------------ Line class ---------------------------------//
////////////////////////////////////////////////////////////////////////////////

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
	// if (f(p) > EPS)
		// cout << p.toString("for pt") << " f(pt) > EPS;\t f(pt) = " << f(p) << endl;
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

point Line::projection_of(point p)
{
	Line l = this->perpendicular(p);
	// cout << l.as_string("vel vec prnd cut line") << endl;
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

std::string Line::as_string(std::string name) const
{
	char str[200];
	sprintf(str, "%s:\t%s; %s; (a=%2f, b=%2f, c=%2f)", name.c_str(), p1.toString("p1").c_str(), p2.toString("p2").c_str(), a, b, c);
	return std::string(str);
}