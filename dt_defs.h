#ifndef DT_DEFS_H
#define DT_DEFS_H

#include "geometry.h"

// TODO расчет относительно широты координат вектора будет более точным
#define EQUATOR_LENGTH 40075686. // [meters]
#define EQUATOR_RADIUS  6378245. // [meters]
#define POLAR_RADIUS    6356863. // [meters]
#define EARTH_FLATTENING (EQUATOR_RADIUS - POLAR_RADIUS) / EQUATOR_RADIUS // [meters] - полярное сжатие Земли

#define GR EQUATOR_LENGTH / 360. 	// метров в одном градусе
#define METERS_IN_ONE_DEG GR

#define D2M(a) a * GR

#define sqr(a) ((a)*(a))

struct movement
{
	vec mv;
	double velocity;
	double error;		// априорная ошибка
	movement(vec v, double vl, double err) : mv(v), velocity(vl), error(err) {} 
};

struct wvector
{
	movement mvn;
	point norm_comp;	//normal component for cut
	point proj;			//start point projection on cut

	wvector(movement _mvn, point nc, point pr) : mvn(_mvn), norm_comp(nc), proj(pr) {}

	vec v() { return mvn.mv; }
};

struct scut // Разрез
{
	point start;
	point end;
	double width;
	double itp_diameter;
	double weight_coef;

	// cut(vec v, double w) : 
		// start(v.start), end(v.end), width(w) {}

	scut(vec v, double w, double id = 0.0, double wc = 0.0) : 
		start(v.start), end(v.end), width(w),
		itp_diameter(id), weight_coef(wc) {}

	vec v()	{ return vec(start, end); }

	string toString(string name = string(""))
	{
		char str[30];
		sprintf(str, "%s(%2f, %2f)", name.c_str(), width, weight_coef);
		return string(str);
	}
};


#endif // DT_DEFS_H