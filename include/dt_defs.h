#ifndef DT_DEFS_H
#define DT_DEFS_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "geometry.h"

// #define EQUATOR_LENGTH 40075686. // [meters]
#define MERIDIAN_LENGTH 40008600. // [meters]
#define EQUATOR_RADIUS  6378245. // [meters]
#define POLAR_RADIUS    6356863. // [meters]
#define EARTH_FLATTENING (EQUATOR_RADIUS - POLAR_RADIUS) / EQUATOR_RADIUS // [meters] - полярное сжатие Земли

#define METERS_IN_ONE_DEG MERIDIAN_LENGTH / 360. 	// метров в одном градусе

// угловая скорость вращения Земли
#define EARTH_OMEGA 7.2921e-5

#define G 9.81

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
	point curvature_center; // центр кривизны потока
	bool curvature_correction = false; // учитывать ли кривизнц потока при расчете ДТ

	// cut(vec v, double w) : 
		// start(v.start), end(v.end), width(w) {}

	scut() {}

	scut(vec v, double w, double id, double wc);

	scut(vec v, double w, double id, double wc, point cc, bool correction = true);

	vec v();

	std::string toString(std::string name = std::string(""));
};

double coriolis_koef(double fi);

// origin - центр локальной декартовой системы координат в глобальных декартовых координатах
point geo2dec(const point &gp, const point &origin = point(0, 0));

point dec2geo(const point &dp, const point &origin = point(0, 0));

void to_cartesian_cs(std::vector <movement> &mvn, std::vector <scut> &station);

#endif // DT_DEFS_H