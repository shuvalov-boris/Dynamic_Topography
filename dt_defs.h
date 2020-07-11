#ifndef DT_DEFS_H
#define DT_DEFS_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "geometry.h"

// #define EQUATOR_LENGTH 40075686. // [meters]
#define MERIDIAN_LENGTH 40008600 // [meters]
#define EQUATOR_RADIUS  6378245. // [meters]
#define POLAR_RADIUS    6356863. // [meters]
#define EARTH_FLATTENING (EQUATOR_RADIUS - POLAR_RADIUS) / EQUATOR_RADIUS // [meters] - полярное сжатие Земли

#define METERS_IN_ONE_DEG MERIDIAN_LENGTH / 360. 	// метров в одном градусе

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

	scut() {}

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

// origin - центр локальной декартовой системы координат в глобальных декартовых координатах
point geo2dec(const point &gp, const point &origin = point(0, 0))
{
	point ggp(gp.x, gp.y);
	ggp.to_dec_cs(origin);
	return ggp;
}

point dec2geo(const point &dp, const point &origin = point(0, 0))
{
	point gdp(dp.x, dp.y);
	gdp.to_geo_cs(origin);
	return gdp;
}

// origin - в географических координатах
void point::to_dec_cs(const point &origin) 
{
	this->to_related_cs(origin);

	// радиус Земли (от её центра к точке на поверхности заданной географической широты)
	double earth_radius_at_latitude = EQUATOR_RADIUS * (1. - EARTH_FLATTENING * sqr(sin(this->y * M_PI / 180)));

	// радиус сечения Земли на заданной географической широте
	double little_radius_at_latitude = earth_radius_at_latitude * cos(this->y * M_PI / 180);

	this->x *= little_radius_at_latitude * M_PI / 180;
	this->y *= METERS_IN_ONE_DEG;

	this->x /= 1000;
	this->y /= 1000;
}

void point::to_geo_cs(const point &origin)
{
	this->y *= 1000;
	this->x *= 1000;

	this->y = this->y / MERIDIAN_LENGTH * 360.;	
	// cout << fixed << setprecision(10) << "to_geo_cs: gdp.y = " << this->y << endl;

	// радиус Земли (от её центра к точке на поверхности заданной географической широты)
	double earth_radius_at_latitude = EQUATOR_RADIUS * 
								(1. - EARTH_FLATTENING * sqr(sin(this->y * M_PI / 180)));
	// cout << "to_geo_cs: earth_radius_at_latitude is " << earth_radius_at_latitude << " m\n";

	// радиус сечения Земли на заданной географической широте
	double little_radius_at_latitude = earth_radius_at_latitude * cos(this->y * M_PI / 180);
	// cout << "to_geo_cs: little_radius_at_latitude is " << little_radius_at_latitude << " m\n";

	this->x *= 180. / M_PI / little_radius_at_latitude;
	// cout << "to_geo_cs: longitude = " << this->x << endl;			
	// cout << "to_geo_cs: latitude = " << this->y << endl;	

	this->to_global_cs(origin);
	// cout << fixed << setprecision(10) << "to_geo_cs: " << this->toString("pt at global DCS") << endl;
}


void to_cartesian_cs(vector <movement> &mvn, vector <scut> &station)
{
	point dcs_dec_origin = station[0].v().middle(); // not dec but geo -> rename variable

	ofstream f_dec_mvn_1;
	f_dec_mvn_1.open("dec_mvn_1.txt");
	for (size_t j = 0; j < mvn.size(); ++j)
	{
		mvn[j].mv.start.to_dec_cs(dcs_dec_origin);
		mvn[j].mv.end.to_dec_cs(dcs_dec_origin);
		f_dec_mvn_1 << mvn[j].mv.start.toString() << "   " << mvn[j].mv.end.toString() << endl;
	}
	f_dec_mvn_1.close();

	for (size_t j = 0; j < station.size(); ++j)
	{
		station[j].start.to_dec_cs(dcs_dec_origin);
		station[j].end.to_dec_cs(dcs_dec_origin);
	}
}


#endif // DT_DEFS_H