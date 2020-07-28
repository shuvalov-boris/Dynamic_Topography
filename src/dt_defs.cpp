#include "dt_defs.h"


////////////////////////////////////////////////////////////////////////////////
// ----------------------------- scut struct ---------------------------------//
////////////////////////////////////////////////////////////////////////////////

scut::scut(vec v, double w, double id, double wc) : 
	scut(v, w, id, wc, point(0., 0.), false) {}

scut::scut(vec v, double w, double id, double wc, point cc, bool correction) :
	start(v.start), end(v.end), width(w), itp_diameter(id), weight_coef(wc / 1000),
	 curvature_center(cc), curvature_correction(correction) {}

vec scut::v(){
	return vec(start, end);
}

std::string scut::toString(std::string name)
{
	char str[200];
	sprintf(str, "%s(%2f, %2f) with %s", name.c_str(), width, weight_coef * 1000, 
		curvature_center.toString("CC").c_str());
	return std::string(str);
}


////////////////////////////////////////////////////////////////////////////////
// ------------------------------ functions ----------------------------------//
////////////////////////////////////////////////////////////////////////////////


double coriolis_koef(double fi)
{
	return 2 * EARTH_OMEGA * sin(fi * M_PI / 180);
}

// origin - центр локальной декартовой системы координат в глобальных декартовых координатах
point geo2dec(const point &gp, const point &origin)
{
	point ggp(gp.x, gp.y);
	ggp.to_dec_cs(origin);
	return ggp;
}

point dec2geo(const point &dp, const point &origin)
{
	point gdp(dp.x, dp.y);
	gdp.to_geo_cs(origin);
	return gdp;
}

// origin - центр декартовой СК в географических координатах
void point::to_dec_cs(const point &origin) 
{
	float latitude = this->y;

	this->to_related_cs(origin);

	// радиус Земли (от её центра к точке на поверхности заданной географической широты)
	double earth_radius_at_latitude = EQUATOR_RADIUS * (1. - EARTH_FLATTENING * sqr(sin(latitude * M_PI / 180)));

	// радиус сечения Земли на заданной географической широте
	double little_radius_at_latitude = earth_radius_at_latitude * cos(latitude * M_PI / 180);

	this->x *= little_radius_at_latitude * M_PI / 180;
	this->y *= METERS_IN_ONE_DEG;

	this->x = M2KM(this->x);
	this->y = M2KM(this->y);
}

void point::to_geo_cs(const point &origin)
{
	this->x = KM2M(this->x);
	this->y = KM2M(this->y);

	this->y = this->y / MERIDIAN_LENGTH * 360.;	

	float latitude = this->y + origin.y;

	// радиус Земли (от её центра к точке на поверхности заданной географической широты)
	double earth_radius_at_latitude = EQUATOR_RADIUS * 
								(1. - EARTH_FLATTENING * sqr(sin(latitude * M_PI / 180)));

	// радиус сечения Земли на заданной географической широте
	double little_radius_at_latitude = earth_radius_at_latitude * cos(latitude * M_PI / 180);

	this->x *= 180. / M_PI / little_radius_at_latitude;

	this->to_global_cs(origin);
}


void to_cartesian_cs(std::vector <movement> &mvn, std::vector <scut> &station)
{
	point dcs_geo_origin = station[0].v().middle();

	for (size_t j = 0; j < mvn.size(); ++j)
	{
		mvn[j].mv.start.to_dec_cs(dcs_geo_origin);
		mvn[j].mv.end.to_dec_cs(dcs_geo_origin);
	}

	for (size_t j = 0; j < station.size(); ++j)
	{
		station[j].start.to_dec_cs(dcs_geo_origin);
		station[j].end.to_dec_cs(dcs_geo_origin);
		station[j].curvature_center.to_dec_cs(dcs_geo_origin);
	}
}