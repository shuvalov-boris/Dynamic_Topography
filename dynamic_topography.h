#ifndef DYNAMIC_TOPOGRAPHY_H
#define DYNAMIC_TOPOGRAPHY_H

// #include "geometry.h"
#include "integration.h"

#include <iomanip>
#include <iostream>
#include <vector>

#define CUT_WIDTH 0.1
#define G 9.81
#define O 7.29e-5

double coriolis(double fi)
{
	return 2 * O * sin(fi * M_PI / 180);
}

double dyn_top_koef(vec v)
{
	double angle_avr = (v.start.y + v.end.y) / 2;
	double crl = coriolis(angle_avr);
	return crl / G;
}

// origin - центр локальной декартовой системы координат в глобальных декартовых координатах
point geo2dec(const point &gp, const point &origin = point(0, 0))
{
	// радиус Земли (от её центра к точке на поверхности заданной географической широты)
	double earth_radius_at_latitude = EQUATOR_RADIUS * (1. - EARTH_FLATTENING * sqr(sin(gp.y * M_PI / 180)));
	// cout << fixed << setprecision(10) << "geo2dec: earth_radius_at_latitude is " << earth_radius_at_latitude << " m\n";

	// радиус сечения Земли на заданной географической широте
	double little_radius_at_latitude = earth_radius_at_latitude * cos(gp.y * M_PI / 180);
	// cout << "geo2dec: little_radius_at_latitude is " << little_radius_at_latitude << " m\n";

	// cout << "geo2dec: x = " << little_radius_at_latitude * gp.x * M_PI / 180 << endl;
	// cout << "geo2dec: y = " << METERS_IN_ONE_DEG * gp.y << endl;

	return point(
		little_radius_at_latitude * gp.x * M_PI / 180,
		METERS_IN_ONE_DEG * gp.y
	).to_origin_related(origin);
}

point dec2geo(const point &dp, const point &origin = point(0, 0))
{
	const point gdp = dp.to_global(origin);

	double latitude = gdp.y / EQUATOR_LENGTH * 360.;	
	// cout << fixed << setprecision(10) << "dec2geo: dp.y = " << dp.y << " ;\t" << "METERS_IN_ONE_DEG = " << METERS_IN_ONE_DEG << endl;
	// cout << "dec2geo: dp.y / EQUATOR_LENGTH = " << dp.y / EQUATOR_LENGTH << endl;

	// радиус Земли (от её центра к точке на поверхности заданной географической широты)
	double earth_radius_at_latitude = EQUATOR_RADIUS * 
								(1. - EARTH_FLATTENING * sqr(sin(latitude * M_PI / 180)));
	// cout << "dec2geo: earth_radius_at_latitude is " << earth_radius_at_latitude << " m\n";

	// радиус сечения Земли на заданной географической широте
	double little_radius_at_latitude = earth_radius_at_latitude * cos(latitude * M_PI / 180);
	// cout << "dec2geo: little_radius_at_latitude is " << little_radius_at_latitude << " m\n";

	double 	longitude = gdp.x * 180. / M_PI / little_radius_at_latitude;	
	// cout << "dec2geo: longitude = " << longitude << endl;			
	// cout << "dec2geo: latitude = " << latitude << endl;	

	return point(longitude, latitude);
}


void to_cartesian_cs(vector <movement> &mvn, point dcs_geo_origin)
{
	point dcs_dec_origin = geo2dec(dcs_geo_origin);

	for (size_t j = 0; j < mvn.size(); ++j)
	{
		mvn[j].mv.start = geo2dec(mvn[j].mv.start, dcs_dec_origin);
		mvn[j].mv.end = geo2dec(mvn[j].mv.end, dcs_dec_origin);
	}
}


struct dt_result
{
	scut cut; 					// разрез
	double dt; 					// динамическая топография
	double interpolation_accuracy; //точность интерполяции
	double integration_error;	// ошибка расчета
	double dt_error;			// разность двух значени ДТ с разным количеством интервалов
	double a_priori_error; 		// апроиорная ошибка
	double ms_deviation;		// среднеквадратичное отклонение
	double step_size;
	int step_count;
	int vector_count;
	double itp_diameter;
	double weight_coef;	

	dt_result(itg_result ir, double cut_width, int vc = 0) : 
		cut(scut(ir.v(), cut_width)), dt(ir.value),
		interpolation_accuracy(ir.interpolation_accuracy),
		integration_error(ir.integration_error), ms_deviation(ir.ms_deviation),
		dt_error(-1.0), a_priori_error(1.0), step_count(ir.step_count), vector_count(vc), 
		step_size(ir.step_size), itp_diameter(ir.itp_diameter), weight_coef(ir.weight_coef) 
	{


		// switch (itg_err)
		// {
		// 	case EIE_SUCCESS:
		// 	{
		// 		error = itg_res.error;
		// 		break;
		// 	}
		// 	case EIE_NOT_ENOUGH_DATA:
		// 	{
		// 		error = 1000.0;
		// 		break;
		// 	}
		// 	case EIE_OTHER:
		// 	{
		// 		error = 1.0;
		// 		break;
		// 	}
		// 	default: 
		// 	{
		// 		break;
		// 	}
					
		// }
	}

	void calc_dt()
	{
		double k = dyn_top_koef(cut.v());
		dt *= k;
		// interpolation_accuracy *= k; 
	}

	string toDTFormat()
	{
		char str[200];
		sprintf(str, "%2f %2f %2f %2f %2f %2f %2f %2f %2f %2f %2f %2f %2f %2f %2f %2f %d %d", 
				cut.start.x, cut.start.y, cut.end.x, cut.end.y, dt, 
				cut.width, itp_diameter, weight_coef, 
				dt_error, interpolation_accuracy, integration_error, ms_deviation, a_priori_error,
				D2M(cut.v().length()), dyn_top_koef(cut.v()),
				step_size, step_count, vector_count);
		return string(str);
	}
};

string get_NV_filename(int ind)
{
	char str[10];
	sprintf(str, "NV%d.vec", ind);
	return string(str);
}

string get_AV_filename(int ind)
{
	char str[10];
	sprintf(str, "AV%d.vec", ind);
	return string(str);
}

class DynamicTopography
{
	scut cut;
	vector <movement> mvn;

	ofstream fNV;
	int file_index;

public:
	DynamicTopography(vector <movement> &m);

	void set_file_index(int index);
	void set_cut(scut c);

	dt_result take();

};

DynamicTopography::DynamicTopography(vector <movement> &m) : mvn(m), cut(vec(point(0, 0), point(0, 0)), 0.0)
{

}

void DynamicTopography::set_file_index(int index)
{
	file_index = index;
}

void DynamicTopography::set_cut(scut c)
{
	cut = c;
	if (cut.width == -1) cut.width = CUT_WIDTH;
}

dt_result DynamicTopography::take()
{

	fNV.open(get_NV_filename(file_index).c_str());

	vector <wvector> wv;

	Line interval(cut.v());

	double to_start = cut.v().length(), to_end = cut.v().length();
	point start = cut.end, end = cut.start;

	double apr_err = 0.0;
	int apr_err_count = 0;

	for (size_t j = 0; j < mvn.size(); ++j)
	{
						
		double f = interval.distance_to(mvn[j].mv.start);
		if (abs(f) < cut.width && mvn[j].velocity > 0)
		{
			point prj = interval.projection(mvn[j].mv.start);
			if (interval.contains(prj))
			{
				Line prl = interval.parallel(mvn[j].mv.end);
				point norm = prl.projection(mvn[j].mv.start);
	
				// flog << "\t\t" << mvn[j].mv.toString("mv") << "\t" << ppt.toString("NV") << endl;

				wv.push_back(wvector(mvn[j], norm, prj));

				if (prj.distance_to(cut.start) < to_start)
				{
					to_start = prj.distance_to(cut.start);
					start = prj;
				}
				if (prj.distance_to(cut.end) < to_end)
				{
					to_end = prj.distance_to(cut.end);
					end = prj;
				}

				apr_err += mvn[j].error;
				++apr_err_count;

				fNV << vec(mvn[j].mv.start, norm).toGlanceFormat();
			}
		}
	}

	cut.start = start, cut.end = end;

	// cout << start.toString("start") << "\t" << end.toString("end") << endl;

	Integral integral(cut, wv);
	integral.set_filename(get_AV_filename(file_index));

	integral.set_partitioning_count(wv.size() * 5);
	struct dt_result dtr(integral.take(), cut.width, wv.size());	

	integral.set_partitioning_count(wv.size() *10);
	struct dt_result dtr2(integral.take(EPM_OFF), cut.width);
	
	if (itg_err == EIE_SUCCESS) dtr.a_priori_error = apr_err / apr_err_count;
	dtr.dt_error = fabs(dtr.dt - dtr2.dt) * dyn_top_koef(dtr.cut.v());

	dtr.calc_dt();

	fNV << dtr.cut.v().toGlanceFormat();
	fNV.close();

	return dtr;
}

#endif //DYNAMIC_TOPOGRAPHY_H