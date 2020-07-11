#ifndef DYNAMIC_TOPOGRAPHY_H
#define DYNAMIC_TOPOGRAPHY_H

#include "geometry.h"
#include "integration.h"

#include <vector>
#include <fstream>
#include <iostream>

#define CUT_WIDTH 10 // [км]
#define G 9.81
// угловая скорость вращения Земли вокруг оси
#define EARTH_OMEGA 7.2921e-5

// коды ошибок при расчете перепада динамических высот
#define EC_DT_SUCCESS 1000
#define EC_DT_FVF_EMPTY 1002

double coriolis_koef(double fi)
{
	return 2 * EARTH_OMEGA * sin(fi * M_PI / 180);
}

double dyn_top_koef(double latitude)
{
	return coriolis_koef(latitude) / G;
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
	double cut_length;
	double dt_coef;
	double step_size;
	int step_count;
	int vector_count;
	double itp_diameter;
	double weight_coef;


	dt_result() {}

	dt_result(itg_result ir, double cut_width, int vc) : 
		cut(scut(ir.v(), cut_width)), dt(ir.value),
		interpolation_accuracy(ir.interpolation_accuracy),
		integration_error(ir.integration_error), dt_error(-1.0), a_priori_error(1.0), 
		ms_deviation(ir.ms_deviation), step_size(ir.step_size), step_count(ir.step_count), 
		vector_count(vc), itp_diameter(ir.itp_diameter), weight_coef(ir.weight_coef) 
	{

		cut_length = ir.v().length();

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

	void calc_dt(double latitude)
	{
		dt_coef = dyn_top_koef(latitude);
		dt *= dt_coef; // dt = integral.value * k
		// interpolation_accuracy *= k; 
	}

	// string toDTFormat()
	// {
	// 	char str[200];
	// 	sprintf(str, "%2f %2f %2f %2f %2f %2f %2f %2f %2e %2f %2f %2f %2f %2f %2f %2f %d %d", 
	// 			cut.start.x, cut.start.y, cut.end.x, cut.end.y, dt, 
	// 			cut.width, itp_diameter, weight_coef, 
	// 			dt_error, interpolation_accuracy, integration_error, ms_deviation, a_priori_error,
	// 			cut.v().length(), dt,
	// 			step_size, step_count, vector_count);
	// 	return string(str);
	// }

	void print_to(ofstream &file)
	{
		file << cut.start.x << " " << cut.start.y << " " << cut.end.x << " " << cut.end.y << " " << dt << " " << 
			cut.width << " " << itp_diameter << " " << weight_coef * 1000 << " " << 
			dt_error << " " << interpolation_accuracy << " " << integration_error << " " << ms_deviation << " " << a_priori_error << " " <<
			cut_length << " " << dt_coef << " " <<
			step_size * 1000 << " " << step_count << " " << vector_count << endl;
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
	point dcs_origin; // начало локальной Декартовой СК в глобальной ДСК
	vector <movement> mvn;

	ofstream fNV;
	int file_index;

public:
	DynamicTopography(vector <movement> &m);

	void set_file_index(int index);
	void set_cut(scut c);
	void set_dcs_origin(const point &dcs_orn);
	int take(struct dt_result &dtr);

};

DynamicTopography::DynamicTopography(vector <movement> &m) : mvn(m)
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

void DynamicTopography::set_dcs_origin(const point &dcs_orn)
{
	dcs_origin = dcs_orn;
}

int DynamicTopography::take(struct dt_result &dtr)
{
	fNV.open(get_NV_filename(file_index).c_str());
	vector <wvector> wv;

	Line cut_line(cut.v());

	double to_start = cut.v().length(), to_end = cut.v().length();
	point start = cut.end, end = cut.start;

	double apr_err = 0.0;
	int apr_err_count = 0;

	int along_line_count = 0;
	// ofstream fUV;
	// fUV.open("UV1.vec");
	for (size_t j = 0; j < mvn.size(); ++j)
	{
						
		double dist_to_vec = cut_line.distance_to(mvn[j].mv.start);

		if (fabs(dist_to_vec) < cut.width && mvn[j].velocity > 0)
		{
			along_line_count++;

			// fUV << mvn[j].mv.at_geo_cs(dcs_origin).toGlanceFormat();
			
			// проекция начала вектора скорости на разрез
			point prj = cut_line.projection_of(mvn[j].mv.start);

			
			if (cut_line.contains(prj))
			{
				// прямая, параллельная разрезу и проходящая через конечную точку вектора скорости
				Line prl = cut_line.parallel(mvn[j].mv.end);
				// проекция начала вектора скорости на прямую prl (!) 
				point norm = prl.projection_of(mvn[j].mv.start);

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

				fNV << vec(mvn[j].mv.start, norm).at_geo_cs(dcs_origin).toGlanceFormat();
			}
		}
	}
	// fUV.close();
	cut.start = start, cut.end = end;

	// cout << "along_line_count is " << along_line_count << endl;
	// cout << "wv size is " << wv.size() << endl;

	if (wv.empty())
	{
		cerr << "Desired flow velocity vectors near the cut are not found\n";
		return EC_DT_FVF_EMPTY;
	}

	// cout << start.toString("start") << "\t" << end.toString("end") << endl;

	Integral integral(cut, wv);
	integral.set_filename(get_AV_filename(file_index));
	integral.set_dcs_origin(dcs_origin);

	integral.set_partitioning_count(wv.size() * 5);
	dtr = dt_result(integral.take(), cut.width, wv.size());	

	integral.set_partitioning_count(wv.size() * 10);
	struct dt_result dtr2(integral.take(EPM_ON), cut.width, wv.size());
	
	if (itg_err == EIE_SUCCESS) dtr.a_priori_error = apr_err / apr_err_count;
	dtr.dt_error = fabs(dtr.dt - dtr2.dt) * dtr.dt_coef;

	dtr.calc_dt(dcs_origin.y);

	dtr.cut.start.to_geo_cs(dcs_origin);
	dtr.cut.end.to_geo_cs(dcs_origin);

	fNV << dtr.cut.v().toGlanceFormat(); // зачем? чтоб выделить разрез?
	fNV.close();

	return EC_DT_SUCCESS;
}

#endif //DYNAMIC_TOPOGRAPHY_H