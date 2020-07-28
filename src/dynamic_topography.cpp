#include "dynamic_topography.h"

////////////////////////////////////////////////////////////////////////////////
// --------------------------- dt_result struct ------------------------------//
////////////////////////////////////////////////////////////////////////////////

dt_result::dt_result() : dt_error(-1.0), a_priori_error(-1.0)
{}

void dt_result::set(scut _cut, int vc) 
{
	cut = _cut;
	vector_count = vc;
	cut_length = _cut.v().length();
}

void dt_result::calc_dt(double latitude)
{
	dt_coef = coriolis_koef(latitude) / G;
	
	dt = (itg_res.sqr_value + itg_res.lin_value) / G;
}

void dt_result::print_to(std::ofstream &file)
{
	file << cut.start.x << " " << cut.start.y << " " << cut.end.x << " " << cut.end.y << " " << dt << " " 
		 << cut.width << " " << itg_res.itp_diameter << " " << itg_res.weight_coef * 1000 << " " << 
		dt_error << " " << itg_res.interpolation_accuracy << " " << itg_res.integration_error << " " 
		<< itg_res.ms_deviation << " " << a_priori_error << " " << cut_length << " " << cr_coef << " " 
		<< KM2M(itg_res.step_size) << " " << itg_res.step_count << " " << vector_count << std::endl;
}

std::string get_NV_filename(int ind)
{
	char str[10];
	sprintf(str, "NV%d.vec", ind);
	return std::string(str);
}

std::string get_AV_filename(int ind)
{
	char str[10];
	sprintf(str, "AV%d.vec", ind);
	return std::string(str);
}

////////////////////////////////////////////////////////////////////////////////
// ----------------------- DynamicTopography class ---------------------------//
////////////////////////////////////////////////////////////////////////////////

DynamicTopography::DynamicTopography(std::vector <movement> &m) : mvn(m)
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

int DynamicTopography::take(struct dt_result &dt_res)
{
	fNV.open(get_NV_filename(file_index).c_str());
	std::vector <wvector> wv;

	Line cut_line(cut.v());

	double to_start = cut.v().length(), to_end = cut.v().length();
	point start = cut.end, end = cut.start;

	double apr_err = 0.0;
	int apr_err_count = 0;

	std::ofstream fNVdec;
	fNVdec.open("NVdec.txt");

	std::ofstream fNVgeo;
	fNVgeo.open("NVgeo.txt");

	for (size_t j = 0; j < mvn.size(); ++j)
	{						
		double dist_to_vec = cut_line.distance_to(mvn[j].mv.start);

		if (fabs(dist_to_vec) < cut.width && mvn[j].velocity > 0)
		{
			// проекция начала вектора скорости на разрез
			point prj = cut_line.projection_of(mvn[j].mv.start);
			
			if (cut_line.contains(prj))
			{
				// прямая, параллельная разрезу и проходящая через конечную точку вектора скорости
				Line prl = cut_line.parallel(mvn[j].mv.end);
				// проекция начала вектора скорости на прямую prl (!) 
				point norm = prl.projection_of(mvn[j].mv.start);

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

				fNVdec << mvn[j].mv.start.x << " " << mvn[j].mv.start.y << " " << 
								norm.x << " " << norm.y << "\n";

				vec v = vec(mvn[j].mv.start, norm).at_geo_cs(dcs_origin);
				fNVgeo << v.start.x << " " << v.start.y << " " << 
								v.end.x << " " << v.end.y << "\n";
			}
		}
	}
	cut.start = start, cut.end = end;

	if (wv.empty())
	{
		std::cerr << "Desired flow velocity vectors near the cut are not found\n";
		return EC_DT_FVF_EMPTY;
	}

	Integral integral(cut, wv);
	integral.set_filename(get_AV_filename(file_index));
	integral.set_dcs_origin(dcs_origin);

	// расчёт интеграла
	integral.set_partitioning_count(wv.size() * 5);	
	int itg_code_error = integral.take(dt_res.itg_res);
	if (itg_code_error != EC_ITG_SUCCESS) 
		return itg_code_error;

	// расчет перепада ДТ по результатам интегрирования
	dt_res.set(cut, wv.size());
	dt_res.calc_dt(dcs_origin.y);	
	dt_res.a_priori_error = apr_err / apr_err_count;
	// расчет ошибки интегрирования
	integral.set_partitioning_count(wv.size() * 10);
	struct itg_result itg_res_2;
	itg_code_error = integral.take(itg_res_2);
	dt_res.dt_error = fabs(dt_res.itg_res.lin_value - itg_res_2.lin_value) * dt_res.dt_coef;

	fNVdec << " " << dt_res.cut.start.x << " " << dt_res.cut.start.y << " " << 
				dt_res.cut.end.x << " " << dt_res.cut.end.y << "\n";
	fNVdec.close();

	dt_res.cut.start.to_geo_cs(dcs_origin);
	dt_res.cut.end.to_geo_cs(dcs_origin);

	fNVgeo << " " << dt_res.cut.start.x << " " << dt_res.cut.start.y << " " << 
				dt_res.cut.end.x << " " << dt_res.cut.end.y << "\n";
	fNVgeo.close();

	fNV << dt_res.cut.v().toGlanceFormat(); // рисуем разрез вектором
	fNV.close();



	return EC_DT_SUCCESS;
}
