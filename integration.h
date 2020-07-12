#ifndef INTEGRATION_H
#define INTEGRATION_H

#include "interpolation.h"
#include "dt_defs.h"
#include "log.h"

#include <string>
#include <cmath>
#include <vector>
#include <algorithm>

#define MIN_POINT_COUNT 10

// коды ошибок при расчете интеграла
#define EC_ITG_SUCCESS 1000
#define EC_ITG_NOT_ENOUGH_DATA 2004

enum E_INTERPOLATION_MODE
{
	EIM_LINEAR,
	ETM_WEIGTH_FUNC
};

enum E_PRINT_MODE
{
	EPM_ON,
	EPM_OFF
};

struct itg_result
{
	double value;
	double interpolation_accuracy; // точность интерполяции
	double integration_error; // ошибка расчета интеграла
	double ms_deviation; // mean-square deviation
	double step_size;
	double step_count;
	double itp_diameter;
	double weight_coef;

	itg_result() : value(0.0), interpolation_accuracy(1.0), integration_error(1.0), 
		ms_deviation(1.0), step_size(0.0), step_count(0), itp_diameter(0.0), weight_coef(0.0)
	{}
};

ofstream fAV;
// E_INTEGRATION_ERROR itg_err;

class Integral
{
	scut cut;
	point dcs_origin;
	vector <wvector> wv;

	int n; // количество интервалов разбиения

	E_PRINT_MODE print_mode;
	ofstream fitp;

	double get_integration_error(vector <double> &val, double h, int n);

public:
	Integral(scut c, vector <wvector> &_wv);
	~Integral();

	void set_cut(scut c);
	void set_dcs_origin(const point &dcs_orn);
	void set_partitioning_count(int _n);
	void set_filename(string filename);

	int take(struct itg_result &itg_res, E_PRINT_MODE pm);
};

/*
double Integral::get_integration_error(vector <double> &val, double h, int n)
{
	vector <double> dv1(val.size() - 1);
	vector <double> dv2(dv1.size() - 1);
	double err = 0.0;
	for (int i = 0; i < val.size() - 1; ++i)
	{
		dv1[i] = fabs(val[i + 1] - val[i]);
		if (i > 0) 
		{
			dv2[i - 1] = fabs(dv1[i] - dv1[i - 1]);
			err = max(err, dv2[i - 1]);
		}
	}
	return err * h * h * h * n / 24;
}*/

Integral::Integral(scut c, vector <wvector> &_wv) : cut(c), wv(_wv)
{

}

Integral::~Integral()
{
	// fitp.close();
}

void Integral::set_cut(scut c)
{
	cut = c;
}

void Integral::set_dcs_origin(const point &dcs_orn)
{
	dcs_origin = dcs_orn;
}

void Integral::set_partitioning_count(int _n)
{
	n = _n;
}

void Integral::set_filename(string filename)
{
	fitp.open(filename.c_str());
}

int Integral::take(struct itg_result &itg_res, E_PRINT_MODE pm = EPM_ON)
{
	if (wv.size() < MIN_POINT_COUNT)
	{
		cerr << "Error: not enough data to calculate the integral\n";
		return EC_ITG_NOT_ENOUGH_DATA;
	}

	// refresh_cut();

	double K = wv[0].mvn.mv.length() / wv[0].mvn.velocity;

	double h = cut.v().length() * 1000 / n; // шаг в метрах
	
	double dx = (cut.end.x - cut.start.x) / n;
	double dy = (cut.end.y - cut.start.y) / n;

	Line cut_line(cut.v());

	// fitg << "dist(interval) = " << interval.length() << "\tdist_metr(interval) = " << dist_metr(interval) << endl;

	vector <double> val(n, 0.0);
	vector <double> tan_avr;
	vector <int> pc(n);

	Interpolation itp(cut.v(), wv);

	if (cut.itp_diameter == -1)
		itp.calc_radius();
	if (cut.itp_diameter >= 0.0) 
		itp.set_radius(cut.itp_diameter / 2);
	if (cut.weight_coef >= 0.0) itp.set_weight_coef(cut.weight_coef);

	// double apr_err = 0.0;
	// int apr_err_count = 0;

	int step_count = 0;

	// fitg << cut.v().toString("station") << endl;
	// fitg << "Point count (for n = " << n << "\th = " << h << "\th_m = " << h_m << "\tdx = " << dx << "\tdy = " << dy << ")\n\n";
	for (int i = 0; i < n; ++i)
	{
		point gr1(cut.start.x + i * dx, cut.start.y + i * dy);
		point gr2(cut.start.x + (i + 1) * dx, cut.start.y + (i + 1) * dy);
		point gr_avr = vec(gr1, gr2).middle();

		// fitg << gr1.toString("gr1") << "\t" << gr2.toString("gr2") << "\t\ti = " << i << endl;

		val[i] = itp.take_for(gr_avr);

		++step_count;

		vec prnd = cut_line.perpendicular(gr_avr);
		prnd.shorten(val[i] * K);

		if (pm == EPM_ON)
			fitp << prnd.at_geo_cs(dcs_origin).toGlanceFormat();

		// fitg << "\t\t" << prnd.toString("avr vec") << endl;

	}
	// fitg << endl;

	vector <double> itp_acr;	// точность интерполяции

	itp.calc_accuracy(itp_acr);

	int count = itp_acr.size();

	double acr_sum = 0.0;
	double acr2_sum = 0.0;
	for (int i = 0; i < count; ++i)
	{
		acr_sum += itp_acr[i];
		acr2_sum += itp_acr[i] * itp_acr[i];

	}
	// cout << "sum of sq is" << acr2_sum << "\t" << count <<  endl;
	itg_res.interpolation_accuracy = acr_sum / count;
	itg_res.integration_error = sqrt(acr2_sum / count) /** interval.length() * GR*/;

	double msd_sum = 0.0;
	for (int i = 0; i < count; ++i)
		msd_sum += (itg_res.interpolation_accuracy - itp_acr[i]) * (itg_res.interpolation_accuracy - itp_acr[i]);
	itg_res.ms_deviation = sqrt(msd_sum / count );

	itg_res.step_size = h;
	itg_res.step_count = step_count;
	itg_res.itp_diameter = itp.get_radius() * 2;
	itg_res.weight_coef = itp.get_weight_coef();	

	// print_array(pc, "point count in every interval is");

	// print_array(val, "calculated avrage value in metres is");

	double sum = 0.0;
	for (size_t i = 0; i < val.size(); ++i)
		sum += val[i] * h;

	itg_res.value = sum;

	if (pm == EPM_ON)
	{
		fitp << cut.v().at_geo_cs(dcs_origin).toGlanceFormat();
		fitp.close();
	}

	return EC_ITG_SUCCESS;
}

/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/

// void print_array(vector <double> ar, string str)
// {
// 	flog << str << endl;
// 	for (size_t i = 0; i < ar.size(); ++i)
// 		flog << i << ":" << ar[i] << "\t";
// 	flog << endl;
// }

// double dist_metr(vec v)
// {
// 	// double GR = 4e7 / 360;
// 	double dx = fabs(v.end.x - v.start.x) * GR/* * cos((v.start.y + v.end.y) / 2)*/, 
// 		   dy = fabs(v.end.y - v.end.x) * GR;
// 	// return sqrt(dx * dx + dy * dy);
// 	return v.length() * GR;
// }


// void interpolate(vector <double> &val) // line (линейное)
// {
// 	int i = 0;
// 	while (i < val.size() && val[i] == 0) ++i;
// 	if (i > 0)
// 	{
// 		for (int j = 0; j < i; ++j)
// 			val[j] = val[i];
// 	}
// 	while (i < val.size())
// 	{
// 		while (i < val.size() && val[i] != 0) ++i;
// 		if (i >= val.size()) break;
// 		int j = i;
// 		while (j < val.size() && val[j] == 0) ++j;
// 		if (j == val.size())
// 		{
// 			for (int l = i; l < j; ++l)
// 				val[l] = val[i - 1];
// 			break;
// 		} else {
// 			int nn = j - i;
// 			for (int l = 0; l < nn; ++l)
// 				val[i + l] = val[i - 1] + (val[j] - val[i - 1]) / (nn + 1) * (l + 1);
// 		}
// 		i = j;
// 	}
// }

/*double integral(vector <movement> list, vector <point> proj, vector <point> tan, vec st, string filename) //list - point index from move point vector; 
																										//points - point on station line
{

	fAV.open(filename.c_str());

	double K = list[0].mv.length() / list[0].velocity;
	// cout << "size of vec is " << vec.size() << endl;

	int n = list.size() * 2;
	double h = st.length() / n, h_m = dist_metr(st) / n;
	double dx = (st.end.x - st.start.x) / n;
	double dy = (st.end.y - st.start.y) / n;

	fitg << "dist(st) = " << st.length() << "\tdist_metr(st) = " << dist_metr(st) << endl;

	vector <double> val(n);
	vector <double> lval(n);
	vector <double> tan_avr;
	vector <double> pc(n, 0);

	Line interval(st);

	fitg << st.toString("station") << endl;
	fitg << "Point count (for n = " << n << "\th = " << h << "\th_m = " << h_m << "\tdx = " << dx << "\tdy = " << dy << ")\n\n"; //in each interval is\n\t";
	for (int i = 0; i < n; ++i)
	{
		int point_count = 0;
		point gr1(st.start.x + i * dx, st.start.y + i * dy);
		point gr2(st.start.x + (i + 1) * dx, st.start.y + (i + 1) * dy);
		point gr_avr = vec(gr1, gr2).middle();;

		double len_sum = 0.0;
		double sp_sum = 0.0;
		fitg << gr1.toString("gr1") << "\t" << gr2.toString("gr2") << "\t\ti = " << i << endl;

		for (int j = 0; j < proj.size(); ++j)
		{
			if (vec(gr1, gr2).contains(proj[j]))
			{
				++point_count;
				double sgn = sign_vv(st, list[j].mv);
				double d = list[j].mv.start.distance_to(tan[j]);
				// double d = dist_metr(vec(list[j].start, tan[j]));
				fitg << "\t" << j << ":\t" << list[j].mv.toString("move vector") << ", " << vec(list[j].mv.start, tan[j]).toString("tangential component") << 
						"\n\t\t\t\t\t\t\t\t\t\ttc:\tlength is " << list[j].mv.start.distance_to(tan[j]) << "\tsign is " << sgn << endl;;
				len_sum += d * sgn;
				sp_sum += list[j].velocity * sgn;
			}
		}

		// cout << point_count << " ";
		pc[i] = point_count;
		if (point_count != 0)
		{
			tan_avr.push_back(len_sum / point_count);
			val[i] = sp_sum / point_count;
			lval[i] = len_sum / point_count;
			
		}
	}
	fitg << endl;

	print_array(val, "init val is");


	interpolate(val);

	for (int i = 0; i < n; ++i)
	{
		point gr1(st.start.x + i * dx, st.start.y + i * dy);
		point gr2(st.start.x + (i + 1) * dx, st.start.y + (i + 1) * dy);
		point gr_avr = vec(gr1, gr2).middle();;

		// point avr_long_end = interval.perpendicular(gr_avr);
		// point avr_end = short_v(vec(gr_avr, avr_long_end), val[i] * K);

		vec prnd = interval.perpendicular(gr_avr);
		prnd.shorten(val[i] * K);
		// favrtan << str_vec_format(vec(gr_avr, avr_long_end));
		fAV << prnd.toGlanceFormat();
		fitg << "\t\t" << prnd.toString("avr vec") << endl;
		// fitg << "\tpoint_count is " << pc[i] << ", vector length sum is " << len_sum << ", vls*dx = " << lval[i] << "vss*dx = " << val[i] << endl; 	
	}

	print_array(pc, "point count in every interval is");

	print_array(tan_avr, "tangential component length is");

	print_array(val, "calculated avrage [speed value * dx] in metres is");

	print_array(lval, "calculated avrage [tan com length * dx] in gr is");

	fAV.close();

	double sum = 0.0;
	for (int i = 0; i < val.size(); ++i)
		sum += val[i] * h_m;

	return sum;
}*/

// double get_integration_error(vector <double> val, double h, int n)
// {
// 	vector <double> dv1(val.size() - 1);
// 	vector <double> dv2(dv1.size() - 1);
// 	double err = 0.0;
// 	for (int i = 0; i < val.size() - 1; ++i)
// 	{
// 		dv1[i] = fabs(val[i + 1] - val[i]);
// 		if (i > 0) 
// 		{
// 			dv2[i - 1] = fabs(dv1[i] - dv1[i - 1]);
// 			err = max(err, dv2[i - 1]);
// 		}
// 	}
// 	return err * h * h * h * n / 24;
// }

#endif // INTEGRATION_H