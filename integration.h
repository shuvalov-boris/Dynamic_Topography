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
	double lin_value; // I(Vt)dn
	double sqr_value; // I(Vt*Vt)dn
	double interpolation_accuracy; // точность интерполяции
	double integration_error; // ошибка расчета интеграла
	double ms_deviation; // mean-square deviation
	double step_size;
	double step_count;
	double itp_diameter;
	double weight_coef;

	itg_result() : lin_value(0.0), sqr_value(0.0), interpolation_accuracy(1.0), integration_error(1.0), 
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

	double len_K = wv[0].mvn.mv.length() / wv[0].mvn.velocity;

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
		prnd.shorten(val[i] * len_K);

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

	double lin_sum = 0.0;
	double sqr_sum = 0.0;
	for (size_t i = 0; i < val.size(); ++i)
	{
		lin_sum += val[i] * h;
		sqr_sum += val[i] * val[i] * h * sign(val[i]);
	}

	itg_res.lin_value = lin_sum;
	itg_res.sqr_value = sqr_sum;

	if (pm == EPM_ON)
	{
		fitp << cut.v().at_geo_cs(dcs_origin).toGlanceFormat();
		fitp.close();
	}

	return EC_ITG_SUCCESS;
}

#endif // INTEGRATION_H