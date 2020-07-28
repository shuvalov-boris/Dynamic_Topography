#include "integration.h"


Integral::Integral(scut c, std::vector <wvector> &_wv) : cut(c), wv(_wv)
{

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

void Integral::set_filename(std::string filename)
{
	fitp.open(filename.c_str());
}

int Integral::take(struct itg_result &itg_res, E_PRINT_MODE pm)
{
	if (wv.size() < MIN_POINT_COUNT)
	{
		std::cerr << "Error: not enough data to calculate the integral\n";
		return EC_ITG_NOT_ENOUGH_DATA;
	}

	// refresh_cut();

	double len_K = wv[0].mvn.mv.length() / wv[0].mvn.velocity;

	double h = KM2M(cut.v().length()) / n; // шаг в метрах
	
	double dx = (cut.end.x - cut.start.x) / n;
	double dy = (cut.end.y - cut.start.y) / n;

	Line cut_line(cut.v());

	// fitg << "dist(interval) = " << interval.length() << "\tdist_metr(interval) = " << dist_metr(interval) << endl;

	Interpolation itp(cut.v(), wv);

	if (cut.itp_diameter == -1)
		itp.calc_radius();
	if (cut.itp_diameter >= 0.0) 
		itp.set_radius(cut.itp_diameter / 2);
	if (cut.weight_coef >= 0.0) itp.set_weight_coef(cut.weight_coef);

	double lin_sum = 0.0;
	double sqr_sum = 0.0;

	int step_count = 0;

	// fitg << cut.v().toString("station") << endl;
	// fitg << "Point count (for n = " << n << "\th = " << h << "\th_m = " << h_m << "\tdx = " << dx << "\tdy = " << dy << ")\n\n";
	for (int i = 0; i < n; ++i)
	{
		point gr1(cut.start.x + i * dx, cut.start.y + i * dy);
		point gr2(cut.start.x + (i + 1) * dx, cut.start.y + (i + 1) * dy);
		point gr_avr = vec(gr1, gr2).middle();

		double velocity = itp.take_for(gr_avr);

		++step_count;

		vec prnd = cut_line.perpendicular(gr_avr);
		prnd.shorten(velocity * len_K);

		if (pm == EPM_ON)
			fitp << prnd.at_geo_cs(dcs_origin).toGlanceFormat();

		// fitg << "\t\t" << prnd.toString("avr vec") << endl;

		double coriolis = coriolis_koef(gr_avr.at_geo_cs(dcs_origin).y);

		double curv_K = 0.;
		// учёт кривизны потока
		if (cut.curvature_correction == true)
			curv_K = 1 / KM2M(cut.curvature_center.distance_to(gr_avr)); 

		lin_sum += coriolis * velocity * h;
		sqr_sum += curv_K * velocity * velocity * h * sign(velocity);
	}
	// fitg << endl;

	std::vector <double> itp_acr;	// точность интерполяции

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

	itg_res.lin_value = lin_sum;
	itg_res.sqr_value = sqr_sum;

	if (pm == EPM_ON)
	{
		fitp << cut.v().at_geo_cs(dcs_origin).toGlanceFormat();
		fitp.close();
	}

	return EC_ITG_SUCCESS;
}

/*
double Integral::get_integration_error(std::vector <double> &val, double h, int n)
{
	std::vector <double> dv1(val.size() - 1);
	std::vector <double> dv2(dv1.size() - 1);
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