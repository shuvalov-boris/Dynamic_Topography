#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "dt_defs.h"
#include "log.h"

#include <algorithm>
#include <vector>

#define WEIGHT_COEF 0.01
#define WEIGHT_COEF_TRANSFORM 1. // потребовалось при переходе от градусов к метрам для сохранение прежней размерности WEIGHT_COEF

class Interpolation
{
	double R;
	double weight_coef;
	vec interval;
	vector <wvector> wv;

	double get_norm_comp(int idx);

	double weight_func(double _r);
	double calc_weight_sum(vector <int> &idx, point pt, int omit_idx);
	double get_interpolation_result(vector <int> &idx, point pt, double sum);

public:
	Interpolation(vec itv, vector <wvector> &_wv);

	void calc_radius();
	void set_radius(double r);
	void set_weight_coef(double coef);

	double get_radius();
	double get_weight_coef();

	double take_for(point pt);

	void calc_accuracy(vector <double> &err);

};

double Interpolation::get_norm_comp(int idx)
{
	Line cut_line(interval);
	double ang = cut_line.angle(wv[idx].mvn.mv);
	return wv[idx].mvn.velocity * -sin(ang * M_PI / 180);
}

double Interpolation::weight_func(double r)
{
	return exp(- weight_coef * r * r) - exp(- weight_coef * R * R);
}

double Interpolation::calc_weight_sum(vector <int> &idx, point pt, int omit_idx = -1)
{
	double S = 0.0;
	for (size_t j = 0; j < wv.size(); ++j)
	{
		double r = pt.distance_to(wv[j].proj);
		if ((size_t)omit_idx != j && r <= R /*/ 2*/)
		{
			idx.push_back(j);
			double wf = weight_func(r);
			// cout << "for r = " << r << " weight func is " << wf << endl;
			S += wf;
			// flog << r << " ";
		}
	}
	// flog << endl;

	return S;
}

double Interpolation::get_interpolation_result(vector <int> &idx, point pt, double sum)
{
	if (sum == 0.0) return 0.0;

	double val = 0.0;
	// flog << "sum elemnts\t";
	for (size_t i = 0; i < idx.size(); ++i)
	{
		double r = pt.distance_to(wv[idx[i]].proj);
		// double ang = cut_line.angle(wv[idx[i]].mvn.mv);
		// cout << "get_norm_comp(idx[i]) = " << get_norm_comp(idx[i]) << endl;
		val += get_norm_comp(idx[i]) * weight_func(r) / sum;
		// flog << d << "," << r << ":" << func(r) << "," << func(r) / S << " ";
		// flog << get_norm_comp(idx[i]) * weight_func(r) / sum << " ";
	}
	return val;
}

Interpolation::Interpolation(vec itv, vector <wvector> &_wv) : weight_coef(WEIGHT_COEF / WEIGHT_COEF_TRANSFORM), interval(itv), wv(_wv)
{
	R = itv.length();
}

void Interpolation::calc_radius()
{
	double res;
	vector <double> dist;
	for (size_t i = 0 ; i < wv.size(); ++i)
	{
		double d = interval.start.distance_to(wv[i].proj);
		dist.push_back(d);
	}
	sort(dist.begin(), dist.end());
	
	res = 0.0;
	for (size_t i = 1; i < dist.size(); ++i)
	{
		res = max(res, dist[i] - dist[i - 1]);
	}
	R = res * 1.5;
}

void Interpolation::set_radius(double _r)
{
	R = _r;
}

void Interpolation::set_weight_coef(double coef)
{
	weight_coef = coef * WEIGHT_COEF_TRANSFORM;
}


double Interpolation::get_radius()
{
	return R;
}

double Interpolation::get_weight_coef()
{
	return weight_coef / WEIGHT_COEF_TRANSFORM;
}

double Interpolation::take_for(point pt)
{

	double val = 0.0;
	vector <int> act_p_ind;

	double S = calc_weight_sum(act_p_ind, pt);

	// cout << " interpolation: take for: act_p_ind vector size is " << act_p_ind.size() << " & weight sum is " << S << endl;

	val = get_interpolation_result(act_p_ind, pt, S);

	return val;
}

void Interpolation::calc_accuracy(vector <double> &err)
{
	vector <int> act_p_ind;

	vector <double> vl;
	vector <double> nc;
	vector <double> sum;

	for (size_t i = 0; i < wv.size(); ++i)
	{
		act_p_ind.clear();

		double S = calc_weight_sum(act_p_ind, wv[i].proj, i);
	
		double val = get_interpolation_result(act_p_ind, wv[i].proj, S);

		err.push_back(val - get_norm_comp(i));

		vl.push_back(val);
		nc.push_back(get_norm_comp(i));
		sum.push_back(S);
	}
/*
	cout << "sum: ";
	for (int i = 0; i < sum.size(); ++i)
		cout << sum[i] << " ";
	cout << endl;

	cout << "interpolation: ";
	for (int i = 0; i < vl.size(); ++i)
		cout << vl[i] << " ";
	cout << endl;

	cout << "norm comp: ";
	for (int i = 0; i < nc.size(); ++i)
		cout << nc[i] << " ";
	cout << endl;	

	double _s = 0.0, _ss = 0.0;
	cout << "accuracy: ";
	for (int i = 0; i < nc.size(); ++i)
	{
		cout << vl[i] - nc[i] << " ";
		_s += vl[i] - nc[i];
		_ss += (vl[i] - nc[i]) * (vl[i] - nc[i]);
	}
	cout << endl;

	double s = 0.0;
	for (int i = 0; i < err.size(); ++i)
		s += err[i];

	cout << "err sum is " << s << endl;
	cout << "err size is " << err.size() << endl << endl;

	cout << "sum " << _s << endl;
	cout << "sq sum " << _ss << endl << endl << endl;*/

}


#endif // INTERPOLATION_H