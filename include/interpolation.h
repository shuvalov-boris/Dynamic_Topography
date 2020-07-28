#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "dt_defs.h"

#include <algorithm>
#include <vector>

#define WEIGHT_COEF 0.01
#define WEIGHT_COEF_TRANSFORM 1. // потребовалось при переходе от градусов к метрам для сохранение прежней размерности WEIGHT_COEF

class Interpolation
{
	double R;
	double weight_coef;
	vec interval;
	std::vector <wvector> wv;

	double get_norm_comp(int idx);

	double weight_func(double _r);
	double calc_weight_sum(std::vector <int> &idx, point pt, int omit_idx);
	double get_interpolation_result(std::vector <int> &idx, point pt, double sum);

public:
	Interpolation(vec itv, std::vector <wvector> &_wv);

	void calc_radius();
	void set_radius(double r);
	void set_weight_coef(double coef);

	double get_radius();
	double get_weight_coef();

	double take_for(point pt);

	void calc_accuracy(std::vector <double> &err);

};

#endif // INTERPOLATION_H