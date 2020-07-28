#ifndef INTEGRATION_H
#define INTEGRATION_H

#include "dt_defs.h"
#include "interpolation.h"

#include <string.h>
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


class Integral
{
	scut cut;
	point dcs_origin;
	std::vector <wvector> wv;

	int n; // количество интервалов разбиения

	E_PRINT_MODE print_mode;
	std::ofstream fitp;

	double get_integration_error(std::vector <double> &val, double h, int n);

public:
	Integral(scut c, std::vector <wvector> &_wv);

	void set_cut(scut c);
	void set_dcs_origin(const point &dcs_orn);
	void set_partitioning_count(int _n);
	void set_filename(std::string filename);

	int take(struct itg_result &itg_res, E_PRINT_MODE pm = EPM_OFF);
};



#endif // INTEGRATION_H