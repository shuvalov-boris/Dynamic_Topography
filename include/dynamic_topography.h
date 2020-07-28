#ifndef DYNAMIC_TOPOGRAPHY_H
#define DYNAMIC_TOPOGRAPHY_H

#include "dt_defs.h"
#include "integration.h"

#include <vector>
#include <fstream>
#include <iostream>

#define CUT_WIDTH 10 // [км]

// коды ошибок при расчете перепада динамических высот
#define EC_DT_SUCCESS 1000
#define EC_DT_FVF_EMPTY 1002

struct dt_result
{
	scut cut; 					// разрез
	double dt; 					// динамическая топография
	struct itg_result itg_res;
	double dt_error;			// разность двух значени ДТ с разным количеством интервалов
	double a_priori_error; 		// апроиорная ошибка
	double cut_length;
	double dt_coef; // 
	double cr_coef;
	int vector_count;

	dt_result();

	void set(scut _cut, int vc);

	void calc_dt(double latitude);

	void print_to(std::ofstream &file);
};

class DynamicTopography
{
	scut cut;
	point dcs_origin; // начало локальной Декартовой СК в географических координатах
	std::vector <movement> mvn;

	std::ofstream fNV;
	int file_index;

public:
	DynamicTopography(std::vector <movement> &m);

	void set_file_index(int index);
	void set_cut(scut c);
	void set_dcs_origin(const point &dcs_orn);
	int take(struct dt_result &dt_res);

};

#endif //DYNAMIC_TOPOGRAPHY_H