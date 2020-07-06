#ifndef DT_TESTS_H
#define DT_TESTS_H

#include <iostream>
#include <fstream>

#include "dt_defs.h"
#include "dynamic_topography.h"
#include "geometry.h"

void test_geo2dec2geo(vector <movement> &mvn, point dcs_geo_origin)
{
	ofstream f_dec_mvn, f_geo_again;
	f_dec_mvn.open("dec_mvn.txt");
	f_geo_again.open("geo_again.txt");

	unsigned int failed_tests_amount = 0;

	point dcs_dec_origin = geo2dec(dcs_geo_origin);

	for (size_t j = 0; j < mvn.size(); ++j)
	{
		// cout << "geo_x_0 = " << mvn[j].mv.start.x << endl;
		// cout << "geo_y_0 = " << mvn[j].mv.start.y << endl;

		point p0 = geo2dec(mvn[j].mv.start, dcs_dec_origin);
		// cout << mvn[j].mv.start.toString("geo start") << "\t->\t" << p0.toString("dec start") << endl;
		point g0 = dec2geo(p0, dcs_dec_origin);
		
		// cout << endl;
		// cout << "geo_x_1 = " << mvn[j].mv.end.x << endl;
		// cout << "geo_y_1 = " << mvn[j].mv.end.y << endl;

		point p1 = geo2dec(mvn[j].mv.end, dcs_dec_origin);
		// cout << mvn[j].mv.end.toString("geo end") << "\t->\t" << p1.toString("dec staendrt") << endl;
		point g1 = dec2geo(p1, dcs_dec_origin);
		

		f_dec_mvn << fixed << p0.x << " " << p0.y << " " << p1.x << " " << p1.y << endl;
		// << setprecision(n)

		f_geo_again << fixed << g0.x << " " << g0.y << " " << g1.x << " " << g1.y << endl;

		if (mvn[j].mv.start.equal_eps(g0) == false || mvn[j].mv.end.equal_eps(g1) == false)
		{
			cout << "j = " << j << ": FAIL\n";
			if (mvn[j].mv.start.equal(g0) == false)
			{
				cout << fixed << setprecision(10) << mvn[j].mv.start.toString("origin geo 0") << endl;
				cout << g0.toString("transformed geo 0") << endl;

			}

			if (mvn[j].mv.end.equal(g1) == false)
			{
				cout << fixed << setprecision(10) << mvn[j].mv.end.toString("origin geo 1") << endl;
				cout << g1.toString("transformed geo 1") << endl;
			}

			failed_tests_amount++;
		}

		// cout << endl;

	}

	if (failed_tests_amount == 0)
		cout << "geo2dec & dec2geo functions tests -- SUCCESS\n";

	f_dec_mvn.close();
	f_geo_again.close();
}

#endif // DT_TESTS_H