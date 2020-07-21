#ifndef DT_TESTS_H
#define DT_TESTS_H

#include <iostream>
#include <fstream>

#include "dt_defs.h"
#include "geometry.h"

void test_to_geo_transforms()
{
	point origin = geo2dec(point(148.382800, 42.621050));
	point pt(148.883000, 43.123100);
	point dec_pt = geo2dec(pt, origin);
	
	point pt2(dec_pt.x, dec_pt.y);

	pt2.to_geo_cs(origin);
	
	point pt1 = dec2geo(dec_pt, origin);

	if (pt1.equal_eps(pt2) == true)
		cout << "to geo CS different transforms tests\t\tSUCCESS";
	else
	{
		cout << "to geo CS transforms tests FAILED";
		cout << pt.toString("input point") << endl;
		cout << pt1.toString("dec2geo result") << endl;
		cout << pt2.toString("to_geo_cs result") << endl;
	}

}

void test_geo2dec2geo(vector <movement> &mvn, point dcs_geo_origin)
{
	ofstream f_dec_mvn, f_geo_again;
	f_dec_mvn.open("at_dec_cs.txt");
	f_geo_again.open("at_geo_cs.txt");

	unsigned int failed_tests_amount = 0;

	for (size_t j = 0; j < mvn.size(); ++j)
	{
		point p0 = mvn[j].mv.start.at_dec_cs(dcs_geo_origin);
		point g0 = p0.at_geo_cs(dcs_geo_origin);
		
		point p1 = mvn[j].mv.end.at_dec_cs(dcs_geo_origin);
		point g1 = p1.at_geo_cs(dcs_geo_origin);		

		f_dec_mvn << fixed << p0.x << " " << p0.y << " " << p1.x << " " << p1.y << endl;
		f_geo_again << fixed << g0.x << " " << g0.y << " " << g1.x << " " << g1.y << endl;

		if (mvn[j].mv.start.equal_eps(g0) == false || mvn[j].mv.end.equal_eps(g1) == false)
		{
			cout << "j = " << j << ": FAIL\n";
			if (mvn[j].mv.start.equal_eps(g0) == false)
			{
				cout << fixed << setprecision(10) << mvn[j].mv.start.toString("origin geo 0 ") << endl;
				cout << g0.toString("transformed geo 0 ") << endl;
			}

			if (mvn[j].mv.end.equal_eps(g1) == false)
			{
				cout << fixed << setprecision(10) << mvn[j].mv.end.toString("origin geo 1 ") << endl;
				cout << g1.toString("transformed geo 1 ") << endl;
			}

			failed_tests_amount++;
			break;
		}
	}

	cout << "geo2dec & dec2geo transformations test -- " << ((failed_tests_amount == 0) ? "SUCCESS" : "FAIL") << "\n";

	f_dec_mvn.close();
	f_geo_again.close();
}

#endif // DT_TESTS_H