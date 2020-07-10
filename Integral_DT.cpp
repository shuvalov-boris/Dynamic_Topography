#include <fstream>
#include <io.h>
#include <iostream>
#include <locale>
#include <math.h>
#include <sstream>
#include <string.h>
#include <vector>

#include "dt_tests.h"
#include "dynamic_topography.h"
#include "log.h"

using namespace std;


void print_eng_usage()
{
	cout << "This Software is for calculating dynamic topography between two points\n\n";
	cout << "Integral_DT release data - 10 July 2020\n\n";
	cout << "USAGE:\tintegral_dt.exe <vp_out_file> <boundary_points_list> <dt_out_file>\n\n";

	cout << "Example: ""integral_DT.exe out_2006-05-04_0730_n27799.m.pro_2006-05-04_1300_n70056.m.pro.txt stations.txt DT_out.txt""\n\n"

		<< "\t<vp_out_file>\t\tVecPlotter output file describing flow velocity field\n"
		<< "\t    string format:\n"
		<< "\t\tgeo longitude start\n"
		<< "\t\tgeo latitude start\n"
		<< "\t\tgeo longitude end\n"
		<< "\t\tgeo latitude end\n"
		<< "\t\tpixel longitude start\n"
		<< "\t\tpixel latitude start\n"
		<< "\t\tpixel longitude end\n"
		<< "\t\tpixel latitude end\n"
		<< "\t\tcorrelation\n"
		<< "\t\tvelocity of movement\n"
		<< "\t\ta priori error\n"

		<< "\t<boundary_points_list>\tList of cut definitions\n"
		<< "\t    string format:\n"
		<< "\t\tgeo longitude start\n"
		<< "\t\tgeo latitude start\n"
		<< "\t\tgeo longitude end\n"
		<< "\t\tgeo latitude end\n"
		<< "\t\tcut width (radius), [km]\tor -1 by default (10 km)\n"
		<< "\t\tinterpolation interval (diameter), [km]\tor -1 by default\n"
		<< "\t\tweight coefficient\tor -1 by default\n"

		<< "\t<dt_out_file>\t\tfile for result output\n"
		<< "\t    string format:\n"
		<< "\t\tgeo longitude start\n"
		<< "\t\tgeo latitude start\n"
		<< "\t\tgeo longitude end\n"
		<< "\t\tgeo latitude end\n"
		<< "\t\tDynamic Topography, [meters]\n"
		<< "\t\tcut width (radius), [km]\n"
		<< "\t\tinterpolation interval (diameter), [km]\n"
		<< "\t\tweight coefficient\n"
		<< "\t\tintegration error (K|I(n)-I(2n)|=dDT)\n"
		<< "\t\tinterpolation accuracy\n"
		<< "\t\tintegration error\n"
		<< "\t\tmean-square deviation\n"
		<< "\t\ta priori error\n"
		<< "\t\tcut length, [km]\n"
		<< "\t\tDT coefficient\n"
		<< "\t\tintegration step size [meters]\n"
		<< "\t\tintegration steps count\n"
		<< "\t\tvector count\n"
		<< "\n";
}

char* move_points_file;
char* station_points_file;
char* out_file;
// char* output_log = (char *)"log.txt";
// char* itg_log = (char *)"itg_log.txt";

bool file_exists(const char *fname)
{
	//return !(ifstream(fname) == NULL);
	return _access(fname, 0) != -1;
}

bool is_options_correct(int argc, char** argv)
{
	if (argc > 3)
	{
		move_points_file = argv[1];
		station_points_file = argv[2];
		out_file = argv[3];
		if (strcmp(move_points_file, station_points_file) && file_exists(move_points_file) && file_exists(station_points_file))
			return true;
	}
	// else
	// 	print_usage();
	return false;
}

void read_cuts(char *file_name, vector <scut> &cut)
{
	fstream fcut;
	fcut.open(file_name);
	double gsx, gsy, gex, gey, cut_width, itp_diameter, weight_coef;

	const char *dm = " \t";

	while (fcut >> gsx >> gsy >> gex >> gey >> cut_width >> itp_diameter >> weight_coef)
	{
		vec v(point(gsx, gsy), point(gex, gey));
		cut.push_back(scut(v, cut_width, itp_diameter, weight_coef / 1000.));
	}
	fcut.close();
}

void read_movement_field(char *file_name, vector <movement> &mvn)
{
	fstream fmoves;
	fmoves.open(file_name);
	double gsx, gsy, gex, gey, crl, vlc, err;
	int psx, psy, pex, pey;
	while (fmoves >> gsx >> gsy >> gex >> gey >> psx >> psy >> pex >> pey >> crl >> vlc >> err)
	{
		vec v(point(gsx, gsy), point(gex, gey));
		mvn.push_back(movement(v, vlc, err));
		// getline(fmoves, s);
	}
	fmoves.close();
}

/*wstring AnsiToWide(const string& in_sAnsi)
{
	wstring wsWide;
	wsWide.resize(in_sAnsi.length(), 0);
	MultiByteToWideChar(
		1251,
		0,
		&in_sAnsi[0],
		(int)in_sAnsi.length(),
		&wsWide[0],
		(int)wsWide.length());

	return wsWide;
}*/

int main(int argc, char** argv)
{

	if (!is_options_correct(argc, argv))
	{
		cout << "Incorrect arguments!\n\n";
		print_eng_usage();
	}

	ofstream fres;

	vector <movement> mvn;
	vector <scut> station;

	read_cuts(station_points_file, station);

	read_movement_field(move_points_file, mvn);
	
	if (station.size() == 0)
	{
		cerr << "Error: Station amount is zero\n";
		return 0;
	}


// todo реализовать выполнение через флаг в командной строке
/* test of geo to dec transformation*/
	// test_geo2dec2geo(mvn, station[0].v().middle());
	// test_to_geo_transforms();
	// return 0;
/* end of test of geo to dec transformation*/


	// flog.open(output_log);
	// fitg.open(itg_log);

	point geo_origin = station[0].v().middle();
	to_cartesian_cs(mvn, station);

	ofstream fDSC;
	fDSC.open("DSC.txt");
	for (size_t i = 0; i < mvn.size(); i++)
		fDSC << mvn[i].mv.start.x << " " << mvn[i].mv.start.y << " " << mvn[i].mv.end.x << " " << mvn[i].mv.end.y << endl;

	fres.open(out_file);

	DynamicTopography dyn_tpg(mvn);

	for (size_t i = 0; i < station.size(); ++i)
	{
		dyn_tpg.set_cut(station[i]);
		dyn_tpg.set_file_index(i + 1);
		dyn_tpg.set_dcs_origin(geo_origin);

		struct dt_result dt_res;
		int ce = dyn_tpg.take(dt_res);

		if (ce == EC_DT_SUCCESS)
			dt_res.print_to(fres);
		else
			cerr << "Error: DT taking: " << ce << endl;

	}

	fres.close();

	// flog.close();
	// fitg.close();

	return 0;
}