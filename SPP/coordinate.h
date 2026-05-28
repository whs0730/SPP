#pragma once
#include<iostream>
#include<cmath>
#include<iomanip>
using namespace std;
#define PI 3.1415926535898 
//define a Cartesian coordinate system structure
struct XYZ
{
	double X = 0.0;//m
	double Y = 0.0;
	double Z = 0.0;
};
//define a Geodetic coordinate system structure
struct BLH
{
	//radian system
	double B = 0.0;
	double L = 0.0;
	double H = 0.0;//m
};
struct ENU
{
	double E = 0.0;
	double N = 0.0;
	double U = 0.0;
};
//define a finction to convert BLH to XYZ
inline XYZ* BLHtoXYZ(BLH& blh, double a /*long axis of ellipsoid in WGS-84*/, double f/*flattening of ellipsoid in WGS-84*/) {
	double e_2 = 2 * f - f * f;
	double sinb = sin(blh.B);
	double cosb = cos(blh.B);
	double sinl = sin(blh.L);
	double cosl = cos(blh.L);
	double N = a / sqrt(1 - e_2 * sinb * sinb);
	XYZ* xyz = new XYZ();
	xyz->X = (N + blh.H) * cosb * cosl;
	xyz->Y = (N + blh.H) * cosb * sinl;
	xyz->Z = (N * (1 - e_2) + blh.H) * sinb;
	return xyz;
}
//define a function to convert XYZ to BLH
//define a function to calculate B
inline double calculate_B(XYZ& xyz, double a /*long axis of ellipsoid in CGCS-2000*/, double f/*flattening of ellipsoid in CGCS-2000*/) {
	double e_2 = 2 * f - f * f;
	double p = sqrt(xyz.X * xyz.X + xyz.Y * xyz.Y);
	double tanb_0 = xyz.Z / p;
	while (true)
	{
		double tanb = (xyz.Z + (a * e_2 * tanb_0) / (sqrt(1 + (1 - e_2) * tanb_0 * tanb_0))) / p;
		if (abs(tanb - tanb_0) < 5e-10) {
			return atan(tanb);
		}
		tanb_0 = tanb;
	}
}
inline BLH* XYZtoBLH(XYZ& xyz, double a /*long axis of ellipsoid in CGCS-2000*/, double f/*flattening of ellipsoid in CGCS-2000*/) {
	double e_2 = 2 * f - f * f;
	BLH* blh = new BLH();
	blh->B = calculate_B(xyz, a, f);
	blh->L = atan2(xyz.Y, xyz.X);
	double N = a / sqrt(1 - e_2 * sin(blh->B) * sin(blh->B));
	blh->H = sqrt(xyz.X * xyz.X + xyz.Y * xyz.Y) / cos(blh->B) - N;
	return blh;
}
inline ENU* xyz2enu(BLH* blh, double* e)
{
	ENU* enu = new ENU;

	double sinp = sin(blh->B), cosp = cos(blh->B);
	double sinl = sin(blh->L), cosl = cos(blh->L);
	double E[9];

	E[0] = -sinl;           E[1] = cosl;           E[2] = 0.0;
	E[3] = -sinp * cosl;    E[4] = -sinp * sinl;    E[5] = cosp;
	E[6] = cosp * cosl;    E[7] = cosp * sinl;    E[8] = sinp;

	enu->E = E[0] * e[0] + E[1] * e[1] + E[2] * e[2];
	enu->N = E[3] * e[0] + E[4] * e[1] + E[5] * e[2];
	enu->U = E[6] * e[0] + E[7] * e[1] + E[8] * e[2];

	return enu;
}
inline double satazel(BLH*blh,double*e, double* azel)
{
	ENU* enu=new ENU;

	enu=xyz2enu(blh, e);
	azel[0] = atan2(enu->E, enu->N);
	if (azel[0] < 0.0) azel[0] += 2.0 * PI;
	azel[1] = asin(enu->U);
	delete enu;
	return azel[1];
}
