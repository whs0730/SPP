#pragma once
#include<cmath>
#include "obs.h"
#include "matrix.h"
#include "define.h"
#define pi 3.1415926535898 
#define F -4.442807633e-10 //s*m-1/2
#define mu_gps 3.986005e14 //m3/s22 WGS84 value of the earth's gravitational constant for GPS user
#define mu_bds 3.986004418e14 //m3/s22 CGCS2000 value of the earth's gravitational constant for BDS user


void CalculateGPS(gtime_t t, eph_t* eph, satpos_t* sat) {
	// GPS星历和观测历元均按GPST组织，直接用GPST计算tk和钟差。
	//计算位置
	double e = eph->e;
	double n0 = sqrt(mu_gps) / pow(eph->sqrtA, 3);
	double tk = timediff(t, eph->toe);//相对于星历参考历元的时间
	if (tk > 302400) {
		tk = tk - 604800;
	}
	else if (tk < -302400) {
		tk = tk + 604800;
	}
	else {
		tk = tk;
	}
	double n = n0 + eph->deln;//平均运动角速度修正
	double Mk = eph->M0 + n * tk;//平近点角
	double Ek = Mk;
	for (int iter = 0; iter < 30; iter++) {
		double E_old = Ek;
		Ek = Mk + e * sin(E_old);

		if (!isfinite(Ek)) {
			sat->pos[0] = sat->pos[1] = sat->pos[2] = 0.0;
			sat->vel[0] = sat->vel[1] = sat->vel[2] = 0.0;
			sat->clk = sat->dclk = 0.0;
			return;
		}

		if (fabs(Ek - E_old) < 1e-12) {
			break;
		}

		if (iter == 29) {
			sat->pos[0] = sat->pos[1] = sat->pos[2] = 0.0;
			sat->vel[0] = sat->vel[1] = sat->vel[2] = 0.0;
			sat->clk = sat->dclk = 0.0;
			return;
		}
	}

	//偏近点角
	double fenZi = sqrt(1 - pow(e, 2)) * sin(Ek) / (1 - e * cos(Ek));
	double fenMu = (cos(Ek) - e) / (1 - e * cos(Ek));
	double vk = atan2(fenZi, fenMu);//真近点角
	double PHIk = vk + eph->omg;//升交角距
	//二阶调和改正数
	double phiuk = eph->cus * sin(2 * PHIk) + eph->cuc * cos(2 * PHIk);
	double phirk = eph->crs * sin(2 * PHIk) + eph->crc * cos(2 * PHIk);
	double phiik = eph->cis * sin(2 * PHIk) + eph->cic * cos(2 * PHIk);
	//改正后的数据
	double uk = PHIk + phiuk;//改正的升交角距
	double rk = eph->A * (1 - e * cos(Ek)) + phirk;//改正的向径
	double ik = eph->i0 + phiik + eph->idot * tk;//改正的轨道倾角
	//轨道平面位置
	double xk_prime = rk * cos(uk);
	double yk_prime = rk * sin(uk);
	double OMGk = eph->OMG0 + (eph->OMGd - OMGed_GPS) * tk - OMGed_GPS * (eph->toe.sec);
	//计算地心地固系下的坐标
	sat->pos[0] = xk_prime * cos(OMGk) - yk_prime * cos(ik) * sin(OMGk);
	sat->pos[1] = xk_prime * sin(OMGk) + yk_prime * cos(ik) * cos(OMGk);
	sat->pos[2] = yk_prime * sin(ik);
	//计算速度
	double Ekd = n / (1 - e * cos(Ek));
	double PHIkd = sqrt((1 + e) / (1 - e)) * pow(cos(vk / 2) / cos(Ek / 2), 2) * Ekd;
	double ukd = 2 * (eph->cus * cos(2 * PHIk) - eph->cuc * sin(2 * PHIk)) * PHIkd + PHIkd;
	double rkd = eph->A * e * sin(Ek) * Ekd + 2 * (eph->crs * cos(2 * PHIk) - eph->crc * sin(2 * PHIk)) * PHIkd;
	double ikd = eph->idot + 2 * (eph->cis * cos(2 * PHIk) - eph->cic * sin(2 * PHIk)) * PHIkd;
	double OMGkd = eph->OMGd - OMGed_GPS;
	Matrix Rd = {
		{cos(OMGk),-sin(OMGk) * cos(ik),-(xk_prime * sin(OMGk) + yk_prime * cos(OMGk) * cos(ik)),yk_prime * sin(OMGk) * sin(ik)},
		{sin(OMGk),cos(OMGk) * cos(ik),xk_prime * cos(OMGk) - yk_prime * sin(OMGk) * cos(ik),yk_prime * cos(OMGk) * sin(ik)},
		{0,sin(ik),0,yk_prime * cos(ik)}
	};
	double xk_primed = rkd * cos(uk) - rk * ukd * sin(uk);
	double yk_primed = rkd * sin(uk) + rk * ukd * cos(uk);
	Matrix M = {
		{xk_primed},
		{yk_primed},
		{OMGkd},
		{ikd}
	};
	Matrix speed = mul(Rd, M);
	sat->vel[0] = speed[0][0];
	sat->vel[1] = speed[1][0];
	sat->vel[2] = speed[2][0];
	//计算钟差和钟速
	double delt = timediff(t, eph->toc);
	double deltr = F * e * eph->sqrtA * sin(Ek);//相对论改正
	sat->clk = eph->af0 + eph->af1 * delt + eph->af2 * delt * delt + deltr;//计算钟差Δt=a0​+a1​(t−toc​)+a2​(t−toc​)2+Δtr​
	double deltrd = F * e * eph->sqrtA * cos(Ek) * Ekd;
	sat->dclk = eph->af1 + 2.0 * eph->af2 * delt + deltrd;
}


void CalculateBDS(gtime_t t, eph_t* eph, satpos_t* sat) {
	// 传入的t为观测历元GPST。
	// BDS星历解码时已把toe/toc由BDT转换到GPST，因此timediff(t, eph->toe/toc)可直接使用；
	// 同时eph->toes保留原始BDT周内秒，用于BDS升交点经度中的地球自转项。
	//计算位置
	double e = eph->e;
	double n0 = sqrt(mu_bds) / pow(eph->sqrtA, 3);
	double tk = timediff(t, eph->toe);//相对于星历参考历元的时间
	if (tk > 302400) {
		tk = tk - 604800;
	}
	else if (tk < -302400) {
		tk = tk + 604800;
	}
	else {
		tk = tk;
	}
	double n = n0 + eph->deln;//平均运动角速度修正
	double Mk = eph->M0 + n * tk;//平近点角
	double Ek = Mk;
	for (int iter = 0; iter < 30; iter++) {
		double E_old = Ek;
		Ek = Mk + e * sin(E_old);

		if (!isfinite(Ek)) {
			sat->pos[0] = sat->pos[1] = sat->pos[2] = 0.0;
			sat->vel[0] = sat->vel[1] = sat->vel[2] = 0.0;
			sat->clk = sat->dclk = 0.0;
			return;
		}

		if (fabs(Ek - E_old) < 1e-12) {
			break;
		}

		if (iter == 29) {
			sat->pos[0] = sat->pos[1] = sat->pos[2] = 0.0;
			sat->vel[0] = sat->vel[1] = sat->vel[2] = 0.0;
			sat->clk = sat->dclk = 0.0;
			return;
		}
	}
	//偏近点角
	double fenZi = sqrt(1 - pow(e, 2)) * sin(Ek) / (1 - e * cos(Ek));
	double fenMu = (cos(Ek) - e) / (1 - e * cos(Ek));
	double vk = atan2(fenZi, fenMu);//真近点角
	double PHIk = vk + eph->omg;//升交角距
	//二阶调和改正数
	double phiuk = eph->cus * sin(2 * PHIk) + eph->cuc * cos(2 * PHIk);
	double phirk = eph->crs * sin(2 * PHIk) + eph->crc * cos(2 * PHIk);
	double phiik = eph->cis * sin(2 * PHIk) + eph->cic * cos(2 * PHIk);
	//改正后的数据
	double uk = PHIk + phiuk;//改正的升交角距
	double rk = eph->A * (1 - e * cos(Ek)) + phirk;//改正的向径
	double ik = eph->i0 + phiik + eph->idot * tk;//改正的轨道倾角
	//轨道平面位置
	double xk_prime = rk * cos(uk);
	double yk_prime = rk * sin(uk);
	//计算速度所需参数
	double Ekd = n / (1 - e * cos(Ek));
	double vkd = sqrt(1 - e * e) * Ekd / (1 - e * cos(Ek));
	double PHIkd = vkd;
	double ukd = 2 * (eph->cus * cos(2 * PHIk) - eph->cuc * sin(2 * PHIk)) * PHIkd + PHIkd;
	double rkd = eph->A * e * sin(Ek) * Ekd + 2 * (eph->crs * cos(2 * PHIk) - eph->crc * sin(2 * PHIk)) * PHIkd;
	double ikd = eph->idot + 2 * (eph->cis * cos(2 * PHIk) - eph->cic * sin(2 * PHIk)) * PHIkd;
	//计算MEO和IGSO卫星的位置和速度
	if (eph->prn > 5 && eph->prn < 59) {
		//位置
		double OMGk = eph->OMG0 + (eph->OMGd - OMGed_BDS) * tk - OMGed_BDS * (eph->toes);
		//计算地心地固系下的坐标
		double Xk = xk_prime * cos(OMGk) - yk_prime * cos(ik) * sin(OMGk);
		double Yk = xk_prime * sin(OMGk) + yk_prime * cos(ik) * cos(OMGk);
		double Zk = yk_prime * sin(ik);
		sat->pos[0] = Xk;
		sat->pos[1] = Yk;
		sat->pos[2] = Zk;
		//速度
		double OMGkd = eph->OMGd - OMGed_BDS;
		double xkd = rkd * cos(uk) - rk * ukd * sin(uk);
		double ykd = rkd * sin(uk) + rk * ukd * cos(uk);
		sat->vel[0] = -Yk * OMGkd - (ykd * cos(ik) - Zk * ikd) * sin(OMGk) + xkd * cos(OMGk);
		sat->vel[1] = Xk * OMGkd + (ykd * cos(ik) - Zk * ikd) * cos(OMGk) + xkd * sin(OMGk);
		sat->vel[2] = ykd * sin(ik) + yk_prime * ikd * cos(ik);

	}
	//计算GEO卫星的位置和速度
	else {
		//计算位置
		double OMGk = eph->OMG0 + eph->OMGd * tk - OMGed_BDS * (eph->toes);
		//计算地心地固系下的坐标
		double XGK = xk_prime * cos(OMGk) - yk_prime * cos(ik) * sin(OMGk);
		double YGK = xk_prime * sin(OMGk) + yk_prime * cos(ik) * cos(OMGk);
		double ZGK= yk_prime * sin(ik);
		double phi1 = OMGed_BDS *tk;
		double phi2 = -5 * pi / 180.0;
		Matrix Rz = {
			{cos(phi1),sin(phi1),0},
			{-sin(phi1),cos(phi1),0},
			{0,0,1}
		};
		Matrix Rx = {
			{1,0,0},
			{0,cos(phi2),sin(phi2)},
			{0,-sin(phi2),cos(phi2)}
		};
		Matrix XYZGK = {
			{XGK},
			{YGK},
			{ZGK}
		};
		Matrix XYZ = mul(mul(Rz, Rx), XYZGK);
		sat->pos[0] = XYZ[0][0];
		sat->pos[1] = XYZ[1][0];
		sat->pos[2] = XYZ[2][0];
		//计算速度
		double OMGkd = eph->OMGd;
		double xkd = rkd * cos(uk) - rk * ukd * sin(uk);
		double ykd = rkd * sin(uk) + rk * ukd * cos(uk);
		double XGKd = -YGK * OMGkd - (ykd * cos(ik) - ZGK * ikd) * sin(OMGk) + xkd * cos(OMGk);
		double YGKd = XGK * OMGkd + (ykd * cos(ik) - ZGK * ikd) * cos(OMGk) + xkd * sin(OMGk);
		double ZGKd = ykd * sin(ik) + yk_prime * ikd * cos(ik);
		Matrix Rzd = {
			{-sin(phi1) * OMGed_BDS,cos(phi1) * OMGed_BDS,0},
			{-cos(phi1) * OMGed_BDS,-sin(phi1) * OMGed_BDS,0},
			{0,0,0}
		};
		Matrix XYZGKd = {
			{XGKd},
			{YGKd},
			{ZGKd}
		};
		Matrix XYZkd = add(mul(mul(Rzd, Rx), XYZGK), mul(mul(Rz, Rx), XYZGKd));
		sat->vel[0] = XYZkd[0][0];
		sat->vel[1] = XYZkd[1][0];
		sat->vel[2] = XYZkd[2][0];

	}
	//计算钟差和钟速
	double delt = timediff(t, eph->toc);
	double deltr = F * e * eph->sqrtA * sin(Ek);//相对论改正
	sat->clk = eph->af0 + eph->af1 * delt + eph->af2 * delt * delt + deltr;//计算钟差Δt=a0​+a1​(t−toc​)+a2​(t−toc​)2+Δtr​
	double deltrd = F * e * eph->sqrtA * cos(Ek) * Ekd;
	sat->dclk = eph->af1 + 2.0 * eph->af2 * delt + deltrd;
}
