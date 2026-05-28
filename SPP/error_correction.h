#pragma once
#include "timeTransform.h"
#include "coordinate.h"
#include<cmath>
#include "define.h"
#include "obs.h"
const double D2R = PI / 180.0;
const double R2D = 180.0 / PI;
//对流层延迟改正误差
inline double Hopfield_delTrop(double H/*测站高度*/, double E/*卫星相对测站的高度角*/) {
    if (H < -1000.0 || H > 50000.0 || E <= 0.0) {
        return 0.0;
    }
    //标准气象元素
    double H0 = 0.0;//m 海平面H
    double T0 = 15.0 + 273.16;//K 温度
    double P0 = 1013.25;//mbar 气压
    double RH0 = 0.5;//相对湿度
   
	double hd = 40136.0 + 148.72 * (T0 - 273.16);//m
	double hw = 11000.0;//m
	double T = T0 - 0.0065 * (H - H0);
	double P = P0 * pow((1 - 0.0000226 * (H - H0)), 5.225);
	double RH = RH0 * exp(-0.0006396 * (H - H0));
	double e = RH * exp(-37.2465 + 0.213166 * T - 0.000256908 * T * T);
	double Kd = 155.2e-7 * P * (hd - H) / T;
	double Kw = 155.2e-7 * 4810 * e * (hw - H)/(T*T);
    // Hopfield 映射函数里的 E 是角度
    double elev_deg =E * 180.0 / PI;

    double md = 1.0 / sin(sqrt(elev_deg * elev_deg + 6.25) * PI / 180.0);
    double mw = 1.0 / sin(sqrt(elev_deg * elev_deg + 2.25) * PI / 180.0);

    return Kd * md + Kw * mw;
}
//电离层改正误差
inline double Klobuchar(gtime_t t, const double* pos,const double* azel, double* ionpara)
{
    
    if (!pos || !azel || !ionpara) return 0.0;
    if (azel[1] <= 0.0) return 0.0;

    // 若参数全0，可视为无效
    int valid = 0;
    for (int i = 0; i < 8; i++) {
        if (fabs(ionpara[i]) > 1E-30) {
            valid = 1;
            break;
        }
    }
    if (!valid) return 0.0;

    double psi, phi, lam, phi_m;
    double tt, f, amp, per, x;
    double az = azel[0];
    double el = azel[1];

    // 半圆单位大地坐标semicircle
    double phi_u = pos[0] / PI;
    double lam_u = pos[1] / PI;
    double E = el / PI;
    double A = az / PI;

    psi = 0.0137 / (E + 0.11) - 0.022;//信号穿过电离层时偏离正上方的角度
    //电离层穿刺点（IPP）
    phi = phi_u + psi * cos(az);
    if (phi > 0.416) phi = 0.416;
    else if (phi < -0.416) phi = -0.416;
    lam = lam_u + psi * sin(az) / cos(phi * PI);
    //地磁纬度
    phi_m = phi + 0.064 * cos((lam - 1.617) * PI);
    //地方时
    int week = 0;
    double sow = time2gpst(t, &week);
    tt = 43200.0 * lam + sow;
    while (tt >= 86400.0) tt -= 86400.0;
    while (tt < 0.0)      tt += 86400.0;
    //振幅
    amp = ionpara[0] + phi_m * (ionpara[1] + phi_m * (ionpara[2] + phi_m * ionpara[3]));
    if (amp < 0.0) amp = 0.0;
    //周期
    per = ionpara[4] + phi_m * (ionpara[5] + phi_m * (ionpara[6] + phi_m * ionpara[7]));
    if (per < 72000.0) per = 72000.0;
    //延迟计算
    x = 2.0 * PI * (tt - 50400.0) / per;
    //斜路径映射
    f = 1.0 + 16.0 * pow(0.53 - E, 3.0);
    double delay_sec;
    if (fabs(x) < 1.57) {
        //白天
        delay_sec = f * (5E-9 + amp * (1.0 - x * x / 2.0 + x * x * x * x / 24.0));
    }
    else {
        //晚上
        delay_sec = f * 5E-9;
    }

    return Clight * delay_sec; // 单位 m
}
// =====================================================
// 地球自转改正：把发射时刻的卫星坐标旋转到接收时刻对应的地固系。
// tau 为信号传播时间，单位 s。
// =====================================================
static void ApplyEarthRotationCorrection(satpos_t& sat, double tau)
{
    double theta = OMGed_GPS * tau;

    double x = sat.pos[0];
    double y = sat.pos[1];

    sat.pos[0] = cos(theta) * x + sin(theta) * y;
    sat.pos[1] = -sin(theta) * x + cos(theta) * y;
    // Z 不变
}

// =====================================================
// SPP内部旧版近似改正：
// 从接收时刻卫星状态近似外推到发射时刻，再做地球自转改正。
// 第4步后主流程已优先在ComputeSatStateAtTransmitTime()中直接计算发射时刻卫星状态。
// 保留此函数用于说明旧流程和可能的调试对照。
// =====================================================
static void ApplySagnacCorrectionForSPP(satpos_t& sat, double tau)
{
    // tau: 信号传播时间，单位 s

    // 1. 近似改到信号发射时刻
    sat.pos[0] -= sat.vel[0] * tau;
    sat.pos[1] -= sat.vel[1] * tau;
    sat.pos[2] -= sat.vel[2] * tau;

    // 卫星钟差也近似改到发射时刻
    sat.clk -= sat.dclk * tau;

    // 2. 地球自转改正
    ApplyEarthRotationCorrection(sat, tau);
}
