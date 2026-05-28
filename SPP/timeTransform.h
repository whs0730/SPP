#pragma once
#include<cmath>
#include<iostream>
#include<iomanip>
using namespace std;
//define a commontime structure
struct CommonTime
{
	unsigned short Year = 0;
	unsigned short Month = 0;
	unsigned short Day = 0;
	unsigned short Hour = 0;
	unsigned short Minute = 0;
	double		   Second = 0.0;
};
//define a MJDTime structure
struct MJDTime
{
	int Days;
	double FracDay;
	MJDTime()
	{
		Days = 0;
		FracDay = 0.0;
	}
};
//define a GPSTime structure 
struct GPSTime
{
	unsigned short Week;
	double SecOfWeek;
	GPSTime(int week = 0, double secOfWeek = 0.0) // 新增的构造函数
	{
		Week = static_cast<unsigned short> (week);
		SecOfWeek = secOfWeek;
	}
};
// ====================== 时间结构体 ======================
// 便于计算统一时间表示：周号和周内秒
struct  gtime_t{
	int week;    // GPS 周
	double sec;  // 周内秒
};

//将GPST转化为统一时间
static gtime_t gpst2time(int week, double sec) {
	gtime_t t = { week, sec };
	return t;
}

// 计算两个 gtime_t 的时间差（秒）
static double timediff(gtime_t t1, gtime_t t2) {
	return (t1.week - t2.week) * 604800.0 + (t1.sec - t2.sec);
}

// 从 gtime_t 中提取周号（简化实现）
static int time2gpst(gtime_t t, int* week) {
	if (week) *week = t.week;
	return t.sec;
}

// 调整 GPS 周号以处理较小的周号（每到1024就会重新计算，到现在至少两次）
static int adjgpsweek(int week) {
	if (week < 300) week += 2048;
	return week;
}

//define a function to convert the commontime to the MJDTime
inline MJDTime* CommonTimeToMJDTime(CommonTime& common_time)
{
	long double JD = 0.0;
	MJDTime* MJD = new MJDTime();
	unsigned short M = common_time.Month;
	unsigned short Y = common_time.Year;
	unsigned short D = common_time.Day;
	unsigned short Hour = common_time.Hour;
	unsigned short Minute = common_time.Minute;
	double		   Second = common_time.Second;
	double UT = Hour + Minute / 60.0 + Second / 3600.0;
	unsigned short y;
	unsigned short m;
	if (M <= 2) {
		y = Y - 1;
		m = M + 12;
	}
	else {
		y = Y;
		m = M;
	}
	JD = (int)(365.25 * y) + (int)(30.6001 * (m + 1)) + D + UT / 24 + 1720981.5;
	MJD->Days = int(JD - 2400000.5);
	MJD->FracDay = JD - 2400000.5 - MJD->Days;
	return MJD;
}
//define a function to convert the MJDTime to the commontime
inline CommonTime* MJDToCommonTime(MJDTime& MJD)
{
	CommonTime* ct = new CommonTime();

	double jd = MJD.Days + MJD.FracDay + 2400000.5;
	int Z = (int)floor(jd + 0.5);
	double F = (jd + 0.5) - Z;

	int A = Z;
	int alpha = (int)((A - 1867216.25) / 36524.25);
	A = A + 1 + alpha - alpha / 4;

	int B = A + 1524;
	int C = (int)((B - 122.1) / 365.25);
	int D = (int)(365.25 * C);
	int E = (int)((B - D) / 30.6001);

	double day = B - D - (int)(30.6001 * E) + F;

	ct->Day = (int)day;

	if (E < 14) ct->Month = E - 1;
	else        ct->Month = E - 13;

	if (ct->Month > 2) ct->Year = C - 4716;
	else               ct->Year = C - 4715;

	double frac = day - ct->Day;
	double sod = frac * 86400.0;

	ct->Hour = (unsigned short)(sod / 3600.0);
	sod -= ct->Hour * 3600.0;
	ct->Minute = (unsigned short)(sod / 60.0);
	sod -= ct->Minute * 60.0;
	ct->Second = sod;

	return ct;
}

//define a function to convert the MJDTime to the GPSTime
inline GPSTime* MJDToGPSTime(MJDTime& MJD) {
	GPSTime* GPS_Time = new GPSTime();
	GPS_Time->Week = (int)((MJD.Days + MJD.FracDay - 44244) / 7);
	GPS_Time->SecOfWeek = (MJD.Days + MJD.FracDay - 44244 - GPS_Time->Week * 7) * 86400.0;
	return GPS_Time;
}
//define a function to convert the GPSTime  to the MJDTime
inline MJDTime* GPSTimeToMJD(GPSTime& GPSTime) {
	MJDTime* MJD_Time = new MJDTime();
	double MJD = 44244 + GPSTime.Week * 7 + GPSTime.SecOfWeek / 86400.0;
	MJD_Time->Days = (int)(MJD);
	MJD_Time->FracDay = MJD - MJD_Time->Days;
	return MJD_Time;
}
//define a function to convert the GPSTime  to the commontime
inline CommonTime* GPSTimeToCommonTime(GPSTime& GPSTime) {
	MJDTime* MJDTime = GPSTimeToMJD(GPSTime);
	return MJDToCommonTime(*MJDTime);
}
//define a function to convert the commontime  to the GPSTime 
inline GPSTime* CommonTimeToGPSTime(CommonTime& CommonTime) {
	MJDTime* MJDTime = CommonTimeToMJDTime(CommonTime);
	return MJDToGPSTime(*MJDTime);
}

#define BDT_GPS_WEEK_OFFSET 1356  // 2006-01-01 到 1980-01-06 的固定周差
// 把BDT周秒转换为gtime_t。这里先只统一周号基准，秒仍保持BDT周内秒。
// 后续如需和GPST观测历元做timediff，需要再通过bdt2gpst()加14秒。
static gtime_t bdt2time(int week, double sec)
{
	gtime_t t;
	t.week = week + BDT_GPS_WEEK_OFFSET;  // 先把 BDT 周号换到 GPS 周号基准
	t.sec = sec;
	return t;
}

// gtime_t 加秒
static gtime_t timeadd(gtime_t t, double dt)
{
	t.sec += dt;

	while (t.sec < 0.0) {
		t.sec += 604800.0;
		t.week--;
	}
	while (t.sec >= 604800.0) {
		t.sec -= 604800.0;
		t.week++;
	}
	return t;
}

// BDT比GPST慢14秒：GPST = BDT + 14s。
// BDS星历的toe/toc来自BDT，转为GPST后才能和观测历元GPST直接相减。
static gtime_t bdt2gpst(gtime_t t)
{
	return timeadd(t, 14.0);
}

// GPST转BDT：BDT = GPST - 14s。
// 当前程序主要把BDS toe/toc转到GPST参与计算；保留此函数用于说明时间关系和后续扩展。
static gtime_t gpst2bdt(gtime_t t)
{
	return timeadd(t, -14.0);
}
