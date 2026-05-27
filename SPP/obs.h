#pragma once
#include "timeTransform.h"
#include "define.h"
#include "matrix.h"
// 观测数据结构
struct obsd_t {
    gtime_t time;                           // 观测时间
    int sat;                                // 卫星编号
    double L[NFREQ + NEXOBS];               // 相位观测（周）
    double P[NFREQ + NEXOBS];               // 伪距观测（米）
    float D[NFREQ + NEXOBS];                // 多普勒/码率
    double SNR[NFREQ + NEXOBS];             // 信噪比（dB-Hz）
    unsigned char LLI[NFREQ + NEXOBS];      // 周跳标记
    unsigned char code[NFREQ + NEXOBS];     // 观测码类型
} ;
// 观测集合，包含多个卫星观测
struct obs_t {
    int n;                  // 当前观测数
    obsd_t data[MAXOBS];    // 观测数组
} ;
//导航电文结构
struct eph_t {
    int sat;        // 内部卫星号
    int prn;        // 系统内PRN号
    int week;       // 周号
    double toes;   // 原始toe秒（BDT周内秒）
    double tocs;   // 原始toc秒 
    gtime_t toe, toc;

    // ===== 轨道参数 =====
    double sqrtA;
    double A;//轨道长半轴
    double e;//偏心率
    double M0;//参考时刻的平近点角
    double OMG0;//参考时刻toe升交点赤经OMGoe与GPS周开始时格尼尼治赤经GAST只差
    double omg;//近地点角距
    double i0;//参考时刻的轨道倾角

    double deln;
    double OMGd;
    double idot;

    // ===== 摄动参数 =====
    double cuc, cus;//升交角距的余弦和正弦调和改正的振幅
    double crc, crs;//轨道半径的余弦和正弦调和改正的振幅
    double cic, cis;//轨道倾角的余弦和正弦调和改正的振幅

    // ===== 钟差参数 =====
    double af0, af1, af2;
    //TGD群延迟
    double tgd[2];

} ;
// 简化的导航数据结构（电文、离子/UTC、GLONASS 频率通道）
struct nav_t {
    eph_t eph[MAXSAT];
    double ion_gps[8];
    double utc_gps[8];
    int glo_fcn[32];
} ;
//卫星计算结果
struct satpos_t{
    int sat;
    gtime_t time;

    double pos[3];//大地坐标系坐标
    double vel[3];//速度

    double clk;//钟差
    double dclk;//钟速

} ;
//定位结果
struct sol_t{
    gtime_t time;      // 当前历元时间
    double XYZ[3];      // 接收机ECEF坐标 (m)
    double dtr[2];     // 接收机钟差 (m)
                         // dtr[0] = GPS
                         // dtr[1] = BDS（如果不用可以只用一个）
    Matrix Q;       // 协方差矩阵（可选）

    double pdop;       // PDOP
    double sigma0;     // 单位权中误差
    int ns;            // 参与解算卫星数
    int stat;          // 0失败 1成功

} ;
struct solvel_t {
    gtime_t time;       // 当前历元时间

    double V[3];        // 接收机ECEF速度，单位 m/s
                        // V[0] = Vx
                        // V[1] = Vy
                        // V[2] = Vz

    double dtrd;        // 接收机钟速，单位 m/s

    Matrix Q;           // 协方差矩阵
    double sigma0;      // 单位权中误差
    int ns;             // 参与测速卫星数
    int stat;       // 0失败 1成功

};
//根据系统类型和内部prn编号
static int satno(int sys, int prn)
{
    if (prn <= 0) return 0;
    switch (sys) {
    case SYS_GPS:
        if (prn < MINPRNGPS || MAXPRNGPS < prn) return 0;
        return prn - MINPRNGPS + 1;
    case SYS_GLO:
        if (prn < MINPRNGLO || MAXPRNGLO < prn) return 0;
        return NSATGPS + prn - MINPRNGLO + 1;
    case SYS_GAL:
        if (prn < MINPRNGAL || MAXPRNGAL < prn) return 0;
        return NSATGPS + NSATGLO + prn - MINPRNGAL + 1;
    case SYS_QZS:
        if (prn < MINPRNQZS || MAXPRNQZS < prn) return 0;
        return NSATGPS + NSATGLO + NSATGAL + prn - MINPRNQZS + 1;
    case SYS_CMP:
        if (prn < MINPRNCMP || MAXPRNCMP < prn) return 0;
        return NSATGPS + NSATGLO + NSATGAL + NSATQZS + prn - MINPRNCMP + 1;
    case SYS_IRN:
        if (prn < MINPRNIRN || MAXPRNIRN < prn) return 0;
        return NSATGPS + NSATGLO + NSATGAL + NSATQZS + NSATCMP + prn - MINPRNIRN + 1;
    case SYS_LEO:
        if (prn < MINPRNLEO || MAXPRNLEO < prn) return 0;
        return NSATGPS + NSATGLO + NSATGAL + NSATQZS + NSATCMP + NSATIRN +
            prn - MINPRNLEO + 1;
    case SYS_SBS:
        if (prn < MINPRNSBS || MAXPRNSBS < prn) return 0;
        return NSATGPS + NSATGLO + NSATGAL + NSATQZS + NSATCMP + NSATIRN + NSATLEO +
            prn - MINPRNSBS + 1;
    }
    return 0;
}
// 根据内部sat编号反推出系统和PRN
static int satsys(int sat, int* prn)
{
    if (sat <= 0) return SYS_NONE;

    if (sat <= NSATGPS) {
        if (prn) *prn = sat;
        return SYS_GPS;
    }
    else if (sat <= NSATGPS + NSATGLO) {
        if (prn) *prn = sat - NSATGPS;
        return SYS_GLO;
    }
    else if (sat <= NSATGPS + NSATGLO + NSATGAL) {
        if (prn) *prn = sat - NSATGPS - NSATGLO;
        return SYS_GAL;
    }
    else if (sat <= NSATGPS + NSATGLO + NSATGAL + NSATQZS) {
        if (prn) *prn = sat - NSATGPS - NSATGLO - NSATGAL;
        return SYS_QZS;
    }
    else if (sat <= NSATGPS + NSATGLO + NSATGAL + NSATQZS + NSATCMP) {
        if (prn) *prn = sat - NSATGPS - NSATGLO - NSATGAL - NSATQZS;
        return SYS_CMP;
    }
    else if (sat <= NSATGPS + NSATGLO + NSATGAL + NSATQZS + NSATCMP + NSATIRN) {
        if (prn) *prn = sat - NSATGPS - NSATGLO - NSATGAL - NSATQZS - NSATCMP;
        return SYS_IRN;
    }
    return SYS_NONE;
}
//生成卫星名称
static string GetSatName(int sat)
{
    int prn = 0;
    int sys = satsys(sat, &prn);

    if (sys == SYS_GPS) {
        return "G" + string(prn < 10 ? "0" : "") + to_string(prn);
    }
    else if (sys == SYS_CMP) {
        return "C" + string(prn < 10 ? "0" : "") + to_string(prn);
    }

    return "";
}
// =====================================================
// 判断卫星位置是否有效
// =====================================================
static bool IsValidSatpos(satpos_t* sat)
{
    if (sat == nullptr) {
        return false;
    }

    if (sat->pos[0] == 0.0 && sat->pos[1] == 0.0 && sat->pos[2] == 0.0) {
        return false;
    }

    return true;
}
