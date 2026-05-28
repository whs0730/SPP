#include "obs.h"
#include <cmath>
#include "matrix.h"
#include "define.h"
#include "error_correction.h"
// 计算双频无电离层组合伪距。
// GPS L1/L2 IF组合不额外处理TGD；
// BDS B1I/B3I的TGD在这里改到伪距观测值上，不放进卫星钟差。
double GetPIF(obsd_t* obs,const eph_t* eph=nullptr)
{
    if (obs == nullptr)
    {
        return 0.0;
    }

    int prn = 0;
    int sys = satsys(obs->sat, &prn);

    if (sys == SYS_GPS)
    {
        // GPS L1C + L2P(Y)
        if (obs->P[0] == 0.0 || obs->P[1] == 0.0)
        {
            return 0.0;
        }

        double f1 = FREQ_GPS_L1;
        double f2 = FREQ_GPS_L2;

        return (f1 * f1 * obs->P[0] - f2 * f2 * obs->P[1]) /
            (f1 * f1 - f2 * f2);
    }
    else if (sys == SYS_CMP)
    {
        // BDS B1I + B3I
        // 注意：B3I 按老师的 Freq=1 存在 P[1]
        if (obs->P[0] == 0.0 || obs->P[1] == 0.0)
        {
            return 0.0;
        }

        double f1 = FREQ_BDS_B1;
        double f3 = FREQ_BDS_B3;
        double P1 = obs->P[0];  // B1I
        double P3 = obs->P[1];  // B3I，按你的代码存在 P[1]
        // BDS B1I 要改 TGD，B3I 不改。
        // TGD 单位是秒，乘光速变成米；这里处理后，后续最小二乘直接使用IF伪距。
        if (eph != nullptr)
        {
            P1 -= Clight * eph->tgd[0];
        }

        return (f1 * f1 * P1 - f3 * f3 * P3) /
            (f1 * f1 - f3 * f3);
    }

    return 0.0;
}
// SNR 当前已按 dB-Hz 存入结构体，单独封装便于以后兼容比例存储。
static double GetSnrDbHz(double snr)
{
    if (snr == 0)
    {
        return 0.0;
    }

    return snr;
}
// GPS/BDS 双频信噪比质量控制，低于 30 dB-Hz 的观测不参与解算。
static bool PassSnrCheck(const obsd_t& obs)
{
    int prn = 0;
    int sys = satsys(obs.sat, &prn);

    double snr1 = 0.0;
    double snr2 = 0.0;

    if (sys == SYS_GPS)
    {
        // GPS L1C + L2P(Y)
        snr1 = GetSnrDbHz(obs.SNR[0]);
        snr2 = GetSnrDbHz(obs.SNR[1]);
    }
    else if (sys == SYS_CMP)
    {
        // BDS B1I + B3I，现在你已经按老师方式存在 P[0]/P[1]
        snr1 = GetSnrDbHz(obs.SNR[0]);
        snr2 = GetSnrDbHz(obs.SNR[1]);
    }
    else
    {
        return false;
    }

    // 有些数据 SNR 可能没存，等于0时先不筛
    if (snr1 > 0.0 && snr1 < 30.0)
    {
        return false;
    }

    if (snr2 > 0.0 && snr2 < 30.0)
    {
        return false;
    }

    return true;
}
// 根据当前接收机近似坐标进行高度角筛选。
static bool PassElevationCheck(const satpos_t& sat,
    const double rec_xyz[3],
    double mask_deg = 10.0)
{
    if (rec_xyz == nullptr)
    {
        return true;
    }

    if (sat.pos[0] == 0.0 && sat.pos[1] == 0.0 && sat.pos[2] == 0.0)
    {
        return false;
    }

    XYZ recXYZ;
    recXYZ.X = rec_xyz[0];
    recXYZ.Y = rec_xyz[1];
    recXYZ.Z = rec_xyz[2];

    BLH* recBLH = XYZtoBLH(recXYZ, 6378137.0, 1.0 / 298.257223563);

    double dx = sat.pos[0] - rec_xyz[0];
    double dy = sat.pos[1] - rec_xyz[1];
    double dz = sat.pos[2] - rec_xyz[2];

    double rho = sqrt(dx * dx + dy * dy + dz * dz);

    if (rho < 1.0)
    {
        delete recBLH;
        return false;
    }

    double los[3];
    los[0] = dx / rho;
    los[1] = dy / rho;
    los[2] = dz / rho;

    double azel[2] = { 0.0, 0.0 };
    double elev_rad = satazel(recBLH, los, azel);

    delete recBLH;

    return elev_rad >= mask_deg * PI / 180.0;
}

static bool IsValidSatposSPP(const satpos_t& sat)
{
    return !(sat.pos[0] == 0.0 &&
        sat.pos[1] == 0.0 &&
        sat.pos[2] == 0.0);
}

// SPP只使用已经成功解码的广播星历。
// 这里集中判断星历是否存在、卫星号是否有效、轨道长半轴是否有效。
static bool HasValidEphSPP(const eph_t* eph)
{
    if (eph == nullptr)
    {
        return false;
    }

    if (eph->sat <= 0)
    {
        return false;
    }

    if (eph->sqrtA <= 0.0)
    {
        return false;
    }

    return true;
}

// 根据内部卫星号从导航数据中取星历。
// 后续PIF、最小二乘和输出统计都通过这个函数取星历，避免各处判断标准不一致。
const eph_t* GetSppEph(const nav_t* nav, int sat)
{
    if (nav == nullptr || sat < 1 || sat > MAXSAT)
    {
        return nullptr;
    }

    const eph_t* eph = &nav->eph[sat - 1];
    return HasValidEphSPP(eph) ? eph : nullptr;
}

// 计算当前接收机位置下的卫星高度角，单位为弧度。
// 高度角筛选、对流层改正和最终用星统计都复用这个函数。
static bool ComputeElevationRad(const satpos_t& sat,
    const double rec_xyz[3],
    double& elev_rad)
{
    elev_rad = 0.0;

    if (rec_xyz == nullptr)
    {
        return false;
    }

    if (!IsValidSatposSPP(sat))
    {
        return false;
    }

    XYZ recXYZ;
    recXYZ.X = rec_xyz[0];
    recXYZ.Y = rec_xyz[1];
    recXYZ.Z = rec_xyz[2];

    BLH* recBLH = XYZtoBLH(recXYZ, 6378137.0, 1.0 / 298.257223563);

    double dx = sat.pos[0] - rec_xyz[0];
    double dy = sat.pos[1] - rec_xyz[1];
    double dz = sat.pos[2] - rec_xyz[2];

    double rho = sqrt(dx * dx + dy * dy + dz * dz);

    if (rho < 1.0)
    {
        delete recBLH;
        return false;
    }

    double los[3];
    los[0] = dx / rho;
    los[1] = dy / rho;
    los[2] = dz / rho;

    double azel[2] = { 0.0, 0.0 };
    elev_rad = satazel(recBLH, los, azel);

    delete recBLH;
    return true;
}


// =====================================================
// SPP统一质量控制：
// 1. GPS/BDS
// 2. 卫星位置有效
// 3. 双频PIF有效
// 4. SNR >= 30 dB-Hz
// 5.高度角 >= 10°
// =====================================================
bool PassSppBasicCheck(const obsd_t& obs,
    const satpos_t& sat,
    const eph_t* eph,
    const double* rec_xyz,
    bool check_elev,
    double* pif_out = nullptr)
{
    int prn = 0;
    int sys = satsys(obs.sat, &prn);

    // 当前SPP只解GPS和BDS；其他系统先不参与定位。
    if (sys != SYS_GPS && sys != SYS_CMP)
    {
        return false;
    }

    // 没有卫星坐标时无法建立观测方程。
    if (!IsValidSatposSPP(sat))
    {
        return false;
    }

    // 统一在这里计算双频无电离层组合。
    // BDS的TGD改正在GetPIF()内部完成，因此预筛和最小二乘使用同一个PIF值。
    double pif = GetPIF((obsd_t*)&obs, eph);
    if (pif == 0.0)
    {
        return false;
    }

    // SNR门限统一放在这里，避免输出统计和真正解算使用不同卫星。
    if (!PassSnrCheck(obs))
    {
        return false;
    }

    // 初始迭代时还没有可靠接收机坐标，所以高度角检查可选。
    // 一旦有接收机坐标，就启用10度截止高度角。
    if (check_elev && rec_xyz != nullptr)
    {
        double elev_rad = 0.0;
        if (!ComputeElevationRad(sat, rec_xyz, elev_rad) ||
            elev_rad < 10.0 * PI / 180.0)
        {
            return false;
        }
    }

    // 调用方如果需要PIF数值，可以直接取走，避免重复计算造成不一致。
    if (pif_out != nullptr)
    {
        *pif_out = pif;
    }

    return true;
}
// 构建 SPP 线性化观测方程并求一次最小二乘改正。
Matrix* LeastSquares(Matrix X, obsd_t* obs, satpos_t* sat, const nav_t* nav, int n, int nv)
{
    Matrix B = zeros(nv, 5);
    Matrix w = zeros(nv, 1);
    int index = 0;
    for (int i = 0; i < n; i++)
    {
        satpos_t s = sat[i];
        // 统一获取星历，并交给PassSppBasicCheck()完成系统、伪距、SNR等预筛。
        const eph_t* eph = GetSppEph(nav, obs[i].sat);

        double pif_check = 0.0;
        if (!PassSppBasicCheck(obs[i], s, eph, nullptr, false, &pif_check))
        {
            continue;
        }

        if (index >= nv)
        {
            return nullptr;
        }
        // sat[]在主流程ComputeSatStateAtTransmitTime()中已经按发射时刻计算，
        // 并完成地球自转改正；这里直接用该状态建立伪距观测方程，避免重复改正。
        double Xs = s.pos[0], Ys = s.pos[1], Zs = s.pos[2];
        double delX = Xs - X[0][0], delY = Ys - X[1][0], delZ = Zs - X[2][0];
        double p = sqrt(delX * delX + delY * delY + delZ * delZ);
        double Trop = 0.0;

        if (!(X[0][0] == 0.0 && X[1][0] == 0.0 && X[2][0] == 0.0))
        {
            // 从第二轮迭代开始已有接收机近似坐标，可以计算高度角并进行10度截止。
            double rec_xyz[3] = {
                X[0][0],
                X[1][0],
                X[2][0]
            };

            double elev_rad = 0.0;
            if (!ComputeElevationRad(s, rec_xyz, elev_rad) ||
                elev_rad < 10.0 * PI / 180.0)
            {
                continue;
            }

            XYZ recXYZ;
            recXYZ.X = X[0][0];
            recXYZ.Y = X[1][0];
            recXYZ.Z = X[2][0];

            BLH* recBLH = XYZtoBLH(recXYZ, 6378137.0, 1.0 / 298.257223563);

            // 对流层改正使用同一个高度角，避免筛选和改正使用不同几何量。
            Trop = Hopfield_delTrop(recBLH->H, elev_rad);

            delete recBLH;
        }
        if (p < 1.0)
        {
            return nullptr;
        }
        double B1 = -delX / p, B2 = -delY / p, B3 = -delZ / p;
        double rec_clk = 0.0;
        int prn = 0;
        int sys = satsys(obs[i].sat, &prn);
        if (sys == SYS_GPS)
        {
            B[index][0] = B1;
            B[index][1] = B2;
            B[index][2] = B3;
            B[index][3] = 1.0;
            B[index][4] = 0.0;
            rec_clk = X[3][0];
        }
        else if (sys == SYS_CMP)
        {
            B[index][0] = B1;
            B[index][1] = B2;
            B[index][2] = B3;
            B[index][3] = 0.0;
            B[index][4] = 1.0;
            rec_clk = X[4][0];
        }
        w[index][0] = pif_check - (p + rec_clk - Clight * s.clk + Trop);
        index++;
    }
    if (index < 5)
    {
        return nullptr;
    }
    // 只保留真正参与解算的行
    B.resize(index);
    w.resize(index);
    Matrix N, dx;
    try
    {
        N = mul(transpose(B), B);
        N = inverse(N);
        dx = mul(N, transpose(B));
        dx = mul(dx, w);
    }
    catch (...)
    {
        return nullptr;
    }
    Matrix *m = new Matrix[4];
    m[0] = add(dx, X);
    m[1] = B;
    m[2] = w;
    m[3] = dx;
    return m;
}
// 迭代执行最小二乘，直到接收机位置改正量足够小或达到最大迭代次数。
Matrix* IterativeSolution(obsd_t* obs, satpos_t* sat, const nav_t* nav, double* PIF, int n, int nv, int times = 1, double Threshold = 1e-6)
{
    Matrix X = zeros(5, 1);
    Matrix *m = nullptr; // X,B,w
    while (times <= 10)
    {
        if (m != nullptr)
        {
            delete[] m;
            m = nullptr;
        }
        m = LeastSquares(X, obs, sat, nav,n, nv);
        if (m == nullptr)
        {
            return nullptr;
        }
        Matrix dX = sub(m[0], X);
        double dpos = sqrt(dX[0][0] * dX[0][0] + dX[1][0] * dX[1][0] + dX[2][0] * dX[2][0]);
        if (dpos < Threshold)
        {
            return m;
        }

        X = m[0];
        times++;
    }
    return m;
}
// 预计算可用卫星的 IF 组合伪距，并统计 GPS/BDS 数量。
void PIF(obsd_t* obs,
    int n,
    satpos_t* sat,
    const nav_t* nav,
    int* nv,
    double* pif,
    int* a = nullptr,
    int* b = nullptr)
{
    int index = 0;

    if (nv) *nv = 0;
    if (a)  *a = 0;
    if (b)  *b = 0;

    for (int i = 0; i < n; i++)
    {
        obsd_t o = obs[i];
        satpos_t s = sat[i];
        // PIF预筛也传入星历，使BDS TGD改正和后续最小二乘保持一致。
        const eph_t* eph = GetSppEph(nav, obs[i].sat);
        double pif_value = 0.0;

        // PIF阶段先不做高度角，因为此时还没有可靠接收机坐标。
        // 高度角筛选会在最小二乘迭代中有近似坐标后再执行。
        if (!PassSppBasicCheck(o, s, eph, nullptr, false, &pif_value))
        {
            continue;
        }

        int prn = 0;
        int sys = satsys(o.sat, &prn);

        if (sys == SYS_GPS)
        {
            if (a) (*a)++;
        }
        else if (sys == SYS_CMP)
        {
            if (b) (*b)++;
        }

        pif[index++] = pif_value;
    }

    if (nv) *nv = index;
}
// GPS/BDS 双系统 SPP 主函数：先筛选观测，再迭代估计 x/y/z/dtG/dtC。
bool SPP(obsd_t *obs, int n, const nav_t *nav, sol_t *sol, satpos_t *sat, int *nv)
{
    if (obs == nullptr || sol == nullptr || sat == nullptr || nv == nullptr)
    {
        return false;
    }
    double X0 = 0.0, Y0 = 0.0, Z0 = 0.0; // 位置初始化
    int times = 0;
    double *pif = new double[n];
    int a = 0;
    int b = 0;
    // 计算双频组合观测值
    PIF(obs, n, sat, nav, nv, pif, &a, &b);
    if (*nv < 5 || a <= 0 || b <= 0)
    {
        sol->stat = 0;
        delete[] pif;
        return false;
    }
    // 迭代解算
    Matrix *m = IterativeSolution(obs, sat, nav,pif, n, *nv);
    if (m == nullptr)
    {
        sol->stat = 0;
        delete[] pif;
        return false;
    }
    Matrix X = m[0]; // 坐标和接收机钟差
    Matrix B = m[1];
    Matrix w = m[2];
    Matrix dx = m[3];  
    *nv = (int)B.size();// 最后一次迭代真实值与参考值的差
    Matrix v = sub(mul(B, dx), w); // 残差
    double vtv = mul(transpose(v), v)[0][0];
    double sigma0 = 0.0;
    if (*nv > 5)
    {
        sigma0 = sqrt(vtv / double(*nv - 5));
    }
    Matrix Qxx = inverse(mul(transpose(B), B));
    Matrix Dxx = scalar_mul(sigma0*sigma0, Qxx);
    double PDOP = sqrt(Qxx[0][0] + Qxx[1][1] + Qxx[2][2]);
    sol->time = obs[0].time;
    sol->XYZ[0] = X[0][0];
    sol->XYZ[1] = X[1][0];
    sol->XYZ[2] = X[2][0];
    sol->dtr[0] = X[3][0]; // GPS
    sol->dtr[1] = X[4][0]; // BDS
    sol->Q = Qxx;
    sol->pdop = PDOP;
    sol->sigma0 = sigma0;
    sol->ns = *nv;
    sol->stat = 1;
    delete[] pif;
    delete[] m;
    return true;
}
// 多普勒转距离率函数。
double DopplerToRangeRate(obsd_t o)
{
    int prn = 0;
    int sys = satsys(o.sat, &prn);

    double f = 0.0;

    if (sys == SYS_GPS)
    {
        f = FREQ_GPS_L1;
    }
    else if (sys == SYS_CMP)
    {
        f = FREQ_BDS_B1;
    }
    else
    {
        return 0.0;
    }

    double lambda = Clight / f;
    // 卫星接近时 Doppler 为正，而距离率为负
    double D_mps = -lambda * o.D[0];

    return D_mps;
}
// 统计可参与测速的观测数量。
int CountSpeedObs(obsd_t *obs, satpos_t *sat, int n, sol_t *sol)
{
    if (obs == nullptr || sat == nullptr || sol == nullptr)
    {
        return 0;
    }

    if (sol->stat != 1)
    {
        return 0;
    }

    int nv = 0;

    for (int i = 0; i < n; i++)
    {
        obsd_t o = obs[i];
        satpos_t s = sat[i];

        int prn = 0;
        int sys = satsys(o.sat, &prn);

        if (sys != SYS_GPS && sys != SYS_CMP)
        {
            continue;
        }

        // 多普勒不能为0
        if (fabs(o.D[0]) < 1e-9)
        {
            continue;
        }

        // 卫星位置不能为0
        if (s.pos[0] == 0.0 && s.pos[1] == 0.0 && s.pos[2] == 0.0)
        {
            continue;
        }

        // 卫星速度不能为0
        if (s.vel[0] == 0.0 && s.vel[1] == 0.0 && s.vel[2] == 0.0)
        {
            continue;
        }

        nv++;
    }

    return nv;
}
// 速度最小二乘求解，未知数为 Vx、Vy、Vz 和接收机钟速。
Matrix *SpeedLeastSquares(obsd_t *obs, satpos_t *sat, int n, sol_t *sol, int nv)
{
    Matrix B = zeros(nv, 4);
    Matrix w = zeros(nv, 1);

    int index = 0;

    for (int i = 0; i < n; i++)
    {
        obsd_t o = obs[i];
        satpos_t s = sat[i];

        int prn = 0;
        int sys = satsys(o.sat, &prn);

        if (sys != SYS_GPS && sys != SYS_CMP)
        {
            continue;
        }

        if (fabs(o.D[0]) < 1e-9)
        {
            continue;
        }

        if (s.pos[0] == 0.0 && s.pos[1] == 0.0 && s.pos[2] == 0.0)
        {
            continue;
        }

        if (s.vel[0] == 0.0 && s.vel[1] == 0.0 && s.vel[2] == 0.0)
        {
            continue;
        }

        double Xs = s.pos[0];
        double Ys = s.pos[1];
        double Zs = s.pos[2];

        double Xr = sol->XYZ[0];
        double Yr = sol->XYZ[1];
        double Zr = sol->XYZ[2];

        double delX = Xs - Xr;
        double delY = Ys - Yr;
        double delZ = Zs - Zr;

        double rho = sqrt(delX * delX + delY * delY + delZ * delZ);

        if (rho < 1.0)
        {
            continue;
        }

        // 方向余弦
        double l = delX / rho;
        double m = delY / rho;
        double nn = delZ / rho;

        // 多普勒转距离率，单位 m/s
        double D_mps = DopplerToRangeRate(o);

        // 卫星速度沿视线方向投影
        double sat_rate = l * s.vel[0] + m * s.vel[1] + nn * s.vel[2];
        if (index >= nv)
        {
            return nullptr;
        }
        // B矩阵
        // D = sat_rate - l*Vx - m*Vy - n*Vz + dtrd - c*dclk
        // 整理：
        // D - sat_rate + c*dclk = -l*Vx -m*Vy -n*Vz + dtrd
        B[index][0] = -l;
        B[index][1] = -m;
        B[index][2] = -nn;
        B[index][3] = 1.0;

        // w矩阵
        w[index][0] = D_mps - sat_rate + Clight * s.dclk;

        index++;
    }
    if (index != nv)
    {
        return nullptr;
    }
    Matrix N, x;

    N = mul(transpose(B), B);
    N = inverse(N);

    x = mul(N, transpose(B));
    x = mul(x, w);

    Matrix v = sub(mul(B, x), w);

    Matrix *result = new Matrix[4];

    result[0] = x; // Vx Vy Vz dtrd
    result[1] = B;
    result[2] = w;
    result[3] = v; // 残差

    return result;
}
// 单点测速入口：定位成功后，用多普勒观测估计速度。
bool SPP_Speed(obsd_t *obs, int n, sol_t *sol, satpos_t *sat, solvel_t *vsol, int *nv)
{
    if (obs == nullptr || sat == nullptr || sol == nullptr || vsol == nullptr || nv == nullptr)
    {
        return false;
    }

    if (sol->stat != 1)
    {
        vsol->stat = 0;
        *nv = 0;
        return false;
    }

    *nv = CountSpeedObs(obs, sat, n, sol);

    // 测速至少4颗卫星
    if (*nv < 4)
    {
        vsol->stat = 0;
        return false;
    }

    Matrix *m = nullptr;

    try
    {
        m = SpeedLeastSquares(obs, sat, n, sol, *nv);
    }
    catch (...)
    {
        vsol->stat = 0;
        if (m != nullptr)
        {
            delete[] m;
        }
        return false;
    }

    if (m == nullptr)
    {
        vsol->stat = 0;
        return false;
    }

    if (m == nullptr)
    {
        vsol->stat = 0;
        return false;
    }
    Matrix x = m[0];
    Matrix B = m[1];
    Matrix w = m[2];
    Matrix v = m[3];

    double vtv = mul(transpose(v), v)[0][0];

    double sigma0 = 0.0;
    if (*nv > 4)
    {
        sigma0 = sqrt(vtv / double(*nv - 4));
    }

    Matrix Qxx = inverse(mul(transpose(B), B));

    vsol->time = sol->time;

    vsol->V[0] = x[0][0];
    vsol->V[1] = x[1][0];
    vsol->V[2] = x[2][0];

    // 接收机钟速，单位 m/s
    vsol->dtrd = x[3][0];

    vsol->Q = Qxx;
    vsol->sigma0 = sigma0;
    vsol->ns = *nv;
    vsol->stat = 1;

    delete[] m;

    return true;
}
