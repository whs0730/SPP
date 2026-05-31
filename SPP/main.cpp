#include "decode.h"
#include "stream_decode.h"
#include "satpos.h"
#include "coordinate.h"
#include "error_correction.h"
#include "spp.h"
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
using namespace std;

static const char* DEFAULT_INPUT_FILE = "NovatelOEM20211114-01.log";
// 实时模式：SPP.exe --stream [ip] [port]。
// 不带 --stream 时，argv[1] 按 OEM4 二进制日志文件处理。
//参考坐标(-2267810.173，5009324.109，3221016.632)
static const char* DEFAULT_STREAM_IP = "8.148.22.229";
static const unsigned short DEFAULT_STREAM_PORT = 7003;

// 离线结果精度分析使用的参考坐标，与 NovatelOEM20211114-01.pos 中的 REF-ECEF 保持一致。
static const double REF_X = -2267804.5260;
static const double REF_Y = 5009342.3720;
static const double REF_Z = 3220991.8630;

void CountSolutionSats(obsd_t* obs,
    satpos_t* sats,
    int n,
    const nav_t* nav,
    const sol_t* sol,
    int& gps_count,
    int& bds_count,
    int& used_n,
    string& sat_list)
{
    gps_count = 0;
    bds_count = 0;
    used_n = 0;
    sat_list.clear();

    if (obs == nullptr || sats == nullptr || sol == nullptr)
    {
        return;
    }

    if (sol->stat != 1)
    {
        return;
    }

    double rec_xyz[3] = {
        sol->XYZ[0],
        sol->XYZ[1],
        sol->XYZ[2]
    };

    for (int i = 0; i < n; i++)
    {
        // =====================================================
        // 统一筛选：
        // 1. GPS/BDS
        // 2. 卫星位置有效
        // 3. 双频PIF有效
        // 4. SNR >= 30 dB-Hz
        // 5. 高度角 >= 10°
        // =====================================================
        // 输出用星统计也走SPP统一筛选，保证GS/BS/n和解算规则尽量一致。
        const eph_t* eph = GetSppEph(nav, obs[i].sat);
        if (!PassSppBasicCheck(obs[i], sats[i], eph, rec_xyz, true))
        {
            continue;
        }

        int prn = 0;
        int sys = satsys(obs[i].sat, &prn);

        char name[16];

        if (sys == SYS_GPS)
        {
            gps_count++;
            sprintf_s(name, "G%02d", prn);
        }
        else if (sys == SYS_CMP)
        {
            bds_count++;
            sprintf_s(name, "C%02d", prn);
        }
        else
        {
            continue;
        }

        used_n++;
        sat_list += name;
    }
}

static bool IsStreamModeArg(const char* arg)
{
    if (arg == nullptr)
    {
        return false;
    }

    return strcmp(arg, "--stream") == 0 ||
        strcmp(arg, "-s") == 0 ||
        strcmp(arg, "stream") == 0;
}

// 解析命令行参数：
// 1. 默认读取离线OEM4日志；
// 2. 遇到 --stream/-s/stream 时切换到实时TCP流；
// 3. 实时模式允许额外指定 ip 和 port。
static void ParseCommandLine(int argc,
    char* argv[],
    bool& realtime,
    string& input_file,
    string& stream_ip,
    unsigned short& stream_port)
{
    if (argc > 1 && IsStreamModeArg(argv[1]))
    {
        realtime = true;
        if (argc > 2)
        {
            stream_ip = argv[2];
        }
        if (argc > 3)
        {
            int port = atoi(argv[3]);
            if (port > 0 && port <= 65535)
            {
                stream_port = static_cast<unsigned short>(port);
            }
        }
    }
    else if (argc > 1)
    {
        input_file = argv[1];
    }
}

// spp_solution.txt 表头，字段宽度与后续数据输出保持一致，便于和 .pos 文件对照。
static void WriteSolutionHeader(ofstream& outSol)
{
    outSol << "  Wk        SOW     ECEF-X/m       ECEF-Y/m       ECEF-Z/m"
        << "    REF-ECEF-X/m    REF-ECEF-Y/m   REF-ECEF-Z/m"
        << "   EAST/m   NORTH/m  UP/m"
        << "         B/deg         L/deg             H/m"
        << "      VX/m     VY/m     VZ/m"
        << "     PDOP    SigmaP   SigmaV"
        << "  GS  BS  n"
        << endl;
}

// 实时流解算摘要：控制台和 stream.txt 复用同一格式，避免两处输出不一致。
static void OutputStreamSolutionEpoch(ostream& out,
    const obsd_t* obs0,
    const sol_t& sol,
    const ENU* enu,
    int used_n)
{
    if (obs0 == nullptr || enu == nullptr)
    {
        return;
    }

    out << fixed
        << "Wk " << obs0->time.week
        << " SOW " << setprecision(3) << obs0->time.sec
        << " XYZ "
        << setprecision(4) << sol.XYZ[0] << " "
        << sol.XYZ[1] << " "
        << sol.XYZ[2]
        << " ENU "
        << setprecision(3) << enu->E << " "
        << enu->N << " "
        << enu->U
        << " n=" << used_n
        << " PDOP=" << sol.pdop
        << endl;
}

// 输出Novatel原始观测值。
// 这里输出第一频点的伪距、载波、多普勒和SNR，主要用于检查解码是否正常。
static void OutputObservationEpoch(ofstream& outObs, const obs_t& obs_set)
{
    if (obs_set.n <= 0)
    {
        return;
    }

    const obsd_t* obs0 = &obs_set.data[0];

    for (int i = 0; i < obs_set.n; i++)
    {
        const obsd_t* obs = &obs_set.data[i];

        outObs << obs0->time.week << " "
            << fixed << setprecision(3) << obs0->time.sec << " "
            << sat2id(obs->sat) << " "
            << fixed << setprecision(3)
            << obs->P[0] << " "
            << fixed << setprecision(4)
            << obs->L[0] << " "
            << fixed << setprecision(3)
            << obs->D[0] << " "
            << fixed << setprecision(3)
            << obs->SNR[0] << endl;
    }
}

// 按卫星系统调用广播星历计算函数。
// 注意：本工程解码时已经把BDS星历toe/toc从BDT转换为GPST，
// 因此这里GPS和BDS都传入GPST时刻；BDS轨道项需要的原始BDT秒保存在eph->toes中。
static bool CalculateBroadcastSatState(gtime_t calc_time,
    eph_t* eph,
    int sat_no,
    satpos_t& sat)
{
    if (!has_valid_eph(eph))
    {
        return false;
    }

    satpos_t tmp{};
    tmp.sat = sat_no;
    tmp.time = calc_time;

    int prn = 0;
    int sys = satsys(sat_no, &prn);

    if (sys == SYS_GPS)
    {
        CalculateGPS(calc_time, eph, &tmp);
    }
    else if (sys == SYS_CMP)
    {
        CalculateBDS(calc_time, eph, &tmp);
    }
    else
    {
        return false;
    }

    if (tmp.pos[0] == 0.0 && tmp.pos[1] == 0.0 && tmp.pos[2] == 0.0)
    {
        return false;
    }

    sat = tmp;
    return true;
}

// 计算信号发射时刻的卫星状态，并把卫星坐标旋转到接收时刻对应的地固系。
// 步骤：
// 1. 用双频IF伪距 PIF/c 得到近似传播时间；
// 2. 得到近似发射时刻，并计算一次卫星钟差；
// 3. 用卫星钟差修正发射时刻，再重新计算卫星位置、速度、钟差、钟速；
// 4. 按传播时间做地球自转改正。
static bool ComputeSatStateAtTransmitTime(const obsd_t& obs,
    eph_t* eph,
    satpos_t& sat)
{
    double pif = GetPIF((obsd_t*)&obs, eph);
    if (pif <= 0.0)
    {
        return false;
    }

    double tau0 = pif / Clight;
    gtime_t tx0 = timeadd(obs.time, -tau0);

    satpos_t sat0{};
    if (!CalculateBroadcastSatState(tx0, eph, obs.sat, sat0))
    {
        return false;
    }

    // 卫星钟差为“卫星钟 - 系统时”，发射时刻需要再扣除该钟差。
    gtime_t tx = timeadd(tx0, -sat0.clk);

    if (!CalculateBroadcastSatState(tx, eph, obs.sat, sat))
    {
        return false;
    }

    double tau = timediff(obs.time, tx);
    ApplyEarthRotationCorrection(sat, tau);
    sat.time = tx;

    return true;
}

// 对当前历元每颗卫星计算位置、速度、钟差和钟速。
// validSat统计的是“有可用星历且成功算出发射时刻卫星状态”的GPS/BDS卫星数，
static int ComputeSatellitePositions(obsd_t* obs,
    int n,
    nav_t* nav,
    satpos_t* sats)
{
    int validSat = 0;

    for (int i = 0; i < n; i++)
    {
        obsd_t* ob = &obs[i];
        eph_t* eph = find_eph(nav, ob->sat);

        if (!has_valid_eph(eph))
        {
            continue;
        }

        int prn = 0;
        int sys = satsys(ob->sat, &prn);
        if (sys != SYS_GPS && sys != SYS_CMP)
        {
            continue;
        }

        satpos_t sat{};
        if (!ComputeSatStateAtTransmitTime(*ob, eph, sat))
        {
            continue;
        }

        sats[i] = sat;
        validSat++;
    }

    return validSat;
}

// 输出单点定位结果。
static void OutputSolutionEpoch(ofstream& outSol,
    obsd_t* obs,
    int n,
    nav_t* nav,
    satpos_t* sats,
    int validSat,
    const sol_t& sol,
    bool realtime,
    ostream* outStream)
{
    obsd_t* obs0 = &obs[0];

    XYZ solXYZ;
    solXYZ.X = sol.XYZ[0];
    solXYZ.Y = sol.XYZ[1];
    solXYZ.Z = sol.XYZ[2];

    BLH* solBLH = XYZtoBLH(solXYZ, 6378137.0, 1.0 / 298.257223563);

    double B_deg = solBLH->B * 180.0 / PI;
    double L_deg = solBLH->L * 180.0 / PI;
    double H_m = solBLH->H;

    XYZ refXYZ;
    refXYZ.X = REF_X;
    refXYZ.Y = REF_Y;
    refXYZ.Z = REF_Z;

    BLH* refBLH = XYZtoBLH(refXYZ, 6378137.0, 1.0 / 298.257223563);

    double dxyz[3];
    dxyz[0] = sol.XYZ[0] - REF_X;
    dxyz[1] = sol.XYZ[1] - REF_Y;
    dxyz[2] = sol.XYZ[2] - REF_Z;

    ENU* enu = xyz2enu(refBLH, dxyz);

    int gps_count = 0;
    int bds_count = 0;
    int used_n = 0;
    string sat_list;

    CountSolutionSats(obs, sats, n, nav, &sol,
        gps_count, bds_count, used_n, sat_list);

    solvel_t vsol{};
    int vel_nv = 0;

    double vx = 0.0;
    double vy = 0.0;
    double vz = 0.0;
    double sigma_v = 999.990;

    bool vel_ok = SPP_Speed(obs, n, (sol_t*)&sol, sats, &vsol, &vel_nv);

    if (vel_ok && vsol.stat == 1)
    {
        vx = vsol.V[0];
        vy = vsol.V[1];
        vz = vsol.V[2];
        sigma_v = vsol.sigma0;
    }

    outSol << fixed;

    outSol << setw(4) << obs0->time.week << " "
        << setw(10) << setprecision(3) << obs0->time.sec << "  "

        << setw(14) << setprecision(4) << sol.XYZ[0] << " "
        << setw(14) << setprecision(4) << sol.XYZ[1] << " "
        << setw(14) << setprecision(4) << sol.XYZ[2] << "  "

        << setw(14) << setprecision(4) << REF_X << " "
        << setw(14) << setprecision(4) << REF_Y << " "
        << setw(14) << setprecision(4) << REF_Z << "  "

        << setw(7) << setprecision(3) << enu->E << " "
        << setw(8) << setprecision(3) << enu->N << " "
        << setw(7) << setprecision(3) << enu->U << "  "

        << setw(13) << setprecision(8) << B_deg << " "
        << setw(14) << setprecision(8) << L_deg << " "
        << setw(14) << setprecision(3) << H_m << "  "

        << setw(7) << setprecision(3) << vx << " "
        << setw(7) << setprecision(3) << vy << " "
        << setw(7) << setprecision(3) << vz << "  "

        << setw(7) << setprecision(3) << sol.pdop << " "
        << setw(8) << setprecision(3) << sol.sigma0 << " "
        << setw(8) << setprecision(3) << sigma_v << " "

        << setw(3) << gps_count << " "
        << setw(3) << bds_count << " "
        << setw(3) << validSat << " "
        << sat_list
        << endl;

    if (realtime)
    {
        OutputStreamSolutionEpoch(cout, obs0, sol, enu, used_n);
        if (outStream != nullptr)
        {
            OutputStreamSolutionEpoch(*outStream, obs0, sol, enu, used_n);
        }
    }

    delete solBLH;
    delete refBLH;
    delete enu;
}

// 控制satpos1.txt的输出频率：只输出整秒，并从指定起始时刻起每10分钟输出一次。
static bool ShouldOutputSatposEpoch(double tow,
    double start_tow,
    int interval,
    int& last_sat_output_tow)
{
    if (fabs(tow - round(tow)) > 1e-3)
    {
        return false;
    }

    int cur_tow = (int)round(tow);

    if (cur_tow < start_tow)
    {
        return false;
    }

    if (last_sat_output_tow < 0)
    {
        last_sat_output_tow = cur_tow;
        return true;
    }

    if (cur_tow - last_sat_output_tow < interval)
    {
        return false;
    }

    last_sat_output_tow = cur_tow;
    return true;
}

// 输出接收时刻的广播星历卫星状态，用于和参考 satpos 文件对照。
// 位置按参考格式输出为 km，速度为 m/s，钟差/钟速转成微秒量级。
static void OutputSatposEpoch(ofstream& outSat,
    obsd_t* obs,
    nav_t* nav,
    int n)
{
    if (obs == nullptr || nav == nullptr || n <= 0)
    {
        return;
    }

    const obsd_t* obs0 = &obs[0];

    outSat << "* "
        << obs0->time.week << " "
        << fixed << setprecision(3) << obs0->time.sec
        << endl;

    for (int i = 0; i < n; i++)
    {
        obsd_t* ob = &obs[i];
        eph_t* eph = find_eph(nav, ob->sat);
        satpos_t sat{};

        // satpos1.txt 用于和参考卫星位置文件对照，输出接收时刻的广播星历结果；
        // SPP 解算用的发射时刻改正后卫星状态仍保存在主循环的 sats[] 中。
        if (!CalculateBroadcastSatState(ob->time, eph, ob->sat, sat))
        {
            continue;
        }

        int prn = 0;
        int sys = satsys(sat.sat, &prn);

        char satname[16];

        if (sys == SYS_GPS)
        {
            sprintf_s(satname, "BG%02d", prn);
        }
        else if (sys == SYS_CMP)
        {
            sprintf_s(satname, "BC%02d", prn);
        }
        else
        {
            continue;
        }

        outSat << satname << " "
            << fixed << setprecision(6)
            << sat.pos[0] / 1000.0 << " "
            << sat.pos[1] / 1000.0 << " "
            << sat.pos[2] / 1000.0 << " "
            << setprecision(4)
            << sat.vel[0] << " "
            << sat.vel[1] << " "
            << sat.vel[2] << " "
            << setprecision(6)
            << sat.clk * 1e6 << " "
            << sat.dclk * 1e6
            << endl;
    }
}

int main(int argc, char* argv[])
{
    // =====================================================
    // 0. 运行参数与输入源设置
    // =====================================================
    FILE* fp = nullptr;
    TcpOem4Stream stream;
    bool realtime = false;
    string input_file = DEFAULT_INPUT_FILE;
    string stream_ip = DEFAULT_STREAM_IP;
    unsigned short stream_port = DEFAULT_STREAM_PORT;

    ParseCommandLine(argc, argv, realtime, input_file, stream_ip, stream_port);

    if (realtime)
    {
        // 进入处理循环前先打开 TCP 实时流；socket 错误会保存在
        // stream.LastError() 中，便于输出具体原因。
        if (!stream.Open(stream_ip.c_str(), stream_port))
        {
            cerr << "Cannot open realtime stream " << stream_ip << ":" << stream_port
                << " (" << stream.LastError() << ")" << endl;
            return -1;
        }

        cout << "Realtime stream opened: " << stream_ip << ":" << stream_port << endl;
    }
    else if (fopen_s(&fp, input_file.c_str(), "rb") != 0 || !fp)
    {
        cerr << "Cannot open file: " << input_file << endl;
        return -1;
    }

    // =====================================================
    // 1. 输出文件初始化
    // =====================================================
    raw_t raw;
    int ret = 0;

    ofstream outObs("observations1.txt");
    ofstream outSat("satpos1.txt");
    ofstream outSol("spp_solution.txt");
    ofstream outRaw;
    ofstream outStream;

    if (realtime)
    {
        // 保存实时接收到的有效 OEM4 观测/星历报文，后续可按文件模式回放。
        outRaw.open("realtime_raw_oem.log", ios::binary);

        // 保存实时解算摘要，内容与控制台实时输出保持一致。
        outStream.open("stream.txt");
        if (!outStream.is_open())
        {
            cerr << "Cannot open stream.txt" << endl;
            return -1;
        }
    }

    WriteSolutionHeader(outSol);

    // satpos1.txt抽样输出
    double start_tow = 7 * 3600 + 25 * 60;  // 07:25:00
    const int interval = 600;               // 10分钟
    int last_sat_output_tow = -1;
    sol_t last_sol{};
    bool has_last_sol = false;

    // =====================================================
    // 2. 主循环：读帧 -> 解码 -> 按历元定位/测速/输出
    // =====================================================
    while (true)
    {
        // ret 含义来自 decode_oem4()：
        // 1 = 观测历元，2 = 星历更新，0 = 跳过消息。
        ret = realtime ? input_oem4s(&raw, &stream) : input_oem4f(&raw, fp);

        if (ret == -2)
        {
            if (realtime && !stream.LastError().empty())
            {
                cerr << "Realtime stream stopped: " << stream.LastError() << endl;
            }
            break;
        }

        if (realtime && outRaw.is_open() && (ret == 1 || ret == 2) &&
            raw.len > 0 && raw.len + 4 <= MAXRAWLEN)
        {
            // 当前 raw.len 不包含尾部 CRC，因此写入 len + 4 字节以保留完整 OEM4 帧。
            outRaw.write(reinterpret_cast<const char*>(raw.buff), raw.len + 4);
        }

        // 只有观测历元才触发定位；星历消息已在 decode_oem4() 中更新 raw.nav，
        // 然后等待下一组观测值参与解算。
        if (ret != 1) continue;
        if (raw.obs.n <= 0) continue;

        obsd_t* obs0 = &raw.obs.data[0];
        double tow = obs0->time.sec;

        // 2.1 输出当前历元观测值，便于检查解码结果。
        OutputObservationEpoch(outObs, raw.obs);

        // 2.2 根据当前历元观测和已解码星历计算卫星状态。
        satpos_t sats[MAXOBS]{};
        int validSat = ComputeSatellitePositions(raw.obs.data, raw.obs.n, &raw.nav, sats);

        // 2.3 单点定位。SPP内部按可用系统自动选择单系统4参数或双系统5参数。
        sol_t sol{};
        int spp_nv = 0;

        bool spp_ok = false;

        if (validSat >= 4)
        {
            const sol_t* init_sol = has_last_sol ? &last_sol : nullptr;
            spp_ok = SPP(raw.obs.data, raw.obs.n, &raw.nav, &sol, sats, &spp_nv, init_sol);
        }

        if (spp_ok && sol.stat == 1)
        {
            last_sol = sol;
            has_last_sol = true;

            // 2.4 输出定位、测速、ENU误差和用星统计。
            OutputSolutionEpoch(outSol, raw.obs.data, raw.obs.n, &raw.nav,
                sats, validSat, sol, realtime,
                realtime && outStream.is_open() ? &outStream : nullptr);
        }

        // 2.5 卫星位置文件只做低频抽样输出，用于和参考satpos格式对照。
        if (!ShouldOutputSatposEpoch(tow, start_tow, interval, last_sat_output_tow))
        {
            continue;
        }

        OutputSatposEpoch(outSat, raw.obs.data, &raw.nav, raw.obs.n);
    }

    if (fp != nullptr)
    {
        fclose(fp);
    }

    outObs.close();
    outSat.close();
    outSol.close();
    if (outRaw.is_open())
    {
        outRaw.close();
    }
    if (outStream.is_open())
    {
        outStream.close();
    }

    cout << (realtime ? "Realtime processing finished." : "Finished processing.") << endl;

    return 0;
}
