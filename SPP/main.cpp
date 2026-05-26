#include "decode.h"
#include "stream_decode.h"
#include "calculate_location_spped.h"
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
static const char* DEFAULT_STREAM_IP = "8.148.22.229";
static const unsigned short DEFAULT_STREAM_PORT = 7003;

static const double REF_X = -2267809.273;
static const double REF_Y = 5009323.033;
static const double REF_Z = 3221015.978;

void CountSolutionSats(obsd_t* obs,
    satpos_t* sats,
    int n,
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
        if (!PassSppBasicCheck(obs[i], sats[i], rec_xyz, true))
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

int main(int argc, char* argv[])
{
    FILE* fp = nullptr;
    TcpOem4Stream stream;
    bool realtime = false;
    string input_file = DEFAULT_INPUT_FILE;
    string stream_ip = DEFAULT_STREAM_IP;
    unsigned short stream_port = DEFAULT_STREAM_PORT;

    // 实时模式和文件模式只是在数据来源上不同；读到完整 OEM4 帧之后，
    // 两者共用 decode_oem4() 和后续 SPP 解算流程。
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

    raw_t raw;
    int ret = 0;

    ofstream outObs("observations1.txt");
    ofstream outSat("satpos1.txt");
    ofstream outSol("spp_solution.txt");
    ofstream outRaw;

    if (realtime)
    {
        // 保存实时接收到的有效 OEM4 观测/星历报文，后续可按文件模式回放。
        outRaw.open("realtime_raw_oem.log", ios::binary);
    }

    outSol << "  Wk        SOW"
        << "     ECEF-X/m       ECEF-Y/m       ECEF-Z/m"
        << "    REF-ECEF-X/m    REF-ECEF-Y/m   REF-ECEF-Z/m"
        << "   EAST/m   NORTH/m  UP/m"
        << "         B/deg         L/deg             H/m"
        << "      VX/m     VY/m     VZ/m"
        << "     PDOP    SigmaP   SigmaV"
        << "  GS  BS  n"
        << endl;
    // ===== 时间控制 =====
    double start_tow = 7 * 3600 + 25 * 60;  // 07:25:00
    const int interval = 600;               // 10分钟
    int last_sat_output_tow = -1;

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

        // =====================================================
        // 1. 输出观测值（全部历元）
        // =====================================================
        for (int i = 0; i < raw.obs.n; i++)
        {
            obsd_t* obs = &raw.obs.data[i];

            outObs << obs0->time.week << " "
                << fixed << setprecision(3) << obs0->time.sec << " "
                << sat2id(obs->sat) << " "
                << obs->P[0] << " "
                << obs->L[0] << " "
                << obs->D[0] << endl;
        }

        // =====================================================
        // 2. 计算卫星位置
        // =====================================================
        satpos_t sats[MAXOBS]{};
        int validSat = 0;

        for (int i = 0; i < raw.obs.n; i++)
        {
            obsd_t* obs = &raw.obs.data[i];
            eph_t* eph = find_eph(&raw.nav, obs->sat);

            // 实时流中观测值可能先于对应星历到达；此时暂时跳过该卫星。
            if (!has_valid_eph(eph)) continue;

            satpos_t sat{};
            sat.sat = obs->sat;
            sat.time = obs->time;

            int prn = 0;
            int sys = satsys(obs->sat, &prn);

            if (sys == SYS_GPS)
                CalculateGPS(obs->time, eph, &sat);
            else if (sys == SYS_CMP)
                CalculateBDS(obs->time, eph, &sat);
            else
                continue;

            if (sat.pos[0] == 0.0 && sat.pos[1] == 0.0 && sat.pos[2] == 0.0)
                continue;

            sats[i] = sat;
            validSat++;
        }

        // =====================================================
        // 3. SPP
        // =====================================================
        sol_t sol{};
        int spp_nv = 0;

        bool spp_ok = false;

        if (validSat >= 5)
        {
            spp_ok = SPP(raw.obs.data, raw.obs.n, &raw.nav, &sol, sats, &spp_nv);
        }

        if (spp_ok && sol.stat == 1)
        {
            // =====================================================
// 输出 SPP solution：对齐 .pos 格式
// =====================================================

            // 解算坐标转 BLH
            XYZ solXYZ;
            solXYZ.X = sol.XYZ[0];
            solXYZ.Y = sol.XYZ[1];
            solXYZ.Z = sol.XYZ[2];

            BLH* solBLH = XYZtoBLH(solXYZ, 6378137.0, 1.0 / 298.257223563);

            double B_deg = solBLH->B * 180.0 / PI;
            double L_deg = solBLH->L * 180.0 / PI;
            double H_m = solBLH->H;

            // 参考坐标转 BLH，用于 ENU
            XYZ refXYZ;
            refXYZ.X = REF_X;
            refXYZ.Y = REF_Y;
            refXYZ.Z = REF_Z;

            BLH* refBLH = XYZtoBLH(refXYZ, 6378137.0, 1.0 / 298.257223563);

            // XYZ误差
            double dxyz[3];
            dxyz[0] = sol.XYZ[0] - REF_X;
            dxyz[1] = sol.XYZ[1] - REF_Y;
            dxyz[2] = sol.XYZ[2] - REF_Z;

            // 转 ENU
            ENU* enu = xyz2enu(refBLH, dxyz);

            int gps_count = 0;
            int bds_count = 0;
            int used_n = 0;
            string sat_list;

            CountSolutionSats(raw.obs.data, sats, raw.obs.n, &sol,
                gps_count, bds_count, used_n, sat_list);

            // =====================================================
 // 单点测速
 // =====================================================
            solvel_t vsol{};
            int vel_nv = 0;

            double vx = 0.0;
            double vy = 0.0;
            double vz = 0.0;
            double sigma_v = 999.990;

            bool vel_ok = SPP_Speed(raw.obs.data, raw.obs.n, &sol, sats, &vsol, &vel_nv);

            if (vel_ok && vsol.stat == 1)
            {
                vx = vsol.V[0];
                vy = vsol.V[1];
                vz = vsol.V[2];
                sigma_v = vsol.sigma0;
            }

            // 如果你已经实现了 SPP_Speed，可以把上面四个变量换成测速结果

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
                << setw(3) << used_n << " "
                << sat_list
                << endl;

            if (realtime)
            {
                cout << fixed
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

            delete solBLH;
            delete refBLH;
            delete enu;
        }

        // =====================================================
        // 4. satpos筛选（只影响卫星输出）
        // =====================================================

        // 只保留整秒
        if (fabs(tow - round(tow)) > 1e-3)
            continue;

        int cur_tow = (int)round(tow);

        if (cur_tow < start_tow)
            continue;

        if (last_sat_output_tow < 0)
        {
            last_sat_output_tow = cur_tow;
        }
        else if (cur_tow - last_sat_output_tow < interval)
        {
            continue;
        }

        last_sat_output_tow = cur_tow;

        // =====================================================
        // 5. 输出satpos（对齐参考格式）
        // =====================================================
        outSat << "* "
            << obs0->time.week << " "
            << fixed << setprecision(3) << obs0->time.sec
            << endl;

        for (int i = 0; i < raw.obs.n; i++)
        {
            satpos_t* sat = &sats[i];

            if (sat->pos[0] == 0.0 && sat->pos[1] == 0.0 && sat->pos[2] == 0.0)
                continue;

            // ===== BG / BC =====
            int prn = 0;
            int sys = satsys(sat->sat, &prn);

            char satname[16];

            if (sys == SYS_GPS)
                sprintf_s(satname, "BG%02d", prn);
            else if (sys == SYS_CMP)
                sprintf_s(satname, "BC%02d", prn);
            else
                continue;

            outSat << satname << " "
                << fixed << setprecision(6)
                << sat->pos[0] / 1000.0 << " "
                << sat->pos[1] / 1000.0 << " "
                << sat->pos[2] / 1000.0 << " "
                << setprecision(4)
                << sat->vel[0] << " "
                << sat->vel[1] << " "
                << sat->vel[2] << " "
                << setprecision(6)
                << sat->clk * 1e6 << " "
                << sat->dclk * 1e6
                << endl;
        }
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

    cout << (realtime ? "Realtime processing finished." : "Finished processing.") << endl;

    return 0;
}
