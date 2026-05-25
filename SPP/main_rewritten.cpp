#include "decode.h"
#include "spp.h"
#include "error_correction.h"
#include "coordinate.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>

using namespace std;

static bool IsValidSatpos(const satpos_t *sat)
{
    if (sat == nullptr) return false;
    return !(sat->pos[0] == 0.0 && sat->pos[1] == 0.0 && sat->pos[2] == 0.0);
}

static string GetSatName(int sat)
{
    int prn = 0;
    int sys = satsys(sat, &prn);

    if (sys == SYS_GPS) {
        return "G" + string(prn < 10 ? "0" : "") + to_string(prn);
    }
    if (sys == SYS_CMP) {
        return "C" + string(prn < 10 ? "0" : "") + to_string(prn);
    }
    return "";
}

static void CountUsedSats(obsd_t *obs, satpos_t *sats, int n, int &gps_count, int &bds_count)
{
    gps_count = 0;
    bds_count = 0;

    for (int i = 0; i < n; i++)
    {
        int prn = 0;
        int sys = satsys(obs[i].sat, &prn);

        if (sys != SYS_GPS && sys != SYS_CMP) continue;
        if (!IsValidSatpos(&sats[i])) continue;
        if (obs[i].P[0] == 0.0 || obs[i].P[1] == 0.0) continue;

        if (sys == SYS_GPS) gps_count++;
        else if (sys == SYS_CMP) bds_count++;
    }
}

static void OutputObservationEpoch(ofstream &outObs, obsd_t *obs, int n)
{
    for (int i = 0; i < n; i++)
    {
        string satid = sat2id(obs[i].sat);
        GPSTime gps(obs[i].time.week, obs[i].time.sec);
        CommonTime *ct = GPSTimeToCommonTime(gps);

        outObs << ct->Year << " "
               << ct->Month << " "
               << ct->Day << " "
               << ct->Hour << " "
               << ct->Minute << " "
               << fixed << setprecision(3) << ct->Second << " "
               << satid << " "
               << fixed << setprecision(3) << obs[i].P[0] << " "
               << fixed << setprecision(4) << obs[i].L[0] << " "
               << fixed << setprecision(3) << obs[i].D[0] << " "
               << fixed << setprecision(3) << obs[i].SNR[0]
               << endl;

        delete ct;
    }
}

static void OutputSatellitePositions(ofstream &outSat,
                                     CommonTime *ct0,
                                     obsd_t *obs0,
                                     obsd_t *obs,
                                     satpos_t *sats,
                                     int n)
{
    outSat << "*  "
           << ct0->Year << " "
           << ct0->Month << " "
           << ct0->Day << "  "
           << ct0->Hour << " "
           << ct0->Minute << "  "
           << fixed << setprecision(8) << ct0->Second
           << endl;

    outSat << "*  "
           << obs0->time.week << "  "
           << fixed << setprecision(8) << obs0->time.sec
           << endl;

    for (int i = 0; i < n; i++)
    {
        int prn = 0;
        int sys = satsys(obs[i].sat, &prn);
        if (sys != SYS_GPS && sys != SYS_CMP) continue;
        if (!IsValidSatpos(&sats[i])) continue;

        string satname = GetSatName(obs[i].sat);
        if (satname.empty()) continue;

        outSat << satname << "  "
               << fixed << setprecision(6)
               << sats[i].pos[0] / 1000.0 << " "
               << sats[i].pos[1] / 1000.0 << " "
               << sats[i].pos[2] / 1000.0 << " "
               << setprecision(4)
               << sats[i].vel[0] << " "
               << sats[i].vel[1] << " "
               << sats[i].vel[2] << " "
               << setprecision(6)
               << sats[i].clk * 1e6 << " "
               << sats[i].dclk * 1e6
               << endl;
    }
}

static void OutputSppDebug(ofstream &outDbg,
                           obsd_t *obs,
                           satpos_t *sats,
                           int n,
                           sol_t *sol)
{
    if (sol == nullptr || sol->stat != 1) return;

    XYZ recXYZ;
    recXYZ.X = sol->XYZ[0];
    recXYZ.Y = sol->XYZ[1];
    recXYZ.Z = sol->XYZ[2];

    BLH *recBLH = XYZtoBLH(recXYZ, 6378137.0, 1.0 / 298.257223563);

    for (int i = 0; i < n; i++)
    {
        int prn = 0;
        int sys = satsys(obs[i].sat, &prn);
        if (sys != SYS_GPS && sys != SYS_CMP) continue;
        if (!IsValidSatpos(&sats[i])) continue;

        double pif = GetPIF(&obs[i]);
        if (pif == 0.0) continue;

        string satname = GetSatName(obs[i].sat);
        if (satname.empty()) continue;

        double dx = sats[i].pos[0] - sol->XYZ[0];
        double dy = sats[i].pos[1] - sol->XYZ[1];
        double dz = sats[i].pos[2] - sol->XYZ[2];
        double rho = sqrt(dx * dx + dy * dy + dz * dz);
        if (rho < 1.0) continue;

        double los[3] = {dx / rho, dy / rho, dz / rho};
        double azel[2] = {0.0, 0.0};
        double elev_rad = satazel(recBLH, los, azel);
        double elev_deg = elev_rad * 180.0 / PI;
        if (elev_rad <= 0.0) continue;

        double trop = Hopfield_delTrop(recBLH->H, elev_rad);
        double rec_clk = (sys == SYS_GPS) ? sol->dtr[0] : sol->dtr[1];
        double resv = pif - (rho + rec_clk - CLIGHT * sats[i].clk + trop);

        outDbg << satname
               << " SOW=" << fixed << setprecision(3) << obs[i].time.sec
               << " X=" << setprecision(3) << sats[i].pos[0]
               << " Y=" << setprecision(3) << sats[i].pos[1]
               << " Z=" << setprecision(3) << sats[i].pos[2]
               << " Clk=" << scientific << uppercase << setprecision(5) << sats[i].clk
               << nouppercase << fixed
               << " Vx=" << setprecision(4) << sats[i].vel[0]
               << " Vy=" << setprecision(4) << sats[i].vel[1]
               << " Vz=" << setprecision(4) << sats[i].vel[2]
               << " Clkd=" << scientific << uppercase << setprecision(5) << sats[i].dclk
               << nouppercase << fixed
               << " PIF=" << setprecision(4) << pif
               << " Trop=" << setprecision(3) << trop
               << " E=" << setprecision(3) << elev_deg << "deg"
               << " ResV=" << setprecision(3) << resv
               << endl;
    }

    delete recBLH;
}

int main()
{
    FILE *fp = nullptr;
    if (fopen_s(&fp, "NovatelOEM20211114-01.log", "rb") != 0 || !fp)
    {
        cerr << "Cannot open file" << endl;
        return -1;
    }

    raw_t raw;
    int ret = 0;
    int obs_msg_count = 0;

    ofstream outObs("observations1.txt");
    ofstream outSat("satpos1.txt");
    ofstream outSol("spp_solution.txt");
    ofstream outDbg("spp_debug.txt");

    outSol << "week tow X Y Z B L H GPS_Clk_m BDS_Clk_m PDOP Sigma GPSSats BDSSats Sats Stat" << endl;

    // 从 07:25:00 开始，每整秒输出一次；如果要每10分钟输出一次，把 interval_sec 改成 600
    const int start_tow = 7 * 3600 + 25 * 60;
    const int interval_sec = 1;
    int last_output_week = -1;
    int last_output_tow = -1;

    while ((ret = input_oem4f(&raw, fp)) != -2)
    {
        if (ret != 1) continue;   // 只处理观测历元

        obs_msg_count++;
        if (raw.obs.n <= 0) continue;

        OutputObservationEpoch(outObs, raw.obs.data, raw.obs.n);

        obsd_t *obs0 = &raw.obs.data[0];
        double tow = obs0->time.sec;
        if (fabs(tow - round(tow)) > 1e-3) continue;

        int cur_tow = (int)round(tow);
        int cur_week = obs0->time.week;

        if (cur_tow < start_tow) continue;
        if ((cur_tow - start_tow) % interval_sec != 0) continue;
        if (last_output_week == cur_week && cur_tow <= last_output_tow) continue;

        last_output_week = cur_week;
        last_output_tow = cur_tow;

        GPSTime gps0(obs0->time.week, obs0->time.sec);
        CommonTime *ct0 = GPSTimeToCommonTime(gps0);

        satpos_t sats[MAXOBS]{};
        int valid_sat = 0;

        for (int i = 0; i < raw.obs.n; i++)
        {
            eph_t *eph = find_eph(&raw.nav, raw.obs.data[i].sat);
            if (!has_valid_eph(eph)) continue;

            satpos_t sat{};
            sat.sat = raw.obs.data[i].sat;
            sat.time = raw.obs.data[i].time;

            int prn = 0;
            int sys = satsys(raw.obs.data[i].sat, &prn);

            if (sys == SYS_GPS) {
                CalculateGPS(raw.obs.data[i].time, eph, &sat);
            }
            else if (sys == SYS_CMP) {
                CalculateBDS(raw.obs.data[i].time, eph, &sat);
            }
            else {
                continue;
            }

            if (!IsValidSatpos(&sat)) continue;

            sats[i] = sat;
            valid_sat++;
        }

        OutputSatellitePositions(outSat, ct0, obs0, raw.obs.data, sats, raw.obs.n);

        sol_t sol{};
        int spp_nv = 0;
        bool spp_ok = false;

        if (valid_sat >= 5) {
            spp_ok = SPP_GPS(raw.obs.data, raw.obs.n, &raw.nav, &sol, sats, &spp_nv);
        }

        int gps_count = 0;
        int bds_count = 0;
        CountUsedSats(raw.obs.data, sats, raw.obs.n, gps_count, bds_count);

        if (spp_ok && sol.stat == 1)
        {
            XYZ xyz;
            xyz.X = sol.XYZ[0];
            xyz.Y = sol.XYZ[1];
            xyz.Z = sol.XYZ[2];
            BLH *blh = XYZtoBLH(xyz, 6378137.0, 1.0 / 298.257223563);

            outSol << obs0->time.week << " "
                   << fixed << setprecision(3) << obs0->time.sec << " "
                   << setprecision(4)
                   << sol.XYZ[0] << " "
                   << sol.XYZ[1] << " "
                   << sol.XYZ[2] << " "
                   << setprecision(8)
                   << blh->B * 180.0 / PI << " "
                   << blh->L * 180.0 / PI << " "
                   << setprecision(3)
                   << blh->H << " "
                   << setprecision(3)
                   << sol.dtr[0] << " "
                   << sol.dtr[1] << " "
                   << setprecision(3)
                   << sol.pdop << " "
                   << sol.sigma0 << " "
                   << gps_count << " "
                   << bds_count << " "
                   << sol.ns << " "
                   << sol.stat
                   << endl;

            OutputSppDebug(outDbg, raw.obs.data, sats, raw.obs.n, &sol);
            delete blh;
        }
        else
        {
            outSol << obs0->time.week << " "
                   << fixed << setprecision(3) << obs0->time.sec << " "
                   << setprecision(4)
                   << 0.0 << " " << 0.0 << " " << 0.0 << " "
                   << setprecision(8)
                   << 0.0 << " " << 0.0 << " "
                   << setprecision(3)
                   << 0.0 << " " << 0.0 << " " << 0.0 << " "
                   << 999.900 << " " << 999.990 << " "
                   << gps_count << " " << bds_count << " " << raw.obs.n << " " << 0
                   << endl;
        }

        delete ct0;
    }

    fclose(fp);
    outObs.close();
    outSat.close();
    outSol.close();
    outDbg.close();

    cout << "Processed " << obs_msg_count << " observation messages" << endl;
    return 0;
}
