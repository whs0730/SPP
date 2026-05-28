#pragma once

#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <string>

#include "obs.h"

// NovAtel OEM4 原始数据状态：
// 保存当前报文缓冲区、已解码观测值、导航星历和跨历元辅助状态。
class raw_t {
public:
    gtime_t time;
    gtime_t tobs[MAXSAT][NFREQ + NEXOBS];
    double lockt[MAXSAT][NFREQ + NEXOBS];
    uint8_t halfc[MAXSAT][NFREQ + NEXOBS];
    obs_t obs;
    nav_t nav;
    uint8_t buff[MAXRAWLEN];
    int nbyte;
    int len;
    int ephsat;
    int ephset;
    char opt[256];

    raw_t() {
        nbyte = 0;
        len = 0;
        ephsat = 0;
        ephset = 0;
        memset(buff, 0, sizeof(buff));
        memset(tobs, 0, sizeof(tobs));
        memset(lockt, 0, sizeof(lockt));
        memset(halfc, 0, sizeof(halfc));
        memset(&obs, 0, sizeof(obs_t));
        memset(&nav, 0, sizeof(nav_t));
        opt[0] = '\0';
    }

    // 当读取到的新时间与当前已缓存观测的时间不同时，重置观测集。
    void reset_obs() {
        if (obs.n > 0 && fabs(timediff(obs.data[0].time, time)) > 1E-9) {
            obs.n = 0;
        }
    }
};

std::string sat2id(int sat);
eph_t* find_eph(nav_t* nav, int sat);
bool has_valid_eph(const eph_t* eph);
uint16_t U2(uint8_t* p);
int decode_oem4(raw_t* raw);
int sync_oem4(uint8_t* buff, uint8_t data);

int input_oem4f(raw_t* raw, FILE* fp);
