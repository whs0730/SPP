#pragma once

#include "obs.h"

void CalculateGPS(gtime_t t, eph_t* eph, satpos_t* sat);
void CalculateBDS(gtime_t t, eph_t* eph, satpos_t* sat);

// 卫星位置计算模块：根据广播星历计算位置、速度、钟差和钟速。
