#pragma once

#include "obs.h"

double GetPIF(obsd_t* obs, const eph_t* eph = nullptr);

const eph_t* GetSppEph(const nav_t* nav, int sat);

bool PassSppBasicCheck(const obsd_t& obs,
    const satpos_t& sat,
    const eph_t* eph,
    const double* rec_xyz,
    bool check_elev,
    double* pif_out = nullptr);

bool PassSppBasicCheck(const obsd_t& obs,
    const satpos_t& sat,
    const double* rec_xyz,
    bool check_elev);

bool SPP(obsd_t* obs,
    int n,
    const nav_t* nav,
    sol_t* sol,
    satpos_t* sat,
    int* nv);

bool SPP_Speed(obsd_t* obs,
    int n,
    sol_t* sol,
    satpos_t* sat,
    solvel_t* vsol,
    int* nv);
