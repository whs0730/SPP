#pragma once

// ====================== NovAtel OEM 报文与观测常量 ======================
#define OEM4SYNC1       0xAA
#define OEM4SYNC2       0x44
#define OEM4SYNC3       0x12
#define OEM3SYNC1       0xAA
#define OEM3SYNC2       0x44
#define OEM3SYNC3       0x11
#define OEM4HLEN        28
#define OEM3HLEN        12
#define MAXRAWLEN       8192
#define MAXOBS          512
#define MAXSAT          220
#define NFREQ           3
#define NEXOBS          3
#define CODE_NONE       0
#define LLI_SLIP        1
#define LLI_HALFC       2
#define LLI_HALFA       4
#define SNR_UNIT        1.0
#define WL1             0.1902936727984
#define WL2             0.2442102134246
#define Clight          299792458.0
#define MAXVAL          8388608.0
#define OFF_FRQNO       -7

// ====================== 卫星系统标识 ======================
#define SYS_NONE    0x00
#define SYS_GPS     0x01
#define SYS_SBS     0x02
#define SYS_GLO     0x04
#define SYS_GAL     0x08
#define SYS_QZS     0x10
#define SYS_CMP     0x20
#define SYS_IRN     0x40
#define SYS_LEO     0x80
#define SYS_ALL     0xFF

// ====================== 观测码类型 ======================
#define CODE_L1C        1
#define CODE_L1P        2
#define CODE_L1W        3
#define CODE_L1Y        4
#define CODE_L1M        5
#define CODE_L1N        6
#define CODE_L1S        7
#define CODE_L1L        8
#define CODE_L1X        9
#define CODE_L2I        40  // BDS B1I 信号
#define CODE_L2C        10
#define CODE_L2D        11
#define CODE_L2S        12
#define CODE_L2L        13
#define CODE_L2X        14
#define CODE_L2P        15
#define CODE_L2W        16
#define CODE_L2Y        17
#define CODE_L2M        18
#define CODE_L2N        19
#define CODE_L5I        20
#define CODE_L5Q        21
#define CODE_L5X        22
#define CODE_L7I        23
#define CODE_L7Q        24
#define CODE_L7X        25
#define CODE_L6A        26
#define CODE_L6B        27
#define CODE_L6C        28
#define CODE_L6X        29
#define CODE_L6Z        30
#define CODE_L6S        31
#define CODE_L6L        32
#define CODE_L6I        60
#define CODE_L8I        33
#define CODE_L8Q        34
#define CODE_L8X        35

// ====================== NovAtel 消息 ID ======================
#define ID_RANGECMP     140
#define ID_RANGE        43
#define ID_RAWEPHEM     41
#define ID_IONUTC       8
#define ID_GLOEPHEMERIS 723
#define ID_GALEPHEMERIS 1122
#define ID_GPSEPHEM     7
#define ID_BDSEPHEMERIS 1696
#define ID_QZSSRAWEPHEM 1331
#define ID_NAVICEPHEMERIS 2123
#define ID_RGEB         15
#define ID_RGED         65
#define ID_REPB         14
#define ID_FRMB         54
#define ID_IONB         16
#define ID_UTCB         17

// ====================== 各系统 PRN 范围与内部编号偏移 ======================
#define MINPRNGPS   1
#define MAXPRNGPS   32
#define NSATGPS     (MAXPRNGPS-MINPRNGPS+1)
#define MINPRNGLO   1
#define MAXPRNGLO   27
#define NSATGLO     (MAXPRNGLO-MINPRNGLO+1)
#define MINPRNGAL   1
#define MAXPRNGAL   36
#define NSATGAL     (MAXPRNGAL-MINPRNGAL+1)
#define MINPRNQZS   193
#define MAXPRNQZS   202
#define MINPRNQZS_S 183
#define MAXPRNQZS_S 191
#define NSATQZS     (MAXPRNQZS-MINPRNQZS+1)
#define MINPRNCMP   1
#define MAXPRNCMP   63
#define NSATCMP     (MAXPRNCMP-MINPRNCMP+1)
#define MINPRNIRN   1
#define MAXPRNIRN   14
#define NSATIRN     (MAXPRNIRN-MINPRNIRN+1)
#define MINPRNLEO   1
#define MAXPRNLEO   10
#define NSATLEO     (MAXPRNLEO-MINPRNLEO+1)
#define MINPRNSBS   120
#define MAXPRNSBS   158
#define NSATSBS     (MAXPRNSBS-MINPRNSBS+1)

// ====================== 常用频率与地球自转角速度 ======================
#define FREQ_GPS_L1 1575.42e6
#define FREQ_GPS_L2 1227.60e6

#define FREQ_BDS_B1 1561.098e6
#define FREQ_BDS_B2 1207.14e6
#define FREQ_BDS_B3 1268.52e6

#define OMGed_GPS 7.2921151467e-5
#define OMGed_BDS 7.2921150e-5
