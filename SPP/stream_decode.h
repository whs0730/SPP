#pragma once

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif

#include "decode.h"

#include <winsock2.h>
#include <ws2tcpip.h>

#include <cstdint>
#include <string>

#pragma comment(lib, "Ws2_32.lib")

// NovAtel OEM4 实时 TCP 流读取器。
// 只负责网络连接和按字节读取，报文同步与解码仍复用 decode.cpp。
class TcpOem4Stream {
public:
    TcpOem4Stream() = default;
    ~TcpOem4Stream();

    // 连接实时数据源，并设置接收超时时间。
    bool Open(const char* ip, unsigned short port, int recv_timeout_ms = 5000);

    // 关闭 socket 并释放 Winsock。
    void Close();

    // 读取单字节或指定长度字节。
    bool ReadByte(uint8_t* data);
    bool ReadExact(uint8_t* data, int len);

    // 返回最近一次网络错误描述。
    const std::string& LastError() const;

private:
    void SetLastSocketError(const char* prefix);

    SOCKET sock_ = INVALID_SOCKET;
    bool wsa_started_ = false;
    std::string last_error_;
};

// 实时流输入：从 TCP 流中拼出完整 OEM4 报文后调用 decode_oem4()。
int input_oem4s(raw_t* raw, TcpOem4Stream* stream);
