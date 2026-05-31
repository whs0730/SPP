#pragma once

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif

#include "decode.h"

#include <winsock2.h>
#include <cstdint>
#include <string>

#pragma comment(lib, "Ws2_32.lib")

class TcpOem4Stream {
public:
    TcpOem4Stream() = default;
    ~TcpOem4Stream();

    bool Open(const char* ip, unsigned short port, int recv_timeout_ms = 5000);
    void Close();

    bool ReadByte(uint8_t* data);
    bool ReadExact(uint8_t* data, int len);
    int RecvPacket(uint8_t* data, int max_len, int sleep_ms = 980);

    const std::string& LastError() const;

private:
    void SetLastSocketError(const char* prefix);

    SOCKET sock_ = INVALID_SOCKET;
    bool wsa_started_ = false;
    std::string last_error_;
};

int input_oem4s(raw_t* raw, TcpOem4Stream* stream);
