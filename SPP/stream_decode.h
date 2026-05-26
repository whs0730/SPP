#pragma once

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif

#include "decode.h"

#include <winsock2.h>
#include <ws2tcpip.h>

#include <cstdint>
#include <iostream>
#include <string>

#pragma comment(lib, "Ws2_32.lib")

// NovAtel OEM4 二进制实时流的 TCP 读取器。
// 这里只负责 socket 管理和按字节读取；OEM4 组帧和解码仍放在
// input_oem4s()/decode_oem4() 中，保证实时流和文件回放走同一套流程。
class TcpOem4Stream {
public:
    TcpOem4Stream() = default;

    ~TcpOem4Stream()
    {
        Close();
    }

    bool Open(const char* ip, unsigned short port, int recv_timeout_ms = 5000)
    {
        Close();

        // 使用 Windows socket API 前必须先初始化 Winsock。
        WSADATA wsa_data{};
        int wsa_ret = WSAStartup(MAKEWORD(2, 2), &wsa_data);
        if (wsa_ret != 0) {
            last_error_ = "WSAStartup failed: " + std::to_string(wsa_ret);
            return false;
        }
        wsa_started_ = true;

        sock_ = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP);
        if (sock_ == INVALID_SOCKET) {
            SetLastSocketError("socket failed");
            Close();
            return false;
        }

        sockaddr_in addr{};
        addr.sin_family = AF_INET;
        addr.sin_port = htons(port);

        if (inet_pton(AF_INET, ip, &addr.sin_addr) != 1) {
            last_error_ = "invalid IPv4 address: ";
            last_error_ += ip;
            Close();
            return false;
        }

        if (connect(sock_, reinterpret_cast<sockaddr*>(&addr), sizeof(addr)) == SOCKET_ERROR) {
            SetLastSocketError("connect failed");
            Close();
            return false;
        }

        // recv() 可能一次读不到请求的全部字节；设置超时可以避免断流时
        // ReadExact() 一直阻塞。
        setsockopt(sock_, SOL_SOCKET, SO_RCVTIMEO,
            reinterpret_cast<const char*>(&recv_timeout_ms), sizeof(recv_timeout_ms));

        last_error_.clear();
        return true;
    }

    void Close()
    {
        if (sock_ != INVALID_SOCKET) {
            closesocket(sock_);
            sock_ = INVALID_SOCKET;
        }

        if (wsa_started_) {
            WSACleanup();
            wsa_started_ = false;
        }
    }

    bool ReadByte(uint8_t* data)
    {
        return ReadExact(data, 1);
    }

    bool ReadExact(uint8_t* data, int len)
    {
        if (data == nullptr || len <= 0 || sock_ == INVALID_SOCKET) {
            last_error_ = "stream is not open";
            return false;
        }

        int got = 0;
        while (got < len) {
            // TCP 是连续字节流，不保留报文边界，因此需要循环 recv，
            // 直到凑齐调用方要求的字节数。
            int n = recv(sock_, reinterpret_cast<char*>(data + got), len - got, 0);
            if (n > 0) {
                got += n;
                continue;
            }

            if (n == 0) {
                last_error_ = "remote host closed the connection";
                return false;
            }

            int err = WSAGetLastError();
            // 超时或中断按临时状态处理；实时数据有短暂间隔时继续等待。
            if (err == WSAEINTR || err == WSAETIMEDOUT) {
                continue;
            }

            last_error_ = "recv failed: " + std::to_string(err);
            return false;
        }

        return true;
    }

    const std::string& LastError() const
    {
        return last_error_;
    }

private:
    void SetLastSocketError(const char* prefix)
    {
        last_error_ = std::string(prefix) + ": " + std::to_string(WSAGetLastError());
    }

    SOCKET sock_ = INVALID_SOCKET;
    bool wsa_started_ = false;
    std::string last_error_;
};

// input_oem4f() 的 socket 输入版本。返回值沿用 decode_oem4()：
// 1 = 观测历元，2 = 星历，0 = 跳过消息，负数 = 错误或结束。
static int input_oem4s(raw_t* raw, TcpOem4Stream* stream)
{
    if (raw == nullptr || stream == nullptr) {
        return -2;
    }

    uint8_t data = 0;

    if (raw->nbyte == 0) {
        // 在连续 TCP 字节流中查找 OEM4 同步头 AA 44 12。
        // 同步头之前的字节视为噪声或残缺报文并丢弃。
        for (int i = 0;; i++) {
            if (!stream->ReadByte(&data)) {
                return -2;
            }
            if (sync_oem4(raw->buff, data)) {
                break;
            }
            if (i >= 4096) {
                return 0;
            }
        }
    }

    // 已经读到 3 字节同步头；再读够头部字节，取得 buff + 8 处的报文长度字段。
    if (!stream->ReadExact(raw->buff + 3, 7)) {
        return -2;
    }
    raw->nbyte = 10;

    raw->len = U2(raw->buff + 8) + OEM4HLEN;
    if (raw->len > MAXRAWLEN - 4) {
        raw->nbyte = 0;
        return -1;
    }

    // 继续读取剩余头部、消息体和 CRC，然后把完整 OEM4 帧交给文件模式共用的解码器。
    if (!stream->ReadExact(raw->buff + 10, raw->len - 6)) {
        return -2;
    }
    raw->nbyte = 0;

    return decode_oem4(raw);
}
