#include "stream_decode.h"

#include <iostream>

TcpOem4Stream::~TcpOem4Stream()
{
    Close();
}

// 打开实时 TCP 数据流，成功后后续按字节读取 OEM4 数据。
bool TcpOem4Stream::Open(const char* ip, unsigned short port, int recv_timeout_ms)
{
    Close();

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

    // 设置接收超时，避免实时流中断时一直阻塞。
    setsockopt(sock_, SOL_SOCKET, SO_RCVTIMEO,
        reinterpret_cast<const char*>(&recv_timeout_ms), sizeof(recv_timeout_ms));

    last_error_.clear();
    return true;
}

// 关闭 socket，并在需要时释放 Winsock 环境。
void TcpOem4Stream::Close()
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

// 读取单个字节，供同步头扫描使用。
bool TcpOem4Stream::ReadByte(uint8_t* data)
{
    return ReadExact(data, 1);
}

// 循环读取指定长度的数据，直到读满或出现网络错误。
bool TcpOem4Stream::ReadExact(uint8_t* data, int len)
{
    if (data == nullptr || len <= 0 || sock_ == INVALID_SOCKET) {
        last_error_ = "stream is not open";
        return false;
    }

    int got = 0;
    while (got < len) {
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
        if (err == WSAEINTR || err == WSAETIMEDOUT) {
            continue;
        }

        last_error_ = "recv failed: " + std::to_string(err);
        return false;
    }

    return true;
}

const std::string& TcpOem4Stream::LastError() const
{
    return last_error_;
}

// 记录最近一次 socket 错误，便于主函数输出具体原因。
void TcpOem4Stream::SetLastSocketError(const char* prefix)
{
    last_error_ = std::string(prefix) + ": " + std::to_string(WSAGetLastError());
}

// 实时输入逻辑：扫描同步头、读取完整报文、交给 decode_oem4() 解码。
int input_oem4s(raw_t* raw, TcpOem4Stream* stream)
{
    if (raw == nullptr || stream == nullptr) {
        return -2;
    }

    uint8_t data = 0;

    if (raw->nbyte == 0) {
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

    if (!stream->ReadExact(raw->buff + 3, 7)) {
        return -2;
    }
    raw->nbyte = 10;

    raw->len = U2(raw->buff + 8) + OEM4HLEN;
    if (raw->len > MAXRAWLEN - 4) {
        raw->nbyte = 0;
        return -1;
    }

    if (!stream->ReadExact(raw->buff + 10, raw->len - 6)) {
        return -2;
    }
    raw->nbyte = 0;

    return decode_oem4(raw);
}
