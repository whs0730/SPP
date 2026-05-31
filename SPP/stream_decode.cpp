#include "stream_decode.h"

#include <cstring>
#include <windows.h>

#pragma warning(disable:4996)

static bool OpenSocket(SOCKET& sock, const char IP[], const unsigned short Port)
{
    WSADATA wsaData;
    SOCKADDR_IN addrSrv;

    if (!WSAStartup(MAKEWORD(1, 1), &wsaData))
    {
        if ((sock = socket(AF_INET, SOCK_STREAM, 0)) != INVALID_SOCKET)
        {
            addrSrv.sin_addr.S_un.S_addr = inet_addr(IP);
            addrSrv.sin_family = AF_INET;
            addrSrv.sin_port = htons(Port);
            if (connect(sock, (SOCKADDR*)&addrSrv, sizeof(SOCKADDR)) == SOCKET_ERROR)
            {
                closesocket(sock);
                sock = INVALID_SOCKET;
                WSACleanup();
                return false;
            }
            return true;
        }
        WSACleanup();
    }
    return false;
}

static void CloseSocket(SOCKET& sock)
{
    closesocket(sock);
    WSACleanup();
}

TcpOem4Stream::~TcpOem4Stream()
{
    Close();
}

bool TcpOem4Stream::Open(const char* ip, unsigned short port, int recv_timeout_ms)
{
    Close();

    if (!OpenSocket(sock_, ip, port)) {
        SetLastSocketError("OpenSocket failed");
        Close();
        return false;
    }
    wsa_started_ = true;

    setsockopt(sock_, SOL_SOCKET, SO_RCVTIMEO,
        reinterpret_cast<const char*>(&recv_timeout_ms), sizeof(recv_timeout_ms));

    last_error_.clear();
    return true;
}

void TcpOem4Stream::Close()
{
    if (sock_ != INVALID_SOCKET || wsa_started_) {
        CloseSocket(sock_);
        sock_ = INVALID_SOCKET;
        wsa_started_ = false;
    }
}

bool TcpOem4Stream::ReadByte(uint8_t* data)
{
    return ReadExact(data, 1);
}

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

int TcpOem4Stream::RecvPacket(uint8_t* data, int max_len, int sleep_ms)
{
    if (data == nullptr || max_len <= 0 || sock_ == INVALID_SOCKET) {
        last_error_ = "stream is not open";
        return -1;
    }

    if (sleep_ms > 0) {
        Sleep(sleep_ms);
    }

    int n = recv(sock_, reinterpret_cast<char*>(data), max_len, 0);
    if (n > 0) {
        return n;
    }

    if (n == 0) {
        last_error_ = "remote host closed the connection";
        return -1;
    }

    int err = WSAGetLastError();
    if (err == WSAEINTR || err == WSAETIMEDOUT || err == WSAEWOULDBLOCK) {
        return 0;
    }

    last_error_ = "recv failed: " + std::to_string(err);
    return -1;
}

const std::string& TcpOem4Stream::LastError() const
{
    return last_error_;
}

void TcpOem4Stream::SetLastSocketError(const char* prefix)
{
    last_error_ = std::string(prefix) + ": " + std::to_string(WSAGetLastError());
}

static bool TryDecodeOem4FromStreamBuff(raw_t* raw, uint8_t* Buff, int& lenD, int& ReadFlag)
{
    if (raw == nullptr) {
        ReadFlag = -2;
        return true;
    }

    while (lenD >= 3) {
        int sync_pos = -1;

        for (int i = 0; i <= lenD - 3; i++) {
            if (Buff[i] == OEM4SYNC1 &&
                Buff[i + 1] == OEM4SYNC2 &&
                Buff[i + 2] == OEM4SYNC3) {
                sync_pos = i;
                break;
            }
        }

        if (sync_pos < 0) {
            Buff[0] = Buff[lenD - 2];
            Buff[1] = Buff[lenD - 1];
            lenD = 2;
            return false;
        }

        if (sync_pos > 0) {
            lenD -= sync_pos;
            memmove(Buff, Buff + sync_pos, lenD);
        }

        if (lenD < 10) {
            return false;
        }

        raw->len = U2(Buff + 8) + OEM4HLEN;
        if (raw->len > MAXRAWLEN - 4) {
            lenD--;
            memmove(Buff, Buff + 1, lenD);
            raw->nbyte = 0;
            ReadFlag = -1;
            return true;
        }

        int frame_len = raw->len + 4;
        if (lenD < frame_len) {
            return false;
        }

        memcpy(raw->buff, Buff, frame_len);
        raw->nbyte = 0;

        lenD -= frame_len;
        if (lenD > 0) {
            memmove(Buff, Buff + frame_len, lenD);
        }

        ReadFlag = decode_oem4(raw);
        return true;
    }

    return false;
}

int input_oem4s(raw_t* raw, TcpOem4Stream* stream)
{
    static uint8_t buff[MAXRAWLEN];
    static uint8_t Buff[2 * MAXRAWLEN];
    static int lenD = 0;

    if (raw == nullptr || stream == nullptr) {
        return -2;
    }

    int lenR = 0;
    int ReadFlag = 0;

    if (TryDecodeOem4FromStreamBuff(raw, Buff, lenD, ReadFlag)) {
        return ReadFlag;
    }

    lenR = stream->RecvPacket(buff, MAXRAWLEN);
    if (lenR < 0) {
        return -2;
    }

    if (lenR <= 0) {
        return 0;
    }

    if ((lenD + lenR) > 2 * MAXRAWLEN) {
        lenD = 0;
    }

    memcpy(Buff + lenD, buff, lenR);
    lenD += lenR;

    if (TryDecodeOem4FromStreamBuff(raw, Buff, lenD, ReadFlag)) {
        return ReadFlag;
    }

    return 0;
}
