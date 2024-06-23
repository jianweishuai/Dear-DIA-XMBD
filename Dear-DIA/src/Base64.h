#ifndef BASE_H
#define BASE_H

#include <math.h>
#include <vector>
#include <stdlib.h>

namespace Base64 
{
    typedef unsigned char byte;

    char charTable[64] =
    {
        'A','B','C','D','E','F','G','H','I','J',
        'K','L','M','N','O','P','Q','R','S','T',
        'U','V','W','X','Y','Z','a','b','c','d',
        'e','f','g','h','i','j','k','l','m','n',
        'o','p','q','r','s','t','u','v','w','x',
        'y','z','0','1','2','3','4','5','6','7',
        '8','9', '+','/'
    };


    char byteTable[256];
    bool byteTableInitialized = false;

    void initializeByteTable()
    {
        for (size_t i = 0; i < 64; i++)
            byteTable[static_cast<int>(charTable[i])] = static_cast<char>(i);

        byteTableInitialized = true;
    }

    size_t binaryToTextSize(size_t byteCount);
    size_t binaryToText(const void* from, size_t byteCount, char* to);
    size_t textToBinarySize(size_t charCount);
    size_t textToBinary(const char* from, size_t charCount, void* to);
}

size_t Base64::binaryToTextSize(size_t byteCount)
{
    return (size_t)ceil(byteCount / 3.) * 4;
}

size_t Base64::binaryToText(const void* from, size_t byteCount, char* to)
{
    byte* it = (byte*)from;
    byte* end = it + byteCount;
    size_t written = 0;

    while (it != end)
    {
        int int24bit = 0;
        int paddingCount = 0;

        // construct 24-bit integer from 3 bytes
        for (int i = 0; i < 3; i++)
        {
            if (it != end)
                int24bit |= *it++ << ((2 - i) * 8);
            else
                paddingCount++;
        }

        // write out 4 characters
        for (int i = 3; i >= 0; i--)
        {
            to[i] = charTable[int24bit & 0x3F];
            int24bit >>= 6;
        }

        // fixup for padding
        if (paddingCount > 0)
            to[3] = '=';
        if (paddingCount > 1)
            to[2] = '=';

        to += 4;
        written += 4;
    }

    return written;
}

size_t Base64::textToBinarySize(size_t charCount)
{
    return (size_t)ceil(charCount / 4.) * 3;
}

size_t Base64::textToBinary(const char* from, size_t charCount, void* to)
{
    if (!byteTableInitialized)
        initializeByteTable();

    byte* it = (byte*)from;
    byte* end = it + charCount;
    byte* result = (byte*)to;
    size_t written = 0;

    while (it != end)
    {
        int int24bit = 0;
        int paddingCount = 0;

        // construct 24-bit integer from 4 characters
        for (int i = 0; i < 4 && it != end; i++, it++)
        {
            if (*it != '=')
            {
                int24bit |= byteTable[*it] << ((3 - i) * 6);
            }
            else
                paddingCount++;
        }

        // write out bytes
        for (int i = 0; i < 3 - paddingCount; i++)
        {
            byte temp = static_cast<byte>(int24bit >> ((2 - i) * 8));
            *result++ = temp;
            int24bit ^= temp << ((2 - i) * 8);
            written++;
        }
    }

    return written;
}

#endif