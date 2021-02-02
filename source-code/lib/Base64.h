#ifndef BASE64_H
#define BASE64_H

#include <math.h>
#include <stdlib.h>

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

size_t textToBinarySize(size_t charCount)
{
    return (size_t)ceil(charCount / 4.0) * 3;
}

void textToBinary(const char* from, size_t charCount, unsigned char* result)
{
    if (!byteTableInitialized)
        initializeByteTable();

    unsigned char* it = (unsigned char*)from;
    unsigned char* end = it + charCount;
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
            unsigned char temp = static_cast<unsigned char>(int24bit >> ((2 - i) * 8));
            *result++ = temp;
            int24bit ^= temp << ((2 - i) * 8);
            written++;
        }
    }

}

#endif