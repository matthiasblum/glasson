#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <search.h>
#include <stdlib.h>
#include "zlib.h"

typedef struct chromsize_t {
    uint16_t n;
    char **chrom;
    uint32_t *size;
};

uint32_t read_uint32(uint8_t *buffer) {
    // LE
    return buffer[0] | (buffer[1]<<8) | (buffer[2]<<16) | (buffer[3]<<24);

    // BE
    //return buffer[3] | (buffer[2]<<8) | (buffer[1]<<16) | (buffer[0]<<24);
}

void write_uint32(uint32_t i, FILE *fp) {
    uint8_t buffer[4];

    // BE
    // buffer[0] = (i >> 24) & 0xFF;
    // buffer[1] = (i >> 16) & 0xFF;
    // buffer[2] = (i >> 8) & 0xFF;
    // buffer[3] = i & 0xFF;

    // LE
    buffer[0] = i & 0xFF;
    buffer[1] = (i >> 8) & 0xFF;
    buffer[2] = (i >> 16) & 0xFF;
    buffer[3] = (i >> 24) & 0xFF;

    fwrite(&buffer, sizeof(uint8_t), 4, fp);
}

void write_uint64(uint64_t i, FILE *fp) {
    uint8_t buffer[8];

    // LE
    buffer[0] = i & 0xFF;
    buffer[1] = (i >> 8) & 0xFF;
    buffer[2] = (i >> 16) & 0xFF;
    buffer[3] = (i >> 24) & 0xFF;
    buffer[4] = (i >> 32) & 0xFF;
    buffer[5] = (i >> 40) & 0xFF;
    buffer[6] = (i >> 48) & 0xFF;
    buffer[7] = (i >> 56) & 0xFF;

    fwrite(&buffer, sizeof(uint8_t), 8, fp);
}

int main(int argc, char const *argv[]) {
    FILE *fp = fopen("test.dat", "w");
    uint32_t i = 1989;

    char *compr;
    uint64_t destLen;

    char string[] = "Hello, world!";
    unsigned long sourceLen = strlen(string);
    destLen = compressBound(sourceLen);
    compr = malloc(destLen * sizeof(char));
    int ret = compress2 (compr, &destLen, string, sourceLen, Z_DEFAULT_COMPRESSION);

    write_uint32(destLen, fp);
    fwrite(compr, 1, destLen, fp);
    free(compr);

    uint8_t *bytes = malloc(10 * sizeof(uint8_t));
    for (i = 0; i < 10; i++) {
        bytes[i] = i * 2;
    }

    sourceLen = sizeof(uint8_t) * 10;
    destLen = compressBound(sourceLen);
    compr = malloc(destLen * sizeof(char));
    ret = compress2 (compr, &destLen, (const Bytef*)bytes, sourceLen, Z_DEFAULT_COMPRESSION);
    free(bytes);
    if (ret != Z_OK) {
        printf("nope\n");
    }

    write_uint32(destLen, fp);
    fwrite(compr, 1, destLen, fp);
    free(compr);

    uint32_t *values = malloc(10 * sizeof(uint32_t));
    for (i = 0; i < 10; i++) {
        values[i] = i * 10;
    }

    bytes = malloc(10 * 4 * sizeof(uint8_t));
    for (i = 0; i < 10; i++) {
        write_uint32(values[i], fp);
        bytes[i*4] = values[i] & 0xFF;
        bytes[i*4+1] = (values[i] >> 8) & 0xFF;
        bytes[i*4+2] = (values[i] >> 16) & 0xFF;
        bytes[i*4+3] = (values[i] >> 24) & 0xFF;
    }
    free(values);

    sourceLen = 10 * 4 * sizeof(uint8_t);
    destLen = compressBound(sourceLen);
    compr = malloc(destLen * sizeof(char));
    printf("%lu\n", compressBound(sourceLen));
    ret = compress2 (compr, &destLen, (const Bytef*)bytes, sourceLen, Z_DEFAULT_COMPRESSION);
    free(bytes);

    write_uint32(destLen, fp);
    fwrite(compr, 1, destLen, fp);

    free(compr);

    fclose(fp);
    return 0;
}
