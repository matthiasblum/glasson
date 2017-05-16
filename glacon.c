#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <search.h>
#include <stdlib.h>
#include "zlib.h"

// typedef struct chromsize_t {
//     uint16_t n;
//     char **chrom;
//     uint32_t *size;
// };

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
    buffer[4] = (i >> 32) & 0xFF;
    buffer[5] = (i >> 40) & 0xFF;
    buffer[6] = (i >> 48) & 0xFF;
    buffer[7] = (i >> 56) & 0xFF;

    fwrite(&buffer, sizeof(uint8_t), 8, fp);
}

int main(int argc, char const *argv[]) {
    FILE *fp = fopen("test.dat", "w");
    uint32_t i = 1989;

    //write_uint32(fp, i);

    char string[] = "Hello, world!";
    unsigned long sourceLen = strlen(string);
    uint64_t destLen = compressBound(sourceLen);
    char *buffer = malloc(destLen * sizeof(char));
    int ret = compress2 (buffer, &destLen, string, sourceLen, Z_DEFAULT_COMPRESSION);

    write_uint64(destLen, fp);
    fwrite(buffer, 1, destLen, fp);

    uint32_t *buffer2 = malloc(10 * sizeof(uint32_t));
    for (i = 0; i < 10; i++) {
        buffer2[i] = i * 10;
    }

    ret = compress2 (buffer, &destLen, (const Bytef*)buffer2, sizeof(buffer2), Z_DEFAULT_COMPRESSION);
    write_uint64(destLen, fp);
    fwrite(buffer, 1, destLen, fp);

    free(buffer);
    fclose(fp);

    /*FILE *fp;
    fp = fopen("test.dat", "r");

    uint8_t buf[4];
    fread(&buf, sizeof(uint8_t), 4, fp);

    uint32_t i = read_uint32(buf);
    printf("%u\n", i);

    char *zbuffer = malloc(i * sizeof(char));

    fread(zbuffer, i, 1, fp);

    free(zbuffer);

    fclose(fp);
    return 0;*/

    // FILE *fp = fopen("test.dat", "w");
    //
    // z_stream strm;
    // int ret, flush;
    // /* allocate deflate state */
    // strm.zalloc = Z_NULL;
    // strm.zfree = Z_NULL;
    // strm.opaque = Z_NULL;
    // ret = deflateInit(&strm, level);
    // if (ret != Z_OK)
    //     return ret;



    // fread(buffer, 11, 1, fp);
    //
    // printf("%s\n", buffer);
    //
    // z_stream defstream;
    // z_stream infstream;
    // infstream.zalloc = Z_NULL;
    // infstream.zfree = Z_NULL;
    // infstream.opaque = Z_NULL;
    // setup "b" as the input and "c" as the compressed output
    // infstream.avail_in = (uInt)((char*)defstream.next_out - buffer); // size of input
    // infstream.next_in = (Bytef *)buffer; // input char array
    // infstream.avail_out = (uInt)sizeof(unc_buffer); // size of output
    // infstream.next_out = (Bytef *)unc_buffer; // output char array
    //
    // // the actual DE-compression work.
    // inflateInit(&infstream);
    // inflate(&infstream, Z_NO_FLUSH);
    // inflateEnd(&infstream);
    //
    // printf("%s\n", unc_buffer);
    // fclose(fp);
    return 0;
}
