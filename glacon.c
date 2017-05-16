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
    return (buffer[0]<<0) | (buffer[1]<<8) | (buffer[2]<<16) | (buffer[3]<<24);

    // BE
    //return (buffer[3]<<0) | (buffer[2]<<8) | (buffer[1]<<16) | (buffer[0]<<24);

}

int main(int argc, char const *argv[]) {
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

    FILE *fp = fopen("test.dat", "w");

    z_stream strm;
    int ret, flush;
    /* allocate deflate state */
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    ret = deflateInit(&strm, level);
    if (ret != Z_OK)
        return ret;



    // fread(buffer, 11, 1, fp);
    //
    // printf("%s\n", buffer);
    //
    // z_stream defstream;
    z_stream infstream;
    infstream.zalloc = Z_NULL;
    infstream.zfree = Z_NULL;
    infstream.opaque = Z_NULL;
    // // setup "b" as the input and "c" as the compressed output
    infstream.avail_in = (uInt)((char*)defstream.next_out - buffer); // size of input
    infstream.next_in = (Bytef *)buffer; // input char array
    infstream.avail_out = (uInt)sizeof(unc_buffer); // size of output
    infstream.next_out = (Bytef *)unc_buffer; // output char array
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
