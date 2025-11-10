#include "../bc_crunch.h"

#define ARENA_IMPLEMENTATION
#include "arena.h"

#include <string.h>
#define STB_DXT_IMPLEMENTATION
#include "stb_dxt.h"

#include <stdbool.h>

Arena g_arena = {};

//----------------------------------------------------------------------------------------------------------------------------
void extract_4x4_block(const uint8_t* rgba, uint32_t width, uint32_t x, uint32_t y, uint8_t block[64])
{
    for (uint32_t j = 0; j < 4; j++)
    {
        const uint8_t* src = &rgba[(y + j) * width * 4 + x * 4];
        memcpy(&block[j * 4 * 4], src, 16);
    }
}

//----------------------------------------------------------------------------------------------------------------------------
// Take an rgba image in input (width and height multiple of 4)
//  - compress it to bc1 using stb_dxt
//  - crunch the bc1 blocks into a stream
//  - compare size of bc1 texture vs crunched texture
//  - decompress the stream and check it's identical to the bc1 texture
bool test_bc1(const uint8_t* rgba, uint32_t width, uint32_t height)
{
    arena_reset(&g_arena);

    uint32_t bc1_size = (width * height) / 2;
    uint8_t* bc1_texture = arena_alloc(&g_arena, bc1_size);

    uint32_t block_index=0;
    for(uint32_t y=0; y<height/4; ++y)
    {
        for(uint32_t x=0; x<width/4; ++x)
        {
            uint8_t block_image[64];
            extract_4x4_block(rgba, width, x*4, y*4, block_image);
            stb_compress_dxt_block(&bc1_texture[8 * block_index], block_image, 0, STB_DXT_HIGHQUAL);
            block_index++;
        }
    }

    // TODO : add profiling
    size_t worst_case = bc1_size * 2;
    void* crunched_texture = arena_alloc(&g_arena, worst_case);
    size_t crunched_size = bc1_crunch(bc1_texture, width, height, crunched_texture, worst_case);

    fprintf(stdout, "BC1 size %u bytes => crunched size %zu bytes\ncompression ratio : %f\n", bc1_size, crunched_size, (double) bc1_size / (double) crunched_size); 

    uint8_t* uncompressed_texture = arena_alloc(&g_arena, bc1_size);
    bc1_decrunch(crunched_texture, crunched_size, width, height, uncompressed_texture);

    fprintf(stdout, "comparing decrushed vs original : ");
    for(uint32_t i=0; i<bc1_size; ++i)
    {
        if (uncompressed_texture[i] != bc1_texture[i])
        {
            fprintf(stdout, "failed\n");
            return false;
        }
    }

    fprintf(stdout, "ok\n");
    return true;
}


void uniform_texture(uint8_t* rgba, uint32_t width, uint32_t height, uint8_t red, uint8_t green, uint8_t blue)
{
    uint32_t row_bytes = width * 4;
    for(uint32_t y=0; y<height; ++y)
    {
        for(uint32_t x=0; x<width; ++x)
        {
            uint32_t index = y * row_bytes + x * 4;
            rgba[index+0] = red;
            rgba[index+1] = green;
            rgba[index+2] = blue;
        }
    }
}

#define TEXTURE_WIDTH (512)
#define TEXTURE_HEIGHT (256)

int main(void)
{
    fprintf(stdout, "bc_crunch test suite\n\n");

    uint8_t* rgba = arena_alloc(&g_arena, TEXTURE_WIDTH * TEXTURE_HEIGHT * 4);

    uniform_texture(rgba, TEXTURE_WIDTH, TEXTURE_HEIGHT, 0xef, 0x7d, 0x57);
    if (!test_bc1(rgba, TEXTURE_WIDTH, TEXTURE_HEIGHT))
        return -1;


    size_t bytes_allocated, byte_used;
    arena_stats(&g_arena, &bytes_allocated, &byte_used);
    fprintf(stdout, "\n%zu kb allocated\n%zu kb used\n\n", bytes_allocated>>10, byte_used>>10);

    arena_free(&g_arena);

    return 0;
}
