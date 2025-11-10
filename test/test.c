#include "../bc_crunch.h"

#define ARENA_IMPLEMENTATION
#include "arena.h"

#include <string.h>
#define STB_DXT_IMPLEMENTATION
#include "stb_dxt.h"

#define STBI_ONLY_PNG
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <stdbool.h>

Arena g_arena = {};
float global_ratio = 0.f;
uint32_t num_ratios = 0;

//----------------------------------------------------------------------------------------------------------------------------
static inline int iq_random( int* seed)
{
    *seed *= 16807;
    return (*seed) >> 9;
}

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
float test_bc1(const uint8_t* rgba, uint32_t width, uint32_t height)
{
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
    size_t worst_case = bc1_size * 10;
    void* crunched_texture = arena_alloc(&g_arena, worst_case);
    size_t crunched_size = bc1_crunch(bc1_texture, width, height, crunched_texture, worst_case);
    float ratio = (float) bc1_size / (float) crunched_size;

    fprintf(stdout, "BC1 size %u bytes => crunched size %zu bytes\ncompression ratio : %f\n", bc1_size, crunched_size, ratio); 

    uint8_t* uncompressed_texture = arena_alloc(&g_arena, bc1_size);
    bc1_decrunch(crunched_texture, crunched_size, width, height, uncompressed_texture);

    fprintf(stdout, "comparing decrushed vs original : ");
    for(uint32_t i=0; i<bc1_size; ++i)
    {
        if (uncompressed_texture[i] != bc1_texture[i])
        {
            fprintf(stdout, "failed, divergence at the %uth bytes\n", i);
            return -1.f;
        }
    }

    fprintf(stdout, "ok\n");
    return ratio;
}

//-----------------------------------------------------------------------------------------------------------------------------
static inline uint32_t lerp_color(uint32_t a, uint32_t b, float t)
{
    int tt = (int)(t * 256.f);
    int oneminust = 256 - tt;

    uint32_t A = (((a >> 24) & 0xFF) * oneminust + ((b >> 24) & 0xFF) * tt) >> 8;
    uint32_t B = (((a >> 16) & 0xFF) * oneminust + ((b >> 16) & 0xFF) * tt) >> 8;
    uint32_t G = (((a >> 8)  & 0xFF) * oneminust + ((b >> 8)  & 0xFF) * tt) >> 8;
    uint32_t R = (((a >> 0)  & 0xFF) * oneminust + ((b >> 0)  & 0xFF) * tt) >> 8;

    return (A << 24) | (B << 16) | (G << 8) | R;
}

//-----------------------------------------------------------------------------------------------------------------------------
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

//-----------------------------------------------------------------------------------------------------------------------------
void make_radial(uint32_t* p, uint32_t w, uint32_t h, uint32_t a, uint32_t b)
{
    float cx = (float)w * 0.5f;
    float cy = (float)h * 0.5f;
    float maxd = sqrtf(cx*cx + cy*cy);

    for (uint32_t y = 0; y < h; y++)
        for (uint32_t x = 0; x < w; x++) 
        {
            float dx = (float)x - cx;
            float dy = (float)y - cy;
            float t = sqrtf(dx*dx + dy*dy) / maxd;
            if (t > 1.f) t = 1.f;
            *p++ = lerp_color(a, b, t);
        }
}

//-----------------------------------------------------------------------------------------------------------------------------
void make_random(uint32_t* p, uint32_t w, uint32_t h)
{
    int seed = 0x12345678;
     for (uint32_t y = 0; y < h; y++)
        for (uint32_t x = 0; x < w; x++) 
            *p++ = iq_random(&seed);
}

//-----------------------------------------------------------------------------------------------------------------------------
void test_multiple(const char* format, uint32_t range_min, uint32_t range_max)
{
    for(uint32_t i=range_min; i<=range_max; ++i)
    {
        char filename[256];
        snprintf(filename, 256, format, i);

        fprintf(stdout, "\nopening %s ", filename);

        int width, height, num_channels;
        unsigned char *data = stbi_load(filename, &width, &height, &num_channels, 4);

        fprintf(stdout, "width : %d height :%d num_channels :%d\n", width, height, num_channels);

        arena_reset(&g_arena);

        float ratio = test_bc1(data, width, height);
        if (ratio >= 0.f)
        {
            global_ratio += ratio;
            num_ratios++;
        }

        stbi_image_free(data);
    }
}

#define TEXTURE_WIDTH (512)
#define TEXTURE_HEIGHT (512)

//-----------------------------------------------------------------------------------------------------------------------------
int main(void)
{
    fprintf(stdout, "bc_crunch test suite\n\n");

    // synthetic textures tests
    fprintf(stdout, "synthetic textures tests\n");
    uint8_t* rgba = arena_alloc(&g_arena, TEXTURE_WIDTH * TEXTURE_HEIGHT * 4);
    uniform_texture(rgba, TEXTURE_WIDTH, TEXTURE_HEIGHT, 0xef, 0x7d, 0x57);
    if (test_bc1(rgba, TEXTURE_WIDTH, TEXTURE_HEIGHT) < 0.f)
        return -1;

    make_radial((uint32_t*)rgba, TEXTURE_WIDTH, TEXTURE_HEIGHT, 0xffff1010, 0xff1010ff);
    if (test_bc1(rgba, TEXTURE_WIDTH, TEXTURE_HEIGHT) < 0.f)
        return -1;

    make_random((uint32_t*)rgba, TEXTURE_WIDTH, TEXTURE_HEIGHT);
    if (test_bc1(rgba, TEXTURE_WIDTH, TEXTURE_HEIGHT) < 0.f)
        return -1;

    arena_reset(&g_arena);

    // kodak photos
    test_multiple("../textures/kodim%02u.png", 1, 5);

    test_multiple("../textures/Elements_%02u-512x512.png", 1, 6);

    test_multiple("../textures/Dirt_%02u-512x512.png", 12, 17);

    test_multiple("../textures/Wood_%02u-512x512.png", 4, 7);

    test_multiple("../textures/Brick_%02d-512x512.png", 10, 14);

    fprintf(stdout, "\nglobal ratio : %f\n", global_ratio / (float) num_ratios);

    size_t bytes_allocated, byte_used;
    arena_stats(&g_arena, &bytes_allocated, &byte_used);
    fprintf(stdout, "\n%zu kb allocated\n%zu kb used\n\n", bytes_allocated>>10, byte_used>>10);

    arena_free(&g_arena);

    return 0;
}
