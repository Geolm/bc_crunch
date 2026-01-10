//----------------------------------------------------------------------------------------------------------------------------
// Disclaimer : this is badly written, the purpose of the file is just to test quickly progress and regression
//              Do not use this file as reference, do not copy this code, it's probably buggy and slow (yes both)
//----------------------------------------------------------------------------------------------------------------------------

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
#include "default_font_atlas.h"

#include "miniz.h"

Arena g_arena = {};
float global_ratio = 0.f;
float global_zlib_ratio = 0.f;
uint32_t num_ratios = 0;
void* cruncher_memory = NULL;

//----------------------------------------------------------------------------------------------------------------------------
static inline int iq_random( int* seed)
{
    *seed *= 16807;
    return (*seed) >> 9;
}

//----------------------------------------------------------------------------------------------------------------------------
static inline void extract_4x4_rgba_block(const uint8_t* rgba, uint32_t width, uint32_t x, uint32_t y, uint8_t block[64])
{
    for (uint32_t j = 0; j < 4; j++)
    {
        const uint8_t* src = &rgba[(y + j) * width * 4 + x * 4];
        memcpy(&block[j * 4 * 4], src, 16);
    }
}

//----------------------------------------------------------------------------------------------------------------------------
static inline void extract_4x4_red_block(const uint8_t* rgba, uint32_t width, uint32_t x, uint32_t y, uint8_t block[16], uint32_t channel)
{
    for (uint32_t j = 0; j < 4; j++)
    {
        for(uint32_t i=0; i<4; ++i)
        {
            block[j * 4 + i] = rgba[((y + j) * width + (x + i)) * 4 + channel];
        }
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
            extract_4x4_rgba_block(rgba, width, x*4, y*4, block_image);
            stb_compress_dxt_block(&bc1_texture[8 * block_index], block_image, 0, STB_DXT_HIGHQUAL);
            block_index++;
        }
    }

    // TODO : add profiling
    size_t worst_case = bc1_size * 10;
    void* crunched_texture = arena_alloc(&g_arena, worst_case);
    size_t crunched_size = bc_crunch(cruncher_memory, bc1_texture, width, height, bc1, crunched_texture, worst_case);
    float ratio = (float) bc1_size / (float) crunched_size;

    

    uint8_t* uncompressed_texture = arena_alloc(&g_arena, bc1_size);
    bc_decrunch(crunched_texture, crunched_size, width, height, bc1, uncompressed_texture);

    for(uint32_t i=0; i<bc1_size; ++i)
    {
        if (uncompressed_texture[i] != bc1_texture[i])
        {
            fprintf(stdout, "comparison vs original failed, divergence at the %uth bytes\n", i);
            return -1.f;
        }
    }

    mz_ulong zlib_size = worst_case;
    mz_compress2(crunched_texture, &zlib_size, bc1_texture, bc1_size, MZ_BEST_COMPRESSION);

    float zlib_ratio = (float) bc1_size / (float) zlib_size;
    global_zlib_ratio += zlib_ratio;

    fprintf(stdout, "BC1 size %u bytes => crunched size %zu bytes\ncompression ratio : %f vs zlib : %f\n", bc1_size, crunched_size, ratio, zlib_ratio); 

    return ratio;
}

//-----------------------------------------------------------------------------------------------------------------------------
bool test_bc3(const char* filename)
{
    uint32_t width, height, num_channels;

    arena_reset(&g_arena);

    fprintf(stdout, "\nopening %s\n", filename);

    unsigned char *data = stbi_load(filename, (int*)&width, (int*)&height, (int*)&num_channels, 4);

    fprintf(stdout, "width : %d height :%d num_channels :%d\n", width, height, num_channels);

    if (data == NULL || num_channels != 4)
        return false;

    uint32_t bc3_size = (width * height);
    uint8_t* bc3_texture = arena_alloc(&g_arena, bc3_size);

    uint32_t block_index=0;
    for(uint32_t y=0; y<height/4; ++y)
    {
        for(uint32_t x=0; x<width/4; ++x)
        {
            uint8_t block_image[128];
            extract_4x4_rgba_block(data, width, x*4, y*4, block_image);
            stb_compress_dxt_block(&bc3_texture[16 * block_index], block_image, 1, STB_DXT_HIGHQUAL);
            block_index++;
        }
    }

    void* crunched_texture = arena_alloc(&g_arena, bc3_size * 2);
    size_t crunched_size = bc_crunch(cruncher_memory, bc3_texture, width, height, bc3, crunched_texture, bc3_size*2);
    float ratio = (float) bc3_size / (float) crunched_size;

    fprintf(stdout, "BC3 size %u bytes => crunched size %zu bytes\ncompression ratio : %f\n", bc3_size, crunched_size, ratio); 

    uint8_t* uncompressed_texture = arena_alloc(&g_arena, bc3_size);
    bc_decrunch(crunched_texture, crunched_size, width, height, bc3, uncompressed_texture);

    fprintf(stdout, "comparing decrushed vs original : ");
    for(uint32_t i=0; i<bc3_size; ++i)
    {
        if (uncompressed_texture[i] != bc3_texture[i])
        {
            fprintf(stdout, "failed, divergence at the %uth bytes\n", i);
            return false;
        }
    }

    fprintf(stdout, "ok\n");
    stbi_image_free(data);

    return true;
}

//-----------------------------------------------------------------------------------------------------------------------------
bool test_bc4(const char* filename)
{
    uint32_t width, height, num_channels;

    arena_reset(&g_arena);

    fprintf(stdout, "\nopening %s\n", filename);

    unsigned char *data = stbi_load(filename, (int*)&width, (int*)&height, (int*)&num_channels, 4);

    fprintf(stdout, "width : %d height :%d num_channels :%d\n", width, height, num_channels);

    if (data == NULL)
        return false;

    uint32_t bc4_size = (width * height) / 2;
    uint8_t* bc4_texture = arena_alloc(&g_arena, bc4_size);

    uint32_t block_index=0;
    for(uint32_t y=0; y<height/4; ++y)
    {
        for(uint32_t x=0; x<width/4; ++x)
        {
            uint8_t block_image[64];
            extract_4x4_red_block(data, width, x*4, y*4, block_image, 0);
            stb_compress_bc4_block(&bc4_texture[8 * block_index], block_image);
            block_index++;
        }
    }

    size_t worst_case = bc4_size * 2;

    void* crunched_texture = arena_alloc(&g_arena, worst_case);
    size_t crunched_size = bc_crunch(cruncher_memory, bc4_texture, width, height, bc4, crunched_texture, worst_case);
    float ratio = (float) bc4_size / (float) crunched_size;

    global_ratio += ratio;
    num_ratios++;

    uint8_t* uncompressed_texture = arena_alloc(&g_arena, bc4_size);
    bc_decrunch(crunched_texture, crunched_size, width, height, bc4, uncompressed_texture);

    for(uint32_t i=0; i<bc4_size; ++i)
    {
        if (uncompressed_texture[i] != bc4_texture[i])
        {
            fprintf(stdout, "comparison vs original failed, divergence at the %uth bytes\n", i);
            return false;
        }
    }

    mz_ulong zlib_size = worst_case;
    mz_compress2(crunched_texture, &zlib_size, bc4_texture, bc4_size, MZ_BEST_COMPRESSION);

    float zlib_ratio = (float) bc4_size / (float) zlib_size;
    global_zlib_ratio += zlib_ratio;

    fprintf(stdout, "BC4 size %u bytes => crunched size %zu bytes\ncompression ratio : %f vs zlib : %f\n", bc4_size, crunched_size, ratio, zlib_ratio); 

    stbi_image_free(data);

    return true;
}

//-----------------------------------------------------------------------------------------------------------------------------
bool test_bc5(const char* filename)
{
    uint32_t width, height, num_channels;

    arena_reset(&g_arena);

    fprintf(stdout, "\nopening %s\n", filename);

    unsigned char *data = stbi_load(filename, (int*)&width, (int*)&height, (int*)&num_channels, 4);

    fprintf(stdout, "width : %d height :%d num_channels :%d\n", width, height, num_channels);

    if (data == NULL)
        return false;

    uint32_t bc5_size = (width * height);
    uint8_t* bc5_texture = arena_alloc(&g_arena, bc5_size);

    uint32_t block_index=0;
    for(uint32_t y=0; y<height/4; ++y)
    {
        for(uint32_t x=0; x<width/4; ++x)
        {
            uint8_t block_image[64];
            extract_4x4_red_block(data, width, x*4, y*4, block_image, 0);       // extract red
            stb_compress_bc4_block(&bc5_texture[16 * block_index], block_image);
            extract_4x4_red_block(data, width, x*4, y*4, block_image, 2);       // extract green
            stb_compress_bc4_block(&bc5_texture[16 * block_index + 8], block_image);
            block_index++;
        }
    }

    void* crunched_texture = arena_alloc(&g_arena, bc5_size * 2);
    size_t crunched_size = bc_crunch(cruncher_memory, bc5_texture, width, height, bc5, crunched_texture, bc5_size*2);
    float ratio = (float) bc5_size / (float) crunched_size;

    global_ratio += ratio;
    num_ratios++;

    fprintf(stdout, "BC5 size %u bytes => crunched size %zu bytes\ncompression ratio : %f\n", bc5_size, crunched_size, ratio); 

    uint8_t* uncompressed_texture = arena_alloc(&g_arena, bc5_size);
    bc_decrunch(crunched_texture, crunched_size, width, height, bc5, uncompressed_texture);

    fprintf(stdout, "comparing decrushed vs original : ");
    for(uint32_t i=0; i<bc5_size; ++i)
    {
        if (uncompressed_texture[i] != bc5_texture[i])
        {
            fprintf(stdout, "failed, divergence at the %uth bytes\n", i);
            return false;
        }
    }

    fprintf(stdout, "ok\n");
    stbi_image_free(data);

    return true;
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
bool test_multiple(const char* format, uint32_t range_min, uint32_t range_max)
{
    bool success = true;
    for(uint32_t i=range_min; i<=range_max; ++i)
    {
        char filename[256];
        snprintf(filename, 256, format, i);

        fprintf(stdout, "\nopening %s ", filename);

        int width, height, num_channels;
        unsigned char *data = stbi_load(filename, &width, &height, &num_channels, 4);

        if (data != NULL)
        {
            fprintf(stdout, "width : %d height :%d num_channels :%d\n", width, height, num_channels);

            arena_reset(&g_arena);

            float ratio = test_bc1(data, width, height);
            if (ratio >= 0.f)
            {
                global_ratio += ratio;
                num_ratios++;
            }
            else 
                success = false;

            stbi_image_free(data);
        }
        else
            fprintf(stdout, "failed\n");
    }
    return success;
}

#define TEXTURE_WIDTH (512)
#define TEXTURE_HEIGHT (512)
#define SMALL_SET

//-----------------------------------------------------------------------------------------------------------------------------
int main(void)
{
    fprintf(stdout, "bc_crunch test suite\n\n");

    cruncher_memory = malloc(crunch_min_size());

    fprintf(stdout, "-----------------------------------\n");
    fprintf(stdout, "| BC1 tests                       |\n");
    fprintf(stdout, "-----------------------------------\n\n");

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

    global_zlib_ratio = 0.f;

#ifdef SMALL_SET
    // kodak photos
    if (!test_multiple("../textures/bc1/kodim%02u.png", 1, 5))
        return -1;

    if (!test_multiple("../textures/bc1/Elements_%02u-512x512.png", 1, 6))
        return -1;

    if (!test_multiple("../textures/bc1/Dirt_%02u-512x512.png", 12, 17))
        return -1;

    if (!test_multiple("../textures/bc1/Wood_%02u-512x512.png", 4, 7))
        return -1;

    if (!test_multiple("../textures/bc1/Brick_%02d-512x512.png", 10, 14))
        return -1;

#else

    if (!test_multiple("../all/kodim%02u.png", 1, 24))
        return -1;

    if (!test_multiple("../all/Dirt_%02u-512x512.png", 1, 20))
        return -1;

    if (!test_multiple("../all/Wood_%02u-512x512.png", 1, 20))
        return -1;

    if (!test_multiple("../all/Elements_%02u-512x512.png", 1, 20))
        return -1;

    if (!test_multiple("../all/Brick_%02u-512x512.png", 1, 20))
        return -1;

    if (!test_multiple("../all/Metal_%02u-512x512.png", 1, 20))
        return -1;

    if (!test_multiple("../all/Plaster_%02u-512x512.png", 1, 20))
        return -1;

    if (!test_multiple("../all/Stone_%02u-512x512.png", 1, 20))
        return -1;
    
    if (!test_multiple("../all/Tile_%02u-512x512.png", 1, 20))
        return -1;

    if (!test_multiple("../all/{blaztree%u.png", 1, 7))
        return -1;

    if (!test_multiple("../all/{SW_Tree%02u.png", 1, 7))
        return -1;

    if (!test_multiple("../all/pkf_concrete%u.png", 1, 4))
        return -1;
    
    if (!test_multiple("../all/pk02_floor%02u_C.png", 1, 13))
        return -1;

    if (!test_multiple("../all/pk02_trim%02u_C.png", 1, 5))
        return -1;

    if (!test_multiple("../all/rock%02u.png", 1, 4))
        return -1;

    if (!test_multiple("../all/brick%u.png", 1, 2))
        return -1;

#endif

    float average_ratio = global_ratio / (float) num_ratios;
    fprintf(stdout, "\n\nBC1 average compression ratio : %f vs zlib ratio : %f\n\n", average_ratio, global_zlib_ratio / (float) num_ratios);

    if (average_ratio < (global_zlib_ratio / (float) num_ratios))
        return -1;

    fprintf(stdout, "-----------------------------------\n");
    fprintf(stdout, "| BC4 tests                       |\n");
    fprintf(stdout, "-----------------------------------\n\n");

    void* crunched = arena_alloc(&g_arena, default_font_atlas_size*2);
    size_t crunched_size = bc_crunch(cruncher_memory, default_font_atlas, 256, 256, bc4, crunched, default_font_atlas_size*2);

    fprintf(stdout, "satoshi font atlas bc4 = %zu ratio : %f\n", crunched_size, (float)default_font_atlas_size / (float)crunched_size);
    fprintf(stdout, "verifying decrunch ");

    uint8_t* decrunched = arena_alloc(&g_arena, default_font_atlas_size);
    bc_decrunch(crunched, crunched_size, 256, 256, bc4, decrunched);

    for(uint32_t i=0; i<default_font_atlas_size; ++i)
        if (decrunched[i] != default_font_atlas[i])
        {
            fprintf(stdout, "failed, bytes %u is different\n", i);
            return -1;
        }

    fprintf(stdout, "ok\n\n");

    global_zlib_ratio = 0.f;
    global_ratio = 0.f;
    num_ratios = 0;

    if (!test_bc4("../textures/bc4/grey_roof_tiles_02_ao_1k.png"))
        return -1;

    if (!test_bc4("../textures/bc4/grey_roof_tiles_02_disp_1k.png"))
        return -1;

    if (!test_bc4("../textures/bc4/patterned_cobblestone_ao_1k.png"))
        return -1;

    if (!test_bc4("../textures/bc4/patterned_cobblestone_disp_1k.png"))
        return -1;

    if (!test_bc4("../textures/bc4/red_brick_ao_1k.png"))
        return -1;

    if (!test_bc4("../textures/bc4/red_brick_disp_1k.png"))
        return -1;

    if (!test_bc4("../textures/bc4/rough_wood_ao_1k.png"))
        return -1;

    if (!test_bc4("../textures/bc4/rough_wood_disp_1k.png"))
        return -1;

    average_ratio = global_ratio / (float) num_ratios;
    fprintf(stdout, "\n\nBC4 average compression ratio : %f vs zlib ratio : %f\n\n", average_ratio, global_zlib_ratio / (float) num_ratios);

    fprintf(stdout, "\n\n-----------------------------------\n");
    fprintf(stdout, "| BC3 tests                       |\n");
    fprintf(stdout, "-----------------------------------\n\n");

    if (!test_bc3("../textures/bc3/apps-internet-web-browser.png"))
        return -1;

    if (!test_bc3("../textures/bc3/plant_05.png"))
        return -1;

    if (!test_bc3("../textures/bc3/plant_60.png"))
        return -1;

    fprintf(stdout, "\n\n-----------------------------------\n");
    fprintf(stdout, "| BC5 tests                       |\n");
    fprintf(stdout, "-----------------------------------\n\n");

    global_ratio = 0.f;
    num_ratios = 0;

    if (!test_bc5("../textures/bc5/grey_roof_tiles_02_nor_dx_1k.png"))
        return -1;

    if (!test_bc5("../textures/bc5/patterned_cobblestone_nor_dx_1k.png"))
        return -1;

    if (!test_bc5("../textures/bc5/red_brick_nor_dx_1k.png"))
        return -1;

    if (!test_bc5("../textures/bc5/rough_wood_nor_dx_1k.png"))
        return -1;

    average_ratio = global_ratio / (float) num_ratios;
    fprintf(stdout, "\n\nBC5 average compression ratio : %f\n\n", average_ratio);

    size_t bytes_allocated, byte_used;
    arena_stats(&g_arena, &bytes_allocated, &byte_used);
    fprintf(stdout, "\n%zu kb allocated\n%zu kb used\n\n", bytes_allocated>>10, byte_used>>10);

    arena_free(&g_arena);
    free(cruncher_memory);

    return 0;
}
