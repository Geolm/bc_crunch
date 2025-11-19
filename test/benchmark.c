#include "../bc_crunch.h"

#define SOKOL_TIME_IMPL
#include "sokol_time.h"

#include <string.h>
#define STB_DXT_IMPLEMENTATION
#include "stb_dxt.h"

#define STBI_ONLY_PNG
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"


#define NUM_ITERATIONS (100U)

void* cruncher_memory = NULL;

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
void profile_bc1(const char* filename)
{
    uint32_t width, height, num_channels;
    unsigned char *rgba = stbi_load(filename, (int*)&width, (int*)&height, (int*)&num_channels, 4);

    if (rgba)
    {
        size_t bc1_size = (width * height) / 2;
        uint8_t* bc1_texture = malloc(bc1_size);

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

        size_t worst_case = bc1_size*2;
        void* crunched_texture = malloc(worst_case);

        uint64_t start = stm_now();

        for(uint32_t i=0; i<NUM_ITERATIONS; ++i)
            bc_crunch(cruncher_memory, bc1_texture, width, height, bc1, crunched_texture, worst_case);

        double duration = stm_sec(stm_now() - start);

        printf("crunch %u times a %ux%u texture in %f seconds\n", NUM_ITERATIONS, width, height, duration);
        printf("output : %f MB/s\n\n", ((double)bc1_size / 1048576.0) * (double) NUM_ITERATIONS / duration);

        start = stm_now();

        for(uint32_t i=0; i<NUM_ITERATIONS; ++i)
            bc_decrunch(crunched_texture, worst_case, width, height, bc1, bc1_texture);

        duration = stm_sec(stm_now() - start);

        printf("decrunch %u times a %ux%u texture in %f seconds\n", NUM_ITERATIONS, width, height, duration);
        printf("output : %f MB/s\n\n", ((double)bc1_size / 1048576.0) * (double) NUM_ITERATIONS / duration);

        free(crunched_texture);
        free(bc1_texture);
        stbi_image_free(rgba);
    }
}


int main(void)
{
    printf("\nbc_crunch profiler\n\n");
    stm_setup();
    cruncher_memory = malloc(crunch_min_size());

    profile_bc1("../textures/bc1/roof_tiles.png");


    printf("done\n\n");

    free(cruncher_memory);
    return 0;
}


/* 

Optimizations history (all on my M1 Pro machine)

Base version
    crunch 50 times a 1024x1024 texture in 1.596338 seconds
    output : 15.660845 MB/s

    decrunch 50 times a 1024x1024 texture in 1.096381 seconds
    output : 22.802282 MB/s

Range encoder decoder table
    crunch 50 times a 1024x1024 texture in 1.601962 seconds
    output : 15.605862 MB/s

    decrunch 50 times a 1024x1024 texture in 0.723421 seconds
    output : 34.558022 MB/s

Neon SIMD dictionary_nearest

    crunch 50 times a 1024x1024 texture in 1.252861 seconds
    output : 19.954335 MB/s

    decrunch 50 times a 1024x1024 texture in 0.720878 seconds
    output : 34.679910 MB/s

Neon SIMD popcount 32 bits for BC1

    crunch 100 times a 1024x1024 texture in 2.233203 seconds
    output : 22.389362 MB/s

    decrunch 100 times a 1024x1024 texture in 1.448726 seconds
    output : 34.513074 MB/s

*/