/*

zlib License

(C) 2025 Geolm

This software is provided 'as-is', without any express or implied
warranty.  In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.


Also as part of the code is based on FastAC, here's a copy of the license:

The only purpose of this program is to demonstrate the basic principles   
of arithmetic coding. It is provided as is, without any express or        
implied warranty, without even the warranty of fitness for any particular 
purpose, or that the implementations are correct.                         
                                                                          
Permission to copy and redistribute this code is hereby granted, provided 
that this warning and copyright notices are not removed or altered.       
                                                                          
Copyright (c) 2004 by Amir Said (said@ieee.org) &                         
                      William A. Pearlman (pearlw@ecse.rpi.edu)   

*/

#include <assert.h>
#include "bc_crunch.h"
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include <stdio.h>

#if defined(__aarch64__) || defined(_M_ARM64) || defined(__ARM_NEON)
    #include <arm_neon.h>
    #define BC_CRUNCH_NEON
    #define BC_TARGET_SIMD 
#elif defined(__x86_64__) || defined(_M_X64)
    #include <immintrin.h>
    #define BC_CRUNCH_SSE

    #if defined(_MSC_VER)
        #define BC_TARGET_SIMD
    #else
        #define BC_TARGET_SIMD __attribute__((target("sse4.1")))
    #endif
#endif

// enable this macro to increase BC1 compression ratio by about 1% but compression is slower (50%)
// decompression speed and output is not affected by this macro
#define BC_CRUNCH_USE_VECTOR_QUANTIZATION

//----------------------------------------------------------------------------------------------------------------------------
// Private structures & functions
//----------------------------------------------------------------------------------------------------------------------------

#define RC__MinLength (0x01000000U)
#define RC__MaxLength (0xFFFFFFFFU)
#define DM_MAX_SYMBOLS (256)
#define DM_MAX_TABLE_BITS 6
#define DM_MAX_TABLE_SIZE (1 << DM_MAX_TABLE_BITS)
#define DM__LengthShift (15)
#define DM__MaxCount    (1 << DM__LengthShift)
#define HASHMAP_SIZE (1 << 20)
#define TABLE_SIZE (128)
#define COLOR_DELTA_NUM_BITS (7)
#define DICTIONARY_SIZE (256)
#define MAKE48(r0, r1, r2) ( (((uint64_t)(r0) << 32) | ((uint64_t)(r1) << 16) | (uint64_t)(r2)) & 0x0000FFFFFFFFFFFFULL )
#define BC4_COLOR_NUM_BITS (8)
#define BC4_INDEX_NUM_BITS (3)
#define BC1_DIFFTABLE_SIZE (16)

static const uint32_t block_zigzag[16] = {0,  1,  2,  3, 7,  6,  5,  4, 8,  9, 10, 11, 15, 14, 13, 12};

#if defined(_MSC_VER)
    #include <intrin.h>
    #pragma intrinsic(__popcnt)
    #pragma intrinsic(__popcnt64)
    #define popcount64(x) ((int)__popcnt64(x))
    #define popcount(x) ((int)__popcnt(x))
#else
    #define popcount64(x) __builtin_popcountll(x)
    #define popcount(x) __builtin_popcount(x)
#endif

//----------------------------------------------------------------------------------------------------------------------
static inline int int_abs(int a) {return (a>=0) ? a : -a;}

//----------------------------------------------------------------------------------------------------------------------
typedef struct bytestream
{
    uint8_t* base;
    uint8_t* ptr;
    uint8_t* end;
} bytestream;

//----------------------------------------------------------------------------------------------------------------------
static inline void bs_init(bytestream* bs, void* buffer, size_t size)
{
    bs->base = (uint8_t*)buffer;
    bs->ptr  = bs->base;
    bs->end  = bs->base + size;
}

//----------------------------------------------------------------------------------------------------------------------
static inline uint32_t bs_remaining(const bytestream* bs)
{
    return (uint32_t)(bs->end - bs->ptr);
}

//----------------------------------------------------------------------------------------------------------------------
static inline uint32_t bs_offset(const bytestream* bs)
{
    return (uint32_t)(bs->ptr - bs->base);
}

//----------------------------------------------------------------------------------------------------------------------
static inline void bs_write_u8(bytestream* bs, uint8_t v)
{
    if (bs->ptr >= bs->end)
        return;

    *bs->ptr++ = v;
}

//----------------------------------------------------------------------------------------------------------------------
static inline void bs_write_buffer(bytestream* bs, const bytestream* other)
{
    uint32_t offset = bs_offset(other);
    if (bs->ptr + offset >= bs->end)
        return;

    memcpy(bs->ptr, other->base, offset);
    bs->ptr += offset;
}

//----------------------------------------------------------------------------------------------------------------------
static inline uint8_t bs_read_u8(bytestream* bs)
{
    if (bs->ptr >= bs->end)
        return 0;

    return *bs->ptr++;
}

//----------------------------------------------------------------------------------------------------------------------
static inline void bs_rewind(bytestream* bs)
{
    bs->ptr = bs->base;
}



// Maps signed:  0, -1,  1, -2,  2
// To unsigned:  0,  1,  2,  3,  4
static inline uint32_t zigzag_encode(int32_t n)
{
    return (uint32_t)((n << 1) ^ (n >> 31));
}

// Maps unsigned: 0,  1,  2,  3,  4
// To signed:    0, -1,  1, -2,  2
static inline int32_t zigzag_decode(uint32_t n)
{
    return (int32_t)((n >> 1) ^ (-(int32_t)(n & 1)));
}


//----------------------------------------------------------------------------------------------------------------------------
typedef struct bc1_block
{
    uint16_t color[2];
    uint32_t indices;
} bc1_block;

//----------------------------------------------------------------------------------------------------------------------------
typedef struct bc4_block
{
    uint8_t color[2];
    uint16_t indices[3];
} bc4_block;

//----------------------------------------------------------------------------------------------------------------------------
typedef struct entry
{
    uint32_t key;
    uint32_t count;
} entry;

//----------------------------------------------------------------------------------------------------------------------------
static inline void* ptr_shift(void *base, size_t shift_bytes)
{
    return (void *)((uint8_t *)base + shift_bytes);
}

//----------------------------------------------------------------------------------------------------------------------------
static inline const void* ptr_shift_const(const void *base, size_t shift_bytes)
{
    return (const void *)((const uint8_t *)base + shift_bytes);
}

//----------------------------------------------------------------------------------------------------------------------------
static inline const void* get_block(const void *base, size_t elem_size, uint32_t width_blocks, uint32_t x, uint32_t y)
{
    size_t row_index = (size_t)y * width_blocks;
    size_t idx = row_index + x;
    return (const uint8_t *)base + idx * elem_size;
}

//----------------------------------------------------------------------------------------------------------------------------
static inline uint8_t delta_encode_wrap(uint8_t prev, uint8_t curr)
{
    return curr - prev; // automatically wraps modulo 256
}

//----------------------------------------------------------------------------------------------------------------------------
static inline uint8_t delta_decode_wrap(uint8_t prev, uint8_t delta_encoded)
{
    return prev + delta_encoded; // automatically wraps modulo 256
}


//----------------------------------------------------------------------------------------------------------------------------
static inline void bc1_extract_565(uint16_t color, uint8_t *r5, uint8_t *g6, uint8_t *b5)
{
    *r5 = (uint8_t)((color >> 11) & 0x1F);
    *g6 = (uint8_t)((color >> 5)  & 0x3F);
    *b5 = (uint8_t)(color & 0x1F);
}

//----------------------------------------------------------------------------------------------------------------------------
static inline uint16_t bc1_pack_565(uint8_t r5, uint8_t g6, uint8_t b5)
{
    return (uint16_t)(((uint16_t)r5 << 11) | ((uint16_t)g6 << 5) | (uint16_t)b5);
}

//----------------------------------------------------------------------------------------------------------------------------
static inline uint32_t hash32(uint32_t x)
{
    x ^= x >> 16;
    x *= 0x7feb352d;
    x ^= x >> 15;
    x *= 0x846ca68b;
    x ^= x >> 16;
    return x;
}

//----------------------------------------------------------------------------------------------------------------------------
BC_TARGET_SIMD
uint32_t nearest32(const uint32_t* table, uint32_t table_size, uint32_t bitfield)
{
    uint32_t scores[TABLE_SIZE];
    uint32_t i = 0;

#ifdef BC_CRUNCH_NEON
    uint32x4_t bf_vec = vdupq_n_u32(bitfield);

    for (; i + 3 < table_size; i += 4)
    {
        uint32x4_t dict_vec = vld1q_u32(&table[i]);
        uint32x4_t delta = veorq_u32(dict_vec, bf_vec);
        uint8x16_t delta_bytes = vreinterpretq_u8_u32(delta);
        uint8x16_t counts = vcntq_u8(delta_bytes);
        uint16x8_t sum16 = vpaddlq_u8(counts);
        uint32x4_t sum32 = vpaddlq_u16(sum16);

        vst1q_u32(&scores[i], sum32);
    }
#elif defined(BC_CRUNCH_SSE)
    const __m128i bf_vec   = _mm_set1_epi32(bitfield);
    const __m128i mask_low = _mm_set1_epi8(0x0F);
    const __m128i lookup   = _mm_setr_epi8(0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4);

    for (; i + 3 < table_size; i += 4) 
    {
        __m128i x = _mm_xor_si128(_mm_loadu_si128((const __m128i*)&table[i]), bf_vec);

        __m128i low  = _mm_and_si128(x, mask_low);
        __m128i high = _mm_and_si128(_mm_srli_epi32(x, 4), mask_low); 
        __m128i cnt  = _mm_add_epi8(_mm_shuffle_epi8(lookup, low), 
                                    _mm_shuffle_epi8(lookup, high));

        __m128i lo_words = _mm_and_si128(cnt, _mm_set1_epi16(0x00FF));
        __m128i hi_words = _mm_srli_epi16(cnt, 8);
        __m128i sums = _mm_add_epi16(lo_words, hi_words);
        __m128i final = _mm_madd_epi16(sums, _mm_set1_epi16(1));

        _mm_storeu_si128((__m128i*)&scores[i], final);
    }
#endif

    // tail (or whole array on x64)
    for (; i < table_size; ++i)
    {
        uint32_t score = popcount(table[i] ^ bitfield);
        if (score == 0) 
            return (0 << 16) | (i & 0xffff);
        scores[i] = score;
    }
    

    // find best
    uint32_t best_index = 0;
    uint32_t best_score = UINT32_MAX;
    for (uint32_t j = 0; j < table_size; ++j)
    {
        uint32_t score = scores[j];
        if (score < best_score || (score == best_score && table[j] > table[best_index]))
        {
            best_score = score;
            best_index = j;
        }
    }
    return ((best_score&0xffff)<<16) | (best_index&0xffff);
}

//----------------------------------------------------------------------------------------------------------------------------
void vq_top_table(const void* input, size_t stride, uint32_t num_blocks, uint32_t* output, uint32_t* num_entries)
{
    // static array of odd steps to avoid aliasing and branches
    // first iteration is always 1 to ensure 100% initial coverage.
    static const uint32_t steps[4] = { 1, 5, 11, 17 };

    uint32_t centroids[TABLE_SIZE];
    struct 
    {
        uint32_t bit_diff_count[32];
        uint32_t count;
    } clusters[TABLE_SIZE];

    // use the top table entries as candidate for the cluster
    for(uint32_t i=0; i<*num_entries; ++i)
        centroids[i] = output[i];

    // multiple iteration to stabilize cluster
    for(uint32_t iteration=0; iteration<4; ++iteration)
    {
        // clear cluster counters
        for(uint32_t i=0; i<*num_entries; ++i)
        {
            for(uint32_t j=0; j<32; ++j)
                clusters[i].bit_diff_count[j] = 0;

            clusters[i].count = 0;
        }

        const uint32_t sample_step = steps[iteration];

        uint32_t bucket = 0;
        const uint32_t threshold = 16;

        // find the best cluster for each block, jittering the first block
        for(uint32_t block_index=(iteration % sample_step); block_index<num_blocks;)
        {
            const bc1_block* b = (const bc1_block*) get_block(input, stride, 0, block_index, 0);

            uint32_t result = nearest32(centroids, *num_entries, b->indices);
            uint32_t score = (result >> 16);
            uint32_t best_entry = result & 0xffff;
            uint32_t diff = b->indices ^ centroids[best_entry];

            for(uint32_t i=0; i<32; ++i)
                clusters[best_entry].bit_diff_count[i] += (diff >> i) & 1;

            clusters[best_entry].count++;
            bucket += score;

            if (bucket >= threshold) 
            {
                // if error is high, we only move 1 block
                block_index++;
                bucket -= threshold; 
            } else 
            {
                // normal skip
                block_index += sample_step;
            }
        }

        // move centroid
        for(uint32_t i=0; i<*num_entries; ++i)
        {
            if (clusters[i].count>0)
            {
                for(uint32_t bit=0; bit<32; ++bit)
                    if (clusters[i].bit_diff_count[bit] > (clusters[i].count/2))
                        centroids[i] ^= (1u << bit);

            }
            
        }
    }

    uint32_t num_clusters = 0;

    // fill the table with centroid
    for(uint32_t i=0; i<*num_entries; ++i)
        if (clusters[i].count > 0)
            output[num_clusters++] = centroids[i];

    // reduce if needed the size of the table to the number of cluster with at least one block
    *num_entries = num_clusters;
}

//----------------------------------------------------------------------------------------------------------------------------
int compare_entries(const void* a, const void* b)
{
    uint32_t entry_a = *(const uint32_t*) a;
    uint32_t entry_b = *(const uint32_t*) b;

    if (entry_a < entry_b)
        return -1;

    if (entry_a > entry_b)
        return 1;

    return 0;
}

//----------------------------------------------------------------------------------------------------------------------------
void build_top_table(entry* hashmap, const void* input, size_t stride, uint32_t num_blocks, uint32_t* output, uint32_t* num_entries)
{
    // clear the hashmap
    for(uint32_t i=0; i<HASHMAP_SIZE; ++i)
        hashmap[i].count = 0;

    // insert all blocks indices in the hashmap
    for(uint32_t i=0; i<num_blocks; ++i)
    {
        const bc1_block* b = (const bc1_block*) get_block(input, stride, 0, i, 0);

        uint32_t h = hash32(b->indices);
        uint32_t index = h & (HASHMAP_SIZE - 1);
        uint32_t first_index = index;

        bool inserted = false;
        while (!inserted)
        {
            if ((hashmap[index].count == 0) || (hashmap[index].key == b->indices))
            {
                hashmap[index].key = b->indices;
                hashmap[index].count++;
                inserted = true;
            }
            else
            {
                index = (index + 1) & (HASHMAP_SIZE - 1);
                assert(index != first_index);
            }
        }
    }

    // clear the table
    entry table[TABLE_SIZE];
    for(uint32_t i=0; i<TABLE_SIZE; ++i)
        table[i].count = 0;

    // fill the table with top most used indices
    for (uint32_t i = 0; i<HASHMAP_SIZE; ++i)
    {
        uint32_t c = hashmap[i].count;
        if (c == 0) 
            continue;

        // if smaller than current min, skip
        if (c <= table[0].count)
            continue;

        // replace min
        table[0] = hashmap[i];

        // bubble new smallest to front
        for (uint32_t j = 1; j < TABLE_SIZE; j++)
        {
            if (table[j-1].count > table[j].count)
            {
                entry tmp = table[j - 1];
                table[j-1] = table[j];
                table[j] = tmp;
            }
            else break;
        }
    }
    
    // reverse the table for output and count
    *num_entries = 0;
    for(uint32_t i=0; i<TABLE_SIZE; ++i)
    {
        output[i] = table[TABLE_SIZE-i-1].key;
        if (table[TABLE_SIZE-i-1].count>0)
            (*num_entries)++;
    }

#ifdef BC_CRUNCH_USE_VECTOR_QUANTIZATION
    // vector quantization
    vq_top_table(input, stride, num_blocks, output, num_entries);
#endif

    // // sort table for compression
    // qsort(output, *num_entries, sizeof(uint32_t), compare_entries);
}

//----------------------------------------------------------------------------------------------------------------------------
uint32_t nearest48(const uint64_t* table, uint32_t table_size, uint64_t bitfield)
{
    uint32_t scores[TABLE_SIZE];
    uint32_t i = 0;

    const uint64_t MASK48 = 0x0000FFFFFFFFFFFFull;

#ifdef BC_CRUNCH_NEON
    uint64x2_t bf_vec   = vdupq_n_u64(bitfield & MASK48);
    uint64x2_t mask_vec = vdupq_n_u64(MASK48);

    for (; i + 1 < table_size; i += 2)
    {
        uint64x2_t dict  = vandq_u64(vld1q_u64(&table[i]), mask_vec);
        uint64x2_t delta = veorq_u64(dict, bf_vec);

        uint8x16_t cnt = vcntq_u8(vreinterpretq_u8_u64(delta));

        // split lanes correctly
        uint8x8_t cnt_lo = vget_low_u8(cnt);
        uint8x8_t cnt_hi = vget_high_u8(cnt);

        // horizontal add per 64-bit lane
        scores[i + 0] = vaddlv_u8(cnt_lo);
        scores[i + 1] = vaddlv_u8(cnt_hi);
    }
#endif

    // tail
    for (; i < table_size; ++i)
    {
        uint64_t delta = (table[i] ^ bitfield) & MASK48;
        scores[i] = (uint32_t)popcount64(delta);
    }

    // find best
    uint32_t best_index = 0;
    uint32_t best_score = UINT32_MAX;

    for (uint32_t j = 0; j < table_size; ++j)
    {
        uint32_t score = scores[j];
        if (score < best_score ||
           (score == best_score && table[j] > table[best_index]))
        {
            best_score = score;
            best_index = j;
        }
    }

    return ((best_score & 0xffff) << 16) | (best_index & 0xffff);
}

//----------------------------------------------------------------------------------------------------------------------------
static inline uint8_t bc4_get_index(const bc4_block* b, uint32_t pixel_index)
{
    uint32_t bit_offset = pixel_index * 3;
    uint64_t bits = ((uint64_t)b->indices[0]) | ((uint64_t)b->indices[1] << 16) | ((uint64_t)b->indices[2] << 32);
    uint8_t index = (uint8_t)((bits >> bit_offset) & 0x7);
    return index;
}

//----------------------------------------------------------------------------------------------------------------------------
static inline void bc4_set_index(bc4_block* b, uint32_t pixel_index, uint8_t data)
{
    uint32_t bit_offset = pixel_index * 3;
    uint32_t word_index = bit_offset >> 4;
    uint32_t bit_in_word = bit_offset & 0xF;

    uint16_t mask = 0x7 << bit_in_word;
    b->indices[word_index] = (b->indices[word_index] & ~mask) | ((data & 0x7) << bit_in_word);

    // if the 3 bits spill into the next word
    if (bit_in_word > 13)  // last 2 or 1 bits spill
    {
        uint16_t spill_bits = (data & 0x7) >> (16 - bit_in_word);
        b->indices[word_index + 1] = (b->indices[word_index + 1] & ~(0x7 >> (16 - bit_in_word))) | spill_bits;
    }
}

//----------------------------------------------------------------------------------------------------------------------------
// static inline range_model* bc4_select_model(const bc4_block* b, range_model* indices)
// {
//     int endpoints_delta = int_abs(b->color[0] - b->color[1]);

//     if (endpoints_delta < 8)
//         return &indices[0];
//     else if (endpoints_delta < 32)
//         return &indices[8];

//     return &indices[16];
// }

//----------------------------------------------------------------------------------------------------------------------------
void bc1_crunch(bytestream* bs, void* cruncher_memory, const void* input, size_t stride, uint32_t width, uint32_t height)
{
    assert((width%4 == 0) && (height%4 == 0));
    assert(((uintptr_t)cruncher_memory)%sizeof(uintptr_t) == 0);

    uint32_t height_blocks = height/4;
    uint32_t width_blocks = width/4;

    // build a histogram and select the TABLE_SIZE block indices which are most used
    entry* hashmap = (entry*) cruncher_memory;
    uint32_t top_table[TABLE_SIZE];
    uint32_t top_table_size;
    build_top_table(hashmap, input, stride, height_blocks*width_blocks, top_table, &top_table_size);

    bs_write_u8(bs, (uint8_t)(top_table_size-1));

    for(uint32_t i=0; i<top_table_size; ++i)
        for(uint32_t j=0; j<4; ++j)
            bs_write_u8(bs, (top_table[i] >> (j*8)) & 0xff);

    bc1_block previous = {0};

    uint32_t block_index = 0;
    uint8_t block_buffer[128];
    bytestream block_stream;
    uint8_t block_indices_raw = 0;
    uint8_t block_color_raw = 0;
    uint8_t diff_table[BC1_DIFFTABLE_SIZE] = {0x00,0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80,0x03,0x05,0x09,0x0A,0x0C,0x30,0xC0};

    uint32_t fail_to_compress = 0;

    bs_init(&block_stream, block_buffer, sizeof(block_buffer));

    for(uint32_t y = 0; y < height_blocks; ++y)
    {
        for(uint32_t x = 0; x < width_blocks; ++x)
        {
            // zig-zag pattern delta compression for colors
            uint32_t zigzag_x = (y&1) ? x : width_blocks - x - 1;
            const bc1_block* current = get_block(input, stride, width_blocks, zigzag_x, y);

            uint8_t dgreen[2], dblue[2], dred[2];
            uint8_t flag_up[2];
            
            for(uint32_t j=0; j<2; ++j)
            {
                uint8_t current_red, current_green, current_blue;
                uint8_t previous_red, previous_green, previous_blue;

                bc1_extract_565(current->color[j], &current_red, &current_green, &current_blue);
                bc1_extract_565(previous.color[j], &previous_red, &previous_green, &previous_blue);

                flag_up[j] = 0;
                if (y>0 && x!=0)
                {
                    const bc1_block* up = get_block(input, stride, width_blocks, zigzag_x, y-1);
                    uint8_t up_red, up_green, up_blue;
                    bc1_extract_565(up->color[j], &up_red, &up_green, &up_blue);

                    int previous_delta = int_abs(current_red - previous_red) + int_abs(current_green-previous_green) + int_abs(current_blue-previous_blue);
                    int up_delta = int_abs(current_red-up_red) + int_abs(current_green-up_green) + int_abs(current_blue-up_blue);

                    flag_up[j] = (up_delta < previous_delta) ? 1 : 0;

                    // overwrite previous value to avoid using a new set of variables
                    if (up_delta < previous_delta)
                    {
                        previous_red = up_red;
                        previous_green = up_green;
                        previous_blue = up_blue;
                    }
                }

                int delta_green = current_green - previous_green;

                dred[j] =  zigzag_encode((int)current_red - (int)previous_red - delta_green/2);
                dgreen[j] =  zigzag_encode(delta_green);
                dblue[j] = zigzag_encode((int)current_blue - (int)previous_blue - delta_green/2);
            }

            bool color_escape = (dgreen[0] > 7 || dgreen[1] > 7 || dblue[0] > 3 || dblue[1] > 3 ||dred[0] > 3 || dred[1] > 3);
            block_color_raw |= (color_escape ? 1 : 0) << block_index;

            uint8_t color_byte0 = flag_up[0] << 7 | flag_up[1] << 6 | (dgreen[0]&0x7) << 3 | (dgreen[1]&0x7);

            if (color_escape)
            {
                // escape case, paying the full price
                bs_write_u8(&block_stream, (current->color[0] >> 8) & 0xff);
                bs_write_u8(&block_stream, current->color[0] & 0xff);
                bs_write_u8(&block_stream, (current->color[1] >> 8) & 0xff);
                bs_write_u8(&block_stream, current->color[1] & 0xff);
            }
            else
            {
                bs_write_u8(&block_stream, color_byte0);
                bs_write_u8(&block_stream, dred[0] << 6 | dred[1] << 4 | dblue[0] << 2 | dblue[1]);
            }

            uint32_t table_index = nearest32(top_table, top_table_size, current->indices) & 0xffff;
            uint32_t difference = current->indices ^ top_table[table_index];
            uint32_t diff_indices[4];

            // try to encode the diff with the difftable
            for(uint32_t j=0; j<4; ++j)
            {
                uint8_t value = (difference>>(j*8)) & 0xff;
                diff_indices[j] = UINT32_MAX;

                for(uint32_t i=0; i<BC1_DIFFTABLE_SIZE && (diff_indices[j] == UINT32_MAX); ++i)
                    if (value == diff_table[i])
                        diff_indices[j] = i;
            }

            bool indices_escape = (diff_indices[0] > 15) || (diff_indices[1] > 15) || 
                               (diff_indices[2] > 15) || (diff_indices[3] > 15);

            block_indices_raw |= (indices_escape ? 1 : 0) << block_index;

            if (indices_escape)
            {
                fail_to_compress++;
                for(uint32_t j=0; j<4; ++j)
                    bs_write_u8(&block_stream, (current->indices >> (j*8)) & 0xff);
            }
            else
            {
                uint8_t mode = ((diff_indices[0] < 4) && (diff_indices[1] < 4) && 
                                (diff_indices[2] < 4) && (diff_indices[3] < 4)) ? 1 : 0;

                bs_write_u8(&block_stream, (mode<<7) | table_index);

                if (mode == 1)
                {
                    bs_write_u8(&block_stream, (diff_indices[0]<<6) | (diff_indices[1] << 4) |
                                               (diff_indices[2]<<2) | diff_indices[3]);
                }
                else
                {
                    bs_write_u8(&block_stream, (diff_indices[0]<<4) | diff_indices[1]);
                    bs_write_u8(&block_stream, (diff_indices[2]<<4) | diff_indices[3]);
                }
            }

            // flush 8 blocks
            if (++block_index==8)
            {
                bs_write_u8(bs, block_color_raw);
                bs_write_u8(bs, block_indices_raw);
                bs_write_buffer(bs, &block_stream);
                bs_rewind(&block_stream);
                block_index = 0;
                block_indices_raw = 0;
                block_color_raw = 0;
            }

            previous = *current;
        }
    }

    bs_write_buffer(bs, &block_stream);

    float num_blocks = height_blocks * width_blocks;
    printf("fail_to_compress ratio : %f\n", (float) fail_to_compress / num_blocks);
}

//----------------------------------------------------------------------------------------------------------------------------
void bc1_decrunch(bytestream* bs, uint32_t width, uint32_t height, void* output, size_t stride)
{
    assert((width % 4 == 0) && (height % 4 == 0));

    uint32_t height_blocks = height/4;
    uint32_t width_blocks = width/4;

    uint32_t top_table[TABLE_SIZE];
    uint32_t top_table_size = bs_read_u8(bs) + 1;

    for(uint32_t i=0; i<top_table_size; ++i)
    {
        top_table[i] = 0;
        for(uint32_t j=0; j<4; ++j)
            top_table[i] |= bs_read_u8(bs) << (j*8);
    }

    bc1_block previous = {0};

    uint32_t block_index = 0;
    uint8_t block_indices_raw = 0;
    uint8_t block_color_raw = 0;
    uint8_t diff_table[BC1_DIFFTABLE_SIZE] = {0x00,0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80,0x03,0x05,0x09,0x0A,0x0C,0x30,0xC0};

    for(uint32_t y = 0; y < height_blocks; ++y)
    {
        for(uint32_t x = 0; x < width_blocks; ++x)
        {
            if (block_index == 0)
            {
                block_color_raw = bs_read_u8(bs);
                block_indices_raw = bs_read_u8(bs);
            }

            // zig-zag pattern color
            uint32_t zigzag_x = (y&1) ? x : width_blocks - x - 1;
            bc1_block* current = (bc1_block*) get_block(output, stride, width_blocks, zigzag_x, y);

            // escape code, full r5g6b5 color encoded :(
            if (block_color_raw & (1<<block_index))
            {
                current->color[0] = bs_read_u8(bs) << 8;
                current->color[0]|= bs_read_u8(bs);
                current->color[1] = bs_read_u8(bs) << 8;
                current->color[1]|= bs_read_u8(bs);
            }
            else
            {
                uint8_t color_byte0 = bs_read_u8(bs);
                uint8_t flag_up[2];
                flag_up[0] = (color_byte0>>7) & 0x1;
                flag_up[1] = (color_byte0>>6) & 0x1;
                uint8_t color_byte1 = bs_read_u8(bs);

                int dgreen[2], dred[2], dblue[2];

                dgreen[0] = zigzag_decode((color_byte0>>3)&0x7);
                dgreen[1] = zigzag_decode(color_byte0&0x7);
                dred[0] = zigzag_decode((color_byte1>>6)&0x3) + dgreen[0] / 2;
                dred[1] = zigzag_decode((color_byte1>>4)&0x3) + dgreen[1] / 2;
                dblue[0] = zigzag_decode((color_byte1>>2)&0x3) + dgreen[0] / 2;
                dblue[1] = zigzag_decode(color_byte1&0x3) + dgreen[1] / 2;

                for (uint32_t j = 0; j < 2; ++j)
                {
                    uint8_t reference_red, reference_green, reference_blue;
                    bc1_extract_565(previous.color[j], &reference_red, &reference_green, &reference_blue);
                    if (flag_up[j])
                    {
                        bc1_block* up = (bc1_block*) get_block(output, stride, width_blocks, zigzag_x, y-1);
                        bc1_extract_565(up->color[j], &reference_red, &reference_green, &reference_blue);
                    }

                    uint8_t current_green_value = (int)reference_green + dgreen[j];
                    uint8_t current_red_value = (int)reference_red + dred[j];
                    uint8_t current_blue_value = (int)reference_blue + dblue[j];

                    current->color[j] = bc1_pack_565(current_red_value, current_green_value, current_blue_value);
                }
            }

            uint32_t difference = 0; 
            uint32_t diff_indices[4];

            if (block_indices_raw & (1<<block_index))
            {
                current->indices = 0;
                for(uint32_t j=0; j<4; ++j)
                    current->indices |= ((uint32_t)bs_read_u8(bs)) << (j*8);

                difference = previous.indices ^ current->indices;
            }
            else
            {
                uint8_t data = bs_read_u8(bs);
                uint8_t table_index = data&0x7f;
                uint8_t mode = (data>>7)&1;

                if (mode == 0)
                {
                    data = bs_read_u8(bs);
                    diff_indices[0] = data >> 4;
                    diff_indices[1] = data & 0xf;

                    data = bs_read_u8(bs);
                    diff_indices[2] = data >> 4;
                    diff_indices[3] = data & 0xf;
                }
                else
                {
                    data = bs_read_u8(bs);
                    diff_indices[0] = (data >> 6) & 0x3;
                    diff_indices[1] = (data >> 4) & 0x3;
                    diff_indices[2] = (data >> 2) & 0x3;
                    diff_indices[3] = (data >> 0) & 0x3;
                }

                for(uint32_t j=0; j<4; ++j)
                    difference |= diff_table[diff_indices[j]] << (j*8);

                current->indices = top_table[table_index] ^ difference;
            }

            if (++block_index == 8)
            {
                block_index = 0;
                block_indices_raw = 0;
                block_color_raw = 0;
            }

            previous = *current;
        }
    }
}

/*

//----------------------------------------------------------------------------------------------------------------------------
void bc4_crunch(range_codec* codec, void* cruncher_memory, const void* input, size_t stride, uint32_t width, uint32_t height)
{
    assert((width%4 == 0) && (height%4 == 0));
    assert(((uintptr_t)cruncher_memory)%sizeof(uintptr_t) == 0);

    uint32_t height_blocks = height/4;
    uint32_t width_blocks = width/4;

    range_model color_delta[2];
    model_init(&color_delta[0], 1<<BC4_COLOR_NUM_BITS);
    model_init(&color_delta[1], 1<<BC4_COLOR_NUM_BITS);

    range_model color_reference, first_index, use_dict, dict_reference; 
    model_init(&color_reference, 2);
    model_init(&first_index, 1<<3);
    model_init(&use_dict, 2);
    model_init(&dict_reference, DICTIONARY_SIZE);

    range_model indices[24];
    for(uint32_t i=0; i<24; ++i)
        model_init(&indices[i], 1<<BC4_INDEX_NUM_BITS);

    range_model dict_delta[16];
    for(uint32_t i=0; i<16; ++i)
        model_init(&dict_delta[i], 1<<3);

    bc4_block empty_block = {.color = {0, 128}};
    const bc4_block* previous = &empty_block;

    // dictionary initialization
    uint64_t dictionary[DICTIONARY_SIZE];
    for(uint32_t i=0; i<DICTIONARY_SIZE; ++i)
        dictionary[i] = UINT64_MAX;

    for(uint32_t y = 0; y < height_blocks; ++y)
    {
        for(uint32_t x = 0; x < width_blocks; ++x)
        {
            const bc4_block* current = get_block(input, stride, width_blocks, x, y);

            int reference = previous->color[0];
            if (y>0)
            {
                const bc4_block* up = get_block(input, stride, width_blocks, x, y-1);
                if (x>0)
                {
                    const bc4_block* up_left = get_block(input, stride, width_blocks, x-1, y-1);
                    reference += up->color[0] - up_left->color[0];
                }
                else
                    reference = up->color[0];
            }

            if (reference < 0) reference = 0;
            if (reference > 255) reference = 255;

            enc_put(codec, &color_delta[0], delta_encode_wrap((uint8_t)reference, current->color[0]));
            enc_put(codec, &color_delta[1], delta_encode_wrap(current->color[0], current->color[1]));

            // search in the dictionary for the current bitfield
            uint64_t bitfield = MAKE48(current->indices[0], current->indices[1], current->indices[2]);
            uint32_t dict_lookup = nearest48(dictionary, DICTIONARY_SIZE, bitfield);
            uint16_t score = dict_lookup>>16;
            uint16_t found_index = dict_lookup&0xffff;
            
            // found or similar? just write the dictionary index
            if (score < 4)
            {
                enc_put(codec, &use_dict, 1);
                enc_put(codec, &dict_reference, found_index);
                
                uint64_t reference = dictionary[found_index];
                uint64_t bitfield_delta = reference ^ bitfield;
                for(uint32_t j=0; j<16; ++j)
                    enc_put(codec, &dict_delta[j], (bitfield_delta>>(j*3))&0x7);

                if(found_index > 0)
                {
                    uint64_t temp = dictionary[found_index];
                    memmove(&dictionary[1], &dictionary[0], found_index * sizeof(uint64_t));
                    dictionary[0] = temp;
                }
            }
            else
            {
                // store the entry in the dictionary
                memmove(&dictionary[1], &dictionary[0], (DICTIONARY_SIZE - 1) * sizeof(uint64_t));
                dictionary[0] = bitfield;

                // write the indices with local difference delta encoded
                enc_put(codec, &use_dict, 0);

                uint8_t block_previous = bc4_get_index(current, 0);
                enc_put(codec, &first_index, block_previous);

                range_model* model = bc4_select_model(current, indices);
                for(uint32_t j=1; j<16; ++j)
                {
                    uint8_t data = bc4_get_index(current, block_zigzag[j]);
                    enc_put(codec, &model[block_previous], block_previous ^ data);
                    block_previous = data;
                }
            }
            previous = current;
        }
    }
}

//----------------------------------------------------------------------------------------------------------------------------
void bc4_decrunch(range_codec* codec, uint32_t width, uint32_t height, void* output, size_t stride)
{
    assert((width % 4 == 0) && (height % 4 == 0));

    uint32_t height_blocks = height/4;
    uint32_t width_blocks = width/4;

    range_model color_delta[2];
    model_init(&color_delta[0], 1<<BC4_COLOR_NUM_BITS);
    model_init(&color_delta[1], 1<<BC4_COLOR_NUM_BITS);

    range_model color_reference, first_index, use_dict, dict_reference;
    model_init(&color_reference, 2);
    model_init(&first_index, 1<<3);
    model_init(&use_dict, 2);
    model_init(&dict_reference, DICTIONARY_SIZE);

    range_model indices[24];
    for(uint32_t i=0; i<24; ++i)
        model_init(&indices[i], 1<<BC4_INDEX_NUM_BITS);

    range_model dict_delta[16];
    for(uint32_t i=0; i<16; ++i)
        model_init(&dict_delta[i], 1<<3);

    bc4_block empty_block = {.color = {0, 128}};
    bc4_block* previous = &empty_block;

    // dictionary initialization
    uint64_t dictionary[DICTIONARY_SIZE];
    for(uint32_t i=0; i<DICTIONARY_SIZE; ++i)
        dictionary[i] = UINT64_MAX;

    for(uint32_t y = 0; y < height_blocks; ++y)
    {
        for(uint32_t x = 0; x < width_blocks; ++x)
        {
            bc4_block* current = (bc4_block*) get_block(output, stride, width_blocks, x, y);
            int reference = previous->color[0];
            if (y>0)
            {
                const bc4_block* up = get_block(output, stride, width_blocks, x, y-1);
                if (x>0)
                {
                    const bc4_block* up_left = get_block(output, stride, width_blocks, x-1, y-1);
                    reference += up->color[0] - up_left->color[0];
                }
                else
                    reference = up->color[0];
            }

            if (reference < 0) reference = 0;
            if (reference > 255) reference = 255;

            current->color[0] = delta_decode_wrap((uint8_t)reference, dec_get(codec, &color_delta[0]));
            current->color[1] = delta_decode_wrap(current->color[0], dec_get(codec, &color_delta[1]));

            if (dec_get(codec, &use_dict))
            {
                // data should be in the dictionary
                uint32_t found_index = dec_get(codec, &dict_reference);
                uint64_t reference = dictionary[found_index];
                uint64_t bitfield = 0;
    
                for(uint32_t j=0; j<16; ++j)
                {
                    uint64_t byte = (uint64_t)dec_get(codec, &dict_delta[j]);
                    bitfield |= (byte << (j*3));
                }

                bitfield ^= reference;
                current->indices[0] = ((bitfield>>32) & 0xffff);
                current->indices[1] = ((bitfield>>16) & 0xffff);
                current->indices[2] = bitfield & 0xffff;

                if(found_index > 0)
                {
                    uint64_t temp = dictionary[found_index];
                    memmove(&dictionary[1], &dictionary[0], found_index * sizeof(uint64_t));
                    dictionary[0] = temp;
                }
            }
            else
            {
                uint8_t block_previous = dec_get(codec, &first_index);
                bc4_set_index(current, 0, block_previous);

                range_model* model = bc4_select_model(current, indices);
                for(uint32_t j=1; j<16; ++j)
                {
                    uint8_t delta = dec_get(codec, &model[block_previous]);
                    uint8_t data = block_previous ^ delta;
                    bc4_set_index(current, block_zigzag[j], data);
                    block_previous = data;
                }

                // store the entry in the dictionary
                memmove(&dictionary[1], &dictionary[0], (DICTIONARY_SIZE - 1) * sizeof(uint64_t));
                dictionary[0] = MAKE48(current->indices[0], current->indices[1], current->indices[2]);
            }
            previous = current;
        }
    }
}

*/

/*
BC4 average compression ratio history

base : 1.051032
zig-zag : 1.090952
up-left-xor : 1.196489
xor endpoints : 1.124062
just left or up : 1.211573
block zig-zag xor : 1.211644
xor-delta dictionary :  : 1.227614
move-to-front dictionary : 1.232521
16x 3bits model for dict xor 1.234481
contextual model for dictionary miss : 1.298764
use endpoints range for context (<32) : 1.301263
endpoints range (<16) :  1.304174
endpoints range (<8) : 1.312253
multiple buckets : 1.313912
second endpoint encoding from first one : 1.3175
removed zig-zag, use left+up-up_left : 1.324589
removed circular dictionnary 1.327702
fixed x=0 block: 1.327747
score < 4 : 1.336385
*/


//----------------------------------------------------------------------------------------------------------------------------
// Public functions
//----------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------
size_t crunch_min_size(void)
{
    return sizeof(entry) * HASHMAP_SIZE;
}

//----------------------------------------------------------------------------------------------------------------------------
size_t bc_crunch(void* cruncher_memory, const void* input, uint32_t width, uint32_t height, enum bc_format format, void* output, size_t length)
{
    assert(cruncher_memory != NULL && "bc_crunch needs memory to run, allocate a buffer of size crunch_min_size()");
    assert(input != NULL);

    bytestream stream;
    bs_init(&stream, output, length);

    switch(format)
    {
    case bc1 : 
        {
            bc1_crunch(&stream, cruncher_memory, input, sizeof(bc1_block), width, height);
            break;
        }
    // case bc3 : 
    //     {
    //         size_t block_size = sizeof(bc1_block) + sizeof(bc4_block);
    //         bc4_crunch(&codec, cruncher_memory, input, block_size, width, height);
    //         bc1_crunch(&codec, cruncher_memory, ptr_shift_const(input, sizeof(bc4_block)), block_size, width, height);
    //         break;
    //     }
    // case bc4 : 
    //     {
    //         bc4_crunch(&codec, cruncher_memory, input, sizeof(bc4_block), width, height);
    //         break;
    //     }
    // case bc5 :
    //     {
    //         size_t block_size = sizeof(bc4_block) * 2;
    //         bc4_crunch(&codec, cruncher_memory, input, block_size, width, height);
    //         bc4_crunch(&codec, cruncher_memory, ptr_shift_const(input, sizeof(bc4_block)), block_size, width, height);
    //         break;
    //     }

    default: break;
    }
    return bs_offset(&stream);
}

//----------------------------------------------------------------------------------------------------------------------------
void bc_decrunch(const void* input, size_t length, uint32_t width, uint32_t height, enum bc_format format, void* output)
{
    bytestream stream;
    bs_init(&stream, (void*)input, length);

    switch(format)
    {
    case bc1 : bc1_decrunch(&stream, width, height, output, sizeof(bc1_block)); break;
    // case bc3 : 
    //     {
    //         size_t block_size = sizeof(bc1_block) + sizeof(bc4_block);
    //         bc4_decrunch(&codec, width, height, output, block_size);
    //         bc1_decrunch(&codec, width, height, ptr_shift(output, sizeof(bc4_block)), block_size);
    //         break;
    //     }
    // case bc4 : bc4_decrunch(&codec, width, height, output, sizeof(bc4_block)); break;
    // case bc5 :
    //     {
    //         size_t block_size = sizeof(bc4_block) * 2;
    //         bc4_decrunch(&codec, width, height, output, block_size);
    //         bc4_decrunch(&codec, width, height, ptr_shift(output, sizeof(bc4_block)), block_size);
    //         break;
    //     }
    default: break;
    }
}