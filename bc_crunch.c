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

#define HASHMAP_SIZE (1 << 20)
#define TABLE_INDEX_NUM_BITS (8)
#define TABLE_SIZE (1<<TABLE_INDEX_NUM_BITS)
#define COLOR_DELTA_NUM_BITS (7)
#define DICTIONARY_SIZE (256)
#define MAKE48(r0, r1, r2) ( (((uint64_t)(r0) << 32) | ((uint64_t)(r1) << 16) | (uint64_t)(r2)) & 0x0000FFFFFFFFFFFFULL )
#define BC4_COLOR_NUM_BITS (8)
#define BC4_INDEX_NUM_BITS (3)

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


#define LE_MODEL_MAX_HOT (36)
#define LE_HISTOGRAM_SIZE (256)

enum le_mode
{
    le_mode_idle,
    le_mode_encode,
    le_mode_decode
};

typedef struct le_stream
{
    uint8_t* buffer;
    uint32_t bit_offset;
    size_t position;
    size_t size;

    uint64_t bit_reservoir;
    uint32_t bits_available;

    enum le_mode mode;
} le_stream;

typedef struct le_model
{
    uint8_t hot_values[LE_MODEL_MAX_HOT];
    uint8_t cold_min;
    uint8_t cold_max;
    uint8_t cold_num_bits;

#ifdef LE_STATS
    uint32_t num_hot_tier0;
    uint32_t num_hot_tier1;
    uint32_t num_hot_tier2;
    uint32_t num_raw;
#endif
} le_model;

typedef struct le_histogram
{
    uint32_t count[LE_HISTOGRAM_SIZE];
    uint32_t num_symbols;
} le_histogram;


// Dibit-delta encoding model, specialized in small delta, use a good predictor up-front to maximize compression
//
// Use dibit to encode delta. 
//      (0, -1, 1)      use 2 bits
//      (-2, +2, -3)    use 4 bits
//      (+3, -4, +4)    use 6 bits
//      (-5, +5, -6)    use 8 bits
//      anything greater use 16 bits

// ----------------------------------------------------------------------------------------------------------------------------
static inline void le_refill(le_stream* s)
{
    // Pull bytes into the reservoir until it's full enough for any standard read
    while (s->bits_available <= 56 && s->position < s->size)
    {
        s->bit_reservoir |= ((uint64_t)s->buffer[s->position]) << s->bits_available;
        s->bits_available += 8;
        s->position++;
    }
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline void le_flush(le_stream* s)
{
    while (s->bits_available >= 8)
    {
#ifdef LE_CHECKS
        assert(s->position < s->size);
#endif
        s->buffer[s->position] = (uint8_t)(s->bit_reservoir & 0xFF);
        s->bit_reservoir >>= 8;
        s->bits_available -= 8;
        s->position++;
    }
}

static inline void le_init(le_stream *s, void* buffer, size_t size)
{
#ifdef LE_CHECKS
    assert(s != NULL);
    assert(buffer != NULL);
#endif
    s->buffer = (uint8_t*)buffer;
    s->size = size;
    s->position = 0;
    s->bit_reservoir = 0;
    s->bits_available = 0;
    s->mode = le_mode_idle;
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline void le_begin_encode(le_stream* s)
{
    s->position = 0;
    s->bit_reservoir = 0;
    s->bits_available = 0;
    s->mode = le_mode_encode;
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline size_t le_end_encode(le_stream* s)
{
    while (s->bits_available > 0)
    {
        if (s->position < s->size)
        {
            s->buffer[s->position] = (uint8_t)(s->bit_reservoir & 0xFF);
            s->bit_reservoir >>= 8;
            s->position++;
        }
        
        if (s->bits_available > 8)
        {
            s->bits_available -= 8;
        }
        else
        {
            s->bits_available = 0;
        }
    }

    s->mode = le_mode_idle;
    return s->position;
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline void le_begin_decode(le_stream* s)
{
    s->position = 0;
    s->bit_reservoir = 0;
    s->bits_available = 0;
    s->mode = le_mode_decode;
    le_refill(s);
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline void le_end_decode(le_stream* s)
{
    s->mode = le_mode_idle;
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline void le_write_dibit(le_stream* s, uint8_t dibit)
{
    s->bit_reservoir |= ((uint64_t)(dibit & 0x03) << s->bits_available);
    s->bits_available += 2;
    if (s->bits_available >= 32)
    {
        le_flush(s);
    }
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline void le_write_nibble(le_stream* s, uint8_t nibble)
{
    s->bit_reservoir |= ((uint64_t)(nibble & 0x0F) << s->bits_available);
    s->bits_available += 4;
    if (s->bits_available >= 32)
    {
        le_flush(s);
    }
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline void le_write_byte(le_stream* s, uint8_t value)
{
    s->bit_reservoir |= ((uint64_t)value << s->bits_available);
    s->bits_available += 8;
    if (s->bits_available >= 32)
    {
        le_flush(s);
    }
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline uint8_t le_read_dibit(le_stream* s)
{
    if (s->bits_available < 2)
    {
        le_refill(s);
    }

    uint8_t value = (uint8_t)(s->bit_reservoir & 0x03);
    s->bit_reservoir >>= 2;
    s->bits_available -= 2;
    return value;
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline uint8_t le_read_nibble(le_stream* s)
{
    if (s->bits_available < 4)
    {
        le_refill(s);
    }

    uint8_t value = (uint8_t)(s->bit_reservoir & 0x0F);
    s->bit_reservoir >>= 4;
    s->bits_available -= 4;
    return value;
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline uint8_t le_read_byte(le_stream* s)
{
    if (s->bits_available < 8)
    {
        le_refill(s);
    }

    uint8_t value = (uint8_t)(s->bit_reservoir & 0xFF);
    s->bit_reservoir >>= 8;
    s->bits_available -= 8;
    return value;
}

//----------------------------------------------------------------------------------------------------------------------------
static inline void histogram_init(le_histogram* h, uint32_t num_symbols)
{
#ifdef LE_CHECKS
    assert(num_symbols > 3 && num_symbols <= 256);
#endif

    h->num_symbols = num_symbols;
    for(uint32_t i=0; i<num_symbols; ++i)
        h->count[i] = 0;
}

#define CONTROL_HOT_TIER0 (0)
#define CONTROL_HOT_TIER1 (1)
#define CONTROL_HOT_TIER2 (2)
#define CONTROL_ESCAPE (3)

// ----------------------------------------------------------------------------------------------------------------------------
void le_model_init(le_model *model, const uint32_t *histogram, uint32_t num_symbols)
{
    memset(model->hot_values, 0, sizeof(model->hot_values));

#ifdef LE_STATS
    model->num_hot_tier0 = 0;
    model->num_hot_tier1 = 0;
    model->num_hot_tier2 = 0;
    model->num_raw = 0;
#endif

    uint32_t selected[LE_MODEL_MAX_HOT];
    for (uint32_t i = 0; i < LE_MODEL_MAX_HOT; i++)
        selected[i] = UINT32_MAX;

    uint32_t hot_used = 0;

    for (uint32_t i = 0; i < LE_MODEL_MAX_HOT; i++)
    {
        uint32_t max_count = 0;
        uint32_t max_index = UINT32_MAX;

        for (uint32_t s = 0; s < num_symbols; s++)
        {
            bool already = false;
            for (uint32_t j = 0; j < i; j++)
            {
                if (selected[j] == s) 
                { 
                    already = true; 
                    break; 
                }
            }

            if (already) 
                continue;

            if (histogram[s] > max_count)
            {
                max_count = histogram[s];
                max_index = s;
            }
        }

        if (max_index == UINT32_MAX || max_count == 0)
            break;

        model->hot_values[i] = (uint8_t)max_index;
        selected[i] = max_index;
        hot_used++;
    }

    model->cold_min = UINT8_MAX;
    model->cold_max = 0;
    for(uint32_t s = 0; s < num_symbols; s++)
    {
        // skip hot symbols
        bool is_hot = false;
        for (uint32_t i = 0; i < hot_used; i++)
        {
            if (model->hot_values[i] == s)
            {
                is_hot = true;
                break;
            }
        }

        if (is_hot || histogram[s] == 0)
            continue;

        if (s < model->cold_min)
            model->cold_min = (uint8_t)s;
        if (s > model->cold_max)
            model->cold_max = (uint8_t)s;
    }

    model->cold_num_bits = 0;
    if (model->cold_max >= model->cold_min)
    {
        uint8_t range = model->cold_max - model->cold_min;

        if (range >= 64)
            model->cold_num_bits = 8;
        else if (range >= 16)
            model->cold_num_bits = 6;
        else if (range >= 4)
            model->cold_num_bits = 4;
        else
            model->cold_num_bits = 2;
    }
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline void le_encode(le_stream *s, le_model *model, uint8_t value)
{
#ifdef LE_CHECKS
    assert(s->mode == le_mode_encode);
#endif

    // hot values
    for (uint32_t i = 0; i < LE_MODEL_MAX_HOT; i++)
    {
        if (model->hot_values[i] == value)
        {
            if (i<4)
            {
                le_write_dibit(s, CONTROL_HOT_TIER0);
                le_write_dibit(s, i);

            #ifdef LE_STATS
                model->num_hot_tier0++;
            #endif
            }
            else if (i<20)
            {
                le_write_dibit(s, CONTROL_HOT_TIER1);
                le_write_nibble(s, i-4);

            #ifdef LE_STATS
                model->num_hot_tier1++;
            #endif
            }
            else
            {
                le_write_dibit(s, CONTROL_HOT_TIER2);
                le_write_nibble(s, i-20);

            #ifdef LE_STATS
                model->num_hot_tier2++;
            #endif
            }
            return;
        }
    }

    // or escape code
    le_write_dibit(s, CONTROL_ESCAPE);
    le_write_byte(s, value);

#ifdef LE_STATS
    model->num_raw++;
#endif
}

static const uint32_t consumption_lut[4] = { 2, 4, 4, 8 };
static const uint32_t offset_lut[4]      = { 0, 4, 20, 0 };

// ----------------------------------------------------------------------------------------------------------------------------
uint8_t static inline le_decode(le_stream *s, le_model *model)
{
    if (s->bits_available < 16) le_refill(s);

    uint64_t res = s->bit_reservoir;
    uint32_t ctrl = (uint32_t)(res & 0x03);
    res >>= 2;

    uint32_t extra_bits = consumption_lut[ctrl];
    uint32_t payload = (uint32_t)(res & ((1 << extra_bits) - 1));

    s->bit_reservoir = res >> extra_bits;
    s->bits_available -= (2 + extra_bits);

    uint8_t hot_val = model->hot_values[payload + offset_lut[ctrl]];
    uint8_t raw_val = (uint8_t)payload;

    return (ctrl == 3) ? raw_val : hot_val;
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline void le_model_save(le_stream *s, const le_model *model)
{
    le_write_nibble(s, model->cold_num_bits);
    le_write_byte(s, model->cold_min);

    for(uint32_t i=0; i<LE_MODEL_MAX_HOT; ++i)
        le_write_byte(s, model->hot_values[i]);
}

// ----------------------------------------------------------------------------------------------------------------------------
void le_model_load(le_stream *s, le_model *model)
{
    model->cold_num_bits = le_read_nibble(s);
    model->cold_min = le_read_byte(s);
    

    for(uint32_t i=0; i<LE_MODEL_MAX_HOT; ++i)
        model->hot_values[i] = le_read_byte(s);
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline uint8_t zigzag8_encode(int8_t v)
{
    return (uint8_t)((v << 1) ^ (v >> 7));
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline int8_t zigzag8_decode(uint8_t v)
{
    return (int8_t)((v >> 1) ^ -(int8_t)(v & 1));
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline void le_encode_delta(le_stream *s, int8_t delta)
{
    uint8_t zz = zigzag8_encode(delta);

    for(uint32_t i = 0; i < 4; i++)
    {
        if(zz < 3)
        {
            le_write_dibit(s, zz);
            return;
        }

        le_write_dibit(s, 3); // escape
        zz -= 3;
    }

    le_write_byte(s, zigzag8_encode(delta));
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline int8_t le_decode_delta(le_stream* s)
{
    if (s->bits_available < 16) 
        le_refill(s);

    uint64_t res = s->bit_reservoir;
    uint64_t flipped = ~res;
    uint32_t first_zero_bit = __builtin_ctzll(flipped); 
    uint32_t escapes = first_zero_bit >> 1;

    if (escapes > 4) 
        escapes = 4;

    uint32_t bits_to_consume;
    uint8_t zz_value;

    if (escapes < 4)
    {
        bits_to_consume = (escapes + 1) << 1;
        uint32_t final_dibit = (uint32_t)((res >> (escapes << 1)) & 0x03);
        zz_value = (uint8_t)(final_dibit + (escapes * 3));
    }
    else
    {
        bits_to_consume = 16;
        zz_value = (uint8_t)((res >> 8) & 0xFF);
    }

    s->bit_reservoir = res >> bits_to_consume;
    s->bits_available -= bits_to_consume;

    return zigzag8_decode(zz_value);
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

// //----------------------------------------------------------------------------------------------------------------------------
// static inline void* ptr_shift(void *base, size_t shift_bytes)
// {
//     return (void *)((uint8_t *)base + shift_bytes);
// }

// //----------------------------------------------------------------------------------------------------------------------------
// static inline const void* ptr_shift_const(const void *base, size_t shift_bytes)
// {
//     return (const void *)((const uint8_t *)base + shift_bytes);
// }

//----------------------------------------------------------------------------------------------------------------------------
static inline const void* get_block(const void *base, size_t elem_size, uint32_t width_blocks, uint32_t x, uint32_t y)
{
    size_t row_index = (size_t)y * width_blocks;
    size_t idx = row_index + x;
    return (const uint8_t *)base + idx * elem_size;
}

//----------------------------------------------------------------------------------------------------------------------------
// static inline uint8_t delta_encode_wrap(uint8_t prev, uint8_t curr)
// {
//     return curr - prev; // automatically wraps modulo 256
// }

// //----------------------------------------------------------------------------------------------------------------------------
// static inline uint8_t delta_decode_wrap(uint8_t prev, uint8_t delta_encoded)
// {
//     return prev + delta_encoded; // automatically wraps modulo 256
// }

//----------------------------------------------------------------------------------------------------------------------
// static inline int int_abs(int a) {return (a>=0) ? a : -a;}

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

    // sort table for compression
    qsort(output, *num_entries, sizeof(uint32_t), compare_entries);
}

//----------------------------------------------------------------------------------------------------------------------------
// uint32_t nearest48(const uint64_t* table, uint32_t table_size, uint64_t bitfield)
// {
//     uint32_t scores[TABLE_SIZE];
//     uint32_t i = 0;

//     const uint64_t MASK48 = 0x0000FFFFFFFFFFFFull;

// #ifdef BC_CRUNCH_NEON
//     uint64x2_t bf_vec   = vdupq_n_u64(bitfield & MASK48);
//     uint64x2_t mask_vec = vdupq_n_u64(MASK48);

//     for (; i + 1 < table_size; i += 2)
//     {
//         uint64x2_t dict  = vandq_u64(vld1q_u64(&table[i]), mask_vec);
//         uint64x2_t delta = veorq_u64(dict, bf_vec);

//         uint8x16_t cnt = vcntq_u8(vreinterpretq_u8_u64(delta));

//         // split lanes correctly
//         uint8x8_t cnt_lo = vget_low_u8(cnt);
//         uint8x8_t cnt_hi = vget_high_u8(cnt);

//         // horizontal add per 64-bit lane
//         scores[i + 0] = vaddlv_u8(cnt_lo);
//         scores[i + 1] = vaddlv_u8(cnt_hi);
//     }
// #endif

//     // tail
//     for (; i < table_size; ++i)
//     {
//         uint64_t delta = (table[i] ^ bitfield) & MASK48;
//         scores[i] = (uint32_t)popcount64(delta);
//     }

//     // find best
//     uint32_t best_index = 0;
//     uint32_t best_score = UINT32_MAX;

//     for (uint32_t j = 0; j < table_size; ++j)
//     {
//         uint32_t score = scores[j];
//         if (score < best_score ||
//            (score == best_score && table[j] > table[best_index]))
//         {
//             best_score = score;
//             best_index = j;
//         }
//     }

//     return ((best_score & 0xffff) << 16) | (best_index & 0xffff);
// }

//----------------------------------------------------------------------------------------------------------------------------
// static inline uint8_t bc4_get_index(const bc4_block* b, uint32_t pixel_index)
// {
//     uint32_t bit_offset = pixel_index * 3;
//     uint64_t bits = ((uint64_t)b->indices[0]) | ((uint64_t)b->indices[1] << 16) | ((uint64_t)b->indices[2] << 32);
//     uint8_t index = (uint8_t)((bits >> bit_offset) & 0x7);
//     return index;
// }

//----------------------------------------------------------------------------------------------------------------------------
// static inline void bc4_set_index(bc4_block* b, uint32_t pixel_index, uint8_t data)
// {
//     uint32_t bit_offset = pixel_index * 3;
//     uint32_t word_index = bit_offset >> 4;
//     uint32_t bit_in_word = bit_offset & 0xF;

//     uint16_t mask = 0x7 << bit_in_word;
//     b->indices[word_index] = (b->indices[word_index] & ~mask) | ((data & 0x7) << bit_in_word);

//     // if the 3 bits spill into the next word
//     if (bit_in_word > 13)  // last 2 or 1 bits spill
//     {
//         uint16_t spill_bits = (data & 0x7) >> (16 - bit_in_word);
//         b->indices[word_index + 1] = (b->indices[word_index + 1] & ~(0x7 >> (16 - bit_in_word))) | spill_bits;
//     }
// }

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

// static const uint32_t block_zigzag[16] = {0,  1,  2,  3, 7,  6,  5,  4, 8,  9, 10, 11, 15, 14, 13, 12};

//----------------------------------------------------------------------------------------------------------------------------
// as we use a static huffman entropy encoder to be the fastest at decompression, we need two passes :
//   - first pass : for each model, build a histogram and compute probabilities.
//   - second pass : save the model in the stream,  then compress the texture
//
// static models are needed for decompression obivously
void bc1_crunch(le_stream* restrict codec, void* restrict cruncher_memory, const void* restrict input, size_t stride, uint32_t width, uint32_t height)
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

    // ----------
    // FIRST PASS
    // ----------
    le_histogram htable_index, hmask, htable_difference;
    histogram_init(&htable_index, 256);
    histogram_init(&hmask, 16);
    histogram_init(&htable_difference, 256);

    le_histogram htable_entry;
    histogram_init(&htable_entry, 256);

    for(uint32_t i=1; i<top_table_size; ++i)
    {
        // table is sorted from small to big, so diff is always positive
        uint32_t diff = top_table[i] - top_table[i-1];

        for(uint32_t j=0; j<4; ++j)
            htable_entry.count[(diff >> (j*8)) & 0xff]++;
    }

    bc1_block previous = {0};
    for(uint32_t y = 0; y < height_blocks; ++y)
    {
        for(uint32_t x = 0; x < width_blocks; ++x)
        {
            uint32_t zigzag_x = (y&1) ? x : width_blocks - x - 1;
            const bc1_block* current = get_block(input, stride, width_blocks, zigzag_x, y);

            uint32_t reference = nearest32(top_table, top_table_size, current->indices) & 0xffff;

            assert(reference < 256);
            htable_index.count[reference]++;

            uint32_t difference = current->indices ^ top_table[reference];

            uint32_t mask = 0;
            if ((difference & 0x000000FF) != 0) mask |= 1;
            if ((difference & 0x0000FF00) != 0) mask |= 2;
            if ((difference & 0x00FF0000) != 0) mask |= 4;
            if ((difference & 0xFF000000) != 0) mask |= 8;

            hmask.count[mask]++;

            for(uint32_t j=0; j<4; ++j)
                if (mask & (1u << j))
                    htable_difference.count[(difference >> (j*8)) & 0xff]++;

            previous = *current;
        }
    }

    le_model table_index, diff_mask, table_entry, table_difference;
    le_model_init(&table_index, htable_index.count, htable_index.num_symbols);
    le_model_init(&diff_mask, hmask.count, hmask.num_symbols);
    le_model_init(&table_entry, htable_entry.count, htable_entry.num_symbols);
    le_model_init(&table_difference, htable_difference.count, htable_difference.num_symbols);

    le_model_save(codec, &table_index);
    le_model_save(codec, &diff_mask);
    le_model_save(codec, &table_entry);
    le_model_save(codec, &table_difference);

    // write the top-table
    le_write_byte(codec, top_table_size-1);
    for(uint32_t j=0; j<4; ++j)
        le_write_byte(codec, (top_table[0] >> (j*8)) & 0xff);    // first entry not compressed

    for(uint32_t i=1; i<top_table_size; ++i)
    {
        // table is sorted from small to big, so diff is always positive
        uint32_t diff = top_table[i] - top_table[i-1];

        for(uint32_t j=0; j<4; ++j)
            le_encode(codec, &table_entry, (diff >> (j*8)) & 0xff);
    }

    // ----------
    // SECOND PASS
    // ----------

    previous = (bc1_block) {0};
    for(uint32_t y = 0; y < height_blocks; ++y)
    {
        for(uint32_t x = 0; x < width_blocks; ++x)
        {
            // zig-zag pattern delta compression for colors
            uint32_t zigzag_x = (y&1) ? x : width_blocks - x - 1;
            const bc1_block* current = get_block(input, stride, width_blocks, zigzag_x, y);
            for(uint32_t j=0; j<2; ++j)
            {
                uint8_t current_red, current_green, current_blue;
                uint8_t previous_red, previous_green, previous_blue;

                bc1_extract_565(current->color[j], &current_red, &current_green, &current_blue);
                bc1_extract_565(previous.color[j], &previous_red, &previous_green, &previous_blue);

                if (y>0 && x!=0)
                {
                    const bc1_block* up = get_block(input, stride, width_blocks, zigzag_x, y-1);
                    uint8_t up_red, up_green, up_blue;
                    bc1_extract_565(up->color[j], &up_red, &up_green, &up_blue);

                    if (y>0 && x!=0)
                    {
                        const bc1_block* up = get_block(input, stride, width_blocks, zigzag_x, y-1);
                        uint8_t up_red, up_green, up_blue;
                        bc1_extract_565(up->color[j], &up_red, &up_green, &up_blue);

                        previous_red = (previous_red + up_red) / 2;
                        previous_green = (previous_green + up_green) / 2;
                        previous_blue = (previous_blue + up_blue) / 2;
                    }
                }

                int dred = current_red - previous_red;
                int dgreen = current_green - previous_green;
                int dblue = current_blue - previous_blue;

                // first encode green delta
                le_encode_delta(codec, (int8_t) dgreen);

                // then encode red and blue delta based on green delta
                // assuming some relation between green and other components
                dgreen /= 2;
                dred -= dgreen;
                dblue -= dgreen;

                le_encode_delta(codec, (int8_t) (dred));
                le_encode_delta(codec, (int8_t) (dblue));
            }

            // for indices, we store the reference to "nearest" indices (can be exactly the same)
            // and the delta with this reference
            uint32_t reference = nearest32(top_table, top_table_size, current->indices) & 0xffff;
            le_encode(codec, &table_index, reference);

            // xor the difference and encode (could be 0 if equal to reference)
            uint32_t difference = current->indices ^ top_table[reference];

            uint32_t mask = 0;
            if ((difference & 0x000000FF) != 0) mask |= 1;
            if ((difference & 0x0000FF00) != 0) mask |= 2;
            if ((difference & 0x00FF0000) != 0) mask |= 4;
            if ((difference & 0xFF000000) != 0) mask |= 8;

            le_write_nibble(codec, mask);

            // only encode the bytes that are actually non-zero
            for(uint32_t j=0; j<4; ++j)
                if (mask & (1u << j))
                    le_encode(codec, &table_difference, (difference >> (j*8)) & 0xff);

            previous = *current;
        }
    }

    // printf("hmask histogram\n");

    // for(uint32_t i=0; i<16; ++i)
    //     printf("h[%u] = %u\t", i, hmask.count[i]);
    
    // printf("\n");
}

//----------------------------------------------------------------------------------------------------------------------------
void bc1_decrunch(le_stream* restrict codec, uint32_t width, uint32_t height, void* restrict output, size_t stride)
{
    assert((width % 4 == 0) && (height % 4 == 0));

    uint32_t height_blocks = height/4;
    uint32_t width_blocks = width/4;
    
    le_model table_index, diff_mask, table_entry, table_difference;
    le_model_load(codec, &table_index);
    le_model_load(codec, &diff_mask);
    le_model_load(codec, &table_entry);
    le_model_load(codec, &table_difference);

    uint32_t top_table[TABLE_SIZE];
    uint32_t top_table_size = le_read_byte(codec)+1;

    top_table[0] = 0;
    for(uint32_t j=0; j<4; ++j)
        top_table[0] |= le_read_byte(codec) << (j*8);

    for(uint32_t i=1; i<top_table_size; ++i)
    {
        uint32_t diff = 0;
        for(uint32_t j=0; j<4; ++j)
            diff |= le_decode(codec, &table_entry) << (j*8);

        top_table[i] = top_table[i-1] + diff;
    }

    bc1_block previous = {0};
    for(uint32_t y = 0; y < height_blocks; ++y)
    {
        for(uint32_t x = 0; x < width_blocks; ++x)
        {
            // zig-zag pattern color
            uint32_t zigzag_x = (y&1) ? x : width_blocks - x - 1;
            bc1_block* current = (bc1_block*) get_block(output, stride, width_blocks, zigzag_x, y);
            for (uint32_t j = 0; j < 2; ++j)
            {
                uint8_t reference_red, reference_green, reference_blue;
                bc1_extract_565(previous.color[j], &reference_red, &reference_green, &reference_blue);

                if (y>0 && x!=0)
                {
                    const bc1_block* up = get_block(output, stride, width_blocks, zigzag_x, y-1);
                    uint8_t up_red, up_green, up_blue;
                    bc1_extract_565(up->color[j], &up_red, &up_green, &up_blue);

                    reference_red = (reference_red + up_red) / 2;
                    reference_green = (reference_green + up_green) / 2;
                    reference_blue = (reference_blue + up_blue) / 2;
                }

                int8_t delta_green = le_decode_delta(codec);
                int8_t delta_red = le_decode_delta(codec);
                int8_t delta_blue = le_decode_delta(codec);

                // red and blue delta are based on green delta
                int dgreen_orig = delta_green;
                int current_green_value = reference_green + dgreen_orig;
                int dgreen_halved = dgreen_orig / 2;

                int dred_orig = delta_red + dgreen_halved;
                int dblue_orig = delta_blue + dgreen_halved;

                int current_red_value = reference_red + dred_orig;
                int current_blue_value = reference_blue + dblue_orig;

                current->color[j] = bc1_pack_565((uint8_t)current_red_value, (uint8_t)current_green_value, (uint8_t)current_blue_value);
            }

            // indices difference with top table
            uint32_t reference = le_decode(codec, &table_index);
            uint32_t mask = le_read_nibble(codec);

            uint32_t difference=0;
            for(uint32_t j=0; j<4; ++j)
                if (mask & (1 << j))
                    difference = difference | (le_decode(codec, &table_difference) << (j*8));

            current->indices =  difference ^ top_table[reference];

            previous = *current;
        }
    }
}

//----------------------------------------------------------------------------------------------------------------------------
// void bc4_crunch(range_codec* codec, void* cruncher_memory, const void* input, size_t stride, uint32_t width, uint32_t height)
// {
//     assert((width%4 == 0) && (height%4 == 0));
//     assert(((uintptr_t)cruncher_memory)%sizeof(uintptr_t) == 0);

//     uint32_t height_blocks = height/4;
//     uint32_t width_blocks = width/4;

//     range_model color_delta[2];
//     model_init(&color_delta[0], 1<<BC4_COLOR_NUM_BITS);
//     model_init(&color_delta[1], 1<<BC4_COLOR_NUM_BITS);

//     range_model color_reference, first_index, use_dict, dict_reference; 
//     model_init(&color_reference, 2);
//     model_init(&first_index, 1<<3);
//     model_init(&use_dict, 2);
//     model_init(&dict_reference, DICTIONARY_SIZE);

//     range_model indices[24];
//     for(uint32_t i=0; i<24; ++i)
//         model_init(&indices[i], 1<<BC4_INDEX_NUM_BITS);

//     range_model dict_delta[16];
//     for(uint32_t i=0; i<16; ++i)
//         model_init(&dict_delta[i], 1<<3);

//     bc4_block previous = {.color = {0, 128}};

//     // dictionary initialization
//     uint64_t dictionary[DICTIONARY_SIZE];
//     for(uint32_t i=0; i<DICTIONARY_SIZE; ++i)
//         dictionary[i] = UINT64_MAX;

//     for(uint32_t y = 0; y < height_blocks; ++y)
//     {
//         for(uint32_t x = 0; x < width_blocks; ++x)
//         {
//             const bc4_block* current = get_block(input, stride, width_blocks, x, y);

//             int reference = previous.color[0];
//             if (y>0)
//             {
//                 const bc4_block* up = get_block(input, stride, width_blocks, x, y-1);
//                 if (x>0)
//                 {
//                     const bc4_block* up_left = get_block(input, stride, width_blocks, x-1, y-1);
//                     reference += up->color[0] - up_left->color[0];
//                 }
//                 else
//                     reference = up->color[0];
//             }

//             if (reference < 0) reference = 0;
//             if (reference > 255) reference = 255;

//             enc_put(codec, &color_delta[0], delta_encode_wrap((uint8_t)reference, current->color[0]));
//             enc_put(codec, &color_delta[1], delta_encode_wrap(current->color[0], current->color[1]));

//             // search in the dictionary for the current bitfield
//             uint64_t bitfield = MAKE48(current->indices[0], current->indices[1], current->indices[2]);
//             uint32_t dict_lookup = nearest48(dictionary, DICTIONARY_SIZE, bitfield);
//             uint16_t score = dict_lookup>>16;
//             uint16_t found_index = dict_lookup&0xffff;
            
//             // found or similar? just write the dictionary index
//             if (score < 5 && ((y*width_blocks) + x > 32))
//             {
//                 enc_put(codec, &use_dict, 1);
//                 enc_put(codec, &dict_reference, found_index);
                
//                 uint64_t reference = dictionary[found_index];
//                 uint64_t bitfield_delta = reference ^ bitfield;
//                 for(uint32_t j=0; j<16; ++j)
//                     enc_put(codec, &dict_delta[j], (bitfield_delta>>(j*3))&0x7);

//                 if(found_index > 0)
//                 {
//                     uint64_t temp = dictionary[found_index];
//                     uint32_t target = found_index / 2;
//                     memmove(&dictionary[target+1], &dictionary[target], (found_index - target) * sizeof(uint64_t));
//                     dictionary[target] = temp;
//                 }
//             }
//             else
//             {
//                 // store the entry in the middle of dictionary
//                 uint32_t middle = DICTIONARY_SIZE/2;
//                 memmove(&dictionary[middle+1], &dictionary[middle], (DICTIONARY_SIZE - middle - 1) * sizeof(uint64_t));
//                 dictionary[middle] = bitfield;

//                 // write the indices with local difference delta encoded
//                 enc_put(codec, &use_dict, 0);

//                 uint8_t block_previous = bc4_get_index(current, 0);
//                 enc_put(codec, &first_index, block_previous);

//                 range_model* model = bc4_select_model(current, indices);
//                 for(uint32_t j=1; j<16; ++j)
//                 {
//                     uint8_t data = bc4_get_index(current, block_zigzag[j]);
//                     enc_put(codec, &model[block_previous], block_previous ^ data);
//                     block_previous = data;
//                 }
//             }
//             previous = *current;
//         }
//     }
// }

//----------------------------------------------------------------------------------------------------------------------------
// void bc4_decrunch(range_codec* codec, uint32_t width, uint32_t height, void* output, size_t stride)
// {
//     assert((width % 4 == 0) && (height % 4 == 0));

//     uint32_t height_blocks = height/4;
//     uint32_t width_blocks = width/4;

//     range_model color_delta[2];
//     model_init(&color_delta[0], 1<<BC4_COLOR_NUM_BITS);
//     model_init(&color_delta[1], 1<<BC4_COLOR_NUM_BITS);

//     range_model color_reference, first_index, use_dict, dict_reference;
//     model_init(&color_reference, 2);
//     model_init(&first_index, 1<<3);
//     model_init(&use_dict, 2);
//     model_init(&dict_reference, DICTIONARY_SIZE);

//     range_model indices[24];
//     for(uint32_t i=0; i<24; ++i)
//         model_init(&indices[i], 1<<BC4_INDEX_NUM_BITS);

//     range_model dict_delta[16];
//     for(uint32_t i=0; i<16; ++i)
//         model_init(&dict_delta[i], 1<<3);

//     bc4_block previous = {.color = {0, 128}};

//     // dictionary initialization
//     uint64_t dictionary[DICTIONARY_SIZE];
//     for(uint32_t i=0; i<DICTIONARY_SIZE; ++i)
//         dictionary[i] = UINT64_MAX;

//     for(uint32_t y = 0; y < height_blocks; ++y)
//     {
//         for(uint32_t x = 0; x < width_blocks; ++x)
//         {
//             bc4_block* current = (bc4_block*) get_block(output, stride, width_blocks, x, y);
//             int reference = previous.color[0];
//             if (y>0)
//             {
//                 const bc4_block* up = get_block(output, stride, width_blocks, x, y-1);
//                 if (x>0)
//                 {
//                     const bc4_block* up_left = get_block(output, stride, width_blocks, x-1, y-1);
//                     reference += up->color[0] - up_left->color[0];
//                 }
//                 else
//                     reference = up->color[0];
//             }

//             if (reference < 0) reference = 0;
//             if (reference > 255) reference = 255;

//             current->color[0] = delta_decode_wrap((uint8_t)reference, dec_get(codec, &color_delta[0]));
//             current->color[1] = delta_decode_wrap(current->color[0], dec_get(codec, &color_delta[1]));

//             if (dec_get(codec, &use_dict))
//             {
//                 // data should be in the dictionary
//                 uint32_t found_index = dec_get(codec, &dict_reference);
//                 uint64_t reference = dictionary[found_index];
//                 uint64_t bitfield = 0;
    
//                 for(uint32_t j=0; j<16; ++j)
//                 {
//                     uint64_t byte = (uint64_t)dec_get(codec, &dict_delta[j]);
//                     bitfield |= (byte << (j*3));
//                 }

//                 bitfield ^= reference;
//                 current->indices[0] = ((bitfield>>32) & 0xffff);
//                 current->indices[1] = ((bitfield>>16) & 0xffff);
//                 current->indices[2] = bitfield & 0xffff;

//                 if(found_index > 0)
//                 {
//                     uint64_t temp = dictionary[found_index];
//                     uint32_t target = found_index / 2;  // bring the hit up but no in front (multiple hit will do that)
//                     memmove(&dictionary[target+1], &dictionary[target], (found_index - target) * sizeof(uint64_t));
//                     dictionary[target] = temp;
//                 }
//             }
//             else
//             {
//                 uint8_t block_previous = dec_get(codec, &first_index);
//                 bc4_set_index(current, 0, block_previous);

//                 range_model* model = bc4_select_model(current, indices);
//                 for(uint32_t j=1; j<16; ++j)
//                 {
//                     uint8_t delta = dec_get(codec, &model[block_previous]);
//                     uint8_t data = block_previous ^ delta;
//                     bc4_set_index(current, block_zigzag[j], data);
//                     block_previous = data;
//                 }

//                 // store the entry in the middle of dictionary
//                 uint32_t middle = DICTIONARY_SIZE/2;
//                 memmove(&dictionary[middle+1], &dictionary[middle], (DICTIONARY_SIZE - middle - 1) * sizeof(uint64_t));
//                 dictionary[middle] = MAKE48(current->indices[0], current->indices[1], current->indices[2]);
//             }
//             previous = *current;
//         }
//     }
// }


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

    le_stream codec;
    le_init(&codec, output, length);
    le_begin_encode(&codec);

    switch(format)
    {
    case bc1 : 
        {
            bc1_crunch(&codec, cruncher_memory, input, sizeof(bc1_block), width, height);
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
    return le_end_encode(&codec);
}

//----------------------------------------------------------------------------------------------------------------------------
void bc_decrunch(const void* input, size_t length, uint32_t width, uint32_t height, enum bc_format format, void* output)
{
    le_stream codec;
    le_init(&codec, (void*)input, length);
    le_begin_decode(&codec);

    switch(format)
    {
    case bc1 : bc1_decrunch(&codec, width, height, output, sizeof(bc1_block)); break;
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

    le_end_decode(&codec);
}
