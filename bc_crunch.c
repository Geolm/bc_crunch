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


#define LE_ALPHABET_SIZE (256)
#define LE_K_TREND_THRESHOLD (12)

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
    uint8_t alphabet[LE_ALPHABET_SIZE];
    uint8_t k;  // rice k-value
    int8_t k_trend;
} le_model;

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

// ----------------------------------------------------------------------------------------------------------------------------
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
static inline void le_write_bits(le_stream* s, uint8_t data, uint8_t num_bits)
{
    s->bit_reservoir |= (uint64_t)(data & ((1 << num_bits) - 1)) << s->bits_available;
    s->bits_available += num_bits;
    if (s->bits_available >= 32)
        le_flush(s);
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline uint8_t le_read_bits(le_stream* s, uint8_t num_bits)
{
    if (s->bits_available < num_bits)
        le_refill(s);

    uint8_t value = (uint8_t)(s->bit_reservoir & ((1U<<num_bits)-1U));
    s->bit_reservoir >>= num_bits;
    s->bits_available -= num_bits;
    return value;
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline void le_write_byte(le_stream* s, uint8_t value)
{
    s->bit_reservoir |= ((uint64_t)value << s->bits_available);
    s->bits_available += 8;
    if (s->bits_available >= 32)
        le_flush(s);
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline uint8_t le_read_byte(le_stream* s)
{
    if (s->bits_available < 8)
        le_refill(s);

    uint8_t value = (uint8_t)(s->bit_reservoir & 0xFF);
    s->bit_reservoir >>= 8;
    s->bits_available -= 8;
    return value;
}


// ----------------------------------------------------------------------------------------------------------------------------
void le_model_init(le_model *model)
{
    for(uint32_t i=0; i<LE_ALPHABET_SIZE; ++i)
        model->alphabet[i] = i;
    model->k = 2;
    model->k_trend = 0;
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline void rice_encode(le_stream *s, uint32_t value, uint8_t k) 
{
    uint32_t q = value >> k;
    uint32_t r = value & ((1U << k) - 1U);

    // write q
    for (uint32_t i = 0; i < q; ++i)
        le_write_bits(s, 1, 1);
    
    // terminator '0'
    le_write_bits(s, 0, 1);

    // remainder
    if (k > 0) 
        le_write_bits(s, (uint8_t)r, k);
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline void le_model_update(le_model* model, uint8_t value)
{
    if (value < (1U << model->k) && model->k > 0) 
        model->k_trend--;
    else if (value > (3U << model->k) && model->k < 6) 
        model->k_trend++;

    // soft adaptation
    if (model->k_trend > LE_K_TREND_THRESHOLD)
    {
        model->k++;
        model->k_trend = 0;
    }
    else if (model->k_trend < -LE_K_TREND_THRESHOLD)
    {
        model->k--;
        model->k_trend = 0;
    }
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline void le_encode(le_stream *s, le_model *model, uint8_t value)
{
#ifdef LE_CHECKS
    assert(s->mode == le_mode_encode);
#endif

    uint32_t index = 0;
    for (; index < 256; index++)
        if (model->alphabet[index] == value)
            break;

    rice_encode(s, index, model->k);

    // move up this value in the alphabet
    if (index > 0) 
    {
        uint8_t temp = model->alphabet[index];
        uint32_t target_index = index / 2;  // lowpass filter, prevent jittering

        for (uint32_t i = index; i > target_index; i--)
            model->alphabet[i] = model->alphabet[i - 1];
        
        model->alphabet[target_index] = temp;
    }

    le_model_update(model, index);
}


// ----------------------------------------------------------------------------------------------------------------------------
static inline uint8_t rice_decode(le_stream *s, uint8_t k) 
{
        uint32_t q = 0;
    while (true) 
    {
        if (s->bits_available == 0) 
            le_refill(s);
        
        if ((s->bit_reservoir & 1ULL) != 0) 
        {
            q++;
            s->bit_reservoir >>= 1;
            s->bits_available--;
        } else 
        {
            s->bit_reservoir >>= 1; 
            s->bits_available--;
            break;
        }
    }

    if (s->bits_available < k)
        le_refill(s);
    
    uint32_t r = 0;
    if (k > 0) 
    {
        r = (uint32_t)(s->bit_reservoir & ((1ULL << k) - 1U));
        s->bit_reservoir >>= k;
        s->bits_available -= k;
    }

    return (uint8_t) ((q << k) | r);
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline uint8_t le_decode(le_stream *restrict s, le_model *restrict model) 
{
#ifdef LE_CHECKS
    assert(index < LE_ALPHABET_SIZE);
#endif

    uint8_t index = rice_decode(s, model->k);
    uint8_t value = model->alphabet[index];

     // move up this value in the alphabet
    if (index > 0) 
    {
        uint8_t temp = model->alphabet[index];
        uint32_t target_index = index / 2;  // lowpass filter, prevent jittering

        for (uint32_t i = index; i > target_index; i--)
            model->alphabet[i] = model->alphabet[i - 1];
        
        model->alphabet[target_index] = temp;
    }

    le_model_update(model, index);

    return value;
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
static inline void le_encode_delta(le_stream *s, le_model* model, int8_t delta)
{
    uint8_t zz = zigzag8_encode(delta);
    rice_encode(s, zz, model->k);
    le_model_update(model, zz);
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline int8_t le_decode_delta(le_stream* s, le_model* model)
{
    uint8_t zz = rice_decode(s, model->k);
    le_model_update(model, zz);
    return zigzag8_decode(zz);
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
    size_t index = row_index + x;
    return (const uint8_t *)base + index * elem_size;
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
// static inline le_model* bc4_select_model(const bc4_block* b, le_model* indices)
// {
//     int endpoints_delta = int_abs(b->color[0] - b->color[1]);

//     if (endpoints_delta < 8)
//         return &indices[0];
//     else if (endpoints_delta < 32)
//         return &indices[8];

//     return &indices[16];
// }

// static const uint32_t block_zigzag[16] = {0,  1,  2,  3, 7,  6,  5,  4, 8,  9, 10, 11, 15, 14, 13, 12};

//----------------------------------------------------------------------------------------------------------------------
static inline int int_abs(int a) {return (a>=0) ? a : -a;}

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

    // write the table
    le_model table_entry;
    le_model_init(&table_entry);

    le_write_byte(codec, top_table_size-1);   // entries count
    for(uint32_t j=0; j<4; ++j)
        le_write_bits(codec, (top_table[0] >> (j*8)) & 0xff, 8);    // first entry not compressed

    for(uint32_t i=1; i<top_table_size; ++i)
    {
        // table is sorted from small to big, so diff is always positive
        uint32_t diff = top_table[i] - top_table[i-1];

        for(uint32_t j=0; j<4; ++j)
            le_encode(codec, &table_entry, (diff >> (j*8)) & 0xff);
    }

    le_model red, green, blue;
    le_model_init(&red);
    le_model_init(&green);
    le_model_init(&blue);

    le_model table_index, table_difference, diff_mask, color_reference;
    le_model_init(&table_index);
    le_model_init(&table_difference);
    le_model_init(&diff_mask);
    le_model_init(&color_reference);

    bc1_block previous = {0};

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

                    int previous_delta = int_abs(current_red - previous_red) + int_abs(current_green-previous_green) + int_abs(current_blue-previous_blue);
                    int up_delta = int_abs(current_red-up_red) + int_abs(current_green-up_green) + int_abs(current_blue-up_blue);

                    le_write_bits(codec, (up_delta < previous_delta) ? 1 : 0, 1);

                    // overwrite previous value to avoid using a new set of variables
                    if (up_delta < previous_delta)
                    {
                        previous_red = up_red;
                        previous_green = up_green;
                        previous_blue = up_blue;
                    }
                }

                int dred = current_red - previous_red;
                int dgreen = current_green - previous_green;
                int dblue = current_blue - previous_blue;

                // first encode green delta
                le_encode_delta(codec, &green, dgreen);

                // then encode red and blue delta based on green delta
                // assuming some relation between green and other components
                dgreen /= 2;
                dred -= dgreen;
                dblue -= dgreen;

                le_encode_delta(codec, &red, dred);
                le_encode_delta(codec, &blue, dblue);
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

            le_encode(codec, &diff_mask, mask);

            // only encode the bytes that are actually non-zero
            for(uint32_t j=0; j<4; ++j)
                if (mask & (1u << j))
                    le_encode(codec, &table_difference, (difference >> (j*8)) & 0xff);

            previous = *current;
        }
    }
}

//----------------------------------------------------------------------------------------------------------------------------
void bc1_decrunch(le_stream* codec, uint32_t width, uint32_t height, void* output, size_t stride)
{
    assert((width % 4 == 0) && (height % 4 == 0));

    uint32_t height_blocks = height/4;
    uint32_t width_blocks = width/4;

    le_model red, green, blue;
    le_model_init(&red);
    le_model_init(&green);
    le_model_init(&blue);

    le_model table_entry;
    le_model_init(&table_entry);

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

    le_model table_index, table_difference, diff_mask;
    le_model_init(&table_index);
    le_model_init(&table_difference);
    le_model_init(&diff_mask);

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
                if (y>0 && x!=0 && le_read_bits(codec, 1) == 1)
                {
                    bc1_block* up = (bc1_block*) get_block(output, stride, width_blocks, zigzag_x, y-1);
                    bc1_extract_565(up->color[j], &reference_red, &reference_green, &reference_blue);
                }

                int delta_green = le_decode_delta(codec, &green);
                int delta_red = le_decode_delta(codec, &red);
                int delta_blue = le_decode_delta(codec, &blue);

                // red and blue delta are based on green delta
                int current_green_value = reference_green + delta_green;
                int dgreen_halved = delta_green / 2;

                int dred_orig = delta_red + dgreen_halved;
                int dblue_orig = delta_blue + dgreen_halved;

                int current_red_value = reference_red + dred_orig;
                int current_blue_value = reference_blue + dblue_orig;

                current->color[j] = bc1_pack_565((uint8_t)current_red_value, (uint8_t)current_green_value, (uint8_t)current_blue_value);
            }

            // indices difference with top table
            uint32_t reference = le_decode(codec, &table_index);
            uint32_t mask = le_decode(codec, &diff_mask);

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

//     le_model color_delta[2];
//     le_model_init(&color_delta[0], 1<<BC4_COLOR_NUM_BITS);
//     le_model_init(&color_delta[1], 1<<BC4_COLOR_NUM_BITS);

//     le_model color_reference, first_index, use_dict, dict_reference; 
//     le_model_init(&color_reference, 2);
//     le_model_init(&first_index, 1<<3);
//     le_model_init(&use_dict, 2);
//     le_model_init(&dict_reference, DICTIONARY_SIZE);

//     le_model indices[24];
//     for(uint32_t i=0; i<24; ++i)
//         le_model_init(&indices[i], 1<<BC4_INDEX_NUM_BITS);

//     le_model dict_delta[16];
//     for(uint32_t i=0; i<16; ++i)
//         le_model_init(&dict_delta[i], 1<<3);

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

//             le_encode(codec, &color_delta[0], delta_encode_wrap((uint8_t)reference, current->color[0]));
//             le_encode(codec, &color_delta[1], delta_encode_wrap(current->color[0], current->color[1]));

//             // search in the dictionary for the current bitfield
//             uint64_t bitfield = MAKE48(current->indices[0], current->indices[1], current->indices[2]);
//             uint32_t dict_lookup = nearest48(dictionary, DICTIONARY_SIZE, bitfield);
//             uint16_t score = dict_lookup>>16;
//             uint16_t found_index = dict_lookup&0xffff;
            
//             // found or similar? just write the dictionary index
//             if (score < 5 && ((y*width_blocks) + x > 32))
//             {
//                 le_encode(codec, &use_dict, 1);
//                 le_encode(codec, &dict_reference, found_index);
                
//                 uint64_t reference = dictionary[found_index];
//                 uint64_t bitfield_delta = reference ^ bitfield;
//                 for(uint32_t j=0; j<16; ++j)
//                     le_encode(codec, &dict_delta[j], (bitfield_delta>>(j*3))&0x7);

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
//                 le_encode(codec, &use_dict, 0);

//                 uint8_t block_previous = bc4_get_index(current, 0);
//                 le_encode(codec, &first_index, block_previous);

//                 le_model* model = bc4_select_model(current, indices);
//                 for(uint32_t j=1; j<16; ++j)
//                 {
//                     uint8_t data = bc4_get_index(current, block_zigzag[j]);
//                     le_encode(codec, &model[block_previous], block_previous ^ data);
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

//     le_model color_delta[2];
//     le_model_init(&color_delta[0], 1<<BC4_COLOR_NUM_BITS);
//     le_model_init(&color_delta[1], 1<<BC4_COLOR_NUM_BITS);

//     le_model color_reference, first_index, use_dict, dict_reference;
//     le_model_init(&color_reference, 2);
//     le_model_init(&first_index, 1<<3);
//     le_model_init(&use_dict, 2);
//     le_model_init(&dict_reference, DICTIONARY_SIZE);

//     le_model indices[24];
//     for(uint32_t i=0; i<24; ++i)
//         le_model_init(&indices[i], 1<<BC4_INDEX_NUM_BITS);

//     le_model dict_delta[16];
//     for(uint32_t i=0; i<16; ++i)
//         le_model_init(&dict_delta[i], 1<<3);

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

//             current->color[0] = delta_decode_wrap((uint8_t)reference, le_decode(codec, &color_delta[0]));
//             current->color[1] = delta_decode_wrap(current->color[0], le_decode(codec, &color_delta[1]));

//             if (le_decode(codec, &use_dict))
//             {
//                 // data should be in the dictionary
//                 uint32_t found_index = le_decode(codec, &dict_reference);
//                 uint64_t reference = dictionary[found_index];
//                 uint64_t bitfield = 0;
    
//                 for(uint32_t j=0; j<16; ++j)
//                 {
//                     uint64_t byte = (uint64_t)le_decode(codec, &dict_delta[j]);
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
//                 uint8_t block_previous = le_decode(codec, &first_index);
//                 bc4_set_index(current, 0, block_previous);

//                 le_model* model = bc4_select_model(current, indices);
//                 for(uint32_t j=1; j<16; ++j)
//                 {
//                     uint8_t delta = le_decode(codec, &model[block_previous]);
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
