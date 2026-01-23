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


//----------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------
enum codec_mode
{
    mode_undefined = 0,
    mode_encoder = 1,
    mode_decoder
};

//----------------------------------------------------------------------------------------------------------------------------
#define HF_MAX_SYMBOLS 256
#define HF_TABLE_BITS  11
#define HF_TABLE_SIZE  (1 << HF_TABLE_BITS)

typedef struct hf_model 
{
    int32_t  decd_table[HF_TABLE_SIZE];
    int32_t  tree[2][HF_MAX_SYMBOLS];
    uint32_t length[HF_MAX_SYMBOLS];
    uint32_t codeword[HF_MAX_SYMBOLS];
    uint32_t symbol_count[HF_MAX_SYMBOLS];
    uint32_t data_symbols;
    uint32_t table_shift;
    uint32_t table_bits;
} hf_model;

//----------------------------------------------------------------------------------------------------------------------------
typedef struct hf_context
{
    uint64_t bit_buffer;
    int32_t bits_in_buffer;
    uint8_t *bc_pointer;
    uint8_t *code_buffer;
    size_t buffer_size;
    enum codec_mode mode;
} hf_context;

//----------------------------------------------------------------------------------------------------------------------------
static void form_huffman_tree(int32_t number_of_symbols, uint32_t original_count[], int32_t left_branch[], int32_t right_branch[])
{
    uint32_t h_count[HF_MAX_SYMBOLS + 1];
    int32_t h_sym[HF_MAX_SYMBOLS + 1];
    int32_t h_size = 0;

    for (int32_t i = 0; i < number_of_symbols; i++) 
    {
        if (original_count[i] > 0) 
        {
            h_size++;
            int32_t c = h_size;
            while (c > 1 && original_count[i] < h_count[c >> 1]) 
            {
                h_count[c] = h_count[c >> 1];
                h_sym[c] = h_sym[c >> 1];
                c >>= 1;
            }
            h_count[c] = original_count[i];
            h_sym[c] = i; 
        }
    }

    if (h_size <= 1) 
    {
        // Force at least a root so the rest of the code doesn't crash
        left_branch[1] = h_size ? h_sym[1] : 0;
        right_branch[1] = h_size ? h_sym[1] : 0;
        return;
    }

    int32_t next_internal = 2; 

    while (h_size > 1) 
    {
        int32_t s[2]; uint32_t c[2];
        for (int k = 0; k < 2; k++) 
        {
            s[k] = h_sym[1]; c[k] = h_count[1];
            uint32_t lc = h_count[h_size];
            int32_t ls = h_sym[h_size--];
            int32_t i = 1, j;
            while ((j = i << 1) <= h_size) 
            {
                if (j < h_size && h_count[j+1] < h_count[j]) j++;
                if (lc <= h_count[j]) break;
                h_count[i] = h_count[j]; h_sym[i] = h_sym[j];
                i = j;
            }
            h_count[i] = lc; h_sym[i] = ls;
        }

        // Determine where to save this internal node
        // The very last merge MUST be stored at index 1 (the root)
        int32_t node_idx = (h_size == 0) ? 1 : next_internal++;
        
        assert(node_idx < HF_MAX_SYMBOLS);

        left_branch[node_idx] = s[0];
        right_branch[node_idx] = s[1];

        // Push internal node back to heap (negative value)
        uint32_t nc = c[0] + c[1];
        int32_t cur = ++h_size;
        while (cur > 1 && nc < h_count[cur >> 1]) 
        {
            h_count[cur] = h_count[cur >> 1];
            h_sym[cur] = h_sym[cur >> 1];
            cur >>= 1;
        }
        h_count[cur] = nc;
        h_sym[cur] = -node_idx;
    }
}

//----------------------------------------------------------------------------------------------------------------------------
static void set_huffman_code(int32_t table_bits, int32_t tree[2][HF_MAX_SYMBOLS], int32_t *decoder_table, uint32_t *codeword, uint32_t *codeword_length)
{
    int32_t stack[32], node, depth = 0, current_code = 0;

    for (stack[0] = 1; depth >= 0;)
    {
        // transverse tree setting codewordsx
        if ((node = tree[current_code & 1][stack[depth]]) < 0)
        {
            do
            {
                if (depth + 1 == table_bits)
                    decoder_table[current_code] = node;
                current_code <<= 1;
                stack[++depth] = -node;
            } while ((node = tree[0][-node]) < 0);
        }

        codeword[node] = current_code; // set codeword to leaf symbol
        int32_t nb = codeword_length[node] = depth + 1;

        if (nb <= table_bits)
        {
            // add entries to decoder table
            if (nb == table_bits)
                decoder_table[current_code] = node;
            else
            {
                int32_t db = table_bits - nb, sc = current_code << db;
                for (int32_t k = 1 << db; k;)
                    decoder_table[sc + --k] = node;
            }
        }

        while (current_code & 1)
        { // backtrack
            current_code >>= 1;
            if (--depth < 0)
                break;
        }
        current_code |= 1;
    }
}

//----------------------------------------------------------------------------------------------------------------------------
static inline void hf_refill(hf_context *c) 
{
    if (c->bits_in_buffer <= 32) 
    {
        // Read 4 bytes at once
        uint32_t next = (uint32_t)c->bc_pointer[0] << 24 | (uint32_t)c->bc_pointer[1] << 16 |
                        (uint32_t)c->bc_pointer[2] << 8  | (uint32_t)c->bc_pointer[3];
        c->bit_buffer |= (uint64_t)next << (32 - c->bits_in_buffer);
        c->bc_pointer += 4;
        c->bits_in_buffer += 32;
    }
}

//----------------------------------------------------------------------------------------------------------------------------
void hf_init(hf_context *c, void *user_buffer, size_t length_in_bytes)
{
    assert(user_buffer != NULL);

    c->mode = mode_undefined;
    c->buffer_size = length_in_bytes;
    c->code_buffer = (uint8_t *)user_buffer;
    c->bc_pointer = (uint8_t *)user_buffer;
    c->bit_buffer = 0;
    c->bits_in_buffer = 0;
}

//----------------------------------------------------------------------------------------------------------------------------
void hf_start_encoder(hf_context *c)
{
    c->mode = mode_encoder;
    c->bit_buffer = 0;
    c->bits_in_buffer = 0;
    c->bc_pointer = c->code_buffer;
}

//----------------------------------------------------------------------------------------------------------------------------
size_t hf_stop_encoder(hf_context *c)
{
    assert(c->mode == mode_encoder);
    c->mode = mode_undefined;

    while (c->bits_in_buffer > 0)
    {
        *c->bc_pointer++ = (uint8_t)(c->bit_buffer >> 56);
        c->bit_buffer <<= 8;
        c->bits_in_buffer -= 8;
    }

    size_t code_bytes = (size_t)(c->bc_pointer - c->code_buffer);
    assert(code_bytes <= c->buffer_size);

    return code_bytes;
}

//----------------------------------------------------------------------------------------------------------------------------
void hf_start_decoder(hf_context *c)
{
    c->mode = mode_decoder;
    c->bit_buffer = 0;
    c->bits_in_buffer = 0;
    c->bc_pointer = c->code_buffer;
    
    // prime buffer
    hf_refill(c);
    hf_refill(c); 
}

//----------------------------------------------------------------------------------------------------------------------------
void hf_stop_decoder(hf_context *c)
{
    assert(c->mode == mode_decoder);
    c->mode = mode_undefined;
}

//----------------------------------------------------------------------------------------------------------------------------
void hf_encode(hf_context *c, const hf_model *model, uint32_t data) 
{
    uint32_t code = model->codeword[data];
    uint32_t len  = model->length[data];

    
    c->bit_buffer |= (uint64_t)code << (64 - c->bits_in_buffer - len);
    c->bits_in_buffer += len;

    // flush 32 bits at a time
    if (c->bits_in_buffer >= 32) 
    {
        uint32_t out = (uint32_t)(c->bit_buffer >> 32);

        c->bc_pointer[0] = (uint8_t)(out >> 24);
        c->bc_pointer[1] = (uint8_t)(out >> 16);
        c->bc_pointer[2] = (uint8_t)(out >> 8);
        c->bc_pointer[3] = (uint8_t)(out);
        c->bc_pointer += 4;
        c->bit_buffer <<= 32;
        c->bits_in_buffer -= 32;
    }
}

//----------------------------------------------------------------------------------------------------------------------------
uint32_t hf_decode(hf_context *c, const hf_model *model) 
{
    hf_refill(c);

    uint32_t idx = (uint32_t)(c->bit_buffer >> (64 - HF_TABLE_BITS));
    int32_t data = model->decd_table[idx];

    if (data >= 0) 
    {
        uint32_t len = model->length[data];
        c->bit_buffer <<= len;
        c->bits_in_buffer -= len;
    } 
    else 
    {
        c->bit_buffer <<= HF_TABLE_BITS;
        c->bits_in_buffer -= HF_TABLE_BITS;
        do 
        {
            data = model->tree[c->bit_buffer >> 63][-data];
            c->bit_buffer <<= 1;
            c->bits_in_buffer--;
        } while (data < 0);
    }
    return (uint32_t)data;
}

//----------------------------------------------------------------------------------------------------------------------------
void hf_put_bits(hf_context *c, uint32_t data, uint32_t len)
{
    assert(c->mode == mode_encoder);
    assert(len > 0 && len <= 32);

    c->bit_buffer |= (uint64_t)data << (64 - c->bits_in_buffer - len);
    c->bits_in_buffer += len;

    // Flush 32-bit chunks to memory
    if (c->bits_in_buffer >= 32)
    {
        uint32_t out = (uint32_t)(c->bit_buffer >> 32);

        c->bc_pointer[0] = (uint8_t)(out >> 24);
        c->bc_pointer[1] = (uint8_t)(out >> 16);
        c->bc_pointer[2] = (uint8_t)(out >>  8);
        c->bc_pointer[3] = (uint8_t)(out);
        c->bc_pointer += 4;
        c->bit_buffer <<= 32;
        c->bits_in_buffer -= 32;
    }
}

//----------------------------------------------------------------------------------------------------------------------------
uint32_t hf_get_bits(hf_context *c, uint32_t len)
{
    assert(c->mode == mode_decoder);
    assert(len > 0 && len <= 32);

    if (c->bits_in_buffer < (int32_t)len) 
        hf_refill(c); 

    uint32_t data = (uint32_t)(c->bit_buffer >> (64 - len));
    c->bit_buffer <<= len;
    c->bits_in_buffer -= len;

    return data;
}

//----------------------------------------------------------------------------------------------------------------------------
void hf_model_init(hf_model *model, uint32_t num_symbols, uint32_t *symbol_count)
{
    assert(num_symbols >= 3 && num_symbols <= UINT8_MAX + 1);

    model->table_bits = HF_TABLE_BITS;
    model->table_shift = 32 - model->table_bits;
    model->data_symbols = num_symbols;

    if (symbol_count != NULL)
    {
        for (uint32_t n = 0; n < model->data_symbols; n++)
            model->symbol_count[n] = symbol_count[n];

        // construct optimal code
        form_huffman_tree(model->data_symbols, model->symbol_count, model->tree[0], model->tree[1]);
        set_huffman_code(model->table_bits, model->tree, model->decd_table, model->codeword, model->length);
    }
}

//----------------------------------------------------------------------------------------------------------------------------
void hf_save_model(hf_context *c, const hf_model *model)
{
    hf_put_bits(c, model->data_symbols - 1, 8);

    uint32_t max_val = 0;
    for (uint32_t i = 0; i < model->data_symbols; i++)
    {
        if (model->symbol_count[i] > max_val)
            max_val = model->symbol_count[i];
    }

    uint32_t bits_needed = 0;
    while ((1U << bits_needed) <= max_val && bits_needed < 32)
        bits_needed++;
    if (bits_needed == 0)
        bits_needed = 1;

    hf_put_bits(c, bits_needed, 6);

    for (uint32_t i = 0; i < model->data_symbols; i++)
        hf_put_bits(c, model->symbol_count[i], bits_needed);
}

//----------------------------------------------------------------------------------------------------------------------------
void hf_load_model(hf_context *c, hf_model *model)
{
    uint32_t num_symbols = hf_get_bits(c, 8) + 1;
    uint32_t bits_per_count = hf_get_bits(c, 6);

    hf_model_init(model, num_symbols, NULL);

    for (uint32_t i = 0; i < num_symbols; i++)
        model->symbol_count[i] = hf_get_bits(c, bits_per_count);

    form_huffman_tree(model->data_symbols, model->symbol_count, model->tree[0], model->tree[1]);
    set_huffman_code(model->table_bits, model->tree, model->decd_table, model->codeword, model->length);
}

//----------------------------------------------------------------------------------------------------------------------
#define HISTOGRAM_SIZE (256)

typedef struct histogram
{
    uint32_t count[HISTOGRAM_SIZE];
    uint32_t num_symbols;
} histogram;

//----------------------------------------------------------------------------------------------------------------------------
void histogram_init(histogram* h, uint32_t num_symbols)
{
    assert(num_symbols > 3 && num_symbols <= 256);
    h->num_symbols = num_symbols;
    for(uint32_t i=0; i<num_symbols; ++i)
        h->count[i] = 0;
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
static inline int int_abs(int a) {return (a>=0) ? a : -a;}

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
void bc1_crunch(hf_context* codec, void* cruncher_memory, const void* input, size_t stride, uint32_t width, uint32_t height)
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
    histogram hred, hgreen, hblue;
    histogram_init(&hred, 1 << COLOR_DELTA_NUM_BITS);
    histogram_init(&hgreen, 1 << COLOR_DELTA_NUM_BITS);
    histogram_init(&hblue, 1 << COLOR_DELTA_NUM_BITS);

    histogram htable_index, hmask, htable_difference;
    histogram_init(&htable_index, 256);
    histogram_init(&hmask, 16);
    histogram_init(&htable_difference, 256);

    histogram htable_entry;
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

                assert((dgreen + 64) < 256);
                hgreen.count[(dgreen + 64)]++;
                
                dgreen /= 2;
                dred -= dgreen;
                dblue -= dgreen;

                assert((dred + 64) < 256);
                assert((dblue + 64) < 256);
                hred.count[dred + 64]++;
                hblue.count[dblue + 64]++;
            }

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

    hf_model red, green, blue;
    hf_model_init(&red, hred.num_symbols, hred.count);
    hf_model_init(&green, hred.num_symbols, hgreen.count);
    hf_model_init(&blue, hred.num_symbols, hblue.count);

    hf_save_model(codec, &red);
    hf_save_model(codec, &green);
    hf_save_model(codec, &blue);

    hf_model table_index, diff_mask, table_entry, table_difference;
    hf_model_init(&table_index, htable_index.num_symbols, htable_index.count);
    hf_model_init(&diff_mask, hmask.num_symbols, hmask.count);
    hf_model_init(&table_entry, htable_entry.num_symbols, htable_entry.count);
    hf_model_init(&table_difference, htable_difference.num_symbols, htable_difference.count);

    hf_save_model(codec, &table_index);
    hf_save_model(codec, &diff_mask);
    hf_save_model(codec, &table_entry);
    hf_save_model(codec, &table_difference);

    // write the top-table
    hf_put_bits(codec, top_table_size-1, TABLE_INDEX_NUM_BITS);   // entries count
    for(uint32_t j=0; j<4; ++j)
        hf_put_bits(codec, (top_table[0] >> (j*8)) & 0xff, 8);    // first entry not compressed

    for(uint32_t i=1; i<top_table_size; ++i)
    {
        // table is sorted from small to big, so diff is always positive
        uint32_t diff = top_table[i] - top_table[i-1];

        for(uint32_t j=0; j<4; ++j)
            hf_encode(codec, &table_entry, (diff >> (j*8)) & 0xff);
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

                    int previous_delta = int_abs(current_red - previous_red) + int_abs(current_green-previous_green) + int_abs(current_blue-previous_blue);
                    int up_delta = int_abs(current_red-up_red) + int_abs(current_green-up_green) + int_abs(current_blue-up_blue);

                    hf_put_bits(codec, (up_delta < previous_delta) ? 1 : 0, 1);

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
                hf_encode(codec, &green, (uint8_t) (dgreen + 64));

                // then encode red and blue delta based on green delta
                // assuming some relation between green and other components
                dgreen /= 2;
                dred -= dgreen;
                dblue -= dgreen;

                hf_encode(codec, &red, (uint8_t) (dred + 64));
                hf_encode(codec, &blue, (uint8_t) (dblue + 64));
            }

            // for indices, we store the reference to "nearest" indices (can be exactly the same)
            // and the delta with this reference
            uint32_t reference = nearest32(top_table, top_table_size, current->indices) & 0xffff;
            hf_encode(codec, &table_index, reference);

            // xor the difference and encode (could be 0 if equal to reference)
            uint32_t difference = current->indices ^ top_table[reference];

            uint32_t mask = 0;
            if ((difference & 0x000000FF) != 0) mask |= 1;
            if ((difference & 0x0000FF00) != 0) mask |= 2;
            if ((difference & 0x00FF0000) != 0) mask |= 4;
            if ((difference & 0xFF000000) != 0) mask |= 8;

            hf_encode(codec, &diff_mask, mask);

            // only encode the bytes that are actually non-zero
            for(uint32_t j=0; j<4; ++j)
                if (mask & (1u << j))
                    hf_encode(codec, &table_difference, (difference >> (j*8)) & 0xff);

            previous = *current;
        }
    }
}

//----------------------------------------------------------------------------------------------------------------------------
void bc1_decrunch(hf_context* codec, uint32_t width, uint32_t height, void* output, size_t stride)
{
    assert((width % 4 == 0) && (height % 4 == 0));

    uint32_t height_blocks = height/4;
    uint32_t width_blocks = width/4;

    hf_model red, green, blue;
    hf_load_model(codec, &red);
    hf_load_model(codec, &green);
    hf_load_model(codec, &blue);
    
    hf_model table_index, diff_mask, table_entry, table_difference;
    hf_load_model(codec, &table_index);
    hf_load_model(codec, &diff_mask);
    hf_load_model(codec, &table_entry);
    hf_load_model(codec, &table_difference);

    uint32_t top_table[TABLE_SIZE];
    uint32_t top_table_size = hf_get_bits(codec, TABLE_INDEX_NUM_BITS)+1;

    top_table[0] = 0;
    for(uint32_t j=0; j<4; ++j)
        top_table[0] |= hf_get_bits(codec, 8) << (j*8);

    for(uint32_t i=1; i<top_table_size; ++i)
    {
        uint32_t diff = 0;
        for(uint32_t j=0; j<4; ++j)
            diff |= hf_decode(codec, &table_entry) << (j*8);

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
                if (y>0 && x!=0 && hf_get_bits(codec, 1))
                {
                    bc1_block* up = (bc1_block*) get_block(output, stride, width_blocks, zigzag_x, y-1);
                    bc1_extract_565(up->color[j], &reference_red, &reference_green, &reference_blue);
                }

                uint8_t delta_green = (uint8_t)hf_decode(codec, &green);
                uint8_t delta_red = (uint8_t)hf_decode(codec, &red);
                uint8_t delta_blue = (uint8_t)hf_decode(codec, &blue);

                // red and blue delta are based on green delta
                int dgreen_orig = (int)delta_green - 64;
                int current_green_value = reference_green + dgreen_orig;
                int dgreen_halved = dgreen_orig / 2;

                int dred_orig = ((int)delta_red - 64) + dgreen_halved;
                int dblue_orig = ((int)delta_blue - 64) + dgreen_halved;

                int current_red_value = reference_red + dred_orig;
                int current_blue_value = reference_blue + dblue_orig;

                current->color[j] = bc1_pack_565((uint8_t)current_red_value, (uint8_t)current_green_value, (uint8_t)current_blue_value);
            }

            // indices difference with top table
            uint32_t reference = hf_decode(codec, &table_index);
            uint32_t mask = hf_decode(codec, &diff_mask);

            uint32_t difference=0;
            for(uint32_t j=0; j<4; ++j)
                if (mask & (1 << j))
                    difference = difference | (hf_decode(codec, &table_difference) << (j*8));

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

    hf_context codec;
    hf_init(&codec, output, length);
    hf_start_encoder(&codec);

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
    return hf_stop_encoder(&codec);
}

//----------------------------------------------------------------------------------------------------------------------------
void bc_decrunch(const void* input, size_t length, uint32_t width, uint32_t height, enum bc_format format, void* output)
{
    hf_context codec;
    hf_init(&codec, (void*)input, length);
    hf_start_decoder(&codec);

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

    hf_stop_decoder(&codec);
}
