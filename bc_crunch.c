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

#if defined(__ARM_NEON) && defined(__ARM_NEON__)
#include <arm_neon.h>
#define BC_CRUNCH_NEON
#endif

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
typedef struct range_model
{
    uint32_t distribution[2 * DM_MAX_SYMBOLS + DM_MAX_TABLE_SIZE + 2]; 
    uint32_t* symbol_count;
    uint32_t total_count, update_cycle, symbols_until_update;
    uint32_t data_symbols, last_symbol, table_size, table_shift;
    uint32_t *decoder_table;
} range_model;

//----------------------------------------------------------------------------------------------------------------------
static inline int int_abs(int a) {return (a>=0) ? a : -a;}

//----------------------------------------------------------------------------------------------------------------------
void adaptive_model_update(range_model* model)
{
    if ((model->total_count += model->update_cycle) > DM__MaxCount) 
    {
        model->total_count = 0;
        for (uint32_t n = 0; n < model->data_symbols; n++)
            model->total_count += (model->symbol_count[n] = (model->symbol_count[n] + 1) >> 1);
    }

    uint32_t k, sum = 0;
    uint32_t scale = 0x80000000U / model->total_count;

    if (model->table_size == 0)
    {
        for (k = 0; k < model->data_symbols; k++) 
        {
            model->distribution[k] = (scale * sum) >> (31 - DM__LengthShift);
            sum += model->symbol_count[k];
        }
    }
    else
    {
        uint32_t k, sum = 0, s = 0;
        for (k = 0; k < model->data_symbols; k++) 
        {
            model->distribution[k] = (scale * sum) >> (31 - DM__LengthShift);
            sum += model->symbol_count[k];
            uint32_t w = model->distribution[k] >> model->table_shift;
            while (s < w) model->decoder_table[++s] = k - 1;
        }
        model->decoder_table[0] = 0;
        while (s <= model->table_size) model->decoder_table[++s] = model->data_symbols - 1;
    }
    

    model->update_cycle = (5 * model->update_cycle) >> 2;
    uint32_t max_cycle = (model->data_symbols + 6) << 3;
    if (model->update_cycle > max_cycle) 
        model->update_cycle = max_cycle;
    model->symbols_until_update = model->update_cycle;
}

//----------------------------------------------------------------------------------------------------------------------
void model_init(range_model* model, uint32_t number_of_symbols)
{
    model->data_symbols = number_of_symbols;
    model->last_symbol = model->data_symbols - 1;
    model->total_count = 0;
    model->update_cycle = model->data_symbols;

    // define size of table for fast decoding
    if (model->data_symbols > 16) 
    {
        uint32_t table_bits = 3;
        while (model->data_symbols > (1U << (table_bits + 2))) ++table_bits;
        model->table_size  = 1 << table_bits;
        model->table_shift = DM__LengthShift - table_bits;
        model->decoder_table = model->distribution + 2 * model->data_symbols;
        assert(model->distribution != NULL);
    }
    else
    {
        // small alphabet: no table needed
        model->decoder_table = 0;
        model->table_size = model->table_shift = 0;
    }

    model->symbol_count = model->distribution + model->data_symbols;

    for (uint32_t k = 0; k < model->data_symbols; k++) 
        model->symbol_count[k] = 1;

    adaptive_model_update(model);
    model->symbols_until_update = model->update_cycle = (model->data_symbols + 6) >> 1;
}

//----------------------------------------------------------------------------------------------------------------------
typedef struct range_codec
{
    uint8_t *code_buffer, *ac_pointer;
    uint32_t base, value, length;
    uint32_t buffer_size;
} range_codec;

//----------------------------------------------------------------------------------------------------------------------
inline static void propagate_carry(range_codec* codec)
{
    uint8_t * p;            
    // carry propagation on compressed data buffer
    for (p = codec->ac_pointer - 1; *p == 0xFFU; p--) 
        *p = 0;
    ++*p;
}

//----------------------------------------------------------------------------------------------------------------------
inline static void renorm_enc_interval(range_codec* codec)
{
    do  // output and discard top byte
    {
        *codec->ac_pointer++ = (uint8_t)(codec->base >> 24);
        codec->base <<= 8;
    } while ((codec->length <<= 8) < RC__MinLength);        // length multiplied by 256
}

//----------------------------------------------------------------------------------------------------------------------
void enc_init(range_codec* codec, uint8_t* output, uint32_t length)
{
    codec->buffer_size = length;
    codec->code_buffer = output;
    codec->base = 0;
    codec->length = RC__MaxLength;
    codec->ac_pointer = codec->code_buffer;
}

//----------------------------------------------------------------------------------------------------------------------
void dec_init(range_codec* codec, const uint8_t* input, uint32_t length)
{
    codec->buffer_size = length;
    codec->code_buffer = (uint8_t*) input;
    codec->length = RC__MaxLength;
    codec->ac_pointer = codec->code_buffer + 3;
    codec->value = ((uint32_t)(codec->code_buffer[0]) << 24) |
                   ((uint32_t)(codec->code_buffer[1]) << 16) |
                   ((uint32_t)(codec->code_buffer[2]) <<  8) |
                    (uint32_t)(codec->code_buffer[3]);
}

//----------------------------------------------------------------------------------------------------------------------
uint32_t enc_done(range_codec* codec)
{
    uint32_t init_base = codec->base;

    if (codec->length > 2 * RC__MinLength) 
    {
        codec->base  += RC__MinLength;
        codec->length = RC__MinLength >> 1;
    }
    else 
    {
        codec->base  += RC__MinLength >> 1;
        codec->length = RC__MinLength >> 9;
    }

    if (init_base > codec->base) 
        propagate_carry(codec);

    renorm_enc_interval(codec);

    uint32_t code_bytes = (uint32_t)(codec->ac_pointer - codec->code_buffer);
    assert(code_bytes <= codec->buffer_size);

    return code_bytes;
}

//----------------------------------------------------------------------------------------------------------------------
void enc_put(range_codec* codec, range_model* model, uint32_t data)
{
    assert(data < model->data_symbols); // invalid data symbols

    uint32_t x;
    uint32_t init_base = codec->base;

    if (data == model->last_symbol) 
    {
        x = model->distribution[data] * (codec->length >> DM__LengthShift);
        codec->base += x;
        codec->length -= x;
    }
    else 
    {
        x = model->distribution[data] * (codec->length >>= DM__LengthShift);
        codec->base += x;
        codec->length  = model->distribution[data+1] * codec->length - x;
    }

    if (init_base > codec->base) 
        propagate_carry(codec);

    if (codec->length < RC__MinLength) 
        renorm_enc_interval(codec);

    ++model->symbol_count[data];
    if (--model->symbols_until_update == 0)
        adaptive_model_update(model);
}

//----------------------------------------------------------------------------------------------------------------------
static inline void enc_put_bits(range_codec* codec, uint32_t data, uint32_t number_of_bits)
{
    assert(data < DM_MAX_SYMBOLS);

    uint32_t init_base = codec->base;
    codec->base += data * (codec->length >>= number_of_bits);            // new interval base and length

    if (init_base > codec->base) 
        propagate_carry(codec);

    if (codec->length < RC__MinLength) 
        renorm_enc_interval(codec);
}

//----------------------------------------------------------------------------------------------------------------------
uint32_t dec_get(range_codec* codec, range_model* model)
{
    uint32_t n, s, x, y = codec->length;

    if (model->decoder_table) 
    {
        // use table look-up for faster decoding
        uint32_t dv = codec->value / (codec->length >>= DM__LengthShift);
        uint32_t t = dv >> model->table_shift;

        s = model->decoder_table[t];         // initial decision based on table look-up
        n = model->decoder_table[t+1] + 1;

        while (n > s + 1) 
        {                        // finish with bisection search
            uint32_t m = (s + n) >> 1;
            if (model->distribution[m] > dv) 
                n = m; 
            else s = m;
        }

        // compute products
        x = model->distribution[s] * codec->length;
        if (s != model->last_symbol) 
            y = model->distribution[s+1] * codec->length;
    }
    else
    {
        // decode using only multiplications
        x = s = 0;
        codec->length >>= DM__LengthShift;
        uint32_t m = (n = model->data_symbols) >> 1;

        // decode via bisection search
        do 
        {
            uint32_t z = codec->length * model->distribution[m];
            if (z > codec->value) 
            {
                n = m;
                y = z;
            }
            else 
            {
                s = m;
                x = z;
            }
        } while ((m = (s + n) >> 1) != s);
    }

    codec->value -= x;
    codec->length = y - x;

    if (codec->length < RC__MinLength)
    {
        do
        {
            codec->value = (codec->value << 8) | (uint32_t)(*++codec->ac_pointer);
        } while ((codec->length <<= 8) < RC__MinLength);
    }

    ++model->symbol_count[s];
    if (--model->symbols_until_update == 0) 
        adaptive_model_update(model);

    return s;
}

//----------------------------------------------------------------------------------------------------------------------
uint32_t dec_get_bits(range_codec* codec, uint32_t number_of_bits)
{
    unsigned s = codec->value / (codec->length >>= number_of_bits);      
    codec->value -= codec->length * s;
    if (codec->length < RC__MinLength)
    {
        do
        {
            codec->value = (codec->value << 8) | (uint32_t)(*++codec->ac_pointer);
        } while ((codec->length <<= 8) < RC__MinLength);
    }

    return s;
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

    // vector quantization
    vq_top_table(input, stride, num_blocks, output, num_entries);

    // sort table for compression
    qsort(output, *num_entries, sizeof(uint32_t), compare_entries);
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
static inline range_model* bc4_select_model(const bc4_block* b, range_model* indices)
{
    int endpoints_delta = int_abs(b->color[0] - b->color[1]);

    if (endpoints_delta < 8)
        return &indices[0];
    else if (endpoints_delta < 32)
        return &indices[8];

    return &indices[16];
}

//----------------------------------------------------------------------------------------------------------------------------
void bc1_crunch(range_codec* codec, void* cruncher_memory, const void* input, size_t stride, uint32_t width, uint32_t height)
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
    range_model table_entry;
    model_init(&table_entry, 256);

    enc_put_bits(codec, top_table_size-1, TABLE_INDEX_NUM_BITS);   // entries count
    for(uint32_t j=0; j<4; ++j)
        enc_put_bits(codec, (top_table[0] >> (j*8)) & 0xff, 8);    // first entry not compressed

    for(uint32_t i=1; i<top_table_size; ++i)
    {
        // table is sorted from small to big, so diff is always positive
        uint32_t diff = top_table[i] - top_table[i-1];

        for(uint32_t j=0; j<4; ++j)
            enc_put(codec, &table_entry, (diff >> (j*8)) & 0xff);
    }

    range_model red, green, blue;
    model_init(&red,   1 << COLOR_DELTA_NUM_BITS);
    model_init(&green, 1 << COLOR_DELTA_NUM_BITS);
    model_init(&blue,  1 << COLOR_DELTA_NUM_BITS);

    range_model table_index, table_difference, diff_mask, color_reference;
    model_init(&table_index, top_table_size);
    model_init(&table_difference, 256);
    model_init(&diff_mask, 16);
    model_init(&color_reference, 2);

    bc1_block empty_block = {0};
    const bc1_block* previous = &empty_block;

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
                bc1_extract_565(previous->color[j], &previous_red, &previous_green, &previous_blue);

                if (y>0 && x!=0)
                {
                    const bc1_block* up = get_block(input, stride, width_blocks, zigzag_x, y-1);
                    uint8_t up_red, up_green, up_blue;
                    bc1_extract_565(up->color[j], &up_red, &up_green, &up_blue);

                    int previous_delta = int_abs(current_red - previous_red) + int_abs(current_green-previous_green) + int_abs(current_blue-previous_blue);
                    int up_delta = int_abs(current_red-up_red) + int_abs(current_green-up_green) + int_abs(current_blue-up_blue);

                    enc_put(codec, &color_reference, (up_delta < previous_delta) ? 1 : 0);

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
                enc_put(codec, &green, (uint8_t) (dgreen + 64));

                // then encode red and blue delta based on green delta
                // assuming some relation between green and other components
                dgreen /= 2;
                dred -= dgreen;
                dblue -= dgreen;

                enc_put(codec, &red, (uint8_t) (dred + 64));
                enc_put(codec, &blue, (uint8_t) (dblue + 64));
            }

            // for indices, we store the reference to "nearest" indices (can be exactly the same)
            // and the delta with this reference
            uint32_t reference = nearest32(top_table, top_table_size, current->indices) & 0xffff;
            enc_put(codec, &table_index, reference);

            // xor the difference and encode (could be 0 if equal to reference)
            uint32_t difference = current->indices ^ top_table[reference];

            uint32_t mask = 0;
            if ((difference & 0x000000FF) != 0) mask |= 1;
            if ((difference & 0x0000FF00) != 0) mask |= 2;
            if ((difference & 0x00FF0000) != 0) mask |= 4;
            if ((difference & 0xFF000000) != 0) mask |= 8;

            enc_put(codec, &diff_mask, mask);

            // only encode the bytes that are actually non-zero
            for(uint32_t j=0; j<4; ++j)
                if (mask & (1u << j))
                    enc_put(codec, &table_difference, (difference >> (j*8)) & 0xff);

            previous = current;
        }
    }
}

//----------------------------------------------------------------------------------------------------------------------------
void bc1_decrunch(range_codec* codec, uint32_t width, uint32_t height, void* output, size_t stride)
{
    assert((width % 4 == 0) && (height % 4 == 0));

    uint32_t height_blocks = height/4;
    uint32_t width_blocks = width/4;

    range_model red, green, blue;
    model_init(&red,   1 << COLOR_DELTA_NUM_BITS);
    model_init(&green, 1 << COLOR_DELTA_NUM_BITS);
    model_init(&blue,  1 << COLOR_DELTA_NUM_BITS);

    range_model table_entry;
    model_init(&table_entry, 256);

    uint32_t top_table[TABLE_SIZE];
    uint32_t top_table_size = dec_get_bits(codec, TABLE_INDEX_NUM_BITS)+1;

    top_table[0] = 0;
    for(uint32_t j=0; j<4; ++j)
        top_table[0] |= dec_get_bits(codec, 8) << (j*8);

    for(uint32_t i=1; i<top_table_size; ++i)
    {
        uint32_t diff = 0;
        for(uint32_t j=0; j<4; ++j)
            diff |= dec_get(codec, &table_entry) << (j*8);

        top_table[i] = top_table[i-1] + diff;
    }

    range_model table_index, table_difference, diff_mask, color_reference;
    model_init(&table_index, top_table_size);
    model_init(&table_difference, 256);
    model_init(&diff_mask, 16);
    model_init(&color_reference, 2);

    bc1_block empty_block = {0};
    bc1_block* previous = &empty_block;

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
                bc1_extract_565(previous->color[j], &reference_red, &reference_green, &reference_blue);
                if (y>0 && x!=0 && dec_get(codec, &color_reference))
                {
                    bc1_block* up = (bc1_block*) get_block(output, stride, width_blocks, zigzag_x, y-1);
                    bc1_extract_565(up->color[j], &reference_red, &reference_green, &reference_blue);
                }

                uint8_t delta_green = (uint8_t)dec_get(codec, &green);
                uint8_t delta_red = (uint8_t)dec_get(codec, &red);
                uint8_t delta_blue = (uint8_t)dec_get(codec, &blue);

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
            uint32_t reference = dec_get(codec, &table_index);
            uint32_t mask = dec_get(codec, &diff_mask);

            uint32_t difference=0;
            for(uint32_t j=0; j<4; ++j)
                if (mask & (1 << j))
                    difference = difference | (dec_get(codec, &table_difference) << (j*8));

            current->indices =  difference ^ top_table[reference];

            previous = current;
        }
    }
}

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

    range_codec codec;
    enc_init(&codec, (uint8_t*) output, (uint32_t)length);

    switch(format)
    {
    case bc1 : 
        {
            bc1_crunch(&codec, cruncher_memory, input, sizeof(bc1_block), width, height);
            break;
        }
    case bc3 : 
        {
            size_t block_size = sizeof(bc1_block) + sizeof(bc4_block);
            bc4_crunch(&codec, cruncher_memory, input, block_size, width, height);
            bc1_crunch(&codec, cruncher_memory, ptr_shift_const(input, sizeof(bc4_block)), block_size, width, height);
            break;
        }
    case bc4 : 
        {
            bc4_crunch(&codec, cruncher_memory, input, sizeof(bc4_block), width, height);
            break;
        }
    case bc5 :
        {
            size_t block_size = sizeof(bc4_block) * 2;
            bc4_crunch(&codec, cruncher_memory, input, block_size, width, height);
            bc4_crunch(&codec, cruncher_memory, ptr_shift_const(input, sizeof(bc4_block)), block_size, width, height);
            break;
        }

    default: break;
    }
    return enc_done(&codec);
}

//----------------------------------------------------------------------------------------------------------------------------
void bc_decrunch(const void* input, size_t length, uint32_t width, uint32_t height, enum bc_format format, void* output)
{
    range_codec codec;
    dec_init(&codec, (const uint8_t*)input, (uint32_t)length);

    switch(format)
    {
    case bc1 : bc1_decrunch(&codec, width, height, output, sizeof(bc1_block)); break;
    case bc3 : 
        {
            size_t block_size = sizeof(bc1_block) + sizeof(bc4_block);
            bc4_decrunch(&codec, width, height, output, block_size);
            bc1_decrunch(&codec, width, height, ptr_shift(output, sizeof(bc4_block)), block_size);
            break;
        }
    case bc4 : bc4_decrunch(&codec, width, height, output, sizeof(bc4_block)); break;
    case bc5 :
        {
            size_t block_size = sizeof(bc4_block) * 2;
            bc4_decrunch(&codec, width, height, output, block_size);
            bc4_decrunch(&codec, width, height, ptr_shift(output, sizeof(bc4_block)), block_size);
            break;
        }
    default: break;
    }
}