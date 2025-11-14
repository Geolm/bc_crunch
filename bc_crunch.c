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

#include <stdio.h>      // TODO : not needed remove it!


//----------------------------------------------------------------------------------------------------------------------------
// Private structures & functions
//----------------------------------------------------------------------------------------------------------------------------

#define RC__MinLength (0x01000000U)
#define RC__MaxLength (0xFFFFFFFFU)
#define MAX_ALPHABET_SIZE (256)
#define DM__LengthShift (15)
#define DM__MaxCount    (1 << DM__LengthShift)
#define HASHMAP_SIZE (1 << 20)
#define TABLE_INDEX_NUM_BITS (8)
#define TABLE_SIZE (1<<TABLE_INDEX_NUM_BITS)
#define RED_DELTA_NUM_BITS (6)
#define GREEN_DELTA_NUM_BITS (7)
#define BLUE_DELTA_NUM_BITS (6)
#define DICTIONARY_SIZE (256)
#define MAKE48(r0,r1,r2) ( ((uint64_t)(r0) << 32) | ((uint64_t)(r1) << 16) | (uint64_t)(r2) )
#define BC4_COLOR_NUM_BITS (8)

#if defined(_MSC_VER)
    #include <intrin.h>
    #pragma intrinsic(__popcnt)
    #define popcount(x) __popcnt(x)
#else
    #define popcount(x) __builtin_popcount(x)
#endif

//----------------------------------------------------------------------------------------------------------------------
typedef struct range_model
{
    uint32_t distribution[MAX_ALPHABET_SIZE]; 
    uint32_t symbol_count[MAX_ALPHABET_SIZE];
    uint32_t total_count, update_cycle, symbols_until_update;
    uint32_t data_symbols, last_symbol;
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

    for (k = 0; k < model->data_symbols; k++) 
    {
        model->distribution[k] = (scale * sum) >> (31 - DM__LengthShift);
        sum += model->symbol_count[k];
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
    assert(model->distribution != NULL); // adaptive model should be initialized
    
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
    assert(data < MAX_ALPHABET_SIZE);

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
    assert(model->distribution != NULL);

    uint32_t n, s, x, y = codec->length;

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
static inline uint8_t delta_encode(uint8_t prev, uint8_t curr, uint8_t num_bits)
{
    int delta = (int)curr - (int)prev;
    return (uint8_t)(delta + (1 << (num_bits - 1)));
}

//----------------------------------------------------------------------------------------------------------------------------
static inline uint8_t delta_decode(uint8_t prev, uint8_t delta_encoded, uint8_t num_bits)
{
    int delta = (int)delta_encoded - (1 << (num_bits - 1));
    int curr = (int)prev + delta;
    return (uint8_t)curr;
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
void build_top_table(entry* hashmap, const void* input, size_t stride, uint32_t num_blocks, entry* output, uint32_t* num_entries)
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
        output[i] = table[TABLE_SIZE - i - 1];
        if (output[i].count>0)
            (*num_entries)++;
    }
}

//----------------------------------------------------------------------------------------------------------------------------
static inline uint32_t top_table_nearest(const entry* table, uint32_t table_size, uint32_t original_indices)
{
    uint32_t best_score = 32;
    uint32_t best_index = 0;
    for(uint32_t i=0; i<table_size; ++i)
    {
        uint32_t delta = table[i].key ^ original_indices;

        // early exit if identical
        if (delta == 0)
            return i;

        uint32_t score = popcount(delta);
        if (score < best_score ||
           ((score == best_score) && (table[i].key > table[best_index].key)))
        {
            best_score = score;
            best_index = i;
        }
    }
    return best_index;
}

//----------------------------------------------------------------------------------------------------------------------------
int compare_entries(const void* a, const void* b)
{
    const entry* entry_a = (const entry*) a;
    const entry* entry_b = (const entry*) b;

    if (entry_a->key < entry_b->key)
        return -1;

    if (entry_a->key > entry_b->key)
        return 1;

    return 0;
}

//----------------------------------------------------------------------------------------------------------------------------
void sort_top_table(entry* table, uint32_t table_size)
{
    qsort(table, table_size, sizeof(entry), compare_entries);
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
static inline void bc4_set_index(bc4_block* b, uint32_t pixel_index, uint8_t index)
{
    uint32_t bit_offset = pixel_index * 3;
    uint32_t word_index = bit_offset >> 4;
    uint32_t bit_in_word = bit_offset & 0xF;

    uint16_t mask = 0x7 << bit_in_word;
    b->indices[word_index] = (b->indices[word_index] & ~mask) | ((index & 0x7) << bit_in_word);

    // if the 3 bits spill into the next word
    if (bit_in_word > 13)  // last 2 or 1 bits spill
    {
        uint16_t spill_bits = (index & 0x7) >> (16 - bit_in_word);
        b->indices[word_index + 1] = (b->indices[word_index + 1] & ~(0x7 >> (16 - bit_in_word))) | spill_bits;
    }
}

//----------------------------------------------------------------------------------------------------------------------------
static const uint32_t morton_inverse[16] = {0, 1, 4, 5, 2, 3, 6, 7, 8, 9, 12, 13, 10, 11, 14, 15};

//----------------------------------------------------------------------------------------------------------------------------
void bc1_crunch(range_codec* codec, void* cruncher_memory, const void* input, size_t stride, uint32_t width, uint32_t height)
{
    assert((width%4 == 0) && (height%4 == 0));
    assert(((uintptr_t)cruncher_memory)%sizeof(uintptr_t) == 0);

    uint32_t height_blocks = height/4;
    uint32_t width_blocks = width/4;

    // build a histogram and select the TABLE_SIZE block indices which are most used
    entry* hashmap = (entry*) cruncher_memory;
    entry top_table[TABLE_SIZE];
    uint32_t top_table_size;
    build_top_table(hashmap, input, stride, height_blocks*width_blocks, top_table, &top_table_size);
    sort_top_table(top_table, top_table_size);

    // write the table
    range_model table_entry;
    model_init(&table_entry, 256);

    enc_put_bits(codec, top_table_size-1, TABLE_INDEX_NUM_BITS);   // entries count
    for(uint32_t j=0; j<4; ++j)
        enc_put_bits(codec, (top_table[0].key >> (j*8)) & 0xff, 8);    // first entry not compressed

    for(uint32_t i=1; i<top_table_size; ++i)
    {
        // table is sorted from small to big, so diff is always positive
        uint32_t diff = top_table[i].key - top_table[i-1].key;

        for(uint32_t j=0; j<4; ++j)
            enc_put(codec, &table_entry, (diff >> (j*8)) & 0xff);
    }

    range_model red, green, blue;
    model_init(&red,   1 << RED_DELTA_NUM_BITS);
    model_init(&green, 1 << GREEN_DELTA_NUM_BITS);
    model_init(&blue,  1 << BLUE_DELTA_NUM_BITS);

    range_model table_index, table_difference, block_mode, color_reference;
    model_init(&table_index, top_table_size);
    model_init(&table_difference, 256);
    model_init(&block_mode, 2);
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

                if (y>0)
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

                enc_put(codec, &red, delta_encode(previous_red, current_red, RED_DELTA_NUM_BITS));
                enc_put(codec, &green, delta_encode(previous_green, current_green, GREEN_DELTA_NUM_BITS));
                enc_put(codec, &blue, delta_encode(previous_blue, current_blue, BLUE_DELTA_NUM_BITS));
            }

            // for indices, we store the reference to "nearest" indices (can be exactly the same)
            // and the delta with this reference
            uint32_t reference = top_table_nearest(top_table, top_table_size, current->indices);
            enc_put(codec, &table_index, reference);

            // xor the difference and encode (could be 0 if equal to reference)
            uint32_t difference = current->indices ^ top_table[reference].key;

            if (difference == 0)
                enc_put(codec, &block_mode, 0);
            else
            {
                enc_put(codec, &block_mode, 1);

                // encode the full 32 bits difference
                for(uint32_t j=0; j<4; ++j)
                    enc_put(codec, &table_difference, (difference>>(j*8))&0xff);
            }

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
    model_init(&red,   1 << RED_DELTA_NUM_BITS);
    model_init(&green, 1 << GREEN_DELTA_NUM_BITS);
    model_init(&blue,  1 << BLUE_DELTA_NUM_BITS);

    range_model table_entry;
    model_init(&table_entry, 256);

    entry top_table[TABLE_SIZE];
    uint32_t top_table_size = dec_get_bits(codec, TABLE_INDEX_NUM_BITS)+1;

    top_table[0].key = 0;
    for(uint32_t j=0; j<4; ++j)
        top_table[0].key |= dec_get_bits(codec, 8) << (j*8);

    for(uint32_t i=1; i<top_table_size; ++i)
    {
        uint32_t diff = 0;
        for(uint32_t j=0; j<4; ++j)
            diff |= dec_get(codec, &table_entry) << (j*8);

        top_table[i].key = top_table[i-1].key + diff;
    }

    range_model table_index, table_difference, block_mode, color_reference;
    model_init(&table_index, top_table_size);
    model_init(&table_difference, 256);
    model_init(&block_mode, 2);
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
                if (y>0 && dec_get(codec, &color_reference))
                {
                    bc1_block* up = (bc1_block*) get_block(output, stride, width_blocks, zigzag_x, y-1);
                    bc1_extract_565(up->color[j], &reference_red, &reference_green, &reference_blue);
                }

                uint8_t delta_red = (uint8_t)dec_get(codec, &red);
                uint8_t delta_green = (uint8_t)dec_get(codec, &green);
                uint8_t delta_blue = (uint8_t)dec_get(codec, &blue);

                uint8_t current_red = delta_decode(reference_red, delta_red, RED_DELTA_NUM_BITS);
                uint8_t current_green = delta_decode(reference_green, delta_green, GREEN_DELTA_NUM_BITS);
                uint8_t current_blue = delta_decode(reference_blue, delta_blue, BLUE_DELTA_NUM_BITS);

                current->color[j] = bc1_pack_565(current_red, current_green, current_blue);
            }

            // indices difference with top table
            uint32_t reference = dec_get(codec, &table_index);
            uint32_t mode = dec_get(codec, &block_mode);

            if (mode == 0)
                current->indices = top_table[reference].key;
            else
            {
                uint32_t difference=0;
                for(uint32_t j=0; j<4; ++j)
                    difference = difference | (dec_get(codec, &table_difference) << (j*8));
                current->indices =  difference ^ top_table[reference].key;
            }
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

    range_model color_reference, first_index, indices, use_dict, dict_reference;
    model_init(&color_reference, 2);
    model_init(&first_index, 1<<3);
    model_init(&indices, 1<<4);
    model_init(&use_dict, 2);
    model_init(&dict_reference, DICTIONARY_SIZE);

    bc4_block empty_block = {.color = {0, 128}};
    const bc4_block* previous = &empty_block;

    // dictionary initialization
    uint64_t dictionary[DICTIONARY_SIZE];
    uint32_t dict_index = 0;
    for(uint32_t i=0; i<DICTIONARY_SIZE; ++i)
        dictionary[i] = UINT64_MAX;

    for(uint32_t y = 0; y < height_blocks; ++y)
    {
        for(uint32_t x = 0; x < width_blocks; ++x)
        {
            // zig-zag pattern
            uint32_t zigzag_x = (y&1) ? x : width_blocks - x - 1;
            const bc4_block* current = get_block(input, stride, width_blocks, zigzag_x, y);
            for(uint32_t j=0; j<2; ++j)
            {
                int previous_delta = int_abs(current->color[j] - previous->color[j]);
                uint8_t candidate = previous->color[j];
                if (y>0)
                {
                    const bc4_block* up = get_block(input, stride, width_blocks, zigzag_x, y-1);
                    int up_delta = int_abs(current->color[j] - up->color[j]);
                    bool up_better = (up_delta < previous_delta);
                    enc_put(codec, &color_reference,  up_better ? 1 : 0);
                    candidate = (up_better) ? up->color[j] : candidate;
                }
                enc_put(codec, &color_delta[j], delta_encode_wrap(candidate, current->color[j]));
            }

            // search in the dictionary for the current bitfield
            uint64_t bitfield = MAKE48(current->indices[0], current->indices[1], current->indices[2]);
            uint32_t found_index = UINT32_MAX;
            for(uint32_t j=0; j<DICTIONARY_SIZE && found_index == UINT32_MAX; ++j)
                if (dictionary[j] == bitfield)
                    found_index = j;
            
            // found? just write the dictionary index
            if (found_index != UINT32_MAX)
            {
                enc_put(codec, &use_dict, 1);
                enc_put(codec, &dict_reference, found_index);
            }
            else
            {
                // store the entry in the dictionary
                dictionary[dict_index] = bitfield;
                dict_index = (dict_index + 1) & (DICTIONARY_SIZE-1);   // num entries has to be a power of two

                // write the indices with local difference delta encoded
                enc_put(codec, &use_dict, 0);

                uint8_t prev_data = bc4_get_index(current, 0);
                enc_put(codec, &first_index, prev_data);

                // morton delta compress the indices
                for(uint32_t j=1; j<16; ++j)
                {
                    uint8_t data = bc4_get_index(current, morton_inverse[j]);
                    enc_put(codec, &indices, delta_encode(prev_data, data, 4));
                    prev_data = data;
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

    range_model color_reference, first_index, indices, use_dict, dict_reference;
    model_init(&color_reference, 2);
    model_init(&first_index, 1<<3);
    model_init(&indices, 1<<4);
    model_init(&use_dict, 2);
    model_init(&dict_reference, DICTIONARY_SIZE);

    bc4_block empty_block = {.color = {0, 128}};
    bc4_block* previous = &empty_block;

    // dictionary initialization
    uint64_t dictionary[DICTIONARY_SIZE];
    uint32_t dict_index = 0;
    for(uint32_t i=0; i<DICTIONARY_SIZE; ++i)
        dictionary[i] = UINT64_MAX;

    for(uint32_t y = 0; y < height_blocks; ++y)
    {
        for(uint32_t x = 0; x < width_blocks; ++x)
        {
            // zig-zag pattern
            uint32_t zigzag_x = (y&1) ? x : width_blocks - x - 1;
            bc4_block* current = (bc4_block*) get_block(output, stride, width_blocks, zigzag_x, y);
            for(uint32_t j=0; j<2; ++j)
            {
                uint8_t reference = previous->color[j];
                if (y>0 && dec_get(codec, &color_reference))
                {
                    bc4_block* up = (bc4_block*) get_block(output, stride, width_blocks, zigzag_x, y-1);
                    reference = up->color[j];
                }

                uint8_t delta = dec_get(codec, &color_delta[j]);
                current->color[j] = delta_decode_wrap(reference, delta);
            }

            if (dec_get(codec, &use_dict))
            {
                // data should be in the dictionary
                uint64_t bitfield = dictionary[dec_get(codec, &dict_reference)];
                current->indices[0] = (bitfield>>32) & 0xffff;
                current->indices[1] = (bitfield>>16) & 0xffff;
                current->indices[2] = bitfield & 0xffff;
            }
            else
            {
                uint8_t prev_data = dec_get(codec, &first_index);
                bc4_set_index(current, 0, prev_data);

                // morton delta compress the indices
                for(uint32_t j=1; j<16; ++j)
                {
                    uint8_t delta = dec_get(codec, &indices);
                    uint8_t data = delta_decode(prev_data, delta, 4);
                    bc4_set_index(current, morton_inverse[j], data);
                    prev_data = data;
                }

                // store the entry in the dictionary
                dictionary[dict_index] = MAKE48(current->indices[0], current->indices[1], current->indices[2]);
                dict_index = (dict_index + 1) & (DICTIONARY_SIZE-1);
            }
            previous = current;
        }
    }
}

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
            return enc_done(&codec);
        }
    case bc3 : 
        {
            size_t block_size = sizeof(bc1_block) + sizeof(bc4_block);
            bc4_crunch(&codec, cruncher_memory, input, block_size, width, height);
            bc1_crunch(&codec, cruncher_memory, ptr_shift_const(input, sizeof(bc4_block)), block_size, width, height);
            return enc_done(&codec);
        }
    case bc4 : 
        {
            bc4_crunch(&codec, cruncher_memory, input, sizeof(bc4_block), width, height);
            return enc_done(&codec);
        }
    default: return 0;
    }
    return 0;
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
    default: break;
    }
}