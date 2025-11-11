#include <assert.h>
#include "bc_crunch.h"

#include <stdbool.h>
#include <stdio.h> // NOT NEEDED

//----------------------------------------------------------------------------------------------------------------------------
// Private structures & functions
//----------------------------------------------------------------------------------------------------------------------------

#define AC__MinLength (0x01000000U)
#define AC__MaxLength (0xFFFFFFFFU)
#define MAX_ALPHABET_SIZE (256)
#define DM__LengthShift (15)
#define DM__MaxCount    (1 << DM__LengthShift)
#define HASHMAP_SIZE (1 << 20)
#define TABLE_INDEX_NUM_BITS (8)
#define TABLE_SIZE (1<<TABLE_INDEX_NUM_BITS)
#define RED_NUM_BITS (6)
#define GREEN_NUM_BITS (7)
#define BLUE_NUM_BITS (6)
#define BC1_ROW_INDICES_NUM_BITS (8)

typedef struct range_model
{
    uint32_t distribution[MAX_ALPHABET_SIZE]; 
    uint32_t symbol_count[MAX_ALPHABET_SIZE];
    uint32_t total_count, update_cycle, symbols_until_update;
    uint32_t data_symbols, last_symbol;
} range_model;

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
    uint32_t base, value, length;                     // arithmetic coding state
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
    } while ((codec->length <<= 8) < AC__MinLength);        // length multiplied by 256
}

//----------------------------------------------------------------------------------------------------------------------
void enc_init(range_codec* codec, uint8_t* output, uint32_t length)
{
    codec->buffer_size = length;
    codec->code_buffer = output;
    codec->base = 0;
    codec->length = AC__MaxLength;
    codec->ac_pointer = codec->code_buffer;
}

//----------------------------------------------------------------------------------------------------------------------
void dec_init(range_codec* codec, const uint8_t* input, uint32_t length)
{
    codec->buffer_size = length;
    codec->code_buffer = (uint8_t*) input;
    codec->length = AC__MaxLength;
    codec->ac_pointer = codec->code_buffer + 3;
    codec->value = ((uint32_t)(codec->code_buffer[0]) << 24) |
                   ((uint32_t)(codec->code_buffer[1]) << 16) |
                   ((uint32_t)(codec->code_buffer[2]) <<  8) |
                    (uint32_t)(codec->code_buffer[3]);
}

//----------------------------------------------------------------------------------------------------------------------
uint32_t enc_done(range_codec* codec)
{
    uint32_t init_base = codec->base;            // done encoding: set final data bytes

    if (codec->length > 2 * AC__MinLength) 
    {
        codec->base  += AC__MinLength;                                     // base offset
        codec->length = AC__MinLength >> 1;             // set new length for 1 more byte
    }
    else 
    {
        codec->base  += AC__MinLength >> 1;                                // base offset
        codec->length = AC__MinLength >> 9;            // set new length for 2 more bytes
    }

    if (init_base > codec->base) 
        propagate_carry(codec);                 // overflow = carry

    renorm_enc_interval(codec);                // renormalization = output last bytes

    uint32_t code_bytes = (uint32_t)(codec->ac_pointer - codec->code_buffer);
    assert(code_bytes <= codec->buffer_size); // code buffer overflow

    return code_bytes;                                   // number of bytes used
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

    if (codec->length < AC__MinLength) 
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

    if (codec->length < AC__MinLength) 
        renorm_enc_interval(codec);
}

//----------------------------------------------------------------------------------------------------------------------
uint32_t dec_get(range_codec* codec, range_model* model)
{
    assert(model->distribution != NULL); // adaptive model should be initialized

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

    if (codec->length < AC__MinLength)
    {
        do
        {
            codec->value = (codec->value << 8) | (uint32_t)(*++codec->ac_pointer);
        } while ((codec->length <<= 8) < AC__MinLength);
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
    if (codec->length < AC__MinLength)
    {
        do
        {
            codec->value = (codec->value << 8) | (uint32_t)(*++codec->ac_pointer);
        } while ((codec->length <<= 8) < AC__MinLength);
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
static inline void bc1_extract_565(uint16_t color, uint8_t *r5, uint8_t *g6, uint8_t *b5)
{
    *r5 = (uint8_t)((color >> 11) & 0x1F);   // 5 bits
    *g6 = (uint8_t)((color >> 5)  & 0x3F);   // 6 bits
    *b5 = (uint8_t)(color & 0x1F);           // 5 bits
}

//----------------------------------------------------------------------------------------------------------------------------
static inline uint16_t bc1_pack_565(uint8_t r5, uint8_t g6, uint8_t b5)
{
    return (uint16_t)(((uint16_t)r5 << 11) | ((uint16_t)g6 << 5) | (uint16_t)b5);
}

//----------------------------------------------------------------------------------------------------------------------------
typedef struct entry
{
    uint32_t key;
    uint32_t count;
} entry;

static entry hashmap[HASHMAP_SIZE];

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
void build_top_table(const bc1_block* blocks, uint32_t num_blocks, entry* output, uint32_t* num_entries)
{
    // clear the hashmap
    for(uint32_t i=0; i<HASHMAP_SIZE; ++i)
        hashmap[i].count = 0;

    // insert all blocks indices in the hashmap
    for(uint32_t i=0; i<num_blocks; ++i)
    {
        const bc1_block* b = &blocks[i];

        uint32_t h = hash32(b->indices);
        uint32_t index = h & (HASHMAP_SIZE - 1);

        bool inserted = false;
        while (!inserted)
        {
            if (hashmap[index].count == 0)
            {
                hashmap[index].key = b->indices;
                hashmap[index].count = 1;
                inserted = true;
            }
            else if (hashmap[index].key == b->indices)
            {
                hashmap[index].count++;
                inserted = true;
            }
            else
                index = (index + 1) & (HASHMAP_SIZE - 1); // linear probe
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
        uint32_t xor_value = table[i].key ^ original_indices;
        uint32_t score = __builtin_popcount(xor_value);
        if (score < best_score)
        {
            best_score = score;
            best_index = i;
        }
    }
    return best_index;
}

//----------------------------------------------------------------------------------------------------------------------------
// Public functions
//----------------------------------------------------------------------------------------------------------------------------

size_t bc1_crunch(const void* input, uint32_t width, uint32_t height, void* output, size_t length)
{
    assert((width%4 == 0) && (height%4 == 0));

    bc1_block* blocks = (bc1_block*) input;
    uint32_t height_blocks = height/4;
    uint32_t width_blocks = width/4;

    range_codec codec;
    enc_init(&codec, (uint8_t*) output, length);

    // build a histogram and select the TABLE_SIZE block indices which are most used
    entry top_table[TABLE_SIZE];
    uint32_t top_table_size;
    build_top_table(blocks, height_blocks*width_blocks, top_table, &top_table_size);

    // write the table (no compression for now)
    enc_put_bits(&codec, top_table_size-1, TABLE_INDEX_NUM_BITS);
    for(uint32_t i=0; i<top_table_size; ++i)
        for(uint32_t j=0; j<4; ++j)
            enc_put_bits(&codec, (top_table[i].key >> (j*8)) & 0xff, 8);

    range_model red[2], green[2], blue[2];
    for(uint32_t i=0; i<2; ++i)
    {
        model_init(&red[i],   1 << RED_NUM_BITS);
        model_init(&green[i], 1 << GREEN_NUM_BITS);
        model_init(&blue[i],  1 << BLUE_NUM_BITS);
    }

    range_model table_index, table_difference, block_mode;
    model_init(&table_index, top_table_size);
    model_init(&table_difference, 256);
    model_init(&block_mode, 4);

    bc1_block empty_block = {0};
    bc1_block* previous = &empty_block;

    for(uint32_t y = 0; y < height_blocks; ++y)
    {
        for(uint32_t x = 0; x < width_blocks; ++x)
        {
            // zig-zag pattern delta compression for colors
            bc1_block* current = (y&1) ? &blocks[y * width_blocks + x] : &blocks[y * width_blocks + width_blocks-x-1];
            for(uint32_t j=0; j<2; ++j)
            {
                uint8_t current_red, current_green, current_blue;
                uint8_t previous_red, previous_green, previous_blue;

                bc1_extract_565(current->color[j], &current_red, &current_green, &current_blue);
                bc1_extract_565(previous->color[j], &previous_red, &previous_green, &previous_blue);

                enc_put(&codec, &red[j], delta_encode(previous_red, current_red, RED_NUM_BITS));
                enc_put(&codec, &green[j], delta_encode(previous_green, current_green, GREEN_NUM_BITS));
                enc_put(&codec, &blue[j], delta_encode(previous_blue, current_blue, BLUE_NUM_BITS));
            }

            // for indices, we store the reference to "nearest" indices (can be exactly the same)
            // and the delta with this reference
            uint32_t reference = top_table_nearest(top_table, top_table_size, current->indices);
            enc_put(&codec, &table_index, reference);

            // xor the difference and encode (could be 0 if equal to reference)
            uint32_t difference = current->indices ^ top_table[reference].key;
            uint32_t hamming = __builtin_popcount(difference);
            
            if (hamming == 0)
            {
                enc_put(&codec, &block_mode, 0);
            }
            else if (hamming == 1)
            {
                enc_put(&codec, &block_mode, 1);
                enc_put_bits(&codec, 31 - __builtin_clz(difference), 5);
            }
            else if (hamming == 2)
            {
                enc_put(&codec, &block_mode, 2);

                uint32_t bit_index = 31 - __builtin_clz(difference);
                enc_put_bits(&codec, bit_index, 5);
                difference &= ~(1<<bit_index);
                bit_index = 31 - __builtin_clz(difference);
                enc_put_bits(&codec, bit_index, 5);
            }
            else
            {
                enc_put(&codec, &block_mode, 3);

                // encode the full 32 bits difference
                for(uint32_t j=0; j<4; ++j)
                    enc_put(&codec, &table_difference, (difference>>(j*8))&0xff);
            }

            previous = current;
        }
    }

    return enc_done(&codec);
}

//----------------------------------------------------------------------------------------------------------------------------
void bc1_decrunch(const void* input, size_t length, uint32_t width, uint32_t height, void* output)
{
    assert((width % 4 == 0) && (height % 4 == 0));

    uint32_t height_blocks = height/4;
    uint32_t width_blocks = width/4;

    range_codec codec;
    dec_init(&codec, (const uint8_t*)input, length);

    range_model red[2], green[2], blue[2];
    for (uint32_t i = 0; i < 2; ++i)
    {
        model_init(&red[i],   1 << RED_NUM_BITS);
        model_init(&green[i], 1 << GREEN_NUM_BITS);
        model_init(&blue[i],  1 << BLUE_NUM_BITS);
    }

    entry top_table[TABLE_SIZE];
    uint32_t top_table_size = dec_get_bits(&codec, TABLE_INDEX_NUM_BITS)+1;

    for(uint32_t i=0; i<top_table_size; ++i)
    {
        top_table[i].key = 0;
        for(uint32_t j=0; j<4; ++j)
            top_table[i].key |= dec_get_bits(&codec, 8) << (j*8);
    }

    range_model table_index, table_difference, block_mode;
    model_init(&table_index, top_table_size);
    model_init(&table_difference, 256);
    model_init(&block_mode, 4);

    bc1_block empty_block = {0};
    bc1_block* blocks = (bc1_block*)output;
    bc1_block* previous = &empty_block;

    for(uint32_t y = 0; y < height_blocks; ++y)
    {
        for(uint32_t x = 0; x < width_blocks; ++x)
        {
            // zig-zag pattern color
            bc1_block* current = (y&1) ? &blocks[y * width_blocks + x] : &blocks[y * width_blocks + width_blocks-x-1];
            for (uint32_t j = 0; j < 2; ++j)
            {
                uint8_t previous_red, previous_green, previous_blue;
                bc1_extract_565(previous->color[j], &previous_red, &previous_green, &previous_blue);

                uint8_t delta_red = (uint8_t)dec_get(&codec, &red[j]);
                uint8_t delta_green = (uint8_t)dec_get(&codec, &green[j]);
                uint8_t delta_blue = (uint8_t)dec_get(&codec, &blue[j]);

                uint8_t current_red = delta_decode(previous_red, delta_red, RED_NUM_BITS);
                uint8_t current_green = delta_decode(previous_green, delta_green, GREEN_NUM_BITS);
                uint8_t current_blue = delta_decode(previous_blue, delta_blue, BLUE_NUM_BITS);

                current->color[j] = bc1_pack_565(current_red, current_green, current_blue);
            }

            // indices difference with top table
            uint32_t reference = dec_get(&codec, &table_index);
            uint32_t mode = dec_get(&codec, &block_mode);

            switch(mode)
            {
            case 0 : current->indices = top_table[reference].key; break;
            case 1 : 
                {
                    current->indices = top_table[reference].key;
                    current->indices ^= (1 << dec_get_bits(&codec, 5));
                    break;
                }
            case 2 : 
                {
                    current->indices = top_table[reference].key;
                    current->indices ^= (1 << dec_get_bits(&codec, 5));
                    current->indices ^= (1 << dec_get_bits(&codec, 5));
                    break;
                }
            case 3 :
                {
                    uint32_t difference=0;
                    for(uint32_t j=0; j<4; ++j)
                        difference = difference | (dec_get(&codec, &table_difference) << (j*8));

                    current->indices =  difference ^ top_table[reference].key;
                }
            }
            previous = current;
        }
    }
}
