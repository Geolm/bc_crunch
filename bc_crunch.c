#include <assert.h>
#include "bc_crunch.h"

#include <stdio.h> // NOT NEEDED

//----------------------------------------------------------------------------------------------------------------------------
// Private structures & functions
//----------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------
#define TOP                 0x01000000u
#define MAX_ALPHABET_SIZE   256
#define LIMIT               ((MAX_ALPHABET_SIZE * 128 < (1 << 16)) ? (MAX_ALPHABET_SIZE * 128) : (1 << 16))

//----------------------------------------------------------------------------------------------------------------------------
typedef struct
{
    uint16_t count[MAX_ALPHABET_SIZE];
    uint32_t total;
    uint32_t alphabet_size;
} range_model;

//----------------------------------------------------------------------------------------------------------------------------
static inline void model_init(range_model *m, uint32_t alphabet_size)
{
    assert(alphabet_size <= MAX_ALPHABET_SIZE);
    m->alphabet_size = alphabet_size;
    for (uint32_t i = 0; i < alphabet_size; ++i)
        m->count[i] = 1;
    m->total = alphabet_size;
}

//----------------------------------------------------------------------------------------------------------------------------
static inline void model_rescale(range_model *m)
{
    m->total = 0;
    for (uint32_t i = 0; i < m->alphabet_size; ++i)
    {
        m->count[i] = (m->count[i] + 1) >> 1;
        m->total += m->count[i];
    }
}

//----------------------------------------------------------------------------------------------------------------------------
static inline void model_update(range_model *m, uint32_t sym)
{
    assert(sym < m->alphabet_size);
    m->count[sym]++;
    m->total++;
    if (m->total >= LIMIT)
        model_rescale(m);
}

//----------------------------------------------------------------------------------------------------------------------------
static inline void model_build(const range_model *m, uint32_t *freq)
{
    freq[0] = 0;
    for (uint32_t i = 0; i < m->alphabet_size; ++i)
        freq[i + 1] = freq[i] + m->count[i];
}

//----------------------------------------------------------------------------------------------------------------------------
typedef struct
{
    uint32_t low, high, code;
    uint8_t *start, *ptr, *end;
} range_codec;

//----------------------------------------------------------------------------------------------------------------------------
static inline void enc_init(range_codec *rc, uint8_t *out, size_t cap)
{
    rc->low  = 0;
    rc->high = 0xFFFFFFFFu;
    rc->start = out;
    rc->ptr  = out;
    rc->end  = out + cap;
}

//----------------------------------------------------------------------------------------------------------------------------
void enc_put(range_codec *rc, range_model *m, uint32_t sym)
{
    assert(sym < m->alphabet_size);

    uint32_t freq[MAX_ALPHABET_SIZE + 1];
    model_build(m, freq);
    uint32_t total = m->total;

    uint64_t range = (uint64_t)rc->high - rc->low + 1ull;
    rc->high = rc->low + (uint32_t)((range * (uint64_t)freq[sym + 1]) / (uint64_t)total) - 1u;
    rc->low  = rc->low + (uint32_t)((range * (uint64_t)freq[sym])     / (uint64_t)total);

    while ((rc->high ^ rc->low) < TOP) 
    {
        if (rc->ptr >= rc->end)
            return;

        *rc->ptr++ = (uint8_t)(rc->low >> 24);
        rc->low <<= 8;
        rc->high = (rc->high << 8) | 0xFFu;
    }

    model_update(m, sym);
}


//----------------------------------------------------------------------------------------------------------------------------
static inline size_t enc_done(range_codec *rc)
{
    for (uint32_t i = 0; i < 4; ++i)
    {
        if (rc->ptr >= rc->end)
            break;

        *rc->ptr++ = rc->low >> 24;
        rc->low <<= 8;
    }
    return (size_t)(rc->ptr - rc->start);
}

//----------------------------------------------------------------------------------------------------------------------------
static inline void dec_init(range_codec *rc, const uint8_t *in, size_t size)
{
    rc->low = 0;
    rc->high = 0xFFFFFFFFu;
    rc->ptr = (uint8_t*)in;
    rc->end = (uint8_t*)in + size;
    rc->code = 0;
    for (uint32_t i = 0; i < 4; ++i)
        rc->code = (rc->code << 8) | (rc->ptr < rc->end ? *rc->ptr++ : 0);
}

//----------------------------------------------------------------------------------------------------------------------------
uint32_t dec_get(range_codec *rc, range_model *m)
{
    uint32_t freq[MAX_ALPHABET_SIZE + 1];
    model_build(m, freq);
    uint32_t total = m->total;

    uint64_t range = (uint64_t)rc->high - rc->low + 1ull;
    uint64_t tmp = (uint64_t)(rc->code - rc->low + 1u);
    uint32_t scaled = (uint32_t)((tmp * (uint64_t)total - 1ull) / range);

    uint32_t sym = 0;
    while (freq[sym + 1] <= scaled)
        ++sym;

    rc->high = rc->low + (uint32_t)((range * (uint64_t)freq[sym + 1]) / (uint64_t)total) - 1u;
    rc->low  = rc->low + (uint32_t)((range * (uint64_t)freq[sym])     / (uint64_t)total);

    while ((rc->high ^ rc->low) < TOP)
    {
        rc->code = (rc->code << 8) | (rc->ptr < rc->end ? *rc->ptr++ : 0);
        rc->low <<= 8;
        rc->high = (rc->high << 8) | 0xFFu;
    }

    model_update(m, sym);
    return sym;
}


//----------------------------------------------------------------------------------------------------------------------------
typedef struct bc1_block
{
    uint16_t color[2];
    uint8_t indices[4];
} bc1_block;

//----------------------------------------------------------------------------------------------------------------------------
static inline uint8_t delta_encode(uint8_t prev, uint8_t curr, uint8_t num_bits) 
{
    uint8_t mask = (1<<num_bits)-1;
    return (curr - prev) & mask;
}

//----------------------------------------------------------------------------------------------------------------------------
static inline uint8_t delta_decode(uint8_t prev, uint8_t delta, uint8_t num_bits)
{
    uint8_t mask = (1<<num_bits)-1;
    return (prev + delta) & mask;
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



#define RED_NUM_BITS (5)
#define GREEN_NUM_BITS (6)
#define BLUE_NUM_BITS (5)
#define BC1_ROW_INDICES_NUM_BITS (8)

//----------------------------------------------------------------------------------------------------------------------------
// Public functions
//----------------------------------------------------------------------------------------------------------------------------

size_t bc1_crunch(const void* input, uint32_t width, uint32_t height, void* output, size_t length)
{
    assert((width%4 == 0) && (height%4 == 0));

    uint32_t num_blocks = (width * height) / 16;

    range_codec codec;
    enc_init(&codec, (uint8_t*) output, length);

    range_model red[2], green[2], blue[2];
    for(uint32_t i=0; i<2; ++i)
    {
        model_init(&red[i],   1 << RED_NUM_BITS);
        model_init(&green[i], 1 << GREEN_NUM_BITS);
        model_init(&blue[i],  1 << BLUE_NUM_BITS);
    }

    range_model rows[4];
    for(uint32_t i=0; i<4; ++i)
        model_init(&rows[i], 1<<BC1_ROW_INDICES_NUM_BITS);

    bc1_block empty_block = {0};
    bc1_block* blocks = (bc1_block*) input;
    bc1_block* previous = &empty_block;

    for(uint32_t i=0; i<num_blocks; ++i)
    {
        bc1_block* current = &blocks[i];

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

        for(uint32_t j=0; j<4; ++j)
            enc_put(&codec, &rows[j], current->indices[j] ^ previous->indices[j]);

        previous = current;
    }

    return enc_done(&codec);
}

//----------------------------------------------------------------------------------------------------------------------------
void bc1_decrunch(const void* input, size_t length, uint32_t width, uint32_t height, void* output)
{
    assert((width % 4 == 0) && (height % 4 == 0));

    uint32_t num_blocks = (width * height) / 16;

    range_codec codec;
    dec_init(&codec, (const uint8_t*)input, length);

    range_model red[2], green[2], blue[2];
    for (uint32_t i = 0; i < 2; ++i)
    {
        model_init(&red[i],   1 << RED_NUM_BITS);
        model_init(&green[i], 1 << GREEN_NUM_BITS);
        model_init(&blue[i],  1 << BLUE_NUM_BITS);
    }

    range_model rows[4];
    for (uint32_t i = 0; i < 4; ++i)
        model_init(&rows[i], 1 << BC1_ROW_INDICES_NUM_BITS);

    bc1_block empty_block = {0};
    bc1_block* blocks = (bc1_block*)output;
    bc1_block* previous = &empty_block;

    for (uint32_t i = 0; i < num_blocks; ++i)
    {
        bc1_block* current = &blocks[i];

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

        for (uint32_t j = 0; j < 4; ++j)
        {
            uint8_t delta = (uint8_t)dec_get(&codec, &rows[j]);
            current->indices[j] = previous->indices[j] ^ delta;
        }

        previous = current;
    }
}
