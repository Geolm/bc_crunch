/*

zlib License

(C) 2026 Geolm

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


*/


#ifndef LITE_ENCODING_H
#define LITE_ENCODING_H

#include <stdint.h>
#include <stddef.h>
#include <assert.h>
#include <string.h>
#include <stdbool.h>

#define LE_STATS
#define LE_CHECKS
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
    uint8_t bit;
    size_t position;
    size_t size;
    enum le_mode mode;
} le_stream;

typedef struct le_model
{
    uint8_t hot_values[LE_MODEL_MAX_HOT];
    uint8_t cold_min;
    uint8_t cold_max;
    uint8_t cold_num_bits;
    bool no_compression;

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

static void le_init(le_stream *s, void* buffer, size_t size);
static void le_begin_encode(le_stream* s);
static size_t le_end_encode(le_stream* s);
static void le_begin_decode(le_stream* s);
static void le_end_decode(le_stream* s);

static void le_write_byte(le_stream *s, uint8_t value);
static void le_write_nibble(le_stream *s, uint8_t nibble);
static void le_write_dibit(le_stream *s, uint8_t dibit);

static uint8_t le_read_byte(le_stream *s);
static uint8_t le_read_nibble(le_stream *s);
static uint8_t le_read_dibit(le_stream *s);

void le_model_init(le_model *model, const uint32_t *histogram, uint32_t num_symbols);

// Nibble-14 Encoding Model
//
// The Nibble-14 model compresses data by focusing on the 14 most frequent byte values in the input:
// The model stores the 14 hot values derived from the input histogram.
// Each byte is encoded as a 4-bit nibble:
//      0–13 → index into the hot values.
//      14 → run-length code (RLE) for repeating the previous byte.
//      15 → escape code, followed by the full byte if it doesn’t match the hot values or RLE.
//
// This approach efficiently compresses data with skewed distributions and repeated bytes while keeping the stream compact.

void le_encode_byte(le_stream *s, le_model *model, uint8_t value);
uint8_t le_decode_byte(le_stream *s, le_model *model);

// Dibit-delta encoding model, specialized in small delta, use a good predictor up-front to maximize compression
//
// Use dibit to encode delta. 
//      (0, -1, 1)      use 2 bits
//      (-2, +2, -3)    use 4 bits
//      (+3, -4, +4)    use 6 bits
//      (-5, +5, -6)    use 8 bits
//      anything greater use 16 bits
void le_encode_delta(le_stream *s, int8_t delta);
int8_t le_decode_delta(le_stream *s);
void le_model_save(le_stream *s, const le_model *model);
void le_model_load(le_stream *s, le_model *model);

static void histogram_init(le_histogram* h, uint32_t num_symbols);


// ----------------------------------------------------------------------------------------------------------------------------
static inline void le_init(le_stream *s, void* buffer, size_t size)
{
#ifdef LE_CHECKS
    assert(s != NULL);
    assert(buffer != NULL);
#endif

    s->buffer = (uint8_t*) buffer;
    s->size = size;
    s->mode = le_mode_idle;
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline void le_begin_encode(le_stream* s)
{
#ifdef LE_CHECKS
    assert(s->mode == le_mode_idle);
#endif

    memset(s->buffer, 0, s->size);

    s->position = 0;
    s->bit = 0;
    s->mode = le_mode_encode;
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline size_t le_end_encode(le_stream *s)
{
#ifdef LE_CHECKS
    assert(s->mode == le_mode_encode);
#endif

    size_t final_size = s->position;

    if (s->bit > 0)
        final_size++;

    s->mode = le_mode_idle;

    return final_size;
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline void le_begin_decode(le_stream* s)
{
#ifdef LE_CHECKS
    assert(s->mode == le_mode_idle);
#endif

    s->position = 0;
    s->bit = 0;
    s->mode = le_mode_decode;
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline void le_end_decode(le_stream* s)
{
#ifdef LE_CHECKS
    assert(s->mode == le_mode_decode);
#endif

    s->mode = le_mode_idle;
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline void le_write_nibble(le_stream *s, uint8_t nibble)
{
#ifdef LE_CHECKS
    assert(s->mode == le_mode_encode);
    assert(nibble <= 0x0F);
    assert(s->position < s->size);
#endif

    switch (s->bit)
    {
        case 0:
        {
            s->buffer[s->position] |= nibble << 4;
            s->bit = 4;
            break;
        }

        case 2:
        {
            s->buffer[s->position] |= nibble << 2;
            s->bit = 6;
            break;
        }

        case 4:
        {
            s->buffer[s->position] |= nibble;
            s->bit = 0;
            s->position++;
            break;
        }

        case 6:
        {
            assert(s->position + 1 < s->size);

            s->buffer[s->position]     |= nibble >> 2;
            s->buffer[s->position + 1] |= nibble << 6;
            s->bit = 2;
            s->position++;
            break;
        }

        default:
            assert(0);
    }
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline uint8_t le_read_nibble(le_stream *s)
{
#ifdef LE_CHECKS
    assert(s->mode == le_mode_decode);
    assert(s->position < s->size);
#endif

    uint8_t nibble;
    switch (s->bit)
    {
        case 0:
        {
            nibble = s->buffer[s->position] >> 4;
            s->bit = 4;
            break;
        }

        case 2:
        {
            nibble = (s->buffer[s->position] >> 2) & 0x0F;
            s->bit = 6;
            break;
        }

        case 4:
        {
            nibble = s->buffer[s->position] & 0x0F;
            s->bit = 0;
            s->position++;
            break;
        }

        case 6:
        {
            assert(s->position + 1 < s->size);

            uint8_t hi = (s->buffer[s->position] & 0x03) << 2;
            uint8_t lo =  s->buffer[s->position + 1] >> 6;
            nibble = hi | lo;

            s->bit = 2;
            s->position++;
            break;
        }

        default:
            assert(0);
    }

    return nibble;
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline void le_write_dibit(le_stream *s, uint8_t dibit)
{
#ifdef LE_CHECKS
    assert(s->mode == le_mode_encode);
    assert(dibit <= 0x03);
    assert(s->position < s->size);
#endif

    uint8_t shift = 6 - s->bit;
    s->buffer[s->position] |= (dibit << shift);
    s->bit += 2;

    if (s->bit == 8)
    {
        s->bit = 0;
        s->position++;
    }
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline void le_write_byte(le_stream *s, uint8_t value)
{
#ifdef LE_CHECKS
    assert(s->mode == le_mode_encode);
#endif

    if (s->bit == 0) 
    {
    #ifdef LE_CHECKS
        assert(s->position < s->size);
    #endif

        s->buffer[s->position++] = value;
    } 
    else 
    {
        // byte is split across two memory locations

    #ifdef LE_CHECKS
        assert(s->position + 1 < s->size);
    #endif
        
        uint8_t shift = s->bit;
        uint8_t remaining = 8 - shift;

        s->buffer[s->position]     |= (value >> shift);
        s->buffer[s->position + 1] |= (value << remaining);
        s->position++;
    }
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline uint8_t le_read_dibit(le_stream *s)
{
#ifdef LE_CHECKS
    assert(s->mode == le_mode_decode);
    assert(s->position < s->size);
#endif

    // Extract 2 bits based on the current bit offset (reading from MSB to LSB)
    uint8_t shift = 6 - s->bit;
    uint8_t dibit = (s->buffer[s->position] >> shift) & 0x03;

    s->bit += 2;

    if (s->bit == 8)
    {
        s->bit = 0;
        s->position++;
    }

    return dibit;
}

// ----------------------------------------------------------------------------------------------------------------------------
static inline uint8_t le_read_byte(le_stream *s)
{
#ifdef LE_CHECKS
    assert(s->mode == le_mode_decode);
#endif

    uint8_t value = 0;
    if (s->bit == 0)
    {
    #ifdef LE_CHECKS
        assert(s->position < s->size);
    #endif
        value = s->buffer[s->position++];
    }
    else
    {
    #ifdef LE_CHECKS
        assert(s->position + 1 < s->size);
    #endif

        uint8_t shift = s->bit;
        uint8_t remaining = 8 - shift;

        uint8_t hi = (s->buffer[s->position] << shift);
        uint8_t lo = (s->buffer[s->position + 1] >> remaining);
        
        value = hi | lo;
        s->position++;
        // s->bit remains unchanged (e.g., if it was 4, it's still 4 in the next byte)
    }

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

#endif

