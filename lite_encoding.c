#include "lite_encoding.h"


#define CONTROL_HOT_TIER0 (0)
#define CONTROL_HOT_TIER1 (1)
#define CONTROL_HOT_TIER2 (2)
#define CONTROL_ESCAPE (3)

// ----------------------------------------------------------------------------------------------------------------------------
void le_model_init(le_model *model, const uint32_t *histogram, uint32_t num_symbols)
{
    memset(model->hot_values, 0, sizeof(model->hot_values));
    model->no_compression = false;

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

    // check if the compression will bring some gain
    // if the LE_MODEL_MAX_HOT hot values are not enough present, escape code will make the compressed data bigger
    // a hot value hit uses 4 bits vs 12 bits for escape, basically we need more than 50% of the bytes to use
    // the hot values.

    uint64_t total_count = 0;
    for (uint32_t s = 0; s < num_symbols; s++)
        total_count += histogram[s];

    uint64_t hot_count = 0;
    for (uint32_t i = 0; i < hot_used; i++)
        hot_count += histogram[model->hot_values[i]];

    if (total_count > 0 && hot_count * 2 < total_count)
        model->no_compression = true;

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
void le_encode_byte(le_stream *s, le_model *model, uint8_t value)
{
#ifdef LE_CHECKS
    assert(s->mode == le_mode_encode);
#endif

    // no compression
    if (model->no_compression)
    {
        le_write_byte(s, value);
        return;
    }

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

    uint8_t residual = value - model->cold_min;

    switch(model->cold_num_bits)
    {
        case 2 : le_write_dibit(s, residual); break;
        case 4 : le_write_nibble(s, residual); break;
        case 6 : 
            {
                le_write_dibit(s, residual >> 4);
                le_write_nibble(s, residual&0xf);
                break;
            }
        case 8 : le_write_byte(s, residual);break;
    }

#ifdef LE_STATS
    model->num_raw++;
#endif
}

// ----------------------------------------------------------------------------------------------------------------------------
uint8_t le_decode_byte(le_stream *s, le_model *model)
{
    if (model->no_compression)
        return le_read_byte(s);

    uint8_t dibit = le_read_dibit(s);

    uint8_t value;
    if (dibit == CONTROL_HOT_TIER0)
        return model->hot_values[le_read_dibit(s)];

    if (dibit == CONTROL_HOT_TIER1)
        return model->hot_values[le_read_nibble(s) + 4];
    
    if (dibit == CONTROL_HOT_TIER2)
        return model->hot_values[le_read_nibble(s) + 20];

    // raw
    switch(model->cold_num_bits)
    {
        case 2 : value = le_read_dibit(s); break;
        case 4 : value = le_read_nibble(s); break;
        case 6 : value = (le_read_dibit(s) << 4) | le_read_nibble(s); break;
        case 8 : value = le_read_byte(s);break;
    }
    value += model->cold_min;

    return value;
}

// ----------------------------------------------------------------------------------------------------------------------------
void le_model_save(le_stream *s, const le_model *model)
{
    le_write_dibit(s, model->no_compression ? 1 : 0);
    if (!model->no_compression)
    {
        le_write_nibble(s, model->cold_num_bits);
        le_write_byte(s, model->cold_min);

        for(uint32_t i=0; i<LE_MODEL_MAX_HOT; ++i)
            le_write_byte(s, model->hot_values[i]);
    }
}

// ----------------------------------------------------------------------------------------------------------------------------
void le_model_load(le_stream *s, le_model *model)
{
    model->no_compression = (le_read_dibit(s) == 1);
    if (!model->no_compression)
    {
        model->cold_num_bits = le_read_nibble(s);
        model->cold_min = le_read_byte(s);
        

        for(uint32_t i=0; i<LE_MODEL_MAX_HOT; ++i)
            model->hot_values[i] = le_read_byte(s);
    }
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
void le_encode_delta(le_stream *s, int8_t delta)
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
int8_t le_decode_delta(le_stream *s)
{
    uint8_t base = 0;

    for(uint32_t i = 0; i < 4; i++)
    {
        uint8_t d = le_read_dibit(s);

        if(d < 3)
            return zigzag8_decode(d + base);

        // escape
        base += 3;
    }

    // fallback
    return zigzag8_decode(le_read_byte(s));
}


