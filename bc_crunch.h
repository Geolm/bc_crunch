#ifndef __BC_CRUNCH_H__
#define __BC_CRUNCH_H__

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif


//----------------------------------------------------------------------------------------------------------------------------
// Returns the size in bytes needed for the cruncher
size_t crunch_min_size(void);

//----------------------------------------------------------------------------------------------------------------------------
// Compresses with lossless technique a BC1 image
//      [cruncher_memory]   Memory for the cruncher process, allocated by the user,
//                          should at least be crunch_min_size() bytes large and aligned on 8 bytes
//      [input]             Pointer to BC1 blocks
//      [width, height]     Size of the image in pixels
//      [output]            Pointer to the pre-allocated buffer for compressed bc1, most of the time smaller than
//                          the original bc1 image but in case of random/weird image can be a bit bigger. To be safe
//                          allocate twice the size of the bc1 image
//      [length]            Length of the compressed image buffer
size_t bc1_crunch(void* cruncher_memory, const void* input, uint32_t width, uint32_t height, void* output, size_t length);

//----------------------------------------------------------------------------------------------------------------------------
// Decompresses 
void bc1_decrunch(const void* input, size_t length, uint32_t width, uint32_t height, void* output);

#ifdef __cplusplus
}
#endif

#endif
