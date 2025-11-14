#ifndef __BC_CRUNCH_H__
#define __BC_CRUNCH_H__

#include <stdint.h>
#include <stddef.h>


enum bc_format
{
    bc1,     // rgb only
    bc3,     // rgba
    bc4,     // red only
    bc5      // red-green (a.k.a normalmap)
};

#ifdef __cplusplus
extern "C" {
#endif




//----------------------------------------------------------------------------------------------------------------------------
// Returns the size in bytes needed for the cruncher
size_t crunch_min_size(void);

//----------------------------------------------------------------------------------------------------------------------------
// Compresses with lossless technique a BC image
//      [cruncher_memory]   Memory for the cruncher process, allocated by the user,
//                          should at least be crunch_min_size() bytes large and aligned on 8 bytes
//      [input]             Pointer to BC blocks
//      [width, height]     Size of the image in pixels
//      [format]            Format of the image (BC1, BC3, BC4 or BC5)
//      [output]            Pointer to the pre-allocated buffer for compressed bc1, most of the time smaller than
//                          the original bc1 image but in case of random/weird image can be a bit bigger. To be safe
//                          allocate twice the size of the bc1 image
//      [length]            Length of the compressed image buffer
size_t bc_crunch(void* cruncher_memory, const void* input, uint32_t width, uint32_t height, enum bc_format format, void* output, size_t length);


//----------------------------------------------------------------------------------------------------------------------------
// Decompresses the crunched image into a BC1 image ready to be used by the GPU
//      [input]             Pointer to the crunched data
//      [length]            Length in bytes of the crunched data
//      [width, height]     Size of the image in pixels
//      [output]            User-allocated memory to receive the BC1 image (should be big enough)
//      [format]            Format of the output, *must* match the input format used in bc_crunch
void bc_decrunch(const void* input, size_t length, uint32_t width, uint32_t height, enum bc_format format, void* output);


#ifdef __cplusplus
}
#endif

#endif
