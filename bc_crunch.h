#ifndef __BC_CRUNCH_H__
#define __BC_CRUNCH_H__

#include <stdint.h>
#include <stddef.h>

size_t bc1_crunch(const void* data, uint32_t width, uint32_t height, void* output, size_t length);

#endif
