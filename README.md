# bc_crunch

![C Standard](https://img.shields.io/badge/C-C99-F4A261)
[![Build Status](https://github.com/Geolm/bc_crunch/actions/workflows/c.yml/badge.svg)](https://github.com/geolm/crunch/actions)


`bc_crunch` is a tiny, dependency-free C99 library for *lossless* compression of GPU-compressed texture blocks **BC1, BC4, BC3, and BC5**.

- ~700 LOC total (single-file encoder/decoder, no build system tricks)
- No malloc, no external libs, no threads, no SIMD requirement  
- Deterministic, bit-exact reconstruction
- Fully unit-tested with cross-validation and fuzz passes  
- Focused on production textures: albedo, masks, normals, heightmaps, etc.
- GPU-ready output: decompressed blocks are written in standard BC1/BC4/BC3/BC5 format, ready to be uploaded directly to GPU memory (or directly written in shared memory)
- No extra memory for decompression: only the encoder needs a temporary buffer; decoding writes straight to the output buffer.



This is *not* another general-purpose compressor.  
`bc_crunch` is specialized for already-compressed GPU formats — it exploits the internal structure of BC1/BC4 blocks, spatial patterns, endpoint deltas, popcount distances, and Morton-ordered indices to achieve significant size reductions with extremely low CPU cost.


---

## Benchmark

Soon to be added.

## Technical details

### BC1
- Zigzag block traversal  
- Endpoint deltas using left/up predictors  
- Top-table of the 256 most frequent index bitfields (prepass histogram)  
- Nearest-match selection via popcount distance  
- Delta-encoded index patches (no raw fallback)

### BC4
- Zigzag traversal  
- Cyclic wrapped endpoint deltas (mod 256)
- 256-entry sliding dictionary
- Morton-order delta path for index fallback  

### Composite Formats
- BC3 is just (BC1 + BC4) — no extra logic  
- BC5 is a dual BC4 channels

