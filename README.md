# bc_crunch

![C Standard](https://img.shields.io/badge/C-C99-F4A261)
[![Build Status](https://github.com/Geolm/bc_crunch/actions/workflows/c.yml/badge.svg)](https://github.com/geolm/crunch/actions)


`bc_crunch` is a tiny, dependency-free C99 library for *lossless* compression of GPU-compressed texture blocks **BC1, BC4, BC3, and BC5**.

- ~800 LOC total (single-file encoder/decoder, no build system tricks)
- No malloc, no external libs
- Deterministic, bit-exact reconstruction, fully tested with bytes-for-bytes comparison
- Focused on production textures: albedo, masks, normals, heightmaps, etc.
- GPU-ready output: decompressed blocks are written in standard BC1/BC4/BC3/BC5 format, ready to be uploaded directly to GPU memory (or directly written in shared memory)
- No extra memory for decompression: only the encoder needs a temporary buffer; decoding writes straight to the output buffer.

This is *not* another general-purpose compressor. `bc_crunch` is specialized for already-compressed GPU formats — it exploits the internal structure of BC1/BC4 blocks, spatial patterns, endpoint deltas, bitfield indices to achieve significant size reductions with very low CPU cost.

---

## Benchmarks

Compression ratio naturally depends on the input content. Repetitive patterns, smooth gradients, and large uniform regions compress very well, while randomness or high-frequency noise significantly reduces efficiency (e.g., dirt textures are typically very noisy).

### BC1 benchmarks

Average compression ratio : 1.493586

| Category       | # Samples | Original BC1 Size (bytes) | Avg Crunched Size (bytes) | Avg Compression Ratio |
|----------------|-----------|---------------------------|---------------------------|---------------------|
| Kodak photos   | 24        | 196,608                   | 136,254                   | 1.44               |
| Dirt           | 20        | 131,072                   | 105,386                   | 1.24               |
| Wood           | 20        | 131,072                   | 79,306                    | 1.65               |
| Elements       | 20        | 131,072                   | 90,412                    | 1.45               |
| Brick          | 20        | 131,072                   | 88,212                    | 1.49               |
| Metal          | 20        | 131,072                   | 87,132                    | 1.51               |
| Plaster        | 20        | 131,072                   | 99,118                    | 1.33               |
| Stone          | 20        | 131,072                   | 97,014                    | 1.35               |
| Tile           | 20        | 131,072                   | 101,087                   | 1.30               |
| Blaz Tree      | 7         | 32,768–32,768             | 16,900                    | 1.84               |
| SW Tree        | 7         | 20,480–32,768             | 14,227                    | 2.08               |
| pkf Concrete   | 4         | 32,768                    | 24,583                    | 1.33               |
| pk02 Floor     | 20        | 32,768–524,288            | 100,604                   | 1.83               |
| pk02 Trim      | 5         | 4,096–131,072             | 11,777                    | 1.89               |
| Rock           | 4         | 524,288                   | 412,365                   | 1.27               |
| Brick Large    | 2         | 524,288                   | 408,746                   | 1.28               |
| **Average**    | -         | -                         | -                         | 1.493586           |

[BC1 textures samples](./textures/bc1/)

### BC4 benchmarks

Average compression ratio : 1.232521

| Category     | Samples | Uncompressed (bytes) | Avg Compressed (bytes) | Avg Compression Ratio |
| ------------ | ------- | -------------------- | ---------------------- | --------------------- |
| AO           | 4       | 524,288              | 472,347                | 1.110685              |
| Displacement | 4       | 524,288              | 390,158                | 1.354357              |
| **Average**  | -       | -                    | -                      | 1.232521           |


[BC4 textures samples](./textures/bc4/)


### BC5 benchmarks

Average compression ratio : 1.150844

| Texture                          | Width | Height | BC5 Size (bytes) | Crunched Size (bytes) | Compression Ratio |
|----------------------------------|-------|--------|-----------------|----------------------|-----------------|
| grey_roof_tiles_02_nor_dx_1k.png | 1024  | 1024   | 1,048,576       | 856,423              | 1.224           |
| patterned_cobblestone_nor_dx_1k.png | 1024 | 1024  | 1,048,576       | 945,970              | 1.108           |
| red_brick_nor_dx_1k.png          | 1024  | 1024   | 1,048,576       | 910,903              | 1.151           |
| rough_wood_nor_dx_1k.png         | 1024  | 1024   | 1,048,576       | 936,729              | 1.119           |
| **Average**                       | -     | -      | -               | -                    | 1.151           |


[BC5 textures samples](./textures/bc5/)

### Performance

Using a precomputed decoder table, decrunching is now significantly faster—up to ~1.5× speedup. Further optimizations are in progress.

Current performance on an M1 Pro MacBook Pro:  

Crunch  50× 1024×1024 texture: 1.602 s  → **15.61 MB/s**  
Decrunch 50× 1024×1024 texture: 0.723 s  → **34.56 MB/s**

[benchmark.c](./test/benchmark.c)

---

## Technical details

### BC1
- Zigzag block traversal  
- Endpoint deltas using left/up predictors  
- Top-table of the 256 most frequent index bitfields (prepass histogram)  
- Nearest-match selection via popcount distance  
- Delta-encoded index patches (no raw fallback)

### BC4
- Zigzag block traversal  
- Cyclic wrapped endpoint deltas (mod 256) using left/up predictors
- 256-entry sliding dictionary
- Move-to-front heuristic when hit
- Nearest-match selection via popcount distance
- Block zigzag traversal with xor-delta encoding for index fallback  

### Composite Formats
- BC3 is just (BC1 + BC4) — no extra logic  
- BC5 is a dual independent BC4 channels

## How to build tests

* open a terminal in the folder
* mkdir build
* cd build
* cmake -DCMAKE_BUILD_TYPE=Release ..
* cmake --build .
* ./test
* ./benchmark

More tests will be added soon.

## FAQ

#### How can I compress a texture with all mipmaps

Right now you need to call bc_crunch separately for each mip level. I may add a helper later, but I want to keep the encoder simple and stateless, without any mip-level awareness or cross-mip data.
