# bc_crunch

![C Standard](https://img.shields.io/badge/C-C99-F4A261)
[![Build Status](https://github.com/Geolm/bc_crunch/actions/workflows/c.yml/badge.svg)](https://github.com/geolm/bc_crunch/actions)


`bc_crunch` is a tiny, dependency-free C99 library for *lossless* compression of GPU-compressed texture blocks **BC1, BC4, BC3, and BC5**.

- ~800 LOC total (single-file encoder/decoder, no build system tricks)
- No malloc, no external libs
- Deterministic, bit-exact reconstruction, fully tested with bytes-for-bytes comparison
- Focused on production textures: albedo, masks, normals, heightmaps, etc.
- GPU-ready output: decompressed blocks are written in standard BC1/BC4/BC3/BC5 format, ready to be uploaded directly to GPU memory (or directly written in shared memory)
- No extra memory for decompression: only the encoder needs a temporary buffer; decoding writes straight to the output buffer.
- Achieves a higher average compression ratio than zlib’s maximum-compression setting

This is *not* another general-purpose compressor. `bc_crunch` is specialized for already-compressed GPU formats — it exploits the internal structure of BC1/BC4 blocks, spatial patterns, endpoint deltas, bitfield indices to achieve significant size reductions with very low CPU cost.

Check out the technical [documentation](./doc.md) for more.

---

## Benchmarks

Compression ratio naturally depends on the input content. Repetitive patterns, smooth gradients, and large uniform regions compress very well, while randomness or high-frequency noise significantly reduces efficiency (e.g., dirt textures are typically very noisy).

### BC1 benchmarks

bc_crunch average compression ratio: **1.593820**  
zlib (best compression) ratio: 1.436144

| Category | Samples | Average bc_crunch Ratio | Average zlib Ratio |
| :--- | :--- | :--- | :--- |
| **kodim** | 24 | 1.5434 | 1.2958 |
| **Wood** | 20 | 1.7068 | 1.3712 |
| **Metal** | 20 | 1.7012 | 1.5511 |
| **Brick** | 20 | 1.6845 | 1.4616 |
| **Elements** | 20 | 1.5540 | 1.3340 |
| **Stone** | 20 | 1.4469 | 1.3061 |
| **Tile** | 20 | 1.4422 | 1.2766 |
| **Plaster** | 20 | 1.3934 | 1.3231 |
| **Dirt** | 20 | 1.3535 | 1.2797 |
| **pk02_floor** | 13 | 2.1121 | 1.9542 |
| **SW_Tree** | 7 | 2.1064 | 2.4498 |
| **blaztree** | 6 | 1.8856 | 2.0893 |
| **pk02_trim** | 5 | 1.8229 | 1.8290 |
| **pkf_concrete** | 4 | 1.4197 | 1.3413 |
| **rock** | 4 | 1.4003 | 1.1317 |

[BC1 textures samples](./textures/bc1/)

### BC4 benchmarks

bc_crunch average compression ratio: **1.336385**  
zlib (best compression) ratio: 1.104125

| Category | Samples | Average bc_crunch Ratio | Average zlib Ratio |
| :--- | :--- | :--- | :--- |
| **grey\_roof\_tiles** | 2 | 1.4623 | 1.1907 |
| **patterned\_cobblestone** | 2 | 1.2525 | 1.0674 |
| **red\_brick** | 2 | 1.3248 | 1.0842 |
| **rough\_wood** | 2 | 1.3059 | 1.0742 |


[BC4 textures samples](./textures/bc4/)


### BC5 benchmarks

Average compression ratio : 1.198832

| Texture | Width | Height | BC5 Size (bytes) | Crunched Size (bytes) | Compression Ratio |
| :--- | ---: | ---: | ---: | ---: | ---: |
| grey\_roof\_tiles\_02\_nor\_dx\_1k.png | 1024 | 1024 | 1,048,576 | 831,233 | 1.261 |
| patterned\_cobblestone\_nor\_dx\_1k.png | 1024 | 1024 | 1,048,576 | 915,387 | 1.146 |
| red\_brick\_nor\_dx\_1k.png | 1024 | 1024 | 1,048,576 | 877,993 | 1.194 |
| rough\_wood\_nor\_dx\_1k.png | 1024 | 1024 | 1,048,576 | 914,356 | 1.147 |
| **Average** | - | - | - | - | **1.187** |


[BC5 textures samples](./textures/bc5/)

### Performance

Using a precomputed decoder table, decrunching is now significantly faster—up to ~1.5× speedup. Some parts of the cruncher have been optimized using NEON instructions (I don't have a x64 available to do the AVX/SSE port).  

Current performance on an M1 Pro MacBook Pro:  

Crunch  100× 1024×1024 texture: 2.26 s  → **22.11 MB/s**  
Decrunch 100× 1024×1024 texture: 1.38 s  → **36.02 MB/s**

[benchmark.c](./test/benchmark.c)

---

## Technical details

### BC1
- Zigzag block traversal  
- Endpoint deltas using left/up predictors  
- Top-table of the 256 most frequent index bitfields (prepass histogram)  
- Nearest-match selection via popcount distance  
- Delta-encoded index patches (no raw fallback)
- Inter-Channel differential coding

### BC4
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

#### How can I compress a texture with all mipmaps?

Right now you need to call bc_crunch separately for each mip level. I may add a helper later, but I want to keep the encoder simple and stateless, without any mip-level awareness or cross-mip data.


#### Is there a plan to improve performance or compression ratio?

Short answer: no. This library was intended as a lightweight testbed for a few ideas, not a state-of-the-art compressor. That said, if you want to experiment with improvements, here are some directions:

* Replace the adaptive range encoder with a static one. This requires gathering statistics in a first pass or defining a generic static model.
* Switch to a cheaper entropy coder (FSE, rANS, Huffman).
* Decompress multiple stream in parallel (multiple texture ou mipmaps).
* For better compression ratio, collecting more data in a first pass could allow a stronger model. Also, BC4 currently has no histogram / first-pass analysis, so there is clear room for improvement there.
