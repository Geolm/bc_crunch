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

This is *not* another general-purpose compressor. `bc_crunch` is specialized for already-compressed GPU formats — it exploits the internal structure of BC1/BC4 blocks, spatial patterns, endpoint deltas, bitfield indices to achieve significant size reductions with very low CPU cost.

---

## Benchmarks

Compression ratio naturally depends on the input content. Repetitive patterns, smooth gradients, and large uniform regions compress very well, while randomness or high-frequency noise significantly reduces efficiency (e.g., dirt textures are typically very noisy).

### BC1 benchmarks

Average compression ratio : 1.580368

| Texture Type | Number of Textures | Total BC1 Size (bytes) | Total Crunched Size (bytes) | Average Compression Ratio |
|:---------------------|---------------------:|-------------------------:|------------------------------:|----------------------------:|
| Brick Textures | 20 | 2621440 | 1644322 | 1.594 |
| Dirt Textures | 20 | 2621440 | 1957069 | 1.339 |
| Elements Textures | 20 | 2621440 | 1718316 | 1.526 |
| Kodim (Photographic) | 24 | 4718592 | 3074531 | 1.535 |
| Metal Textures | 20 | 2621440 | 1661363 | 1.578 |
| Other Misc | 2 | 1048576 | 734665 | 1.427 |
| PK02 Floor (Misc) | 13 | 1900544 | 1147816 | 1.656 |
| PK02 Trim (Misc) | 5 | 188416 | 122745 | 1.535 |
| PKF Concrete (Misc) | 4 | 131072 | 93659 | 1.399 |
| Plaster Textures | 20 | 2621440 | 1853888 | 1.414 |
| Rock (Misc) | 4 | 2097152 | 1504195 | 1.394 |
| Stone Textures | 20 | 2621440 | 1807316 | 1.450 |
| Tile Textures | 20 | 2621440 | 1823480 | 1.438 |
| Wood Textures | 20 | 2621440 | 1580632 | 1.658 |

[BC1 textures samples](./textures/bc1/)

### BC4 benchmarks

Average compression ratio : 1.2988

| Category | Samples | Uncompressed (bytes) | Avg Compressed (bytes) | Avg Compression Ratio (Arithmetic Mean) |
| :--- | :--- | ---: | ---: | ---: |
| AO | 4 | 524,288 | 451,051 | 1.1624 |
| Displacement | 4 | 524,288 | 367,660 | 1.4260 |
| **Arithmetic Average** | - | - | - | **1.2988** |


[BC4 textures samples](./textures/bc4/)


### BC5 benchmarks

Average compression ratio : 1.187013

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

Crunch  100× 1024×1024 texture: 2.23 s  → **21.97 MB/s**  
Decrunch 100× 1024×1024 texture: 1.44 s  → **36.02 MB/s**

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
