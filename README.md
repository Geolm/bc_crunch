# bc_crunch

![C Standard](https://img.shields.io/badge/C-C99-F4A261)
[![Build Status](https://github.com/Geolm/bc_crunch/actions/workflows/c.yml/badge.svg)](https://github.com/geolm/crunch/actions)


`bc_crunch` is a tiny, dependency-free C99 library for *lossless* compression of GPU-compressed texture blocks **BC1, BC4, BC3, and BC5**.

- ~700 LOC total (single-file encoder/decoder, no build system tricks)
- No malloc, no external libs, no threads, no SIMD requirement  
- Deterministic, bit-exact reconstruction
- Fully tested with bytes-per-bytes comparison
- Focused on production textures: albedo, masks, normals, heightmaps, etc.
- GPU-ready output: decompressed blocks are written in standard BC1/BC4/BC3/BC5 format, ready to be uploaded directly to GPU memory (or directly written in shared memory)
- No extra memory for decompression: only the encoder needs a temporary buffer; decoding writes straight to the output buffer.

This is *not* another general-purpose compressor. `bc_crunch` is specialized for already-compressed GPU formats — it exploits the internal structure of BC1/BC4 blocks, spatial patterns, endpoint deltas, popcount distances, and Morton-ordered indices to achieve significant size reductions with extremely low CPU cost.

---

## Benchmarks

Compression ratio naturally depends on the input content. Repetitive patterns, smooth gradients, and large uniform regions compress very well, while randomness or high-frequency noise significantly reduces efficiency (e.g., dirt textures are typically very noisy).

### BC1 benchmarks

Average compression ratio : 1.493586

| Category       | # Samples | Original BC1 Size (bytes) | Avg Crunched Size (bytes) | Avg Compression Ratio |
|----------------|-----------|---------------------------|---------------------------|---------------------|
| Kodak photos   | 24        | 196,608                   | 136,254                   | 1.44×               |
| Dirt           | 20        | 131,072                   | 105,386                   | 1.24×               |
| Wood           | 20        | 131,072                   | 79,306                    | 1.65×               |
| Elements       | 20        | 131,072                   | 90,412                    | 1.45×               |
| Brick          | 20        | 131,072                   | 88,212                    | 1.49×               |
| Metal          | 20        | 131,072                   | 87,132                    | 1.51×               |
| Plaster        | 20        | 131,072                   | 99,118                    | 1.33×               |
| Stone          | 20        | 131,072                   | 97,014                    | 1.35×               |
| Tile           | 20        | 131,072                   | 101,087                   | 1.30×               |
| Blaz Tree      | 7         | 32,768–32,768             | 16,900                    | 1.84×               |
| SW Tree        | 7         | 20,480–32,768             | 14,227                    | 2.08×               |
| pkf Concrete   | 4         | 32,768                    | 24,583                    | 1.33×               |
| pk02 Floor     | 20        | 32,768–524,288            | 100,604                   | 1.83×               |
| pk02 Trim      | 5         | 4,096–131,072             | 11,777                    | 1.89×               |
| Rock           | 4         | 524,288                   | 412,365                   | 1.27×               |
| Brick Large    | 2         | 524,288                   | 408,746                   | 1.28×               |

[BC1 textures samples](./textures/bc1/)

### BC4 benchmarks

Average compression ratio : 1.211644

| Category      | Samples | Uncompressed (bytes) | Avg Compressed (bytes) | Avg Compression Ratio |
|---------------|----------|-----------------------|--------------------------|------------------------|
| AO            | 4        | 524,288               | 472,594                 | 1.109384               |
| Displacement  | 4        | 524,288               | 400,749                 | 1.308269               |

[BC4 textures samples](./textures/bc4/)

---

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


## How to build tests

* open a terminal in the folder
* mkdir build
* cd build
* cmake -DCMAKE_BUILD_TYPE=Release ..
* cmake --build .
* ./test

More tests will be added soon.

## FAQ

#### How can I compress a texture with all mipmaps

Right now you need to call bc_crunch separately for each mip level. I may add a helper later, but I want to keep the encoder simple and stateless, without any mip-level awareness or cross-mip data.
