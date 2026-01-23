# bc_crunch

![C Standard](https://img.shields.io/badge/C-C99-F4A261)
[![Build Status](https://github.com/Geolm/bc_crunch/actions/workflows/c.yml/badge.svg)](https://github.com/geolm/bc_crunch/actions)


`bc_crunch` is a tiny, dependency-free C99 library for *lossless* compression of GPU-compressed texture blocks **BC1, BC4, BC3, and BC5**.

## Huffman Branch

This branch introduces a Huffman-based backend aimed at improving decompression speed while preserving a strong compression ratio.  
At the moment, it is implemented for **BC1 only**.

Despite the focus on fast decoding, the average compression ratio remains **better than zlib**.

### Compression Ratio (BC1)

- bc_crunch: **1.505671**
- zlib: 1.436144

### Performance

Benchmarks were run on my M1 Pro by compressing and decompressing a 1024×1024 BC1 texture 100 times.

Compression
- Time: 5.829384 s
- Throughput: **8.577236 MB/s**

Decompression
- Time: 0.353572 s
- Throughput: **141.413896 MB/s**

Decompression throughput is **four times (4x)** higher than the arithmetic encoding implementation (35MB/s).
