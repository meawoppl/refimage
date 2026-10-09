# `refimage`
[![crates.io](https://img.shields.io/crates/v/refimage)](https://crates.io/crates/refimage)
[![Documentation](https://docs.rs/refimage/badge.svg)](https://docs.rs/refimage)

## A Serializable Image Container

This crate provides a type-erased image container (`GenericImageRef`), backed by a contiguous slice (owned or referenced) of primitive (`u8`, `u16`, `f32`) pixels, with arbitrary color space (grayscale, Bayer pattern, RGB, ...) and color channels support (max. 255).
Image sizes are limited to 65536 &times; 65536 for practical reasons.

`GenericImageRef` and `GenericImageOwned` always carry a `Metadata` block: a mandatory typed core (`timestamp: chrono::DateTime<Utc>`, `exposure: Duration` — `Duration::ZERO` meaning "not applicable") plus an insertion-ordered map of extra (`key`, `value`) pairs with optional comments. Metadata keys are case-insensitive strings up to 80 characters (`TIMESTAMP` / `EXPOSURE` are reserved); values are integers (stored as `i64`), 32-/64-bit floats (stored as `f64`), strings, `ColorSpace`, `std::time::Duration`, or a UTC `chrono::DateTime`. `chrono` is re-exported as `refimage::chrono`.

`GenericImageRef` supports serialization, and can be saved to the open [Flexible Image Transport System (FITS)](https://fits.gsfc.nasa.gov/fits_standard.html) through the `FitsWrite` trait. It supports uncompressed output (`FitsCompression::NONE`) and the tile-compression convention via the builder structs `Gzip` (lossless, any pixel type), `Rice` and `Hcompress` (both lossless for `u8`/`u16`; `f32` is quantized with `SUBTRACTIVE_DITHER_1`).
Each builder exposes its settings, e.g. tile shape (`tile_rows` / `tile_dims`), quantization level, HCOMPRESS `scale` / `smooth`.

`GenericImageRef` serializes to an internal representation which allows <i>de</i>serialization to `GenericImageOwned`.
`GenericImageOwned` contains owned data, and is freely <i>ser</i>ialized-<i>de</i>serialized. 

## The path to a `GenericImageRef`

A `GenericImageRef` is obtained from a `ImageRef` object, created with the
appropriate, contiguous, backing storage and image format:
```rust
use refimage::{BayerPattern, ImageRef, DemosaicMethod, DynamicImageRef, GenericImageOwned};
use refimage::pipeline::Pipeline;
use refimage::chrono::DateTime;
use std::time::Duration;

let mut data = vec![0u8; 256]; // this is the backing store
// acquire(&mut data); // this function populates the backing store with the image pixels
let img = ImageRef::new(&mut data, 16, 16, BayerPattern::Grbg.into()).unwrap(); // 16x16 image backed by the vector
let img = DynamicImageRef::from(img); // convert the `ImageRef` object to `DynamicImageRef`
// All pixel conversions go through a `Pipeline`; `apply` runs it once and owns the result.
let img = Pipeline::new().debayer(DemosaicMethod::Nearest).apply(&img).expect("Could not debayer");
// A UTC timestamp and an exposure (`Duration::ZERO` if not applicable) are mandatory.
// In application code the timestamp is usually `chrono::Utc::now()`.
let ts = DateTime::from_timestamp(1_700_000_000, 0).unwrap();
let mut img = GenericImageOwned::new(ts, Duration::from_millis(20), img);
img.insert_key("CAMERA", ("Rust Test Program", "Name of the camera used to capture the image")).unwrap();
let json = serde_json::to_string(&img).unwrap(); // serialize the image to JSON
let rimg: GenericImageOwned = serde_json::from_str(&json).unwrap(); // deserialize to GenericImageOwned
assert_eq!(&img, &rimg); // Confirm that deserialized image matches the original
```
`*ImageRef` to minimizes unnecessary allocations. Running a
`GenericImageRef` through `Pipeline::apply` returns a `GenericImageOwned`, with its
metadata carried across unchanged.

## `GenericImageOwned` and other `ImageOwned` types
An image can be loaded using the [`image`](https://crates.io/crates/image) crate from disk, by enabling the `image` feature:
```rust,ignore
use refimage::DynamicImageOwned;
use image::open;

let img = open("/path/to/image.png").expect("Could not load image");
let img = DynamicImageOwned::try_from(img).expect("Could not convert image");
```

## Loading and storing a `GenericImageRef`
A `GenericImageRef` can be stored in a standard format such as [bincode](https://crates.io/crates/bincode) — which follows trivially from its `serde` implementation — or written to a FITS file through the `FitsWrite` trait. The FITS writer is pure Rust and always available.

```rust,no_run
use refimage::{FitsCompression, FitsWrite, GenericImageRef, Gzip, Rice, Quantize};
use std::path::Path;
let img: GenericImageRef = { todo!() }; // obtain a GenericImageRef

// RICE_1 tile compression; overwrite if it exists.
img.write_fits(Path::new("/path/to/image.fits"), Rice::new(), true)
    .expect("Could not write FITS file.");

// Each algorithm is a builder: tile shape, quantization level, HCOMPRESS scale, ...
img.write_fits(
    Path::new("/path/to/image.fits"),
    Rice::new().tile_rows(16).quantize(Quantize::new().level(16.0)),
    true,
).unwrap();

// Or straight to bytes — handy on wasm, where there is no filesystem.
let bytes: Vec<u8> = img.fits_bytes(Gzip::new()).unwrap();
let plain: Vec<u8> = img.fits_bytes(FitsCompression::NONE).unwrap();
```

## Conversions

Debayering, luminance, pixel-type conversion, affine pixel scaling, cropping, ROI
extraction, flips, 90° rotations and aspect-preserving resize (`resize_to_fit`,
with a `Bilinear` / `Bicubic` / `Lanczos3` filter) are all `Op`s on a
[`Pipeline`](#reusable-processing-pipelines).
Build a `Pipeline` and call `apply` (a `GenericImage*` input keeps its metadata on the result), or compile a `Runner` for a frame stream.

## Reusable processing pipelines
Every conversion is a step in a `Pipeline`.
For a one-off conversion, `Pipeline::apply(&img)` compiles, runs once, and returns an owned image.
For acquisition loops that run the same chain over every frame, compile a `Runner` once and reuse its buffers:
```rust
use refimage::{BayerPattern, ColorSpace, DemosaicMethod, DynamicImageRef, ImageRef, PixelType};
use refimage::pipeline::{ImageSpec, Pipeline, ResizeFilter, Strategy};

// A declarative, cloneable, serializable recipe.
let recipe = Pipeline::new()
    .crop(0, 0, 64, 64)
    .debayer(DemosaicMethod::Linear)
    .to_luma()
    .scale(1.2, 4.0)
    .convert(PixelType::U8)
    .resize_to_fit(48, 48, ResizeFilter::Lanczos3)
    .flip_vertical();

// Compile against a concrete frame format: validates the chain and
// pre-allocates every buffer it will ever need.
let spec = ImageSpec::new(64, 64, ColorSpace::Bayer(BayerPattern::Rggb), PixelType::U16);
let mut runner = recipe.compile(spec, Strategy::tiled_parallel(16)).unwrap();

// Run it per frame; the result borrows the runner's internal buffer.
let mut frame = vec![0u16; 64 * 64];
let img = DynamicImageRef::from(
    ImageRef::new(&mut frame, 64, 64, ColorSpace::Bayer(BayerPattern::Rggb)).unwrap(),
);
let luma = runner.run(&img).unwrap();
assert_eq!(luma.color_space(), ColorSpace::Gray);
assert_eq!((luma.width(), luma.height()), (48, 48));
```
`Strategy::Sequential` uses one ping-pong buffer pair; `Strategy::tiled`/`tiled_parallel`/`tiled_2d`/`tiled_2d_parallel` cut the frame into cache-sized tiles. A serial-tiled `Runner` does **zero** heap allocation per frame once compiled. A geometric op (rotation, mid-chain crop, `resize_to_fit`) splits the chain into segments that each tile independently — `debayer luma → resize → scale convert` tiles on both sides of the resize (`Runner::tiled_pass_count`). See the module docs for the details and for `Runner::recompile` / the `grow` feature.

# FITS
`GenericImageRef` / `GenericImageOwned` can be stored to FITS via the `FitsWrite` trait.
It writes uncompressed images (`FitsCompression::NONE`) and the tile-compression convention through the builder structs `Gzip`, `Rice` and `Hcompress` (anything `Into<FitsCompression>`).
`Gzip` is lossless for any pixel type; `Rice` and `Hcompress` are lossless for `u8`/`u16` and quantize `f32` with `SUBTRACTIVE_DITHER_1`.
`Gzip` and `Rice` tiles are encoded by the pure-Rust [`fitsio-pure`](https://crates.io/crates/fitsio-pure) crate, whose encoders are checked against cfitsio.
Tile shape is set with `tile_rows(n)` or `tile_dims([nx, ny])`, defaulting to one image row per tile for `Rice`/`Gzip` and one whole channel plane for `Hcompress` (the H-transform + quadtree coder, min. tile 4×4); `Rice`/`Hcompress` also take a `Quantize` (level / dither seed) and `Hcompress` a `scale` / `smooth`.
Multi-channel images are stored planar, `u16` carries `BZERO = 32768`, and metadata is stored in header cards (`HIERARCH` for long keys, `CONTINUE` for long strings).
`FitsWrite::fits_bytes` returns a single image as a `Vec<u8>` (in-memory), and `create_fits_to(sink, compression)` builds a multi-HDU file into a `Write` sink (`create_fits` is the on-disk equivalent).

# Optional Features
Features are available to extend the functionalities of the core `refimage` data types:
- `rayon`: Parallelizes luma conversion / demosaic / cast inside the `pipeline`, and enables the parallel `pipeline` strategies (<b>enabled</b> by default).
- `image`: Enables `TryFrom` conversions between `image::DynamicImage` and `refimage::DynamicImageRef`, `refimage::DynamicImageOwned` (<b>disabled</b> by default).
- `grow`: Lets a compiled pipeline `Runner` reallocate its buffers on a frame-shape change instead of erroring (<b>enabled</b> by default).

# WebAssembly
The core of `refimage` — image construction, `pipeline`, metadata, and `serde`, runs on `wasm32-unknown-unknown`.
The pipeline falls back to its serial kernels when `rayon` is disabled. `tests/wasm.rs` is a `wasm-bindgen-test` suite covering this implementation.
To run the tests, install `wasm-bindgen-cli` and execute:

```sh
cargo install wasm-bindgen-cli
cargo test --target wasm32-unknown-unknown --no-default-features --tests
```

(`--tests` skips the doc-tests, which rustdoc does not cross-compile.)

FITS export is available on WASM targets.
`FitsWrite::fits_bytes` returns a single image as a `Vec<u8>`, and `create_fits_to(Vec::new(), compression)` builds a multi-HDU file in memory.