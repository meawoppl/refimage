# Changelog for `refimage`

## Changes

### Unreleased

- `GZIP_1` and `RICE_1` tiles are now encoded by
  [`fitsio-pure`](https://crates.io/crates/fitsio-pure)'s `compress` module, a
  pure-Rust (no_std, wasm32) port of cfitsio's encoders that is checked against
  cfitsio in its CI. The crate's own `rice.rs` and `gzip.rs` are removed, and so is
  the `miniz_oxide` dependency. `HCOMPRESS_1`, the tile grid and the global `f32`
  quantization step are still computed here, so every builder setting behaves as
  before. Output is the same size and decodes to the same values in astropy,
  quantized `f32` included. Lossless `GZIP_1` floats now carry
  `ZQUANTIZ = 'NONE'`, as cfitsio writes them.
- `FitsError::Compression` reports a `fitsio-pure` encoding failure.

### 1.0.0-pre6 (2026-09-03)

- **Breaking:** `PixelType` now has only the six variants the crate can actually
  store and process — `U8`, `U10`, `U12`, `U14`, `U16`, `F32`. The signed and
  wide-integer variants (`I8`/`I16`/`I32`/`I64`/`U32`/`U64`/`F64`) and the
  `PixelStor` impls for `i8`/`i16`/`i32`/`u32`/`f64` are gone; they were never
  representable in a `DynamicImage*` or a pipeline. `PixelType::UnsupportedPixelType`
  can still surface if a future variant is added ahead of kernel support.
- Added `PixelStor` for three new types, `U10` / `U12` / `U14`: `#[repr(transparent)]`
  newtypes over `u16` whose `DEFAULT_MIN/MAX_VALUE` is the sensor's true range
  (`0..=1023` / `4095` / `16383`) instead of the full `u16` one. A `Pipeline` now
  dispatches every per-pixel kernel (`Convert`, `Scale`, `ScalePixels`, luma,
  debayer, `ResizeToFit`) to the matching sub-range type whenever the source
  carries an [`ImageRef::with_bit_depth`] tag, so e.g. converting a `U12`-tagged
  image to `u8` scales against `4095`, not `65535` — previously the tag was
  read back by `pixel_type()` but silently ignored by every kernel, so tagged
  data was scaled and clamped as if it were full-range `u16`. `Op::Convert`'s
  target must now be a real storage type (`U8`/`U16`/`F32`); requesting `U10`/
  `U12`/`U14` as a target is a `PipelineError::ConvertTargetNotStorage` — tag
  the result with `with_bit_depth` instead.
- **Breaking:** the raw-buffer accessors (`as_raw_u8`, `as_raw_u8_checked`,
  `as_slice_u8` / `_u16` / `_f32`, and the `_mut` variants) moved from inherent
  methods to a new `PixelData` trait, implemented for `DynamicImageRef` /
  `DynamicImageOwned` / `GenericImageRef` / `GenericImageOwned` / `GenericImage`
  and the new `DynamicImageView`. `use refimage::PixelData;` to call them.
- **Breaking:** `get_`-prefixed accessors are renamed to drop the prefix, matching
  `Metadata`'s existing style: `GenericImageRef` / `GenericImageOwned` /
  `GenericImage` gain `timestamp()`, `exposure()`, `frame_id()`, `image()`,
  `image_mut()`, `metadata()`, `metadata_mut()`, `key()`; `GenericLineItem` gains
  `value()` / `comment()`; `GenericValue` gains `value_u16()` / `value_f64()` /
  `value_timestamp()` / `value_string()` etc. `set_*` names are unchanged.
- Added `ImageView<'a, T>` and `DynamicImageView<'a>` — a shared *read-only*
  borrow of image samples, produced by `ImageRef::view()` / `ImageOwned::view()` /
  `DynamicImageRef::view()` / `DynamicImageOwned::view()`. `GenericImage` is now
  fully wired: `image()` returns a `DynamicImageView` over whichever half it
  holds, and it gains `metadata_mut()` and `into_owned()`.

### 1.0.0-pre5 (2026-08-28)

- Updated Rust edition to 2024.
- `Gzip` gained a `level` setter (DEFLATE effort `0`..=`9`, default `6`), matching
  the builder style of the other compression types. The FITS writer now streams
  each HDU straight to the `Write` sink — header, then the pixel data one row at
  a time, then padding — instead of building the whole file in a `Vec<u8>` first,
  so writing a large uncompressed image no longer needs a second full-size buffer.
- Added an optional frame ID to `Metadata`, alongside `timestamp` and
  `exposure`. `set_frame_id` takes a `u32`; `frame_id` returns `Option<u32>`.
  `GenericImage`, `GenericImageRef` and `GenericImageOwned` gain matching
  `get_frame_id` / `set_frame_id`. `FRAMEID` is reserved against the extra-item
  map, and a set frame ID is written to FITS as a `FRAMEID` card. `Metadata`
  serialized before this release deserializes with the frame ID unset.
- Added the `pipeline` module: a declarative, cloneable, serializable `Pipeline`
  of `Op`s (debayer, luminance, affine pixel scale, pixel-type conversion, crop,
  ROI, flip, 90°/180°/270° rotation). `Pipeline::compile` validates the chain
  against an `ImageSpec` and pre-allocates all buffers, yielding a `Runner` that
  reuses them across frames. The serial-tiled `Strategy` does zero per-frame heap
  allocation; tiled strategies produce bit-identical output to `Sequential`.
  `Pipeline::apply` runs the chain once and returns an owned image: it accepts any
  of the four image types (via `pipeline::ApplyInput`), and a `GenericImageRef` /
  `GenericImageOwned` input yields a `GenericImageOwned` with its metadata carried
  across unchanged. `Runner::run` / `run_into` accept any `pipeline::Frame`
  (`DynamicImageRef` or `DynamicImageOwned`); `run_into` renders into a
  caller-supplied buffer.
- Added `Op::ResizeToFit` / `Pipeline::resize_to_fit(max_width, max_height, filter)`:
  aspect-preserving resampling to the largest size fitting within a bounding box
  (enlarging when the source is smaller). `pipeline::ResizeFilter` selects
  `Bilinear`, `Bicubic` (Catmull-Rom), or `Lanczos3`; on downscale the kernel is
  widened by the ratio so the result is area-averaged (anti-aliased). It is a
  geometric op — it runs as a whole-frame pass between tiled segments, fused and
  cache-blocked into bands of output rows that fan out over the `rayon` pool
  (each band allocates a small strip, not a full intermediate plane; output is
  independent of band size and thread count). Not valid on a Bayer image
  (`PipelineError::ResizeOnBayer`); debayer first.
- Added `Pipeline::optimize`: an output-preserving peephole cleanup of the op
  list, to run before `compile`. It drops `Op::Nop` stages, folds any run of
  flips and 90°/180°/270° rotations into the shortest equivalent sequence (at
  most two ops — `rotate_90().rotate_90()` becomes `rotate_180()`,
  `flip_horizontal().flip_horizontal()` vanishes), merges nested `Op::Crop`
  stages into one, and drops a redundant second luminance stage or a repeated
  `Op::Convert` to the same type. Arithmetic stages (`Scale`, `ScalePixels`,
  distinct `Convert`s) are left unfused — the rounding and saturation between
  them is observable.
- Tiled plans are now **segmented** at every geometric op (resize, rotation,
  mid-chain crop/ROI, flip): the pixel-op run on the far side of one is retiled
  against its own post-transform dimensions instead of dropping to a single
  whole-frame pass, so `debayer luma → resize → scale convert` tiles on both
  sides. The two full-frame buffers this needs were already allocated for the
  old geo tail, so per-frame allocation is unchanged. `Runner::tiled_pass_count`
  reports how many independently-tiled passes a plan runs. Output stays
  bit-identical to `Sequential`.
- **Breaking:** the pipeline is now the _only_ API for pixel conversions. The
  `Debayer`, `ToLuma`, `ConvertPixelType`, `SelectRoi` and `CopyRoi` traits and
  every per-type impl are removed. `SelectRoi` becomes `Op::Roi`; `CopyRoi`
  becomes `Runner::run_into` with a leading `Op::Roi`.
- **Breaking:** fallible operations now return concrete `thiserror` enums instead
  of `&'static str` / `String`: `ImageError` (construction / byte casts),
  `MetadataError` (key/comment/value validation, reserved & missing keys),
  `SerdeError` (wire-format decode), `ExposureError` (`OptimumExposureBuilder` /
  `CalcOptExp`), and `InteropError` (the `image` feature). `BayerError::InvalidColorSpace`
  now carries a `ColorSpace`.
- **Breaking:** `GenericImageOwned::operate` is removed — `Pipeline::apply` on a
  `GenericImageOwned` now covers it (transform the pixels, keep the metadata). For a
  non-pipeline transform, `GenericImageOwned::with_metadata(img.get_metadata().clone(),
  f(img.get_image())?)`.
- Added the `grow` feature: lets a compiled `Runner` reallocate its buffers on a
  frame-shape change (via `Runner::recompile`) instead of erroring.
- **Breaking:** metadata is now a typed `Metadata` value, not a bare map.
  `GenericImageRef` / `GenericImageOwned` carry
  `Metadata { timestamp: chrono::DateTime<Utc>, exposure: Duration, extra: IndexMap<String, GenericLineItem> }`: `timestamp` and `exposure` are
  mandatory typed fields (both `new` constructors take `(timestamp, exposure, image)`; `image` is now `impl Into<DynamicImage…>`), so `get_timestamp()` is
  infallible and returns a UTC `DateTime` (was `SystemTime`), and `get_exposure()`
  returns `Duration` (`Duration::ZERO` = not applicable) instead of
  `Option<Duration>`. `TIMESTAMP` and `EXPOSURE` are reserved keys.
  `get_metadata()` returns `&Metadata`; the `extra` map is insertion-ordered
  (`indexmap`), so FITS header cards and the serialized metadata are
  deterministic. `GenericValue::SystemTime` becomes
  `GenericValue::Timestamp(DateTime<Utc>)` (`get_value_systemtime` →
  `get_value_timestamp`); pre-1970 timestamps now round-trip. New: `Metadata`,
  `with_metadata` constructors, `set_exposure`, `get_metadata_mut`, and
  `GenericImage*::optimum_exposure(eval, bin)` (sources the reference exposure
  from metadata). `FitsError::{MissingTimestamp, TimestampBeforeEpoch}` removed.
  The `serde`/bincode wire format for metadata changes.
- **Breaking:** `CalcOptExp` / `OptimumExposure` reworked.
  `calc_opt_exp` now takes `&mut self` (was `self`) and both it and
  `OptimumExposure::calculate` return an `OptimumExposureResult { exposure, bin, measured, within_target, clamped }` instead of a bare `(Duration, u16)`.
  `bin` is `u16` in and out (was `u8` in). `OptimumExposure::calculate` drops its
  redundant `len` argument and its `T: Ord` bound, so **`f32` images are now
  supported** (`ExposureError::FloatUnsupported` removed); it uses an `O(n)`
  partial selection instead of a full sort. Empty images now yield
  `ExposureError::EmptyImage` rather than panicking, and `pixel_exclusion`
  near the pixel count no longer underflows. The `pixel_exclusion` hot-pixel
  cutoff is now always applied (previously only when the percentile index fell
  below it).
- `PixelStor::from_f64` / `from_f32` now clamp into the type's range (NaN maps to
  the minimum) _before_ narrowing, so out-of-range / non-finite inputs saturate
  instead of panicking. This makes `Op::Scale` with a large gain or negative
  offset well-defined.
- `DemosaicMethod` now derives `Serialize` / `Deserialize`.
- Added `ColorSpace::channels()`. `ImageRef` / `ImageOwned` no longer cache a
  separate channel count — it is derived from the `ColorSpace` — so the two can
  never disagree.
- Added the `pipeline` example; `test_debayer` / `test_flow` rewritten to the
  pipeline API.
- `refimage` now builds and is tested on `wasm32-unknown-unknown`: added
  `tests/wasm.rs` (a `wasm-bindgen-test` suite over image construction, the
  pipeline, metadata, and serde) and a GitHub Actions job that runs it. `rand`
  is now a non-wasm dev-dependency.
- **Breaking:** the FITS writer is now a pure-Rust implementation with **no
  `cfitsio` linkage**. It is always available (the `fitsio` Cargo feature is
  removed) and compiles on every target, `wasm32` included. `FitsError` is now a
  concrete `#[non_exhaustive]` enum (was a re-export of `fitsio::errors::Error`).
  Compression is chosen with **builder structs** — `Gzip` (`GZIP_1`), `Rice`
  (`RICE_1`), `Hcompress` (`HCOMPRESS_1`), or `FitsCompression::NONE` — passed to the
  write methods as `impl Into<FitsCompression>` (`FitsCompression` is now an opaque
  type; the old unit enum is gone). Each builder exposes settings: `tile_rows(n)` /
  `tile_dims([nx, ny])` for the tile grid (mutually exclusive at compile time via a
  type-state parameter; default one image row per tile, or one plane for `Hcompress`),
  `Quantize` (`level` = cfitsio's `q` / "noise bits", and dither `seed`) for `f32` on
  `Rice` / `Hcompress`, and `scale` / `smooth` for `Hcompress`. `Rice` and `Hcompress`
  are lossless for `u8` / `u16`; `f32` is quantized with
  `ZQUANTIZ = 'SUBTRACTIVE_DITHER_1'` (the dither RNG is reproduced exactly, so
  `astropy` / `cfitsio` reconstruct the same values). `Hcompress` is a pure-Rust port
  of cfitsio's H-transform + binary-quadtree encoder and requires every tile to be at
  least 4×4. `BZIP2_1` was dropped from the plan — neither `astropy.io.fits` nor the
  public `cfitsio` release can read it. New: `FitsWrite::write_fits_to` / `fits_bytes`
  (serialise a single image without a filesystem), `FitsResult`, `FitsCompressionKind`
  (the algorithm tag behind a `FitsCompression`).
  `FitsWriter` is now generic over any `Write` sink (`FitsWriter<W = BufWriter<File>>`);
  the new `create_fits_to(sink, compress)` builds a multi-HDU file in memory (e.g. into
  a `Vec<u8>`), and `FitsWriter::finish` now returns the sink. Multi-channel images are
  now written **planar** (`NAXIS3 = channels`); `write_fits` no longer rewrites the
  path extension. Output is checked against `astropy.io.fits` and `cfitsio` in
  the test suite.
- **Breaking:** `GenericValue`'s numeric variants are collapsed — every integer
  width folds into `GenericValue::Integer(i64)` and both floats into
  `GenericValue::Real(f64)`. `u64` is no longer a supported metadata value type
  (it cannot be promoted to `i64` losslessly). The public `insert_key` /
  `get_value_*` API is unchanged apart from the dropped `get_value_u64`. The
  `serde`/bincode wire format for metadata changes accordingly.
- Dependency cleanup: dropped the unused `itertools` dependency and the
  `color-backtrace` dev-dependency; dropped the unused `aarch64_simd`,
  `wasm_simd` and `extern_crate_alloc` `bytemuck` features. Added `miniz_oxide`
  (pure Rust, for `GZIP_1`), `indexmap` (ordered metadata), a non-wasm `fitsio`
  dev-dependency (cross-checks), and `chrono` (`default-features = false`,
  `["alloc", "serde"]` only — no `clock`, so the caller supplies the timestamp;
  re-exported as `refimage::chrono`).
  Bumped `thiserror` to 2.0 and refreshed the remaining version requirements.
- Internal changes to `fitsio_interop` function signatures.
- Fixed a bug in `CopyRoi` function.

### 1.0.0-pre4 (2025-01-21)

- Added `create_fits` function to library to create an empty FITS file.
- Changed type of the `path` input in `FitsWrite::write_fits` to `AsRef<Path>`.

### 1.0.0-pre3 (2025-01-15)

- Clippy: Elided lifetimes where possible.
- Added `append_fits` function to `FitsWrite` trait to create multi-image FITS files.

### 1.0.0-pre2 (2024-12-18)

- Fixed PixelStor declarations for rare types.

### 1.0.0-pre1 (2024-12-18)

- Added traits `SelectRoi` and `CopyRoi` for image ROI extraction.

### 1.0.0-pre0 (2024-10-23)

- Updated documentation for `GenericLineItem`.
- Bugfix: ImageRef to ImageOwned did not respect data size.
- Bugfix: GenericImageRef FitsWrite did not respect data size.
- Bugfix: DATE-OBS correctly reports time up to microseconds.

### 0.13.1 (2024-10-20)

- Save timestamp and duration with a float second value with the original key,
  and additional integer values of seconds and nanoseconds, with keys suffixed
  by 's' and 'ns'.

### 0.13.0 (2024-10-16)

- Removed default features from `image` crate.
- Removed default features from `chrono` crate.
- Removed feature `serde_flate`. Serialized image does not compress data.

### 0.12.2 (2024-10-10)

- Added `is_bayer` method to `ColorSpace`.

### 0.12.1 (2024-10-09)

- Bugfix: `rayon` imports in `src/coreimpls.rs` not hidden behind the rayon feature flag.

### 0.12.0 (2024-10-08)

- Removed alpha channel support to allow in-place luminance calculations.
- `ToLuma` trait operations are now in-place, and require mutable borrows.

### 0.11.0 (2024-09-24)

- Added `flip_*` methods to `BayerShift` trait.

### 0.10.0 (2024-09-24)

- Complete API overhaul. `DataStor` type has been removed.
- `ImageData` and similar are renamed to `*ImageRef`, and only hold data backed by some storage.
- `ImageOwned` and derivatives own the image data.
- `GenericImageRef` and `DynamicImageRef` are only serializable. `*ImageRef` do not implement clone.
- All methods that mutate an image (e.g. luminance calculation, alpha channel addition/removal, debayering) either fail when mutation is not necessary, or return an owned image on success.

### 0.8.2 (2024-09-22)

- Added missing `InsertKey` implementation for `ColorSpace`.
- Added [crates.io](https://crates.io/crates/refimage) documentation support for features.

### 0.8.0 (2024-09-22)

- Added `from_u8` methods to `ImageData` and `ImageOwned` to get these images from a slice of `u8` data in machine-endian order.

### 0.7.0 (2024-09-21)

- Alpha channel support has been added.
- API change: Bayer patterns moved to `BayerPattern` enum, `ColorSpace` now contains a variant `ColorSpace::Bayer(BayerPattern)`.
- Enhancement: `BayerPattern` implements the `BayerShift` trait to allow `BayerPattern` selection for smaller, shifted ROIs.

### 0.6.2 (2024-09-21)

- Bugfix: `into_luma()` errored on a grayscale image.

### 0.6.1 (2024-09-21)

- Added `EXPOSURE_KEY` to standardize setting of exposure time.

### 0.6.0 (2024-09-21)

- Added `CalcOptExp` trait to calculate optimum exposure for a given `DynamicImageData` or `ImageData`.
  This trait is also implemented for the equivalent Owned types.
- Changed the signature of `OptimumExposure::calculate` function to accept a mutable reference to an array instead of a `Vec<T>`.

### 0.5.1 (2024-09-14)

- Refactored `rayon` code, reduced code duplication.
- Added feature documentation.

### 0.5.0 (2024-09-13)

- Optimum exposure calculation internally casts to `f32`.
- Added `ColorSpace` metadata, removed `repr(u8)` requirement since `serde` is now required.
- Added `to_f32` and `from_f32` functions for `PixelStor` trait.
- Added `cast_u8` for `PixelStor`, allowing casting any value to [0, 255] range.
- Added owned data counterparts to `GenericImage`, `DynamicImageData` and `ImageData`, along with
  accompanying `To` and `From` implementations.
- Added `into_u8` (parallelized with `rayon`) functions to convert all images to owned images with `u8` pixels.

### 0.4.0 (2024-08-27)

- Exposed a `GenericValue` enum, with `From` and `TryInto` implementations to concrete types.
  `GenericLineItem` now does not include functions to access the underlying concrete types.
- Added optimum exposure calculation support.
- Added luminance conversion support.

### 0.3.0 (2024-08-23)

- Changed the signature of the closure taken by the `GenericImage::operate` function, and
  changed the behavior of the function to return a copy of the original image. Updated the
  relevant test.
- Updated the example to show the use of `DynamicImageData` and `GenericImage`.

### 0.2.0 (2024-08-21)

- Moved the `debayer` method to its own trait, and removed public `debayer` methods
  from `DynamicImageData` and `ImageData`.
- Renamed enum `Demosaic` to `DemosaicMethod`.
- Implemented `Debayer::debayer` for `GenericImage`.
- Implemented `GenericImage::operate()` that allows replacement of the underlying
  image while conserving the metadata.
- Moved the `GenericImage::write_fits` function (enabled with `fitsio` feature) to the `FitsWrite` trait.
- Added custom FITS compression support through `FitsCompression::Custom`.
- Formatted all files and documentation.
- Added this change log.

### 0.1.1 (2024-08-20)

- Updated documentation.

## Pre-publish steps

### Check documentation

First you add this to your `Cargo.toml`:

```toml
[package.metadata.docs.rs]
all-features = true
rustdoc-args = ["--cfg", "docsrs"]
```

This will make docs.rs compile your documentation with `RUSTFLAGS="--cfg docsrs"` set. We can use this to conditionally enable nightly features when building documentation, but not when building it normally.

Now you add this to your `lib.rs` file:

`#![cfg_attr(docsrs, feature(doc_cfg))]`
Then, on anything that needs such an annotation, you add the following

`#[cfg_attr(docsrs, doc(cfg(feature = "rt")))]`
This supports all of the syntax that `#[cfg(...)]` does, including stuff like platform-specifiers and any/all. Tokio often sets this annotation using the macros in [`src/macros/cfg.rs`](https://github.com/tokio-rs/tokio/blob/master/tokio/src/macros/cfg.rs). Note that the annotation works recursively, so putting it on a module also affects the stuff inside (sometimes, it's a bit flaky).

Finally, when building the documentation locally, you run the following command:

```sh
$ RUSTDOCFLAGS="--cfg docsrs" cargo +nightly doc --all-features
```

### Run `cargo clippy`

```sh
$ cargo clippy --all-features
```

### Run `cargo fmt`
