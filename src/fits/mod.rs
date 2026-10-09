//! Writing images to the Flexible Image Transport System (FITS) format.
//!
//! It produces uncompressed image HDUs and the FITS
//! tile-compression convention (`GZIP_1`, `RICE_1`, `HCOMPRESS_1`).
//! Output is verified to be compatible with `astropy.io.fits`
//! and `cfitsio` in the test suite.
//!
//! # Entry points
//!
//! [`FitsWrite`] is implemented for [`GenericImageRef`](crate::GenericImageRef),
//! [`GenericImageOwned`](crate::GenericImageOwned) and
//! [`GenericImage`](crate::GenericImage). It writes one image to a path
//! ([`write_fits`](FitsWrite::write_fits)), to any [`Write`] sink
//! ([`write_fits_to`](FitsWrite::write_fits_to)), or to an owned buffer
//! ([`fits_bytes`](FitsWrite::fits_bytes)). A [`FitsWriter`], created by
//! [`create_fits`] or [`create_fits_to`], collects several images into one multi-HDU
//! file with [`append_fits`](FitsWrite::append_fits).
//!
//! # Compression
//!
//! The `compress` argument of every write method accepts any type that converts into
//! [`FitsCompression`]: [`FitsCompression::NONE`] for an uncompressed HDU, or one of
//! the builders [`Gzip`], [`Rice`], [`Hcompress`]. Each builder carries its own
//! settings:
//!
//! | Setting | Applies to | Default |
//! |---|---|---|
//! | [`tile_rows`](Rice::tile_rows) / [`tile_dims`](Rice::tile_dims) — the tile grid | all | one image row per tile (`Gzip`, `Rice`); one channel plane (`Hcompress`) |
//! | [`level`](Gzip::level) — DEFLATE effort, `0`..=`9` | `Gzip` | `6` (zlib / cfitsio default) |
//! | [`Quantize`] — `f32` quantization step and dither seed | `Rice`, `Hcompress` | step from an image-noise estimate ([`level`](Quantize::level) 4.0) |
//! | [`scale`](Hcompress::scale) — H-transform divisor | `Hcompress` | `0`, lossless |
//! | [`smooth`](Hcompress::smooth) — decompression smoothing flag | `Hcompress` | off |
//!
//! [`tile_rows`](Rice::tile_rows) and [`tile_dims`](Rice::tile_dims) are mutually
//! exclusive. Calling either changes the builder's type so the other is no longer in
//! scope.
//!
//! `Gzip` is lossless for every pixel type. `Rice` and `Hcompress` are lossless for
//! `u8` and `u16`; `f32` images are quantized (`ZQUANTIZ = 'SUBTRACTIVE_DITHER_1'`),
//! a lossy step whose error is bounded by [`Quantize::level`].
//!
//! # File structure
//!
//! - An uncompressed single image is the primary HDU. A compressed single image, and
//!   every file produced through a [`FitsWriter`], begins with an empty primary HDU
//!   followed by one image extension per image.
//! - Multi-channel images are stored one plane after another (`NAXIS3 = channels`).
//! - `u16` data is offset to signed values and written with `BZERO = 32768`, so a
//!   reader reports it as unsigned 16-bit.
//!
//! # Metadata
//!
//! The image timestamp is written as `DATE-OBS`. A non-zero exposure is written as the
//! `EXPOSURE_S` / `EXPOSURE_NS` integer pair together with `EXPOSURE` in seconds. A
//! frame ID, when set, is written as `FRAMEID`. Each
//! entry of the image's [`Metadata`](crate::Metadata) then becomes a header card, in
//! insertion order, using the `HIERARCH` convention for keywords longer than eight
//! characters and `CONTINUE` for string values longer than one card. A
//! [`Duration`](std::time::Duration) or timestamp value is written as a `<KEY>_S` /
//! `<KEY>_NS` integer pair with a convenience card in seconds or ISO-8601. A metadata
//! keyword that would collide with a structural FITS keyword is rejected with
//! [`FitsError::ReservedKeyword`].

mod card;
mod compress;
mod config;
mod datetime;
mod hcompress;
mod hdu;
mod quantize;

use chrono::{DateTime, Utc};
use std::fs::OpenOptions;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use thiserror::Error;

use crate::{GenericLineItem, GenericValue, Metadata, EXPOSURE_KEY, FRAMEID_KEY};

use card::Header;
use config::Method;
use fitsio_pure::compress::{compress_image_hdu, TileCompression};
use hdu::{bayer_pattern, colorspace_str, ImageView, Planar};

pub use config::{
    AutoTile, DitherSeed, FitsCompression, FitsCompressionKind, FixedTile, Gzip, Hcompress,
    Quantize, Rice,
};

/// Errors produced while writing a FITS file.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum FitsError {
    /// Underlying I/O failure.
    #[error(transparent)]
    Io(#[from] std::io::Error),
    /// A metadata key collides with a structural FITS keyword.
    #[error("metadata key {0:?} is a reserved FITS keyword")]
    ReservedKeyword(String),
    /// A metadata key does not fit in a FITS keyword, even as `HIERARCH`.
    #[error("metadata key {0:?} is too long for a FITS keyword")]
    KeywordTooLong(String),
    /// A metadata string value is too long to encode.
    #[error("metadata string value is too long")]
    MetadataValueTooLong,
    /// `HCOMPRESS` was requested with a tile (or whole image) smaller than 4×4.
    #[error("HCOMPRESS needs every tile to be at least 4x4 pixels")]
    HcompressTooSmall,
    /// A tile specification could not be applied to the image.
    #[error("invalid tile specification: {0}")]
    InvalidTiling(String),
    /// `fitsio-pure` could not tile-compress the image.
    #[error("tile compression failed: {0}")]
    Compression(String),
}

/// `Result` alias for [`FitsError`].
pub type FitsResult<T> = Result<T, FitsError>;

/// A multi-HDU FITS file open for writing, over any [`Write`] sink (a file by default,
/// or an in-memory buffer — handy on `wasm32`, which has no filesystem).
///
/// Append image HDUs with [`FitsWrite::append_fits`], then call [`FitsWriter::finish`]
/// to flush and recover the sink.
pub struct FitsWriter<W: Write = BufWriter<std::fs::File>> {
    out: W,
    compress: FitsCompression,
    n_hdus: usize,
}

impl<W: Write> FitsWriter<W> {
    /// Flush the sink and return it (a `Vec<u8>` now holds the whole file, a file
    /// handle is fully written).
    pub fn finish(mut self) -> FitsResult<W> {
        self.out.flush()?;
        Ok(self.out)
    }

    /// Number of HDUs written so far (the empty primary counts as one).
    pub fn hdu_count(&self) -> usize {
        self.n_hdus
    }
}

/// Create a FITS file on disk with an empty primary HDU, ready for
/// [`FitsWrite::append_fits`]. `compress` is the compression each appended image uses.
pub fn create_fits<P: AsRef<Path>, C: Into<FitsCompression>>(
    path: P,
    compress: C,
    overwrite: bool,
) -> FitsResult<FitsWriter> {
    let out = BufWriter::new(open(path.as_ref(), overwrite)?);
    create_fits_to(out, compress)
}

/// Like [`create_fits`], but writes to any [`Write`] sink instead of a path — e.g. a
/// `Vec<u8>` for an in-memory multi-HDU FITS file. The empty primary HDU is written
/// immediately.
pub fn create_fits_to<W: Write, C: Into<FitsCompression>>(
    mut sink: W,
    compress: C,
) -> FitsResult<FitsWriter<W>> {
    sink.write_all(&primary_empty())?;
    Ok(FitsWriter {
        out: sink,
        compress: compress.into(),
        n_hdus: 1,
    })
}

fn open(path: &Path, overwrite: bool) -> FitsResult<std::fs::File> {
    if path.is_dir() {
        return Err(FitsError::Io(std::io::Error::new(
            std::io::ErrorKind::IsADirectory,
            "path is a directory",
        )));
    }
    Ok(OpenOptions::new()
        .write(true)
        .create(true)
        .create_new(!overwrite)
        .truncate(true)
        .open(path)?)
}

/// Write an image, with its metadata, to a FITS file.
///
/// `compress` accepts anything that converts into [`FitsCompression`]: the bare builders
/// [`Gzip`] / [`Rice`] / [`Hcompress`], a built [`FitsCompression`] value, or
/// [`FitsCompression::NONE`] for an uncompressed HDU.
pub trait FitsWrite {
    /// Write the image to `path` (created, or truncated if `overwrite`). Returns `path`.
    fn write_fits<P: AsRef<Path>, C: Into<FitsCompression>>(
        &self,
        path: P,
        compress: C,
        overwrite: bool,
    ) -> FitsResult<PathBuf>;

    /// Write the image to any [`Write`] sink (no filesystem needed).
    fn write_fits_to<W: Write, C: Into<FitsCompression>>(
        &self,
        sink: W,
        compress: C,
    ) -> FitsResult<()>;

    /// Serialise the image to an in-memory FITS file.
    fn fits_bytes<C: Into<FitsCompression>>(&self, compress: C) -> FitsResult<Vec<u8>> {
        let mut buf = Vec::new();
        self.write_fits_to(&mut buf, compress)?;
        Ok(buf)
    }

    /// Append the image as a new HDU of an open [`FitsWriter`] (file- or memory-backed).
    fn append_fits<W: Write>(&self, writer: &mut FitsWriter<W>) -> FitsResult<()>;
}

fn write_one<W: Write>(
    view: &ImageView<'_>,
    meta: &Metadata,
    compress: &FitsCompression,
    mut sink: W,
    lead_primary: bool,
) -> FitsResult<()> {
    match &compress.0 {
        Method::None => write_uncompressed(view, meta, &mut sink, lead_primary)?,
        method => {
            if lead_primary {
                sink.write_all(&primary_empty())?;
            }
            match method {
                Method::Hcompress {
                    tiling,
                    scale,
                    smooth,
                    quantize,
                } => {
                    let c = compress::hcompress(view, tiling, *scale, *smooth, quantize)?;
                    write_hcompress(view, meta, &c, &mut sink)?;
                }
                _ => write_tiled(view, meta, method, &mut sink)?,
            }
        }
    }
    Ok(())
}

macro_rules! impl_fitswrite {
    ($t:ty) => {
        impl FitsWrite for $t {
            fn write_fits<P: AsRef<Path>, C: Into<FitsCompression>>(
                &self,
                path: P,
                compress: C,
                overwrite: bool,
            ) -> FitsResult<PathBuf> {
                let path = path.as_ref().to_path_buf();
                let file = open(&path, overwrite)?;
                self.write_fits_to(BufWriter::new(file), compress)?;
                Ok(path)
            }

            fn write_fits_to<W: Write, C: Into<FitsCompression>>(
                &self,
                sink: W,
                compress: C,
            ) -> FitsResult<()> {
                let compress = compress.into();
                let view = ImageView::from_ref_like(self.image());
                write_one(&view, self.metadata(), &compress, sink, true)
            }

            fn append_fits<Sink: Write>(&self, writer: &mut FitsWriter<Sink>) -> FitsResult<()> {
                let view = ImageView::from_ref_like(self.image());
                write_one(
                    &view,
                    self.metadata(),
                    &writer.compress,
                    &mut writer.out,
                    false,
                )?;
                writer.n_hdus += 1;
                Ok(())
            }
        }
    };
}

impl_fitswrite!(crate::GenericImageRef<'_>);
impl_fitswrite!(crate::GenericImageOwned);

impl FitsWrite for crate::GenericImage<'_> {
    fn write_fits<P: AsRef<Path>, C: Into<FitsCompression>>(
        &self,
        path: P,
        compress: C,
        overwrite: bool,
    ) -> FitsResult<PathBuf> {
        let compress = compress.into();
        match self {
            crate::GenericImage::Ref(i) => i.write_fits(path, compress, overwrite),
            crate::GenericImage::Own(i) => i.write_fits(path, compress, overwrite),
        }
    }

    fn write_fits_to<W: Write, C: Into<FitsCompression>>(
        &self,
        sink: W,
        compress: C,
    ) -> FitsResult<()> {
        let compress = compress.into();
        match self {
            crate::GenericImage::Ref(i) => i.write_fits_to(sink, compress),
            crate::GenericImage::Own(i) => i.write_fits_to(sink, compress),
        }
    }

    fn append_fits<W: Write>(&self, writer: &mut FitsWriter<W>) -> FitsResult<()> {
        match self {
            crate::GenericImage::Ref(i) => i.append_fits(writer),
            crate::GenericImage::Own(i) => i.append_fits(writer),
        }
    }
}

/// An empty primary HDU (`NAXIS = 0`), block-padded.
fn primary_empty() -> Vec<u8> {
    let mut h = Header::new();
    h.logical("SIMPLE", true, Some("conforms to FITS standard"))
        .unwrap();
    h.integer("BITPIX", 8, None).unwrap();
    h.integer("NAXIS", 0, None).unwrap();
    h.logical("EXTEND", true, None).unwrap();
    h.finish()
}

/// A plain (uncompressed) image HDU, streamed to `sink` (header, then the pixel
/// data one row at a time, then block padding — the image is never buffered
/// whole).
fn write_uncompressed<W: Write>(
    view: &ImageView<'_>,
    meta: &Metadata,
    sink: &mut W,
    primary: bool,
) -> FitsResult<()> {
    let axes = view.axes();
    let mut h = Header::new();

    if primary {
        h.logical("SIMPLE", true, Some("conforms to FITS standard"))?;
        h.integer("BITPIX", view.bitpix(), None)?;
        h.integer("NAXIS", axes.len() as i64, None)?;
        for (i, &n) in axes.iter().enumerate() {
            h.integer(&format!("NAXIS{}", i + 1), n as i64, None)?;
        }
        h.logical("EXTEND", true, None)?;
    } else {
        h.string("XTENSION", "IMAGE", Some("image extension"))?;
        h.integer("BITPIX", view.bitpix(), None)?;
        h.integer("NAXIS", axes.len() as i64, None)?;
        for (i, &n) in axes.iter().enumerate() {
            h.integer(&format!("NAXIS{}", i + 1), n as i64, None)?;
        }
        h.integer("PCOUNT", 0, None)?;
        h.integer("GCOUNT", 1, None)?;
    }

    if let Some(bz) = view.bzero() {
        h.integer("BZERO", bz, Some("offset for unsigned integers"))?;
        h.integer("BSCALE", 1, None)?;
    }
    if let Some(bits) = view.adc_bits {
        h.integer("BITADC", bits as i64, Some("meaningful ADC bits per sample"))?;
    }

    write_common_cards(&mut h, view, meta.timestamp())?;
    write_metadata(&mut h, meta)?;

    sink.write_all(&h.finish())?;
    view.write_native_be(sink)?;
    card::pad_writer(sink, view.native_be_len())?;
    Ok(())
}

/// A `GZIP_1` or `RICE_1` tile-compressed image HDU (`BINTABLE` with `ZIMAGE = T`),
/// written to `sink`. `fitsio_pure::compress` writes the table and its compression
/// keywords; this crate's own cards follow them, before `END`.
fn write_tiled<W: Write>(
    view: &ImageView<'_>,
    meta: &Metadata,
    method: &Method,
    sink: &mut W,
) -> FitsResult<()> {
    let (tiling, opts, quantize) = match method {
        Method::Gzip { tiling, level } => (
            tiling,
            TileCompression::gzip().gzip_level(*level).lossless(),
            None,
        ),
        Method::Rice { tiling, quantize } => (tiling, TileCompression::rice(), Some(quantize)),
        Method::None | Method::Hcompress { .. } => unreachable!("written elsewhere"),
    };
    let grid = compress::resolve_grid(view, tiling, false)?;
    let mut opts = opts.tile_dims(&grid.ztile);
    if let Some(q) = quantize.filter(|_| view.is_float()) {
        let (delta, seed) = compress::global_quant(view, &grid, q);
        // One global step in every tile, as `HCOMPRESS_1` uses.
        opts = opts.quantize(
            fitsio_pure::compress::Quantize::new()
                .step(delta)
                .seed(fitsio_pure::compress::DitherSeed::Fixed(seed as u16)),
        );
    }
    let naxes = view.axes();
    let hdu = match view.planar() {
        Planar::U8(p) => compress_image_hdu(&naxes, &p, &opts, &[]),
        Planar::U16(p) => compress_image_hdu(&naxes, &p, &opts, &[]),
        Planar::F32(p) => compress_image_hdu(&naxes, &p, &opts, &[]),
    }
    .map_err(|e| FitsError::Compression(e.to_string()))?;

    // `fitsio-pure` writes `BZERO`/`BSCALE` for `u16` itself.
    let mut h = Header::new();
    if let Some(bits) = view.adc_bits {
        h.integer("BITADC", bits as i64, Some("meaningful ADC bits per sample"))?;
    }
    write_common_cards(&mut h, view, meta.timestamp())?;
    write_metadata(&mut h, meta)?;
    let ours = h.finish();

    let theirs = header_cards(&hdu);
    let data = &hdu[(theirs.len() + 80).div_ceil(card::BLOCK) * card::BLOCK..];
    let mut header = theirs.to_vec();
    header.extend_from_slice(header_cards(&ours));
    sink.write_all(&card::end_block(header))?;
    sink.write_all(data)?;
    Ok(())
}

/// The 80-byte cards of the header at the start of `bytes`, `END` excluded.
fn header_cards(bytes: &[u8]) -> &[u8] {
    let end = bytes
        .chunks(80)
        .position(|c| c.starts_with(b"END     "))
        .expect("a serialised header ends with END");
    &bytes[..end * 80]
}

/// An `HCOMPRESS_1` tile-compressed image HDU (`BINTABLE` with `ZIMAGE = T`),
/// streamed to `sink`.
fn write_hcompress<W: Write>(
    view: &ImageView<'_>,
    meta: &Metadata,
    c: &compress::Compressed,
    sink: &mut W,
) -> FitsResult<()> {
    let quantized = c.quant.is_some();
    let mut h = Header::new();

    h.string("XTENSION", "BINTABLE", Some("binary table extension"))?;
    h.integer("BITPIX", 8, None)?;
    h.integer("NAXIS", 2, None)?;
    h.integer(
        "NAXIS1",
        c.row_bytes as i64,
        Some("width of table row (bytes)"),
    )?;
    h.integer("NAXIS2", c.naxis2 as i64, Some("number of tiles"))?;
    h.integer("PCOUNT", c.pcount as i64, Some("heap size (bytes)"))?;
    h.integer("GCOUNT", 1, None)?;
    h.integer("TFIELDS", if quantized { 3 } else { 1 }, None)?;
    h.string("TTYPE1", "COMPRESSED_DATA", None)?;
    h.string("TFORM1", "1PB", Some("variable-length array of bytes"))?;
    if quantized {
        h.string("TTYPE2", "ZSCALE", None)?;
        h.string("TFORM2", "1D", None)?;
        h.string("TTYPE3", "ZZERO", None)?;
        h.string("TFORM3", "1D", None)?;
    }

    h.logical("ZIMAGE", true, Some("tile-compressed image"))?;
    h.string("ZCMPTYPE", "HCOMPRESS_1", Some("compression algorithm"))?;
    h.integer("ZBITPIX", c.zbitpix, Some("data type of original image"))?;
    h.integer("ZNAXIS", c.zaxes.len() as i64, None)?;
    for (i, &n) in c.zaxes.iter().enumerate() {
        h.integer(&format!("ZNAXIS{}", i + 1), n as i64, None)?;
    }
    for (i, &n) in c.ztiles.iter().enumerate() {
        h.integer(&format!("ZTILE{}", i + 1), n as i64, None)?;
    }
    h.string("ZNAME1", "SCALE", Some("HCOMPRESS scale factor"))?;
    h.real("ZVAL1", c.hscale as f64, Some("HCOMPRESS scale factor"))?;
    h.string("ZNAME2", "SMOOTH", Some("HCOMPRESS smooth option"))?;
    h.integer("ZVAL2", c.hsmooth as i64, Some("HCOMPRESS smooth option"))?;
    if let Some(zdither0) = c.quant {
        h.string(
            "ZQUANTIZ",
            "SUBTRACTIVE_DITHER_1",
            Some("lossy float quantization"),
        )?;
        h.integer("ZDITHER0", zdither0 as i64, Some("dithering offset"))?;
    }

    if let Some(bz) = view.bzero() {
        h.integer("BZERO", bz, Some("offset for unsigned integers"))?;
        h.integer("BSCALE", 1, None)?;
    }
    if let Some(bits) = view.adc_bits {
        h.integer("BITADC", bits as i64, Some("meaningful ADC bits per sample"))?;
    }

    write_common_cards(&mut h, view, meta.timestamp())?;
    write_metadata(&mut h, meta)?;

    sink.write_all(&h.finish())?;
    sink.write_all(&c.data)?;
    card::pad_writer(sink, c.data.len())?;
    Ok(())
}

/// `DATE-OBS`, `COLORSPC`, `BAYERPAT` — the image's own descriptive cards.
fn write_common_cards(h: &mut Header, view: &ImageView<'_>, ts: DateTime<Utc>) -> FitsResult<()> {
    h.string(
        "DATE-OBS",
        &datetime::to_iso8601(ts),
        Some("UTC of exposure start"),
    )?;
    h.string(
        "COLORSPC",
        &colorspace_str(&view.cspace),
        Some("refimage colour space"),
    )?;
    if let Some(p) = bayer_pattern(&view.cspace) {
        h.string("BAYERPAT", p, Some("Bayer mosaic pattern"))?;
    }
    Ok(())
}

/// Emit the typed exposure plus one or more cards for every extra metadata entry,
/// in insertion order.
fn write_metadata(h: &mut Header, meta: &Metadata) -> FitsResult<()> {
    let exp = meta.exposure();
    if !exp.is_zero() {
        h.integer(
            &format!("{EXPOSURE_KEY}_S"),
            exp.as_secs() as i64,
            Some("[s]"),
        )?;
        h.integer(
            &format!("{EXPOSURE_KEY}_NS"),
            exp.subsec_nanos() as i64,
            Some("[ns]"),
        )?;
        h.real(EXPOSURE_KEY, exp.as_secs_f64(), Some("[s] exposure time"))?;
    }
    if let Some(frame_id) = meta.frame_id() {
        h.integer(
            FRAMEID_KEY,
            i64::from(frame_id),
            Some("acquisition frame id"),
        )?;
    }
    for (key, item) in meta.iter() {
        if card::is_reserved(key) {
            return Err(FitsError::ReservedKeyword(key.clone()));
        }
        write_item(h, key, item)?;
    }
    Ok(())
}

fn write_item(h: &mut Header, key: &str, item: &GenericLineItem) -> FitsResult<()> {
    let comment = item.comment();
    match item.value() {
        GenericValue::Integer(v) => h.integer(key, *v, comment)?,
        GenericValue::Real(v) => h.real(key, *v, comment)?,
        GenericValue::String(s) => {
            if s.len() > 65_536 {
                return Err(FitsError::MetadataValueTooLong);
            }
            h.string(key, s, comment)?;
        }
        GenericValue::ColorSpace(cs) => h.string(key, &colorspace_str(cs), comment)?,
        GenericValue::Duration(d) => {
            h.integer(&format!("{key}_S"), d.as_secs() as i64, Some("[s]"))?;
            h.integer(&format!("{key}_NS"), d.subsec_nanos() as i64, Some("[ns]"))?;
            h.real(key, d.as_secs_f64(), comment.or(Some("[s]")))?;
        }
        GenericValue::Timestamp(t) => {
            // The image's own timestamp is emitted as `DATE-OBS` by
            // `write_common_cards`; this arm is for any *other* timestamp value.
            let (secs, nanos) = datetime::epoch_parts(*t);
            h.integer(&format!("{key}_S"), secs, Some("[s] since 1970-01-01"))?;
            h.integer(&format!("{key}_NS"), nanos as i64, Some("[ns]"))?;
            h.string(key, &datetime::to_iso8601(*t), comment)?;
        }
    }
    Ok(())
}

// Bridge so the macro can accept either a `&DynamicImageRef` or `&DynamicImageOwned`.
impl<'a> ImageView<'a> {
    fn from_ref_like<T: AsImageView>(v: &'a T) -> Self {
        v.as_image_view()
    }
}

trait AsImageView {
    fn as_image_view(&self) -> ImageView<'_>;
}

impl AsImageView for crate::DynamicImageRef<'_> {
    fn as_image_view(&self) -> ImageView<'_> {
        ImageView::from_ref(self)
    }
}

impl AsImageView for crate::DynamicImageOwned {
    fn as_image_view(&self) -> ImageView<'_> {
        ImageView::from_owned(self)
    }
}

#[cfg(test)]
mod tests;
