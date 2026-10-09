//! Tile grids and float quantization for the compressed writers, and the
//! `HCOMPRESS_1` tile-compressed `BINTABLE` HDU (`ZIMAGE` convention).
//!
//! `GZIP_1` and `RICE_1` tiles are compressed by `fitsio_pure::compress`
//! (`mod.rs`); `HCOMPRESS_1` is built here. Both take their tile grid and their
//! float quantization from this module, so the three algorithms tile and quantize
//! alike.
//!
//! Tiling is resolved from the caller's [`Tiling`] choice and the concrete image
//! size: the default is one image row per tile for `GZIP_1` / `RICE_1` and one whole
//! channel plane for `HCOMPRESS_1`, but [`Rice::tile_rows`](super::Rice::tile_rows) /
//! [`Rice::tile_dims`](super::Rice::tile_dims) (and the same on the other builders) pick
//! any rectangular grid. Tiles are numbered fastest-FITS-axis-first; edge tiles are
//! clipped, not padded.
//!
//! `f32` + `Rice` / `Hcompress` is quantized to 32-bit integers with one global step
//! for the whole image ([`quantize::global_delta`]); `ZSCALE` is constant and `ZZERO`
//! is per-tile, both written as `1D` columns alongside `COMPRESSED_DATA`.

use super::config::{DitherSeed, Quantize, Tiling};
use super::hdu::ImageView;
use super::{hcompress, quantize, FitsError, FitsResult};

/// Everything `mod.rs` needs to emit an `HCOMPRESS_1` image HDU.
pub(super) struct Compressed {
    pub zbitpix: i64,
    pub zaxes: Vec<usize>,
    pub ztiles: Vec<usize>,
    /// HCOMPRESS scale factor (`ZVAL1`); `0` = lossless.
    pub hscale: i64,
    /// HCOMPRESS `SMOOTH` flag (`ZVAL2`).
    pub hsmooth: bool,
    /// `Some` when `f32` data was quantized: carries `ZDITHER0`.
    pub quant: Option<u32>,
    pub naxis2: usize,
    pub pcount: usize,
    /// `NAXIS1` — table row width in bytes (8, or 24 when quantized).
    pub row_bytes: usize,
    /// Main table (fixed columns + descriptors) followed by the heap; not block-padded.
    pub data: Vec<u8>,
}

/// The tile grid, resolved against the image size.
pub(super) struct Grid {
    tx: usize,
    ty: usize,
    ntx: usize,
    nty: usize,
    /// `ZTILE1..ZTILEn`.
    pub ztile: Vec<usize>,
}

pub(super) fn resolve_grid(
    view: &ImageView<'_>,
    tiling: &Tiling,
    hcompress: bool,
) -> FitsResult<Grid> {
    let axes = view.axes();
    let (mut tx, mut ty) = match tiling {
        Tiling::Default => {
            if hcompress {
                (view.w, view.h)
            } else {
                (view.w, 1)
            }
        }
        Tiling::Rows(n) => (view.w, (*n).max(1)),
        Tiling::Dims(d) => {
            if d.is_empty() || d.iter().take(2).any(|&v| v == 0) {
                return Err(FitsError::InvalidTiling(
                    "tile dimensions must be non-zero".into(),
                ));
            }
            if d.len() > axes.len() || d.iter().skip(2).any(|&v| v != 1) {
                return Err(FitsError::InvalidTiling(
                    "only the first two tile dimensions may exceed 1".into(),
                ));
            }
            (d[0], d.get(1).copied().unwrap_or(1))
        }
    };
    tx = tx.min(view.w).max(1);
    ty = ty.min(view.h).max(1);

    let ntx = view.w.div_ceil(tx);
    let nty = view.h.div_ceil(ty);

    if hcompress {
        let last_w = view.w - (ntx - 1) * tx;
        let last_h = view.h - (nty - 1) * ty;
        if tx < 4 || ty < 4 || last_w < 4 || last_h < 4 {
            return Err(FitsError::HcompressTooSmall);
        }
    }

    let mut ztile = vec![1usize; axes.len()];
    ztile[0] = tx;
    ztile[1] = ty;

    Ok(Grid {
        tx,
        ty,
        ntx,
        nty,
        ztile,
    })
}

/// The global quantization step and the dither seed (`ZDITHER0`) for an `f32` image.
pub(super) fn global_quant(view: &ImageView<'_>, grid: &Grid, q: &Quantize) -> (f64, u32) {
    let planar = view.planar_f32();
    let delta = quantize::global_delta(&planar, view.w, view.h * view.ch, q.level);
    let seed = match q.seed {
        DitherSeed::Auto => {
            let tw = grid.tx.min(view.w);
            let th = grid.ty.min(view.h);
            quantize::dither_seed(&view.rect_f32(0, 0, 0, tw, th))
        }
        DitherSeed::Fixed(n) => n.clamp(1, 10_000),
    };
    (delta, seed)
}

pub(super) fn hcompress(
    view: &ImageView<'_>,
    tiling: &Tiling,
    scale: i32,
    smooth: bool,
    q: &Quantize,
) -> FitsResult<Compressed> {
    let grid = resolve_grid(view, tiling, true)?;
    let quant = view.is_float().then(|| global_quant(view, &grid, q));

    let quantized = quant.is_some();
    let row_bytes = if quantized { 24 } else { 8 };
    let n_tiles = grid.ntx * grid.nty * view.ch;

    let mut heap = Vec::new();
    let mut table = Vec::with_capacity(n_tiles * row_bytes);

    let mut tile_index = 0usize;
    for plane in 0..view.ch {
        for iy in 0..grid.nty {
            for ix in 0..grid.ntx {
                let x0 = ix * grid.tx;
                let y0 = iy * grid.ty;
                let tw = grid.tx.min(view.w - x0);
                let th = grid.ty.min(view.h - y0);

                let (compressed, zzero) = if let Some((delta, seed)) = quant {
                    let f = view.rect_f32(plane, x0, y0, tw, th);
                    let mut q = quantize::quantize_tile(&f, delta, tile_index, seed);
                    let bytes = hcompress::compress(&mut q.idata, tw, th, scale);
                    (bytes, Some(q.zzero))
                } else {
                    let mut d = view.rect_i32(plane, x0, y0, tw, th);
                    (hcompress::compress(&mut d, tw, th, scale), None)
                };

                let offset = heap.len() as i32;
                let nbytes = compressed.len() as i32;
                table.extend_from_slice(&nbytes.to_be_bytes());
                table.extend_from_slice(&offset.to_be_bytes());
                if let Some((delta, _)) = quant {
                    table.extend_from_slice(&delta.to_be_bytes());
                    table.extend_from_slice(&zzero.unwrap_or(0.0).to_be_bytes());
                }
                heap.extend_from_slice(&compressed);
                tile_index += 1;
            }
        }
    }

    let pcount = heap.len();
    let mut data = table;
    data.extend_from_slice(&heap);

    Ok(Compressed {
        zbitpix: view.bitpix(),
        zaxes: view.axes(),
        ztiles: grid.ztile,
        hscale: scale as i64,
        hsmooth: smooth,
        quant: quant.map(|(_, seed)| seed),
        naxis2: n_tiles,
        pcount,
        row_bytes,
        data,
    })
}
