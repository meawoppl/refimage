//! Unit tests for the FITS writer. Structural checks only — `astropy` round-trips live
//! in `tests/fits_astropy.rs`, cfitsio cross-checks in `tests/fits_roundtrip.rs`.

use std::time::Duration;

use chrono::{DateTime, Utc};

use super::*;
use crate::{ColorSpace, DynamicImageOwned, GenericImageOwned, ImageOwned};

/// A fixed, deterministic UTC timestamp for tests.
fn ts() -> DateTime<Utc> {
    DateTime::from_timestamp(1_700_000_000, 0).unwrap()
}

fn gray_u16(w: usize, h: usize) -> GenericImageOwned {
    let data: Vec<u16> = (0..w * h).map(|i| (i as u16).wrapping_mul(97)).collect();
    let img =
        DynamicImageOwned::from(ImageOwned::from_owned(data, w, h, ColorSpace::Gray).unwrap());
    let mut g = GenericImageOwned::new(ts(), Duration::from_millis(1500), img);
    g.insert_key("CAMERA", "Test Cam").unwrap();
    g.insert_key("GAIN", 3u16).unwrap();
    g
}

fn is_block_aligned(b: &[u8]) -> bool {
    b.len().is_multiple_of(2880)
}

fn find_card<'a>(bytes: &'a [u8], key: &str) -> Option<&'a str> {
    // Search only header blocks (stop at first block with no '=' cards is overkill;
    // just scan everything in 80-byte steps until END).
    for card in bytes.chunks(80) {
        if card.len() < 80 {
            break;
        }
        let Ok(text) = std::str::from_utf8(card) else {
            continue;
        };
        if text.starts_with(&format!("{key:<8}=")) || text.starts_with(&format!("HIERARCH {key} ="))
        {
            return Some(text);
        }
    }
    None
}

#[test]
fn uncompressed_primary_structure() {
    let g = gray_u16(8, 4);
    let bytes = g.fits_bytes(FitsCompression::NONE).unwrap();
    assert!(is_block_aligned(&bytes));
    assert!(bytes.starts_with(b"SIMPLE  =                    T"));
    assert!(find_card(&bytes, "BITPIX").unwrap().contains("16"));
    assert!(find_card(&bytes, "NAXIS1").unwrap().contains("8"));
    assert!(find_card(&bytes, "NAXIS2").unwrap().contains("4"));
    assert!(find_card(&bytes, "BZERO").unwrap().contains("32768"));
    assert!(find_card(&bytes, "DATE-OBS").is_some());
    assert!(find_card(&bytes, "COLORSPC").unwrap().contains("GRAY"));
    // Duration split into _S / _NS plus a float base card.
    assert!(find_card(&bytes, "EXPOSURE_S").unwrap().contains("1"));
    assert!(find_card(&bytes, "EXPOSURE_NS")
        .unwrap()
        .contains("500000000"));
    // header + 8*4*2 = 64 data bytes -> 2 blocks
    assert_eq!(bytes.len(), 2880 * 2);
}

#[test]
fn frame_id_card() {
    // An unset frame ID is not written.
    let g = gray_u16(8, 4);
    assert_eq!(g.frame_id(), None);
    let bytes = g.fits_bytes(FitsCompression::NONE).unwrap();
    assert!(find_card(&bytes, "FRAMEID").is_none());

    // A set frame ID is written as an integer card.
    let mut g = gray_u16(8, 4);
    g.set_frame_id(12345);
    let bytes = g.fits_bytes(FitsCompression::NONE).unwrap();
    assert!(find_card(&bytes, "FRAMEID").unwrap().contains("12345"));
}

#[test]
fn compressed_is_bintable_extension() {
    let g = gray_u16(16, 16);
    for (comp, want) in [
        (FitsCompression::from(Gzip::new()), "GZIP_1"),
        (FitsCompression::from(Rice::new()), "RICE_1"),
    ] {
        let bytes = g.fits_bytes(&comp).unwrap();
        assert!(is_block_aligned(&bytes), "{comp:?}");
        // primary first
        assert!(bytes.starts_with(b"SIMPLE  =                    T"));
        assert!(find_card(&bytes, "XTENSION").unwrap().contains("BINTABLE"));
        assert!(find_card(&bytes, "ZIMAGE").unwrap().contains('T'));
        let zc = find_card(&bytes, "ZCMPTYPE").unwrap();
        assert!(zc.contains(want));
        assert!(find_card(&bytes, "ZNAXIS1").unwrap().contains("16"));
    }
}

#[test]
fn gzip_level_is_wired_through() {
    // A low-entropy plane so DEFLATE effort visibly matters.
    let data = vec![4321u16; 64 * 64];
    let img =
        DynamicImageOwned::from(ImageOwned::from_owned(data, 64, 64, ColorSpace::Gray).unwrap());
    let g = GenericImageOwned::new(ts(), Duration::from_millis(1), img);

    // One tile, so the effort has the whole plane to work with.
    let stored = g
        .fits_bytes(Gzip::new().level(0).tile_dims([64, 64]))
        .unwrap();
    let deflated = g.fits_bytes(Gzip::new().tile_dims([64, 64])).unwrap();
    // Effort 0 stores the bytes verbatim; the default effort (6) must beat it.
    assert!(deflated.len() < stored.len());
    // `level` is still callable after the tile shape is fixed (type-state).
    let _ = g.fits_bytes(Gzip::new().tile_rows(16).level(9)).unwrap();
}

#[test]
fn hcompress_is_bintable_with_scale_cards() {
    let g = gray_u16(20, 16);
    let bytes = g.fits_bytes(Hcompress::new()).unwrap();
    assert!(is_block_aligned(&bytes));
    assert!(bytes.starts_with(b"SIMPLE  =                    T"));
    assert!(find_card(&bytes, "ZCMPTYPE")
        .unwrap()
        .contains("HCOMPRESS_1"));
    assert!(find_card(&bytes, "ZNAME1").unwrap().contains("SCALE"));
    assert!(find_card(&bytes, "ZNAME2").unwrap().contains("SMOOTH"));
    // Whole image is one tile: ZTILE1 == ZNAXIS1.
    assert!(find_card(&bytes, "ZTILE1").unwrap().contains("20"));
    assert!(find_card(&bytes, "ZTILE2").unwrap().contains("16"));
}

#[test]
fn hcompress_rejects_tiny_images() {
    let data = vec![0u16; 9];
    let img =
        DynamicImageOwned::from(ImageOwned::from_owned(data, 3, 3, ColorSpace::Gray).unwrap());
    let g = GenericImageOwned::new(ts(), Duration::ZERO, img);
    assert!(matches!(
        g.fits_bytes(Hcompress::new()),
        Err(FitsError::HcompressTooSmall)
    ));
}

#[test]
fn float_compresses_both_ways() {
    let data: Vec<f32> = (0..64 * 8).map(|i| (i as f32 * 0.1).sin()).collect();
    let img =
        DynamicImageOwned::from(ImageOwned::from_owned(data, 64, 8, ColorSpace::Gray).unwrap());
    let g = GenericImageOwned::new(ts(), Duration::ZERO, img);

    // Gzip: lossless raw float bytes, marked `ZQUANTIZ = 'NONE'` as cfitsio marks them.
    let gz = g.fits_bytes(Gzip::new()).unwrap();
    assert!(find_card(&gz, "ZQUANTIZ").unwrap().contains("'NONE "));
    assert!(find_card(&gz, "ZDITHER0").is_none());

    // Rice: quantized, so ZQUANTIZ / ZDITHER0 / ZSCALE column appear.
    let rc = g.fits_bytes(Rice::new()).unwrap();
    assert!(find_card(&rc, "ZQUANTIZ")
        .unwrap()
        .contains("SUBTRACTIVE_DITHER_1"));
    assert!(find_card(&rc, "ZDITHER0").is_some());
    assert!(find_card(&rc, "TTYPE2").unwrap().contains("ZSCALE"));
    assert!(find_card(&rc, "ZBITPIX").unwrap().contains("-32"));
}

#[test]
fn rgb_is_planar_cube() {
    let data: Vec<u8> = (0..3 * 4 * 5).map(|i| i as u8).collect();
    let img = DynamicImageOwned::from(ImageOwned::from_owned(data, 5, 4, ColorSpace::Rgb).unwrap());
    let g = GenericImageOwned::new(ts(), Duration::ZERO, img);
    let bytes = g.fits_bytes(FitsCompression::NONE).unwrap();
    assert!(find_card(&bytes, "NAXIS").unwrap().contains('3'));
    assert!(find_card(&bytes, "NAXIS3").unwrap().contains('3'));
}

#[test]
fn reserved_metadata_key_errors() {
    let data = vec![0u8; 4];
    let img =
        DynamicImageOwned::from(ImageOwned::from_owned(data, 2, 2, ColorSpace::Gray).unwrap());
    let mut g = GenericImageOwned::new(ts(), Duration::ZERO, img);
    g.insert_key("NAXIS1", 5u16).unwrap();
    assert!(matches!(
        g.fits_bytes(FitsCompression::NONE),
        Err(FitsError::ReservedKeyword(_))
    ));
}

#[test]
fn multi_hdu_file() {
    let dir = std::env::temp_dir();
    let path = dir.join(format!("refimage_fits_test_{}.fits", std::process::id()));
    let g = gray_u16(8, 8);
    {
        let mut w = create_fits(&path, Rice::new(), true).unwrap();
        g.append_fits(&mut w).unwrap();
        g.append_fits(&mut w).unwrap();
        w.finish().unwrap();
    }
    let bytes = std::fs::read(&path).unwrap();
    assert!(is_block_aligned(&bytes));
    let n_xtension = bytes
        .chunks(80)
        .filter(|c| c.starts_with(b"XTENSION= 'BINTABLE'"))
        .count();
    assert_eq!(n_xtension, 2);
    std::fs::remove_file(&path).ok();
}

#[test]
fn in_memory_multi_hdu_file() {
    let g = gray_u16(8, 8);

    // Owned Vec sink, recovered from finish().
    let mut w = create_fits_to(Vec::new(), Rice::new()).unwrap();
    g.append_fits(&mut w).unwrap();
    g.append_fits(&mut w).unwrap();
    assert_eq!(w.hdu_count(), 3);
    let bytes = w.finish().unwrap();

    assert!(is_block_aligned(&bytes));
    assert!(bytes.starts_with(b"SIMPLE  =                    T"));
    let n_ext = bytes
        .chunks(80)
        .filter(|c| c.starts_with(b"XTENSION= 'BINTABLE'"))
        .count();
    assert_eq!(n_ext, 2);

    // Borrowed Vec sink: caller keeps ownership, same bytes.
    let mut buf = Vec::new();
    {
        let mut w = create_fits_to(&mut buf, Rice::new()).unwrap();
        g.append_fits(&mut w).unwrap();
        g.append_fits(&mut w).unwrap();
        w.finish().unwrap();
    }
    assert_eq!(buf, bytes);
}

#[test]
fn tile_rows_sets_ztile2_and_tile_count() {
    let g = gray_u16(16, 12);
    let bytes = g.fits_bytes(Rice::new().tile_rows(4)).unwrap();
    assert!(is_block_aligned(&bytes));
    assert!(find_card(&bytes, "ZTILE1").unwrap().contains("16"));
    assert!(find_card(&bytes, "ZTILE2").unwrap().contains('4'));
    // 12 rows / 4 per tile = 3 tiles.
    assert!(find_card(&bytes, "NAXIS2").unwrap().contains('3'));
}

#[test]
fn tile_dims_makes_a_rectangular_grid() {
    let g = gray_u16(20, 20);
    let bytes = g.fits_bytes(Gzip::new().tile_dims([8, 8])).unwrap();
    assert!(is_block_aligned(&bytes));
    assert!(find_card(&bytes, "ZTILE1").unwrap().contains('8'));
    assert!(find_card(&bytes, "ZTILE2").unwrap().contains('8'));
    // ceil(20/8) == 3 in each direction -> 9 tiles.
    assert!(find_card(&bytes, "NAXIS2").unwrap().contains('9'));
}

#[test]
fn quantize_level_is_honoured_for_floats() {
    // Needs real noise: `global_delta` only scales by the level when it finds any.
    let mut rng = 0x9e37_79b9u32;
    let data: Vec<f32> = (0..32 * 32)
        .map(|i| {
            rng = rng.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            let noise = (rng >> 27) as f32 - 16.0;
            (i as f32 * 0.05).sin() * 500.0 + noise
        })
        .collect();
    let img =
        DynamicImageOwned::from(ImageOwned::from_owned(data, 32, 32, ColorSpace::Gray).unwrap());
    let g = GenericImageOwned::new(ts(), Duration::ZERO, img);

    let pcount = |q: Quantize| -> i64 {
        let b = g.fits_bytes(Rice::new().quantize(q)).unwrap();
        find_card(&b, "PCOUNT")
            .unwrap()
            .split('=')
            .nth(1)
            .unwrap()
            .split('/')
            .next()
            .unwrap()
            .trim()
            .parse()
            .unwrap()
    };
    // A finer step keeps more bits, so the compressed heap is larger.
    let coarse = pcount(Quantize::new().level(0.5));
    let fine = pcount(Quantize::new().level(256.0));
    assert!(fine > coarse, "{fine} !> {coarse}");
}

#[test]
fn fixed_dither_seed_is_reproducible() {
    let data: Vec<f32> = (0..16 * 16)
        .map(|i| (i as f32 * 0.1).cos() * 20.0)
        .collect();
    let img =
        DynamicImageOwned::from(ImageOwned::from_owned(data, 16, 16, ColorSpace::Gray).unwrap());
    let g = GenericImageOwned::new(ts(), Duration::ZERO, img);
    let a = g
        .fits_bytes(Rice::new().quantize(Quantize::new().seed(DitherSeed::Fixed(1234))))
        .unwrap();
    let b = g
        .fits_bytes(Rice::new().quantize(Quantize::new().seed(DitherSeed::Fixed(1234))))
        .unwrap();
    assert_eq!(a, b);
    assert!(find_card(&a, "ZDITHER0").unwrap().contains("1234"));
}

#[test]
fn hcompress_rejects_tiny_tiles() {
    let g = gray_u16(16, 16);
    assert!(matches!(
        g.fits_bytes(Hcompress::new().tile_dims([2, 2])),
        Err(FitsError::HcompressTooSmall)
    ));
}

#[test]
fn bincode_roundtrip_still_works() {
    // Part A regression: numeric-collapsed GenericValue serialises.
    let g = gray_u16(4, 4);
    let ser = bincode::serialize(&g).unwrap();
    let de: GenericImageOwned = bincode::deserialize(&ser).unwrap();
    assert_eq!(g.metadata(), de.metadata());
}
