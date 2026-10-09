//! FITS header cards: 80-byte records, 2880-byte blocks.
//!
//! Everything here writes *fixed-format* cards where it can and falls back to the
//! `HIERARCH` (long keyword) and `CONTINUE` (long string) conventions — both of which
//! `astropy.io.fits` reads.

use super::{FitsError, FitsResult};

const CARD: usize = 80;
pub(super) const BLOCK: usize = 2880;

/// Structural keywords a caller must not shadow with metadata.
pub(super) const RESERVED: &[&str] = &[
    "SIMPLE", "BITPIX", "NAXIS", "EXTEND", "END", "XTENSION", "PCOUNT", "GCOUNT", "BSCALE",
    "BZERO", "TFIELDS", "THEAP", "GROUPS", "BLOCKED", "CONTINUE", "COMMENT", "HISTORY", "DATE-OBS",
    "COLORSPC", "BAYERPAT",
];

/// Is `key` a valid 8-character fixed-format FITS keyword?
fn is_fixed(key: &str) -> bool {
    !key.is_empty()
        && key.len() <= 8
        && key
            .bytes()
            .all(|b| b.is_ascii_uppercase() || b.is_ascii_digit() || b == b'-' || b == b'_')
}

/// Accumulates header cards and pads the finished header to a block boundary.
pub(super) struct Header {
    buf: Vec<u8>,
    longstrn: bool,
}

impl Header {
    pub(super) fn new() -> Self {
        Self {
            buf: Vec::with_capacity(BLOCK),
            longstrn: false,
        }
    }

    /// Append a verbatim card, space-padded to 80 columns.
    fn raw(&mut self, text: &str) {
        debug_assert!(text.len() <= CARD, "card over 80 cols: {text:?}");
        let start = self.buf.len();
        self.buf.extend_from_slice(text.as_bytes());
        self.buf.resize(start + CARD, b' ');
    }

    /// `KEYWORD = value` prefix for a fixed keyword (10 cols).
    fn fixed_prefix(key: &str) -> String {
        format!("{key:<8}= ")
    }

    pub(super) fn logical(
        &mut self,
        key: &str,
        val: bool,
        comment: Option<&str>,
    ) -> FitsResult<()> {
        self.valued(key, if val { "T" } else { "F" }, comment, true)
    }

    pub(super) fn integer(&mut self, key: &str, val: i64, comment: Option<&str>) -> FitsResult<()> {
        self.valued(key, &val.to_string(), comment, true)
    }

    pub(super) fn real(&mut self, key: &str, val: f64, comment: Option<&str>) -> FitsResult<()> {
        if !val.is_finite() {
            self.commentary("COMMENT", &format!("{key} omitted: non-finite value"));
            return Ok(());
        }
        // `{:?}` keeps a decimal point (so it never reads as an integer) and is
        // round-trip-shortest; uppercase the exponent for portability.
        let s = format!("{val:?}").replace('e', "E");
        self.valued(key, &s, comment, true)
    }

    /// Write a numeric/logical value: right-justified to col 30 for fixed keywords,
    /// free-format for HIERARCH.
    fn valued(
        &mut self,
        key: &str,
        val: &str,
        comment: Option<&str>,
        numeric: bool,
    ) -> FitsResult<()> {
        let mut card = if is_fixed(key) {
            let mut c = Self::fixed_prefix(key);
            if numeric && val.len() <= 20 {
                c.push_str(&format!("{val:>20}"));
            } else {
                c.push_str(val);
            }
            c
        } else {
            let c = format!("HIERARCH {key} = {val}");
            if c.len() > CARD {
                return Err(FitsError::KeywordTooLong(key.to_string()));
            }
            c
        };
        if let Some(cmt) = comment
            && card.len() + 3 < CARD {
                let room = CARD - card.len() - 3;
                card.push_str(" / ");
                card.push_str(&cmt[..cmt.len().min(room)]);
            }
        self.raw(&card);
        Ok(())
    }

    /// A `COMMENT` / `HISTORY` style card (keyword, then free text, no `=`).
    pub(super) fn commentary(&mut self, key: &str, text: &str) {
        for chunk in text.as_bytes().chunks(CARD - 8) {
            let line = std::str::from_utf8(chunk).unwrap_or("");
            self.raw(&format!("{key:<8}{line}"));
        }
    }

    /// A string-valued card, using `CONTINUE` if the value does not fit one card.
    pub(super) fn string(&mut self, key: &str, val: &str, comment: Option<&str>) -> FitsResult<()> {
        let escaped = val.replace('\'', "''");

        let prefix = if is_fixed(key) {
            Self::fixed_prefix(key)
        } else {
            let p = format!("HIERARCH {key} = ");
            if p.len() + 10 > CARD {
                return Err(FitsError::KeywordTooLong(key.to_string()));
            }
            p
        };

        // Columns available for quoted content on the first card (leave 2 for quotes).
        let first_room = CARD - prefix.len() - 2;

        if escaped.len() <= first_room {
            let body = pad8(&escaped);
            let mut card = format!("{prefix}'{body}'");
            if let Some(cmt) = comment
                && card.len() + 3 < CARD {
                    let room = CARD - card.len() - 3;
                    card.push_str(" / ");
                    card.push_str(&cmt[..cmt.len().min(room)]);
                }
            self.raw(&card);
            return Ok(());
        }

        // Long string: OGIP CONTINUE convention.
        if !self.longstrn {
            // Must appear before the first long string; header order is otherwise free.
            let card = format!("{:<8}= 'OGIP 1.0'", "LONGSTRN");
            self.raw(&card);
            self.longstrn = true;
        }

        let bytes = escaped.as_bytes();
        let mut pos = 0;
        // First segment.
        let take = first_room.saturating_sub(1).min(bytes.len());
        let seg = std::str::from_utf8(&bytes[..take]).unwrap_or("");
        self.raw(&format!("{prefix}'{seg}&'"));
        pos += take;

        while pos < bytes.len() {
            let room = CARD - 10 - 2; // "CONTINUE  " + quotes
            let last = pos + room >= bytes.len();
            let n = if last { bytes.len() - pos } else { room - 1 };
            let seg = std::str::from_utf8(&bytes[pos..pos + n]).unwrap_or("");
            pos += n;
            let tail = if last { "'" } else { "&'" };
            let mut card = format!("CONTINUE  '{seg}{tail}");
            if last
                && let Some(cmt) = comment
                    && card.len() + 3 < CARD {
                        let room = CARD - card.len() - 3;
                        card.push_str(" / ");
                        card.push_str(&cmt[..cmt.len().min(room)]);
                    }
            self.raw(&card);
        }
        Ok(())
    }

    /// Emit `END` and pad the header to a 2880-byte boundary with spaces.
    pub(super) fn finish(self) -> Vec<u8> {
        end_block(self.buf)
    }
}

/// Whole 80-byte `cards` followed by `END`, padded to a 2880-byte boundary with
/// spaces.
pub(super) fn end_block(mut cards: Vec<u8>) -> Vec<u8> {
    let start = cards.len();
    cards.extend_from_slice(b"END");
    cards.resize(start + CARD, b' ');
    cards.resize(cards.len().div_ceil(BLOCK) * BLOCK, b' ');
    cards
}

/// Pad a string to at least 8 characters with trailing spaces (FITS minimum string
/// value width).
fn pad8(s: &str) -> String {
    if s.len() >= 8 {
        s.to_string()
    } else {
        format!("{s:<8}")
    }
}

/// Write zero bytes to `w` to pad a data section of `data_len` bytes out to a
/// 2880-byte boundary.
pub(super) fn pad_writer<W: std::io::Write>(w: &mut W, data_len: usize) -> std::io::Result<()> {
    let rem = data_len % BLOCK;
    if rem != 0 {
        w.write_all(&vec![0u8; BLOCK - rem])?;
    }
    Ok(())
}

/// True if `key` (upper-cased) shadows a structural keyword.
pub(super) fn is_reserved(key: &str) -> bool {
    let k = key.to_ascii_uppercase();
    RESERVED.contains(&k.as_str())
        || k.starts_with("NAXIS")
        || k.starts_with('Z')
        || k.starts_with('T') && k[1..].chars().all(|c| c.is_ascii_digit()) && k.len() > 1
}
