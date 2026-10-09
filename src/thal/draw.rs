//! ASCII structure output, matching primer3-py's `ascii_structure`.
//!
//! Rows are built from the uppercased input bytes, as thal.c does; the
//! result is decoded lossily only for non-UTF-8 input.

use std::iter::repeat_n;

/// Four-row duplex drawing (thal.c `drawDimer`). `ps1`/`ps2` give each
/// base's 1-based partner, 0 if unpaired; `oligo2` is already reversed.
pub(super) fn dimer(ps1: &[usize], ps2: &[usize], oligo1: &[u8], oligo2: &[u8]) -> String {
    let (len1, len2) = (oligo1.len(), oligo2.len());
    let mut rows: [Vec<u8>; 4] = Default::default();
    let [top, top_pair, bottom_pair, bottom] = &mut rows;

    // Unpaired 5' overhangs, right-aligned against each other.
    let ss1 = ps1.iter().take_while(|&&p| p == 0).count();
    let ss2 = ps2.iter().take_while(|&&p| p == 0).count();
    if ss1 >= ss2 {
        top.extend_from_slice(&oligo1[..ss1]);
        top_pair.extend(repeat_n(b' ', ss1));
        bottom_pair.extend(repeat_n(b' ', ss1));
        bottom.extend(repeat_n(b' ', ss1 - ss2));
        bottom.extend_from_slice(&oligo2[..ss2]);
    } else {
        bottom.extend_from_slice(&oligo2[..ss2]);
        top_pair.extend(repeat_n(b' ', ss2));
        bottom_pair.extend(repeat_n(b' ', ss2));
        top.extend(repeat_n(b' ', ss2 - ss1));
        top.extend_from_slice(&oligo1[..ss1]);
    }

    let (mut i, mut j) = (ss1, ss2); // 0-based
    while i < len1 {
        // Paired stretch.
        while i < len1 && j < len2 && ps1[i] != 0 && ps2[j] != 0 {
            top.push(b' ');
            top_pair.push(oligo1[i]);
            bottom_pair.push(oligo2[j]);
            bottom.push(b' ');
            i += 1;
            j += 1;
        }
        // Unpaired bases on each strand, padded with '-' to equal length.
        let start1 = i;
        while i < len1 && ps1[i] == 0 {
            top.push(oligo1[i]);
            top_pair.push(b' ');
            i += 1;
        }
        let start2 = j;
        while j < len2 && ps2[j] == 0 {
            bottom_pair.push(b' ');
            bottom.push(oligo2[j]);
            j += 1;
        }
        let (gap1, gap2) = (i - start1, j - start2);
        if gap1 < gap2 {
            top.extend(repeat_n(b'-', gap2 - gap1));
            top_pair.extend(repeat_n(b' ', gap2 - gap1));
        } else if gap1 > gap2 {
            bottom_pair.extend(repeat_n(b' ', gap1 - gap2));
            bottom.extend(repeat_n(b'-', gap1 - gap2));
        }
    }

    let [top, top_pair, bottom_pair, bottom] = rows.map(|mut row| {
        while row.last().is_some_and(u8::is_ascii_whitespace) {
            row.pop();
        }
        String::from_utf8_lossy(&row).into_owned()
    });
    format!("SEQ\t{top}\nSEQ\t{top_pair}\nSTR\t{bottom_pair}\nSTR\t{bottom}\n")
}

/// Two-row hairpin drawing (thal.c `drawHairpin`): '/' and '\' mark the
/// 5' and 3' sides of each pair.
pub(super) fn hairpin(bp: &[usize], oligo: &[u8]) -> String {
    let mut row = vec![b'-'; bp.len()];
    for (i, &partner) in bp.iter().enumerate() {
        // As in thal.c, the mark goes at the partner's position.
        if partner > 0 {
            row[partner - 1] = if partner > i { b'\\' } else { b'/' };
        }
    }
    format!(
        "SEQ\t{}\nSTR\t{}\n",
        String::from_utf8_lossy(&row),
        String::from_utf8_lossy(oligo)
    )
}
