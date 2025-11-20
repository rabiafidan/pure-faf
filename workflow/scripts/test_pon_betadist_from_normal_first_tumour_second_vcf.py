#!/usr/bin/env python3
"""
Add p_ADfit to a gzipped VCF using a Beta-Binomial fit.

For each variant:
- Read BETA=alpha,beta from INFO (key configurable; default "BETA").
- Read AD (allelic depths as "ref,alt") and DP (depth) from the tumor sample.
- Compute p_ADfit = 1 - CDF(AD_alt | DP, alpha, beta) using scipy.stats.betabinom.
- Append p_ADfit to INFO.

Header:
- Injects the INFO definition for p_ADfit exactly once (unless it already exists).

Assumptions / Notes:
- Input VCF is gzipped (.vcf.gz) and bgzipped is fine; script reads via gzip.
- Output is plain-text .vcf (you can bgzip + index afterwards with bcftools).
- Tumor sample is assumed to be the last sample column unless --tumor-name is given.
- If DP is missing but AD exists, falls back to DP = sum(AD parts).
- Lines missing required fields are passed through unchanged (no hard failure).
"""

import sys
import io
import gzip
import argparse
from typing import Dict, Optional, Tuple

# SciPy provides the Beta-Binomial distribution
from scipy.stats import betabinom


NEW_INFO_LINE = (
    '##INFO=<ID=p_ADfit,Number=1,Type=Float,'
    'Description="Probability that the alt AF fits the PON beta-binomial distribution '
    '(1 - CDF(AD_alt | DP, alpha, beta))">'
)

def parse_args():
    """Parse CLI arguments."""
    p = argparse.ArgumentParser(
        description="Append p_ADfit to VCF INFO using Beta-Binomial parameters from INFO:BETA."
    )
    p.add_argument("vcf_gz", help="Input gzipped VCF (.vcf.gz)")
    p.add_argument("out_vcf", help="Output plain-text VCF (.vcf)")
    p.add_argument(
        "--beta-key", default="BETA",
        help="INFO key carrying 'alpha,beta' (default: BETA)"
    )
    p.add_argument(
        "--tumor-name", default=None,
        help="Exact sample name to treat as tumor (overrides last-column default)."
    )
    p.add_argument(
        "--precision", type=int, default=6,
        help="Decimal digits for p_ADfit formatting (default: 6)."
    )
    return p.parse_args()


def parse_info(info_str: str) -> Dict[str, str]:
    """
    Parse INFO field into a dict.
    - 'A=B;C;D=E' -> {'A': 'B', 'C': True, 'D': 'E'}
    """
    d: Dict[str, str] = {}
    if info_str == "." or info_str == "":
        return d
    for kv in info_str.split(";"):
        if not kv:
            continue
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
        else:
            d[kv] = True  # flag without a value
    return d


def ensure_info_declared(header_seen: bool, fout, inserted_flag: bool) -> bool:
    """
    Ensure our NEW_INFO_LINE is present once, right before #CHROM.
    If header already has it, do nothing.
    Return True if we inserted it in this call; otherwise return the prior flag.
    """
    if not header_seen and not inserted_flag:
        fout.write(NEW_INFO_LINE + "\n")
        return True
    return inserted_flag


def locate_tumor_index(samples: list, tumor_name: Optional[str]) -> int:
    """
    Decide which sample column is tumor.
    - If tumor_name provided, find exact match; else use last sample.
    Raise ValueError if tumor_name given but not found.
    """
    if tumor_name:
        try:
            return samples.index(tumor_name)
        except ValueError:
            raise ValueError(f"Tumor sample name '{tumor_name}' not found in header samples: {samples}")
    # default: last sample
    return len(samples) - 1


def extract_fmt_values(format_str: str, sample_str: str) -> Dict[str, str]:
    """
    Map FORMAT keys to the corresponding sample values.
    Example:
      format_str = "GT:AD:DP:AF"
      sample_str = "0/1:12,8:20:0.4"
      -> {"GT":"0/1","AD":"12,8","DP":"20","AF":"0.4"}
    """
    keys = format_str.split(":")
    vals = sample_str.split(":")
    m = {}
    # keys and vals can differ in length; guard with min
    for i in range(min(len(keys), len(vals))):
        m[keys[i]] = vals[i]
    return m


def get_alpha_beta(info_dict: Dict[str, str], beta_key: str) -> Optional[Tuple[float, float]]:
    """
    From INFO dict, extract alpha,beta from key beta_key (e.g., 'BETA' -> '3.1,2.9').
    Return (alpha, beta) as floats or None if missing/malformed.
    """
    tag = info_dict.get(beta_key)
    if not tag:
        return None
    try:
        a_str, b_str = tag.split(",", 1)
        return float(a_str), float(b_str)
    except Exception:
        return None


def compute_p_adfit(AD: Optional[str], DP: Optional[str], alpha: float, beta: float) -> Optional[float]:
    """
    Compute p_ADfit = 1 - CDF(AD_alt | DP, alpha, beta) using Beta-Binomial.
    - AD is expected as "ref,alt".
    - DP should be an integer; if missing, we try sum(AD parts).
    Returns None on any issue.
    """
    if AD is None and DP is None:
        return None

    try:
        ad_ref_alt = AD.split(",") if AD is not None else []
        alt = int(ad_ref_alt[1]) if len(ad_ref_alt) >= 2 else None

        # Determine DP: prefer explicit DP; fallback to sum of AD parts if sensible
        if DP is not None and DP != "." and DP != "":
            depth = int(DP)
        else:
            if AD is None:
                return None
            depth = sum(int(x) for x in ad_ref_alt if x not in (".", ""))

        if alt is None:
            return None
        if depth < 0 or alt < 0:
            return None
        if depth < alt:
            # Depth less than alt count is suspicious; still compute but guard input
            # Here we fail safe: return None to avoid nonsensical stats
            return None

        cdf = betabinom.cdf(alt, depth, alpha, beta)
        p = 1.0 - float(cdf)
        return p
    except Exception:
        return None


def main():
    args = parse_args()

    with io.TextIOWrapper(gzip.open(args.vcf_gz, "rb")) as fin, open(args.out_vcf, "w") as fout:
        inserted_info = False     # Did we inject our INFO meta line already?
        have_info_already = False # Did the header already declare p_ADfit?
        header_seen = False       # Have we seen the #CHROM header line yet?

        tumor_idx = None          # Index (within samples) of the tumor sample
        samples = []              # Sample names from header

        for raw in fin:
            # -----------------------------
            # METADATA LINES (start with ##)
            # -----------------------------
            if raw.startswith("##"):
                if raw.startswith("##INFO=<ID=p_ADfit,"):
                    have_info_already = True
                fout.write(raw)
                continue

            # -----------------------------
            # HEADER LINE (#CHROM ...)
            # -----------------------------
            if raw.startswith("#"):
                # Inject our INFO definition right before #CHROM, unless it's already present
                if not have_info_already and not inserted_info:
                    fout.write(NEW_INFO_LINE + "\n")
                    inserted_info = True

                # Write the #CHROM line
                fout.write(raw)
                header_seen = True

                # Parse sample names from #CHROM header (columns 9+)
                hdr = raw.rstrip("\n").split("\t")
                if len(hdr) > 9:
                    samples = hdr[9:]
                    # Decide tumor index now (if we have samples)
                    try:
                        tumor_idx = locate_tumor_index(samples, args.tumor_name)
                    except ValueError as e:
                        # Fail fast with a clear message (bad tumor name)
                        sys.exit(str(e))
                continue

            # -----------------------------
            # BODY LINES (actual variants)
            # -----------------------------
            fields = raw.rstrip("\n").split("\t")

            # If malformed (<10 columns), pass through unchanged
            if len(fields) < 10:
                fout.write("\t".join(fields) + "\n")
                continue

            # INFO as dict
            info_dict = parse_info(fields[7])

            # Extract (alpha, beta) from INFO
            ab = get_alpha_beta(info_dict, args.beta_key)
            if ab is None:
                # No usable BETA present -> write unchanged
                fout.write("\t".join(fields) + "\n")
                continue
            alpha, beta = ab

            # FORMAT → sample map; tumor sample selection
            format_str = fields[8]
            sample_cols = fields[9:]
            if not sample_cols:
                fout.write("\t".join(fields) + "\n")
                continue

            # Default to last sample if we didn't compute tumor_idx (e.g. weird header)
            t_idx = tumor_idx if tumor_idx is not None else (len(sample_cols) - 1)
            tumor_sample_str = sample_cols[t_idx]

            fmt_map = extract_fmt_values(format_str, tumor_sample_str)
            AD = fmt_map.get("AD")
            DP = fmt_map.get("DP")

            # Compute p_ADfit
            pval = compute_p_adfit(AD, DP, alpha, beta)
            if pval is None:
                # If computation not possible, pass record unchanged
                fout.write("\t".join(fields) + "\n")
                continue

            # Append to INFO (avoid leading ';' when INFO is '.' or empty)
            if fields[7] in (".", ""):
                fields[7] = f"p_ADfit={pval:.{args.precision}g}"
            else:
                fields[7] = fields[7] + f";p_ADfit={pval:.{args.precision}g}"

            # Write the complete line with a single newline AT THE END
            fout.write("\t".join(fields) + "\n")


if __name__ == "__main__":
    main()
