"""Read Table 1 of FRDM2012 (thesis/FRLDM_2012.pdf) into data/frdm2012_table.dat.

Usage:
    pdftotext -layout thesis/FRLDM_2012.pdf frdm2012.txt
    python3 tools/parse_frdm2012_table.py frdm2012.txt data/frdm2012_table.dat

Each row of Table 1 has 18 numbers when the experimental mass is known and 16 when
it is not (M_exp and sigma_exp missing):
    N A eps2 eps3 eps4 eps6 beta2 beta3 beta4 beta6 E_s+p E_mic E_bind M_th [M_exp sigma_exp] E_mic^FL M_th^FL
The proton number comes from the "Z = .." line above the rows. The output has one
nucleus per line; a missing M_exp, sigma_exp is written as nan.

【日本語】FRDM2012 の Table 1 を pdftotext で文字にしたものから，各核の
Z, N, A, eps2, eps3, eps4, eps6, beta2..beta6, E_s+p, E_mic, E_bind, M_th, M_exp, sigma_exp,
E_mic^FL, M_th^FL を読み取り，1 行 1 核のデータファイルにする．
実験値のない核は M_exp, sigma_exp を nan にする．9318 核すべてが読めることを確認している．
"""
import re
import sys

COLUMNS = ["Z", "N", "A", "eps2", "eps3", "eps4", "eps6", "beta2", "beta3", "beta4", "beta6",
           "E_sp", "E_mic", "E_bind", "M_th", "M_exp", "sigma_exp", "E_mic_FL", "M_th_FL"]


def parse(lines):
    rows = []
    Z = None
    for line in lines:
        m = re.match(r"\s*Z\s*=\s*(\d+)\s*(\(|$)", line)
        if m:
            Z = int(m.group(1))
            continue
        if Z is None:
            continue
        try:
            v = [float(x) for x in line.replace("−", "-").split()]
        except ValueError:
            continue
        if len(v) not in (16, 18):
            continue
        N, A = int(v[0]), int(v[1])
        if N + Z != A:
            continue
        if len(v) == 16:
            v = v[:14] + [float("nan"), float("nan")] + v[14:]
        rows.append([Z, N, A] + v[2:])
    return rows


def main():
    text, out = sys.argv[1], sys.argv[2]
    with open(text, encoding="utf-8") as f:
        rows = parse(f)
    keys = {(r[0], r[1]) for r in rows}
    if len(keys) != len(rows):
        sys.exit("duplicate nuclei in the table")
    with open(out, "w") as f:
        f.write("# FRDM2012 Table 1 (P. Moller et al., ADNDT 109-110 (2016) 1), %d nuclei\n" % len(rows))
        f.write("# energies in MeV; M_exp, sigma_exp from AME2003 (nan: not measured)\n")
        f.write("# " + " ".join(COLUMNS) + "\n")
        for r in rows:
            f.write("%4d %4d %4d" % tuple(r[:3])
                    + "".join(" %6.2f" % x for x in r[3:7])
                    + "".join(" %7.3f" % x for x in r[7:11])
                    + "".join(" %9.2f" % x for x in r[11:15])
                    + " %9.2f %6.3f" % (r[15], r[16])
                    + " %7.2f %9.2f\n" % (r[17], r[18]))
    print("%d nuclei written to %s" % (len(rows), out))


if __name__ == "__main__":
    main()
