"""Compute M_th with main.exe at the FRDM2012 ground-state shapes (eps2, eps3, eps4, eps6).

Usage (from the repository root, after bash build.sh):
    python3 tools/run_frdm2012_shapes.py --zmax 20 --measured --jobs 4
    python3 tools/run_frdm2012_shapes.py --zmin 82 --zmax 82 --no-zpe

For each nucleus of data/frdm2012_table.dat in the range, an input file with the table's
shape (gamma = 0) is written to the work directory, main.exe is run there, and E_mac,
E_s+p, E_zp and M_th are read from its output. The result file has one line per nucleus,
next to the table's M_th^FL (FRLDM2012), M_th (FRDM2012) and M_exp (AME2003). Nuclei
already in the result file are skipped, so an interrupted run can be continued.

【日本語】FRDM2012 の表の基底状態の形 (eps2, eps3, eps4, eps6, gamma = 0) で main.exe を
走らせ，E_mac, E_s+p, E_zp, M_th を表の M_th^FL (FRLDM2012), M_th (FRDM2012), M_exp と
並べて 1 行 1 核で書き出す．--jobs 個の核を同時に計算する (各 --threads スレッド)．
結果ファイルにすでにある核は飛ばすので，途中で止めても続きから計算できる．
"""
import argparse
import math
import os
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

HEADER = ("# Z N A eps2 eps3 eps4 eps6 E_mac E_sp E_zp M_th M_th_FL(table) M_th(FRDM2012) M_exp "
          "M_th-M_th_FL M_exp-M_th flags\n")


def read_table(path):
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            v = line.split()
            rows.append(dict(Z=int(v[0]), N=int(v[1]), A=int(v[2]), eps2=float(v[3]), eps3=float(v[4]),
                             eps4=float(v[5]), eps6=float(v[6]), M_th=float(v[14]), M_exp=float(v[15]),
                             M_th_FL=float(v[18])))
    return rows


def last(pattern, text):
    m = re.findall(pattern + r"\s*:\s*(-?[0-9.]+)", text)
    return float(m[-1]) if m else math.nan


def first(pattern, text):
    m = re.findall(pattern + r"\s*:\s*(-?[0-9.]+)", text)
    return float(m[0]) if m else math.nan


def run_one(r, args):
    name = "%03d_%03d" % (r["Z"], r["N"])
    nml = os.path.join(args.workdir, name + ".nml")
    log = os.path.join(args.workdir, name + ".log")
    with open(nml, "w") as f:
        f.write("&frldm_input\n")
        f.write("  Z = %d, N = %d,\n" % (r["Z"], r["N"]))
        f.write("  eps2 = %.2f, gamma_deg = 0.0, eps3 = %.2f, eps4 = %.2f, eps6 = %.2f,\n"
                % (r["eps2"], r["eps3"], r["eps4"], r["eps6"]))
        f.write("  pairing = '%s', n_gh = %d, zpe = %s\n/\n" % (args.pairing, args.n_gh,
                                                             ".true." if args.zpe else ".false."))
    env = dict(os.environ, OMP_NUM_THREADS=str(args.threads))
    p = subprocess.run([args.exe, name + ".nml"], cwd=args.workdir, env=env,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    out = p.stdout
    with open(log, "w") as f:
        f.write(out)
    flags = []
    if p.returncode != 0:
        flags.append("failed")
    if "no Lipkin-Nogami solution" in out:
        flags.append("LN->BCS")
    if "C_eps2 <= 0" in out:
        flags.append("C_eps2<=0")
    if "C_gamma <= 0" in out:
        flags.append("C_gamma<=0")
    e_zp = last(r"E_zp", out) if args.zpe else 0.0
    return dict(r, E_mac=first(r"E_mac", out), E_sp=first(r"E_s\+p", out), E_zp=e_zp,
                M_calc=last(r"mass excess M_th", out), flags=",".join(flags) or "-")


def line(r):
    return ("%4d %4d %4d %6.2f %6.2f %6.2f %6.2f %10.4f %9.4f %8.4f %10.4f %9.2f %9.2f %9.2f %8.3f %8.3f %s\n"
            % (r["Z"], r["N"], r["A"], r["eps2"], r["eps3"], r["eps4"], r["eps6"], r["E_mac"], r["E_sp"],
               r["E_zp"], r["M_calc"], r["M_th_FL"], r["M_th"], r["M_exp"], r["M_calc"] - r["M_th_FL"],
               r["M_exp"] - r["M_calc"], r["flags"]))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--table", default="data/frdm2012_table.dat")
    ap.add_argument("--exe", default="./main.exe")
    ap.add_argument("--out", default="data/frldm_at_frdm2012_shapes.dat")
    ap.add_argument("--workdir", default="runs_frdm2012")
    ap.add_argument("--zmin", type=int, default=0)
    ap.add_argument("--zmax", type=int, default=200)
    ap.add_argument("--nmin", type=int, default=0)
    ap.add_argument("--nmax", type=int, default=300)
    ap.add_argument("--measured", action="store_true", help="only nuclei with a measured mass")
    ap.add_argument("--no-zpe", dest="zpe", action="store_false", help="without the zero-point energy")
    ap.add_argument("--pairing", default="LN")
    ap.add_argument("--n_gh", type=int, default=44)
    ap.add_argument("--jobs", type=int, default=4, help="nuclei computed at the same time")
    ap.add_argument("--threads", type=int, default=1, help="OpenMP threads per nucleus")
    args = ap.parse_args()
    args.exe = os.path.abspath(args.exe)
    os.makedirs(os.path.join(args.workdir, "data"), exist_ok=True)

    done = set()
    if os.path.exists(args.out):
        with open(args.out) as f:
            for l in f:
                if not l.startswith("#"):
                    done.add((int(l.split()[0]), int(l.split()[1])))
    else:
        with open(args.out, "w") as f:
            f.write("# main.exe at the FRDM2012 ground-state shapes (gamma = 0), pairing %s, zpe %s\n"
                    % (args.pairing, args.zpe))
            f.write("# energies in MeV; M_th_FL and M_th from the FRDM2012 table, M_exp from AME2003\n")
            f.write(HEADER)

    todo = [r for r in read_table(args.table)
            if args.zmin <= r["Z"] <= args.zmax and args.nmin <= r["N"] <= args.nmax
            and (not args.measured or not math.isnan(r["M_exp"])) and (r["Z"], r["N"]) not in done]
    print("%d nuclei to compute" % len(todo))
    with ThreadPoolExecutor(max_workers=args.jobs) as pool, open(args.out, "a") as f:
        futures = [pool.submit(run_one, r, args) for r in todo]
        for i, fut in enumerate(as_completed(futures), 1):
            r = fut.result()
            f.write(line(r))
            f.flush()
            print("%5d/%d  Z = %3d N = %3d  M_th = %10.4f  (table FRLDM %9.2f)  %s"
                  % (i, len(todo), r["Z"], r["N"], r["M_calc"], r["M_th_FL"], r["flags"]))


if __name__ == "__main__":
    main()
