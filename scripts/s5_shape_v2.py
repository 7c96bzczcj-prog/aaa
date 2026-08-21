#!/usr/bin/env python3
"""
S-5 §2 / §3, CORRECTED design mapping.

A design error in the first attempt is corrected here and recorded rather
than quietly fixed. That attempt treated GSE228629's "skin" compartment as
carrying the same contrast as E-MTAB-10176. It does not. The methods of
that paper state that the photoconverted skin was "the lower back and
proximal parts of the legs", while "the rest of the body, including the
distal parts of the legs and hind paws, was shielded with aluminum foil";
the sequenced skin is HIND-PAW skin. So in GSE228629:
    RED   = cells that MIGRATED IN from the photoconverted skin
    GREEN = local hind-paw-skin cells
which is "migrant versus local in a destination tissue", with the polarity
INVERTED relative to E-MTAB-10176 (green = new arrival there).

Auditing every candidate the same way leaves exactly one public dataset
with the E-MTAB-10176 contrast (photoconvert tissue X, sample tissue X,
compare arrived-since versus retained):

  PRIMARY   GSE235702 -- mesenteric lymph node of Vav-H2B-Dendra2 mice
            photoconverted, and 7 days later D-Red+ (retained) and D-Red-
            (arrived since) CD19-CD8-CD4+CD62Llo cells sorted from the SAME
            node. Polarity matches: D-Red- == "green", D-Red+ == "red".

  DESCRIPTIVE ONLY  GSE228629 skin and joint -- different contrast, kept in
            the output with the polarity written out, NOT used for the
            §2 judgements.

Tumour side, for lineage matching: E-MTAB-10176 NK (the S-4 population),
CD8T, CD4T and Treg.

Pipeline identical to S-4 throughout: baseline-side gate only, binomial
thinning null at a level-matched c, residual rank position compared. Only
ranks are compared, never absolute levels.
"""
import warnings
import numpy as np
import pandas as pd
import anndata as ad
import scipy.io as sio
from scipy import stats

warnings.filterwarnings("ignore")
pd.set_option("display.width", 250)

OUT = "results"
RNG = np.random.default_rng(20260821)
MIN_P, MIN_POS, NSIM = 0.05, 10, 25
TARGET_SUM = 4648.0
EGRESS = ["Itgam", "S1pr5", "S1pr1", "Cxcr4", "Klf2", "Klf3"]
CHEMO = ["Ccl3", "Ccl4", "Ccl5", "Xcl1"]
S4 = {"Itgam": .9994, "S1pr5": .9993, "S1pr1": .9988, "Cxcr4": .9975,
      "Klf2": .9859, "Klf3": .9823, "Ccl3": .9986, "Ccl4": .9959,
      "Ccl5": .9956, "Xcl1": .1823}


def thin(c, m, rng):
    o = np.zeros(m.shape, np.int32)
    nz = m > 0
    o[nz] = rng.binomial(m[nz], c)
    return o


def analyse(G, R, genes, label, fh):
    pg, pr = (G > 0).mean(axis=0), (R > 0).mean(axis=0)
    npos = (G > 0).sum(axis=0)
    idx = np.flatnonzero((pg >= MIN_P) & (npos >= MIN_POS))
    g, r = pg[idx], pr[idx]
    drop = (g - r) / g
    med = np.median(drop)
    Gi = np.ascontiguousarray(G[:, idx])

    def nd(c, n):
        acc = np.zeros((n, len(idx)))
        for b in range(n):
            acc[b] = (thin(c, Gi, RNG) > 0).mean(axis=0)
        return (g - acc.mean(axis=0)) / g

    best, bg = 0.9, 1e9
    for c in np.round(np.arange(0.10, 1.00, 0.05), 2):
        gap = abs(np.median(nd(c, 8)) - med)
        if gap < bg:
            bg, best = gap, c
    bg, bd, bc = 1e9, None, best
    for c in np.round(np.arange(max(.02, best - .05), min(.99, best + .051), .01), 2):
        d = nd(c, NSIM)
        gap = abs(np.median(d) - med)
        if gap < bg or bd is None:
            bg, bc, bd = gap, c, d
    P = pd.DataFrame(dict(gene=genes[idx], p_green=g, p_red=r, drop_ratio=drop,
                          null_drop=bd, resid=drop - bd, n_pos_green=npos[idx]))
    P["resid_pct"] = P.resid.rank(pct=True)
    out = [f"\n### {label}",
           f"  A-side {G.shape[0]} cells, B-side {R.shape[0]} cells, "
           f"{len(idx)} genes gated, null c = {bc:.2f}",
           f"  drop ratio median {med:+.4f}, "
           f"Spearman(baseline, residual) = "
           f"{stats.spearmanr(g, P.resid).statistic:+.4f}",
           f"  {'gene':8s} {'set':7s} {'p_A':>7} {'p_B':>7} {'drop':>8} "
           f"{'resid':>8} {'resid pct':>10} | {'S-4 tumour':>10}"]
    for gene in EGRESS + CHEMO:
        row = P[P.gene == gene]
        st = "egress" if gene in EGRESS else "chemo"
        if not len(row):
            out.append(f"  {gene:8s} {st:7s}  -- not through the baseline gate")
            continue
        x = row.iloc[0]
        out.append(f"  {gene:8s} {st:7s} {x.p_green:7.4f} {x.p_red:7.4f} "
                   f"{x.drop_ratio:8.4f} {x.resid:8.4f} {x.resid_pct:10.4f} | "
                   f"{S4.get(gene, np.nan):10.4f}")
    for ln in out:
        print(ln); fh.write(ln + "\n")
    return P


def mtx(pref):
    feat = pd.read_csv(pref + "_features.tsv.gz", sep="\t", header=None,
                       names=["id", "name", "type"])
    M = sio.mmread(pref + "_matrix.mtx.gz").tocsc()
    return M, feat.name.values, feat.type.values


def main():
    Ps = {}
    with open(f"{OUT}/s5_shape_v2.log", "w") as fh:
        # ---------- tumour side ----------
        head = "=" * 78 + "\nTUMOUR SIDE: E-MTAB-10176 (green = arrived since, red = retained)\n" + "=" * 78
        print(head); fh.write(head + "\n")
        a = ad.read_h5ad("data/raw/DW_T_NK.h5ad")
        names = np.array(a.raw.var_names)
        X = a.raw.X.tocsr()
        nc = a.obs.n_counts.values.astype(np.float64)
        ct = a.obs.main_celltype.astype(str).values
        col = a.obs.colour.astype(str).values

        def cts(mask):
            i = np.flatnonzero(mask)
            return np.rint(np.expm1(X[i].toarray()) *
                           (nc[i][:, None] / TARGET_SUM)).astype(np.int32)
        for lin in ("NK", "CD8T", "CD4T", "Treg"):
            m = ct == lin
            G, R = cts(m & (col == "Green")), cts(m & (col == "Red"))
            if min(G.shape[0], R.shape[0]) < 60:
                ln = f"\n### tumour {lin}: too few cells ({G.shape[0]}/{R.shape[0]}) -- skipped"
                print(ln); fh.write(ln + "\n"); continue
            Ps[f"tumour_{lin}"] = analyse(G, R, names, f"tumour {lin}", fh)

        # ---------- PRIMARY non-tumour: GSE235702 ----------
        head = ("\n" + "=" * 78 +
                "\nNON-TUMOUR PRIMARY: GSE235702, mesenteric LN, SAME contrast"
                "\n  A-side = D-Red- (arrived in the node since photoconversion) == 'green'"
                "\n  B-side = D-Red+ (retained in the node)                      == 'red'"
                "\n  interval 7 days (E-MTAB-10176: 24 and 72 h) -- see §3\n" + "=" * 78)
        print(head); fh.write(head + "\n")
        Mg, gn, _ = mtx("data/raw/s5/GSM7508111_sample1")
        Mr, rn, _ = mtx("data/raw/s5/GSM7508112_sample2")
        assert (gn == rn).all()
        G = np.asarray(Mg.todense()).T.astype(np.int32)
        R = np.asarray(Mr.todense()).T.astype(np.int32)
        Ps["LN_CD4T"] = analyse(G, R, gn, "non-tumour LN CD4 T (GSE235702)", fh)
        Ps["LN_CD4T"].to_csv(f"{OUT}/s5_LN_CD4T.csv", index=False)

        # ---------- DESCRIPTIVE ONLY: GSE228629 ----------
        head = ("\n" + "=" * 78 +
                "\nDESCRIPTIVE ONLY: GSE228629 -- DIFFERENT CONTRAST, NOT USED FOR §2"
                "\n  photoconverted site = lower back + proximal legs; hind paws SHIELDED"
                "\n  A-side = Kaede GREEN = LOCAL hind-paw cells"
                "\n  B-side = Kaede RED   = cells that MIGRATED IN from the converted skin"
                "\n  polarity is INVERTED relative to the tumour data\n" + "=" * 78)
        print(head); fh.write(head + "\n")
        for lib, pref in (("B6", "data/raw/s5/GSM7134137_B6"),
                          ("BALBc", "data/raw/s5/GSM7134139_C")):
            M, nm, ty = mtx(pref)
            hto = np.flatnonzero(ty == "Antibody Capture")
            gex = np.flatnonzero(ty == "Gene Expression")
            H = np.asarray(M[hto].todense()).astype(float)
            tot, top = H.sum(0), H.max(0)
            sec = np.sort(H, 0)[-2]
            call = np.where((tot >= 20) & (top >= 3 * np.maximum(sec, 1)),
                            np.array(nm[hto])[H.argmax(0)], "UNASSIGNED")
            Xg = M[gex]
            for comp in ("skin", "joint"):
                A = np.asarray(Xg[:, np.flatnonzero(call == f"{comp}_CD45pos_GREEN")]
                               .todense()).T.astype(np.int32)
                B = np.asarray(Xg[:, np.flatnonzero(call == f"{comp}_CD45pos_RED")]
                               .todense()).T.astype(np.int32)
                Ps[f"{comp}_{lib}"] = analyse(
                    A, B, nm[gex],
                    f"GSE228629 {comp} / {lib}  [A=LOCAL, B=MIGRANT -- inverted]", fh)

        # ---------- side by side ----------
        head = ("\n" + "=" * 78 + "\nSIDE BY SIDE (residual percentile)\n"
                "  comparable pair: tumour lineages  vs  LN_CD4T  (same contrast)\n"
                "  GSE228629 columns are DIFFERENT CONTRAST, shown for completeness\n"
                + "=" * 78)
        print(head); fh.write(head + "\n")
        keys = [k for k in ("tumour_NK", "tumour_CD8T", "tumour_CD4T", "tumour_Treg",
                            "LN_CD4T", "skin_B6", "skin_BALBc",
                            "joint_B6", "joint_BALBc") if k in Ps]
        hdr = f"{'gene':8s} {'set':7s} " + " ".join(f"{k:>13s}" for k in keys)
        print(hdr); fh.write(hdr + "\n")
        for gene in EGRESS + CHEMO:
            cells = []
            for k in keys:
                row = Ps[k][Ps[k].gene == gene]
                cells.append(f"{row.resid_pct.iloc[0]:.4f}" if len(row) else "gated-out")
            st = "egress" if gene in EGRESS else "chemo"
            ln = f"{gene:8s} {st:7s} " + " ".join(f"{c:>13s}" for c in cells)
            print(ln); fh.write(ln + "\n")
    print("\nwrote s5_shape_v2.log")


if __name__ == "__main__":
    main()
