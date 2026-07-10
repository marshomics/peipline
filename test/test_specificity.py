#!/usr/bin/env python3
"""Verify the target-specificity analyses against a dataset with a known answer.

What is planted, and what must therefore be recovered:

  * a PeiP-like structure with the catalytic triad and two extra "specificity"
    residues clustered in space, and the rest of the chain far away, so the
    groove shell is unambiguous;
  * SDP columns whose residue is determined by a HOST GENUS label, which the
    active-site barcode never saw;
  * decoy columns whose residue is determined by the barcode subgroup and
    nothing else -- these are the circularity trap, and they must NOT be
    reported as replicated SDPs;
  * two PMBR architectures (2 repeats vs 4 repeats) that are independent of the
    catalytic sequence, so the module trees must come out incongruent;
  * the genus-driven SDP columns placed INSIDE the groove and the barcode-driven
    decoys placed far outside it, so the Fisher enrichment is a real test;
  * two trees, a star and a balanced one, so both the degenerate-null fallback
    and the real within-clade null are exercised.

Run:  python test/test_specificity.py /tmp/c71_spec
"""
from __future__ import annotations

import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd
import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(ROOT, "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, os.path.join(ROOT, "test"))

FAILURES = []
AA = "ACDEFGHIKLMNPQRSTVWY"

L = 300                       # match columns
TRIAD = [112, 147, 171]
SPEC_SEEDS = [150, 168]       # the "Y174/V252" analogues
# Wang et al.'s class positions. position_1 sits 2 residues after the catalytic
# Asp, exactly as V252 sits after D250 -- i.e. INSIDE the +/-5 triad-flank
# barcode. That is the circularity hazard pei_class introduces, and the test
# checks it is excluded from the SDP column set.
CLASS_POS = [173, 186]        # (auth 174, 187)
CLASS_RES = {"I": ("V", "C"), "II": ("V", "V"), "III": ("T", "I"), "IV": ("A", "M")}
# A divalent-cation site placed FAR from the substrate groove. PeiW and PeiP both
# need a metal and disagree about which (Subedi et al. 2015), so the metal shell
# is a second hypothesis space. Planting it 40 A from the groove is what makes
# "the two shells are independent" a claim the test can check rather than assert.
METAL_POS = [230, 234, 238]
METAL_CENTRE = np.array([40.0, 0.0, 0.0])
# columns whose residue is a function of host genus  -> true SDPs
SDP_COLS = [145, 149, 152, 170, 173]
# columns whose residue is a function of the barcode subgroup -> circular decoys
CIRCULAR_COLS = [40, 60, 80, 200, 220]
GENERA = ["g__Methanobrevibacter", "g__Methanosphaera", "g__Methanothermobacter"]
SDP_AA = {0: "GSTNQ", 1: "AVILM", 2: "DEKRH"}      # genus -> residue pool
SUB_AA = {0: "W", 1: "F", 2: "Y"}                  # subgroup -> residue


def check(cond, msg):
    if not cond:
        print(f"!! {msg}")
        FAILURES.append(msg)
    return cond


def run(cmd, env=None):
    print("\n$ " + " ".join(map(str, cmd)), flush=True)
    r = subprocess.run(list(map(str, cmd)), env=env)
    if r.returncode != 0:
        sys.exit(f"FAILED: {' '.join(map(str, cmd))}")


# ---------------------------------------------------------------------------
def make_structure(path, seq):
    """A synthetic PeiP: the triad and the two specificity seeds sit within 6 A
    of each other; every other residue is >20 A away. So the groove shell is
    exactly {triad} u {seeds} u nothing else, and any leakage is a bug."""
    rng = np.random.default_rng(3)
    coords, names = {}, {}
    # The genus-driven SDP columns sit IN the groove; the barcode-driven decoys
    # sit far away. So a groove-enrichment test that passes is a test that means
    # something.
    site = TRIAD + SPEC_SEEDS + SDP_COLS + CLASS_POS
    # place the site residues in a tight cluster around the origin
    for k, c in enumerate(site):
        coords[c] = np.array([0.0, 0.0, 0.0]) + rng.normal(0, 1.6, 3)
    # A separate metal-binding cluster, 40 A away, with the ion at its centre.
    # Placed deterministically at 2.4 A, inside the 3.2 A coordination radius: a
    # random jitter would put a donor outside the first shell on some seeds and
    # the test would fail for reasons that have nothing to do with the code.
    for c, axis in zip(METAL_POS, ([2.4, 0, 0], [0, 2.4, 0], [0, 0, 2.4])):
        coords[c] = METAL_CENTRE + np.array(axis, dtype=float)
    # everything else on a distant sphere
    for c in range(L):
        if c in coords:
            continue
        v = rng.normal(size=3)
        v /= np.linalg.norm(v)
        coords[c] = v * rng.uniform(30, 60)
        # keep the filler away from the metal, or the shell is not the cluster
        while np.linalg.norm(coords[c] - METAL_CENTRE) < 12.0:
            v = rng.normal(size=3)
            v /= np.linalg.norm(v)
            coords[c] = v * rng.uniform(30, 60)
    # The chain must carry a REAL sequence, not poly-Ala: hmmalign --trim would
    # throw away a poly-Ala chain as "outside the profile envelope", and the
    # groove would silently come back empty. Use one of the actual sequences.
    ONE_TO_THREE = {v: k for k, v in
                    {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
                     "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
                     "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
                     "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y",
                     "VAL": "V"}.items()}
    fixed = {112: "C", 147: "H", 171: "D", 150: "Y", 168: "V",
             CLASS_POS[0]: "V", CLASS_POS[1]: "C"}   # reference is class I
    # give the metal three O/N donors, so `n_plausible_ligands` >= 2 and the site
    # is not flagged as a modelling artefact
    fixed.update({METAL_POS[0]: "D", METAL_POS[1]: "E", METAL_POS[2]: "N"})
    for c in range(L):
        names[c] = ONE_TO_THREE[fixed.get(c, seq[c])]

    with open(path, "w") as fh:
        fh.write("data_TEST\n#\nloop_\n")
        for f in ("group_PDB", "id", "type_symbol", "label_atom_id",
                  "label_comp_id", "label_asym_id", "label_seq_id",
                  "auth_asym_id", "auth_seq_id", "Cartn_x", "Cartn_y", "Cartn_z"):
            fh.write(f"_atom_site.{f}\n")
        i = 0
        for c in range(L):
            # auth numbering is 1-based, and the pipeline's match columns are
            # 0-based, so residue c+1 corresponds to match column c
            for atom in ("N", "CA", "C"):
                i += 1
                x, y, z = coords[c] + (0.0 if atom == "CA" else 0.4)
                fh.write(f"ATOM {i} {atom[0]} {atom} {names[c]} A {c+1} A {c+1} "
                         f"{x:.3f} {y:.3f} {z:.3f}\n")
        # The calcium. Note the trap this exercises: its label_atom_id is "CA"
        # and so is the alpha carbon's. Only label_comp_id distinguishes them,
        # and a reader that keys on atom name will pull 300 alpha carbons into
        # the ion list.
        i += 1
        x, y, z = METAL_CENTRE
        fh.write(f"HETATM {i} CA CA CA B 1 B 901 {x:.3f} {y:.3f} {z:.3f}\n")
        fh.write("#\n")
    return coords


def make_sequences(n=240):
    """Barcode subgroup and host genus are assigned INDEPENDENTLY, so a column
    driven by one carries no information about the other. That independence is
    what makes the circularity trap a real trap."""
    rng = np.random.default_rng(11)
    template = "".join(rng.choice(list(AA), L))
    rows, meta = [], []
    for i in range(n):
        genus = i % 3
        subgroup = (i // 3) % 3          # independent of genus
        s = list(template)
        for j in range(L):
            if rng.random() < 0.25:
                s[j] = rng.choice(list(AA))
        for c, r in zip(TRIAD, "CHD"):
            s[c] = r
        for c in SDP_COLS:
            s[c] = rng.choice(list(SDP_AA[genus]))
        for c in CIRCULAR_COLS:
            s[c] = SUB_AA[subgroup]
        cls = ["I", "I", "I", "II", "III", "IV"][i % 6]   # class I is the majority
        s[CLASS_POS[0]], s[CLASS_POS[1]] = CLASS_RES[cls]
        sid = f"c71_{i:07d}"
        rows.append((sid, "".join(s)))
        meta.append((sid, f"SAMPLE_{i:03d}", GENERA[genus], subgroup,
                     4 if i % 2 == 0 else 2, cls))
    return rows, meta


def main() -> None:
    base = sys.argv[1] if len(sys.argv) > 1 else "/tmp/c71_spec"
    work = os.path.join(base, "work")
    tab = os.path.join(base, "outputs", "tables")
    fig = os.path.join(base, "outputs", "figures")
    for d in (work, tab, fig):
        os.makedirs(d, exist_ok=True)

    # groove_map.py shells out to hmmalign. Install the stubs here rather than
    # assuming a seeded PATH: a test that only passes when someone has already
    # run something else is not a test.
    import stubs
    stubs.install(os.path.join(base, "bin"))

    # --- build the world ----------------------------------------------------
    rows, meta = make_sequences()
    cif = os.path.join(base, "peip.cif")
    make_structure(cif, rows[0][1])

    afa = os.path.join(work, "triad_pass_matchcols.afa")
    with open(afa, "w") as fh:
        for sid, s in rows:
            fh.write(f">{sid}\n{s}\n")

    # a profile that is all match states, so match column == residue index
    import pyhmmer
    from pyhmmer.easel import Alphabet, TextMSA, TextSequence
    abc = Alphabet.amino()
    msa = TextMSA(name=b"Peptidase_C71",
                  sequences=[TextSequence(name=f"t{i}".encode(), sequence=s)
                             for i, (_, s) in enumerate(rows[:60])]).digitize(abc)
    hmm, _, _ = pyhmmer.plan7.Builder(abc).build_msa(msa, pyhmmer.plan7.Background(abc))
    hmm.name = b"Peptidase_C71"
    hmm_path = os.path.join(base, "PF12386.hmm")
    with open(hmm_path, "wb") as fh:
        hmm.write(fh)

    with open(os.path.join(work, "triad_columns.json"), "w") as fh:
        json.dump({"residues": ["C", "H", "D"], "match_columns": TRIAD,
                   "n_match_columns": L, "n_input_sequences": len(rows),
                   "n_triad_positive": len(rows)}, fh)

    md = pd.DataFrame(meta, columns=["seq_id", "sample", "gtdb_genus",
                                     "subgroup", "n_pmbr", "true_class"])
    md.assign(cluster=md["seq_id"], weight=1.0)[
        ["seq_id", "cluster", "weight"]].assign(cluster_size=1).to_csv(
        os.path.join(tab, "sequence_weights.tsv"), sep="\t", index=False)
    md[["seq_id", "sample"]].assign(protein_id=md["seq_id"], faa="x",
                                    evidence="specific").to_csv(
        os.path.join(work, "hits_idmap.tsv.gz"), sep="\t", index=False)
    md[["seq_id", "subgroup"]].assign(barcode="X", weight=1.0,
                                      cluster=md["seq_id"]).to_csv(
        os.path.join(tab, "subgroup_assignments.tsv"), sep="\t", index=False)
    md[["sample", "gtdb_genus"]].to_csv(os.path.join(base, "merged.tsv"),
                                        sep="\t", index=False)
    # The planted architectures are 2 and 4 repeats, which sit either side of the
    # 3-motif binding cliff. Derive the binding columns from pmbr_reference so the
    # test cannot drift from the rule it is testing.
    sys.path.insert(0, SCRIPTS)
    from pmbr_reference import predict_binding  # noqa: E402
    _bind = md["n_pmbr"].map(lambda n: predict_binding(n))
    md[["seq_id"]].assign(
        n_pmbr=md["n_pmbr"],
        n_pmbr_permissive=md["n_pmbr"],
        pmbr_count_fragile=0,
        pmbr_binding_competent=_bind.map(lambda t: int(t[1])),
        predicted_binding=_bind.map(lambda t: t[0]),
        pmbr_pi=7.0,
        accessory_binding_domains="",
        architecture_class=md["n_pmbr"].map(
            lambda n: "canonical_pei" if n >= 3 else "reduced_pmbr"),
        pmbr_span="", catalytic_span="").to_csv(
        os.path.join(tab, "domain_architecture.tsv"), sep="\t", index=False)
    # SSN clusters: correlate with genus, so it is a second, non-barcode partition
    md[["seq_id"]].assign(ssn_cluster=md["gtdb_genus"].map(
        {g: i for i, g in enumerate(GENERA)})).to_csv(
        os.path.join(tab, "ssn_clusters.tsv"), sep="\t", index=False)

    # Two trees. The star tree has no clade structure, so the within-clade null
    # is degenerate and the code must detect that and fall back. The balanced
    # tree has real clades, so the within-clade null must run and still recover
    # the planted SDPs. Both paths get tested.
    with open(os.path.join(base, "star.treefile"), "w") as fh:
        fh.write("(" + ",".join(f"{sid}:0.1" for sid, _ in rows) + ");\n")

    def balanced(ids, depth=0):
        if len(ids) == 1:
            return f"{ids[0]}:0.05"
        h = len(ids) // 2
        return f"({balanced(ids[:h], depth+1)},{balanced(ids[h:], depth+1)}):0.05"
    with open(os.path.join(base, "balanced.treefile"), "w") as fh:
        fh.write(balanced([sid for sid, _ in rows]) + ";\n")

    cfg = yaml.safe_load(open(os.path.join(ROOT, "config.yaml")))
    cfg["specificity"].update(
        structure=cif, structure_chain="A", structure_numbering="test",
        groove_seed_residues=[c + 1 for c in SPEC_SEEDS],
        structure_expect={113: "C", 148: "H", 172: "D", 151: "Y", 169: "V",
                          CLASS_POS[0] + 1: "V", CLASS_POS[1] + 1: "C"},
        groove_literature_residues=[c + 1 for c in TRIAD + SPEC_SEEDS + SDP_COLS],
        groove_radius_a=8.0, contact_radius_a=8.0, structure_offset=0,
        sdp_partitions=["ssn_cluster", "host_genus", "pmbr_architecture", "pei_class"],
        sdp_min_partitions=2, sdp_permutations=150, sdp_min_group_size=10,
        panel_size=6)
    cfg["specificity"]["pei_class"] = {
        "enabled": True, "position_1": CLASS_POS[0] + 1, "position_2": CLASS_POS[1] + 1,
        "position_2_alt": None,
        "classes": {"I": {"p1": ["V"], "p2": ["C"]}, "II": {"p1": ["V"], "p2": ["V"]},
                    "III": {"p1": ["T", "S"], "p2": ["I"]},
                    "IV": {"p1": ["A"], "p2": ["M", "W"]}}}
    cfg["outputs"]["prodigal_dir"] = base
    tcfg = os.path.join(base, "config.test.yaml")
    yaml.safe_dump(cfg, open(tcfg, "w"))

    # --- groove_map ----------------------------------------------------------
    run([sys.executable, f"{SCRIPTS}/groove_map.py", "--config", tcfg,
         "--align-hmm", hmm_path, "--chosen", f"{work}/triad_columns.json",
         "--workdir", f"{work}/groove", "--out-columns", f"{tab}/groove_columns.tsv",
         "--out-contacts", f"{work}/groove_contacts.tsv",
         "--out-json", f"{tab}/groove_definition.json"],
        env=dict(os.environ, PATH=f"{base}/bin:{os.environ['PATH']}"))

    gr = pd.read_csv(f"{tab}/groove_columns.tsv", sep="\t")
    gdef = json.load(open(f"{tab}/groove_definition.json"))
    got = set(gr.loc[gr["in_groove"] == 1, "match_col"])
    want = set(TRIAD) | set(SPEC_SEEDS) | set(SDP_COLS) | set(CLASS_POS)
    print("\n" + "=" * 72)
    print(f"groove columns recovered: {sorted(got)}")
    print(f"planted site            : {sorted(want)}")
    check(got == want, f"groove shell wrong: got {sorted(got)}, want {sorted(want)}")
    check(gdef["triad_maps_to_CHD"], "triad columns did not map to C/H/D in the structure")

    D = np.loadtxt(f"{work}/groove_contacts.tsv", delimiter="\t")
    check(D[TRIAD[0], TRIAD[1]] < 8.0, "triad residues not in contact in the Ca matrix")
    check(D[TRIAD[0], 5] > 20.0, "a distant residue is reported as nearby")

    # --- the divalent metal site --------------------------------------------
    import copy as _c, subprocess as _sp, yaml as _y  # noqa: E401

    # The ion's label_atom_id is "CA", exactly like every alpha carbon. If the
    # reader keys on atom name instead of residue name, `ions_found` picks up
    # 300 backbone atoms and the "metal shell" becomes the whole protein.
    met = gdef["metal"]
    check(len(met["ions_found"]) == 1,
          f"expected exactly 1 ion, found {len(met['ions_found'])} "
          f"-- a CA alpha-carbon was probably parsed as calcium")
    if met["ions_found"]:
        ion = met["ions_found"][0]
        check(ion["ion"] == "CA", f"ion misread as {ion['ion']}")
        check(ion["n_plausible_ligands"] >= 2,
              "the planted D/E/N donors were not recognised as ligands")
    mshell = set(gr.loc[gr["in_metal_shell"] == 1, "match_col"])
    mcoord = set(gr.loc[gr["metal_coordinating"] == 1, "match_col"])
    check(set(METAL_POS) <= mshell,
          f"metal shell {sorted(mshell)} misses planted {METAL_POS}")
    check(set(METAL_POS) <= mcoord,
          f"coordinating set {sorted(mcoord)} misses planted {METAL_POS}")
    check(not (mshell & got),
          f"metal shell and substrate groove overlap at {sorted(mshell & got)}; "
          f"they were planted 40 A apart, so this is a coordinate-frame bug")
    check(met["shell_overlaps_groove"] == 0, "shell_overlaps_groove should be 0")
    print(f"metal shell: {sorted(mshell)}  (groove is {sorted(got)}, disjoint)")

    # ...and it must degrade cleanly when no ion is modelled, which is the real
    # case: neither 8JX4 nor 8Z4F resolves the cation.
    nometal = _c.deepcopy(cfg)
    nometal["specificity"]["metal_site"] = {"enabled": True, "ions": [],
                                            "coordination_radius_a": 3.2,
                                            "shell_radius_a": 8.0,
                                            "coordinating_residues": list("DENQSTH")}
    np_ = os.path.join(base, "config.nometal.yaml")
    _y.safe_dump(nometal, open(np_, "w"))
    r = _sp.run([sys.executable, f"{SCRIPTS}/groove_map.py", "--config", np_,
                 "--align-hmm", hmm_path, "--chosen", f"{work}/triad_columns.json",
                 "--workdir", f"{work}/groove_nm",
                 "--out-columns", f"{work}/groove_columns_nm.tsv",
                 "--out-contacts", f"{work}/groove_contacts_nm.tsv",
                 "--out-json", f"{work}/groove_definition_nm.json"],
                capture_output=True, text=True,
                env=dict(os.environ, PATH=f"{base}/bin:{os.environ['PATH']}"))
    check(r.returncode == 0, "groove_map must not fail when no metal is modelled")
    gnm = pd.read_csv(f"{work}/groove_columns_nm.tsv", sep="\t")
    check("in_metal_shell" in gnm.columns and gnm["in_metal_shell"].sum() == 0,
          "with no ion, in_metal_shell must exist and be all zero, so downstream "
          "code is uniform and a missing site cannot masquerade as a found one")
    check("no" in r.stderr and "EDTA" in r.stderr,
          "the no-metal path must say why the absence matters, not fail silently")
    print("no-metal path: columns emitted as zeros, absence explained")

    # --- structure_expect must REJECT a wrong structure ----------------------
    badcfg = _c.deepcopy(cfg)
    badcfg["specificity"]["structure_expect"] = dict(cfg["specificity"]["structure_expect"])
    badcfg["specificity"]["structure_expect"][113] = "A"   # claim Ala where Cys sits
    bp = os.path.join(base, "config.bad.yaml")
    _y.safe_dump(badcfg, open(bp, "w"))
    r = _sp.run([sys.executable, f"{SCRIPTS}/groove_map.py", "--config", bp,
                 "--align-hmm", hmm_path, "--chosen", f"{work}/triad_columns.json",
                 "--workdir", f"{work}/groove_bad", "--out-columns", "/dev/null",
                 "--out-contacts", "/dev/null", "--out-json", "/dev/null"],
                capture_output=True, text=True,
                env=dict(os.environ, PATH=f"{base}/bin:{os.environ['PATH']}"))
    print("\n" + "=" * 72)
    check(r.returncode != 0, "groove_map accepted a structure whose residues do not "
                             "match structure_expect")
    check("8JX4" in r.stderr and "8Z4F" in r.stderr,
          "the mismatch error does not name the two PDB entries")
    print("structure_expect guard fires:", r.stderr.strip().splitlines()[0][:88])

    # --- pei_class ------------------------------------------------------------
    run([sys.executable, f"{SCRIPTS}/pei_class.py", "--afa", afa,
         "--groove-json", f"{tab}/groove_definition.json",
         "--assign", f"{tab}/subgroup_assignments.tsv", "--ssn", f"{tab}/ssn_clusters.tsv",
         "--weights", f"{tab}/sequence_weights.tsv", "--config", tcfg,
         "--figdir", fig, "--out", f"{tab}/pei_class.tsv",
         "--out-agreement", f"{tab}/pei_class_vs_subgroup.tsv"])
    pc = pd.read_csv(f"{tab}/pei_class.tsv", sep="\t")
    truth_cls = dict(zip(md["seq_id"], md["true_class"]))
    pc["truth"] = pc["seq_id"].map(truth_cls)
    acc = (pc["pei_class"] == pc["truth"]).mean()
    print("\n" + "=" * 72)
    print(f"pei_class accuracy vs planted classes: {acc:.3f}")
    print(pc["pei_class"].value_counts().to_string())
    check(acc == 1.0, f"pei_class mis-assigned {(1-acc)*len(pc):.0f} sequences")
    check(set(pc["pei_class"]) == {"I", "II", "III", "IV"},
          "not all four published classes were recovered")

    # --- sdp must refuse an architecture table with no binding call ----------
    # Silently falling back to a repeat-count partition would hide the fact that
    # the 3-motif cliff was never applied.
    stale = pd.read_csv(f"{tab}/domain_architecture.tsv", sep="\t").drop(
        columns=["pmbr_binding_competent"])
    stale_p = os.path.join(base, "arch_stale.tsv")
    stale.to_csv(stale_p, sep="\t", index=False)
    r = _sp.run([sys.executable, f"{SCRIPTS}/sdp.py", "--afa", afa,
                 "--assign", f"{tab}/subgroup_assignments.tsv", "--arch", stale_p,
                 "--ssn", f"{tab}/ssn_clusters.tsv",
                 "--tree", f"{base}/balanced.treefile",
                 "--merged", f"{base}/merged.tsv",
                 "--idmap", f"{work}/hits_idmap.tsv.gz",
                 "--groove", f"{tab}/groove_columns.tsv",
                 "--groove-json", f"{tab}/groove_definition.json",
                 "--pei-class", f"{tab}/pei_class.tsv",
                 "--weights", f"{tab}/sequence_weights.tsv", "--config", tcfg,
                 "--figdir", fig, "--out", os.path.join(base, "sdp_stale.tsv"),
                 "--out-concordance", os.path.join(base, "conc_stale.tsv")],
                capture_output=True, text=True)
    check(r.returncode != 0,
          "sdp.py must exit non-zero when domain_architecture.tsv lacks "
          "pmbr_binding_competent, not quietly partition on repeat counts")
    check("pmbr_binding_competent" in r.stderr,
          "the error must name the missing column")
    print("sdp guard: refuses an architecture table with no 3-motif binding call")

    # a nested/unknown partition mode must also be rejected
    badmode = _c.deepcopy(cfg)
    badmode["specificity"]["pmbr_partition_mode"] = "both"
    bmp = os.path.join(base, "config.badmode.yaml")
    _y.safe_dump(badmode, open(bmp, "w"))
    r = _sp.run([sys.executable, f"{SCRIPTS}/sdp.py", "--afa", afa,
                 "--assign", f"{tab}/subgroup_assignments.tsv",
                 "--arch", f"{tab}/domain_architecture.tsv",
                 "--ssn", f"{tab}/ssn_clusters.tsv",
                 "--tree", f"{base}/balanced.treefile",
                 "--merged", f"{base}/merged.tsv",
                 "--idmap", f"{work}/hits_idmap.tsv.gz",
                 "--groove", f"{tab}/groove_columns.tsv",
                 "--groove-json", f"{tab}/groove_definition.json",
                 "--pei-class", f"{tab}/pei_class.tsv",
                 "--weights", f"{tab}/sequence_weights.tsv", "--config", bmp,
                 "--figdir", fig, "--out", os.path.join(base, "sdp_bad.tsv"),
                 "--out-concordance", os.path.join(base, "conc_bad.tsv")],
                capture_output=True, text=True)
    check(r.returncode != 0 and "pmbr_partition_mode" in r.stderr,
          "sdp.py must reject pmbr_partition_mode='both': count and "
          "binding_competent are nested and would double-count one observation")
    print("sdp guard: refuses pmbr_partition_mode='both' (nested partitions)")

    # --- sdp, on both trees --------------------------------------------------
    for tree_kind, want_null in (("star", "global_fallback"),
                                 ("balanced", "within_clade")):
        run([sys.executable, f"{SCRIPTS}/sdp.py", "--afa", afa,
             "--assign", f"{tab}/subgroup_assignments.tsv",
             "--arch", f"{tab}/domain_architecture.tsv",
             "--ssn", f"{tab}/ssn_clusters.tsv",
             "--tree", f"{base}/{tree_kind}.treefile",
             "--merged", f"{base}/merged.tsv", "--idmap", f"{work}/hits_idmap.tsv.gz",
             "--groove", f"{tab}/groove_columns.tsv",
             "--groove-json", f"{tab}/groove_definition.json",
             "--pei-class", f"{tab}/pei_class.tsv",
             "--weights", f"{tab}/sequence_weights.tsv", "--config", tcfg,
             "--figdir", fig, "--out", f"{tab}/sdp_replicated.tsv",
             "--out-concordance", f"{tab}/sdp_concordance.tsv"])

        sdp = pd.read_csv(f"{tab}/sdp_replicated.tsv", sep="\t")
        pp = pd.read_csv(f"{tab}/sdp_replicated_per_partition.tsv", sep="\t")
        rep = set(sdp.loc[sdp["replicated"], "match_col"])
        nulls = set(pp.loc[pp["partition"] != "tree_clade", "null_mode"])
        print("\n" + "=" * 72)
        print(f"--- {tree_kind} tree ---")
        print(f"null modes used      : {sorted(nulls)}")
        print(f"replicated SDPs      : {sorted(rep)}")
        print(f"planted (genus-driven): {sorted(SDP_COLS)}")
        print(f"circular decoys      : {sorted(CIRCULAR_COLS)}")

        check(want_null in nulls,
              f"{tree_kind} tree: expected null '{want_null}', got {sorted(nulls)}")

        # the pei_class partition must not be allowed to test its own columns
        pcl = pp[pp["partition"] == "pei_class"]
        check(len(pcl) > 0, f"{tree_kind}: pei_class partition never ran")
        tested = set(pcl["match_col"])
        check(not (tested & set(CLASS_POS)),
              f"{tree_kind}: CIRCULARITY: pei_class tested its own defining columns "
              f"{sorted(tested & set(CLASS_POS))}")
        others = pp[pp["partition"] != "pei_class"]
        check(set(CLASS_POS) <= set(others["match_col"]),
              f"{tree_kind}: the class columns were excluded from the OTHER partitions "
              f"too; only pei_class should exclude them")
        print(f"pei_class excluded its defining columns {sorted(CLASS_POS)}; "
              f"other partitions still test them")
        found = rep & set(SDP_COLS)
        circ = rep & set(CIRCULAR_COLS)
        check(len(found) >= 4, f"{tree_kind}: only {len(found)}/5 planted SDPs recovered")
        check(not circ, f"{tree_kind}: CIRCULARITY: barcode-driven decoys {sorted(circ)} "
                        f"were called replicated SDPs. The partitions leaked.")
        fp = rep - set(SDP_COLS) - set(CIRCULAR_COLS)
        check(len(fp) <= 3, f"{tree_kind}: {len(fp)} false-positive SDPs: {sorted(fp)[:5]}")

        summ = pd.read_csv(f"{tab}/sdp_replicated_summary.tsv", sep="\t",
                           index_col="metric")["value"]
        orr = float(summ.get("groove_enrichment_or", np.nan))
        pv = float(summ.get("groove_enrichment_p", np.nan))
        print(f"groove enrichment    : OR = {orr:.2f}, p = {pv:.3g}")
        check(pv < 0.05, f"{tree_kind}: replicated SDPs are not enriched in the "
                         f"groove (p = {pv:.3g}); the groove/SDP intersection is broken")

    # --- assay panel ---------------------------------------------------------
    # Three wall chemistries, drawn from Kandler & Koenig 1978 via
    # cellwall_reference.py, so the panel must span both the P1 and P1' axes.
    import cellwall_reference as cwr
    spp = ["s__Methanobrevibacter smithii",      # Ala / Lys-Orn
           "s__Methanobrevibacter ruminantium",  # Thr / Lys
           "s__Methanothermobacter thermautotrophicus"]  # Ala / Lys
    cw = md.copy()
    cw["species"] = [spp[i % 3] for i in range(len(cw))]
    chem = cwr.annotate(cw["species"])
    cw = pd.concat([cw, chem], axis=1)
    cw.assign(pathway_call="pseudomurein", completeness=99.0, contamination=0.0,
              domain="Archaea", has_c71=1).to_csv(
        f"{tab}/cellwall_genotype.tsv", sep="\t", index=False)
    run([sys.executable, f"{SCRIPTS}/assay_panel.py", "--afa", afa,
         "--sdp", f"{tab}/sdp_replicated.tsv", "--arch", f"{tab}/domain_architecture.tsv",
         "--cellwall", f"{tab}/cellwall_genotype.tsv", "--idmap", f"{work}/hits_idmap.tsv.gz",
         "--assign", f"{tab}/subgroup_assignments.tsv",
         "--weights", f"{tab}/sequence_weights.tsv", "--config", tcfg,
         "--figdir", fig, "--out", f"{tab}/assay_panel.tsv"])
    panel = pd.read_csv(f"{tab}/assay_panel.tsv", sep="\t")
    check(len(panel) == 6, f"panel has {len(panel)} members, wanted 6")
    check(panel["sdp_residues"].nunique() == len(panel),
          "panel members are not distinct in SDP space")
    check(panel["prediction_p1"].notna().all(), "panel rows lack a P1 prediction")
    check("p1_prime_residue" in panel.columns,
          "panel does not carry the acyl-acceptor (P1') axis")
    check(panel["p1_residue"].nunique() >= 2,
          "panel does not span more than one P1 chemistry")
    # Glu-gamma-Thr-pNA EXISTS (Subedi et al. 2015 bought it from JPT) and
    # neither enzyme cleaves it, while PeiW does cleave the Thr isopeptide. So a
    # Thr-host protein must be routed to the isopeptide, and the panel must warn
    # against the chromogen by name rather than pretend it was never made.
    thr = panel[panel["p1_residue"] == "Thr"]
    check(len(thr) and thr["prediction_p1"].str.contains("ISOPEPTIDE").all(),
          "a Thr-host protein is not asked to be assayed on the Thr isopeptide")
    check(thr["substrate_format"].str.contains("isopeptide").all(),
          "substrate_format must send a Thr-host protein to the isopeptide assay")
    check(thr["prediction_p1"].str.contains("Do NOT use EgammaT-pNA").all(),
          "the panel must warn against EgammaT-pNA explicitly: it exists, it was "
          "assayed, and neither characterised enzyme cleaves it")
    ala = panel[panel["p1_residue"].isin(["Ala", "Ser"])]
    if len(ala):
        check(ala["substrate_format"].str.startswith("pNA").all(),
              "Ala and Ser at P1 are cleaved in the pNA series; the cheap "
              "continuous assay is valid and the panel should say so")
    check(panel["assay_metal"].str.contains("PeiP by Ca only").all(),
          "every panel member must be assayed against the cation series: PeiW is "
          "rescued by five cations, PeiP by one")
    # No published substrate carries an Orn acceptor.
    orn = panel[panel["p1_prime_residue"] == "Lys/Orn"]
    check(len(orn) == 0 or orn["prediction_p1_prime"].str.contains("Orn").all(),
          "an Orn-host protein does not flag that no published substrate tests P1'")
    print(f"\npanel: {len(panel)} proteins, "
          f"{panel['sdp_residues'].nunique()} distinct SDP strings")

    print("=" * 72)
    print("\nRESULT:", "PASS" if not FAILURES else f"FAIL ({len(FAILURES)})")
    for f in FAILURES:
        print("  -", f)
    sys.exit(0 if not FAILURES else 1)


if __name__ == "__main__":
    main()
