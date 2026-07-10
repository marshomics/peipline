#!/usr/bin/env python3
"""Static checks on the Snakefile and the conda environments.

Not a substitute for `snakemake -n`, which needs snakemake installed and the
real input files present. This runs anywhere and catches the wiring mistakes
that a dry-run would catch: a rule pointing at a script that does not exist, a
shell command that never activates an environment, an env YAML that would let
conda fall back to the `defaults` channel.

    python scripts/check_snakefile.py
"""
from __future__ import annotations

import os
import re
import sys

import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ENVS = ["hmmer", "py", "phylo", "network", "prodigal", "r", "selection"]
BANNED_CHANNELS = {"defaults", "anaconda", "main", "r", "free", "pro", "msys2"}

problems: list[str] = []
notes: list[str] = []


def bad(msg):
    problems.append(msg)


def parse_rules(text):
    """Yield (name, body) for every rule/checkpoint block."""
    starts = [(m.start(), m.group(1))
              for m in re.finditer(r"^(?:checkpoint|rule)\s+(\w+):", text, re.M)]
    for i, (pos, name) in enumerate(starts):
        end = starts[i + 1][0] if i + 1 < len(starts) else len(text)
        yield name, text[pos:end]


# Walltime caps from sge_probe_report. standard.q's soft limit lands first.
QUEUE_CAP_SEC = {"standard.q": 23 * 3600 + 55 * 60,
                 "long.q": 670 * 3600,
                 "cryo-em.q": 670 * 3600,
                 "test.q": 10 ** 9}
SMALLEST_NODE_CORES = 32


def _sec(hms):
    h, m, s = (list(map(int, str(hms).split(":"))) + [0, 0])[:3]
    return h * 3600 + m * 60 + s


def check_sge_profile(prof, snakefile_text):
    """The profile is where a wrong number costs 700 pending jobs."""
    sub = str(prof.get("cluster-generic-submit-cmd", ""))
    if "submit.sh" not in sub:
        bad("profiles/sge: submit-cmd is a bare qsub. h_vmem is CONSUMABLE on "
            "this cluster, so it is requested per slot; a bare `-l h_vmem={mem_mb}M` "
            "over-requests by a factor of `threads`. Use profiles/sge/submit.sh.")
    if "-l h_vmem" in sub or "-pe " in sub:
        bad("profiles/sge: submit-cmd still passes -l h_vmem / -pe directly; "
            "submit.sh must own those (per-slot arithmetic, no -pe when threads=1)")
    if "openmpi" in sub:
        bad("profiles/sge: openmpi is a $round_robin PE and scatters slots across "
            "hosts; every threaded tool here needs a single host. Use `parallel`.")

    threads = prof.get("set-threads") or {}
    res = prof.get("set-resources") or {}
    for rule, t in threads.items():
        if int(t) > SMALLEST_NODE_CORES:
            notes.append(f"set-threads {rule}={t} exceeds the smallest node "
                         f"({SMALLEST_NODE_CORES} cores); submit.sh will cap it")

    rules = {m.group(1) for m in
             re.finditer(r"^(?:checkpoint|rule)\s+(\w+):", snakefile_text, re.M)} - {"all"}

    for rule, r in res.items():
        if rule not in rules:
            bad(f"profiles/sge set-resources names '{rule}', which is not a rule "
                f"in the Snakefile")
        q = str(r.get("queue", "standard.q")).strip("'\"")
        if q not in QUEUE_CAP_SEC:
            bad(f"profiles/sge: rule {rule} uses unknown queue '{q}'")
            continue
        if "h_rt" not in r:
            bad(f"profiles/sge: rule {rule} has no h_rt")
            continue
        if _sec(r["h_rt"]) > QUEUE_CAP_SEC[q]:
            bad(f"profiles/sge: rule {rule} wants h_rt={r['h_rt']} on {q}, which "
                f"caps at {QUEUE_CAP_SEC[q] // 3600}h. It will never be scheduled.")
        if "mem_mb" not in r:
            bad(f"profiles/sge: rule {rule} has no mem_mb")

    for rule in sorted(rules):
        if rule not in res:
            notes.append(f"rule {rule} has no set-resources entry; it gets "
                         f"default-resources")

    dr = prof.get("default-resources") or []
    joined = " ".join(str(x) for x in dr)
    for k in ("mem_mb", "h_rt", "queue"):
        if k not in joined:
            bad(f"profiles/sge default-resources is missing {k}; submit.sh needs it")

    notes.append(f"SGE profile: {len(res)} rules with explicit resources, "
                 f"{sum(1 for r in res.values() if 'queue' in r)} routed off standard.q")


def main():
    sf = os.path.join(ROOT, "Snakefile")
    text = open(sf).read()

    # --- env YAMLs ----------------------------------------------------------
    for e in ENVS:
        p = os.path.join(ROOT, "envs", f"{e}.yaml")
        if not os.path.exists(p):
            bad(f"envs/{e}.yaml missing")
            continue
        y = yaml.safe_load(open(p))
        ch = [str(c) for c in (y.get("channels") or [])]
        if "nodefaults" not in ch:
            bad(f"envs/{e}.yaml has no `nodefaults` channel: conda will silently "
                f"append `defaults`")
        for c in ch:
            if c.lower() in BANNED_CHANNELS:
                bad(f"envs/{e}.yaml lists banned channel '{c}'")
        if not y.get("dependencies"):
            bad(f"envs/{e}.yaml has no dependencies")
    notes.append(f"{len(ENVS)} env YAMLs checked")

    # --- envs_root plumbing -------------------------------------------------
    cfg = yaml.safe_load(open(os.path.join(ROOT, "config.yaml")))
    if "envs_root" not in cfg:
        bad("config.yaml has no `envs_root` key")
    for e in ENVS:
        if f'_act("{e}")' not in text:
            bad(f"Snakefile never activates the '{e}' environment")

    prof = yaml.safe_load(open(os.path.join(ROOT, "profiles", "sge", "config.yaml")))
    if "software-deployment-method" in prof:
        bad("profiles/sge/config.yaml sets software-deployment-method; it would "
            "override run.sh and make offline jobs try to download")

    check_sge_profile(prof, text)

    # --- per-rule checks ----------------------------------------------------
    env_vars = {"ENV_HMMER", "ENV_PY", "ENV_PHYLO", "ENV_NET", "ENV_PROD", "ENV_R",
                "ENV_SEL"}
    conda_to_var = {"hmmer": "ENV_HMMER", "py": "ENV_PY", "phylo": "ENV_PHYLO",
                    "network": "ENV_NET", "prodigal": "ENV_PROD", "r": "ENV_R",
                    "selection": "ENV_SEL"}
    n_rules = n_shell = 0

    for name, body in parse_rules(text):
        n_rules += 1
        if name == "all":
            continue
        has_shell = re.search(r"^\s*shell:", body, re.M)
        inherits = f"use rule" in body
        cm = re.search(r'conda:\s*"envs/(\w+)\.yaml"', body)

        if cm:
            envname = cm.group(1)
            if envname not in ENVS:
                bad(f"rule {name}: conda env '{envname}' is not one of {ENVS}")
            if not os.path.exists(os.path.join(ROOT, "envs", f"{envname}.yaml")):
                bad(f"rule {name}: envs/{envname}.yaml does not exist")
        elif has_shell and not inherits:
            bad(f"rule {name}: has a shell but no conda: directive")

        if has_shell:
            n_shell += 1
            want = conda_to_var.get(cm.group(1)) if cm else None
            used = set(re.findall(r"\{(ENV_\w+)\}", body))
            # search_batch delegates to a wrapper script that activates several
            # environments in sequence; it must be handed ENVS_ROOT instead.
            delegates = "ENVS_ROOT=" in body
            if not used and not delegates:
                bad(f"rule {name}: shell never activates an environment "
                    f"(no {{ENV_*}} token, no ENVS_ROOT passed); it will run with "
                    f"whatever is on PATH")
            elif used and want and want not in used:
                bad(f"rule {name}: conda env is '{cm.group(1)}' but shell activates "
                    f"{sorted(used)}, expected {{{want}}}")
            for u in used - env_vars:
                bad(f"rule {name}: shell references undefined {{{u}}}")
            for m in re.finditer(r"bash \{SCRIPTS\}/([\w.]+)", body):
                if not delegates and not used:
                    bad(f"rule {name} runs bash {m.group(1)} with neither an "
                        f"{{ENV_*}} prefix nor ENVS_ROOT")

        # scripts referenced must exist
        for m in re.finditer(r"\{SCRIPTS\}/([\w.]+)", body):
            p = os.path.join(ROOT, "scripts", m.group(1))
            if not os.path.exists(p):
                bad(f"rule {name}: scripts/{m.group(1)} does not exist")

    notes.append(f"{n_rules} rules, {n_shell} with a shell block")

    # --- setup scripts ------------------------------------------------------
    for s in ("setup/build_envs.sh", "setup/install_envs.sh"):
        p = os.path.join(ROOT, s)
        if not os.path.exists(p):
            bad(f"{s} missing")
        elif not os.access(p, os.X_OK):
            notes.append(f"{s} is not executable (chmod +x)")

    # --- report -------------------------------------------------------------
    for n in notes:
        print(f"[note]  {n}")
    if problems:
        print()
        for p in problems:
            print(f"[FAIL]  {p}")
        print(f"\n{len(problems)} problem(s)")
        sys.exit(1)
    print("\nSnakefile wiring OK")


if __name__ == "__main__":
    main()
