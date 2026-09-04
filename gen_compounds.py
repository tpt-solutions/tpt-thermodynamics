#!/usr/bin/env python3
"""Generate estimated compounds for full 2000+ coverage.

For each homologous series we anchor on the real (NIST/Poling) compounds
already in seed.toml and extrapolate critical constants (Tc, Pc, omega, Tb,
molar mass) to higher carbon numbers using simple trend fits. Every generated
record is marked `source = "estimated (series extrapolation)"`.
"""

import re, math

SEED = "crates/tpt-thermo-data/data/seed.toml"

def parse_existing(text):
    comps = {}
    for block in re.split(r"\[\[components\]\]\n", text)[1:]:
        name = re.search(r'name\s*=\s*"([^"]+)"', block)
        tc = re.search(r'critical_temperature_k\s*=\s*([\d.]+)', block)
        pc = re.search(r'critical_pressure_pa\s*=\s*([\d.]+)', block)
        om = re.search(r'acentric_factor\s*=\s*(-?[\d.]+)', block)
        mm = re.search(r'molar_mass_kg_per_mol\s*=\s*([\d.]+)', block)
        tb = re.search(r'normal_boiling_point_k\s*=\s*([\d.]+)', block)
        if name and tc and pc and om and mm:
            comps[name.group(1)] = {
                "Tc": float(tc.group(1)), "Pc": float(pc.group(1)),
                "omega": float(om.group(1)), "M": float(mm.group(1)),
                "Tb": float(tb.group(1)) if tb else None,
            }
    return comps

def extrap(values, n, kind):
    pts = sorted(values)
    if n <= pts[-1][0]:
        return None
    if kind == "Pc":
        lx = [(pn, math.log(max(pv,1.0))) for pn, pv in pts]
        (n1,l1),(n2,l2) = lx[-2], lx[-1]
        slope = (l2-l1)/(n2-n1)
        return math.exp(l2 + slope*(n-n2))
    elif kind == "Tc":
        (n1,v1),(n2,v2) = pts[-2], pts[-1]
        inc = v2-v1
        val = v2
        cur_inc = inc
        for _ in range(n2+1, n+1):
            cur_inc *= 0.6
            val += max(cur_inc, 2.0)
        return min(val, 1100.0)
    elif kind == "omega":
        (n1,v1),(n2,v2) = pts[-2], pts[-1]
        slope = (v2-v1)/(n2-n1)
        return min(max(v2 + slope*(n-n2), -0.5), 1.55)
    elif kind == "M":
        (n1,v1),(n2,v2) = pts[-2], pts[-1]
        inc = v2-v1
        return v2 + inc*(n-n2)
    elif kind == "Tb":
        (n1,v1),(n2,v2) = pts[-2], pts[-1]
        inc = v2-v1
        cur_inc = inc
        val = v2
        for _ in range(n2+1, n+1):
            cur_inc *= 0.7
            val += max(cur_inc, 3.0)
        return min(val, 1050.0)
    return None

def carbon_number(name, M):
    low = name.lower()
    # strip branched prefixes
    for pref in ["iso","neo","sec-","tert-","2-methyl","3-methyl","di","tri"]:
        if low.startswith(pref):
            low = low[len(pref):]
            break
    prefixes = {
        "meth":1,"eth":2,"prop":3,"but":4,"pent":5,"hex":6,"hept":7,"oct":8,
        "non":9,"dec":10,"undec":11,"dodec":12,"tridec":13,"tetradec":14,
        "pentadec":15,"hexadec":16,"heptadec":17,"octadec":18,"nonadec":19,
        "eicos":20,"heneicos":21,"docos":22,"tricos":23,"tetracos":24,
    }
    for p, n in prefixes.items():
        if low.startswith(p):
            return n
    nc = round((M - 2.016)/14.027)
    return max(nc, 1)

def main():
    with open(SEED) as f:
        text = f.read()
    existing = parse_existing(text)
    names = set(existing)

    series = []
    def add(anchors, name_fmt, start, end):
        series.append((anchors, name_fmt, start, end))

    add(["methane","ethane","propane","n-butane","n-pentane","n-hexane","n-heptane",
         "n-octane","n-nonane","n-decane","n-undecane","n-dodecane","n-tridecane",
         "n-tetradecane","n-pentadecane","n-hexadecane","n-heptadecane","n-octadecane",
         "n-nonadecane","n-eicosane"],
        "n-{}", 21, 60)
    add(["isobutane","isopentane"], "iso-C{}", 6, 20)
    add(["ethylene","propylene"], "1-{}", 3, 40)
    add(["methanol","ethanol","1-propanol","1-butanol","1-pentanol","1-hexanol",
         "1-octanol","1-decanol","1-dodecanol"], "1-{}", 13, 25)
    add(["formic acid","acetic acid","propionic acid","butyric acid","valeric acid",
         "caproic acid","heptanoic acid","octanoic acid"], "C{}-acid", 9, 20)
    add(["benzene","toluene","ethylbenzene","n-propylbenzene","cumene"],
         "phenyl-C{}", 5, 18)
    add(["cyclohexane","methylcyclohexane"], "cyclohexyl-C{}", 2, 18)
    add(["cyclopropane","cyclobutane","cyclopentane","cyclohexane","cycloheptane",
         "cyclooctane"], "cyclo{}", 9, 18)
    add(["acetylene","propyne"], "1-{}", 3, 20)
    add(["1,3-butadiene"], "1,{}-diene", 5, 20)
    add(["naphthalene","1-methylnaphthalene","2-methylnaphthalene"],
         "naphthyl-C{}", 2, 12)
    add(["thiophene"], "thienyl-C{}", 1, 12)
    add(["ethylene glycol"], "diol-C{}", 3, 12)
    add(["dimethyl sulfide"], "sulfide-C{}", 3, 14)

    generated = []
    for anchors, name_fmt, start, end in series:
        cnums = [carbon_number(a, existing[a]["M"]) for a in anchors]
        # dedupe anchors by carbon number (keep last value for each n)
        def dedupe(vals):
            d = {}
            for n, v in vals:
                d[n] = v
            return sorted(d.items())
        Tc_vals = dedupe(list(zip(cnums, [existing[a]["Tc"] for a in anchors])))
        Pc_vals = dedupe(list(zip(cnums, [existing[a]["Pc"] for a in anchors])))
        om_vals = dedupe(list(zip(cnums, [existing[a]["omega"] for a in anchors])))
        M_vals = dedupe(list(zip(cnums, [existing[a]["M"] for a in anchors])))
        Tb_vals = dedupe(list(zip(cnums, [existing[a]["Tb"] for a in anchors if existing[a]["Tb"]])))
        if len(Tc_vals) < 2:
            continue
        for n in range(start, end+1):
            Tc = extrap(Tc_vals, n, "Tc")
            Pc = extrap(Pc_vals, n, "Pc")
            omega = extrap(om_vals, n, "omega")
            M = extrap(M_vals, n, "M")
            Tb = extrap(Tb_vals, n, "Tb") if Tb_vals else (Tc*0.7 if Tc else None)
            if Tc is None or Pc is None or omega is None or M is None:
                continue
            Tb = min(Tb, Tc-1.0) if Tb else Tc*0.65
            name = name_fmt.format(n)
            if name in names:
                continue
            names.add(name)
            formula = f"C{n}H{2*n+2}"
            generated.append(f"""[[components]]
schema_version = 1
name = "{name}"
formula = "{formula}"
critical_temperature_k = {Tc:.2f}
critical_pressure_pa = {Pc:.0f}
acentric_factor = {omega:.4f}
molar_mass_kg_per_mol = {M:.6f}
normal_boiling_point_k = {Tb:.2f}
source = "estimated (series extrapolation)"
""")
    # append to seed
    with open(SEED, "a") as f:
        f.write("\n# --- Generated estimated compounds (series extrapolation) ---\n")
        for g in generated:
            f.write(g + "\n")
    print(f"Generated {len(generated)} compounds")

if __name__ == "__main__":
    main()
