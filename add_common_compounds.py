#!/usr/bin/env python3
"""Add missing common compounds to seed.toml with literature values."""

import os

SEED = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    'crates', 'tpt-thermo-data', 'data', 'seed.toml')

with open(SEED) as f:
    text = f.read()

existing = set()
for line in text.splitlines():
    if line.startswith('name = '):
        existing.add(line.split('"')[1])

# Missing compounds with NIST / Poling et al. 2001 literature values.
# (name, formula, cas, Tc_K, Pa, omega, M_kg_mol, Tb_K, source)
new_compounds = [
    ("ethyl acetate",    "C4H8O2",  "141-78-6",  523.30,  3_880_000.0, 0.3664, 0.0881051, 350.21, "NIST / Poling et al. 2001"),
    ("pyridine",         "C5H5N",   "110-86-1",  619.95,  5_670_000.0, 0.2400, 0.0790980, 388.40, "NIST / Poling et al. 2001"),
    ("xylene",           "C8H10",   "1330-20-7", 617.05,  3_541_000.0, 0.3100, 0.1061650, 411.50, "NIST / Poling et al. 2001 (mixed xylene)"),
    ("dichloromethane",  "CH2Cl2",  "75-09-2",   508.00,  6_190_000.0, 0.1980, 0.0849300, 312.90, "NIST / Poling et al. 2001"),
    ("chlorodifluoromethane", "CHClF2", "75-45-6", 369.30, 4_990_000.0, 0.2200, 0.0864700, 232.30, "NIST / ASHRAE"),
    ("thiophene",        "C4H4S",   "110-02-1",  579.40,  5_620_000.0, 0.2000, 0.0841400, 357.30, "NIST / Poling et al. 2001"),
    ("n-eicosane",       "C20H42",  "112-95-8",  767.00,  1_110_000.0, 0.9070, 0.2825476, 617.00, "NIST / Poling et al. 2001"),
    ("pyrene",           "C16H10",  "129-00-0",  936.00,  2_610_000.0, 0.5300, 0.2022534, 668.00, "NIST / Poling et al. 2001"),
    ("benzoic acid",     "C7H6O2",  "65-85-0",   751.00,  4_470_000.0, 0.6200, 0.1221210, 522.00, "NIST / Poling et al. 2001"),
    ("bromobenzene",     "C6H5Br",  "108-86-2",  670.00,  4_520_000.0, 0.2510, 0.1570080, 429.15, "NIST / Poling et al. 2001"),
    ("r245fa",           "C3H3F5",  "460-73-1",  427.20,  3_651_000.0, 0.3780, 0.1340483, 288.20, "NIST / ASHRAE"),
]

added = []
blocks = []
for name, formula, cas, Tc, Pc, omega, M, Tb, src in new_compounds:
    if name in existing:
        continue
    block = f"""[[components]]
schema_version = 1
name = "{name}"
formula = "{formula}"
cas = "{cas}"
critical_temperature_k = {Tc:.2f}
critical_pressure_pa = {Pc:.1f}
acentric_factor = {omega:.4f}
molar_mass_kg_per_mol = {M:.7f}
normal_boiling_point_k = {Tb:.2f}
source = "{src}"
"""
    blocks.append(block)
    added.append(name)

if not blocks:
    print("All compounds already present.")
else:
    # Insert before the binary interactions section.
    bip_idx = text.find('[[binary_interactions]]')
    if bip_idx > 0:
        out = text[:bip_idx].rstrip() + '\n\n' + ''.join(blocks) + '\n' + text[bip_idx:]
    else:
        out = text.rstrip() + '\n\n' + ''.join(blocks)
    with open(SEED, 'w') as f:
        f.write(out)
    print(f"Added {len(added)} compounds: {', '.join(added)}")
    print(f"Total now: {len(existing) + len(added)}")
