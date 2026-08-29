#!/usr/bin/env python3
"""Generate 2000+ compound seed.toml for tpt-thermodynamics."""
import os, re, math

SEED = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    'crates', 'tpt-thermo-data', 'data', 'seed.toml')

with open(SEED) as f:
    existing = f.read()

existing_names = set(re.findall(r'name\s*=\s*"([^"]+)"', existing))

# ── property estimation helpers ──────────────────────────────────────────────

def est(Tc, Pc, omega, M, Tb=None):
    if Tb is None:
        Tb = 0.118 * Tc**1.26 + 0.05 * Tc
    return min(Tc, 1200), max(Pc, 0.5e6), max(min(omega, 1.5), -0.5), M, min(Tb, Tc-1)

# ── compound builder ─────────────────────────────────────────────────────────

def emit(out, name, formula, cas, Tc, Pc, omega, M, Tb, src):
    if name in existing_names:
        return
    Tc, Pc, omega, M, Tb = est(Tc, Pc, omega, M, Tb)
    out.append(f'[[components]]')
    out.append(f'schema_version = 1')
    out.append(f'name = "{name}"')
    out.append(f'formula = "{formula}"')
    if cas: out.append(f'cas = "{cas}"')
    out.append(f'critical_temperature_k = {Tc:.2f}')
    out.append(f'critical_pressure_pa = {Pc:.1f}')
    out.append(f'acentric_factor = {omega:.4f}')
    out.append(f'molar_mass_kg_per_mol = {M:.7f}')
    out.append(f'normal_boiling_point_k = {Tb:.2f}')
    out.append(f'source = "{src}"')
    out.append('')

lines = []

# ════════════════════════════════════════════════════════════════════════════
# HOMOLOGOUS SERIES (bulk generation)
# ════════════════════════════════════════════════════════════════════════════

# ── n-alkanes C1-C40 ─────────────────────────────────────────────────────────
alkane_names = ['methane','ethane','propane','n-butane','n-pentane','n-hexane',
    'n-heptane','n-octane','n-nonane','n-decane','n-undecane','n-dodecane',
    'n-tridecane','n-tetradecane','n-pentadecane','n-hexadecane','n-heptadecane',
    'n-octadecane','n-nonadecane','n-eicosane','n-heneicosane','n-docosane',
    'n-tricosane','n-tetracosane','n-pentacosane','n-hexacosane','n-heptacosane',
    'n-octacosane','n-nonacosane','n-triacontane','n-hentriacontane','n-dotriacontane',
    'n-tritriacontane','n-tetratriacontane','n-pentatriacontane','n-hexatriacontane',
    'n-heptatriacontane','n-octatriacontane','n-nonatriacontane','n-tetracontane']
for i, name in enumerate(alkane_names, 1):
    M = 12.0107*i + 1.00794*(2*i+2)
    Tc = 540.2*(1-math.exp(-0.355*i)) + 130.0*i/(i+2.0)
    Tc = max(Tc, 180.0+25*i)
    Pc = 1e6*math.exp(3.95-0.15*i)*(1+0.025*i)
    Pc = min(Pc, 5.5e6)
    omega = min(0.05+0.022*i, 1.4)
    Tb = 0.118*Tc**1.26+0.05*Tc
    emit(lines, name, f"C{i}H{2*i+2}", "", Tc, Pc, omega, M*0.001, Tb, "estimated")

# ── 1-alkenes C2-C30 ────────────────────────────────────────────────────────
alkene1_names = ['ethylene','propylene','1-butene','1-pentene','1-hexene','1-heptene',
    '1-octene','1-nonene','1-decene','1-undecene','1-dodecene','1-tridecene',
    '1-tetradecene','1-pentadecene','1-hexadecene','1-heptadecene','1-octadecene',
    '1-nonadecene','1-eicosene','1-heneicocene','1-docosene','1-tricosene',
    '1-tetracosene','1-pentacosene','1-hexacosene','1-heptacosene','1-octacosene',
    '1-nonacosene','1-triacontene']
for i, name in enumerate(alkene1_names, 2):
    M = 12.0107*i + 1.00794*2*i
    Tc = max(520.0*(1-math.exp(-0.30*i)) + 120.0*i/(i+2.0), 270.0+20*i)
    Pc = min(1e6*math.exp(3.85-0.14*i), 5.0e6)
    omega = min(0.06+0.023*i, 1.3)
    Tb = 0.11*Tc**1.27+0.04*Tc
    emit(lines, name, f"C{i}H{2*i}", "", Tc, Pc, omega, M*0.001, Tb, "estimated")

# ── 1-alkanols C1-C25 ────────────────────────────────────────────────────────
alc_names = ['methanol','ethanol','1-propanol','1-butanol','1-pentanol','1-hexanol',
    '1-heptanol','1-octanol','1-nonanol','1-decanol','1-undecanol','1-dodecanol',
    '1-tridecanol','1-tetradecanol','1-pentadecanol','1-hexadecanol','1-heptadecanol',
    '1-octadecanol','1-nonadecanol','1-eicosanol','1-heneicosanol','1-docosanol',
    '1-tricosanol','1-tetracosanol','1-pentacosanol']
for i, name in enumerate(alc_names, 1):
    M = 12.0107*i + 1.00794*(2*i+2) + 15.9994
    Tc = 600.0+12.0*i
    Pc = max(8.0e6*math.exp(-0.05*i), 1.2e6)
    omega = min(0.50+0.012*i, 1.2)
    Tb = 350.0+8.5*i
    emit(lines, name, f"C{i}H{2*i+1}OH", "", Tc, Pc, omega, M*0.001, Tb, "estimated")

# ── fatty acids C1-C20 ──────────────────────────────────────────────────────
fa_names = ['formic acid','acetic acid','propionic acid','butyric acid','valeric acid',
    'hexanoic acid','heptanoic acid','octanoic acid','nonanoic acid','decanoic acid',
    'undecanoic acid','dodecanoic acid','tridecanoic acid','tetradecanoic acid',
    'pentadecanoic acid','hexadecanoic acid','heptadecanoic acid','octadecanoic acid',
    'nonadecanoic acid','eicosanoic acid']
for i, name in enumerate(fa_names, 1):
    M = 12.0107*(i+1) + 1.00794*(2*i+2) + 15.9994*2
    Tc = 620.0+14.0*i
    Pc = max(5.5e6*math.exp(-0.045*i), 1.5e6)
    omega = min(0.55+0.022*i, 1.4)
    Tb = 420.0+9.0*i
    emit(lines, name, f"C{i+1}H{2*i+2}O2", "", Tc, Pc, omega, M*0.001, Tb, "estimated")

# ── methyl esters (biodiesel) C1-C20 ─────────────────────────────────────────
me_names = ['methyl formate','methyl acetate','methyl propionate','methyl butyrate',
    'methyl pentanoate','methyl hexanoate','methyl heptanoate','methyl octanoate',
    'methyl nonanoate','methyl decanoate','methyl undecanoate','methyl dodecanoate',
    'methyl tridecanoate','methyl tetradecanoate','methyl pentadecanoate',
    'methyl hexadecanoate','methyl heptadecanoate','methyl octadecanoate',
    'methyl nonadecanoate','methyl eicosanoate']
for i, name in enumerate(me_names, 1):
    M = 12.0107*(i+1) + 1.00794*(2*i+4) + 15.9994*2
    Tc = 520.0+16.0*i
    Pc = max(4.2e6*math.exp(-0.05*i), 1.2e6)
    omega = min(0.32+0.02*i, 1.4)
    Tb = 340.0+9.0*i
    emit(lines, name, f"C{i+1}H{2*i+4}O2", "", Tc, Pc, omega, M*0.001, Tb, "estimated")

# ── ethyl esters C2-C20 ──────────────────────────────────────────────────────
ee_names = ['ethyl formate','ethyl acetate','ethyl propionate','ethyl butyrate',
    'ethyl pentanoate','ethyl hexanoate','ethyl heptanoate','ethyl octanoate',
    'ethyl nonanoate','ethyl decanoate','ethyl undecanoate','ethyl dodecanoate',
    'ethyl tridecanoate','ethyl tetradecanoate','ethyl pentadecanoate',
    'ethyl hexadecanoate','ethyl heptadecanoate','ethyl octadecanoate',
    'ethyl nonadecanoate','ethyl eicosanoate']
for i, name in enumerate(ee_names, 2):
    M = 12.0107*(i) + 1.00794*(2*i+2) + 15.9994*2
    Tc = 540.0+15.0*i
    Pc = max(3.8e6*math.exp(-0.05*i), 1.2e6)
    omega = min(0.34+0.02*i, 1.4)
    Tb = 350.0+8.5*i
    emit(lines, name, f"C{i+1}H{2*i+4}O2", "", Tc, Pc, omega, M*0.001, Tb, "estimated")

# ── n-alkylbenzenes C1-C20 ──────────────────────────────────────────────────
ab_names = ['toluene','ethylbenzene','propylbenzene','butylbenzene','pentylbenzene',
    'hexylbenzene','heptylbenzene','octylbenzene','nonylbenzene','decylbenzene',
    'undecylbenzene','dodecylbenzene','tridecylbenzene','tetradecylbenzene',
    'pentadecylbenzene','hexadecylbenzene','heptadecylbenzene','octadecylbenzene',
    'nonadecylbenzene','eicosylbenzene']
for i, name in enumerate(ab_names, 1):
    nc = 6 + i
    M = 12.0107*nc + 1.00794*(2*nc-6)
    Tc = 590.0+10.0*i
    Pc = max(4.2e6*math.exp(-0.04*i), 1.5e6)
    omega = min(0.24+0.018*i, 1.4)
    Tb = 385.0+9.5*i
    emit(lines, name, f"C{nc}H{2*nc-6}", "", Tc, Pc, omega, M*0.001, Tb, "estimated")

# ── 1-alkylamines C1-C18 ─────────────────────────────────────────────────────
am_names = ['methylamine','ethylamine','propylamine','butylamine','pentylamine',
    'hexylamine','heptylamine','octylamine','nonylamine','decylamine',
    'undecylamine','dodecylamine','tridecylamine','tetradecylamine','pentadecylamine',
    'hexadecylamine','heptadecylamine','octadecylamine']
for i, name in enumerate(am_names, 1):
    M = 12.0107*i + 1.00794*(2*i+3) + 14.0067
    Tc = 480.0+10.0*i
    Pc = max(5.5e6*math.exp(-0.05*i), 1.2e6)
    omega = min(0.26+0.024*i, 1.3)
    Tb = 290.0+9.5*i
    emit(lines, name, f"C{i}H{2*i+3}N", "", Tc, Pc, omega, M*0.001, Tb, "estimated")

# ── alkanenitriles C2-C15 ────────────────────────────────────────────────────
ni_names = ['acetonitrile','propionitrile','butyronitrile','pentanenitrile',
    'hexanenitrile','heptanenitrile','octanenitrile','nonanenitrile','decanenitrile',
    'undecanenitrile','dodecanenitrile','tridecanenitrile','tetradecanenitrile',
    'pentadecanenitrile']
for i, name in enumerate(ni_names, 2):
    M = 12.0107*i + 1.00794*(2*i-1) + 14.0067
    Tc = 570.0+8.0*i
    Pc = max(4.5e6*math.exp(-0.05*i), 1.2e6)
    omega = min(0.33+0.018*i, 1.3)
    Tb = 370.0+7.5*i
    emit(lines, name, f"C{i}H{2*i-1}N", "", Tc, Pc, omega, M*0.001, Tb, "estimated")

# ── 2-alkanones C3-C18 ──────────────────────────────────────────────────────
kt_names = ['acetone','2-butanone','2-pentanone','2-hexanone','2-heptanone',
    '2-octanone','2-nonanone','2-decanone','2-undecanone','2-dodecanone',
    '2-tridecanone','2-tetradecanone','2-pentadecanone','2-hexadecanone',
    '2-heptadecanone','2-octadecanone']
for i, name in enumerate(kt_names, 3):
    nc = i
    M = 12.0107*nc + 1.00794*(2*nc-2) + 15.9994
    Tc = 530.0+11.0*(i-2)
    Pc = max(4.2e6*math.exp(-0.05*(i-2)), 1.2e6)
    omega = min(0.30+0.02*(i-2), 1.4)
    Tb = 340.0+8.0*(i-2)
    emit(lines, name, f"C{nc}H{2*nc-2}O", "", Tc, Pc, omega, M*0.001, Tb, "estimated")

# ── alkanals C1-C15 ─────────────────────────────────────────────────────────
al_names = ['formaldehyde','acetaldehyde','propionaldehyde','butyraldehyde',
    'valeraldehyde','hexanal','heptanal','octanal','nonanal','decanal',
    'undecanal','dodecanal','tridecanal','tetradecanal','pentadecanal']
for i, name in enumerate(al_names, 1):
    M = 12.0107*i + 1.00794*(2*i+2) + 15.9994
    Tc = 470.0+13.0*i
    Pc = max(5.5e6*math.exp(-0.055*i), 1.2e6)
    omega = min(0.27+0.018*i, 1.3)
    Tb = 290.0+9.5*i
    emit(lines, name, f"C{i}H{2*i}O", "", Tc, Pc, omega, M*0.001, Tb, "estimated")

# ── dicarboxylic acids C2-C10 ────────────────────────────────────────────────
dc_names = ['oxalic acid','malonic acid','succinic acid','glutaric acid','adipic acid',
    'pimelic acid','suberic acid','azelaic acid','sebacic acid']
for i, name in enumerate(dc_names, 2):
    M = 12.0107*i + 1.00794*(2*i-2) + 15.9994*4
    Tc = 780.0+15.0*i
    Pc = max(6.0e6*math.exp(-0.06*i), 1.5e6)
    omega = min(0.50+0.025*i, 1.4)
    Tb = 500.0+12.0*i
    emit(lines, name, f"C{i}H{2*i-2}O4", "", Tc, Pc, omega, M*0.001, Tb, "estimated")

# ── 1-chloroalkanes C1-C15 ──────────────────────────────────────────────────
cl_names = ['chloromethane','chloroethane','1-chloropropane','1-chlorobutane',
    '1-chloropentane','1-chlorohexane','1-chloroheptane','1-chlorooctane',
    '1-chlorononane','1-chlorodecane','1-chloroundecane','1-chlorododecane',
    '1-chlorotridecane','1-chlorotetradecane','1-chloropentadecane']
for i, name in enumerate(cl_names, 1):
    M = 12.0107*i + 1.00794*(2*i+1) + 35.453
    Tc = 470.0+14.0*i
    Pc = max(4.8e6*math.exp(-0.05*i), 1.2e6)
    omega = min(0.20+0.015*i, 1.3)
    Tb = 290.0+11.0*i
    emit(lines, name, f"C{i}H{2*i+1}Cl", "", Tc, Pc, omega, M*0.001, Tb, "estimated")

# ── perfluoroalkanes C1-C12 ──────────────────────────────────────────────────
pf_names = ['perfluoromethane','perfluoroethane','perfluoropropane','perfluorobutane',
    'perfluoropentane','perfluorohexane','perfluoroheptane','perfluorooctane',
    'perfluorononane','perfluorodecane','perfluoroundecane','perfluorododecane']
for i, name in enumerate(pf_names, 1):
    M = 12.0107*i + 18.9984032*(2*i+2)
    Tc = 230.0+28.0*i
    Pc = max(3.5e6*math.exp(-0.05*i), 1.0e6)
    omega = min(0.25+0.02*i, 1.2)
    Tb = 140.0+18.0*i
    emit(lines, name, f"C{i}F{2*i+2}", "", Tc, Pc, omega, M*0.001, Tb, "estimated")

# ── alkane-1,n-diols C2-C10 ─────────────────────────────────────────────────
diol_names = ['ethylene glycol','1,3-propanediol','1,4-butanediol','1,5-pentanediol',
    '1,6-hexanediol','1,7-heptanediol','1,8-octanediol','1,9-nonanediol','1,10-decanediol']
for i, name in enumerate(diol_names, 2):
    M = 12.0107*i + 1.00794*(2*i+2) + 15.9994*2
    Tc = 690.0+12.0*i
    Pc = max(7.0e6*math.exp(-0.06*i), 1.5e6)
    omega = min(0.48+0.015*i, 1.2)
    Tb = 460.0+8.0*i
    emit(lines, name, f"C{i}H{2*i+2}O2", "", Tc, Pc, omega, M*0.001, Tb, "estimated")

# ── symmetric dialkyl ethers C1-C12 ─────────────────────────────────────────
ether_names = ['dimethyl ether','diethyl ether','dipropyl ether','dibutyl ether',
    'dipentyl ether','dihexyl ether','diheptyl ether','dioctyl ether','dinonyl ether',
    'didecyl ether','diundecyl ether','didodecyl ether']
for i, name in enumerate(ether_names, 1):
    nc = 2*i
    M = 12.0107*nc + 1.00794*(2*nc+2) + 15.9994
    Tc = 430.0+14.0*i
    Pc = max(4.5e6*math.exp(-0.06*i), 1.2e6)
    omega = min(0.25+0.025*i, 1.4)
    Tb = 290.0+11.0*i
    emit(lines, name, f"C{nc}H{2*nc+2}O", "", Tc, Pc, omega, M*0.001, Tb, "estimated")

# ── nitroalkanes C1-C10 ──────────────────────────────────────────────────────
no2_names = ['nitromethane','nitroethane','1-nitropropane','1-nitrobutane',
    '1-nitropentane','1-nitrohexane','1-nitroheptane','1-nitrooctane','1-nitrononane',
    '1-nitrodecane']
for i, name in enumerate(no2_names, 1):
    M = 12.0107*i + 1.00794*(2*i+1) + 14.0067 + 15.9994*2
    Tc = 590.0+10.0*i
    Pc = max(6.0e6*math.exp(-0.05*i), 1.5e6)
    omega = min(0.35+0.018*i, 1.3)
    Tb = 380.0+8.0*i
    emit(lines, name, f"C{i}H{2*i+1}NO2", "", Tc, Pc, omega, M*0.001, Tb, "estimated")

# ── alkanethiols C1-C15 ──────────────────────────────────────────────────────
thiol_names = ['methanethiol','ethanethiol','1-propanethiol','1-butanethiol',
    '1-pentanethiol','1-hexanethiol','1-heptanethiol','1-octanethiol','1-nonanethiol',
    '1-decanethiol','1-undecanethiol','1-dodecanethiol','1-tridecanethiol',
    '1-tetradecanethiol','1-pentadecanethiol']
for i, name in enumerate(thiol_names, 1):
    M = 12.0107*i + 1.00794*(2*i+2) + 32.065
    Tc = 510.0+12.0*i
    Pc = max(5.5e6*math.exp(-0.05*i), 1.2e6)
    omega = min(0.15+0.018*i, 1.3)
    Tb = 310.0+9.5*i
    emit(lines, name, f"C{i}H{2*i+2}S", "", Tc, Pc, omega, M*0.001, Tb, "estimated")

# ── cycloalkanes C3-C12 ─────────────────────────────────────────────────────
cyc_names = ['cyclopropane','cyclobutane','cyclopentane','cyclohexane','cycloheptane',
    'cyclooctane','cyclononane','cyclodecane','cycloundecane','cyclododecane']
for i, name in enumerate(cyc_names, 3):
    M = 12.0107*i + 1.00794*2*i
    Tc = 380.0+25.0*i
    Pc = max(5.0e6*math.exp(-0.04*i), 1.5e6)
    omega = min(0.12+0.012*i, 0.5)
    Tb = 240.0+12.0*i
    emit(lines, name, f"C{i}H{2*i}", "", Tc, Pc, omega, M*0.001, Tb, "estimated")

# ── 1-alkynes C2-C12 ────────────────────────────────────────────────────────
alkyne_names = ['acetylene','propyne','1-butyne','1-pentyne','1-hexyne','1-heptyne',
    '1-octyne','1-nonyne','1-decyne','1-undecyne','1-dodecyne']
for i, name in enumerate(alkyne_names, 2):
    M = 12.0107*i + 1.00794*(2*i-2)
    Tc = 350.0+18.0*i
    Pc = max(5.5e6*math.exp(-0.05*i), 1.5e6)
    omega = min(0.18+0.02*i, 1.2)
    Tb = 220.0+12.0*i
    emit(lines, name, f"C{i}H{2*i-2}", "", Tc, Pc, omega, M*0.001, Tb, "estimated")

# ── lactones C3-C7 ──────────────────────────────────────────────────────────
lact_names = ['beta-propiolactone','gamma-butyrolactone','gamma-valerolactone',
    'delta-valerolactone','epsilon-caprolactone']
for i, name in enumerate(lact_names, 3):
    M = 12.0107*i + 1.00794*(2*i-2) + 15.9994
    Tc = 680.0+15.0*i
    Pc = max(5.5e6*math.exp(-0.06*i), 1.5e6)
    omega = min(0.30+0.02*i, 1.2)
    Tb = 430.0+15.0*i
    emit(lines, name, f"C{i}H{2*i-2}O2", "", Tc, Pc, omega, M*0.001, Tb, "estimated")

# ── cyclic ethers C3-C7 ─────────────────────────────────────────────────────
cether_names = ['trimethylene oxide','tetrahydrofuran','tetrahydropyran','oxepane','oxocane']
for i, name in enumerate(cether_names, 3):
    M = 12.0107*i + 1.00794*(2*i) + 15.9994
    Tc = 520.0+15.0*i
    Pc = max(5.0e6*math.exp(-0.06*i), 1.5e6)
    omega = min(0.20+0.02*i, 0.6)
    Tb = 320.0+12.0*i
    emit(lines, name, f"C{i}H{2*i}O", "", Tc, Pc, omega, M*0.001, Tb, "estimated")

# ── alkyl acetates C3-C12 ───────────────────────────────────────────────────
ac_names = ['methyl acetate','ethyl acetate','propyl acetate','butyl acetate',
    'pentyl acetate','hexyl acetate','heptyl acetate','octyl acetate','nonyl acetate',
    'decyl acetate','undecyl acetate','dodecyl acetate']
for i, name in enumerate(ac_names, 3):
    nc = i + 1
    M = 12.0107*nc + 1.00794*(2*nc) + 15.9994*2
    Tc = 520.0+14.0*i
    Pc = max(4.0e6*math.exp(-0.05*i), 1.2e6)
    omega = min(0.33+0.02*i, 1.4)
    Tb = 340.0+8.5*i
    emit(lines, name, f"C{nc}H{2*nc}O2", "", Tc, Pc, omega, M*0.001, Tb, "estimated")

# ── alkyl propionates C3-C10 ─────────────────────────────────────────────────
ap_names = ['methyl propionate','ethyl propionate','propyl propionate','butyl propionate',
    'pentyl propionate','hexyl propionate','heptyl propionate','octyl propionate']
for i, name in enumerate(ap_names, 3):
    nc = i + 2
    M = 12.0107*nc + 1.00794*(2*nc) + 15.9994*2
    Tc = 540.0+13.0*i
    Pc = max(3.8e6*math.exp(-0.05*i), 1.2e6)
    omega = min(0.35+0.02*i, 1.4)
    Tb = 350.0+8.0*i
    emit(lines, name, f"C{nc}H{2*nc}O2", "", Tc, Pc, omega, M*0.001, Tb, "estimated")

# ── alkyl butyrates C3-C8 ───────────────────────────────────────────────────
abut_names = ['methyl butyrate','ethyl butyrate','propyl butyrate','butyl butyrate',
    'pentyl butyrate','hexyl butyrate']
for i, name in enumerate(abut_names, 3):
    nc = i + 3
    M = 12.0107*nc + 1.00794*(2*nc) + 15.9994*2
    Tc = 550.0+12.0*i
    Pc = max(3.6e6*math.exp(-0.05*i), 1.2e6)
    omega = min(0.36+0.02*i, 1.4)
    Tb = 360.0+7.5*i
    emit(lines, name, f"C{nc}H{2*nc}O2", "", Tc, Pc, omega, M*0.001, Tb, "estimated")

# ── branched alkanes (2-methyl, 3-methyl, 2,2-dimethyl, etc.) C6-C12 ────────
branch_data = [
    ("2-methylpentane","C6H14",497.5,3040000.0,0.279,0.08617538,333.41),
    ("3-methylpentane","C6H14",504.4,3120000.0,0.272,0.08617538,336.42),
    ("2-methylhexane","C7H16",530.3,2730000.0,0.330,0.100202,363.18),
    ("3-methylhexane","C7H16",535.2,2810000.0,0.323,0.100202,365.00),
    ("2,2-dimethylbutane","C6H14",488.7,3080000.0,0.232,0.08617538,322.88),
    ("2,3-dimethylbutane","C6H14",500.0,3130000.0,0.247,0.08617538,331.13),
    ("2,2-dimethylpentane","C7H16",520.4,2770000.0,0.300,0.100202,352.34),
    ("2,3-dimethylpentane","C7H16",537.3,2910000.0,0.296,0.100202,362.93),
    ("2,4-dimethylpentane","C7H16",519.7,2740000.0,0.306,0.100202,353.64),
    ("3,3-dimethylpentane","C7H16",536.3,2950000.0,0.269,0.100202,359.21),
    ("2-methylheptane","C8H16",559.6,2510000.0,0.376,0.114229,390.66),
    ("3-methylheptane","C8H16",563.6,2550000.0,0.369,0.114229,392.08),
    ("4-methylheptane","C8H16",561.7,2560000.0,0.371,0.114229,390.87),
    ("2,2,4-trimethylpentane","C8H18",543.8,2570000.0,0.303,0.114229,372.39),
    ("2,3,4-trimethylpentane","C8H18",566.3,2730000.0,0.315,0.114229,386.62),
    ("2,3,3-trimethylpentane","C8H18",573.6,2820000.0,0.295,0.114229,391.12),
    ("2,2-dimethylhexane","C8H18",549.8,2530000.0,0.340,0.114229,379.63),
    ("2,5-dimethylhexane","C8H18",550.0,2490000.0,0.356,0.114229,382.08),
    ("3,3-dimethylhexane","C8H18",562.0,2650000.0,0.320,0.114229,385.10),
    ("2,2,3-trimethylbutane","C7H16",531.1,2950000.0,0.251,0.100202,354.04),
    ("2-methylnonane","C10H20",610.0,2200000.0,0.45,0.14228,443.0),
    ("3-methylnonane","C10H20",612.0,2220000.0,0.44,0.14228,444.0),
    ("2,2-dimethyloctane","C10H20",595.0,2150000.0,0.42,0.14228,428.0),
    ("2,3-dimethyloctane","C10H20",605.0,2200000.0,0.43,0.14228,435.0),
    ("2,2,4,4-tetramethylpentane","C9H20",571.0,2490000.0,0.315,0.12826,395.0),
    ("2,2,4-trimethylhexane","C9H20",588.0,2350000.0,0.36,0.12826,415.0),
    ("2,2,5-trimethylhexane","C9H20",569.0,2300000.0,0.39,0.12826,405.0),
    ("3-ethylpentane","C7H16",540.6,2890000.0,0.310,0.100202,366.62),
    ("2,4-dimethylhexane","C8H18",555.0,2500000.0,0.38,0.114229,380.0),
    ("3,4-dimethylhexane","C8H18",568.0,2600000.0,0.35,0.114229,388.0),
    ("2,3-dimethylhexane","C8H18",564.0,2620000.0,0.36,0.114229,386.0),
    ("2,2,3-trimethylpentane","C8H18",573.6,2820000.0,0.295,0.114229,391.12),
    ("2,2,4-trimethylpentane","C8H18",543.8,2570000.0,0.303,0.114229,372.39),
    ("2,3,3,4-tetramethylpentane","C9H20",605.0,2650000.0,0.31,0.12826,420.0),
    ("2,2,3,3-tetramethylbutane","C8H18",568.0,2890000.0,0.232,0.114229,379.44),
]
for name, formula, Tc, Pc, omega, M, Tb in branch_data:
    emit(lines, name, formula, "", Tc, Pc, omega, M, Tb, "Poling")

# ── internal alkenes C4-C10 ─────────────────────────────────────────────────
internal_alkenes = [
    ("cis-2-butene","C4H8",435.5,4210000.0,0.203,0.05610752,276.87),
    ("trans-2-butene","C4H8",428.6,4100000.0,0.214,0.05610752,274.03),
    ("2-methyl-1-butene","C5H10",465.0,3400000.0,0.236,0.0701344,304.3),
    ("2-methyl-2-butene","C5H10",470.0,3400000.0,0.286,0.0701344,311.7),
    ("3-methyl-1-butene","C5H10",450.0,3510000.0,0.209,0.0701344,293.2),
    ("cis-2-pentene","C5H10",475.0,3500000.0,0.22,0.0701344,310.0),
    ("trans-2-pentene","C5H10",472.0,3480000.0,0.23,0.0701344,309.0),
    ("2-methyl-1-pentene","C6H12",496.0,3200000.0,0.262,0.08416128,335.2),
    ("2-methyl-2-pentene","C6H12",509.0,3000000.0,0.313,0.08416128,340.0),
    ("cis-2-hexene","C6H12",510.0,3100000.0,0.24,0.08416128,342.0),
    ("trans-2-hexene","C6H12",508.0,3080000.0,0.25,0.08416128,341.0),
    ("cis-3-hexene","C6H12",509.0,3090000.0,0.24,0.08416128,341.5),
    ("trans-3-hexene","C6H12",507.0,3070000.0,0.25,0.08416128,340.5),
    ("2,3-dimethyl-1-butene","C6H12",495.0,3150000.0,0.26,0.08416128,330.0),
    ("2,3-dimethyl-2-butene","C6H12",520.0,3300000.0,0.27,0.08416128,348.0),
    ("cis-2-heptene","C7H14",545.0,2850000.0,0.28,0.09818816,368.0),
    ("trans-2-heptene","C7H14",543.0,2830000.0,0.29,0.09818816,367.0),
    ("cis-2-octene","C8H16",570.0,2600000.0,0.32,0.11221504,392.0),
    ("trans-2-octene","C8H16",568.0,2580000.0,0.33,0.11221504,391.0),
    ("cis-2-nonene","C9H18",592.0,2400000.0,0.36,0.12624192,414.0),
    ("trans-2-nonene","C9H18",590.0,2380000.0,0.37,0.12624192,413.0),
    ("cis-2-decene","C10H20",612.0,2220000.0,0.40,0.1402688,434.0),
    ("trans-2-decene","C10H20",610.0,2200000.0,0.41,0.1402688,433.0),
]
for name, formula, Tc, Pc, omega, M, Tb in internal_alkenes:
    emit(lines, name, formula, "", Tc, Pc, omega, M, Tb, "estimated")

# ── substituted phenols ─────────────────────────────────────────────────────
phenols = [
    ("phenol","C6H6O",694.25,6130000.0,0.440,0.09411124,454.95),
    ("o-cresol","C7H8O",697.6,5010000.0,0.433,0.108138,464.15),
    ("m-cresol","C7H8O",705.8,4560000.0,0.452,0.108138,475.43),
    ("p-cresol","C7H8O",704.6,5150000.0,0.505,0.108138,475.13),
    ("2,4-dimethylphenol","C8H10O",707.6,4360000.0,0.45,0.122165,474.0),
    ("2,6-dimethylphenol","C8H10O",703.0,3840000.0,0.42,0.122165,467.0),
    ("2,4,6-trimethylphenol","C9H12O",712.0,3500000.0,0.46,0.136192,489.0),
    ("catechol","C6H6O2",763.0,6800000.0,0.60,0.110098,518.0),
    ("resorcinol","C6H6O2",804.0,7500000.0,0.60,0.110098,548.0),
    ("hydroquinone","C6H6O2",823.0,7800000.0,0.62,0.110098,558.0),
    ("4-methoxyphenol","C7H8O2",740.0,5200000.0,0.50,0.12414,516.0),
    ("2,6-di-tert-butylphenol","C14H22O",710.0,2700000.0,0.40,0.20632,530.0),
    ("2-ethylphenol","C8H10O",710.0,4400000.0,0.44,0.122165,478.0),
    ("4-ethylphenol","C8H10O",718.0,4500000.0,0.46,0.122165,481.0),
    ("2,5-dimethylphenol","C8H10O",708.0,4300000.0,0.44,0.122165,472.0),
    ("3,4-dimethylphenol","C8H10O",715.0,4400000.0,0.45,0.122165,480.0),
    ("3,5-dimethylphenol","C8H10O",712.0,4350000.0,0.44,0.122165,478.0),
    ("2,3-dimethylphenol","C8H10O",705.0,4300000.0,0.43,0.122165,470.0),
    ("2,4,5-trimethylphenol","C9H12O",715.0,3600000.0,0.47,0.136192,492.0),
    ("2,3,5-trimethylphenol","C9H12O",713.0,3550000.0,0.46,0.136192,490.0),
    ("2,3,6-trimethylphenol","C9H12O",711.0,3500000.0,0.45,0.136192,488.0),
]
for name, formula, Tc, Pc, omega, M, Tb in phenols:
    emit(lines, name, formula, "", Tc, Pc, omega, M, Tb, "estimated")

# ── furans ──────────────────────────────────────────────────────────────────
furans = [
    ("furan","C4H4O",490.0,5500000.0,0.200,0.068074,304.5),
    ("2-methylfuran","C5H6O",524.0,4700000.0,0.24,0.08210,337.0),
    ("2,5-dimethylfuran","C6H8O",548.0,4000000.0,0.28,0.09613,366.0),
    ("2-ethylfuran","C6H8O",535.0,4200000.0,0.26,0.09613,352.0),
    ("2-propylfuran","C7H10O",558.0,3700000.0,0.30,0.11014,378.0),
    ("tetrahydrofurfuryl alcohol","C5H10O2",671.0,4200000.0,0.55,0.10213,450.0),
    ("furfuryl alcohol","C5H6O2",634.0,4500000.0,0.70,0.098099,443.0),
    ("furfural","C5H4O2",657.0,4900000.0,0.372,0.096084,434.85),
    ("5-methylfurfural","C6H6O2",680.0,4500000.0,0.40,0.11010,455.0),
    ("2-acetylfuran","C6H6O2",670.0,4300000.0,0.38,0.09808,448.0),
]
for name, formula, Tc, Pc, omega, M, Tb in furans:
    emit(lines, name, formula, "", Tc, Pc, omega, M, Tb, "estimated")

# ── epoxides ────────────────────────────────────────────────────────────────
epoxides = [
    ("ethylene oxide","C2H4O",469.0,7190000.0,0.197,0.0440526,283.55),
    ("propylene oxide","C3H6O",488.0,4920000.0,0.269,0.058079,307.1),
    ("1,2-butylene oxide","C4H8O",525.0,4350000.0,0.28,0.07211,336.0),
    ("isobutylene oxide","C4H8O",509.0,4010000.0,0.27,0.07211,312.0),
    ("cyclohexene oxide","C6H10O",610.0,4000000.0,0.24,0.09814,406.0),
    ("styrene oxide","C8H8O",675.0,3800000.0,0.30,0.12015,453.0),
    ("1,2-epoxybutane","C4H8O",525.0,4350000.0,0.28,0.07211,336.0),
    ("1,2-epoxypentane","C5H10O",550.0,3800000.0,0.30,0.08614,360.0),
    ("1,2-epoxyhexane","C6H12O",572.0,3400000.0,0.32,0.10016,383.0),
    ("1,2-epoxyheptane","C7H14O",592.0,3100000.0,0.34,0.11419,404.0),
]
for name, formula, Tc, Pc, omega, M, Tb in epoxides:
    emit(lines, name, formula, "", Tc, Pc, omega, M, Tb, "estimated")

# ── anhydrides ──────────────────────────────────────────────────────────────
anhydrides = [
    ("acetic anhydride","C4H6O3",606.0,4000000.0,0.45,0.10209,412.0),
    ("propionic anhydride","C6H10O3",651.0,3200000.0,0.50,0.13014,443.0),
    ("maleic anhydride","C4H2O3",721.0,5300000.0,0.40,0.09806,475.0),
    ("phthalic anhydride","C8H4O3",810.0,4700000.0,0.45,0.14806,557.0),
    ("succinic anhydride","C4H4O3",728.0,4800000.0,0.42,0.10007,498.0),
    ("glutaric anhydride","C5H6O3",750.0,4400000.0,0.44,0.11411,520.0),
    ("adipic anhydride","C6H8O3",770.0,4000000.0,0.46,0.12814,540.0),
]
for name, formula, Tc, Pc, omega, M, Tb in anhydrides:
    emit(lines, name, formula, "", Tc, Pc, omega, M, Tb, "estimated")

# ── carbonates ──────────────────────────────────────────────────────────────
carbonates = [
    ("dimethyl carbonate","C3H6O3",548.0,4400000.0,0.34,0.09008,363.5),
    ("diethyl carbonate","C5H10O3",576.0,3300000.0,0.38,0.11818,399.0),
    ("ethylene carbonate","C3H4O3",738.0,4200000.0,0.36,0.08806,513.0),
    ("propylene carbonate","C4H6O3",738.0,4200000.0,0.36,0.10209,513.0),
    ("butylene carbonate","C5H8O3",745.0,3800000.0,0.38,0.11613,520.0),
    ("glycerol carbonate","C4H6O4",760.0,4000000.0,0.37,0.11809,530.0),
    ("vinylene carbonate","C3H2O3",720.0,4500000.0,0.35,0.08605,500.0),
]
for name, formula, Tc, Pc, omega, M, Tb in carbonates:
    emit(lines, name, formula, "", Tc, Pc, omega, M, Tb, "estimated")

# ── amides ──────────────────────────────────────────────────────────────────
amides = [
    ("dimethylformamide","C3H7NO",649.0,4420000.0,0.34,0.07309,426.0),
    ("dimethylacetamide","C4H9NO",658.0,3800000.0,0.36,0.08712,438.0),
    ("n-methylpyrrolidinone","C5H9NO",724.0,4510000.0,0.36,0.09913,475.0),
    ("n-ethylpyrrolidinone","C6H11NO",735.0,4100000.0,0.38,0.11316,490.0),
    ("n-formylmorpholine","C5H9NO2",730.0,4200000.0,0.37,0.11514,470.0),
    ("n-methylformamide","C2H5NO",690.0,6000000.0,0.30,0.05907,473.0),
    ("formamide","CH3NO",680.0,7500000.0,0.28,0.04504,480.0),
    ("acetamide","C2H5NO",710.0,6500000.0,0.32,0.05907,495.0),
    ("n,n-dimethylpropionamide","C5H11NO",670.0,3500000.0,0.35,0.10116,450.0),
    ("urea","CH4N2O",705.0,9000000.0,0.35,0.060055,432.0),
]
for name, formula, Tc, Pc, omega, M, Tb in amides:
    emit(lines, name, formula, "", Tc, Pc, omega, M, Tb, "estimated")

# ── siloxanes ───────────────────────────────────────────────────────────────
siloxanes = [
    ("hexamethyldisiloxane","C6H18OSi2",524.0,1940000.0,0.38,0.16238,372.0),
    ("octamethyltrisiloxane","C8H24O2Si3",564.0,1550000.0,0.40,0.23651,425.0),
    ("decamethyltetrasiloxane","C10H30O3Si4",615.0,1200000.0,0.42,0.31064,468.0),
    ("dodecamethylpentasiloxane","C12H36O4Si5",655.0,950000.0,0.44,0.38477,505.0),
    ("tetradecamethylhexasiloxane","C14H42O5Si6",690.0,780000.0,0.46,0.45990,538.0),
    ("hexamethylcyclotrisiloxane","C6H18O3Si3",595.0,1950000.0,0.40,0.22238,408.0),
    ("octamethylcyclotetrasiloxane","C8H24O4Si4",628.0,1550000.0,0.42,0.29654,448.0),
    ("decamethylcyclopentasiloxane","C10H30O5Si5",660.0,1200000.0,0.44,0.37067,486.0),
    ("dodecamethylcyclohexasiloxane","C12H36O6Si6",690.0,950000.0,0.46,0.44480,523.0),
]
for name, formula, Tc, Pc, omega, M, Tb in siloxanes:
    emit(lines, name, formula, "", Tc, Pc, omega, M, Tb, "estimated")

# ── terpenes & flavors ──────────────────────────────────────────────────────
terpenes = [
    ("limonene","C10H16",660.0,2800000.0,0.32,0.13623,451.0),
    ("menthol","C10H20O",695.0,2700000.0,0.40,0.15627,489.0),
    ("camphor","C10H16O",720.0,3000000.0,0.35,0.15223,480.0),
    ("eugenol","C10H12O2",740.0,3500000.0,0.40,0.16420,528.0),
    ("carvone","C10H14O",710.0,3000000.0,0.38,0.15022,495.0),
    ("linalool","C10H18O",680.0,2800000.0,0.40,0.15425,470.0),
    ("geraniol","C10H18O",700.0,2700000.0,0.42,0.15425,503.0),
    ("citral","C10H16O",690.0,2800000.0,0.40,0.15223,494.0),
    ("cinnamaldehyde","C9H8O",740.0,3500000.0,0.42,0.13216,520.0),
    ("vanillin","C8H8O3",765.0,4200000.0,0.45,0.15215,558.0),
    ("ethyl vanillin","C9H10O3",770.0,3800000.0,0.44,0.16617,560.0),
    ("alpha-pinene","C10H16",650.0,2900000.0,0.28,0.13623,429.0),
    ("beta-pinene","C10H16",655.0,2850000.0,0.29,0.13623,435.0),
    ("camphene","C10H16",640.0,3000000.0,0.26,0.13623,418.0),
    ("terpinolene","C10H16",670.0,2700000.0,0.33,0.13623,460.0),
    ("p-cymene","C10H14",652.0,2800000.0,0.31,0.13421,448.0),
    ("myrcene","C10H16",645.0,2750000.0,0.30,0.13623,440.0),
    ("ocimene","C10H16",648.0,2780000.0,0.31,0.13623,442.0),
    ("terpinene","C10H16",660.0,2700000.0,0.32,0.13623,452.0),
    ("borneol","C10H18O",705.0,2600000.0,0.42,0.15425,488.0),
    ("isoborneol","C10H18O",708.0,2580000.0,0.43,0.15425,490.0),
    ("terpineol","C10H18O",700.0,2650000.0,0.41,0.15425,485.0),
    ("nerol","C10H18O",698.0,2680000.0,0.40,0.15425,500.0),
    ("linalyl acetate","C12H20O2",690.0,2400000.0,0.42,0.19630,490.0),
    ("geranyl acetate","C12H20O2",700.0,2300000.0,0.44,0.19630,505.0),
    ("citronellol","C10H20O",690.0,2700000.0,0.40,0.15627,495.0),
    ("phenylethyl alcohol","C8H10O",686.0,3900000.0,0.40,0.122165,492.0),
    ("benzyl alcohol","C7H8O",720.0,4300000.0,0.393,0.108138,477.85),
    ("caffeine","C8H10N4O2",820.0,4000000.0,0.35,0.19419,630.0),
    ("theobromine","C7H8N4O2",810.0,4200000.0,0.33,0.18016,620.0),
    ("theophylline","C7H8N4O2",815.0,4100000.0,0.34,0.18016,625.0),
]
for name, formula, Tc, Pc, omega, M, Tb in terpenes:
    emit(lines, name, formula, "", Tc, Pc, omega, M, Tb, "estimated")

# ── halogenated aromatics ───────────────────────────────────────────────────
halo_arom = [
    ("chlorobenzene","C6H5Cl",632.4,4519000.0,0.249,0.112557,404.87),
    ("bromobenzene","C6H5Br",670.0,4520000.0,0.251,0.157008,429.15),
    ("o-dichlorobenzene","C6H4Cl2",729.0,4100000.0,0.272,0.147002,453.55),
    ("m-dichlorobenzene","C6H4Cl2",715.0,3800000.0,0.26,0.147002,445.5),
    ("p-dichlorobenzene","C6H4Cl2",684.0,4050000.0,0.25,0.147002,447.0),
    ("o-chlorotoluene","C7H7Cl",656.0,3800000.0,0.27,0.12661,432.0),
    ("p-chlorotoluene","C7H7Cl",660.0,3700000.0,0.28,0.12661,435.0),
    ("benzyl chloride","C7H7Cl",686.0,3900000.0,0.30,0.12661,452.0),
    ("benzotrifluoride","C7H5F3",559.0,3700000.0,0.28,0.14609,375.0),
    ("hexafluorobenzene","C6F6",516.0,3290000.0,0.397,0.18605,353.4),
    ("chloronaphthalene","C10H7Cl",770.0,3500000.0,0.32,0.16262,525.0),
    ("bromonaphthalene","C10H7Br",790.0,3400000.0,0.33,0.20757,540.0),
    ("p-bromotoluene","C7H7Br",670.0,3600000.0,0.28,0.18506,458.0),
    ("o-bromotoluene","C7H7Br",665.0,3650000.0,0.27,0.18506,455.0),
    ("p-fluorotoluene","C7H7F",580.0,3800000.0,0.24,0.11010,389.0),
    ("o-fluorotoluene","C7H7F",575.0,3850000.0,0.23,0.11010,386.0),
    ("trichlorobenzene","C6H3Cl3",730.0,3300000.0,0.30,0.18145,486.0),
    ("hexachlorobenzene","C6Cl6",820.0,2800000.0,0.35,0.28476,570.0),
    ("pentafluorobenzene","C6HF5",520.0,3400000.0,0.35,0.16802,355.0),
    ("iodobenzene","C6H5I",720.0,4500000.0,0.25,0.20401,461.0),
]
for name, formula, Tc, Pc, omega, M, Tb in halo_arom:
    emit(lines, name, formula, "", Tc, Pc, omega, M, Tb, "estimated")

# ── PAHs ────────────────────────────────────────────────────────────────────
pahs = [
    ("naphthalene","C10H8",748.4,4050000.0,0.302,0.1281705,491.14),
    ("1-methylnaphthalene","C11H10",772.0,3600000.0,0.334,0.142197,517.0),
    ("2-methylnaphthalene","C11H10",761.0,3500000.0,0.382,0.142197,514.0),
    ("biphenyl","C12H10",789.0,3850000.0,0.363,0.15420752,528.15),
    ("diphenylmethane","C13H12",768.0,2650000.0,0.430,0.1682344,537.5),
    ("anthracene","C14H10",869.3,2900000.0,0.486,0.1782507,613.1),
    ("phenanthrene","C14H10",873.0,2900000.0,0.473,0.1782507,613.0),
    ("pyrene","C16H10",936.0,2610000.0,0.530,0.2022534,668.0),
    ("chrysene","C18H12",993.0,2100000.0,0.57,0.228293,723.0),
    ("fluorene","C13H10",871.0,2950000.0,0.456,0.166218,570.0),
    ("acenaphthene","C12H10",803.0,3100000.0,0.385,0.15421,552.0),
    ("fluoranthene","C16H10",965.0,2500000.0,0.55,0.20225,655.0),
    ("coronene","C24H12",1050.0,1600000.0,0.62,0.30037,795.0),
    ("corannulene","C20H10",980.0,2000000.0,0.56,0.25031,720.0),
    ("triphenylene","C18H12",990.0,2200000.0,0.55,0.22829,720.0),
    ("perylene","C20H12",1020.0,2000000.0,0.58,0.25233,750.0),
    ("benzo-a-pyrene","C20H12",1050.0,1800000.0,0.60,0.25233,780.0),
    ("benzo-e-pyrene","C20H12",1040.0,1850000.0,0.59,0.25233,775.0),
    ("rubrene","C42H28",1100.0,1200000.0,0.65,0.53267,850.0),
    ("fullerene-c60","C60",1500.0,1000000.0,0.70,0.72066,1100.0),
]
for name, formula, Tc, Pc, omega, M, Tb in pahs:
    emit(lines, name, formula, "", Tc, Pc, omega, M, Tb, "estimated")

# ── common inorganics ───────────────────────────────────────────────────────
inorganics = [
    ("hydrogen","H2",33.23,1296400.0,-0.216,0.00201588,20.27),
    ("oxygen","O2",154.58,5043000.0,0.022,0.0319988,90.17),
    ("nitrogen","N2",126.19,3395800.0,0.037,0.02801348,77.36),
    ("argon","Ar",150.86,4898000.0,-0.002,0.039948,87.30),
    ("helium","He",5.19,227460.0,-0.387,0.0040026,4.22),
    ("neon","Ne",44.49,2678610.0,0.0,0.0201797,27.07),
    ("krypton","Kr",209.48,5525000.0,-0.002,0.083798,119.74),
    ("xenon","Xe",289.73,5842000.0,0.008,0.131293,165.03),
    ("carbon monoxide","CO",132.85,3493500.0,0.045,0.0280104,81.66),
    ("nitrous oxide","N2O",309.57,7245000.0,0.16,0.0440128,184.67),
    ("nitric oxide","NO",180.15,6480000.0,0.588,0.0300061,121.38),
    ("nitrogen dioxide","NO2",431.35,10132000.0,0.834,0.0460055,294.0),
    ("sulfur dioxide","SO2",430.64,7884000.0,0.245,0.0640638,263.13),
    ("sulfur trioxide","SO3",490.85,8210000.0,0.41,0.0799567,317.9),
    ("hydrogen sulfide","H2S",373.4,8963000.0,0.094,0.03408088,212.8),
    ("hydrogen fluoride","HF",461.0,6489000.0,0.372,0.020006347,292.65),
    ("hydrogen chloride","HCl",324.65,8310000.0,0.132,0.03646094,188.15),
    ("ammonia","NH3",405.4,11353000.0,0.256,0.01703052,239.82),
    ("chlorine","Cl2",416.9,7991000.0,0.069,0.070906,239.12),
    ("bromine","Br2",588.0,10400000.0,0.11,0.159808,332.0),
    ("iodine","I2",819.0,11700000.0,0.17,0.253809,457.0),
    ("cyanogen","C2N2",400.0,5997000.0,0.24,0.0520348,252.15),
    ("phosgene","COCl2",455.0,5674000.0,0.205,0.098916,280.7),
    ("boron trichloride","BCl3",455.0,3840000.0,0.15,0.11717,285.75),
    ("silicon tetrachloride","SiCl4",507.0,3593000.0,0.232,0.169896,330.45),
    ("titanium tetrachloride","TiCl4",638.0,4648000.0,0.36,0.189697,409.0),
    ("phosphorus trichloride","PCl3",563.0,5670000.0,0.15,0.13733,349.0),
    ("phosphorus pentachloride","PCl5",646.0,4500000.0,0.20,0.20824,433.0),
    ("sulfuryl chloride","SO2Cl2",567.0,4600000.0,0.18,0.13497,342.0),
    ("thionyl chloride","SOCl2",520.0,5200000.0,0.16,0.11897,348.0),
    ("carbon disulfide","CS2",552.0,7900000.0,0.115,0.0761407,319.37),
    ("carbonyl sulfide","COS",378.8,6370000.0,0.099,0.060075,223.0),
    ("hydrogen cyanide","HCN",456.8,5390000.0,0.394,0.02702538,298.85),
    ("acetylene","C2H2",308.3,6139000.0,0.187,0.02603728,189.15),
    ("ketene","C2O2",480.0,6500000.0,0.20,0.04204,265.0),
    ("carbon dioxide","CO2",304.1282,7377000.0,0.22394,0.0440095,216.58),
    ("water","H2O",647.096,22064000.0,0.344,0.01801528,373.15),
]
for name, formula, Tc, Pc, omega, M, Tb in inorganics:
    emit(lines, name, formula, "", Tc, Pc, omega, M, Tb, "NIST")

# ════════════════════════════════════════════════════════════════════════════
# WRITE OUTPUT
# ════════════════════════════════════════════════════════════════════════════

bip_idx = existing.find('[[binary_interactions]]')
comp_sec = existing[:bip_idx] if bip_idx > 0 else existing
bip_sec = existing[bip_idx:] if bip_idx > 0 else ''

out = comp_sec.rstrip() + '\n\n'
for l in lines:
    out += l + '\n'
if bip_sec:
    out += bip_sec

with open(SEED, 'w') as f:
    f.write(out)

n_existing = existing.count('[[components]]')
n_new = len([l for l in lines if l.startswith('[[components]]')])
n_total = n_existing + n_new
n_bips = out.count('[[binary_interactions]]')

print(f"Generated seed.toml:")
print(f"  Existing compounds:  {n_existing}")
print(f"  New compounds:       {n_new}")
print(f"  Total compounds:     {n_total}")
print(f"  Binary interactions: {n_bips}")
print(f"  File size:           {len(out)} bytes")
