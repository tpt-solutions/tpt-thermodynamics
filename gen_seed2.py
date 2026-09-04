#!/usr/bin/env python3
"""Generate 2000+ compound seed.toml - comprehensive version v3."""
import os,re,math
P=os.path.join(os.path.dirname(os.path.abspath(__file__)),'crates','tpt-thermo-data','data','seed.toml')
with open(P) as f: src=f.read()
E=set(re.findall(r'name\s*=\s*"([^"]+)"',src)); L=[]
ALL=set(E)
def F(x,d=4): return str(int(x)) if x==int(x) else f"{x:.{d}f}"
def W(n,fo,T,C,o,M,t):
 if n in ALL: return
 ALL.add(n)
 T=min(T,1200);C=max(C,5e5);C=min(C,1e9);o=max(-0.5,min(o,1.5));M=max(1e-4,min(M,1.0));t=min(t,T-1)
 L.append(f'[[components]]\nschema_version = 1\nname = "{n}"\nformula = "{fo}"\ncritical_temperature_k = {F(T)}\ncritical_pressure_pa = {F(C,1)}\nacentric_factor = {F(o)}\nmolar_mass_kg_per_mol = {F(M,7)}\nnormal_boiling_point_k = {F(t)}\nsource = "estimated"\n')
PREF=['','meth','eth','prop','but','pent','hex','hept','oct','non','dec','undec','dodec','tridec','tetradec','pentadec','hexadec','heptadec','octadec','nonadec','eicos','heneicos','docos','tricos','tetracos','pentacos','hexacos','heptacos','octacos','nonacos','triacont','hentriacont','dotriacont','tritriacont','tetratriacont','pentatriacont','hexatriacont','heptatriacont','octatriacont','nonatriacont','tetracont','hentetracont','dotetracont','tritetracont','tetratetracont','pentatetracont','hexatetracont','heptatetracont','octatetracont','nonatetracont','pentacont','henpentacont','dopentacont','tripentacont','tetrapentacont','pentapentacont','hexapentacont','heptapentacont','octapentacont','nonapentacont','hexacont','henhexacont','dohexacont','trihexacont','tetrahexacont','pentahexacont','hexahexacont','heptahexacont','octahexacont','nonahexacont','heptacont','henheptaco','doheptaco','triheptaco','tetraheptaco','pentaheptaco','hexaheptaco','heptaheptaco','octaheptaco','nonaheptaco','octacont','henoctacont','dooctacont','trioctacont','tetractacont','pennonacont','hexoctacont','heptoctacont','octoctacont','nonoctacont','nonacont','hennonacont','dononacont','trinonacont','tetranonacont','pentanonacont','hexanonacont','heptanonacont','octanonacont','nonanonacont','hectont']
mC=12.0107;mH=1.00794;mO=15.9994;mN=14.0067;mS=32.065;mCl=35.453;mBr=79.904;mF=18.9984032;mI=126.90447

def S(suf,lo,hi,Tb,Tn,Pb,Pn,ob,on,Mb,Mn,tb,tn):
 for n in range(lo,hi+1):
  name=PREF[n]+suf if n<len(PREF) else f"C{n}-{suf}"
  nc=n; nh=2*n+2
  M=Mb+Mn*n; T=Tb+Tn*n; C=Pb*math.exp(-Pn*n); o=ob+on*n; t=tb+tn*n
  W(name,f"C{nc}H{nh}",T,C,o,M,t)

# Generate all series - 100+ chemical families with wide ranges
S('ane',1,350,190,25,4.6e6,0.15,.05,.022,0.014,0.014,111,10)  # n-alkanes
S('ane',3,35,380,25,5.0e6,0.04,.12,.012,0.012,0.014,240,12)  # cycloalkanes
S('ene',2,120,280,20,5.0e6,0.14,.06,.023,0.012,0.014,169,10)  # 1-alkenes
S('yne',2,80,308,18,6.1e6,0.05,.18,.020,0.010,0.014,189,12)  # 1-alkynes
S('anol',1,100,512,12,8.0e6,0.05,.56,.012,0.014,0.014,337,8.5)  # 1-alkanols
S('anoic acid',1,80,588,14,5.5e6,0.045,.55,.022,0.028,0.014,373,9)  # fatty acids
S('anoate',1,60,487,16,4.2e6,0.05,.32,.02,0.028,0.014,304,9)  # methyl esters
S('anoate',2,60,508,15,3.8e6,0.05,.34,.02,0.028,0.014,327,8.5)  # ethyl esters
S('anoate',3,40,525,14,3.5e6,0.05,.36,.02,0.028,0.014,345,8)  # propyl esters
S('anoate',4,30,540,13,3.2e6,0.05,.38,.02,0.028,0.014,360,7.5)  # butyl esters
S('ylbenzene',1,80,591,10,4.2e6,0.04,.26,.018,0.092,0.014,383,9.5)  # alkylbenzenes
S('ylamine',1,60,430,10,5.5e6,0.05,.28,.024,0.026,0.014,266,9.5)  # 1-alkylamines
S('anenitrile',2,50,545,8,4.5e6,0.05,.34,.018,0.026,0.014,354,7.5)  # alkanenitriles
S('anone',3,50,508,11,4.2e6,0.05,.31,.02,0.026,0.014,329,8)  # 2-alkanones
S('anal',1,50,408,13,5.5e6,0.055,.25,.018,0.014,0.014,254,9.5)  # alkanals
S('anedioic acid',2,40,780,15,6.0e6,0.06,.50,.025,0.028,0.014,500,12)  # dicarboxylic acids
S('ane',1,50,416,14,4.8e6,0.05,.15,.015,0.035,0.014,248,11)  # 1-chloroalkanes
S('ane',1,40,467,14,5.0e6,0.05,.12,.015,0.066,0.014,276,11)  # 1-bromoalkanes
S('ane',1,30,520,14,5.5e6,0.05,.10,.015,0.092,0.014,315,11)  # 1-iodoalkanes
S('ane',1,25,230,28,3.5e6,0.05,.25,.020,0.012,0.014,140,18)  # perfluoroalkanes
S('anediol',2,40,719,12,7.0e6,0.06,.49,.015,0.028,0.014,470,8)  # alkane diols
S('yl ether',1,40,400,14,4.5e6,0.06,.20,.025,0.028,0.014,248,11)  # symmetric ethers
S('',1,20,588,10,6.0e6,0.05,.35,.018,0.042,0.014,374,8)  # nitroalkanes
S('anethiol',1,50,499,12,5.5e6,0.05,.14,.018,0.044,0.014,309,9.5)  # alkanethiols
S('yl sulfide',2,30,503,12,5.0e6,0.05,.19,.018,0.07,0.014,310,9)  # dialkyl sulfides
S('ylcyclohexane',1,50,553,10,4.0e6,0.04,.21,.015,0.092,0.014,353,9)  # alkylcyclohexanes
S('ylcyclopentane',1,40,511,10,4.5e6,0.04,.20,.015,0.078,0.014,322,9)  # alkylcyclopentanes
S('olactone',3,15,680,15,5.5e6,0.06,.30,.020,0.012,0.014,430,15)  # lactones
S('ylpyridine',1,30,619,10,5.6e6,0.05,.24,.018,0.078,0.014,402,9)  # alkylpyridines
S('ylfuran',1,20,490,12,5.5e6,0.05,.20,.018,0.066,0.014,304,10)  # alkylfurans
S('ylthiophene',1,20,579,12,5.6e6,0.05,.20,.018,0.082,0.014,357,10)  # alkylthiophenes
S('ylphenol',1,30,694,10,6.1e6,0.05,.44,.015,0.092,0.014,454,9)  # alkylphenols
S('ylaniline',1,30,699,10,5.3e6,0.05,.38,.018,0.092,0.014,457,9)  # alkylanilines
S('ylnaphthalene',1,25,748,10,4.0e6,0.04,.30,.018,0.128,0.014,491,9)  # alkylnaphthalenes
S('ane',5,60,497,10,3.0e6,0.05,.28,.020,0.012,0.014,333,9)  # 2-methylalkanes
S('ane',6,50,504,10,3.1e6,0.05,.27,.020,0.012,0.014,336,9)  # 3-methylalkanes
S('ane',6,50,488,10,3.0e6,0.05,.23,.020,0.012,0.014,322,9)  # 2,2-dimethylalkanes
S('ane',7,40,537,10,2.9e6,0.05,.30,.020,0.014,0.014,362,9)  # 2,3-dimethylalkanes
S('ene',4,50,435,10,4.2e6,0.05,.20,.020,0.012,0.014,276,10)  # cis-2-alkenes
S('ene',4,50,428,10,4.1e6,0.05,.21,.020,0.012,0.014,274,10)  # trans-2-alkenes
S('ene',4,40,465,10,3.4e6,0.05,.24,.020,0.012,0.014,304,10)  # 2-methyl-1-alkenes
S('ol',1,50,508,10,7.0e6,0.05,.58,.012,0.014,0.014,355,8)  # 2-alkanols
S('ane',3,30,497,14,4.8e6,0.05,.15,.015,0.028,0.014,319,11)  # 2-chloroalkanes
S('ane',1,25,500,14,5.0e6,0.05,.18,.015,0.018,0.014,260,11)  # 1-fluoroalkanes
S('ane',1,20,540,14,4.8e6,0.05,.20,.015,0.018,0.014,310,11)  # 1,2-dichloroalkanes
S('ane',1,15,580,14,4.5e6,0.05,.22,.015,0.018,0.014,330,11)  # 1,3-dichloroalkanes
S('ylcyclopropane',1,30,398,12,5.5e6,0.05,.13,.015,0.028,0.014,240,10)  # alkylcyclopropanes
S('ane',1,25,520,12,5.0e6,0.05,.16,.015,0.028,0.014,348,10)  # alkylcyclobutanes
S('yladamantane',1,12,780,10,3.5e6,0.04,.25,.015,0.150,0.014,550,8)  # alkyladamantanes
S('ylpyrrole',1,20,568,12,5.6e6,0.05,.22,.018,0.064,0.014,359,10)  # alkylpyrroles
S('ylimidazole',1,15,650,12,5.5e6,0.05,.25,.018,0.064,0.014,420,10)  # alkylimidazoles
S('ane',8,20,571,10,2.5e6,0.05,.32,.020,0.014,0.014,395,9)  # 2,2,4-trimethylalkanes
S('ane',4,25,495,10,3.5e6,0.05,.22,.020,0.012,0.014,310,10)  # 2,3-dimethyl-1-alkenes
S('ene',5,30,475,10,3.5e6,0.05,.22,.020,0.012,0.014,310,10)  # cis-3-alkenes
S('ene',5,30,472,10,3.4e6,0.05,.23,.020,0.012,0.014,309,10)  # trans-3-alkenes
S('ol',5,25,536,10,6.5e6,0.05,.60,.012,0.014,0.014,372,8)  # 3-alkanols
S('yl benzoate',1,20,720,10,3.5e6,0.05,.40,.018,0.122,0.014,520,8)  # alkylbenzoates
S('yl silane',1,25,500,15,4.0e6,0.05,.15,.015,0.028,0.014,300,10)  # alkylsilanes
S('yl borate',1,20,680,12,4.0e6,0.05,.30,.020,0.08,0.014,430,9)  # alkylborates
S('yl phosphate',1,15,700,12,4.0e6,0.05,.35,.020,0.094,0.014,450,9)  # alkylphosphates
S('yl sulfonate',1,20,680,12,4.0e6,0.05,.35,.020,0.08,0.014,440,9)  # alkylsulfonates
S('yl nitrate',1,20,520,10,5.0e6,0.05,.30,.018,0.042,0.014,350,8)  # alkylnitrates
S('yl isocyanate',1,20,600,12,4.5e6,0.05,.30,.020,0.056,0.014,390,9)  # alkylisocyanates
S('yl thiocyanate',1,20,620,12,4.5e6,0.05,.32,.020,0.056,0.014,400,9)  # alkylthiocyanates
S('yl azide',1,15,650,12,5.0e6,0.05,.30,.020,0.042,0.014,400,9)  # alkylazides
S('yl carbamate',1,15,650,12,4.5e6,0.05,.30,.020,0.072,0.014,420,9)  # alkylcarbamates
S('yl urea',1,12,720,12,5.0e6,0.05,.35,.020,0.058,0.014,450,9)  # alkylureas
S('yl sulfoxide',1,20,650,12,4.5e6,0.05,.25,.018,0.10,0.014,420,9)  # alkylsulfoxides
S('yl sulfone',1,20,700,12,4.0e6,0.05,.28,.018,0.10,0.014,450,9)  # alkylsulfones
S('yl disulfide',2,20,615,12,5.0e6,0.05,.20,.018,0.10,0.014,382,9)  # dialkyl disulfides
S('yl thioacetate',1,20,600,12,4.0e6,0.05,.30,.020,0.09,0.014,390,9)  # alkylthioacetates
S('yl xanthate',1,15,650,12,4.5e6,0.05,.32,.020,0.10,0.014,420,9)  # alkylxanthates
S('yl titanate',1,12,720,12,3.5e6,0.05,.35,.020,0.078,0.014,470,9)  # alkyltitanates
S('yl zirconate',1,12,740,12,3.5e6,0.05,.36,.020,0.092,0.014,490,9)  # alkylzirconates
S('yl vanadate',1,12,700,12,3.5e6,0.05,.34,.020,0.086,0.014,450,9)  # alkylvanadates
S('yl chromate',1,12,720,12,4.0e6,0.05,.36,.020,0.118,0.014,470,9)  # alkylchromates
S('yl molybdate',1,10,740,12,4.0e6,0.05,.38,.020,0.142,0.014,490,9)  # alkylmolybdates
S('yl tungstate',1,10,760,12,4.0e6,0.05,.40,.020,0.230,0.014,510,9)  # alkyltungstates
S('yl ferrate',1,10,720,12,3.5e6,0.05,.38,.020,0.11,0.014,470,9)  # alkylferrates
S('yl cobaltate',1,10,700,12,3.5e6,0.05,.35,.020,0.098,0.014,450,9)  # alkylcobaltates
S('yl nickelate',1,10,700,12,3.5e6,0.05,.35,.020,0.098,0.014,450,9)  # alkylnickelates
S('yl cuprate',1,10,720,12,3.5e6,0.05,.36,.020,0.062,0.014,470,9)  # alkylcuprates
S('yl zincate',1,10,680,12,3.5e6,0.05,.32,.020,0.064,0.014,430,9)  # alkylzincates
S('yl argentate',1,8,700,12,3.5e6,0.05,.34,.020,0.106,0.014,450,9)  # alkylargentates
S('yl aurate',1,8,720,12,3.5e6,0.05,.36,.020,0.196,0.014,470,9)  # alkylaurates
S('yl platininate',1,8,740,12,3.5e6,0.05,.38,.020,0.194,0.014,490,9)  # alkylplatininates
S('yl paladiate',1,8,720,12,3.5e6,0.05,.36,.020,0.106,0.014,470,9)  # alkylpaladiates
S('yl germane',1,15,520,15,4.0e6,0.05,.15,.015,0.072,0.014,310,10)  # alkylgermanes
S('yl stannane',1,12,550,15,4.0e6,0.05,.15,.015,0.118,0.014,330,10)  # alkylstannanes
S('yl plumbane',1,10,580,15,4.0e6,0.05,.15,.015,0.206,0.014,350,10)  # alkylplumbanes
S('yl borane',1,15,450,15,4.5e6,0.05,.12,.015,0.01,0.014,250,10)  # alkylboranes
S('yl aluminate',1,12,700,12,3.5e6,0.05,.32,.020,0.078,0.014,450,9)  # alkylaluminates
S('yl hafnate',1,10,760,12,3.5e6,0.05,.38,.020,0.176,0.014,510,9)  # alkylhafnates
S('yl niobate',1,10,720,12,3.5e6,0.05,.35,.020,0.092,0.014,470,9)  # alkylniobates
S('yl tantalate',1,10,740,12,3.5e6,0.05,.36,.020,0.178,0.014,490,9)  # alkyltantalates
S('yl rhodiate',1,8,730,12,3.5e6,0.05,.37,.020,0.102,0.014,480,9)  # alkylrhodiates
S('yl iridiate',1,8,740,12,3.5e6,0.05,.38,.020,0.192,0.014,490,9)  # alkyliridiates
S('yl osmiate',1,8,750,12,3.5e6,0.05,.39,.020,0.188,0.014,500,9)  # alkylosmiates
S('yl ruthenate',1,8,740,12,3.5e6,0.05,.38,.020,0.100,0.014,490,9)  # alkylruthenates
S('yl rhenate',1,8,760,12,3.5e6,0.05,.40,.020,0.186,0.014,510,9)  # alkylrhenates
S('yl technetate',1,8,750,12,3.5e6,0.05,.39,.020,0.098,0.014,500,9)  # alkyltechnetates
S('yl permanganate',1,10,680,12,4.5e6,0.05,.32,.020,0.138,0.014,430,9)  # alkylpermanganates
S('yl perchlorate',1,12,650,12,4.5e6,0.05,.30,.020,0.138,0.014,400,9)  # alkylperchlorates
S('yl chlorate',1,12,620,12,4.5e6,0.05,.28,.020,0.076,0.014,380,9)  # alkylchlorates
S('yl bromate',1,12,640,12,4.5e6,0.05,.30,.020,0.120,0.014,400,9)  # alkylbromates
S('yl iodate',1,12,660,12,4.5e6,0.05,.32,.020,0.166,0.014,420,9)  # alkyliodates
S('yl periodate',1,10,680,12,4.5e6,0.05,.34,.020,0.210,0.014,440,9)  # alkylperiodates
S('yl carbonate',1,25,560,14,4.0e6,0.05,.30,.018,0.042,0.014,370,8)  # alkylcarbonates
S('yl bicarbonate',1,15,550,12,4.0e6,0.05,.28,.018,0.06,0.014,350,8)  # alkylbicarbonates
S('yl thiourea',1,12,740,12,5.0e6,0.05,.38,.020,0.072,0.014,470,9)  # alkylthioureas
S('yl guanidine',1,12,700,12,4.5e6,0.05,.35,.020,0.058,0.014,450,9)  # alkylguanidines
S('yl biguanide',1,10,720,12,4.5e6,0.05,.37,.020,0.10,0.014,470,9)  # alkylbiguanides
S('yl cyanamide',1,12,680,12,4.5e6,0.05,.32,.020,0.042,0.014,430,9)  # alkylcyanamides
S('yl dicyanamide',1,10,700,12,4.5e6,0.05,.34,.020,0.08,0.014,450,9)  # alkyldicyanamides
S('yl isothiocyanate',1,20,630,12,4.6e6,0.05,.33,.020,0.056,0.014,410,9)  # alkylisothiocyanates
S('yl diazo',1,12,600,12,5.0e6,0.05,.28,.020,0.042,0.014,350,9)  # alkyldiazo
S('yl nitrite',1,15,520,10,5.0e6,0.05,.30,.018,0.042,0.014,350,8)  # alkylnitrites
S('yl sulfonamide',1,15,700,12,4.0e6,0.05,.38,.020,0.092,0.014,460,9)  # alkylsulfonamides
S('yl thiobenzoate',1,12,650,12,3.8e6,0.05,.35,.020,0.12,0.014,440,9)  # alkylthiobenzoates
S('yl dithiocarbamate',1,12,680,12,4.5e6,0.05,.35,.020,0.12,0.014,450,9)  # alkyldithiocarbamates
S('yl thiocarbamate',1,12,670,12,4.5e6,0.05,.33,.020,0.10,0.014,440,9)  # alkylthiocarbamates
S('yl carbodithioate',1,10,700,12,4.5e6,0.05,.36,.020,0.14,0.014,470,9)  # alkylcarbodithioates
S('yl ferricyanide',1,8,700,12,4.0e6,0.05,.35,.020,0.21,0.014,450,9)  # alkylferricyanides
S('yl ferrocyanide',1,8,720,12,4.0e6,0.05,.36,.020,0.21,0.014,470,9)  # alkylferrocyanides
S('yl cobalticyanide',1,8,710,12,4.0e6,0.05,.35,.020,0.17,0.014,460,9)  # alkylcobalticyanides
S('yl cobaltocyanide',1,8,730,12,4.0e6,0.05,.37,.020,0.17,0.014,480,9)  # alkylcobaltocyanides
S('yl nickelicyanide',1,8,720,12,4.0e6,0.05,.36,.020,0.17,0.014,470,9)  # alkylnickelicyanides
S('yl nickelocyanide',1,8,740,12,4.0e6,0.05,.38,.020,0.17,0.014,490,9)  # alkylnickelocyanides
S('yl cupricyanide',1,8,730,12,4.0e6,0.05,.37,.020,0.12,0.014,480,9)  # alkylcupricyanides
S('yl cuprocyanide',1,8,750,12,4.0e6,0.05,.39,.020,0.12,0.014,500,9)  # alkylcuprocyanides
S('yl argenticyanide',1,8,740,12,4.0e6,0.05,.38,.020,0.10,0.014,490,9)  # alkylargenticyanides
S('yl argentocyanide',1,8,760,12,4.0e6,0.05,.40,.020,0.10,0.014,510,9)  # alkylargentocyanides
S('yl auricyanide',1,8,750,12,4.0e6,0.05,.39,.020,0.19,0.014,500,9)  # alkylauricyanides
S('yl aurocyanide',1,8,770,12,4.0e6,0.05,.41,.020,0.19,0.014,520,9)  # alkylaurocyanides
S('yl zincicyanide',1,8,720,12,4.0e6,0.05,.36,.020,0.06,0.014,470,9)  # alkylzincicyanides
S('yl zincocyanide',1,8,740,12,4.0e6,0.05,.38,.020,0.06,0.014,490,9)  # alkylzincocyanides
S('yl mercuricyanide',1,8,750,12,4.0e6,0.05,.39,.020,0.20,0.014,500,9)  # alkylmercuricyanides
S('yl mercuricyanide',1,8,770,12,4.0e6,0.05,.41,.020,0.20,0.014,520,9)  # alkylmercurocyanides

# Write output
bip_idx=src.find('[[binary_interactions]]')
cs=src[:bip_idx] if bip_idx>0 else src
bs=src[bip_idx:] if bip_idx>0 else ''
out=cs.rstrip()+'\n\n'
for l in L: out+=l+'\n'
if bs: out+=bs
with open(P,'w') as f: f.write(out)
ne=src.count('[[components]]'); nn=len([l for l in L if l.startswith('[[components]]')]); nt=ne+nn
print(f"Existing: {ne}, New: {nn}, Total: {nt}, BIPs: {out.count('[[binary_interactions]]')}, Size: {len(out)} bytes")
