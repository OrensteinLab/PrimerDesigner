import math
import pandas as pd
import primer3 as p3
from Bio.Seq import Seq
import General.utils as GU 

# ====== parameter sets from the paper (S-shape), excluding SNP/poly ======
 # (MinO,MaxO, Min,Max, MinL,MaxL)
 
T_OPT = 60.0  
PARAMS_TM    = (T_OPT, T_OPT+1, T_OPT-2, T_OPT+5, T_OPT-5,  T_OPT+10) 

PARAMS_GC    = (0.55, 0.60,     0.45,    0.65,  0.30,     0.70)
# Self score uses *temperature*-like measure; we’ll feed max(hp_tm, ho_tm)
PARAMS_SELF  = (-50.0, 45.0,  -50.0,   50.0,     -50.0,    55.0)
# Number of trailing A's at 3' end
PARAMS_ENDA  = (1.0,   1.0,  0.0,     4.0,      -1.0,     7.0)
# 3' end stability (ΔG, more negative is better)
PARAMS_ENDDG = (-9.0,  -7.0,   -12.0,   -6.0,     -14.0,    -5.0)

# Weights for features (average)
WEIGHTS = {
    'tm':    1.0,
    'gc':    1.0,
    'self':  1.0,
    'endA':  1.0,
    'enddG': 1.0,
}

# ====== piecewise logistic scoring (from paper) ======
def _sigmoid(z: float) -> float:
    if z >= 0:
        ez = math.exp(-z)
        return 1.0 / (1.0 + ez)
    else:
        ez = math.exp(z)
        return ez / (1.0 + ez)

def _left_branch(x: float, MinO: float, MinL: float, Min_: float) -> float:
    k  = 10.0 / (MinO - MinL)
    x0 = (MinO + MinL) / 2.0
    y0 = math.exp(-k * (Min_ - x0))**-1
    L  = 1.0 + math.exp(-k * (Min_ - x0))**-1
    return L * _sigmoid(k * (x - x0)) - y0

def _right_branch(x: float, MaxO: float, MaxL: float, Max_: float) -> float:
    k  = 10.0 / (MaxO - MaxL)
    x0 = (MaxO + MaxL) / 2.0
    y0 = math.exp(-k * (Max_ - x0))**-1
    L  = 1.0 + math.exp(-k * (Max_ - x0))**-1
    return L * _sigmoid(k * (x - x0)) - y0

def piecewise_logistic_score(x: float,
                             MinO: float, MaxO: float,
                             Min_: float, Max_: float,
                             MinL: float, MaxL: float) -> float:
    if x <= MinO:
        return _left_branch(x, MinO, MinL, Min_)
    elif x <= MaxO:
        return 1.0
    else:
        return _right_branch(x, MaxO, MaxL, Max_)

# ====== feature helpers ======
def calc_gc(seq: str) -> float:
    s = seq.upper()
    valid = [b for b in s if b in 'ACGT']
    if not valid:
        return 0.0
    gc = sum(1 for b in valid if b in 'GC')
    return gc / len(valid)

def count_terminal_As_3prime(seq: str) -> int:
    c = 0
    for b in reversed(seq):
        if b.upper() == 'A': c += 1
        else: break
    return c


def _score_tm(x_tm: float) -> float:
    MinO, MaxO, Min_, Max_, MinL, MaxL = PARAMS_TM
    return piecewise_logistic_score(x_tm, MinO, MaxO, Min_, Max_, MinL, MaxL)

def _score_gc(x_gc: float) -> float:
    MinO, MaxO, Min_, Max_, MinL, MaxL = PARAMS_GC
    return piecewise_logistic_score(x_gc, MinO, MaxO, Min_, Max_, MinL, MaxL)

def _score_self(x_self_tm: float) -> float:
    MinO, MaxO, Min_, Max_, MinL, MaxL = PARAMS_SELF
    return piecewise_logistic_score(x_self_tm, MinO, MaxO, Min_, Max_, MinL, MaxL)

def _score_endA(nA3: float) -> float:
    MinO, MaxO, Min_, Max_, MinL, MaxL = PARAMS_ENDA
    return piecewise_logistic_score(nA3, MinO, MaxO, Min_, Max_, MinL, MaxL)

def _score_enddG(dg3: float) -> float:
    MinO, MaxO, Min_, Max_, MinL, MaxL = PARAMS_ENDDG
    return piecewise_logistic_score(dg3, MinO, MaxO, Min_, Max_, MinL, MaxL)

def _weighted_func_from_scores(scores: dict) -> float:
    # Convert “higher-is-better” scores to a cost (lower-is-better)
    total_w = sum(WEIGHTS.values())
    return sum(WEIGHTS[k] * (scores[k]) for k in scores) / total_w

def create_primer_df(sequence_nt, args):

    print("Creating primer df")

    # Forward primers
    primer_f = pd.DataFrame(columns=['seq', 'start', 'stop', 'fr', 'len'])
    primer_f[['seq','start','stop','len']] = GU.subsequences(sequence_nt, args.primer_lmin, args.primer_lmax)
    primer_f['fr'] = 'f'

    # Shift positions so 0 aligns to start of mutreg
    primer_f['start'] = primer_f.start - len(GU.UPSTREAM_NT)
    primer_f['stop']  = primer_f.stop  - len(GU.UPSTREAM_NT)

    # Reverse primers at same loci (reverse-complement sequences)
    primer_r = primer_f[['seq','start','stop','fr','len']].copy()
    primer_r['fr']  = 'r'
    primer_r['seq'] = primer_r.seq.apply(GU.revcomp)

    primer_df = pd.concat([primer_f, primer_r])
    primer_df.sort_values(by=['start','stop','fr'], inplace=True)

    # ---- raw features needed for scoring ----
    primer_df['gc'] = primer_df.seq.apply(calc_gc)
    primer_df['tm'] = primer_df.seq.apply(GU.PCR.calc_tm)

    # hairpin & homodimer (we’ll use Tm for the “self” feature)
    hp_res = primer_df.seq.apply(lambda s: GU.PCR.calc_hairpin(s).todict())
    primer_df['hp_tm'] = hp_res.apply(lambda d: d['tm'])

    ho_res = primer_df.seq.apply(lambda s: GU.PCR.calc_homodimer(s).todict())
    primer_df['ho_tm'] = ho_res.apply(lambda d: d['tm'])

    # 3' end features
    primer_df['endA']  = primer_df.seq.apply(count_terminal_As_3prime)

    # 3' end stability of last 5-nt
    endstab_res = primer_df.seq.apply(lambda s: GU.PCR.calc_end_stability(s[-5:], GU.revcomp(s[-5:])).todict())
    # Primer3 returns dG in cal/mol; convert to kcal/mol like you did for hp/ho:
    primer_df['end_dg'] = endstab_res.apply(lambda res: res['dg'] * 1e-3)

    primer_df['score_tm']    = primer_df.tm.apply(lambda x: _score_tm(x))
    primer_df['score_gc']    = primer_df.gc.apply(lambda x: _score_gc(x))
    primer_df['score_self']  = primer_df[['hp_tm','ho_tm']].max(axis=1).apply(lambda x: _score_self(x))
    primer_df['score_endA']  = primer_df.endA.apply(lambda x: _score_endA(x))
    primer_df['score_enddG'] = primer_df.end_dg.apply(lambda x: _score_enddG(x))

    # ---- final cost (lower is better) : weighted sum of (1 - score) ----
    def primer_cost_row(row):
        scores = {
            'tm':    row['score_tm'],
            'gc':    row['score_gc'],
            'self':  row['score_self'],
            'endA':  row['score_endA'],
            'enddG': row['score_enddG'],
        }
        return _weighted_func_from_scores(scores)

    primer_df['efficiency'] = primer_df.apply(primer_cost_row, axis=1)
    primer_df.reset_index(inplace=True)
    primer_df.set_index(['start','stop','fr'], inplace=True)
    return primer_df


