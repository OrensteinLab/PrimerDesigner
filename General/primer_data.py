import pandas as pd
import primer3 as p3
from Bio.Seq import Seq
from Bio.SeqUtils import GC, seq1, seq3
from General.constants import *


def revcomp(seq):
  return str(Seq(seq).reverse_complement())

def subsequences(sequence,primer_lmin,primer_lmax): #Generates all subsequences w/ all poss. start-stop pairs
  ls = []
  for j in range(primer_lmin, primer_lmax+1): #length
    for i in range(len(sequence)-j+1): #starting index
      start = i
      stop = i+j
      ls.append([sequence[start:stop], start, stop, stop-start])
  return ls

def create_primer_df(sequence_nt,args):
  # sets up pcr reaction with parameters
  pcr = p3.thermoanalysis.ThermoAnalysis(dna_conc=250,
                                         mv_conc=50,
                                         dv_conc=0,
                                         dntp_conc=0,
                                         tm_method='santalucia',
                                         salt_correction_method='owczarzy',
                                         temp_c=25)

  # convention: start index of r-primers will be 3' (i.e. start < stop)
  primer_f = pd.DataFrame(columns=['seq','start','stop','fr','len'])
  primer_f[['seq','start','stop','len']] = subsequences(sequence_nt, args.primer_lmin, args.primer_lmax)
  primer_f['fr'] = 'f'

  #Shifting so that 0 is at the start of mutreg (upstream has negative values)
  primer_f['start'] = primer_f.start - len(upstream_nt)
  primer_f['stop'] = primer_f.stop - len(upstream_nt)

  #Creating reverse primers at same locations
  primer_r = primer_f[['seq','start','stop','fr','len']].copy()
  primer_r['fr'] = 'r'
  primer_r['seq'] = primer_r.seq.apply(revcomp)

  #Concatenating Forward & Reverse
  primer_df = pd.concat([primer_f,primer_r])
  primer_df.sort_values(by=['start','stop','fr'], inplace=True)

  #Calculating "Cost" Values
  primer_df['gc'] = primer_df.seq.apply(GC)
  primer_df['tm'] = primer_df.seq.apply(pcr.calc_tm)
  res = primer_df.seq.apply(lambda s: pcr.calc_hairpin(s).todict())
  primer_df['hp_tm'] = res.apply(lambda res: res['tm'])
  primer_df['hp_dg'] = res.apply(lambda res: res['dg']*1e-3)
  res = primer_df.seq.apply(lambda s: pcr.calc_homodimer(s).todict())
  primer_df['ho_tm'] = res.apply(lambda res: res['tm'])
  primer_df['ho_dg'] = res.apply(lambda res: res['dg']*1e-3)

  def primer_cost(primer):
    # calculates primer cost based on homodimer an haipin delta G, tm cost, and len cost
    hp_dg_max = -4
    ho_dg_max = -8
    tm_min = 58

    tm_cost = max(0, tm_min-primer.tm)**1.7
    hp_cost = max(0, hp_dg_max - primer.hp_dg)**1.2
    ho_cost = max(0, ho_dg_max - primer.ho_dg)**1.2
    len_cost = primer.len*1e-5

    cost = hp_cost + ho_cost + len_cost + tm_cost

    return cost

  primer_df['cost'] = primer_df.apply(primer_cost, axis=1)

  primer_df.reset_index(inplace=True)
  primer_f = primer_df.query('fr=="f"').reset_index(drop=True)

  primer_df.set_index(['start','stop','fr'], inplace=True)

  return primer_f, primer_df

def check_threshold(tm,gc,args):
  # returns false only if the threshold flag is on and primers did not pass threshold
  if not args.apply_threshold:
    return True
  else:
    if args.min_gc <= gc <= args.max_gc and  args.min_tm <= tm <=args.max_tm:
      return True
    else:
      return False