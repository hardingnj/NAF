# To handle inactive coaches, give them a "-" instead of a rank. As long as rank not index, no problem.
# Read from config.
fn = snakemake.input.csv
outfn = snakemake.output.csv

import pandas as pd
import numpy as np
from collections import Counter

counter = Counter()

df = pd.read_csv(fn, index_col=0)
# apply phi cutoff
df = df.loc[df.phi < snakemake.params.phi_limit]

rank = 1
url_path = "https://member.thenaf.net/index.php?module=NAF&type=tournamentinfo&uid={nafnum}||{value}"

rank_list = []
race_rank_list = []

for phi, raceid in zip(df.phi.values, df.race.values):
  if phi > snakemake.params.phi_active:
    rank_list.append(999999)
    race_rank_list.append(999999)
  else:
    rank_list.append(rank)
    rank += 1
    
    counter.update([raceid])
    race_rank_list.append(counter[raceid])

df["rank"] = rank_list
df["racerank"] = race_rank_list

df["coachname"] = df.apply(lambda y: url_path.format(nafnum=y.naf_number, value=y.coach), axis=1)
df["rating"] = df["value"]

dfq = df[["rank", "racerank", "coachname", "naf_number", "race", "nation", "mu", "phi", "rating"]]

dfq.to_csv(outfn, index=False, header=True, float_format="%.1f")
