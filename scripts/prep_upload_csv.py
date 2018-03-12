# To handle inactive coaches, give them a "-" instead of a rank. As long as rank not index, no problem.
# Read from config.
fn = snakemake.input.csv
outfn = snakemake.output.csv

import pandas as pd
import numpy as np

df = pd.read_csv(fn, index_col=0)

rank = 1
url_path = "https://member.thenaf.net/index.php?module=NAF&type=tournamentinfo&uid={nafnum}||{value}"

rank_list = []
for phi in df.phi:
  if phi > snakemake.params.phi_limit:
    rank_list.append(999999)
  else:
    rank_list.append(str(rank))
    rank += 1

df["active"] = np.where(df.phi > 100, "Inactive", "Active") 
df["rank"] = rank_list

df["coachname"] = df.apply(lambda y: url_path.format(nafnum=y.naf_number, value=y.coach), axis=1)
df["rating"] = df["value"]

dfq = df[["rank", "coachname", "naf_number", "race", "nation", "mu", "phi", "rating", "active"]]

dfq.to_csv(outfn, index=False, header=True, float_format="%.1f")
