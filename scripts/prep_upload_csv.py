# To handle inactive coaches, give them a "-" instead of a rank. As long as rank not index, no problem.
# Read from config.
fn = snakemake.input.csv
outfn = snakemake.output.csv

import pandas as pd
df = pd.read_csv(fn, index_col=0)

rank = 1
rank_list = []
for phi in df.phi:
  if phi > snakemake.params.phi_limit:
    rank_list.append("-")
  else:
    rank_list.append(str(rank))
    rank += 1

df["rank"] = rank_list

dfq = df[["rank", "coach", "naf_number", "race", "nation", "mu", "phi", "value"]]

dfq.to_csv(outfn, index=False, header=True, float_format="%.1f")
