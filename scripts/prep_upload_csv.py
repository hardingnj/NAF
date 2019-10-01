import pandas as pd
import numpy as np

# Read from config.
fn = snakemake.input.csv
phi_cut = snakemake.params.phi_active

df = pd.read_csv(fn, index_col=0)

# drop records where phi has expired.
df = df.loc[df.curr_phi < snakemake.params.phi_limit]

rank = 1
url_path = "https://member.thenaf.net/index.php?module=NAF&type=tournamentinfo&uid={nafnum}||{value}"

# Used to compute current ranks
df["rating"] = df.curr_rating.copy()
df.loc[df.curr_phi > 100, "rating"] = np.nan

df["old_rating"] = df.last_rating.copy()
df.loc[df.last_phi > 100, "old_rating"] = np.nan

# Compute_rank change
df["decay"] = (df.last_phi - df.curr_phi) * snakemake.params.phi_penalty
df["change"] = df.curr_mu - df.last_mu

# Global ranks.
df["rank"] = df.rating.rank(ascending=False, na_option="bottom", method="min").astype("int")
df["old_rank"] = df.old_rating.rank(ascending=False, na_option="bottom", method="min").astype("int")
df.loc[~(df.curr_phi <= phi_cut), "rank"] = 999999
df.loc[~(df.last_phi <= phi_cut), "old_rank"] = 999999

rank_str = "{0} ({1})"
df["rank_rep"] = df[["rank", "old_rank"]].apply(
  lambda y: rank_str.format(y["rank"], y["old_rank"]).replace("999999", "-"), axis=1)

# Compute race ranks
df["race_rank"] = df.groupby("race").rating.transform(
  lambda y: y.rank(ascending=False, na_option="bottom", method="min")).astype("int")
df["old_race_rank"] = df.groupby("race").old_rating.transform(
  lambda y: y.rank(ascending=False, na_option="bottom", method="min")).astype("int")

# Need to set ranks for global and race to 999999. if phi > X
df.loc[~(df.curr_phi <= phi_cut), "race_rank"] = 999999
df.loc[~(df.last_phi <= phi_cut), "old_race_rank"] = 999999
df["race_rank_rep"] = df[["race_rank", "old_race_rank"]].apply(
  lambda y: rank_str.format(y["race_rank"], y["old_race_rank"]).replace("999999", "-"), axis=1)

# top movers
# if either set rank change to 0.
df["rank_change"] = df["old_rank"] - df["rank"]
df.loc[(df.curr_phi > phi_cut) | (df.last_phi > phi_cut), "rank_change"] = 0

# rename variables
df["mu"] = df.curr_mu
df["phi"] = df.curr_phi
df["rating"] = df.curr_rating
df["qrank"] = df["rank"]
df["rank"] = df["rank_rep"]
df["qrace_rank"] = df["race_rank"]
df["race_rank"] = df["race_rank_rep"]

#if naf_number is missing... 
df["naf_number"] = df.naf_number.fillna(0).astype("int")

# create URL
df["coachname"] = df.apply(lambda y: url_path.format(nafnum=y.naf_number, value=y.coach), axis=1)
core_cols = ["rank", "race_rank", "coachname", "naf_number", "race", "nation", "rating"]
extra_cols = snakemake.params.extra_cols
printcols = core_cols + extra_cols

if snakemake.params.globalmode:
    dfq = df.sort_values("rating", ascending=False)[printcols[2:]]
    g = dfq.groupby("nation")
    for nat, _df in g:
      _df.reset_index(drop=True).to_csv(snakemake.output.upload.replace("csv", nat + ".csv"), index=True, header=True, float_format="%.2f", quoting=2)
    dfq.to_csv(snakemake.output.upload, index=False, header=True, float_format="%.1f", quoting=2)

else:
    dfq = df.sort_values("qrank", ascending=True)[printcols]
    dfq.to_csv(snakemake.output.upload, index=False, header=True, float_format="%.1f", quoting=2)

# top with each race
df.query("qrace_rank == 1")[printcols].to_csv(snakemake.output.races, index=False, header=True, float_format="%.1f")

# print winners
rc = df.sort_values("rank_change")[["rank", "race_rank", "rank_change", "naf_name", "naf_number", "race", "mu", "phi", "rating"]]
rc.head(10).to_csv(snakemake.output.losers, index=False)
rc.tail(10)[::-1].to_csv(snakemake.output.winners, index=False)
