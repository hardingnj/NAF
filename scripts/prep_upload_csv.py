import pandas as pd
import numpy as np

# Read from config.
fn = snakemake.input.csv
phi_cut = snakemake.params.phi_active

df = pd.read_csv(fn, index_col=0).drop(["rank", "race_rank"], axis=1)
df = df.loc[df.curr_phi < snakemake.params.phi_limit]

rank = 1
url_path = "https://member.thenaf.net/index.php?module=NAF&type=tournamentinfo&uid={nafnum}||{value}"

df["rating"] = df.curr_rating.copy()
df.loc[df.curr_phi > 100, "rating"] = np.nan

df["rank"] = df.rating.rank(ascending=False, na_option="bottom", method="min").astype("int")
df["race_rank"] = df.groupby("race").rating.transform(
  lambda y: y.rank(ascending=False, na_option="bottom", method="min")).astype("int")

# rank change
df["old_rating"] = df.last_rating.copy()
df.loc[df.last_phi > 100, "old_rating"] = np.nan
df["old_rank"] = df.old_rating.rank(ascending=False, na_option="bottom", method="min").astype("int")

df["rank_change"] = df["old_rank"] - df["rank"]

# Need to set ranks for global and race to 999999. if phi > X
df.loc[df.curr_phi > phi_cut, "rank"] = 999999
df.loc[df.curr_phi > phi_cut, "race_rank"] = 999999
df.loc[(df.curr_phi > phi_cut) | (df.last_phi > phi_cut), "rank_change"] = 0

df["mu"] = df.curr_mu
df["phi"] = df.curr_phi

df["coachname"] = df.apply(lambda y: url_path.format(nafnum=y.naf_number, value=y.coach), axis=1)
dfq = df[["rank", "race_rank", "coachname", "naf_number", "race", "nation", "mu", "phi", "rating"]]

dfq.sort_values("rank", ascending=True).to_csv(snakemake.output.upload, index=False, header=True, float_format="%.1f")

# top with each race
dfq.query("race_rank == 1").to_csv(snakemake.output.races, index=False, header=True, float_format="%.1f")

# top movers
rc = df.sort_values("rank_change")[["rank", "race_rank", "rank_change", "naf_name", "naf_number", "race", "mu", "phi", "rating"]]
rc.head(10).to_csv(snakemake.output.losers, index=False)
rc.tail(10)[::-1].to_csv(snakemake.output.winners, index=False)
