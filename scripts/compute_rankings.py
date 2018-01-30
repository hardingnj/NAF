import h5py
import pandas as pd
def cleanup(x):         
    if isinstance(x, float):
        return "Unknown"
    y = x.lower()
    y = y.rstrip()
    return y.capitalize()


PHI_PENALTY = snakemake.params.phi_penalty
coach_info = pd.read_csv(snakemake.input.csv)
coach_info["nation"] = coach_info.nation.apply(cleanup)

with h5py.File(snakemake.input.hdf5, "r") as hf:
    
    coaches = list(hf["coaches"].keys())
    race_ids = hf["race_ids"][:].astype("<U16")
    
    df_mu = pd.DataFrame(index=coaches, columns=race_ids)
    df_mu.index.name = "coach"
    df_phi = pd.DataFrame(index=coaches, columns=race_ids)
    df_phi.index.name = "coach"

    for pid in coaches:
        df_mu.loc[pid] = hf["coaches"][pid]["mu"][0]
        df_phi.loc[pid] = hf["coaches"][pid]["phi"][0]
    mu_melt = df_mu.reset_index().melt(var_name="race", id_vars="coach", value_name="mu")
    phi_melt = df_phi.reset_index().melt(var_name="race", id_vars="coach", value_name="phi")

    rank_df = pd.merge(mu_melt, phi_melt, on=["coach", "race"])
    rank_df["value"] = (rank_df.mu - (PHI_PENALTY * rank_df.phi))

    # merge with existing coach info.
    merged = pd.merge(rank_df, coach_info, left_on=["coach", "race"], right_on=["naf_name", "race"])
    merged = merged.dropna(subset=["value"]).sort_values("value", ascending=False)
    merged.reset_index(inplace=True, drop=True)
    merged.index = merged.index.values + 1
    merged.to_csv(snakemake.output.csv, float_format='%.3f')
