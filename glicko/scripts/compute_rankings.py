import h5py
import pandas as pd

PHI_PENALTY = snakemake.params.phi_penalty

with h5py.File(snakemake.input.hdf5, "r") as hf:
    
    coaches = list(hf["coaches"].keys())
    race_ids = hf["race_ids"][:].astype("<U16")
    
    df = pd.DataFrame(index=coaches, columns=race_ids)
    df.index.name = "coach"

    for pid in coaches:
        df.loc[pid] = hf["coaches"][pid]["mu"][0] * (PHI_PENALTY * hf["coaches"][pid]["phi"][0])
    rank_df = df.reset_index().melt(var_name="race", id_vars="coach").sort_values("value", ascending=False)
    rank_df.dropna(subset=["value"]).to_csv(snakemake.output.txt, sep="\t")
