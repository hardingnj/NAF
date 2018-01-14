configfile: "config.yaml"

rule all:
  input:
    "output/player_ranks.txt"

rule process_data:
  input:
    csv="../data/all_matches.csv.gz"
  output:
    txt="processed/cleaned_games.txt.gz"
  script:
    "scripts/process_naf_data.py"

rule compute_rankings:
  input:
    txt=rules.process_data.output.txt,
  output:
    hdf5="output/rankings.h5"
  params:
    mu=config["mu"],
    phi=config["phi"],
    tau=config["tau"],
    sigma=config["sigma"],
    update_freq=config["update_freq"]
  script:
    "scripts/run_glicko.py"

rule get_ranks:
  input:
    hdf5=rules.compute_rankings.output.hdf5
  output:
    txt="output/player_ranks.txt"
  params:
    phi_penalty=config["phi_penalty"]
  script:
    "scripts/compute_rankings.py"
