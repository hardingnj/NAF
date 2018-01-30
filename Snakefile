configfile: "config.yaml"

rule all:
  input:
    "rankings.html"

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

rule dl_gdata:
  output:
    "data/all_matches.csv.gz"
  shell:
    "curl -L www.dropbox.com/s/t8cq8zixe99fl8f/all_matches.csv?dl=1 | gzip > {output}"

rule dl_pdata:
  output:
    "data/all_coaches.csv.gz"
  shell:
    "curl -L https://www.dropbox.com/s/xus9dhoytlhofqs/all_coaches.csv?dl=1 | gzip > {output}"

rule process_data:
  input:
    csv="data/all_matches.csv.gz"
  output:
    txt="output/cleaned_games.txt.gz"
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
    csv="output/player_ranks.csv"
  params:
    phi_penalty=config["phi_penalty"]
  script:
    "scripts/compute_rankings.py"

rule makehtml:
  input:
    txt=rules.get_ranks.output.csv,
    nb="produce_html.ipynb"
  output:
    "rankings.html",
    temp("produce_html.nbconvert.ipynb")
  shell:
    "jupyter nbconvert --to notebook --execute {input.nb}"

