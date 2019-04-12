"""
This script cleans the NAF game data, removing where there are missing values or coaches do not 
meet naming requirements
"""

import pandas as pd
import re
import numpy as np

# define dtypes
datatypes = {
    'tournament_id': int,
    'match_id': int,
    'tournament_name': "|S80", 
    'home_score': int,
    'away_score': int,
    'home_cas': int,
    'away_cas': int,
    'mirror': bool,
    'home_tr': int,
    'away_tr': int,
    'swiss': bool}

naf_data = pd.read_csv(
    snakemake.input.csv, sep=",", index_col=0, dtype=datatypes, parse_dates=[0])

for col in naf_data.columns:
    print(col, naf_data[col].dtype)

#naf_data["date"] = pd.to_datetime(naf_data.date, format="%Y/%m/%d")

#naf_data.set_index("date", inplace=True)
naf_data.sort_index(inplace=True)

allowed_name = re.compile("\w+")

keep = np.ones(naf_data.shape[0], dtype="bool")

# Fix for Marco
for f in "home_coach", "away_coach":
    naf_data[f] = naf_data[f].str.replace("badstorm", "BadStorm")

# VALIDATION
excluded_ids = []

for i, (hid, aid) in enumerate(zip(naf_data.home_coach.values, naf_data.away_coach.values)):
    for xid in (hid, aid):
        if not isinstance(xid, str):
            if xid not in excluded_ids:
                print(xid)
                excluded_ids.append(xid)
            keep[i] = False
            
        else:
            result = allowed_name.match(xid)
            if not result:
                if xid not in excluded_ids:
                    print(xid)
                    excluded_ids.append(xid)
                keep[i] = False

keep = keep & (naf_data.variant == 'Blood Bowl').values

print(naf_data.shape[0], "games")
naf_data = naf_data.loc[keep]
print(naf_data.shape[0], "games after qc")

cols_home = ["home_coach", "home_race", "home_score"]
cols_away = ["away_coach", "away_race", "away_score"]

tmp = naf_data[cols_home + cols_away].copy()
tmp2 = naf_data[cols_away + cols_home].copy()
tmp2.columns = tmp.columns

rank_data = pd.concat([tmp, tmp2])

coach_names = rank_data.home_coach.unique()
diff = (rank_data.home_score - rank_data.away_score)
win = diff > 0
draw = diff == 0
rank_data["result"] = 0 + (0.5 * draw) + (1.0 * win)

rank_data.to_csv(snakemake.output.txt, sep="\t", compression="gzip")
