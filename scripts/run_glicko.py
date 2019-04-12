import pandas as pd
import glicko2 as Glicko
import numpy as np 
import h5py
from player import Player

# improbable results are likely so set tau low
TAU = snakemake.params.tau      # system constant
MU = snakemake.params.mu
PHI = snakemake.params.phi      # Starting rating deviation
SIGMA = snakemake.params.sigma  # starting volatility
UPDATE_FREQ = snakemake.params.update_freq

glck = Glicko.Glicko2(mu=MU, tau=TAU, phi=PHI, sigma=SIGMA)

rank_data = pd.read_table(snakemake.input.txt, index_col=0, parse_dates=[0])
rank_data.sort_index(inplace=True)

uniq_races = rank_data.home_race.unique().tolist()

today = pd.Timestamp.today()
cutoff = pd.Timestamp(year=today.year, month=today.month, day=1) - pd.Timedelta('1 days')
print(cutoff)

# Trim data to cutoff
rank_data = rank_data[:cutoff]

grouped_games = rank_data.groupby(pd.Grouper(freq=UPDATE_FREQ))

rank_periods = [p for p, _ in grouped_games]
ranking_data = dict()

for period, x in grouped_games:
    if period.month == 1:
        print(period)

    coaches = x.home_coach.unique()
    for ix, xid in enumerate(coaches):
        if xid not in ranking_data:
            ranking_data[xid] = Player(xid, rank_periods, uniq_races, glck)
    
    # group by player
    grped = x.groupby(["home_coach", "home_race"])
    
    # first run through all players who have played in this period
    for (player, race), data in grped:
        
        player_rank = ranking_data[player]
        
        # if new race
        if race not in player_rank.rankings:
            player_rank.init_rating(race)
            
        series = list()

        for opp_id, opp_race, result in zip(
            data.away_coach, data.away_race, data.result):
            
            opp_rank = ranking_data[opp_id]
            
            # opponent ranking?
            if opp_race not in opp_rank.rankings:
                opp_rank.init_rating(opp_race)
            
            series.append(
                (result, opp_rank.rankings[opp_race]))
        
        player_rank.new_rankings[race] = glck.rate(
            player_rank.rankings[race], series)
        
    # end all ranking periods.
    for k, v in ranking_data.items():
        v.end_ranking_period(period)

# Save all historical data.
datestrings = np.array([x.strftime('%Y-%m-%d') for x in rank_periods[::-1]], dtype="|S10")

with h5py.File(snakemake.output.hdf5, "w") as hf:
    hf.create_dataset("date", data=datestrings)
    hf.create_dataset("race_ids", data=np.array(uniq_races, dtype="|S16"))
    player_g = hf.create_group("coaches")

    for key, value in ranking_data.items():
        g = player_g.create_group(key)
        g.create_dataset("mu", data=value.hist_mu[::-1].round(3).values, compression="gzip")
        g.create_dataset("phi", data=value.hist_phi[::-1].round(3).values, compression="gzip")
