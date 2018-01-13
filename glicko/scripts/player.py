import pandas as pd
import math
import numpy as np

class Player(object):
    
    # need to maintain only historical rankings, not phi etc.
    hist_mu = None
    hist_phi = None
    
    def __init__(self, player_id, periods, race_ids, glicko):
        
        assert isinstance(player_id, str), "Player ID must be a string"
        self.pid = player_id
        self.hist_mu = pd.DataFrame(
            columns=race_ids, index=periods, dtype=float)
        self.hist_phi = pd.DataFrame(
            columns=race_ids, index=periods, dtype=float)
        self.rankings = {}
        self.new_rankings = {}
        self.glicko = glicko
              
    def init_rating(self, race_id, method="median"):
        
        if (method == "default") or (len(self.rankings) < 2):
            
            self.rankings[race_id] = self.glicko.create_rating()
            
        elif method == "median":
            # other rankings. Possibly exclude stunties? 
            # divide PHI by the number of other ranks?
            # Or always start stunties at 1500?
            mu_vals = [v.mu for v in self.rankings.values()]
            phi_vals = [v.phi for v in self.rankings.values()]
            _mu = np.median(mu_vals)
            _phi = np.max(phi_vals)

            self.rankings[race_id] = self.glicko.create_rating(
                mu=_mu, phi=_phi)

        return self.rankings[race_id]
    
    def end_ranking_period(self, date):
        
        # copy new rankings
        for race in self.hist_mu.columns:
            if (race in self.rankings) and (race not in self.new_rankings):
                self.rankings[race] = self.decay(race)
            elif race in self.new_rankings:
                self.rankings[race] = self.new_rankings[race]
                
        # delete temp ranks
        self.new_rankings = {}
        
        # fill historical with rankings.
        for rid, rank in self.rankings.items():
            self.hist_mu.at[date, rid] = rank.mu
            self.hist_phi.at[date, rid] = rank.phi
            
    def decay(self, race):
        
        # transform into glicko space...
        dnp = self.glicko.scale_down(self.rankings[race])
        
        # increment phi
        phi_star = math.sqrt(min(self.glicko.phi, dnp.phi ** 2 + dnp.sigma ** 2))
                
        return self.glicko.scale_up(
            self.glicko.create_rating(mu=dnp.mu, phi=phi_star, sigma=dnp.sigma))
