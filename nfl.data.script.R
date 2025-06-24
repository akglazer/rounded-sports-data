library(tidyverse)
library(nflreadr) # Get NFL data
library(geosphere) # Get distances from lat/long

# Set season value
season_value <- 2023

## Note: the following is for offensive players
# change snap counts and player stats to get defense, ST

## Load snap counts for the season ##
snap_counts <- load_snap_counts(seasons = season_value)
# Get players that played in each game (those with nonzero snap counts)
players_in_games <- snap_counts %>%
  filter(offense_snaps > 0)

## Get weekly roster ##
weekly_roster <- load_rosters_weekly(season = season_value) 

## Get player IDs for matching ##
# Note: lots of missing IDs which makes this problematic to match on
# Ultimately will match on name instead
player_ids <- load_ff_playerids()

## Game/schedule data ##
games <- load_schedules(seasons = season_value)
# Create long format with 2 rows for each game
# Home games
home_games <- games %>%
  dplyr::select(game_id, season, game_type, week, gameday, weekday,
                team = home_team, opponent = away_team, 
                location, rest = home_rest,
                roof, surface, temp, wind, stadium, stadium_id) %>%
  mutate(is_home = 1)

# Away games
away_games <- games %>%
  dplyr::select(game_id, season, game_type, week, gameday, weekday,
                team = away_team, opponent = home_team, 
                location, rest = away_rest,
                roof, surface, temp, wind, stadium, stadium_id) %>%
  mutate(is_home = 0)

# Bind home and away games
games_long <- home_games %>%
  bind_rows(away_games)

## Weekly player stats for offensive players ##
# NOTE: only player-games where a player records 
# an official offensive stat will appear, 
# even if they played in that game, need to merge with snap counts
player_stats <- load_player_stats(seasons = season_value, 
                                  stat_type = "offense")

## Stadium coordinates ##
stadium_coords <- tibble::tribble(
  ~stadium_name, ~latitude, ~longitude,
  "MetLife Stadium", 40.8135, -74.0744,
  "Soldier Field", 41.8623, -87.6167,
  "Ford Field", 42.3400, -83.0456,
  "SoFi Stadium", 33.9534, -118.3392,
  "FirstEnergy Stadium", 41.5061, -81.6996,
  "NRG Stadium", 29.6847, -95.4107,
  "Lumen Field", 47.5952, -122.3316,
  "Paycor Stadium", 39.0954, -84.5161,
  "Lucas Oil Stadium", 39.7601, -86.1639,
  "AT&T Stadium", 32.7473, -97.0945,
  "State Farm Stadium", 33.5277, -112.2637,
  "M&T Bank Stadium", 39.2780, -76.6227,
  "Allegiant Stadium", 36.0909, -115.1838,
  "New Era Field", 42.7737, -78.7869,
  "Lambeau Field", 44.5013, -88.0622,
  "Gillette Stadium", 42.0909, -71.2643,
  "Mercedes-Benz Superdome", 29.9511, -90.0812,
  "Lincoln Financial Field", 39.9008, -75.1675,
  "FedExField", 38.9077, -76.8644,
  "GEHA Field at Arrowhead Stadium", 39.0489, -94.4849,
  "Raymond James Stadium", 27.9759, -82.5033,
  "Empower Field at Mile High", 39.7439, -105.0201,
  "Hard Rock Stadium", 25.9580, -80.2389,
  "Levi's Stadium", 37.4030, -121.9692,
  "Nissan Stadium", 36.1665, -86.7713,
  "U.S. Bank Stadium", 44.9740, -93.2581,
  "Bank of America Stadium", 35.2258, -80.8528,
  "Tottenham Stadium", 51.6043, -0.0665,
  "TIAA Bank Stadium", 30.3239, -81.6370,
  "Mercedes-Benz Stadium", 33.7554, -84.4008,
  "Acrisure Stadium", 40.4468, -80.0158,
  "Deutsche Bank Park", 50.0682, 8.6458,
  "Wembley Stadium", 51.5558, -0.2796
)

# Mapping team abbreviation to home stadium
team_home_stadiums <- tibble::tribble(
  ~team, ~home_stadium,
  "ARI", "State Farm Stadium",
  "ATL", "Mercedes-Benz Stadium",
  "BAL", "M&T Bank Stadium",
  "BUF", "New Era Field",
  "CAR", "Bank of America Stadium",
  "CHI", "Soldier Field",
  "CIN", "Paycor Stadium",
  "CLE", "FirstEnergy Stadium",
  "DAL", "AT&T Stadium",
  "DEN", "Empower Field at Mile High",
  "DET", "Ford Field",
  "GB",  "Lambeau Field",
  "HOU", "NRG Stadium",
  "IND", "Lucas Oil Stadium",
  "JAX", "TIAA Bank Stadium",
  "KC",  "GEHA Field at Arrowhead Stadium",
  "LA",  "SoFi Stadium",
  "LAC", "SoFi Stadium",
  "LV",  "Allegiant Stadium",
  "MIA", "Hard Rock Stadium",
  "MIN", "U.S. Bank Stadium",
  "NE",  "Gillette Stadium",
  "NO",  "Mercedes-Benz Superdome",
  "NYG", "MetLife Stadium",
  "NYJ", "MetLife Stadium",
  "PHI", "Lincoln Financial Field",
  "PIT", "Acrisure Stadium",
  "SEA", "Lumen Field",
  "SF",  "Levi's Stadium",
  "TB",  "Raymond James Stadium",
  "TEN", "Nissan Stadium",
  "WAS", "FedExField"
)

# join stadium coordinates to data
games_long <- games_long %>%
  left_join(stadium_coords, by = c("stadium" = "stadium_name"))

# Join home stadium coordinates
games_long <- games_long %>%
  left_join(team_home_stadiums, by = "team") %>%
  left_join(stadium_coords %>%
              rename(home_lat = latitude, home_lon = longitude),
            by = c("home_stadium" = "stadium_name"))

# Function to compute Haversine distance (in miles)
calculate_distance <- function(lat1, lon1, lat2, lon2) {
  # Ensure vectorized operations using ifelse()
  missing_values <- is.na(lat1) | is.na(lon1) | is.na(lat2) | is.na(lon2)
  ifelse(missing_values, 0, distHaversine(cbind(lon1, lat1), cbind(lon2, lat2)) / 1609.34)
}

# Compute travel distances between consecutive games for each team
games_long <- games_long %>%
  group_by(team) %>%
  arrange(team, gameday) %>%
  mutate(
    prev_lat = if_else(week == 1, home_lat, lag(latitude)),
    prev_lon = if_else(week == 1, home_lon, lag(longitude))
  ) %>%
  ungroup() %>%
  mutate(travel_distance = calculate_distance(prev_lat, prev_lon, 
                                              latitude, longitude))

# get offensive/defensive EPA per week and season
off_epa_week <- nfl_pbp %>%
  group_by(week, posteam) %>%
  filter(is.na(epa) == F, week <= 18, is.na(posteam) == F) %>%
  summarise(off_epa_per_play = mean(epa, na.rm = TRUE),
            off_plays = n()) %>%
  arrange(desc(off_epa_per_play))

def_epa_week <- nfl_pbp %>%
  group_by(week, defteam) %>%
  filter(is.na(epa) == F, week <= 18, is.na(defteam) == F) %>%
  summarise(def_epa_per_play = mean(epa, na.rm = TRUE),
            def_plays = n()) %>%
  arrange(def_epa_per_play)

rush_defense <- nfl_pbp %>%
  filter(play_type == "run", !is.na(epa)) %>%
  group_by(week, defteam) %>%
  summarise(
    rush_epa_allowed = mean(epa, na.rm = TRUE),
    #success_rate_allowed = mean(success, na.rm = TRUE),
    total_rush_plays = n()
  ) %>%
  arrange(defteam, week) %>%
  ungroup() %>%
  group_by(defteam) %>%
  mutate(cum_total_plays = cumsum(total_rush_plays),
            cum_weighted_epa = cumsum(rush_epa_allowed * total_rush_plays),
            cum_rush_epa_allowed = cum_weighted_epa / cum_total_plays) %>%
  select(week, defteam, cum_rush_epa_allowed)

# Create full grid of weeks and teams
all_weeks_teams <- expand_grid(
  week = 1:18,
  defteam = unique(rush_defense$defteam)
)

# Join and fill missing values
rush_defense_filled <- all_weeks_teams %>%
  left_join(rush_defense, by = c("week", "defteam")) %>%
  arrange(defteam, week) %>%
  group_by(defteam) %>%
  tidyr::fill(cum_rush_epa_allowed, .direction = "down") %>%
  ungroup()

epa_week <- off_epa_week %>%
  inner_join(def_epa_week, 
             by = c("posteam" = "defteam", "week")) %>%
  inner_join(rush_defense, 
             by = c("posteam" = "defteam", "week")) %>%
  filter(week <= 18) %>%
  rename(team = posteam)

odds <- games %>%
  dplyr::select(game_id, week, weekday, home_team,
                away_team, home_score, away_score, 
                home_moneyline, away_moneyline, spread_line)

# write.csv(epa_week, "~/Downloads/epa_week.csv")

## Join relevant offensive player stats and game data ##
nfl_rb_df <- weekly_roster %>%
  dplyr::select(season, week, team, full_name, 
                birth_date, status, position) %>%
  mutate(full_name = clean_player_names(full_name)) %>%
  filter(position == "RB") %>%
  left_join(players_in_games %>%
              mutate(player = clean_player_names(player)),
            by = c("season", "week", "team", "full_name" = "player")) %>%
  left_join(player_stats %>%
  dplyr::select(player_id, player_display_name, position,
                season, week, recent_team,
                rushing_yards, receiving_yards, passing_yards) %>%
    mutate(player_display_name = clean_player_names(player_display_name)), 
                by = c("full_name" = "player_display_name",
                                       "season", "week",
                                       "team" = "recent_team")) %>%
  left_join(
    games_long,
    by = c("season", "week", "team" = "team")
  ) %>%
  mutate(rushing_yards = ifelse(is.na(rushing_yards), 0, rushing_yards),
         receiving_yards = ifelse(is.na(receiving_yards), 0, receiving_yards),
         passing_yards = ifelse(is.na(passing_yards), 0, passing_yards),
         played = ifelse(is.na(offense_snaps), 0, 1),
         age = round(time_length(interval(birth_date, gameday), "years"), 2)) %>%
  left_join(rush_defense_filled %>% mutate(next_week = week + 1), 
            by = c("week" = "next_week", "opponent.y" = "defteam")) %>%
  dplyr::select(season, week, team, full_name, player_id,
                birth_date, age, status, 
                game_id = game_id.x, game_type = game_type.x, 
                position = position.x, played,
                opponent = opponent.x, offense_snaps, offense_pct,
                rushing_yards, receiving_yards, passing_yards,
                opponent_rush_epa_allowed = cum_rush_epa_allowed,
                gameday, weekday, location, rest, roof, surface,
                temp, wind, stadium, stadium_id, is_home, 
                latitude, longitude, prev_lat, prev_lon, travel_distance) 

# write csv
#write.csv(nfl_rb_df, "~/Downloads/nfl-rb-game-yds.csv")

## Get play-by-play data
nfl_pbp <- nflreadr::load_pbp(seasons = season_value)
# Just get rushing yards
nfl_rb_rush_df <- nfl_pbp %>%
  filter(rush_attempt == 1) %>%
  dplyr::select(game_id, game_date, season, week, 
                play_id, team = posteam,
                opponent = defteam,
                score_differential,
                qtr,
                game_seconds_remaining,
                rusher_player_id, rusher_id,
                rusher_player_name, rushing_yards, down, ydstogo, 
                yardline_100, 
                yards_gained) %>%
  left_join(nfl_rb_df %>%
              select(player_id, season, week, team,
                     age, opponent_rush_epa_allowed, 
                     gameday, weekday, location, rest, 
                     roof, surface, temp, wind, stadium, 
                     stadium_id, is_home, latitude,
                     longitude, travel_distance),
            by = c("rusher_player_id" = "player_id", "season", 
                   "week", "team")) %>%
  filter(is.na(age) == F)

# write csv
#write.csv(nfl_rb_rush_df, "~/Downloads/nfl-rb-pbp-rushing.csv")
