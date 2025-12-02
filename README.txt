R Code for Glazer et al. (2026) paper in The American Statistician on modeling rounded NFL yardage data

single.script.R:  Fits model to single player data and creates inference results and graphics
	single.rd.norm.R:  MCMC algorithm for single player model with rounded normal 
	single.rd.ald.R:  MCMC algorithm for single player model with rounded ALD 
	single.rd.tald.R:  MCMC algorithm for single player model with truncated rounded ALD (in appendix)
	single.ppc.R:  Model checking code for posterior predictive checks using skewness and kurtosis 

multi.script.R:  Fits model to multi player data and creates inference results and graphics
	multi.rd.norm.R:  MCMC algorithm for multi player model with rounded normal 
	multi.rd.ald.R:  MCMC algorithm for multi player model with rounded ALD 

nfl.data.script.R: Access and clean NFL data to generate nfl-rb-pbp-rushing.csv

nfl-rb-pbp-rushing.csv:  Data set based on NFL running backs for 2023 season.  Each script above performs additional subsetting of the data to restrict it to regular season and omit missing values.   

Note:  The multiplayer MCMC algorithms require hours to fit to this data set.   

