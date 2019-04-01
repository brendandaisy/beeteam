source('beeteam.R')

### to run the model with the 'current' parameter settings, use quick_run
quick_run()

### the function plot_time_series has several arguments for visualizing a run

plot_time_series() # default run is a quick_run
plot_time_series(pop=c('B'))

### to run with different parameters, can call ode explicitely
### EX: naive recreation of the HI experiment
params_hi = c(p = 0.17,
               q = 0.3,
               d = 1/7,
               z = 0, # equiv. to reinnoculating flowers each day
               Nb = 45,
               Nh = 0, # equiv. to having no HB present
               Nf = 30)

state_hi = c(B = 0, H = 0, FL = 30)

hi_out = ode(y=state_hi, times=0:3, func=beeteam, parms=params_hi)

plot_time_series(ode=hi_out, pop=c('B', 'FL')) +
    geom_hline(yintercept=15, linetype='dashed', col=grey(.5)) # output matches observed trend of 1/3 BB infected after 3 days


### another way to change parameters is to just change the default settings
### EX2: effect of increasing Nf

params['Nf'] = 700 # add more flowers
state['H'] = 50 # start with more infections for aesthetic reasons

plot_time_series() # now, infection is declining


