require(deSolve)
require(tidyverse)
require(ggthemes)
require(wesanderson)

## core model: bumblebees and honeybees share transmission route via flowers

# number of days to run the model
maxT = 100

# model parameters
params = c(p = 0.17, # intrinsic flower to bee infection rate
           q = 0.3, # intrinsic bee to flower infection rate
           d = 1/7, # death rate (lifespan = 1/d days)
           z = 1/10, # flower recovery/death
           Nb = 100, # number bumble bees
           Nh = 100, # num honey bees
           Nf = 100) # num flowers

# model start state. by default, a single honeybee
state = c(B = 0, # B := num infected bumblebees
          H = 1, # H := num infected honeybees
          FL = 0) # FL := num infected flowers

# timesteps for which model output should be recorded
times = seq(0, maxT, by=1)

beeteam = function(t, state, parameters) {
  # model difference equations
  with(as.list(c(state, parameters)), {
    dB = p*(Nb - B)*FL/Nf - d*B
    dH = p*(Nh - H)*FL/Nf - d*H
    dFL = q*((Nf - FL)*(H + B))/Nf - z*FL
    return(list(c(dB, dH, dFL)))
  })
}

quick_run = function() {
  # run the ode with "default" values
  # wrapper to ease notation and ensure above global values are used
  ode(y=state, times=times, func=beeteam, parms=params)
}

## additional functions for aiding analysis

h_fixed_point = function(nf = params[7], nh = params[6], nb = params[5]) {
  # calculate the fixed point for honey bees
  with(as.list(params), {
    num = (nh + nb) * p * q - d * nf * z
    denom = (1 + nb / nh) * (d + p) * q
    return(pmax(num / denom, 0))
  })
}

b_fixed_point = function(nh = params[6], nb=params[5], nf=params[7]) {
  # calculate the fixed point for bumble bees
  with(as.list(params), {
    return(pmax(nb / nh * h_fixed_point(nf=nf, nh=nh, nb=nb), 0))
  })
}

fl_fixed_point = function(nh = params[6], nb=params[5], nf=params[7]) {
  # calculate the fixed point for flowers
  with(as.list(params), {
    num = nf * ((nh + nb) * p * q - d * nf * z)
    denom = p * ((nh + nb) * q + nf * z)
    return(pmax(num / denom, 0))
  })
}

h_crit_point = function(n) {
  # finds the critical value of H, past which the infection will survive
  with(as.list(params), {
    num = d * z * Nf
    s = sum(bb(0:n))
    denom = s * p * q
    return(num / denom)
  })
}

bb = function(i) {
  # compute the ith term in the series for the reproductive number
  with(as.list(params), {
    num = p * q * Nb
    denom = z * Nf * d
    return((num / denom)^i)
  })
}

repr_number = function(n=20) {
  # calculate the nth approximation for the reproductive number R_0
  # increasing n results in slight increase in R_0 (see definition in writeup)
  # assumes the perspective of the honeybee
  hb = (fl_visit^2 * bee_inf * fl_inf * num_honeybee) /
    (bee_mort * fl_rec * num_flower)
  bb = sum(sapply(0:n, bb))
  return(hb * bb)
}

## model visualizations: several helper functions for viewing model with ggplot

# default colors for each population
color_pal = wes_palette(name = "Darjeeling1", 3)
names(color_pal) = c('B', 'H', 'FL')
color_pal2 = wes_palette(name = "Darjeeling1", 3)
names(color_pal2) = c('Nb', 'Nh', 'Nf')

plot_time_series = function(ode = quick_run(),
                            pop = c('B', 'H', 'FL'),
                            colors = color_pal) {
  # plot population infection levels as function of time
  # input: 
    # ode: a matrix of class deSolve
    # pop: columns of ode to select, representing infection levels
    # colors: a color pallette. Must have names corresponding to names in 'pop'
  # output:
    # a ggplot object
    as.data.frame(ode) %>%
    select(pop, time) %>%
    gather(Population, inf, -time) %>%
    ggplot(aes(x=time, y=inf, col=Population)) +
      geom_line(size=1.2) +
      labs(y=paste(str_c(pop, collapse=", "), "Infections"), x="Time (Days)") +
      theme_few() + scale_color_manual(values=colors)
}

plot_state_space = function(ode = quick_run()) {
    ggplot(as_tibble(ode), aes(x=log(H+1), y=log(B+1))) +
        geom_line() +
        theme_few() +
        scale_color_manual(values=colors)
}

plot_fixed_points = function(n = 1:100,
                             pop = c('Nb', 'Nh', 'Nf'),
                             func = b_fixed_point,
                             colors = color_pal2) {
  # plot a fixed points function as a function of provided population levels,
  # with one line for each provided population level
  # input: 
    # n: the range of population levels to vary
    # pop: population levels
    # colors: a color pallette. Must have names corresponding to names in 'pop'
  # output:
    # a ggplot object
    val = c()
    for (p in pop) {
        if (p == 'Nb') val = c(val, func(nb=n))
        if (p == 'Nh') val = c(val, func(nh=n))
        if (p == 'Nf') val = c(val, func(nf=n))
    }
    tibble(Population =  rep(pop, each=length(n)),
           n = rep(n, length(pop)),                
           val = val) %>%
        ggplot(aes(x=n, y=val, col=Population)) +
        geom_line(size=1.2) +
        labs(## title=str_to_title(str_replace_all(deparse(substitute(func)), "_", " ")),
            x="Population Size",
            y=bquote(.(str_to_upper(str_sub(deparse(substitute(func)), 1, 1)))^'*')) +
        theme_few() + scale_color_manual(values=colors)
}

## plot_2d_state_space = function(ode = quick_run(),
##                                pop = c('B', 'H'),
##                                colors = color_pal) {
##   # plot population infection levels as function of time
##   # input: 
##   # ode: a matrix of class deSolve
##   # pop: columns of ode to select, representing infection levels
##   # colors: a color pallette. Must have names corresponding to names in 'pop'
##   # output:
##   # a ggplot object
##   if (length(pop) != 2) {
##     stop('only 2 dimensional state space supported; pop must be of length 2')
##   }
##   as.data.frame(ode) %>%
##     select(pop) %>%
##     ggplot(aes(x=eval(pop[1]), y=eval(pop[2]))) +
##       geom_line(size=1.2) +
##       labs() +
##       theme_few()
## }


## additional simple SIS model for intercolony dynamics between bees

params_colony = c(b = .18,
                  d = 1/7,
                  Nh = 10^3)

state_colony = c(H=1)

inter_colony = function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dH = b*(Nh - H)*H / Nh - d*H
    return(list(c(dH)))
  })
}

# model parameters
params_mites = c(p = 0.112, # intrinsic flower to bee infection rate
                 pm = .2,
                 q = 0.15, # intrinsic bee to flower infection rate
                 qm = .3,
                 d = 1/7, # death rate (lifespan = 1/d days)
                 dm = 1/5,
                 z = 1/5, # flower recovery/death
                 Nb = 100, # number bumble bees
                 Nh = 400, # num honey bees
                 Nf = 100,
                 Nm = 70) # num flowers         

# model start state. by default, a single honeybee
state_mites = c(B = 0, # B := num infected bumblebees
          H = 1, # H := num infected honeybees
          FL = 0,
          M = 0) # FL := num infected flowers

beeteam_mites = function(t, state, parameters) {
  # model difference equations
  with(as.list(c(state, parameters)), {
    dB = p*(Nb - B)*FL/Nf - d*B
    dH = .05 * (Nh - H)*FL/Nf - d*H +
        pm * (Nh - H) / Nh * M
    dM = qm * (H/Nh) * (Nm - M) - dm * M
    dFL = q*((Nf - FL)*(H + B))/Nf - z*FL
    return(list(c(dB, dH, dFL, dM)))
  })
}
