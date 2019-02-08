require(deSolve)
require(tidyverse)
require(ggthemes)
require(wesanderson)

## core model: bumblebees and honeybees share transmission route via flowers

# model parameters
params = c(k = 100, # number of flowers visited per day
           p = .05, # prob. infection picked up from flower
           q = .01, # prob infection deposited on flower
           d = 1/7, # death rate (lifespan = 1/d days)
           z = 1/5, # flower recovery/death
           Nb = 10^2, # number bumble bees
           Nh = 10^2, # num honey bees
           Nf = 10^3) # num flowers 

# model start state. by default, a single honeybee
state = c(B = 0, # B := num infected bumblebees
          H = 1, # H := num infected honeybees
          FL = 0) # FL := num infected flowers

# timesteps for which model output should be recorded
times = seq(0, 100, by=1)

beeteam = function(t, state, parameters) {
  # model difference equations
  with(as.list(c(state, parameters)), {
    dB = k*p*(Nb - B)*FL/Nf - d*B
    dH = k*p*(Nh - H)*FL/Nf - d*H
    dF = k*q*((Nf - FL)*(H + B))/Nf - z*FL
    return(list(c(dB, dH, dF)))
  })
}

run_model = function() {
  # wrapper to ease notation and ensure above global values are used
  ode(y=state, times=times, func=beeteam, parms=params)
}

## additional functions for aiding in analysis

h_fixed_point = function(nf = params[8], nh = params[7], nb=params[6]) {
  # calculate the fixed point for honey bees
  with(as.list(params), {
    num = (nh + nb) * k^2 * p * q - d * nf * z
    denom = (1 + nb / nh) * k * (d + k * p) * q
    return(num / denom)
  })
}

b_fixed_point = function(nh = params[7], nb=params[6], nf=params[8]) {
  # calculate the fixed point for bumble bees
  with(as.list(params), {
    return(nb / nh * h_fixed_point(nf=nf, nh=nh, nb=nb))
  })
}

fl_fixed_point = function(nh = params[7], nb=params[6], nf=params[8]) {
  # calculate the fixed point for flowers
  with(as.list(params), {
    num = nf * ((nh + nb) * k^2 * p * q - d * nf * z)
    denom = k * p * ((nh + nb) * k * q + nf * z)
    return(num/denom)
  })
}

h_crit_point = function(n) {
  # finds the critical value of H, past which the infection will survive
  with(as.list(params), {
    num = d * z * Nf
    s = sum(bb(0:n))
    denom = s * k^2 * p * q
    return(num / denom)
  })
}

bb = function(i) {
  # compute the ith term in the series for the reproductive number
  with(as.list(params), {
    num = k^2 * p * q * Nb
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

plot_time_series = function(ode = NULL, 
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
