require(deSolve)
require(tidyverse)
require(ggthemes)
require(wesanderson)

## core model: bumblebees and honeybees share transmission route via flowers

# model parameters
params_base = c(p = 0.17, # intrinsic flower to bee infection rate (mean)
                q = 0.3, # intrinsic bee to flower infection rate (mean)
                d = 1/7, # death rate (lifespan = 1/d days)
                z = 1/10, # flower recovery/death
                Nb = 100, # number bumble bees
                Nh = 100, # num honey bees
                Nf = 100) # num flowers

# model start state. by default, a single honeybee
## state = c(B = 0, # B := num infected bumblebees
##           H = 1, # H := num infected honeybees
##           FL = 0) # FL := num infected flowers

## beeteam = function(t, state, parameters) {
##   # model difference equations
##   with(as.list(c(state, parameters)), {
##     dB = p*(Nb - B)*FL/Nf - d*B
##     dH = p*(Nh - H)*FL/Nf - d*H
##     dFL = q*((Nf - FL)*(H + B))/Nf - z*FL
##     return(list(c(dB, dH, dFL)))
##   })
## }

## quick_run = function() {
##   # run the ode with "default" values
##   # wrapper to ease notation and ensure above global values are used
##   ode(y=state, times=times, func=beeteam, parms=params)
## }

run_model = function(start = set_state(1, 1),
                     max_t = 100,
                     parms = set_params(1)) {
    ode(y=start, times=0:max_t, func=beeteam, parms=parms) %>%
        as_tibble
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

### default colors for each population
color_pal = wes_palette('Darjeeling2', 3)
names(color_pal) = c('B', 'H', 'F')
color_pal2 = wes_palette('Darjeeling2', 3)
names(color_pal2) = c('Nb', 'Nh', 'Nf')

## color_pal = function(n=fl_types) {
##     ret = plasma(n+2)
##     names(ret) = c('B', 'H', str_c('F', 1:n))
##     ret
## }

## color_pal2 = function(n=fl_types) {
##     ret = viridis(n+2)
##     names(ret) = c('Nb', 'Nh', str_c('Nf', 1:n))
##     ret
## }

plot_time_series = function(out = run_model(),
                            pop = colnames(out[,-1]),
                            colors = color_pal) {
  # plot population infection levels as function of time
  # input: 
    # ode: a matrix of class deSolve
    # pop: columns of ode to select, representing infection levels
    # colors: a color pallette. Must have names corresponding to names in 'pop'
  # output:
    # a ggplot object
    out %>%
        gather(type, inf, starts_with('F')) %>%
        ggplot(aes(x=time, y=inf)) +
        stat_summary(fun.y = sum, col='red', geom='smooth', linetype='dotted', alpha=.7) +
        stat_summary(fun.data = median_hilow, fun.args = list(conf.int = .5), fill=grey(.5), geom='ribbon', alpha=.5) +
        stat_summary(fun.y = median, col='red', geom='smooth', size=1.2) +
        geom_line(mapping=aes(x=time, y=B), col='blue', size=1.2) +
        geom_line(mapping=aes(x=time, y=H), col='green', size=1.2) +
        labs(y="Infections", x="Time (Days)") +
        theme_few()
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

### heterogenus flower type model

set_params = function(fl_types = 1,
                       pd = pmax(rnorm(fl_types, params['p'], ifelse(fl_types > 1, .1, 0)), 0),
                       qd = pmax(rnorm(fl_types, params['q'], ifelse(fl_types > 1, .15, 0)), 0),
                       zd = pmax(rnorm(fl_types, params['z'], ifelse(fl_types > 1, .1, 0)), 0)) {
    with(as.list(c(params_base)), {
        ## pd = if (is.null(pd)) rep(p, fl_types) else pd
        ## qd = if (is.null(qd)) rep(q, fl_types) else qd
        ## zd = if (is.null(zd)) rep(z, fl_types) else zd
        ret = c(pd, qd, d, zd, Nb, Nh, rep(Nf %/% fl_types, fl_types))
        names(ret) = c(str_c('p', 1:fl_types), str_c('q', 1:fl_types), 'd', str_c('z', 1:fl_types), 'Nb', 'Nh', str_c('Nf', 1:fl_types))
        return(ret)
    })
}

set_state = function(fl_types=1, H=1, B=0) {
    with(as.list(c(params_base)), {
        ret = c(B, H, rep(0, fl_types))
        names(ret) = c('B', 'H', str_c('F', 1:fl_types))
        return(ret)
    })
}

## params_flowers = make_params(pmax(rnorm(fl_types, params['p'], .1), 0),
##                              pmax(rnorm(fl_types, params['q'], .15), 0),
##                              pmax(rnorm(fl_types, params['z'], .1), 0))
## state_flowers = make_state(10)

get_fl_params = function(i, parameters) {
    c(parameters[str_c('p', i)],
      parameters[str_c('q', i)],
      parameters[str_c('z', i)],
      parameters[str_c('Nf', i)])
}

de_flowers = function(i, state, parameters) {
    with(as.list(c(state, parameters)), {
        fp = get_fl_params(i, parameters)
        Nf = sum(parameters[str_detect(names(parameters), 'Nf')])
        return(fp[2] * ((fp[4] - state[str_c('F', i)]) * (H + B)) / Nf - fp[3] * state[str_c('F', i)])
    })
}

de_bumble = function(i, state, parameters) {
    with(as.list(c(state, parameters)), {
        fp = get_fl_params(i, parameters)
        Nf = sum(parameters[str_detect(names(parameters), 'Nf')])
        return(fp[1] * (Nb - B) * state[str_c('F', i)] / Nf)
    })
}

de_honey = function(i, state, parameters) {
    with(as.list(c(state, parameters)), {
        fp = get_fl_params(i, parameters)
        Nf = sum(parameters[str_detect(names(parameters), 'Nf')])
        return(fp[1] * (Nh - H) * state[str_c('F', i)] / Nf)
    })
}

beeteam = function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
        fl_types = length(state) - 2
        dB = sum(map_dbl(1:fl_types, de_bumble, state, parameters)) - d * B
        dH = sum(map_dbl(1:fl_types, de_honey, state, parameters)) - d * H
        dFs = 1:fl_types %>% map_dbl(de_flowers, state, parameters)
        return(list(c(dB, dH, dFs)))
  })
}

### single flower model, with mites!

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
                 Nm = 70) 

state_mites = c(B = 0, 
          H = 1,
          FL = 0,
          M = 0)

beeteam_mites = function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dB = p*(Nb - B)*FL/Nf - d*B
    dH = .05 * (Nh - H)*FL/Nf - d*H +
        pm * (Nh - H) / Nh * M
    dM = qm * (H/Nh) * (Nm - M) - dm * M
    dFL = q*((Nf - FL)*(H + B))/Nf - z*FL
    return(list(c(dB, dH, dFL, dM)))
  })
}

## ## additional simple SIS model for intercolony dynamics between bees

## params_colony = c(b = .18,
##                   d = 1/7,
##                   Nh = 10^3)

## state_colony = c(H=1)

## inter_colony = function(t, state, parameters) {
##   with(as.list(c(state, parameters)), {
##     dH = b*(Nh - H)*H / Nh - d*H
##     return(list(c(dH)))
##   })
## }
