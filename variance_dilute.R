source('beeteam.R')

### single flower case: use only mean values, deterministic
plot_time_series()


### multi-flower case: use only mean values, should be same as above
out = run_model(start=set_state(10),
                parms=set_params(10, rep(params_base['p'], 10), rep(params_base['q'], 10), rep(params_base['z'], 10)))
plot_time_series(out)

### multi-flower case: sample from rnorm using base parameters
run_model(start=set_state(10), parms=set_params(10)) %>%
    plot_time_series
