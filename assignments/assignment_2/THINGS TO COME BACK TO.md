THINGS TO COME BACK TO 
----------------------
(or notes on things I need elsewhere)

- 2a: choice of likelihood model (Poisson) and relation to the progression of bin count

- 4.a.ii: for investigation, maybe just plot some of the MCMC chains above each other on one plot? kind of like I did with the flux ramp for Ali (LATER, IF TIME). 

- come back to corner plot (I HATE THIS) and just use code from assignment 1 to generate the new plot

- come back to the analysis of the first corner plot, talk about correlation (I think they are all highly correlated )


- go into the posterior predictive check function and take out the "Poisson Replicates" in title


Between the previous plot and this one, what has changed? 
- the simulated profiles follow a shallower descent from the maximum, and follow a shape that resembles the real data 
- the real data points are actually all contained within the spreads of the simulated data (however the peak $N_j$ at 5 arcminutes has an error bar which extends beyond the simulated spread) but overall the new profile fits quite nicely! 
- even the smaller-scale distributions which I was worried about with my binning, at very small R, are actually captured (the increase from the first to second bin, so peak density occurs technically outside the very center)
- the variance of the simulated data also is roughly adequate to hold the error bars of the observed data as well (whereas before I felt that the simulated data was not disperse enough) 














import emcee 
# parameter setting 
n_dimensions = 3 # for the three theta values
n_walkers_init = 20 # number of chains
n_steps_init = 5000 # steps per chain 
n_burned_init = 500 # burn in steps 
total_steps = n_steps_init * n_walkers_init
print('Initial Parameters for MCMC: \n----------------------------')
print(f'Walkers: {n_walkers_init}')
print(f'Steps per chain: {n_steps_init}')
print(f'Burned steps per chain: {n_burned_init}')
print(f'Total Steps: {total_steps}')

# choosing starting position

sigma_30 = Counts_30 / Areas_30
sigma_30_bg = np.mean(sigma_30[-5:])
sigma_0_guess = np.max(sigma_30 - sigma_30_bg)
print(sigma_0_guess)

# do first position guesses: 
first_pos_guesses = np.array([sigma_0_guess, 0.999999999999999, 9.999999999999999]) 

rng = np.random.default_rng(33)

def get_randomized_inf_check(starting_theta, n_walkers, n_dims):
    '''Just making this in case I need to re-generate the random sample with different emcee params
    COME BACK AND ANNOTATE THIS FUNCTION
    '''
    # randomize position
    position = starting_theta + 0.1*rng.standard_normal(size=(n_walkers, n_dims)) # variance of 10^-2

    # this is a little while bit that makes sure things are within the bounds of the prior
    # and resamples if it has an infinite output 
    out_of_bounds = ~np.isfinite([joint_prior_prob(pos, 10) for pos in position])
    while np.any(out_of_bounds): # this just makes sure 
        position[out_of_bounds] = starting_theta + 0.1 * rng.standard_normal(size = (out_of_bounds.sum(), n_dims))
        out_of_bounds = ~np.isfinite([joint_prior_prob(pos, 10) for pos in position])
    
    return position

# try out randomizer

pos1 = get_randomized_inf_check(first_pos_guesses, n_walkers_init, n_dimensions)

# this block I am using with the help of my friend and classmate, Ana Mare

sampler = emcee.EnsembleSampler(n_walkers_init, n_dimensions, log_post)
print("Sampler Running")
sampler.run_mcmc(pos1, n_steps_init, progress = True)
samples = sampler.get_chain(discard = 0, flat = False)

# first run: 00:01:02


titles = ["Central surface density", "Core radius", "Tidal radius"]
variables = ["$\Sigma_0$", "$R_c$", "$R_t$"]
dimensions = ["stars/arcmin$^2$", "arcmin", "arcmin"]
for i in range(0, 3):
    plt.plot(samples[:, :, i], alpha = 1, linewidth = 0.5)
    plt.title(f"{titles[i]} MCMC")
    plt.xlabel("step number")
    plt.ylabel(f"{variables[i]} [{dimensions[i]}]")
    plt.axvspan(0, n_burned_init, alpha=0.5, color='orange', label = "burn area")
    plt.legend()
    plt.show()

# re-running
burn_in = 300
sampler = emcee.EnsembleSampler(n_walkers_init, n_dimensions, log_post)
print("Sampler Running")
sampler.run_mcmc(pos1, n_steps_init, progress = True)
samples_burned = sampler.get_chain(discard = burn_in, flat = False) 

titles = ["Central surface density", "Core radius", "Tidal radius"]
variables = ["$\Sigma_0$", "$R_c$", "$R_t$"]
dimensions = ["stars/arcmin$^2$", "arcmin", "arcmin"]
print('MCMC Chains, Burn-out:')
for i in range(0, 3):
    plt.plot(samples_burned[:, :, i], alpha = 1, linewidth = 0.5)
    plt.title(f"{titles[i]} MCMC")
    plt.xlabel("step number")
    plt.ylabel(f"{variables[i]} [{dimensions[i]}]")
    plt.ylim(np.mean(samples_burned[:,:,i]) - 15*np.std(samples_burned[:,:,i]), 
             np.mean(samples_burned[:,:,i]) + 15*np.std(samples_burned[:,:,i])) # just to see a bit better
    plt.show()


# Computing integrated autocorrelation time for each parameter
n_steps_burnout = n_steps_init - burn_in
taus = []
N_effs = []
taus_bigsample = []
for i in range(3): 
    tau = emcee.autocorr.integrated_time(samples[:,:,i], c = 5, tol = 50, quiet = False)
    N_eff = n_steps_burnout/tau
    tau_bigsample = tau*n_steps_burnout
    # these will then be an array of times in order Sigma_0, R_c, R_t
    taus.append(tau) 
    N_effs.append(N_eff)
    taus_bigsample.append(tau_bigsample)

print(f'Autocorrelation Times:\n---------------------- \nSigma_0: {taus[0]} \nR_c: {taus[1]} \nR_t: {taus[2]}\n')
print(f'Effective Sample Sizes:\n----------------------- \nSigma_0: {N_effs[0]} \nR_c: {N_effs[1]} \nR_t: {N_effs[2]}\n')
print(f'Times to achieve N_eff:\n----------------------- \nSigma_0: {taus_bigsample[0]} \nR_c: {taus_bigsample[1]} \nR_t: {taus_bigsample[2]}')
