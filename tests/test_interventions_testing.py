'''
Testing the effect of testing interventions in Covasim
'''

#%% Imports and settings
import matplotlib
matplotlib.use('Agg')
import sciris as sc
import covasim as cv

do_plot   = 1
do_show   = 0
do_save   = 1
debug     = 1
keep_sims = 0
fig_paths = [f'results/testing_scen_{i}.png' for i in range(4)]


def test_interventions(do_plot=False, do_show=True, do_save=False, fig_path=None):
    sc.heading('Test of testing interventions')


    sc.heading('Setting up...')

    sc.tic()

    n_runs = 3
    verbose = 1
    base_pars = {
      'pop_size': 1000,
      'use_layers': True,
      }

    base_sim = cv.Sim(base_pars) # create sim object
    n_people = base_sim['pop_size']
    npts = base_sim.npts

    # Define overall testing assumptions
    # As the most optimistic case, we assume countries could get to South Korea's testing levels. S Korea has tested
    # an average of 10000 people/day over March, or 270,000 in total. This is ~200 people per million every day (0.02%).
    max_optimistic_testing = 0.0002
    optimistic_daily_tests = [max_optimistic_testing*n_people]*npts # Very best-case scenario for asymptomatic testing

    # Define the scenarios
    scenarios = {
        'baseline': {
          'name':'Status quo, no testing',
          'pars': {
              'interventions': None,
              }
          },
        'test_skorea': {
          'name':'Assuming South Korea testing levels of 0.02% daily (untargeted); isolate positives',
          'pars': {
              'interventions': cv.test_num(daily_tests=optimistic_daily_tests)
              }
          },
        'tracing': {
          'name':'Assuming South Korea testing levels of 0.02% daily (with contact tracing); isolate positives',
          'pars': {
              'interventions': [cv.test_num(daily_tests=optimistic_daily_tests),
                                cv.dynamic_pars({'quar_trans_factor':{'days':20, 'vals':0.1}})] # This means that people who've been in contact with known positives isolate with 90% effectiveness
              }
          },
        'floating': {
            'name': 'Test with constant probability based on symptoms',
            'pars': {
                'interventions': cv.test_prob(symptomatic_prob=max_optimistic_testing, asymptomatic_prob=0.0)
                }
        },
        'historical': {
            'name': 'Test a known number of positive cases',
            'pars': {
                'interventions': cv.test_historical(n_tests=[100]*npts, n_positive = [1]*npts)
            }
        },
        'sequence': {
            'name': 'Historical switching to probability',
            'pars': {
                'interventions': cv.sequence(days=[10, 51], interventions=[
                    cv.test_historical(n_tests=[100] * npts, n_positive=[1] * npts),
                    cv.test_prob(symptomatic_prob=0.2, asymptomatic_prob=0.002),
                ])
            }
        },

    }

    metapars = {'n_runs': n_runs}

    scens = cv.Scenarios(sim=base_sim, metapars=metapars, scenarios=scenarios)
    scens.run(verbose=verbose, debug=debug)

    if do_plot:
        scens.plot(do_save=do_save, do_show=do_show, fig_path=fig_path)

    return scens


def test_turnaround(do_plot=False, do_show=True, do_save=False, fig_path=None):
    sc.heading('Test impact of reducing delay time for getting test results')

    sc.heading('Setting up...')

    sc.tic()

    n_runs = 3
    verbose = 1
    base_pars = {
      'pop_size': 5000,
      'use_layers': True,
      }

    base_sim = cv.Sim(base_pars) # create sim object
    n_people = base_sim['pop_size']
    npts = base_sim.npts

    # Define overall testing assumptions
    testing_prop = 0.1 # Assumes we could test 10% of the population daily (!!)
    daily_tests = [testing_prop*n_people]*npts # Number of daily tests

    # Define the scenarios
    scenarios = {
        f'{d}dayturnaround': {
            'name':f'Symptomatic testing with {d} days to get results',
            'pars': {
                'interventions': cv.test_num(daily_tests=daily_tests, test_delay=d)
            }
        } for d in range(1, 7+1)
    }

    metapars = {'n_runs': n_runs}

    scens = cv.Scenarios(sim=base_sim, metapars=metapars, scenarios=scenarios)
    scens.run(verbose=verbose, debug=debug)

    if do_plot:
        scens.plot(do_save=do_save, do_show=do_show, fig_path=fig_path)

    return scens


def test_sanitycheck_tracedelay(do_plot=False, do_show=True, do_save=False, fig_path=None):
    sc.heading('Sanity check for contact tracing: removing the intervention vs setting effects to zero')

    sc.tic()

    n_runs = 3
    verbose = 1
    base_pars = {
      'pop_size': 5000,
      'use_layers': True,
      }

    base_sim = cv.Sim(base_pars) # create sim object
    n_people = base_sim['pop_size']
    npts = base_sim.npts

    # Define overall testing assumptions
    testing_prop = 0.01 # Assumes we could test 1% of the population daily
    daily_tests = [testing_prop*n_people]*npts # Number of daily tests

    # Define the scenarios
    scenarios = {
        'none': {
            'name': 'No contact tracing',
            'pars': {
                'interventions': [
                    cv.test_num(daily_tests=daily_tests, trace_test=1.0)]
            }
        },
        'ineffectivetracing': {
            'name': "Contact tracing intervention in place, but it doesn't work",
            'pars': {
                'quar_trans_factor': 1.,
#                'quar_acq_factor': 0,
                'interventions': [
                    cv.test_num(daily_tests=daily_tests, trace_test=1.0),
                    cv.contact_tracing(trace_probs = {'h': 0, 's': 0, 'w': 0, 'c': 0},
                                       trace_time  = {'h': 0, 's': 7,   'w': 7,   'c': 10})]
            }
        },
    }

    metapars = {'n_runs': n_runs}

    scens = cv.Scenarios(sim=base_sim, metapars=metapars, scenarios=scenarios)
    scens.run(verbose=verbose, debug=debug)

    if do_plot:
        scens.plot(do_save=do_save, do_show=do_show, fig_path=fig_path)

    return scens


def test_tracedelay(do_plot=False, do_show=True, do_save=False, fig_path=None):
    sc.heading('Test impact of reducing delay time for finding contacts of positives')

    sc.tic()

    n_runs = 3
    verbose = 1
    base_pars = {
      'pop_size': 20000,
      'use_layers': True,
      }

    base_sim = cv.Sim(base_pars) # create sim object
    n_people = base_sim['pop_size']
    npts = base_sim.npts

    # Define overall testing assumptions
    testing_prop = 0.01 # Assumes we could test 1% of the population daily
    daily_tests = [testing_prop*n_people]*npts # Number of daily tests

    # Define the scenarios
    scenarios = {
        'ineffectivetracing': {
            'name': "Very poor contact tracing, 10% of contacts isolate for 7 days",
            'pars': {
                'quar_trans_factor': 1.,
                'quar_acq_factor':  0.,
                'quarantine_period': 7,
                'interventions': [
                    cv.test_num(daily_tests=daily_tests),
                    cv.contact_tracing(start_day=21,
                                       trace_probs = {'h': 1, 's': 0.8, 'w': 0.5, 'c': 0.1},
                                       trace_time  = {'h': 0, 's': 7,   'w': 7,   'c': 10})]
            }
        },
        'lowtrace': {
            'name': 'Poor contact tracing, 50% of contacts self-isolate for 7 days',
            'pars': {
                'quar_trans_factor': 0.5,
                'quar_acq_factor':  0.5,
                'quarantine_period': 7,
                'interventions': [
                    cv.test_num(daily_tests=daily_tests),
                    cv.contact_tracing(start_day=21,
                                       trace_probs = {'h': 1, 's': 0.8, 'w': 0.5, 'c': 0.1},
                                       trace_time  = {'h': 0, 's': 7,   'w': 7,   'c': 10})]
            }
        },
        'modtrace': {
            'name': 'Moderate contact tracing, 75% of contacts self-isolate for 10 days',
            'pars': {
                'quar_trans_factor': 0.25,
                'quar_acq_factor': 0.75,
                'quarantine_period': 10,
                'interventions': [
                    cv.test_num(daily_tests=daily_tests),
                    cv.contact_tracing(start_day=21,
                                       trace_probs = {'h': 1, 's': 0.8, 'w': 0.5, 'c': 0.1},
                                       trace_time  = {'h': 0,  's': 3,  'w': 3,   'c': 8})]
            }
        },
        'hightrace': {
            'name': 'Fast contact tracing, 90% of contacts self-isolate for 14 days',
            'pars': {
                'quar_trans_factor': 0.1,
                'quar_acq_factor': 0.9,
                'quarantine_period': 14,
                'interventions': [
                    cv.test_num(daily_tests=daily_tests),
                    cv.contact_tracing(start_day=21,
                                       trace_probs = {'h': 1, 's': 0.8, 'w': 0.8, 'c': 0.2},
                                       trace_time  = {'h': 0, 's': 1,   'w': 1,   'c': 5})]
            }
        },
        'crazy': {
            'name': 'Same-day contact tracing, 100% of contacts self-isolate for 21 days',
            'pars': {
                'quar_trans_factor': 0,
                'quar_acq_factor': 1,
                'quarantine_period': 21,
                'interventions': [
                    cv.test_num(daily_tests=daily_tests),
                    cv.contact_tracing(start_day=21,
                                       trace_probs = {'h': 1, 's': 1, 'w': 1, 'c': 1},
                                       trace_time  = {'h': 0, 's': 0, 'w': 0, 'c': 0})]
            }
        },
    }

    metapars = {'n_runs': n_runs}

    scens = cv.Scenarios(sim=base_sim, metapars=metapars, scenarios=scenarios)
    scens.run(verbose=verbose, debug=debug)

    if do_plot:
        scens.plot(do_save=do_save, do_show=do_show, fig_path=fig_path)

    return scens




#%% Run as a script
if __name__ == '__main__':
    sc.tic()

#    scens1 = test_interventions(do_plot=do_plot, do_save=do_save, do_show=do_show, fig_path=fig_paths[0])
#    scens2 = test_turnaround(do_plot=do_plot, do_save=do_save, do_show=do_show, fig_path=fig_paths[1])
    scens3 = test_sanitycheck_tracedelay(do_plot=do_plot, do_save=do_save, do_show=do_show, fig_path=fig_paths[2])
#    scens4 = test_tracedelay(do_plot=do_plot, do_save=do_save, do_show=do_show, fig_path=fig_paths[3])

    sc.toc()


print('Done.')
