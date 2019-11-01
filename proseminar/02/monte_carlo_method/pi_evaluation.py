'''
Created on 22.10.2019

@author: marc_bussjaeger
'''
import numpy as np
import matplotlib.pyplot as plt

def get_trials_from_experiment_data(experiment_data):
    data = experiment_data.split(",")[0].split("trials = ")[1]
    return int(data)

def get_estimate_from_experiment_data(experiment_data):
    data = experiment_data.split("pi is ")[1].split("\nThe")[0]
    return float(data)

def get_seconds_from_experiment_data(experiment_data):
    data = experiment_data.split("took ")[1].split(" seconds")[0]
    return float(data)

def outfile_parser(outfile_path):
    # parse the txt into a list of runs that contains experiment dicts
    f = open(outfile_path, "r")
    data = f.read()
    f.close()
    data = data.replace("\r","").split("finish.")
    
    runs = [] # a run is consisting of 9 experiments
    run = []
    for i, experiment_data in enumerate(data[:-1]):
        if i%9==0:
            if len(run) is not 0:
                runs.append(run)
            run = []
        experiment = {
            "trials":get_trials_from_experiment_data(experiment_data),
            "estimate":get_estimate_from_experiment_data(experiment_data),
            "seconds":get_seconds_from_experiment_data(experiment_data)
            }
        run.append(experiment)
    runs.append(run)
    return runs

def accumulate_data(runs):
    # return the mean and stddev of the runs
    result = runs[0]
    for i, result_exp in enumerate(result):
        trials = result_exp["trials"]
        estimate_list = []
        second_list = []
        for run in runs:
            for exp in run:
                if exp['trials']==trials:
                    estimate_list.append(exp["estimate"])
                    second_list.append(exp["seconds"])
        result[i]["estimate"] = np.mean(estimate_list)
        result[i]["estimate_uncert"] = np.std(estimate_list)
        result[i]["seconds"] = np.mean(second_list)
        result[i]["seconds_uncert"] = np.std(second_list)
    return result

def sort_for_plotting(result):
    # sort the accumulated data for plotting
    trials = []
    estimates = []
    estimate_uncerts = []
    seconds = []
    seconds_uncerts = []
    
    for exp in result:
        trials.append(exp["trials"])
        estimates.append(exp["estimate"])
        estimate_uncerts.append(exp["estimate_uncert"])
        seconds.append(exp["seconds"])
        seconds_uncerts.append(exp["seconds_uncert"])
    return [trials, estimates, estimate_uncerts, seconds, seconds_uncerts]


# get the data
seq_runs = outfile_parser("./output_seq")
mpi_runs_16 = outfile_parser("./output_mpi")
mpi_runs_8 = outfile_parser("./output_mpi_8")
mpi_runs_4 = outfile_parser("./output_mpi_4")

seq_result = accumulate_data(seq_runs)
mpi_result_16 = accumulate_data(mpi_runs_16)
mpi_result_8 = accumulate_data(mpi_runs_8)
mpi_result_4 = accumulate_data(mpi_runs_4)

seq_plot_data = sort_for_plotting(seq_result)
mpi_plot_data_16 = sort_for_plotting(mpi_result_16)
mpi_plot_data_8 = sort_for_plotting(mpi_result_8)
mpi_plot_data_4 = sort_for_plotting(mpi_result_4)


# plot the estimate data
ax = plt.gca()
ax.errorbar(seq_plot_data[0], np.array(seq_plot_data[1])-np.pi, yerr=seq_plot_data[2],fmt=".", label="sequential")
ax.errorbar(mpi_plot_data_16[0], np.array(mpi_plot_data_16[1])-np.pi, yerr=mpi_plot_data_16[2],fmt=".", label="16 ranks")
ax.errorbar(mpi_plot_data_8[0], np.array(mpi_plot_data_8[1])-np.pi, yerr=mpi_plot_data_8[2],fmt=".", label="8 ranks")
ax.errorbar(mpi_plot_data_4[0], np.array(mpi_plot_data_4[1])-np.pi, yerr=mpi_plot_data_4[2],fmt=".", label="4 ranks")
ax.set_xlabel("trials")
ax.set_ylabel(r"$\pi_{\mathrm{estimate}}-\pi$")
plt.xscale('log')
legend = ax.legend()
plt.show()


# plot the seconds data
ax = plt.gca()
ax.errorbar(seq_plot_data[0], np.array(seq_plot_data[3]), yerr=seq_plot_data[4],fmt=".", label="sequential")
ax.errorbar(mpi_plot_data_16[0], np.array(mpi_plot_data_16[3]), yerr=mpi_plot_data_16[4],fmt=".", label="16 ranks")
ax.errorbar(mpi_plot_data_8[0], np.array(mpi_plot_data_8[3]), yerr=mpi_plot_data_8[4],fmt=".", label="8 ranks")
ax.errorbar(mpi_plot_data_4[0], np.array(mpi_plot_data_4[3]), yerr=mpi_plot_data_4[4],fmt=".", label="4 ranks")
ax.set_xlabel("trials")
ax.set_ylabel("time / s")
plt.xscale('log')
plt.yscale('log')
legend = ax.legend()
plt.show()



