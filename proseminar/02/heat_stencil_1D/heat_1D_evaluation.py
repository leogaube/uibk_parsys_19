'''
Created on 22.10.2019

@author: marc_bussjaeger
'''
def get_room_size_from_experiment_data(experiment_data):
    data = experiment_data.split("size N=")[1].split(" for T")[0]
    return int(data)

def get_time_steps_from_experiment_data(experiment_data):
    data = experiment_data.split(" for T=")[1].split(" timesteps")[0]
    return int(data)

def get_data_from_experiment_data(experiment_data):
    data = experiment_data.split("Initial")[1].split("Final")[0]
    if "Verification" in data:
        data = data.split("Verification")[0]
    return data

def outfile_parser(outfile_path):
    # parse the txt into a list of runs that contains experiment dicts
    f = open(outfile_path, "r")
    data = f.read()
    f.close()
    data = data.replace("\r","").split("Computing")
    
    runs = [] # a run is consisting of 9 experiments
    run = []
    for i, experiment_data in enumerate(data[1:-1]):
        if i%1==0:
            if len(run) is not 0:
                runs.append(run)
            run = []
        experiment = {
            "room_size":get_room_size_from_experiment_data(experiment_data),
            "time_steps":get_time_steps_from_experiment_data(experiment_data),
            "data":get_data_from_experiment_data(experiment_data)
            }
        run.append(experiment)
    runs.append(run)
    return runs

# get the data
seq_runs = outfile_parser("./output_seq")
mpi_runs_16 = outfile_parser("./output_mpi_16")
mpi_runs_8 = outfile_parser("./output_mpi_8")
mpi_runs_4 = outfile_parser("./output_mpi_4")

# check if runs are the same
for i, seq_run in enumerate(seq_runs):
    print("seq and mpi_16 check:")
    for j in range(len(seq_run[0]["data"])):
        if seq_run[0]["data"][j] != mpi_runs_16[i][0]["data"][j]:
            print("They missmatch in character", j)
            print(seq_run[0]["data"][max([0, j-50]):min([len(seq_run[0]["data"]), j+50])])
            print(mpi_runs_16[i][0]["data"][max([0, j-50]):min([len(mpi_runs_16[i][0]["data"]), j+50])])
    print("lengths:", len(seq_run[0]["data"]), len(mpi_runs_16[i][0]["data"]))

    print("seq and mpi_8 check:")
    for j in range(len(seq_run[0]["data"])):
        if seq_run[0]["data"][j] != mpi_runs_8[i][0]["data"][j]:
            print("They missmatch in character", j)
            print(seq_run[0]["data"][max([0, j-50]):min([len(seq_run[0]["data"]), j+50])])
            print(mpi_runs_8[i][0]["data"][max([0, j-50]):min([len(mpi_runs_8[i][0]["data"]), j+50])])
    print("lengths:", len(seq_run[0]["data"]), len(mpi_runs_8[i][0]["data"]))

    print("seq and mpi_4 check:")
    for j in range(len(seq_run[0]["data"])):
        if seq_run[0]["data"][j] != mpi_runs_4[i][0]["data"][j]:
            print("They missmatch in character", j)
            print(seq_run[0]["data"][max([0, j-50]):min([len(seq_run[0]["data"]), j+50])])
            print(mpi_runs_4[i][0]["data"][max([0, j-50]):min([len(mpi_runs_4[i][0]["data"]), j+50])])
    print("lengths:", len(seq_run[0]["data"]), len(mpi_runs_4[i][0]["data"]))

print("Analysis done. (No mismatches if no other comments.)")