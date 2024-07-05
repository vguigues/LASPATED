import pandas as pd
import geopandas as gpd

import matplotlib.pyplot as plt
import numpy as np

import laspated as spated
from shapely.geometry import Polygon, MultiPolygon, Point
import sys

x_max = 10
y_max = 10
n_x = 10
n_y = 10
R = n_x * n_y
C = 1
T = 4
D = 7

nb_weeks = 52
nb_obs = nb_weeks * 7
durations = np.ones(T, dtype=float)
nb_land_types = 2
nb_regressors = 1 + nb_land_types
theoretical_beta = np.zeros((C,D,T,nb_regressors), dtype=float)
regressors = np.zeros((nb_regressors, R), dtype=float)

for d in range(D):
    theoretical_beta[0, d, 1, 0] = 0.05
    theoretical_beta[0, d, 3, 0] = 0.05

    theoretical_beta[0, d, 0, 1] = 6
    theoretical_beta[0, d, 1, 1] = 18
    theoretical_beta[0, d, 2, 1] = 6
    theoretical_beta[0, d, 3, 1] = 18

    theoretical_beta[0, d, 0, 2] = 3
    theoretical_beta[0, d, 1, 2] = 6
    theoretical_beta[0, d, 2, 2] = 3
    theoretical_beta[0, d, 3, 2] = 6

is_red = [False for _ in range(R)]
type_region = [-1 for r in range(R)]
for r in range(R):
    x = r % 10
    y = r // 10
    cx = x + 0.5
    cy = y + 0.5
    mid_x = n_x // 2
    mid_y = n_y // 2
    if (x < mid_x and y < mid_y) or (x >= mid_x and y >= mid_y):
        is_red[r] = True
        type_region[r] = 0
        regressors[0, r] = 1
        regressors[1, r] = 0.5
        regressors[2, r] = 0.25
    elif (x >= mid_x and y < mid_y) or (x < mid_x and y >= mid_y):
        is_red[r] = False
        type_region[r] = 1
        regressors[0, r] = 1
        regressors[1, r] = 0.25
        regressors[2, r] = 0.5

sample = np.zeros((C,D,T,R, nb_weeks), dtype=int)
nb_observations = np.zeros((C,D,T,R), dtype=int)
for index in range(nb_obs):
    day = index % 7
    for c in range(C):
        for t in range(T):
            for r in range(R):
                rate = 0.0
                for j in range(nb_regressors):
                    rate += theoretical_beta[c,day,t,j] * regressors[j,r]
                sample[c,day,t,r, nb_observations[c,day,t,r]] = np.random.poisson(rate)
                nb_observations[c,day,t,r] += 1

def write_info(path, T, R, C, D, nb_regressors, nb_obs):
    info_file = open(path, "w")
    nb_holidays = 0
    info_file.write(f"{T} {D} {R} {C} {nb_regressors} {nb_holidays}\n")
    for _ in range(D):
        info_file.write(f"{nb_obs}\n")
    info_file.close()

def write_arrivals(path, sample):
    C,D,T,R,nb_obs = sample.shape
    arrivals_file = open(path, "w")
    for c in range(C):
        for d in range(D):
            for t in range(T):
                for r in range(R):
                    for j in range(nb_obs):
                        arrivals_file.write(f"{t} {d} {r} {c} {j} {sample[c,d,t,r,j]}\n")                               
    arrivals_file.write("END")
    arrivals_file.close()

def write_neighbors(path,  R, regressors, type_region):
    neighbors_file = open(path, "w")
    for r in range(R):
        x = r % 10
        y = r // 10
        cx = x + 0.5
        cy = y + 0.5
        neighbors_file.write(f"{r} {cx} {cy} {type_region[r]} ")
        for j in range(regressors.shape[0]):
            neighbors_file.write(f"{regressors[j,r]} ")
        neighbors_file.write("\n")
    neighbors_file.write("END")
    neighbors_file.close()

def write_durations(path,durations):
    times = open(path,"w")
    times.write(f"{len(durations)}\n")
    for t in range(len(durations)):
        times.write(f"{durations[t]}\n")
    times.close()
    
base_dir = "calib_data_covariates"    
app_params = {"EPS": 1e-5,
              "sigma": 0.5,
              "gap": 0.01,
              "max_iter": 30,
              "lower_lambda": 1e-6,
              "upper_lambda": 1.0,
              "beta_bar": 2.0,
              "cv_proportion": 0.2,
              "model_type": "reg",
              "method": "calibration",
              "algorithm": "feasible",
              "info_file": base_dir + "/info.txt",
              "arrivals_file": base_dir + "/arrivals.txt",
              "neighbors_file": base_dir + "/neighbors.txt",
              "durations_file": base_dir + "/durations.txt",
              "output_file": base_dir + "/output_calib.txt"}


write_info(app_params["info_file"], T, R, C, D, nb_regressors, nb_weeks)
write_arrivals(app_params["arrivals_file"], sample)
write_neighbors(app_params["neighbors_file"], R, regressors, type_region)
write_durations(app_params["durations_file"], durations)

from subprocess import PIPE, Popen


cfg_filename = "test_covariates.cfg"
cfg_file = open(cfg_filename, "w")
for k,v in app_params.items():
    cfg_file.write(f"{k} = {v}\n")
cfg_file.close()


print("Params:")
for k,v in app_params.items():
    print(f"{k}\t\t\t\t{v}")

# Running calibration with unitary weights
try:
    with Popen(
        ["../Model_Calibration/Cpp/laspated", "-f", cfg_filename],
        bufsize=1,
        stdout=PIPE,
        stderr=PIPE,
        text=True,
    ) as p:
        for line in p.stdout:
            print(line, end="")

        for line in p.stderr:
            print(line, end="")
except:
    print(f"Error with C++ experiments")
    exit(1)
    
print("Finished Calibration")