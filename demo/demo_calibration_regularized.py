import numpy as np

x_max = 10
y_max = 10
n_x = 10
n_y = 10
C = 1
R = n_x * n_y
T = 4 * 7
nb_groups = 2
nb_weeks = 10
groups = [[] for _ in range(nb_groups)]
which_group = [-1 for _ in range(T)]

for t in range(T):
    which_group[t] = t % nb_groups
    groups[which_group[t]].append(t)

is_red = [False for r in range(R)]
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
    elif (x >= mid_x and y < mid_y) or (x < mid_x and y >= mid_y):
        is_red[r] = False
        type_region[r] = 1
        
neighbors = [[] for _ in range(R)]
distance = np.ones((R,R), dtype=float)

for r in range(R):
    xi = (r + 1) % n_x
    yi = -1
    if xi == 0:
        yi = (r+1) // n_x
    else:
        yi = 1 + (r+1 - xi) // n_x
    
    if xi == 0:
        xi = n_x
    
    if yi - 1 >= 1:
        neighbors[r].append((yi - 2) * n_x + xi - 1)
    if yi + 1 <= n_y:
        neighbors[r].append(yi * n_x + xi - 1) 
    if xi - 1 >= 1:
        neighbors[r].append((yi - 1) * n_x + xi - 1 - 1)
    if xi + 1 <= n_x:
        neighbors[r].append((yi - 1) * n_x + xi + 1 - 1)   

nb_observations_total = nb_weeks
sample = np.zeros((T,R,C,nb_observations_total), dtype=int)
durations = np.ones(T,dtype=float)
theoretical_lambda = np.zeros((C,R,T), dtype=float)
for t in range(T):
    for r in range(R):
        x = r % 10
        y = r // 10
        cx = x + 0.5
        cy = y + 0.5
        mid_x = n_x // 2
        mid_y = n_y // 2
        for c in range(C):
            if (not is_red[r] and t % 2 == 0) or (is_red[r] and t % 2 != 0):
                theoretical_lambda[c,r,t] = 0.1
            elif (not is_red[r] and t % 2 != 0) or (is_red[r] and t % 2 == 0):
                theoretical_lambda[c,r,t] = 0.5
                
            for j in range(nb_observations_total):
                sample[t,r,c,j] = np.random.poisson(theoretical_lambda[c,r,t]*durations[t])             

weights = np.ones((nb_groups), dtype=float)
alpha = distance

cv_weights = [0, 0.0001, 0.0003, 0.0004, 0.0008, 0.001, 1]

def write_groups(path, groups, which_group, weights):
    time_groups_file = open(path,"w")
    time_groups_file.write(f"{len(groups)}\n")
    for t in range(len(which_group)):
        time_groups_file.write(f"{which_group[t]}\n")
    for g in range(len(groups)):
        time_groups_file.write(f"{weights[g]}\n")
    time_groups_file.close()
    
def write_durations(path,durations):
    times = open(path,"w")
    times.write(f"{len(durations)}\n")
    for t in range(len(durations)):
        times.write(f"{durations[t]}\n")
    times.close()

def write_alpha(path, alpha):
    alpha_regions_file = open(path, "w")
    for r in range(alpha.shape[0]):
        for s in range(alpha.shape[1]):
            alpha_regions_file.write(f"{alpha[r,s]} ")
        alpha_regions_file.write("\n")
    alpha_regions_file.close()

def write_info(path, T, R, C, nb_observations_total):
    info_file = open(path, "w")
    D = 7
    T_file = T // D
    nb_regressors = 0
    nb_holidays = 0
    info_file.write(f"{T_file} {D} {R} {C} {nb_regressors} {nb_holidays}\n")
    # Number of observations per day
    for _ in range(D):
        info_file.write(f"{nb_observations_total}\n")
    info_file.close()

def write_arrivals(path, sample):
    T,R,C,nb_observation_total = sample.shape
    arrivals_file = open(path, "w")
    D = 7
    for t in range(T):
            for r in range(R):
                for c in range(C):
                    for j in range(nb_observation_total):
                        #extract week day and period from t
                        d = t // (T // D) 
                        file_t = t % (T // D)
                        arrivals_file.write(f"{file_t} {d} {r} {c} {j} {sample[t,r,c,j]}\n")                    
    arrivals_file.write("END")
    arrivals_file.close()
    

def write_neighbors(path, R, neighbors, distance, type_region):
    neighbors_file = open(path, "w")
    for r in range(R):
        x = r % 10
        y = r // 10
        cx = x + 0.5
        cy = y + 0.5
        neighbors_file.write(f"{r} {cx} {cy} {type_region[r]} ")
        for s in neighbors[r]:
            neighbors_file.write(f"{s} {distance[r,s]} ")
        neighbors_file.write("\n")
    neighbors_file.write("END")
    neighbors_file.close()
    
def write_cv_weights(path, cv_weights):
    cv_weights_file = open(path,"w")
    for w in cv_weights:
        cv_weights_file.write(f"{w} ")
    cv_weights_file.write("\n")
    cv_weights_file.close()

base_dir = "calib_data_regularized"
app_params = {"EPS": 1e-5,
              "sigma": 0.5,
              "gap": 0.01,
              "max_iter": 30,
              "lower_lambda": 1e-6,
              "upper_lambda": 1.0,
              "beta_bar": 2.0,
              "cv_proportion": 0.2,
              "model_type": "no_reg",
              "method": "calibration",
              "algorithm": "feasible",
              "info_file": base_dir + "/info.txt",
              "arrivals_file": base_dir + "/arrivals.txt",
              "neighbors_file": base_dir + "/neighbors.txt",
              "alpha_regions_file": base_dir + "/alpha_regions.txt",
              "time_groups_file": base_dir + "/time_groups.txt",
              "cv_weights_file": base_dir + "/cv_weights.txt",
              "durations_file": base_dir + "/durations.txt",
              "output_file": base_dir + "/output_calib.txt"}



write_groups(app_params["time_groups_file"], groups, which_group, weights)
write_alpha(app_params["alpha_regions_file"], alpha)
write_info(app_params["info_file"], T, R, C, nb_observations_total)
write_arrivals(app_params["arrivals_file"], sample)
write_neighbors(app_params["neighbors_file"], R, neighbors, distance, type_region)
write_cv_weights(app_params["cv_weights_file"], cv_weights)
write_durations(app_params["durations_file"], durations)


from subprocess import PIPE, Popen


cfg_filename = "test_regularized.cfg"
cfg_file = open(cfg_filename, "w")
for k,v in app_params.items():
    cfg_file.write(f"{k} = {v}\n")
cfg_file.close()


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
# Running cross validation with selected weights
try:
    with Popen(
        ["../Model_Calibration/Cpp/laspated", "-f", cfg_filename, "--method", "cross_validation", "--output_file", "calib_data_regularized/output_cv.txt"],
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




