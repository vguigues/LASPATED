import numpy as np

R = 76
S = 30
D = 7
T = 48
C = 3

nb_obs = [104, 104, 104, 104, 105, 105, 105]
max_obs = np.max(nb_obs)


def read_sample(path="Rect10x10/missing_calls.dat"):
    arq = open(path, "r")
    sample_missing_calls = np.zeros((C, D, T, max_obs))
    for line in arq.readlines():
        if line == "END":
            break
        tokens = line.split()
        t = int(tokens[0])
        d = int(float(tokens[1]))
        c = int(tokens[3])
        n = int(tokens[4])
        val = int(tokens[5])
        sample_missing_calls[c, d, t, n] = val
    arq.close()

    return sample_missing_calls


def read_neighbors(path="Rect10x10/pop.dat"):
    arq = open(path, "r")
    pops = []
    for line in arq.readlines():
        if line == "END":
            break
        tokens = line.split()
        pops.append(float(tokens[4]))
    arq.close()

    return pops


sample_arrivals = read_sample()
pops = read_neighbors()
total_pop = np.sum(pops)
probs = [x / total_pop for x in pops]

mn_samples = np.zeros((C, D, T, max_obs, S, R))
output_file = f"Rect10x10/mn_samples.dat"
arq = open(output_file, "w")
for c in range(C):
    for d in range(D):
        for t in range(T):
            for n in range(nb_obs[d]):
                print(c, d, t, n)
                mn = np.random.multinomial(sample_arrivals[c, d, t, n], probs, size=S)
                for s in range(S):
                    for r in range(R):
                        arq.write(f"{c} {d} {t} {n} {s} {r} {mn[s,r]}\n")
arq.write("END\n")
arq.close()
