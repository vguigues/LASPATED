#!/usr/bin/env python
#-*- encoding: utf-8 -*-


import matplotlib.pyplot as plt
import sys
from operator import itemgetter

def get_values(path):
    return [float(x) for x in open(path,"r").readlines()]


def plot_values(path):
    values = get_values(path)
    plt.plot(values)
    plt.show()



def plot_many_vals(nb_obs):
    # Constant lambdas:
    cv_no_reg = {2: -1, 20: 0.196737, 100: 0.164541, 1000: 0.105242}
    # Linear Lambdas:
    # cv_no_reg = {2: -1, 20: 0.247164, 100: 0.174264, 1000: 0.105847}
    files =  ["err_no_reg_g4_obs%d_n2", "err_no_reg_g4_obs%d_n2_a0", "err_no_reg_g2_obs%d_n2", "err_no_reg_g2_obs%d_n2_a0"]
    legends = ["4 groups, neighbors", "4 groups, no neighbors", "2 groups, neighbors", "2 groups, no neighbors"]

    # x_axis = [0,0.01,0.05,0.01,0.1,1,5,10,20,30,40,50,60,70,80,90,100,
	# 		110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,
	# 		270,280,290,300,310,320,330,340,350,360,370,380,390,400]
    # x_axis = [0,0.01,0.05,0.1,1,5,10,20,30,40,50,60,70,80,90,100]
    x_axis = [0,0.3,0.6, 0.9,1.2,1.5,1.8,2.1,2.4,2.7,3]
    if nb_obs == 20:
        x_axis = [x / 10 for x in x_axis]
    elif nb_obs == 100:
        x_axis = [x / 50 for x in x_axis]
    elif nb_obs == 1000:
        x_axis = [x / 500 for x in x_axis]
    colors = ["black", "red", "blue", "black"]
    lines  = ["-", (0, (3, 5, 1, 5)), ":", (0, (5, 5))]
    
    fig, ax = plt.subplots(figsize=(12, 12))
    max_value = 0
    if not nb_obs in [1,2]:
        ax.axhline(y = cv_no_reg[nb_obs], color = "red", linestyle = "-", label="Cross validation")
        ax.legend(fontsize=15)
    for i,file in enumerate(files):
        fp = file % (nb_obs)
        try:
            values = get_values(fp+".txt")
            max_value = max([max_value] + values)
            f = [(x_axis[i], values[i]) for i in range(len(x_axis))]
            min_point = min(f, key=itemgetter(1))
            plt.scatter(min_point[0], min_point[1], marker="d", facecolors="none", color=colors[i], s=160)
            plt.plot(x_axis, values, linestyle=lines[i],label=legends[i], linewidth=3, color=colors[i])
            ax.set_xlabel("Penalties", fontsize=18)
            ax.set_ylabel("Relative error", fontsize=18)
            ax.legend(fontsize=18)
        except ValueError:
            print("Error at %s: x =" % fp,x_axis, get_values(fp+".txt"))
            input()
    
        
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlim([0, max(x_axis)])
    plt.ylim([0, max_value])
    plt.savefig("%dweeks_Color_Unknown_v3.pdf" % (nb_obs*2 if nb_obs in [1,10] else nb_obs), bbox_inches="tight")
    # plt.show()


def plot_with_obs():
    all_nbs = [2,20, 100,1000]

    for nb_obs in all_nbs:
        plot_many_vals(nb_obs)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Wrong usage! Correct usage: python plot_vals.py path|many|obs")
    elif sys.argv[1] == "many":
        plot_many_vals()
    elif sys.argv[1] == "obs":
        plot_with_obs()
    else:
        plot_values(sys.argv[1])