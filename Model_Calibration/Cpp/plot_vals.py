#!/usr/bin/env python
#-*- encoding: utf-8 -*-


import matplotlib.pyplot as plt
import sys

def get_values(path):
    return [float(x) for x in open(path,"r").readlines()]


def plot_values(path):
    values = get_values(path)
    plt.plot(values)
    plt.show()



def plot_many_vals(nb_obs):
    cv_no_reg = [-1, 0.7577456, 0.5297131, 0.0883604]
    files =  ["err_no_reg_g4_obs%d_n2_a1", "err_no_reg_g4_obs%d_n2_a0", "err_no_reg_g2_obs%d_n2_a1", "err_no_reg_g2_obs%d_n2_a0"]
    legends = ["4 groups, neighbors", "4 groups, no neighbors", "2 groups, neighbors", "2 groups, no neighbors"]

    x_axis = [0,0.01,0.05,0.01,0.1,1,5,10,20,30,40,50,60,70,80,90,100,
			110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,
			270,280,290,300,310,320,330,340,350,360,370,380,390,400]
    
    colors = ["black", "red", "blue", "black"]
    lines  = ["-", (0, (3, 5, 1, 5)), ":", (0, (5, 10))]
    
    fig, ax = plt.subplots(figsize=(12, 6))
    for i,file in enumerate(files):
        fp = file % (nb_obs)
        try:
            values = get_values(fp+".txt")
            plt.plot(x_axis, values, linestyle=lines[i],label=legends[i], linewidth=3, color=colors[i])
            ax.set_xlabel("Penalties", fontsize=18)
            ax.set_ylabel("Relative error", fontsize=18)
            if cv_no_reg[i] >= 0:
                plt.axhline(y=cv_no_reg[i], color="red", linestyle="-")
            if nb_obs == 1000:
                ax.legend(loc="lower left",fontsize=15)
            else:
                ax.legend(fontsize=15)
        except ValueError:
            print("Error at %s: x =" % fp,x_axis, get_values(fp+".txt"))
            input()
    

    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.savefig("results_n%d_crop.pdf" % (nb_obs), bbox_inches="tight")
    # plt.show()


def plot_with_obs():
    all_nbs = [1,10,50,1000]

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