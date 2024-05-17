#!/usr/bin/env python
#-*- encoding: utf-8 -*-


import matplotlib.pyplot as plt
import sys
import numpy as np 

C = 1
T = 28
G = 2


def read_file(filename):
    arq = open(filename,"r")
    rates_true = np.zeros((C,T))
    rates_emp = np.zeros((C,T))
    rates_reg = np.zeros((C,T))
    for line in arq.readlines():
        tokens = line.split()
        c,t = [int(x) for x in tokens[:2]]
        rate_true, rate_emp, rate_reg = [float(x) for x in tokens[2:]]
        rates_true[c,t] = rate_true
        rates_emp[c,t] = rate_emp
        rates_reg[c,t] = rate_reg
    arq.close()

    return rates_true, rates_emp,rates_reg

def plot_rates(nb_obs,r):
    rates_true, rates_emp, rates_reg = read_file("x_no_reg_r%d_%d_%d.txt" % (r,G,nb_obs))
    fig, ax = plt.subplots(figsize=(12, 12))
    labels = [""]
    lines  = ["solid", "dashed", "dashdot"]
    colors = ["black", "blue", "red"]
    legends = ["True intensity", "Regularized estimator", "Empirical estimator"]
    values = [rates_true[0,:], rates_reg[0,:],rates_emp[0,:]]
    x_axis = [x for x in range(28)]
    # ax.legend(fontsize=15)
    ax.set_xlabel("Time periods", fontsize=18)
    ax.set_ylabel("Intensities", fontsize=18)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlim([0, max(x_axis)])
    plt.ylim([0, max([max(x) for x in values])])
    for i in range(3):
        plt.step(x_axis,values[i], linestyle=lines[i],label=legends[i],linewidth=3, color=colors[i])

    ax.legend(fontsize=18)
    plt.savefig("rates%dsite%dv3.pdf" % (nb_obs,r), bbox_inches="tight")
    

plot_rates(2,1)
plot_rates(2,6)
plot_rates(20,1)
plot_rates(20,6)






