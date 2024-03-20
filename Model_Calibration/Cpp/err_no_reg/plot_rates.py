#!/usr/bin/env python
#-*- encoding: utf-8 -*-


import matplotlib.pyplot as plt
import sys
from operator import itemgetter
import numpy as np


C = 3
T = 7*48

# NOREG_FILES = ["rates_no_reg_real76_g336_obs105_n1.txt",
#                "rates_no_reg_real297_g336_obs105_n1.txt",
#                "rates_no_reg_real160_g336_obs105_n1.txt",
#                "rates_no_reg_real1811_g336_obs105_n1.txt"]

NOREG_FILES = ["rates_no_reg_real297_g336_obs105_n1.txt"]

# NOREG_LEGENDS = ["Rectangular 10x10",
#                  "Hexagonal 7",
#                  "District",
#                  "Hexagonal 8"]

NOREG_LEGENDS = ["Hexagonal 7"]

# NOREG_PLOTS = ["AllstepsRectangular%s_v2.pdf",
#                "AllstepsHexagonal7%s_v2.pdf",
#                "AllstepsDistrict%s_v2.pdf",
#                "AllstepsHexagonal%s_v2.pdf"]

NOREG_PLOTS = ["AllstepsHexagonal7%s_v2.pdf"]

# REG_FILES = ["rates_reg_real76_g336_obs105_n1.txt",
#              "rates_reg_real297_g336_obs105_n1.txt",
#                "rates_reg_real160_g336_obs105_n1.txt",
#                "rates_reg_real1811_g336_obs105_n1.txt"]

REG_FILES = ["rates_reg_real297_g336_obs105_n1.txt"]

# REG_LEGENDS = ["Covariates, rectangular 10x10",
#                "Covariates, hexagonal 7",
#                  "Covariates, district",
#                  "Covariates, hexagonal 8"]

REG_LEGENDS = ["Covariates, hexagonal 7"]


def get_reg_values(path):
    reg_values = np.zeros((C,T))
    est_values = np.zeros((C,T))
    arq = open(path, "r") 
    for line in arq.readlines():
        c,t,cval, est_val = [float(x) for x in line.split()]
        c,t = int(c), int(t)
        reg_values[c,t] = cval
        est_values[c,t] = est_val
        print("REG "+line)
    arq.close()
    input()
    return reg_values, est_values

def get_values(path):
    estimated = np.zeros((C,T))
    cv = np.zeros((C,T))

    arq = open(path, "r")
    for line in arq.readlines():
        c,t,cval,est = [float(x) for x in line.split()]
        c,t = int(c), int(t)
        estimated[c,t] = est
        cv[c,t] = cval
        print("NOREG "+line)
    arq.close()
    input()
    return estimated,cv


def get_total(values):
    total_values = np.zeros(T)
    for t in range(T):
        total_values[t] = np.sum(values[:,t])
    return total_values

def sub_plot(i, c, est_values, cv_values, cov_values):
    fig_total, ax_total = plt.subplots(figsize=(12, 10))
    plt.plot([t for t in range(T)], est_values, linestyle="-", 
             label="Empirical", linewidth=3, color="black")
    plt.plot([t for t in range(T)], cv_values, 
             linestyle= (0, (5, 5)), label=NOREG_LEGENDS[i],  
             linewidth=3, color="red")
    plt.plot([t for t in range(T)], cov_values, 
             linestyle= ":", label=REG_LEGENDS[i],  
             linewidth=3, color="violet")
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    ax_total.set_xlabel("Time steps", fontsize=20)
    ax_total.set_ylabel("Intensities", fontsize=20)
    ax_total.legend(loc="upper right", fontsize=20, framealpha=0.5)
    plt.savefig(NOREG_PLOTS[i] % ("P"+str(c+1)+"b"), bbox_inches="tight")

def plot():

    for i,rates_file in enumerate(NOREG_FILES):
        estimated, cv = get_values(rates_file)
        cov_values, estimated = get_reg_values(REG_FILES[i])
        total_est = get_total(estimated)
        total_cv  = get_total(0.333*cv)
        total_cov = get_total(cov_values)

        print(f"File: {rates_file}:")
        for c in range(C):
            for t in range(T):
                print(f"c = {c}, t = {t}, emp = {estimated[c,t]}, no_reg = {0.333*cv[c,t]}, reg = {cov_values[c,t]}")
            input()

        fig_total, ax_total = plt.subplots(figsize=(12, 10))
        plt.plot([t for t in range(T)], total_est, 
                 linestyle="-", label="Empirical", 
                 linewidth=3, color="black")
        plt.plot([t for t in range(T)], total_cv, 
                 linestyle= (0, (5, 5)), label=NOREG_LEGENDS[i],  
                 linewidth=3, color="red")
        plt.plot([t for t in range(T)], total_cov, 
                 linestyle=":",label=REG_LEGENDS[i], 
                 linewidth=3, color="violet")
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        ax_total.set_xlabel("Time steps", fontsize=20)
        ax_total.set_ylabel("Intensities", fontsize=20)
        ax_total.legend(loc="upper right", fontsize=20, framealpha=0.5)
        plt.savefig(NOREG_PLOTS[i] % ("b"), bbox_inches="tight")

        for c in range(C):
            sub_plot(i, c, estimated[c,:], 0.333*cv[c,:], cov_values[c,:])

plot()
