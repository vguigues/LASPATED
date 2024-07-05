import pandas as pd

pd.options.mode.chained_assignment = "raise"
import geopandas as gpd

import matplotlib.pyplot as plt
import numpy as np
import laspated as spated
from shapely.geometry import Polygon, MultiPolygon, Point

from subprocess import run, CalledProcessError, STDOUT, PIPE, Popen
import os
from shutil import copyfile
import time


def process_land_use():
    land_use = gpd.read_file(r"regressores/uso_do_solo/")
    land_types = [
        "Afloramentos rochosos e depósitos sedimentares",
        "Áreas agrícolas",
        "Áreas de comércio e serviços",
        "Áreas de educação e saúde",
        "Áreas de exploração mineral",
        "Áreas de lazer",
        "Áreas de transporte",
        "Áreas industriais",
        "Áreas institucionais e de infraestrutura pública",
        "Áreas não edificadas",
        "Áreas residenciais",
        "Áreas sujeitas à inundação",
        "Cobertura arbórea e arbustiva",
        "Cobertura gramíneo lenhosa",
        "Corpos hídricos",
        "Favela",
    ]

    groups = {
        # Urbanized Areas. Public Activities
        0: [
            "Áreas de comércio e serviços",
            "Áreas de educação e saúde",
            "Áreas de lazer",
            "Áreas institucionais e de infraestrutura pública",
            "Áreas de transporte",
            "Áreas industriais",
        ],
        # Non Urbanized Areas. Non-Populational
        1: [
            "Afloramentos rochosos e depósitos sedimentares",
            "Áreas agrícolas",
            "Áreas de exploração mineral",
            "Áreas sujeitas à inundação",
            "Cobertura arbórea e arbustiva",
            "Cobertura gramíneo lenhosa",
            "Corpos hídricos",
        ],
        # Urbanized Areas. Residential
        2: ["Áreas residenciais", "Favela"],
    }

    for i in range(3):
        land_use[f"subgroup_{i}"] = 0
        for j, row in land_use.iterrows():
            if row["usoagregad"] in groups[i]:
                land_use.at[j, f"subgroup_{i}"] = 1

    land_use = land_use[["subgroup_0", "subgroup_1", "subgroup_2", "geometry"]].copy()
    return land_use


DISCS = {}


def generate_disc(disc_type):
    app = spated.DataAggregator(crs="epsg:4326")  # initializes data aggregator
    max_borders = gpd.read_file(r"rj/rj.shp")

    events = pd.read_csv(r"sorted_events.csv", encoding="ISO-8859-1", sep=",")
    app.add_events_data(
        events,
        datetime_col="data_hora",
        lat_col="lat",
        lon_col="long",
        feature_cols=["prioridade"],
        datetime_format="%m/%d/%y %H:%M:%S",
    )  # %m/%d/%y %H:%M:%S

    app.add_max_borders(max_borders)
    app.add_time_discretization("m", 30, 60 * 24, column_name="hhs")
    app.add_time_discretization("D", 1, 7, column_name="dow")

    if disc_type == "rect":
        app.add_geo_discretization(
            discr_type="R", rect_discr_param_x=10, rect_discr_param_y=10
        )
    elif disc_type == "hex":
        app.add_geo_discretization(discr_type="H", hex_discr_param=7)
    elif disc_type == "district":
        custom_map = gpd.read_file(
            r"rio_de_janeiro_neighborhoods/rio_neighborhoods.shp"
        )
        custom_map = custom_map.set_crs("epsg:29193")
        app.add_geo_discretization("C", custom_data=custom_map)

    land_use = process_land_use()
    app.add_geo_variable(land_use, type_geo_variable="area")

    population = gpd.read_file(r"regressores/populacao/")
    population = population[["population", "geometry"]].copy()
    app.add_geo_variable(population)

    # scaling the population
    app.geo_discretization["population"] /= 10**4
    R = int(np.max(app.events_data["gdiscr"]) + 1)
    app.plot_discretization(to_file=f"replication_results/disc_r{R}.pdf")

    app.write_arrivals(f"../Data/discretizations/{disc_type}/arrivals.dat")
    app.write_regions(f"../Data/discretizations/{disc_type}/neighbors.dat")
    app.write_info("dow", path=f"../Data/discretizations/{disc_type}/info.dat")

    num_regions = int(np.max(app.events_data["gdiscr"]) + 1)
    DISCS[num_regions] = app
    print(f"Wrote discretization files for disc {disc_type}. {num_regions} regions")
    return num_regions


def generate_discretizations():
    disc_types = ["rect", "hex", "district"]
    for disc_type in disc_types:
        generate_disc(disc_type)


def run_cpp_experiments():
    try:
        with Popen(
            ["cpp_tests/laspated", "-e", "all" , "--data_dir", "../Data"],
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


class Ex1Result(object):
    def __init__(
        self, nb_weeks: int, nb_groups: int, neighbor_factor: int, constant_lambda: int
    ) -> None:
        ex = 1 if constant_lambda == 1 else 2
        self.filename = f"cpp_tests/results/ex{ex}/err_by_weight_w{nb_weeks}_g{nb_groups}_a{neighbor_factor}_m{constant_lambda}.txt"
        fmt_neighbors = "no neighbors" if neighbor_factor == 0 else "neighbors"
        self.description = f"{nb_groups} groups, {fmt_neighbors}"
        file_handler = open(self.filename, "r")
        err_by_w = []
        first_line = file_handler.readline().split()
        self.min_error = float(first_line[0])
        self.min_w = int(first_line[1])
        n_weights = int(first_line[2])
        self.mean_emp = float(first_line[3])
        for i in range(n_weights):
            err_by_w.append(float(file_handler.readline().strip()))

        self.has_cv = nb_groups == 2 and neighbor_factor == 1 and nb_weeks > 1
        if self.has_cv:
            self.err_cv, self.w_cv = [float(x) for x in file_handler.readline().split()]
        file_handler.close()

        self.nb_weeks = nb_weeks
        self.nb_groups = nb_groups
        self.neighbor_factor = neighbor_factor
        self.constant_lambda = constant_lambda
        self.n_weights = n_weights
        self.err_by_w = err_by_w

        if nb_groups == 4 and neighbor_factor == 1:
            self.color = "black"
            self.line = "-"
        elif nb_groups == 4 and neighbor_factor == 0:
            self.color = "red"
            self.line = (0, (3, 5, 1, 5))
        elif nb_groups == 2 and neighbor_factor == 1:
            self.color = "blue"
            self.line = ":"
        elif nb_groups == 2 and neighbor_factor == 0:
            self.color = "black"
            self.line = (0, (5, 5))


PLOT_BASE_DIR = "replication_results/plots/"
TABLES_BASE_DIR = "replication_results/"


def experiment_1():
    num_weights = 13
    initial_weights = [[0] * num_weights, [0] * num_weights]
    initial_weights = [
        [0, 0.2, 0.4, 0.6, 0.8, 1, 1.2,1.4,1.6,1.8,2,50,100],
        [0, 0.2, 0.4, 0.6, 0.8, 1, 1.2,1.4,1.6,1.8,2,50,100]
    ]

    test_nb_weeks = [1, 10, 50, 500]
    test_nb_groups = [2, 4]
    test_neighbor_factors = [0, 1]
    test_constant_lambdas = [0,1]

    table_results = {}
    for nb_weeks in test_nb_weeks:
        for constant_lambda in test_constant_lambdas:
            x_axis = [w / nb_weeks for w in initial_weights[constant_lambda]]
            model = 1 if constant_lambda == 1 else 2
            fig, ax = plt.subplots(figsize=(12, 12))
            for nb_groups in test_nb_groups:
                for neighbor_factor in test_neighbor_factors:
                    result = Ex1Result(
                        nb_weeks, nb_groups, neighbor_factor, constant_lambda
                    )
                    if constant_lambda == 1:
                        table_results[nb_weeks, nb_groups, neighbor_factor] = (
                            result.err_by_w[result.min_w],
                            result.min_w,
                            result.mean_emp,
                            result.err_cv if result.has_cv else -1,
                            result.w_cv if result.has_cv else -1,
                        )
                    y_axis = result.err_by_w
                    plt.scatter(
                        x_axis[result.min_w],
                        result.err_by_w[result.min_w],
                        marker="d",
                        facecolors="none",
                        color=result.color,
                        s=160,
                    )
                    plt.plot(
                        x_axis,
                        y_axis,
                        linestyle=result.line,
                        label=result.description,
                        color=result.color,
                        linewidth=3,
                    )

                    if result.has_cv:
                        ax.axhline(
                            y=result.err_cv,
                            color="red",
                            linestyle="-",
                            label="Cross validation",
                        )
            ax.set_xlabel("Penalties", fontsize=18)
            ax.set_ylabel("Relative error", fontsize=18)
            ax.legend(fontsize=18)
            plot_filename = (
                PLOT_BASE_DIR + f"plot_err_by_weight_m{model}_obs{nb_weeks}.pdf"
            )
            plt.savefig(plot_filename, bbox_inches="tight")
            print(
                f"Plot of err by weights for model {model} with {nb_weeks} weeks saved at {plot_filename}"
            )
            plt.close()
    table_filename = "replication_results/tables/table_no_covariates_results.txt"
    table_file = open(table_filename, "w")
    table_file.write("Nb_weeks & Reg 1 & Reg 2 & Reg 3 & Reg 4 & CV & Emp\n")
    for nb_week in test_nb_weeks:
        table_file.write(f"{nb_week} & ")
        index_w = table_results[nb_week, 4, 1][1]
        val = initial_weights[1][index_w] / nb_week
        table_file.write(f"{table_results[nb_week,4,1][0]:.2f}/{val} & ")
        index_w = table_results[nb_week, 4, 0][1]
        val = initial_weights[1][index_w] / nb_week
        table_file.write(f"{table_results[nb_week,4,0][0]:.2f}/{val} & ")
        index_w = table_results[nb_week, 2, 1][1]
        val = initial_weights[1][index_w] / nb_week
        table_file.write(f"{table_results[nb_week,2,1][0]:.2f}/{val} & ")
        index_w = table_results[nb_week, 2, 0][1]
        val = initial_weights[1][index_w] / nb_week
        table_file.write(f"{table_results[nb_week,2,0][0]:.2f}/{val} & ")
        if nb_week > 1:
            table_file.write(
                f"{table_results[nb_week,2,1][3]:.2f}/{table_results[nb_week,2,1][4]} & "
            )
        else:
            table_file.write("- & ")
        table_file.write(f"{table_results[nb_week,2,1][2]}/0 \\\\ \hline\n")

        table_file.write("\n")
    table_file.close()
    print(f"No covariates experiments table saved at {table_filename}")


def experiment_2():
    regions = [1, 6]
    test_nb_weeks = [1, 10]
    for r in regions:
        for nb_weeks in test_nb_weeks:
            fig, ax = plt.subplots(figsize=(12, 12))
            periods = []
            theos = []
            emps = []
            ests = []
            legends = ["True intensity", "Empirical estimator", "Regularized estimator"]
            colors = ["black", "red", "blue"]
            lines = ["-", "dashed", "dashdot"]
            filename = f"cpp_tests/results/ex1/lambdas_r{r}_w{nb_weeks}.txt"
            file_handler = open(filename, "r")
            for line in file_handler.readlines():
                tokens = line.split()
                t = int(tokens[0])
                periods.append(t)
                theo, emp, est = [float(x) for x in tokens[1:]]
                theos.append(theo)
                emps.append(emp)
                ests.append(est)
            file_handler.close()
            plt.step(
                periods, theos, linestyle=lines[0], color=colors[0], label=legends[0]
            )
            plt.step(
                periods, emps, linestyle=lines[1], color=colors[1], label=legends[1]
            )
            plt.step(
                periods, ests, linestyle=lines[2], color=colors[2], label=legends[2]
            )
            ax.set_xlabel("Time periods", fontsize=18)
            ax.set_ylabel("Intensities", fontsize=18)
            ax.legend(fontsize=18)
            plt.savefig(f"replication_results/plots/art_rates_by_t_w{nb_weeks}_r{r}.pdf", bbox_inches="tight")
            plt.close()

def experiment_3():
    copyfile(
        "cpp_tests/results/ex3/table_ex3.txt",
        "replication_results/tables/table_covariates_results.txt",
    )

    table_file_in = open("cpp_tests/results/ex3/table_ex3.txt", "r")
    contents = []
    for line in table_file_in.readlines():
        contents.append(line.split())
    table_file_in.close()
    table_file_out = open("replication_results/tables/table_covariates_results.txt")

    table_file_out.close()
    print(
        "Table with results for artificial data with covariates saved at replication_results/tables/table_covariates_results.txt"
    )


def experiment_4():
    C = 3
    T = 7 * 48
    regions = DISCS.keys()
    # regions = [76]
    description = {}
    for r in regions:
        if r == min(regions):
            description[r] = "Rectangular 10x10"
        elif r == max(regions):
            description[r] = "Hexagonal 7"
        else:
            description[r] = "District"
    for R in regions:
        tr_filename = f"cpp_tests/results/real_data/rates_by_t_r{R}.txt"
        tr_file_handler = open(tr_filename, "r")
        emp, reg, cov = np.zeros((C, T)), np.zeros((C, T)), np.zeros((C, T))
        for line in tr_file_handler.readlines():
            tokens = line.split()
            t = int(tokens[0])
            i = 1
            for c in range(C):
                emp[c, t] = float(tokens[i])
                reg[c, t] = float(tokens[i + 1])
                cov[c, t] = float(tokens[i + 2])
                i += 3
        tr_file_handler.close()
        x_axis = [x for x in range(T)]
        labels = ["Empirical", description[R], f"Covariates, {description[R].lower()}"]
        colors = ["black", "red", "pink"]
        lines = ["-", (0, (5, 5)), ":"]
        # Total
        y_emp_axis = [np.sum(emp[:, t]) for t in range(T)]
        y_reg_axis = [np.sum(reg[:, t]) for t in range(T)]
        y_cov_axis = [np.sum(cov[:, t]) for t in range(T)]
        fig_total, ax_total = plt.subplots(figsize=(12, 10))
        plt.plot(
            x_axis,
            y_emp_axis,
            label=labels[0],
            linewidth=3,
            color=colors[0],
            linestyle=lines[0],
        )
        plt.plot(
            x_axis,
            y_reg_axis,
            label=labels[1],
            linewidth=3,
            color=colors[1],
            linestyle=lines[1],
        )
        plt.plot(
            x_axis,
            y_cov_axis,
            label=labels[2],
            linewidth=3,
            color=colors[2],
            linestyle=lines[2],
        )
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        ax_total.set_xlabel("Time steps", fontsize=20)
        ax_total.set_ylabel("Intensities", fontsize=20)
        ax_total.legend(loc="upper right", fontsize=20, framealpha=0.5)
        plot_filename = PLOT_BASE_DIR + f"plot_rates_by_t_R{R}_total.pdf"
        plt.savefig(plot_filename, bbox_inches="tight")
        print(f"Saved plot by time rates at {plot_filename}")
        plt.close()
        # c = 0
        y_emp_axis = [emp[0, t] for t in range(T)]
        y_reg_axis = [reg[0, t] for t in range(T)]
        y_cov_axis = [cov[0, t] for t in range(T)]
        plt.clf()
        plt.cla()
        fig_total, ax_total = plt.subplots(figsize=(12, 10))
        plt.plot(
            x_axis,
            y_emp_axis,
            label=labels[0],
            linewidth=3,
            color=colors[0],
            linestyle=lines[0],
        )
        plt.plot(
            x_axis,
            y_reg_axis,
            label=labels[1],
            linewidth=3,
            color=colors[1],
            linestyle=lines[1],
        )
        plt.plot(
            x_axis,
            y_cov_axis,
            label=labels[2],
            linewidth=3,
            color=colors[2],
            linestyle=lines[2],
        )
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        ax_total.set_xlabel("Time steps", fontsize=20)
        ax_total.set_ylabel("Intensities", fontsize=20)
        ax_total.legend(loc="upper right", fontsize=20, framealpha=0.5)
        plot_filename = PLOT_BASE_DIR + f"plot_rates_by_t_R{R}_c0.pdf"
        plt.savefig(plot_filename, bbox_inches="tight")
        print(f"Saved plot by time rates at {plot_filename}")
        plt.close()
        # c = 1
        y_emp_axis = [emp[1, t] for t in range(T)]
        y_reg_axis = [reg[1, t] for t in range(T)]
        y_cov_axis = [cov[1, t] for t in range(T)]

        plt.clf()
        plt.cla()
        fig_total, ax_total = plt.subplots(figsize=(12, 10))
        plt.plot(
            x_axis,
            y_emp_axis,
            label=labels[0],
            linewidth=3,
            color=colors[0],
            linestyle=lines[0],
        )
        plt.plot(
            x_axis,
            y_reg_axis,
            label=labels[1],
            linewidth=3,
            color=colors[1],
            linestyle=lines[1],
        )
        plt.plot(
            x_axis,
            y_cov_axis,
            label=labels[2],
            linewidth=3,
            color=colors[2],
            linestyle=lines[2],
        )
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        ax_total.set_xlabel("Time steps", fontsize=20)
        ax_total.set_ylabel("Intensities", fontsize=20)
        ax_total.legend(loc="upper right", fontsize=20, framealpha=0.5)
        plot_filename = PLOT_BASE_DIR + f"plot_rates_by_t_R{R}_c1.pdf"
        plt.savefig(plot_filename, bbox_inches="tight")
        print(f"Saved plot by time rates at {plot_filename}")
        plt.close()
        # c = 2
        y_emp_axis = [emp[2, t] for t in range(T)]
        y_reg_axis = [reg[2, t] for t in range(T)]
        y_cov_axis = [cov[2, t] for t in range(T)]
        plt.close()
        fig_total, ax_total = plt.subplots(figsize=(12, 10))
        plt.plot(
            x_axis,
            y_emp_axis,
            label=labels[0],
            linewidth=3,
            color=colors[0],
            linestyle=lines[0],
        )
        plt.plot(
            x_axis,
            y_reg_axis,
            label=labels[1],
            linewidth=3,
            color=colors[1],
            linestyle=lines[1],
        )
        plt.plot(
            x_axis,
            y_cov_axis,
            label=labels[2],
            linewidth=3,
            color=colors[2],
            linestyle=lines[2],
        )
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        ax_total.set_xlabel("Time steps", fontsize=20)
        ax_total.set_ylabel("Intensities", fontsize=20)
        ax_total.legend(loc="upper right", fontsize=20, framealpha=0.5)
        plot_filename = PLOT_BASE_DIR + f"plot_rates_by_t_R{R}_c2.pdf"
        plt.savefig(plot_filename, bbox_inches="tight")
        print(f"Saved plot by time rates at {plot_filename}")
        plt.close()

        if R != min(regions):
            continue

        fig_total, ax_total = plt.subplots(figsize=(12, 10))
        # Heatmaps
        rr_filename = f"cpp_tests/results/real_data/rates_by_region_r{R}.txt"
        rr_file_handler = open(rr_filename, "r")
        emp, reg, cov = np.zeros((C, R)), np.zeros((C, R)), np.zeros((C, R))
        for line in rr_file_handler.readlines():
            tokens = line.split()
            r = int(tokens[0])
            i = 1
            for c in range(C):
                emp[c, r] = float(tokens[i])
                reg[c, r] = float(tokens[i + 1])
                cov[c, r] = float(tokens[i + 2])
                i += 3
        rr_file_handler.close()
        gd = DISCS[R].geo_discretization
        gd["emp_c0"] = emp[0, :]
        gd["reg_c0"] = reg[0, :]
        gd["cov_c0"] = cov[0, :]
        gd["emp_c1"] = emp[1, :]
        gd["reg_c1"] = reg[1, :]
        gd["cov_c1"] = cov[1, :]
        gd["emp_c2"] = emp[2, :]
        gd["reg_c2"] = reg[2, :]
        gd["cov_c2"] = cov[2, :]
        gd["emp_total"] = [np.sum(emp[:, r]) for r in range(R)]
        gd["reg_total"] = [np.sum(reg[:, r]) for r in range(R)]
        gd["cov_total"] = [np.sum(cov[:, r]) for r in range(R)]

        plt.clf()
        plt.cla()
        fig_size = (12, 10)
        # Total
        heat_filename = PLOT_BASE_DIR + f"heat_total_r{R}_emp.pdf"
        gd.plot(column="emp_total", legend=True, cmap="OrRd", figsize=fig_size)
        plt.savefig(heat_filename, bbox_inches="tight")
        print(f"Heatmap of total intensities model 1 saved at {heat_filename}")
        heat_filename = PLOT_BASE_DIR + f"heat_total_r{R}_m1.pdf"
        plt.close()
        gd.plot(column="reg_total", legend=True, cmap="OrRd", figsize=fig_size)
        plt.savefig(heat_filename, bbox_inches="tight")
        print(f"Heatmap of total intensities model 1 saved at {heat_filename}")
        plt.close()
        heat_filename = PLOT_BASE_DIR + f"heat_total_r{R}_m2.pdf"
        gd.plot(column="cov_total", legend=True, cmap="OrRd", figsize=fig_size)
        plt.savefig(heat_filename, bbox_inches="tight")
        print(f"Heatmap of total intensities model 2 saved at {heat_filename}")
        plt.close()
        # c = 0
        heat_filename = PLOT_BASE_DIR + f"heat_c0_r{R}.pdf"
        gd.plot(column="reg_c0", legend=True, cmap="OrRd", figsize=fig_size)
        plt.savefig(heat_filename, bbox_inches="tight")
        print(f"Heatmap c = 0, model 1 saved at {heat_filename}")
        heat_filename = PLOT_BASE_DIR + f"heat_c0_r{R}.pdf"
        plt.close()
        gd.plot(column="cov_c0", legend=True, cmap="OrRd", figsize=fig_size)
        plt.savefig(heat_filename, bbox_inches="tight")
        print(f"Heatmap c = 0, model 2 saved at {heat_filename}")
        plt.close()
        # c = 1
        heat_filename = PLOT_BASE_DIR + f"heat_c1_r{R}.pdf"
        gd.plot(column="reg_c1", legend=True, cmap="OrRd", figsize=fig_size)
        plt.savefig(heat_filename, bbox_inches="tight")
        print(f"Heatmap c = 1, model 1 saved at {heat_filename}")
        plt.close()
        heat_filename = PLOT_BASE_DIR + f"heat_c1_r{R}.pdf"
        gd.plot(column="cov_c1", legend=True, cmap="OrRd", figsize=fig_size)
        plt.savefig(heat_filename, bbox_inches="tight")
        print(f"Heatmap c = 1, model 2 saved at {heat_filename}")
        plt.close()
        # c = 2
        heat_filename = PLOT_BASE_DIR + f"heat_c2_r{R}.pdf"
        gd.plot(column="reg_c2", legend=True, cmap="OrRd", figsize=fig_size)
        plt.savefig(heat_filename, bbox_inches="tight")
        print(f"Heatmap c = 2, model 1 saved at {heat_filename}")
        plt.close()
        heat_filename = PLOT_BASE_DIR + f"heat_c2_r{R}.pdf"
        gd.plot(column="cov_c2", legend=True, cmap="OrRd", figsize=fig_size)
        plt.savefig(heat_filename, bbox_inches="tight")
        print(f"Heatmap c = 2, model 2 saved at {heat_filename}")
        plt.close()


def plot_results():
    experiment_1()
    experiment_2()
    experiment_3()
    experiment_4()


def main():
    tic = time.perf_counter()
    print("Generating discretizations")
    generate_discretizations()
    # print(DISCS.keys())
    print("Running experiments")
    run_cpp_experiments()

    print("Plotting results")
    plot_results()
    toc = time.perf_counter()
    print(f"Finished replication script in {toc - tic:0.4f} seconds")


if __name__ == "__main__":
    main()
