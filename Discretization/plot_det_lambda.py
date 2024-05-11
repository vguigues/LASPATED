import numpy as np
import matplotlib.pyplot as plt


class Ex1Result(object):
    def __init__(
        self, nb_weeks: int, nb_groups: int, neighbor_factor: int, constant_lambda: int
    ) -> None:
        ex = 1 if constant_lambda == 1 else 2
        self.filename = f"cpp_tests/results/ex{ex}/det_err_by_weight_w{nb_weeks}_g{nb_groups}_a{neighbor_factor}_m{constant_lambda}.txt"
        fmt_neighbors = "no neighbors" if neighbor_factor == 0 else "neighbors"
        self.description = f"{nb_groups} groups, {fmt_neighbors}"
        file_handler = open(self.filename, "r")
        err_by_w = []
        first_line = file_handler.readline().split()
        self.min_error = float(first_line[0])
        self.min_w = int(first_line[1])
        n_weights = int(first_line[2])
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


PLOT_BASE_DIR = ""
TABLES_BASE_DIR = ""


def experiment_1():
    # num_weights = 10
    # initial_weights = [[0] * num_weights, [0] * num_weights]
    # for w in range(num_weights):
    #     # initial_weights[0][w] = 0.001 * w
    #     initial_weights[1][w] = 0.2 * w
    initial_weights = [
        [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0],
        [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0],
    ]
    test_nb_weeks = [1]
    test_nb_groups = [2, 4]
    test_neighbor_factors = [0, 1]
    test_constant_lambdas = [1]

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
                    if constant_lambda == 0:
                        table_results[nb_weeks, nb_groups, neighbor_factor] = (
                            result.err_by_w[result.min_w],
                            result.min_w,
                            result.err_by_w[0],
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
                PLOT_BASE_DIR + f"plot_det_err_by_weight_m{model}_obs{nb_weeks}.pdf"
            )
            plt.savefig(plot_filename, bbox_inches="tight")
            print(
                f"Plot of err by weights for model {model} with {nb_weeks} weeks saved at {plot_filename}"
            )
    # table_filename = "table_det_no_covariates_results.txt"
    # table_file = open(table_filename, "w")
    # table_file.write("Nb_weeks\tReg 1\tReg 2\tReg 3\tReg 4\t CV\tEmp\n")
    # for nb_week in test_nb_weeks:
    #     table_file.write(f"{nb_week}\t")
    #     table_file.write(
    #         f"{table_results[nb_week,4,1][0]}/{table_results[nb_week,4,1][1]}\t"
    #     )
    #     table_file.write(
    #         f"{table_results[nb_week,4,0][0]}/{table_results[nb_week,4,0][1]}\t"
    #     )
    #     table_file.write(
    #         f"{table_results[nb_week,2,1][0]}/{table_results[nb_week,2,1][1]}\t"
    #     )
    #     table_file.write(
    #         f"{table_results[nb_week,2,0][0]}/{table_results[nb_week,2,0][1]}\t"
    #     )
    #     if nb_week > 1:
    #         table_file.write(
    #             f"{table_results[nb_week,2,1][3]}/{table_results[nb_week,2,1][4]}\t"
    #         )
    #     else:
    #         table_file.write("-\t")
    #     table_file.write(f"{table_results[nb_week,2,1][2]}/0\t")

    #     table_file.write("\n")
    # table_file.close()
    # print(f"No covariates experiments table saved at {table_filename}")


def experiment_2():
    regions = [1, 6]
    test_nb_weeks = [1]
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
            filename = f"cpp_tests/results/ex1/det_lambdas_r{r}_w{nb_weeks}.txt"
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
            plt.savefig(f"art_det_rates_by_t_w{nb_weeks}_r{r}.pdf", bbox_inches="tight")


experiment_1()
experiment_2()
