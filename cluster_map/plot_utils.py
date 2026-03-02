import numpy as np
import matplotlib.pyplot as plt


def plot_cumulative_expression(expr, gene_name):
    expr_sorted = np.sort(expr)
    # plt.figure()
    plt.plot(expr_sorted, label=gene_name)
    plt.legend(bbox_to_anchor=(1.1, 1.05), ncol=3)
    plt.title("Cumulative normalized distribution")


#     plt.show()
