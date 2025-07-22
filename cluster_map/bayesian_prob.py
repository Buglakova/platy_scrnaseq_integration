import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd
import logging
from pathlib import Path
import numpy as np


def get_prob(expr, profile):
    positive_profile = profile
    negative_profile = 1 - profile

    positive_expr = expr
    negative_expr = 1 - expr

    prob = (positive_expr * positive_profile + negative_expr * negative_profile)
    
    # logging.info("______")
    # logging.info(profile)
    # logging.info(prob)
    # logging.info("______")
    return prob.prod()


def expr_to_cluster_probability(expr, cluster_prior_prob, cluster_profile, total_profile):
    cluster_prob = get_prob(expr, cluster_profile)
    total_prob = get_prob(expr, total_profile)
    bayesian_prob = cluster_prob * cluster_prior_prob
    return bayesian_prob


def assign_clusters(smfish_df, gene_list, profiles, total_profile, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(exist_ok=True)

    cluster_list = list(profiles.keys())
    result_df = pd.DataFrame(index=smfish_df.index, columns=cluster_list)
    

    # Gene order in smfish_df and in cluster expression profiles are assumed to be the same

    # Prior probability: each cluster has same prior probability
    cluster_prior_prob = 1 / len(cluster_list)

    for cluster_name in cluster_list:
        cluster_profile = profiles[cluster_name]

        bayesian_prob = [expr_to_cluster_probability(expr.to_numpy().flatten(), cluster_prior_prob, cluster_profile, total_profile) for _, expr in smfish_df.iterrows()]
        result_df[cluster_name] = bayesian_prob

    # result_df = result_df / (result_df.sum(axis=1) + 1e-7)

    logging.info(f"Highest cluster probaility: {result_df.max()}")
    logging.info(f"Probability sum: {result_df.sum(axis=1)}")
    logging.info(f"Probability sum min: {result_df.sum(axis=1).min()}")

    result_df = result_df.div(result_df.sum(axis=1), axis=0)

    logging.info(f"Highest cluster probaility: {result_df.max()}")
    logging.info(f"Probability sum min: {result_df.sum(axis=1).min()}")

    result_df["most probable cluster"] = result_df.idxmax(axis=1)

    result_df.to_csv(out_dir / "cluster_probability.tsv", sep="\t")
