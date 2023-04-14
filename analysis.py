import pandas as pd

from hbip.utils.analysis_utils import plot_correlation, add_max_identity


def main(add_max=False):
    pd.options.display.width = 100
    pd.options.display.max_columns = 10

    # path = "hbip_scores/alpha_beta_scores.csv"
    path = "hbip_scores/alpha_beta_new_scores.csv"
    if add_max:
        df_max = add_max_identity(path)
    else:
        df_max = pd.read_csv(path)
    print(df_max.head())
    # plot_correlation(df_max, fig_name="ab_correlation_annotated", size="full")
    plot_correlation(df_max, fig_name="ab_correlation_annotated, updated binding", size="full")


if __name__ == "__main__":
    main(add_max=True)
