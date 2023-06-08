import pandas as pd
import sys

from hbip.utils.analysis_utils import plot_correlation, add_max_identity


def main(scores_path, title=None, add_max=False):
    pd.options.display.width = 100
    pd.options.display.max_columns = 10

    if add_max:
        df_max = add_max_identity(scores_path)
    else:
        df_max = pd.read_csv(scores_path)
    print(df_max.head())
    plot_correlation(df_max, fig_name=title, size="full", print_correlation=True)


if __name__ == "__main__":
    
    """
    Usage CLI: python3 analysis.py
    Usage CLI: python3 analysis.py path_to_scores title 
    Usage CLI: python3 analysis.py hbip_scores/alpha_beta_scores.csv "ab_new_correlation_annotated" 
    If no parameters are provided it will generate the scatterplot from alpha_beta_scores.csv

    """

    if len(sys.argv) == 1 :
        path = "hbip_scores/alpha_beta_scores.csv"
        title= "ab_correlation_annotated"
    elif len(sys.argv) > 1:
        path = sys.argv[1]
        title= sys.argv[2]

    boolean = False
    if pd.read_csv(path).shape[1] == 5:
        boolean = True
    main(scores_path = path, title=title, add_max=boolean)