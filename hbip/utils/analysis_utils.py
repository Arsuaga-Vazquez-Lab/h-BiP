import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt, patches as patches

from binding import bind_human, no_bind_human


def scores_binding(df):
    subset = df[df.Accession.isin(bind_human)].copy()
    subset["Binds"] = "Yes"
    temp = df[df.Accession.isin(no_bind_human)].copy()
    temp["Binds"] = "No"
    return pd.concat([subset, temp])


def add_max_identity(path):
    df = pd.read_csv("./data/alpha_beta_identity_vs_hcov.csv")
    scores = pd.read_csv(path)
    if "MaxSim" not in set(scores.columns):
        df = df.merge(scores, on="Accession")
        df.to_csv(path, index=False)
    else:
        print("MaxSim already in dataset")
    return df


def plot_correlation(df, fig_name=None, size="half", print_correlation=False):
    """size refers to letter
    fig_name: string
        file name to save the plot. if None, it won't be saved"""
    if size == "full":
        figsize = (6, 4)
    elif size == "half":
        figsize = (3.3, 2.5)
    else:  # best for display
        figsize = (10, 6)

    fig, ax = plt.subplots(figsize=figsize)
    plt.rc("font", size=7)  # ab_correlation_annotated, text inside plot

    sns.scatterplot(
        x="MaxSim",
        y="Binds_prob",
        data=df,
        hue="Binds",
        palette=["dimgray", "indianred"],
        style="Binds",
        s=5,
    )  # smaller
    lgnd = plt.legend(labels=["Binds", "Doesn't bind"], bbox_to_anchor=(1, 1.1), ncol=2)
    lgnd.legendHandles[0]._sizes = [40]
    lgnd.legendHandles[1]._sizes = [30]
    plt.xlabel("Max % identity hCoV")
    plt.ylabel("Human Binding Potential score")
    plt.plot(97.46, 0.999, color="indianred", marker="*", markersize=5)  # RaTG13
    plt.axhline(y=0.5, linestyle="--")
    plt.axvline(x=97, linestyle="--")
    plt.text(61.5, 0.78, "Bt133")  # 67.18, 0.78
    plt.text(84, 0.82, "LYRa3")  # 89.7, 0.82
    rect = patches.Rectangle(
        (98, 0.67), 2.3, 0.18, alpha=0.5, edgecolor="indianred", facecolor="indianred"
    )
    ax.add_patch(rect)

    if fig_name != None:
        plt.savefig("./outputs/" + fig_name + ".pdf", bbox_inches="tight")

    plt.show()

    if print_correlation:
        print("\nCorrelation {:.2f}\n".format(df.MaxSim.corr(df.Binds_prob)))
