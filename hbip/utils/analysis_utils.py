import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt, patches as patches
from sklearn.metrics import classification_report, roc_curve, roc_auc_score

# Coordinates are best display for the label (not the exact values)
FP_orig =  {'Bt133': (61.5, 0.78),
            'LYRa3': (84, 0.82)}
FP_mers_mod =  {'Bt133': (61.5, 0.78),
                'LYRa3': (84, 0.79),
                'HKU5r': (58, 0.99)}


def add_max_identity(path):
    df = pd.read_csv("./data/alpha_beta_identity_vs_hcov.csv")
    scores = pd.read_csv(path)
    if "MaxSim" not in set(scores.columns):
        df = df.merge(scores, on="Accession")
        df.to_csv(path, index=False)
    else:
        print("MaxSim already in dataset")
    return df


def plot_correlation(df, fig_name=None, size="half", print_correlation=False, 
                     africa=False, show_viruses=FP_mers_mod):
    """size refers to letter
    fig_name: string
        file name to save the plot. if None, it won't be saved
    size: string
        Plot size: "full" for full letter, "half" for half letter (default), 
        "display" best for display
    africa: boolean
        It will add a red rectangle to signal MERS Africa for the original model
    show_viruses: dict
        Viruses to display the label in the plot {'virus_name': (x, y)}
    """
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
    plt.axvline(x=93, linestyle="--")
    if show_viruses != None:
        for virus, coord in show_viruses.items():
            plt.text(coord[0], coord[1], virus)
    if africa:
        rect = patches.Rectangle(
            (98, 0.67), 2.3, 0.18, alpha=0.5, edgecolor="indianred", facecolor="indianred"
        )
        ax.add_patch(rect)

    if fig_name != None:
        plt.savefig("./outputs/" + fig_name + ".pdf", bbox_inches="tight")

    plt.show()

    if print_correlation:
        print("\nCorrelation {:.2f}\n".format(df.MaxSim.corr(df.Binds_prob)))

def metrics(df, label, plot_roc_curve=False):
    def plot_curve(fp, tp):
        plt.figure()
        plt.plot(fp, tp, color="darkorange", lw=2)
        plt.plot([0, 1], [0, 1], color="navy", lw=2, linestyle="--")
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate")
        plt.show()

    _prob = label + "_prob"
    _pred = label + "_predicted"
    print(classification_report(df[label], df[_pred]))
    print("Note: Sensitivity=Recall(1), Specificity=Recall(0), FPR=1-Recall(0)")

    if plot_roc_curve:
        fpr, tpr, _ = roc_curve(df[label], df[_prob])
        plot_curve(fpr, tpr)
    print("\nROC_AUC score {:0.3f}".format(roc_auc_score(df[label], df[_prob])))


middle_east = ["Jordan", "Oman", "Qatar", "Saudi Arabia", "United Arab Emirates"]
bind_human = {
    "KT444582": "Bat WIV16",
    "MN996532": "RaTG13",
    "KF367457": "Bat WIV1",
    "MZ190137": "Khosta-1",
    "MZ190138": "Khosta-2",
    "KC881005": "Bat RsSHC014",
    "KY417151": "Rs7327",
    "KC881006": "Rs3367",
    "KY417152": "Rs9401",
    "KY417146": "Rs4231",
    "KY417150": "Rs4874",
    "KY417144": "Rs4084",
    "ASL68941": "HKU25",
    "ASL68953": "HKU25",
    "KF569996": "Bat LYRa11",
    "MK211376": "YN2018B",
    "AWH65877": "HKU4",
    "ABN10848": "HKU4",
    "ABN10857": "HKU4",
    "ABN10866": "HKU4",
    "AWH65888": "HKU4",
    "AWH65899": "HKU4",
    "QHA24678": "HKU4",
    "YP_001039953": "HKU4",
    "AY304486": "SZ3",
    "AY686863": "A022",
    "AY274119": "Tor2",
    "YP_009724390": "Wuhan",
    "ACT11030": "HEC",
    "ACJ35486": "HEC",
    "AHN64774": "HKU23",
    "AHN64783": "HKU23",
    "QEY10625": "HKU23",
    "QEY10641": "HKU23",
    "QEY10649": "HKU23",
    "QEY10657": "HKU23",
    "QEY10633": "HKU23",
    "ALA50080": "HKU23",
    "JX163927": "SARS1",
    "JX163926": "SARS1",
    "AY545919": "SARS1",
}


# ab = pd.read_csv("./data/alpha_beta.csv")
# camel_bind_human = ab[(ab.Host_agg == 'Camels') & \
#                       (ab.Species_agg == 'MERSr') & \
#                       (ab.Country.isin(middle_east))]
