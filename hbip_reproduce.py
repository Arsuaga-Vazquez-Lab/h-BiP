from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os

import yaml
from copy import deepcopy

from hbip.predict import predict
from hbip.utils.data_prep import aggregate_triplicates, metrics, prepare_data
from hbip.preprocess import *


def modeling(config, train, test, continuous_trimer=True, save_model=False):
    """
    Generate the model for train and test using the arguments from the config file.
    Produce the scores for the human infection potential

    :param config: Object
        Object containing the information from the config yaml file
    :param train: DataFrame
        Train dataset
    :param test: DataFrame
        Test dataset
    :param continuous_trimer: Boolean
        If True will create traditional trimers from the sequence.
        False will create all three possible trimers (triplicates)
    :param save_model: Boolean
        If True, will save the model in ./models/config['name'].
        Model consist of 5 parts: scale, LR85,  mat_vect, mean_vect, vocab
    :return Xall: ndarray
        Embedding (before normalizing)
    :return df_meta: DataFrame
        Dataframe with the h-BiP scores and meta data
    """
    #  target is a df with columns: [Accession, Slide, Human]
    if continuous_trimer:
        print("\nPREPROCESSING APPROACH: classic trimer")
        # X85, target85, X_test, target_test = embedding_continuous(train, test, label=config['target'])

        X85, target85, X_test, target_test = embedding_continuous(
            config["name"], train, test, config["target"], save_embed=save_model
        )
    else:
        print("\nPREPROCESSING APPROACH: triplicates")
        X85, target85, X_test, target_test = embedding(
            train, test, label=config["target"]
        )

    # Normalizing
    scale = StandardScaler().fit(X85)
    X85_norm = scale.transform(X85)
    # X85_norm = StandardScaler().fit_transform(X85)
    # Xtest_norm = StandardScaler().fit_transform(X_test)

    # Full dataset (with embeddings)
    Xall = deepcopy(X85)
    Xall = np.concatenate((Xall, X_test), axis=0)
    df_meta = pd.concat([target85, target_test], axis=0)

    # Modeling
    LR85 = LogisticRegression(
        C=1.5, penalty="l2", random_state=0, solver="lbfgs", max_iter=400
    )
    LR85 = LR85.fit(X85_norm, target85[config["target"]])
    if save_model:
        save_model_here(config["name"], scale, "scale")
        save_model_here(config["name"], LR85, "LR85")
    print(
        "\nScore train: {:0.2f}".format(
            LR85.score(X85_norm, target85[config["target"]])
        )
    )

    # Predict for all data and remove triplicates
    prob = config["target"] + "_prob"
    # Xall_norm = StandardScaler().fit_transform(Xall)  # old
    Xall_norm = scale.transform(Xall)  # old
    df_meta[prob] = predict(LR85, Xall_norm)
    # df_meta[prob] = predict(LR85, Xall)   # old
    pred = config["target"] + "_predicted"
    df_meta[pred] = 0
    df_meta.loc[df_meta[prob] >= 0.5, pred] = 1
    if not continuous_trimer:
        df_meta = aggregate_triplicates(df_meta, agg="max", label=config["target"])
        print("\nMetrics for TEST set (no triplicates):")
    else:
        print("\nMetrics for TEST set:")
    print("----------------------------------------")
    metrics(
        df_meta[df_meta.Source == "Test"], label=config["target"], plot_roc_curve=False
    )

    if not continuous_trimer:
        print("\nMetrics for FULL set (no triplicates):")
    else:
        print("\nMetrics for FULL set:")
    print("----------------------------------------")
    metrics(df_meta, label=config["target"])

    print(
        "\nPrediction for wuhan is: {} \n".format(
            df_meta.loc[df_meta.Accession == "YP_009724390", prob]
        )
    )
    return Xall, df_meta


def plot_embedding(
    output_name, x, meta, resolution="half", by="Binds", horizontal=True, save=False
):
    """
    2D scatterplot of the embedding after applying TSNE color coded using 'by'

    :param output_name: str
        if save=True, the plot will be save as ./outputs/output_name_tsne.pdf'
    :param x: ndarray
        Embedding
    :param meta: DataFrame
        Scores dataframe with the target feature
    :param resolution: {'display', 'full', 'half}
        Different sizes for the plot. Best for display, full letter size, half letter size
    :param by: str
        Plot the embedding with different colors according to this column (hue)
    :param horizontal: Boolean
        if True, horizontal legend
    :param save: Boolean
        if True, save plot as './outputs/' + output_name + '_tsne.pdf'
    :return: None
    """
    # resolution: display, letter, half letter
    size = {"display": (8, 6), "full": (6, 4), "half": (3.3, 2.5)}
    tsne_all = TSNE(
        n_components=2, verbose=1, perplexity=40, n_iter=300, random_state=0
    ).fit_transform(x)
    tsne_all_df = pd.DataFrame(
        {
            "tsne1": tsne_all[:, 0],
            "tsne2": tsne_all[:, 1],
        }
    )
    print(tsne_all_df.shape)
    print(meta.shape)
    tsne_all_df = pd.concat(
        [tsne_all_df.reset_index(drop=True), meta.reset_index(drop=True)], axis=1
    )

    plt.figure(figsize=(size[resolution]))  # for display
    plt.rc("font", size=6.5)

    sns.scatterplot(
        x="tsne1",
        y="tsne2",
        hue=by,
        style=by,
        data=tsne_all_df,
        s=7,
        palette=["dimgray", "indianred"],
        markers=["o", "P"],
    )  # markers = '.'
    plt.ylabel("TSNE2", labelpad=0)
    plt.xlabel("TSNE1")
    if horizontal:
        legend = plt.legend(bbox_to_anchor=(1, 1.12), ncol=2, labels=["Yes", "No"])
    else:
        legend = plt.legend(bbox_to_anchor=(0.65, 0.12), title=by, labels=["Yes", "No"])
    legend.legendHandles[1].set_sizes([2])
    if save:
        file_name = "./outputs/" + output_name + "_tsne.pdf"
        plt.savefig(file_name, bbox_inches="tight")  # or foo.png
    plt.show()


def main(
    config_path,
    exclude_sars2=False,
    downsample_sars2=False,
    classic_trimer=True,
    plot=False,
    save_model=False,
):
    """
    This function trains a model using the arguments defined by the user at the config file.
    It can plot the embedding after applying dimensionality reduction from TSNE.
    :param config_path: str
        Path to config file related to the data
    :param exclude_sars2: Boolean
        If True will exclude all SARS2 viruses (Species_agg != 'SARS-CoV-2')
    :param downsample_sars2: Boolean
        If True will take only a fraction of SARS-CoV-2 viruses (default=20%)
    :param classic_trimer: Boolean
        If True will create traditional trimers from the sequence.
        False will create all three possible trimers (triplicates)
    :param plot: Boolean
        If True, generate 2D scatterplot after applying TSNE to the embedding
    :param save_model: Boolean
        If True, save the model (binary)
    :return embedding_all: ndarray
        Embedding for full data
    :return final_scores: DataFrame
        DataFrame with final scores stores in a column with name ending in '_prob'
    """

    def save_train_test():
        train_output = config.get("train_output", "")
        if train_output == "":
            train_output = "./data/" + config["name"] + "_train.csv"
        test_output = config.get("test_output", "")
        if test_output == "":
            test_output = "./data/" + config["name"] + "_test.csv"
        df_train.to_csv(train_output, index=False)
        print("Saving train set to: ", train_output)
        df_test.to_csv(test_output, index=False)
        print("Saving test set to: ", test_output)

        # Uncomment if you want to save full dataset (after removals)
        # full_output = './data/' + config['name'] + '_cleaned.csv'
        # df_full.to_csv(full_output, index=False)
        # print('Saving full set to: ', full_output)

    # Create outputs directory
    if not os.path.exists("./outputs"):
        os.makedirs("./outputs")

    with open(config_path, "r") as file:
        config = yaml.safe_load(file)
    target = config["target"]
    name = config["name"]
    print("\nGenerating model for :", name)
    print("Description: ", config["description"], "\n")

    if config["fixed_train_test"]:
        print("Reading user provided train and test sets...")
        df_train = pd.read_csv(config["train"])
        df_test = pd.read_csv(config["test"])

        if len(config["acc_del"]) != 0:
            print("Train length: ", df_train.shape[0])
            to_remove = df_train.loc[
                df_train.Accession.isin(config["acc_del"]), ["Accession", "Virus"]
            ]
            print(
                "\n Removing the following {} viruses from the training set:".format(
                    len(to_remove)
                )
            )
            print(to_remove)
            df_train = df_train[~df_train.Accession.isin(config["acc_del"])]
            print("Final train length: ", df_train.shape[0])

    else:
        print("Generating train and test sets...")
        df_full, df_train, df_test = prepare_data(
            config, exclude_sars2, downsample_sars2
        )
        save_train_test()

    embedding_all = None
    final_scores = None
    embedding_all, final_scores = modeling(
        config, df_train, df_test, classic_trimer, save_model
    )

    if plot:
        plot_embedding(
            name, embedding_all, final_scores, resolution="half", by=target, save=True
        )

    if config["save_scores"]:
        scores_output = "./outputs/" + config["name"] + "_scores.csv"
        final_scores.to_csv(scores_output, index=False)

    return embedding_all, final_scores


if __name__ == "__main__":
    """
    Usage: python3 hbip_reproduce.py config_file_path
    If no parameters are provided, it will use the full dataset for alpha and beta coronaviruses
    """

    if len(sys.argv) > 1:
        config_file_path = sys.argv[1]
    else:
        config_file_path = "./data/alpha_beta_config.yml"

    embed, scores = main(config_file_path, plot=False, save_model=False)
