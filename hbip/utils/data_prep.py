import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.metrics import classification_report, roc_curve, roc_auc_score

from hbip.preprocess import train_test

wuhan = {
    "Accession": "YP_009724390",
    "Sequence": "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAI\
HVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKN\
NKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVD\
LPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTV\
EKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLN\
DLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDI\
STEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTG\
TGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQ\
LTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYS\
NNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQI\
YKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDE\
MIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKL\
QDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATK\
MSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRN\
FYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNE\
VAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT",
    "Species": "Severe acute respiratory syndrome coronavirus 2",
    "Virus": "Wuhan-Hu-1",
    "Note": np.nan,
    "Host": "Homo sapiens",
    "Host_agg": "Homo sapiens",
    "Species_agg": "SARS-CoV-2",
    "Human": 1,
    "Binds": 1,
}


def aggregate_triplicates(df, agg="max", label="Human"):
    # Remove triplicates aggregating by accession (max or average)
    _prob = label + "_prob"
    _pred = label + "_predicted"
    keep = ["Accession", "Source", label, _prob]  # Remove Slide to avoid taking the max
    if agg == "average":
        df_unique = df[keep].groupby(["Accession"], as_index=False).mean()
    else:
        print("Using the max to aggregate")
        df_unique = df[keep].groupby(["Accession"], as_index=False).max()

    df_unique[_pred] = 0
    df_unique.loc[df_unique[_prob] >= 0.5, _pred] = 1
    return df_unique


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


def prepare_data(config, remove_sars2=False, sample_sars2=False):
    data_path = config["path"]
    data = pd.read_csv(data_path)
    print("\nProcessing data:", config["name"])
    print("Shape of initial input data is:", data.shape)

    # Remove undesired viruses
    undesired_viruses = len(data[data.Accession.isin(config["acc_del"])])
    print("\nDataset has {} undesired viruses".format(undesired_viruses))
    print(data.loc[data.Accession.isin(config["acc_del"]), ["Accession", "Virus"]])
    if undesired_viruses != 0:
        data = data[~data.Accession.isin(config["acc_del"])]
        print("Removed {} viruses".format(undesired_viruses))
        print("Shape of data is now: {}".format(data.shape))

    # Remove SARS-CoV-2
    if remove_sars2:
        print("Removing SARS-CoV-2 from dataset")
        data = data[data.Species_agg != "SARS-CoV-2"]
        print("Shape of data is now: {}".format(data.shape))

    if sample_sars2:
        print("Keeping 20% of original SARS-CoV-2 samples")
        small_s2 = data[data.Species_agg == "SARS-CoV-2"].sample(
            frac=0.2, random_state=10
        )
        data = data[data.Species_agg != "SARS-CoV-2"]
        data = pd.concat([data, small_s2], axis=0)
        print("Shape of data is now: {}".format(data.shape))

    # Remove viruses that we want at training set
    add_train = data[data.Accession.isin(config["acc_train"])]
    print(
        "\nThere are {} viruses specifically selected to be at training data".format(
            add_train.shape[0]
        )
    )
    for_split = data[~data.Accession.isin(config["acc_train"])]

    # Split data into train and test set
    print("\nSplitting the dataset...")
    train, test = train_test(
        for_split, strat=config["strat"], test_frac=config["test_fraction"]
    )

    # Add viruses required at training set
    train = pd.concat([train, add_train], axis=0)
    print("\nTraining set after adding selected viruses: ", train.shape[0])

    # Include Wuhan in test dataset
    if config["wuhan_at_test"]:
        print("\nAdded Wuhan virus to test dataset")
        # datasets missing columns such as Species_agg will have nans in the rest of records
        test = test.append(wuhan, ignore_index=True)
        if config["rbd"]:
            test.loc[test.Accession == "YP_009724390", "Sequence"] = wuhan["Sequence"][
                333:526
            ]
    print("Final full dataset: {} \n".format(train.shape[0] + test.shape[0]))
    return data, train, test
