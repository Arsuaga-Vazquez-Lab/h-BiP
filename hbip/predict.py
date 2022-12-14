import pickle

import pandas as pd

from hbip import preprocess as pp


def read_seqs(_file_path):
    """Read text file (a fasta file) with sequences and description

    Parameters
    ----------
    _file_path: str
        Path to fasta file

    Returns
    -------
    df : dataframe
        Dataframe with columns: Description, Sequence
    """

    print("\nReading Sequences...\n")
    seq = []
    meta = []
    read = ""

    with open(_file_path) as f:
        lines = f.readlines()
    last_line = lines[-1]
    for line in lines:
        if not line.startswith(">"):
            read += line.rstrip("\n")
        if line.startswith(">") or line == last_line:
            if line != last_line:
                clean_meta = line[1:].strip("\n")
                meta.append(clean_meta)
            if read != "":
                seq.append(read)
                read = ""
    df = pd.DataFrame({"Description": meta, "Sequence": seq})
    return df


def predict_new(name, sequence, models_dir):
    """
    Use the model to predict a list of new sequences

    :param name: str
        Will read the model from ./models/name
    :param sequence: Series
        Sequence list
    :param models_dir: str
        Path to models directory
    :return: float
        Score for the sequence
    """
    one_x = pp.single_embedding(name, sequence)
    directory = models_dir + name
    LR_path = directory + "/LR85.pkl"
    with open(LR_path, "rb") as f:
        LR = pickle.load(f)
    scale_path = directory + "/scale.pkl"
    with open(scale_path, "rb") as f:
        scale = pickle.load(f)
    one_x_norm = scale.transform(one_x)
    return predict(LR, one_x_norm)


def predict(model, x_scaled):
    """
    Load the logistic regression model and predict the score for a normalized sequence

    :param model: Object
        Logistic regression model
    :param x_scaled:
        Embedding previously scaled by the trained scaler
    :return: Series
        Probability scores from the logistic regression model
    """
    # predict from standardized embedding
    prob = model.predict_proba(x_scaled)
    human_prob = [sample[1] for sample in prob]
    return human_prob
