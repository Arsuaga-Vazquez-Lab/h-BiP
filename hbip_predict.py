import pickle
import sys
import pandas as pd

# HIP libraries
import preprocess as pp
import hip_reproduce as mn


def read_seqs(_file_path):
    """ Read text file (a fasta file) with sequences and description

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
        if not line.startswith('>'):
            read += line.rstrip('\n')
        if line.startswith('>') or line == last_line:
            if line != last_line:
                clean_meta = line[1:].strip('\n')
                meta.append(clean_meta)
            if read != "":
                seq.append(read)
                read = ""
    df = pd.DataFrame({'Description': meta, 'Sequence': seq})
    return df


def predict_new(name, sequence):
    """
    Use the model to predict a list of new sequences

    :param name: str
        Will read the model from ./models/name
    :param sequence: Series
        Sequence list
    :return: float
        Score for the sequence
    """
    one_x = pp.single_embedding(name, sequence)
    directory = './models/' + name
    LR_path = directory + '/LR85.pkl'
    with open(LR_path, 'rb') as f:
        LR = pickle.load(f)
    scale_path = directory + '/scale.pkl'
    with open(scale_path, 'rb') as f:
        scale = pickle.load(f)
    one_x_norm = scale.transform(one_x)
    return mn.predict(LR, one_x_norm)


if __name__ == '__main__':
    """
    Usage CLI: python3 hbip_predict.py path_to_fasta_file model_name
    Usage CLI: python3 hbip_predict.py ./data/virus_predict.fasta alpha_beta
    If no parameters are provided it will compute the scores for SARS-CoV-2 for alpha_beta
    """

    model_name = 'alpha_beta'
    sequences_file_path = './data/sars_cov2.fasta'
    if len(sys.argv) > 1:
        sequences_file_path = sys.argv[1]
        if len(sys.argv) == 3:
            model_name = sys.argv[2]
    df_from_fasta = read_seqs(sequences_file_path)
    sequences_list = df_from_fasta.iloc[:, 1]

    scores = predict_new(model_name, sequences_list)
    df_scores = pd.DataFrame({'Description': df_from_fasta.iloc[:, 0], 'Scores': scores})
    print('\nScores from {} model:'.format(model_name))
    print(df_scores)


