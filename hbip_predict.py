import sys
import pandas as pd

from hbip.predict import read_seqs, predict_new

import warnings
warnings.filterwarnings("ignore")

def main(sequences_file_path, model_name):
    df_from_fasta = read_seqs(sequences_file_path)
    sequences_list = df_from_fasta.iloc[:, 1]

    scores = predict_new(model_name, sequences_list, models_dir="./models/")
    df_scores = pd.DataFrame(
        {"Description": df_from_fasta.iloc[:, 0], "Scores": scores}
    )
    print("\nScores from {} model:".format(model_name))
    print(df_scores)


if __name__ == "__main__":
    """
    Usage CLI: python3 hbip_predict.py path_to_fasta_file model_name
    Usage CLI: python3 hbip_predict.py ./data/virus_predict.fasta alpha_beta
    If no parameters are provided it will compute the scores for SARS-CoV-2 for alpha_beta
    """

    model_name = "alpha_beta"
    sequences_file_path = "./data/sars_cov2.fasta"

    if len(sys.argv) > 1:
        sequences_file_path = sys.argv[1]
    if len(sys.argv) == 3:
        model_name = sys.argv[2]
    main(sequences_file_path, model_name)
