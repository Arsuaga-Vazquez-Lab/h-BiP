from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, classification_report
import yaml
import pandas as pd
from copy import deepcopy

# sPredCov libraries
from preprocess import *
from analysis import *


def continuous_trimer(seqs):
    """Creates all possible trimers for a sequence padding the beginning and end
    (e.g. 'ABCD' = '--A, -AB, ABC, BCD, CD-, D--'
     seqs: Series"""
    def continuous_trimer_one(seq):
        ct = ['--' + seq[0], '--' + seq[:2]]
        last = len(seq) - 1
        for i in range(last - 1):
            ct.append(seq[i:i + 3])
        ct.append(seq[-2:] + '-')
        ct.append(seq[-1] + '--')
        return ct

    trigram = []
    for s in seqs:
    # for i in range(df.shape[0]):
        cont3_seq = continuous_trimer_one(s)
        trigram.append(cont3_seq)
    return trigram


def embedding_continuous(train_df, test_df):
    """Creates trigrams from the training set and feeds them to a skip-gram model with
    negative sampling (word2vec) to create the sequence embedding. The test set uses
    the resulting model to create an embedding. Unknown words are replaced with a mean vector

    Parameters
    ----------
    train_df, test_df : DataFrame
        Dataframes for the train and test set with columns [ 'Sequence', 'Accession', 'Human']

    Returns
    --------
    train_triple, test_triple: np.ndarray
        Embedding for the Sequence
    train_target, test_target : DataFrame
        Details about the target and sequence after data augmentation: ['Accession', 'Slide', 'Human']
    """
    df_3w_train = continuous_trimer(train_df.Sequence)
    train_target = train_df.loc[:, ['Accession', 'Human']].copy()
    train_target['Source'] = ['Train'] * (train_target.shape[0])
    df_3w_test = continuous_trimer(test_df.Sequence)
    test_target = test_df.loc[:, ['Accession', 'Human']].copy()
    test_target['Source'] = ['Test'] * (test_target.shape[0])
    mat_vect, vocab, mean_vec = create_word_vectors(df_3w_train)
    print('Training set:')
    train_embedding = create_seq_embedding(df_3w_train, mat_vect, vocab, mean_vec)
    print('\nTest set:')
    test_embedding = create_seq_embedding(df_3w_test, mat_vect, vocab, mean_vec)
    return train_embedding, train_target, test_embedding, test_target


wuhan = {'Accession': 'YP_009724390',
         'Sequence': 'MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAI\
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
         VAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT',
         'Species': 'Severe acute respiratory syndrome coronavirus 2',
         'Virus': 'Wuhan-Hu-1',
         'Note': np.nan,
         'Host': 'Homo sapiens',
         'Host_agg': 'Homo sapiens',
         'Species_agg': 'SARS-CoV-2',
         'Human': 1}


def predict(model, x):
    prob = model.predict_proba(StandardScaler().fit_transform(x))
    human_prob = [sample[1] for sample in prob]
    return human_prob


def aggregate_triplicates(df, agg='max'):
    # Remove triplicates aggregating by accession (max or average)
    keep = ['Accession', 'Source', 'Human', 'Human_prob']  # Remove Slide to avoid taking the max
    if agg == 'average':
        df_unique = df[keep].groupby(['Accession'], as_index=False).mean()
    else:
        print('Using the max to aggregate')
        df_unique = df[keep].groupby(['Accession'], as_index=False).max()

    df_unique['Predicted'] = 0
    df_unique.loc[df_unique.Human_prob >= 0.5, 'Predicted'] = 1
    return df_unique


def metrics(df):
    print(classification_report(df.Human, df.Predicted))
    print('Note: Sensitivity=Recall(1), Specificity=Recall(0), FPR=1-Recall(0)')

    print('\nROC_AUC score {:0.3f}'.format(roc_auc_score(df.Human, df.Human_prob)))


def main_continuous(data_name, remove_sars2=False, sample_sars2=False):
    config_path = './data/' + data_name + '_config.yml'
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)

    data_path = config['path']
    data = pd.read_csv(data_path)
    print('\nProcessing data:', config['name'])
    print('Shape of initial input data is:', data.shape)

    # Remove undesired viruses
    undesired_viruses = len(data[data.Accession.isin(config['acc_del'])])
    print('\nDataset has {} viruses to be removed'.format(undesired_viruses))
    if undesired_viruses != 0:
        data = data[~data.Accession.isin(config['acc_del'])]
        print('Removed {} viruses'.format(undesired_viruses))
        print('Shape of data is now: {}'.format(data.shape))

    # Remove SARS-CoV-2
    if remove_sars2:
        print('Removing SARS-CoV-2 from dataset')
        data = data[data.Species_agg != 'SARS-CoV-2']
        print('Shape of data is now: {}'.format(data.shape))

    if sample_sars2:
        print('Keeping 20% of original SARS-CoV-2 samples')
        small_s2 = data[data.Species_agg == 'SARS-CoV-2'].sample(frac=.2, random_state=10)
        data = data[data.Species_agg != 'SARS-CoV-2']
        data = pd.concat([data, small_s2], axis=0)
        print('Shape of data is now: {}'.format(data.shape))

    # Remove viruses that we want at training set
    add_train = data[data.Accession.isin(config['acc_train'])]
    print('\nThere are {} viruses specifically selected to be at training data'.format(add_train.shape[0]))
    data = data[~data.Accession.isin(config['acc_train'])]

    # Split data into train and test set
    print('\nSplitting the dataset...')
    train, test = train_test(data, strat=config['strat'], test_frac=config['test_fraction'])

    # Add viruses required at training set
    train = pd.concat([train, add_train], axis=0)
    print('\nTraining set after adding selected viruses: ', train.shape[0])

    # Include Wuhan in test dataset
    if config['wuhan_at_test']:
        print('\nAdded Wuhan virus to test dataset')
        # datasets missing columns such as Species_agg will have nans in the rest of records
        test = test.append(wuhan, ignore_index=True)
        if config['rbd']:
            print('\nExtracting Wuhan RBD')
            test.loc[test.Accession == 'YP_009724390', 'Sequence'] = wuhan['Sequence'][333:526]
    print('Final full dataset: {} \n'.format(train.shape[0] + test.shape[0]))

    #  target is a df with columns: [Accession, Slide, Human]
    X85, target85, X_test, target_test = embedding_continuous(train, test)

    # Normalizing
    X85_norm = StandardScaler().fit_transform(X85)
    Xtest_norm = StandardScaler().fit_transform(X_test)

    # Full dataset (with embeddings)
    Xall = deepcopy(X85)
    Xall = np.concatenate((Xall, X_test), axis=0)
    df_all = pd.concat([target85, target_test], axis=0)

    # Modeling
    LR85 = LogisticRegression(C=1.5, penalty='l2', random_state=0, solver='lbfgs', max_iter=400)
    LR85 = LR85.fit(X85_norm, target85.Human)
    print('\nPREPROCESSING APPROACH: continuous')
    print('\nScore train: {:0.2f}'.format(LR85.score(X85_norm, target85.Human)))

    # Predict for all data and remove triplicates
    df_all['Human_prob'] = predict(LR85, Xall)
    df_all['Predicted'] = 0
    df_all.loc[df_all.Human_prob >= 0.5, 'Predicted'] = 1
    # df_all = aggregate_triplicates(df_all, agg='max')

    print('\nMetrics for TEST set:')
    print('----------------------------------------')
    print('\nFirst 5 sequences and their scores...')
    print(df_all[df_all.Source == 'Test'].head())
    print('')
    metrics(df_all[df_all.Source == 'Test'])

    print('\nMetrics for FULL set:')
    print('----------------------------------------')
    metrics(df_all)

    print('Prediction for wuhan is: {}'.format(df_all.loc[df_all.Accession == 'YP_009724390'].Human_prob))

    print(all_important_scores(df_all))
    return data, df_all


if __name__ == '__main__':
    dataset, scores = main_continuous('ab_same_sars1', remove_sars2=False, sample_sars2=True)
    # dataset, scores = main_continuous('sars1_195')
    # dataset, scores = main_continuous('rbd_189')
    scores.to_csv('./data/ab_same_sars1_balancedS2_scores_cont.csv', index=False)
    # scores.to_csv('./data/rbd_189_scores_cont.csv', index=False)
