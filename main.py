from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, classification_report, roc_curve
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import seaborn as sns

import yaml
import pandas as pd
from copy import deepcopy
import numpy as np

# sPredCov libraries
from preprocess import *

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
         'Human': 1,
         'Binds': 1}


def predict(model, x):
    prob = model.predict_proba(StandardScaler().fit_transform(x))
    human_prob = [sample[1] for sample in prob]
    return human_prob


def aggregate_triplicates(df, agg='max', label='Human'):
    # Remove triplicates aggregating by accession (max or average)
    _prob = label + '_prob'
    _pred = label + '_predicted'
    keep = ['Accession', 'Source', label, _prob]  # Remove Slide to avoid taking the max
    if agg == 'average':
        df_unique = df[keep].groupby(['Accession'], as_index=False).mean()
    else:
        print('Using the max to aggregate')
        df_unique = df[keep].groupby(['Accession'], as_index=False).max()

    df_unique[_pred] = 0
    df_unique.loc[df_unique[_prob] >= 0.5, _pred] = 1
    return df_unique


def metrics(df, label, plot_roc_curve=False):
    def plot_curve(fp, tp):
        plt.figure()
        plt.plot(fp, tp,
            color="darkorange", lw=2)
        plt.plot([0, 1], [0, 1], color="navy", lw=2, linestyle="--")
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate")
        plt.show()

    _prob = label + '_prob'
    _pred = label + '_predicted'
    print(classification_report(df[label], df[_pred]))
    print('Note: Sensitivity=Recall(1), Specificity=Recall(0), FPR=1-Recall(0)')

    if plot_roc_curve:
        fpr, tpr, _ = roc_curve(df[label], df[_prob])
        plot_curve(fpr, tpr)
    print('\nROC_AUC score {:0.3f}'.format(roc_auc_score(df[label], df[_prob])))


def prepare_data(data_name, remove_sars2=False, sample_sars2=False):
    config_path = './data/' + data_name + '_config.yml'
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)

    data_path = config['path']
    data = pd.read_csv(data_path)
    print('\nProcessing data:', config['name'])
    print('Shape of initial input data is:', data.shape)

    # Remove undesired viruses
    undesired_viruses = len(data[data.Accession.isin(config['acc_del'])])
    print('\nDataset has {} undesired viruses'.format(undesired_viruses))
    print(data.loc[data.Accession.isin(config['acc_del']), ['Accession', 'Virus']])
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
            test.loc[test.Accession == 'YP_009724390', 'Sequence'] = wuhan['Sequence'][333:526]
    print('Final full dataset: {} \n'.format(train.shape[0] + test.shape[0]))
    return data, train, test


def modeling(train, test, continuous_trimer, label):
    #  target is a df with columns: [Accession, Slide, Human]
    if continuous_trimer:
        print('\nPREPROCESSING APPROACH: continuous')
        X85, target85, X_test, target_test = embedding_continuous(train, test, label=label)
    else:
        print('\nPREPROCESSING APPROACH: triplicates')
        X85, target85, X_test, target_test = embedding(train, test, label=label)

    # Normalizing
    X85_norm = StandardScaler().fit_transform(X85)
    # Xtest_norm = StandardScaler().fit_transform(X_test)

    # Full dataset (with embeddings)
    Xall = deepcopy(X85)
    Xall = np.concatenate((Xall, X_test), axis=0)
    df_meta = pd.concat([target85, target_test], axis=0)

    # Modeling
    LR85 = LogisticRegression(C=1.5, penalty='l2', random_state=0, solver='lbfgs', max_iter=400)
    LR85 = LR85.fit(X85_norm, target85[label])
    print('\nScore train: {:0.2f}'.format(LR85.score(X85_norm, target85[label])))

    # Predict for all data and remove triplicates
    prob = label + '_prob'
    df_meta[prob] = predict(LR85, Xall)
    pred = label + '_predicted'
    df_meta[pred] = 0
    df_meta.loc[df_meta[prob] >= 0.5, pred] = 1
    if not continuous_trimer:
        df_meta = aggregate_triplicates(df_meta, agg='max', label=label)
        print('\nMetrics for TEST set (no triplicates):')
    else:
        print('\nMetrics for TEST set:')
    print('----------------------------------------')
    metrics(df_meta[df_meta.Source == 'Test'], label=label, plot_roc_curve=False)

    if not continuous_trimer:
        print('\nMetrics for FULL set (no triplicates):')
    else:
        print('\nMetrics for FULL set:')
    print('----------------------------------------')
    metrics(df_meta, label=label)

    print('Prediction for wuhan is: {}'.format(df_meta.loc[df_meta.Accession == 'YP_009724390', prob]))
    return Xall, df_meta


def plot_embedding(output_name, x, meta, resolution='half', by='Binds', horizontal=True, save=False):
    # resolution: display, letter, half letter
    size = {'display': (8,6), 'full': (6,4), 'half': (3.3, 2.5)}
    tsne_all = TSNE(n_components=2, verbose=1, perplexity=40,
                    n_iter=300, random_state=0).fit_transform(x)
    tsne_all_df = pd.DataFrame({'tsne1': tsne_all[:, 0], 'tsne2': tsne_all[:, 1], })
    print(tsne_all_df.shape)
    print(meta.shape)
    tsne_all_df = pd.concat([tsne_all_df.reset_index(drop=True), meta.reset_index(drop=True)], axis=1)

    plt.figure(figsize=(size[resolution]))  # for display
    plt.rc('font', size=6.5)

    sns.scatterplot(x='tsne1', y='tsne2', hue=by, style=by, data=tsne_all_df, s=7,
                    palette=['dimgray', 'indianred'], markers=['o', 'P'])   # markers = '.'
    plt.ylabel('TSNE2', labelpad=0)
    plt.xlabel('TSNE1')
    if horizontal:
        legend = plt.legend(bbox_to_anchor=(1, 1.12), ncol=2,  labels=['Yes', 'No'])
    else:
        legend = plt.legend(bbox_to_anchor=(.65, .12), title=by, labels=['Yes', 'No'])
    legend.legendHandles[1].set_sizes([2])
    if save:
        file_name = './data/' + output_name + '_tsne.pdf'
        plt.savefig(file_name, bbox_inches='tight')   # or foo.png
    plt.show()


def main(name, data_from_file=True, exclude_sars2=False, downsample_sars2=False,
         classic_trimer=True, target='Binds', plot=False):
    """
    :param name: str
        Database name
    :param data_from_file: bool
        Determine if the dataset will be generated or read from a file
    :param exclude_sars2: bool
        It removes all SARS-CoV-2 viruses from the dataset
    :param downsample_sars2: bool
        It randomly selects 20% of SARS-CoV-2 viruses
    :param classic_trimer: bool
        Approach to generate the trimers: traditional (continuous) or triplicates
    :param target: str
        The name of the target label for classification
    :return: none
    """
    config_path = './data/' + name + '_config.yml'
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    print('Generating model for :', config['name'])

    if data_from_file:
        print('Reading provided train and test sets...')
        df_train = pd.read_csv(config['train'])
        df_test = pd.read_csv(config['test'])
    else:
        print('Generating train and test sets...')
        df_full, df_train, df_test = prepare_data(name,
                                                  exclude_sars2, downsample_sars2)
        # Uncomment if you want to save the files
        # df_train.to_csv(config['train'])
        # df_test.to_csv(config['test'])

    embedding_all = None
    final_scores = None
    for exp in config['experiments']:
        print('\n=================')
        print(exp, ': ', config[exp]['description'])
        exp_train = df_train[~df_train.Accession.isin(config[exp]['remove'])]
        embedding_all, final_scores = modeling(exp_train, df_test, classic_trimer, target)

    if plot:
        if len(config['experiments']) != 0:
            # Will only return the last experiment
            plot_embedding(name, embedding_all, final_scores, resolution='half', by=target, save=True)

    return embedding_all, final_scores


if __name__ == '__main__':

    # dataset, scores = main('ab_same_sars1_binds', remove_sars2=True, sample_sars2=False,
    #                        trimer_continuous=True, label='Binds')  # DONE

    # nos2
    # df_full, df_train, df_test = prepare_data('ab_same_sars1_binds',
    #                                           remove_sars2=True, sample_sars2=False)

    # main('nos2', data_from_file=True)
    # main('rbd_189_binds', data_from_file=False)




    # Save training and test set
    # df_full.to_csv('./data/ab_same_sars1_del_6_muts_nos2_binds.csv', index=False)
    # df_train.to_csv('./data/ab_same_sars1_del_6_muts_nos2_binds_train.csv', index=False)
    # df_test.to_csv('./data/ab_same_sars1_del_6_muts_nos2_binds_test.csv', index=False)


    # dataset, scores = main('sars1_195', remove_sars2=False, sample_sars2=False,
    #                       trimer_continuous=False, label='Human')
    # scores.to_csv('./data/sars1_195_scores_cont.csv', index=False)

    # dataset, scores = main('rbd_189', remove_sars2=False, sample_sars2=False,
    #                       trimer_continuous=False, label='Human')  # DONE
    # scores.to_csv('./data/rbd_189_scores.csv', index=False)

    # dataset, scores = main('rbd_189', remove_sars2=False, sample_sars2=False,
    #                       trimer_continuous=True, label='Human')  # DONE
    # scores.to_csv('./data/rbd_189_scores_cont.csv', index=False)

    # dataset, scores = main('ab_same_sars1', remove_sars2=False, sample_sars2=False,
    #                       trimer_continuous=False, label='Human')
    # scores.to_csv('./data/ab_same_sars1_scores.csv', index=False)

    # dataset, scores = main('ab_same_sars1', remove_sars2=True, sample_sars2=False,
    #                       trimer_continuous=False, label='Human')
    # scores.to_csv('./data/ab_same_sars1_nos2_scores.csv', index=False)

    # dataset, scores = main('ab_same_sars1', remove_sars2=True, sample_sars2=False,
    #                       trimer_continuous=True, label='Human')
    # scores.to_csv('./data/ab_same_sars1_nos2_scores_cont.csv', index=False)

    # dataset, scores = main('ab_same_sars1', remove_sars2=False, sample_sars2=True,
    #                       trimer_continuous=False, label='Human')
    # scores.to_csv('./data/ab_same_sars1_balancedS2_scores.csv', index=False)

    # =================
    # RUN BINDING TARGET, TRIPLICATES (OLD)

    # dataset, scores = main('ab_same_sars1_binds', remove_sars2=False, sample_sars2=False,
    #                       trimer_continuous=False, label='Binds')
    # scores.to_csv('./data/ab_same_sars1_binds_scores.csv', index=False)

    # dataset, scores = main('ab_same_sars1_binds', remove_sars2=True, sample_sars2=False,
    #                       trimer_continuous=False, label='Binds')   # DONE
    # scores.to_csv('./data/ab_same_sars1_nos2_binds_scores.csv', index=False)

    # dataset, scores = main('ab_same_sars1_binds', remove_sars2=True, sample_sars2=False,
    #                       trimer_continuous=False, label='Binds')   # DONE
    # scores.to_csv('./data/ab_same_sars1_nos2_binds_del_muts_scores.csv', index=False)

    # =================
    # RUN BINDING TARGET, CONTINUOUS

    # dataset, scores = main('ab_same_sars1_binds', remove_sars2=False, sample_sars2=False,
    #                        trimer_continuous=True, label='Binds')   # DONE
    # scores.to_csv('./data/ab_same_sars1_binds_scores_cont.csv', index=False)

    # dataset, scores = main('ab_same_sars1_binds', remove_sars2=False, sample_sars2=False,
    #                        trimer_continuous=True, label='Binds')   # HERE
    # scores.to_csv('./data/ab_same_sars1_6_muts_at_train_binds_scores_cont.csv', index=False)

    embed, scores = main('ab_same_sars1_binds', data_from_file=True, exclude_sars2=False,
                         downsample_sars2=False, classic_trimer=True, target='Binds', plot=True)   # DONE
    # scores.to_csv('./data/ab_same_sars1_binds_del_muts_scores_cont.csv', index=False)

    # dataset, scores = main('ab_same_sars1_binds', remove_sars2=True, sample_sars2=False,
    #                       trimer_continuous=True, label='Binds')   # DONE
    # scores.to_csv('./data/ab_same_sars1_nos2_binds_scores_cont.csv', index=False)
    #
    # dataset, scores = main('ab_same_sars1_binds', remove_sars2=True, sample_sars2=False,
    #                       trimer_continuous=True, label='Binds')   # DONE
    # scores.to_csv('./data/ab_same_sars1_nos2_6_muts_at_train_binds_scores_cont.csv', index=False)
    #
    # dataset, scores = main('ab_same_sars1_binds', remove_sars2=True, sample_sars2=False,
    #                       trimer_continuous=True, label='Binds')   # DONE
    # scores.to_csv('./data/ab_same_sars1_nos2_del_muts_binds_scores_cont.csv', index=False)
    #
    # dataset, scores = main('ab_same_sars1_binds', remove_sars2=True, sample_sars2=False,
    #                       trimer_continuous=True, label='Binds')   # DONE
    # scores.to_csv('./data/ab_same_sars1_nos2_noratg13_binds_scores_cont.csv', index=False)
    #
    # dataset, scores = main('ab_same_sars1_binds', remove_sars2=True, sample_sars2=False,
    #                       trimer_continuous=True, label='Binds')   # DONE
    # scores.to_csv('./data/ab_same_sars1_nos2_nopangolins_binds_scores_cont.csv', index=False)

    # dataset, scores = main('ab_same_sars1_binds', remove_sars2=True, sample_sars2=False,
    #                       trimer_continuous=True, label='Binds')   # DONE
    # scores.to_csv('./data/ab_same_sars1_nos2_nopangolins_binds_scores_cont.csv', index=False)

    # dataset, scores = main('ab_same_sars1_binds', remove_sars2=True, sample_sars2=False,
    #                       trimer_continuous=True, label='Binds')   # DONE



    #===
    # RBD
    # dataset, scores = main('rbd_189_binds', remove_sars2=False, sample_sars2=False,
    #                       trimer_continuous=True, label='Binds')  # DONE
    # scores.to_csv('./data/rbd_189_binds_scores_cont.csv', index=False)

    # dataset, scores = main('rbd_189_binds', remove_sars2=False, sample_sars2=False,
    #                       trimer_continuous=True, label='Binds')  # DONE
    # scores.to_csv('./data/rbd_189_6_muts_at_train_binds_scores_cont.csv', index=False)

    # dataset, scores = main('rbd_189_binds', remove_sars2=False, sample_sars2=False,
    #                       trimer_continuous=True, label='Binds')  # HERE
    # scores.to_csv('./data/rbd_189_del_muts_binds_scores_cont.csv', index=False)

    # dataset, scores = main('rbd_sarbe_1261_binds', remove_sars2=False, sample_sars2=True,
    #                       trimer_continuous=True, label='Binds')  # DONE
    # scores.to_csv('./data/rbd_sarbe_1261_binds_balancedS2_scores_cont.csv', index=False)

    # dataset, scores = main('rbd_sarbe_1261_binds', remove_sars2=False, sample_sars2=True,
    #                       trimer_continuous=True, label='Binds')  # DONE
    # scores.to_csv('./data/rbd_sarbe_1261_balancedS2_del_muts_binds_scores_cont.csv', index=False)