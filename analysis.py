import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as patches

from binding import *


pd.options.display.width = 100
pd.options.display.max_columns = 10




def scores_binding(df):
    subset = df[df.Accession.isin(bind_human)].copy()
    subset['Binds'] = 'Yes'
    temp = df[df.Accession.isin(no_bind_human)].copy()
    temp['Binds'] = 'No'
    return pd.concat([subset, temp])


def all_important_scores(df):
    refs = {'KT444582': 'WIV16','MN996532': 'RaTG13', 'AY278741': 'Urbani',
             'AY274119': 'Tor2', 'YP_009724390': 'Wuhan', 'KY417150': 'Rs4874'}
    pangolins = {'MT121216':' PCoV MP789', 'MT799524': ' PCoV cDNA18', 'MT799523':' PCoV cDNA16',
                'MT072864':' PCoV GX-P2V', 'MT040336':'PCoV GX-P5E', 'MT040335': ' PCoV GX-P5L',
                'MT040334':' PCoV GX-P1E','MT040333':' PCoV GX-P4L', 'MT799521':' PCoV cDNA8'}
    other = {'ARC95227': 'PHEV','AVN88332' : 'Water deer'}
    all_imp = refs | pangolins | other
    df_all_imp = pd.DataFrame({'Accession':all_imp.keys(), 'Virus':all_imp.values()})
    imp_scores = df_all_imp.merge(df[df.Accession.isin(all_imp)], on='Accession')
    return imp_scores


def add_max_identity(path):
    id = pd.read_csv('./data/ab_same_sars1_identity_vs_hcov.csv')
    scores = pd.read_csv(path)
    id = id.merge(scores, on='Accession')
    return id


def plot_correlation(df, fig_name=None, size='half'):
    """size refers to letter
    fig_name: string
        file name to save the plot. if None, it won't be saved"""
    if size == 'full':
        figsize = (6, 4)
    elif size == 'half':
        figsize = (3.3, 2.5)
    else:   # best for display
        figsize = (10, 6)

    plt.subplots(figsize=figsize)
    plt.rc('font', size=8)
    sns.scatterplot(x='MaxSim', y='Binds_prob', data=df, hue='Binds',
                   palette=['dimgray', 'indianred'], style='Binds', s=5)  # smaller
    lgnd = plt.legend(labels=['Binds', "Doesn't bind"], bbox_to_anchor=(.9, 1.2), ncol=2)
    lgnd.legendHandles[0]._sizes = [40]
    lgnd.legendHandles[1]._sizes = [30]
    plt.xlabel('Max % identity hCoV')
    plt.ylabel('Human-infectivity score')
    plt.axhline(y=.5, linestyle='--')
    plt.axvline(x=97, linestyle='--')
    plt.plot(97.46, 0.999, color='steelblue', marker='x', markersize=2)  # RaTG13

    if fig_name != None:
        plt.savefig('./data/' + fig_name + '.pdf', bbox_inches='tight')

    plt.show()


def plot_score_distribution_by_target(df):
    exclude = ['MT084071', 'YP_009724390']  # spurious pangolin and wuhan
    fig, ax = plt.subplots(1, 1, figsize=(3.3, 2.5))
    plt.rc('font', size=6)
    histbins = [0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1]
    subset_nh = (df.Human == 0) & (~df.Accession.isin(exclude))
    subset_h = (df.Human == 1) & (~df.Accession.isin(exclude))
    df.loc[subset_nh, 'Human_prob'].hist(bins=histbins, hatch='///', color=['white'],
                                                   edgecolor='black', grid=False, label='Non-human')
    df.loc[subset_nh, 'Human_prob'].hist(bins=histbins, color=['white'],
                                                   edgecolor='black', grid=False, alpha=0.5)
    df.loc[subset_h, 'Human_prob'].hist(bins=histbins, alpha=0.5, color=['indianred'],
                                                  edgecolor='black', grid=False, label='Human')
    plt.xlabel('Human-infection score')
    plt.ylabel('Count')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # rect = patches.Rectangle((.5, .5), .61, 50.5, linewidth=1.5, linestyle='dashed',
    #                          edgecolor='darkred', facecolor='none')
    # ax.add_patch(rect)
    plt.legend(title='Host', bbox_to_anchor=(.61, .7))

    # plt.savefig('../MyPaper/sars1_v3_hist_all.pdf', bbox_inches='tight')

    plt.show()


def plot_selected_binding_trip_vs_cont(triplicates, continuous, dataset='ab_same_sars1'):
    binding_trip = scores_binding(triplicates)
    binding_cont = scores_binding(continuous)
    sns.scatterplot(binding_cont.Human_prob, binding_trip.Human_prob, hue=binding_trip.Binds, s=8)
    plt.text(.27, .99, '*', fontsize=7)   # QHA24710 HKU10
    plt.text(.01, .99, '*', fontsize=7)   # ABG47052 BtCoV/133
    plt.text(.001, .85, '*', fontsize=7)   # AHY61337 BtVs
    plt.xlabel('Score from continuous')
    plt.ylabel('Score from triplicates')
    plt.title('Selected viruses ' + dataset)
    plt.legend(title='Binds to human receptor', loc='lower right')
    plt.show()


if __name__ == '__main__':

     #  HUMAN TRIPLICATES
    # path = './data/ab_same_sars1_scores.csv'
    # path = './data/ab_same_sars1_balancedS2_scores.csv'
    # path = './data/ab_same_sars1_nos2_scores.csv'
    # path = './data/rbd_189_scores.csv'

    # HUMAN CONTINUOUS
    # path = './data/ab_same_sars1_balancedS2_scores_cont.csv'
    # path = './data/ab_same_sars1_nos2_scores_cont.csv'
    # path = './data/rbd_189_scores_cont.csv'

    # BINDS TRIPLICATES
    # path = './data/ab_same_sars1_binds_scores.csv'
    # path = './data/ab_same_sars1_nos2_binds_scores.csv'


    # BINDS CONTINUOUS
    # path = './data/ab_same_sars1_binds_scores_cont.csv'
    # path = './data/ab_same_sars1_nos2_binds_scores_cont.csv'
    # path = './data/ab_same_sars1_6_muts_at_train_binds_scores_cont.csv'
    path = './data/ab_same_sars1_binds_del_muts_scores_cont.csv'
    # path = './data/ab_same_sars1_nos2_6_muts_at_train_binds_scores_cont.csv'
    # path = './data/ab_same_sars1_nos2_del_muts_binds_scores_cont.csv'
    # path = './data/rbd_189_binds.csv'
    # path = './data/rbd_189_binds_scores_cont.csv'
    # path = './data/rbd_189_6_muts_at_train_binds_scores_cont.csv'
    # path = './data/rbd_189_del_muts_binds_scores_cont.csv'
    # path = './data/rbd_sarbe_1261_binds_scores_cont.csv'
    # path = './/data/rbd_sarbe_1261_binds_balancedS2_scores_cont.csv'
    # path = './data/rbd_sarbe_1261_balancedS2_del_muts_binds_scores_cont.csv'


    # add_max_identity(path).to_csv(path, index=False)

    # trip = pd.read_csv('./data/ab_same_sars1_scores.csv')
    # cont = pd.read_csv('./data/ab_same_sars1_scores_cont.csv')
    # trip = pd.read_csv('./data/ab_same_sars1_balancedS2_scores.csv')
    # cont = pd.read_csv('./data/ab_same_sars1_balancedS2_scores_cont.csv')
    # trip = pd.read_csv('./data/rbd_189_scores.csv')
    # cont = pd.read_csv('./data/rbd_189_scores_cont.csv')
    cont = pd.read_csv(path)


    # trip_scores = scores_binding(trip)
    # print('\nScore from selected viruses: Triplicates approach:\n',
    # trip_scores[['Binds', 'Accession', 'Virus', 'Human_prob']].sort_values(by=['Binds', 'Human_prob']))

    # cont_scores = scores_binding(cont)
    # print('\nScore from selected viruses: Continuous approach:\n',
    # cont_scores[['Binds', 'Accession', 'Virus', 'Human_prob']].sort_values(by=['Binds', 'Human_prob']))

    # plot_selected_binding_trip_vs_cont(trip, cont, dataset='ab_same_sars1')
    # plot_selected_binding_trip_vs_cont(trip, cont, dataset='ab_same_sars1_balancedS2')
    # plot_selected_binding_trip_vs_cont(trip, cont, dataset='RBD_189')

    # print(binds.Human.value_counts(),'\n')
    # print(cont.Host_agg.value_counts(), '\n')
    # print(cont.Species_agg.value_counts(), '\n')
    # print(pd.crosstab(binds.Species_agg, binds.Human),'\n')
    # # print(binds[(binds.Human == 1) & (binds.Binds == 0)])
    # print('Non-humans:')
    # print(binds[(binds.Host_agg != 'Homo sapiens')].Human.value_counts(),'\n')
    # print('Humans:')
    # print(binds[(binds.Host_agg == 'Homo sapiens')].Human.value_counts(), '\n')

    # rbd = cont.loc[cont.Accession == 'MN996532', 'Sequence']  # series
    # rbd = cont.loc[cont.Accession == 'MN996532', 'Sequence'].values   # np
    # rbd = str(cont.loc[cont.Accession == 'MN996532', 'Sequence'].values)
    # print(type(rbd))
    # rbd = 'ITNLCPFGEVFNATTFASVYAWNRKRISNCVADYSVLYNSTSFSTFKCYGVSPTKL\
    # NDLCFTNVYADSFVITGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSKHIDAKEGGNFNYLYRLFRKANLKPFER\
    # DISTEIYQAGSKPCNGQTGLNCYYPLYRYGFYPTDGVGHQPYRVVVLSFELLNAPATV'
    # print(cont.loc[cont.Accession == 'MN996532'])
    # print(type(cont.Sequence[0]))
    # print(cont[cont.Sequence == rbd])
    # print(str(cont.loc[cont.Accession == 'MN996532', 'Sequence']))

    print('\nCorrelation {:.2f}\n'.format(cont.MaxSim.corr(cont.Binds_prob)))

    plot_correlation(cont, fig_name='ab_same_sars1_binds_del_muts_correlation', size='half')
    # plot_correlation(cont, size='half')

     # plot_score_distribution_by_target(ab_scores)




