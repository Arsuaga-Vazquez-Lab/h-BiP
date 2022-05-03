import pandas as pd

from binding import *

#  IMPORTANT: Human columns will be dropped. A new target for binding (yes or no)
# will be created and renamed as "Human" to facilitate integration with the rest
# of the code.   ALREADY STARTING FIXING THIS
#
# ab = pd.read_csv('./data/ab_same_sars1.csv')
# camels_77 = pd.read_csv('./data/camels_mersr_77.csv')
#
# hCoV = ['MERS', 'hCoV229E', 'hCoVHKU1', 'hCoVNL63', 'hCoVOC43',
#         'SARS-CoV-1', 'SARS-CoV-2']
#
# add_country_mers = {'AVN89365': 'Burkina Faso', 'AVN89376': 'Nigeria',
#                     'ALA49352': 'Saudi Arabia'}
#
# middle_east = ['Jordan', 'Oman', 'Qatar', 'Saudi Arabia',
#                 'United Arab Emirates']
#
# ab = ab.merge(camels_77[['Accession', 'Country']], on='Accession', how='left')
#
# for acc in add_country_mers:
#     ab.loc[ab.Accession == acc, 'Country'] = add_country_mers[acc]
#
# # Create column 'Binds' with 1 if human coronavirus, or binding evidence
# camel_bind_human = ((ab.Host_agg == 'Camels') &
#                     (ab.Species_agg == 'MERSr') &
#                     (ab.Country.isin(middle_east)))
#
# binds = ((ab.Species_agg.isin(hCoV)) |
#          (ab.Accession.isin(bind_human)) |
#          camel_bind_human |
#          (ab.Host_agg == 'Civet') |
#          (ab.Host_agg == 'Pangolin') |
#          (ab.Accession == 'ABI93999'))   # Bovine coronavirus isolate Alpaca
#
# ab['Binds'] = 0
# ab.loc[binds, 'Binds'] = 1
#
# print(pd.crosstab(ab.Human, ab.Binds))
#
# print('\n', pd.crosstab(ab.Host_agg, ab.Binds))
#
# print(ab[(ab.Host_agg == 'Bovinae other') & (ab.Binds == 1)])
#
# # ab.drop(columns='Human', inplace=True)
# # ab.rename(columns={'Binds': 'Human'}, inplace=True)
# print(ab.columns)
# print(pd.crosstab(ab.Human, ab.Binds))
#
# # ab.to_csv('./data/ab_same_sars1_binds.csv', index=False)

# GENERATE RBD_189 SUBSET
# ab = pd.read_csv('./data/ab_same_sars1_binds.csv')
# rbd = pd.read_csv('./data/rbd_189.csv')
# rbd_binds = rbd.merge(ab[['Accession', 'Binds']], on='Accession', how='left')
# rbd_binds.loc[rbd_binds.Binds.isnull(), 'Binds'] = 1
# print('nans', rbd_binds.Binds.isnull().sum())
# print(pd.crosstab(rbd_binds.Human, rbd_binds.Binds))
# rbd_binds.to_csv('./data/rbd_189_binds.csv', index=False)

# GENERATE RBD SARBE SUBSET
ab = pd.read_csv('./data/ab_same_sars1_binds.csv')
ab.drop(columns=['Sequence'], inplace=True)

s1 = pd.read_csv('./data/rbd_189.csv')
s1 = s1[['Accession', 'Sequence']]
s2 = pd.read_csv('./data/found_rbd_sars2_1073.csv')
s2 = s2[['name', 'rbd']]
s2.rename(columns={'name': 'Accession', 'rbd': 'Sequence'}, inplace=True)
#
sarbe = pd.concat([s1, s2], axis=0)
rbd_binds = sarbe.merge(ab, on='Accession', how='left')
rbd_binds = rbd_binds[rbd_binds.Accession != 'MT084071']   #1261
# rbd_binds.to_csv('./data/rbd_sarbe_1261_binds.csv', index=False)

