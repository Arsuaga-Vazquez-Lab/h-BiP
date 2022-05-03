import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



# Position 356 (start=0) in the alignment corresponds to position 319 from urbani
# minimal Tor2-RBD: 319 to 510
rbd_start = 356  # 319
rbd_end = 356 + (510 - 319)  # 510

urbani = 'AY278741'
tor2 = 'AY274119'
a022 = 'AY686863'   #binds 50% tor2
rp3 = 'DQ071615'   # no binding
As6526 = 'KY417142'

bind_human = {'KF569996': 'Bat LYRa11', 'MN996532': 'RaTG13', 'KF367457': 'Bat WIV1',
              'KC881005': 'Bat RsSHC014', 'KY417150': 'Rs4874', 'KY417151': 'Rs7327',
              'KY417146': 'Rs4231', 'MZ190137': 'Khosta-1', 'MZ190138': 'Khosta-2',
              'ASL68941': 'HKU25', 'ASL68953': 'HKU25',
              'KT444582': 'Bat WIV16',
              'AWH65877': 'HKU4', 'ABN10848': 'HKU4', 'ABN10857': 'HKU4', 'ABN10866': 'HKU4',
              'AWH65888': 'HKU4', 'AWH65899': 'HKU4', 'QHA24678': 'HKU4', 'YP_001039953': 'HKU4',
              'MT040334': 'Pangolin PCoV_GX-P1E', 'MT040333': 'Pangolin PCoV_GX-P4L',
              'MT040336': 'Pangolin PCoV_GX-P5E', 'MT040335': 'Pangolin PCoV_GX-P5L',
              'AY304486': 'SZ3', 'AY686863': 'A022',
              'AY274119': 'Tor2', 'YP_009724390': 'Wuhan'}


no_bind_human = {'MG772934': 'Bat CoVZXC21', 'DQ022305': 'Bat HKU3-1',
                 'DQ412042': 'Bat Rf1', 'KY417145': 'Bat Rf4092',
                 'NC_014470': 'Bat BM48-31', 'YP_003858584': 'Bat BM48-31',
                 'YP_001039962': 'HKU5', 'AWH65932': 'HKU5', 'AWH65910': 'HKU5',
                 'AWH65921': 'HKU5', 'LC556375': 'Bat SARSr Rc-o319',
                 'DQ071615': 'bat SARS CoV Rp3', 'KY417142': 'bat SL CoV As6526'}


def plot_scores(data, start, end):
    plt.plot(data.attention.values[0])
    plt.axvline(x=start, color='black', linestyle=':')
    plt.axvline(x=end, color='black', linestyle=':')
    virus = str(data.Virus.values).replace("['",'').replace("']",'')
    plt.suptitle(virus)
    plt.title('Average attention scores -RBD between lines-')
    plt.show()


path = './data/sars1_195_attention.csv'
df = pd.read_csv(path)
df = df[['Accession', 'Virus', 'Host_agg', 'Human', 'Sequence_Aligned', 'predict', 'attention']]
print(df.dtypes)
# attention column is a string (doesn't recognize the list)
df.loc[df.Accession == 'MT084071', 'attention'] = np.nan
df["attention"] = df["attention"].fillna("[]").apply(lambda x: eval(x))   # after this is a list of floats

plot_scores(df.loc[df.Accession == As6526], rbd_start, rbd_end)
