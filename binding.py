
middle_east = ['Jordan', 'Oman', 'Qatar', 'Saudi Arabia',
                'United Arab Emirates']

# camel_bind_human = ab[(ab.Host_agg == 'Camels') & \
#                       (ab.Species_agg == 'MERSr') & \
#                       (ab.Country.isin(middle_east))]

bind_human = {'KT444582': 'Bat WIV16', 'MN996532': 'RaTG13', 'KF367457': 'Bat WIV1',
              'MZ190137': 'Khosta-1', 'MZ190138': 'Khosta-2',
              'KC881005': 'Bat RsSHC014', 'KY417151': 'Rs7327', 'KC881006': 'Rs3367',
              'KY417152': 'Rs9401', 'KY417146': 'Rs4231',  'KY417150': 'Rs4874',
              'KY417144': 'Rs4084',
              'ASL68941': 'HKU25', 'ASL68953': 'HKU25',
              'KF569996': 'Bat LYRa11', 'MK211376': 'YN2018B',
              'AWH65877': 'HKU4', 'ABN10848': 'HKU4', 'ABN10857': 'HKU4', 'ABN10866': 'HKU4',
              'AWH65888': 'HKU4', 'AWH65899': 'HKU4', 'QHA24678': 'HKU4', 'YP_001039953': 'HKU4',
              'AY304486': 'SZ3', 'AY686863': 'A022',
              'AY274119': 'Tor2', 'YP_009724390': 'Wuhan',
              'ACT11030': 'HEC', 'ACJ35486': 'HEC',
              'AHN64774': 'HKU23', 'AHN64783': 'HKU23', 'QEY10625': 'HKU23',
              'QEY10641': 'HKU23', 'QEY10649': 'HKU23', 'QEY10657': 'HKU23',
              'QEY10633': 'HKU23', 'ALA50080': 'HKU23',
              'JX163927': 'SARS1', 'JX163926': 'SARS1', 'AY545919': 'SARS1'
              }

no_bind_human = {'MG772934': 'Bat CoVZXC21', 'DQ022305': 'Bat HKU3-1',
                 'DQ412042': 'Bat Rf1', 'KY417145': 'Bat Rf4092',
                 'NC_014470': 'Bat BM48-31', 'YP_003858584': 'Bat BM48-31',
                 'YP_001039962': 'HKU5', 'AWH65932': 'HKU5', 'AWH65910': 'HKU5',
                 'AWH65921': 'HKU5', 'LC556375': 'Bat SARSr Rc-o319',
                 'DQ071615': 'bat SARS CoV Rp3', 'KY417142': 'bat SL CoV As6526'}