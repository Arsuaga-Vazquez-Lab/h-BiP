import pandas as pd
import numpy as np
import pickle
import random

def read_seqs(_file_path):
    # todo: would be nice to save output to a file
    """ Read text file (a fasta file) with sequences and description containing:
            Accession.version, GenBank_Title, Species and Host
        Classify host in three groups: bat, human u other
        Provides a list with description for every sequence, a 2D array with sequences and a list with the group of host

    Parameters
    ----------
    _file_path: str
        Path to fasta file

    Returns
    -------
    _df : dataframe
        Dataframe with columns: Accession, GenBank_Title, Species, Host, Sequence, Target and Virus.
    _df_unique : dataframe
        Same dataframe as _df after removing duplicates for Sequence
    """

    print("\nReading Sequences...\n")
    seq = []
    meta = []
    read = ""
    target = []
    virus = []

    def detect_host(_description):
        """ Recognize three different groups from the host: bat, human and other. If host is missing will assign NaN

        Parameters
        ----------
        _description : list
            Description for the sequence from NCBI: Accession.version, GenBank_Title, Species, Host

        Returns
        -------
        _group : str
            Either "Bat", "Human", "Other". NaN if missing.
        """

        bat = ['Aselliscus stoliczkanus', 'Chaerephon plicatus', 'Chiroptera', 'Cynopterus sphinx', 'Desmodus rotundus',
               'Eonycteris spelaea', 'Hipposideros', 'Hipposideros larvatus', 'Hipposideros pomona',
               'Hipposideros pratti', 'Hypsugo pulveratus', 'Hypsugo savii', 'Macronycteris vittata',
               'Miniopterus fuliginosus', 'Miniopterus pusillus', 'Miniopterus schreibersii', 'Myotis brandtii',
               'Myotis lucifugus', 'Myotis ricketti', 'Myotis sp.', 'Neoromicia capensis', 'Nyctalus velutinus',
               'Pipistrellus abramus', 'Pipistrellus kuhlii', 'Rhinolophus', 'Rhinolophus affinis',
               'Rhinolophus blasii', 'Rhinolophus ferrumequinum', 'Rhinolophus macrotis', 'Rhinolophus pusillus',
               'Rhinolophus sinicus', 'Rousettus aegyptiacus', 'Rousettus leschenaultii', 'Rousettus sp.',
               'Scotophilus kuhlii', 'Suncus murinus', 'Triaenops afer', 'Tylonycteris pachypus',
               'Tylonycteris robustula', 'Vespertilio sinensis', 'Eidolon helvum', 'Pipistrellus',
               'Cynopterus brachyotis', 'Myotis siligorensis', 'Myotis davidii', 'Myotis daubentonii']
        if _description[3] in bat:
            _group = "Bat"
        elif _description[3] == 'Homo sapiens':
            _group = 'Human'
        elif _description[3] == '':
            if 'bat' in _description[1].lower():
                _group = "Bat"
            elif 'human' in _description[1].lower():
                _group = 'Human'
            elif 'murine' in _description[1].lower():
                _group = "Other"
            elif 'bovine' in _description[1].lower():
                _group = "Other"
            elif 'equine' in _description[1].lower():
                _group = "Other"
            else:
                _group = 'NaN'
        else:
            _group = 'Other'
        return _group

    def detect_virus_type(_description, _host_group):
        """ Classify alphacoronavirus and betacoronavirus into:
            Human:     1) HCoV-NL63',  2) 'HCoV-229E', 3) 'SARS-CoV-2', 4) 'HCoV-HKU1', 5) 'MERS',
                       6) 'HCoV-OC43', 7) 'SARS-CoV-1', 8) 'Other-human-CoV'
            Non-human: 1) 'SARS-related-non-human', 2) 'MERS-related-non-human', 3) 'Other-non-human-CoV'
            Other:     Other virus from unknown host

            Parameters
            ----------
            _description : list
                NCBI sequence's description: Accession.version, GenBank_Title, Species, Host
            _host_group : str
                Host output from detect_host() with: Human, Bat, Other

            Returns
            -------
            _virus : str
                Returns 8 different types of human coronavirus and 3 types of non-human coronavirus
            _host_adjusted : str
                If original _host_group was missing and _virus turns to be either SAS-CoV-1 or SARS-CoV-2,
                _host_adjusted will be imputed as 'Human'
            """
        _host_adjusted = _host_group
        _GenBank = _description[1].replace('[', ' ').replace(']', ' ').replace('/', ' ')
        if _host_group == 'Human':
            if 'NL63' in _GenBank:
                _virus = 'HCoV-NL63'
            elif '229E' in _GenBank:
                _virus = 'HCoV-229E'
            # elif 'SARS-CoV-2' in _description[1].replace('/', ' '):
            elif 'Severe acute respiratory syndrome coronavirus 2' in _GenBank:
                _virus = 'SARS-CoV-2'
            elif 'HKU1' in _GenBank:
                _virus = 'HCoV-HKU1'
            elif 'Middle' in _GenBank:
                _virus = 'MERS'
            elif 'OC43' in _GenBank:
                _virus = 'HCoV-OC43'
            # elif 'SARS coronavirus' in _GenBank:
            #     _virus = 'SARS-CoV-1'
            elif 'Severe acute respiratory syndrome coronavirus' in _GenBank:
                _virus = 'SARS-CoV-1'
            else:
                _virus = 'Other-human-CoV'
        elif _host_group == 'NaN':
            if 'Severe acute respiratory syndrome coronavirus 2' in _GenBank:
                _virus = 'SARS-CoV-2'
                _host_adjusted = 'Human'
            elif 'SARS coronavirus' in _GenBank:
                _virus = 'SARS-CoV-1'
                _host_adjusted = 'Human'
            else:
                _virus = 'Other'
        else:
            # _description[2] is Species
            if 'Severe' in _description[2]:
                _virus = 'SARS-related-non-human'
            elif 'Middle' in _description[2]:
                _virus = 'MERS-related-non-human'
            else:
                _virus = 'Other-non-human-CoV'
        return _virus, _host_adjusted

    def stats(_target, _virus):
        """ Print absolute and relative frequencies for target (host type) and virus type.

        Parameters
        ----------
        _target : list
            Target variable: Human, Non-human, Bat, Other
        _virus : list
            Virus type: 8 types o human coronavirus, 3 types of non-human coronavirus, and "Other"
            """
        # Stats for target (Host)
        print('---- FREQUENCIES HOST -----\n')
        _df = pd.DataFrame(_target)
        print('Absolute:')
        _abs = _df.value_counts()
        _tot = _df.shape[0]
        print('Total sequences = ', _tot, '\n')
        print(_abs)
        print('\nRelative (%): \nIs the dataset balanced?\n')
        print((_abs / _tot) * 100)
        print('-----------------------------\n')
        # Stats for virus type
        print('---- FREQUENCIES VIRUS -----\n')
        _df = pd.DataFrame(_virus)
        print('Absolute:')
        _abs = _df.value_counts()
        print(_abs)
        print('-----------------------------\n')

    def create_df(_meta, _seq, _target, _virus):
        """" Create a dataframe from all the data

        Parameters
        ----------
        _meta : list
            Description for every sequence: Accession.version, GenBank_Title, Species and Host
        _seq: list
            list of strings containing all sequences
        _target: list
            The group of the host: Bat, Human u Other. NaN if missing
        _virus: list
            8 different types of human coronavirus and 3 types of non-human coronavirus

        Returns
        -------
        dataframe
            Dataframe with columns: Accession, GenBank_Title, Species, Host, Sequence, Target and Virus
            """
        # convert array to data frame with headers
        _meta_df = pd.DataFrame(_meta, columns=['Accession', 'GenBank_Title', 'Species', 'Host'])
        # removing trailing spaces in Accession
        _meta_df.Accession = _meta_df.Accession.str.rstrip()
        _data_df = pd.DataFrame({'Sequence': _seq, 'Target': _target, 'Virus': _virus})
        return pd.concat([_meta_df, _data_df], axis=1)

    with open(_file_path) as f:
        lines = f.readlines()
    last_line = lines[-1]
    not_this_meta = 0
    not_this_seq = 0
    for line in lines:
        if not line.startswith('>'):
            read += line.rstrip('\n')
        if line.startswith('>') or line == last_line:
            if line != last_line:
                # Accession.version, GenBank_Title, Species, Host
                clean_meta = line[1:].strip('\n').split('|')
                meta.append(clean_meta)
                host = detect_host(clean_meta)
                # host gets imputed for some missing values
                virus_type, host = detect_virus_type(clean_meta, host)
                target.append(host)
                virus.append(virus_type)
            if read != "":
                seq.append(read)
                read = ""
    stats(target, virus)
    _df = create_df(meta, seq, target, virus)
    return _df


def read_pickled_all_raw():
    """ Reading alpha_beta_all_raw.pkl

    Parameters
    ----------
    None

    Returns
    -------
    dataframe
        Dataframe with raw sequences from alphacoronavirus and betacoronavirus.
        Columns are: Accession, GenBank_Title, Species, Host, Sequence, Target and Virus.
        """
    print('Reading pickle data : alpha_beta_all_raw.pkl ...')
    with open('alpha_beta_all_raw.pkl', 'rb') as f:
        return pickle.load(f)


def clean_raw_data(_df_raw):
    """ Remove duplicates for Sequence and remove samples with missing data for Target

    Parameters
    ----------
    _df_raw : dataframe
        Initial dataframe with at least columns Sequence and Target
    Returns
    -------
    _df_clean : dataframe
        Same dataframe as _df_raw after removing duplicates and missing values
    """

    print('\n------ UNIQUENESS -------\n')
    print(_df_raw.nunique())
    _duplicates = _df_raw.shape[0] - _df_raw.Sequence.nunique()
    print('\n*** There are ', _duplicates, 'duplicated sequences ***\n')
    _df_clean = _df_raw.drop_duplicates(subset="Sequence")
    print('\n------ TARGET -------\n')
    print(_df_clean.Target.value_counts())
    _df_clean = _df_clean.replace('NaN', np.NaN)
    print('\n*** ', _df_clean.Target.isnull().sum(), 'missing values from Target variable will be removed ***\n')
    _df_clean = _df_clean.dropna()
    print('\n Final number of samples in dataset: ', _df_clean.shape[0], '\n')
    return _df_clean


def down_sampling_cov2(_df_unique):
    """ Down-sampling Virus==SARS-CoV-2 samples by 1150 and removing two incomplete sequences.

    Parameters
    ----------
    _df_unique : dataframe
        Dataframe with columns: 'Accession', 'GenBank_Title', 'Species', 'Host', 'Sequence', 'Target', 'Virus'.

    Returns
    -------
    _df_subset : dataframe
        Same dataframe as _df_unique after removing 1150 samples where Virus==SARS-CoV-2.
    """
    print('Down-sampling SARS-CoV-2 samples by 1150 ...')
    _cov2_accessions = _df_unique.Accession[_df_unique.Virus == 'SARS-CoV-2']
    print('SARS-CoV-2 samples: ', _cov2_accessions.shape[0])
    random.seed(0)  # for reproducibility purposes
    _cov2_subset_accessions = random.sample(list(_cov2_accessions), 1150)
    # Convert dataframe to list
    _df_unique_to_list = _df_unique.values.tolist()
    # sample[0] corresponds to Accession
    _unique_subset_list = [sample for sample in _df_unique_to_list if sample[0] not in _cov2_subset_accessions]
    print('Original data set number of samples: ', _df_unique.shape[0])
    _partial_data = ['AYR18479.1', 'AYR18511.1']
    _unique_subset_no_partial_list = [sample for sample in _unique_subset_list if sample[0] not in _partial_data]
    _df_subset = pd.DataFrame(_unique_subset_no_partial_list,
                             columns=['Accession', 'GenBank_Title', 'Species', 'Host', 'Sequence', 'Target', 'Virus'])
    print('Number of samples after down-sampling: ', _df_subset.shape[0])
    return _df_subset


if __name__ == '__main__':
   # NCBI complete amino acid sequences from alphacoronavirus and betacoronavirus for the spike
    file_path = "../NCBI/new/alpha-beta-spike-110520.fasta"

    # df_all_raw = read_seqs(file_path)

    # # PICKLE IT:
    # # ===============================================
    # with open('alpha_beta_all_raw.pkl', 'wb') as f:
    #     pickle.dump(df_all_raw, f)
    # #===============================================

    df_all_raw = read_pickled_all_raw()
    df_raw_unique = clean_raw_data(df_all_raw)
    df_final_raw = down_sampling_cov2(df_raw_unique)

    # PICKLE IT:
    # ===============================================
    with open('alpha_beta_final_raw.pkl', 'wb') as f:
        pickle.dump(df_final_raw, f)
    #===============================================



