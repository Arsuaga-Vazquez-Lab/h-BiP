from Bio.Seq import Seq
from Bio import SeqIO
import pandas as pd
import numpy as np
import re

    
def extract_spike(_seq):
    """Extract spike protein from a coronavirus using the first 9 residues from conserved domain from S2 
    (NCBI accession cl40439) allowing for one substitution. After the domain is located, the sequence is truncated
    to the left and to the right at the stop codon and trimmed on the left 
    to the start codon M
    
    Parameters
    ----------
    _seq: Bio.Seq.Seq
        A string with the nucleotide sequence from where to extract the spike
    
    Returns
    -------
    spike: Bio.Seq.Seq
        A string with the aminoacid sequence for the spike. Empty Seq if the spike wasn't found
    found: int
        Will return 1 or 0 if the spike was found (or not)
    """
    def find_domain(s):
#         s2d = "SFIEDLLFNKVTLADAGF"   # around position 798-815
#         Using half of the domain resulted in 9 more records s2d = "SFIEDLLFN"
        s2d_1mutation = [".FIEDLLFN", "S.IEDLLFN", "SF.EDLLFN","SFI.DLLFN","SFIE.LLFN",
                         "SFIED.LFN","SFIEDL.FN","SFIEDLL.N","SFIEDLLF."]

        for p in s2d_1mutation:
            try:
                res= re.finditer(p, s)
                return next(res).start()
            except:
                pass
        return None
  
    # if sequence not multiple of 3 add either 'n' or 'nn'
    excess = len(_seq) % 3
    if excess == 0:
        adj_seq = _seq + 'nn'
    elif excess == 1:
        adj_seq = _seq + 'n'
    else:
        adj_seq = _seq
    # Translation for the 3 possible ORFs
    found = 0
    ORFs = [adj_seq[:-2].translate(), adj_seq[1:-1].translate(), adj_seq[2:].translate()]
    for ORF in ORFs:
        loc = find_domain(str(ORF))
        if loc != None:
            protein = ORF
            found += 1
            break
    if found == 0:
        spike = Seq('')
    else:
        loc_stop_right = protein.find('*',start=loc)
        loc_stop_left = protein.rfind('*',end=loc)
        # Trim the protein in between stops, but take full length when no stop
        if loc_stop_right == -1:
            if loc_stop_left == -1:
                spike = protein
            else:
                spike = protein[loc_stop_left :]
        elif loc_stop_left == -1:
            spike = protein[:loc_stop_right]
        else:
            spike = protein[loc_stop_left : loc_stop_right]
        # If we want to trim the beginning up to M (start residue)
        start = spike.find('M')
        if start != -1:
            spike = spike[start:]
    return spike, found


def create_spike_fasta(input_path, output_path):
    record_iterator = SeqIO.parse(input_path,"fasta")
    # for seq_record in SeqIO.parse(kan_path,"fasta"):  # Another option for parsing
    total_seqs = 0
    found = 0
    spike_list = []
#     for i in range(2):  # for testing and debugging
    while True:
        try:
            seq_record = next(record_iterator)
        except:
            break
        total_seqs += 1
        spike_seq, flag = extract_spike(seq_record.seq)
        if flag == 0:
            continue
        # Create new SeqRecord
        found += 1
        spike_rec = seq_record
        spike_rec.description = spike_rec.description + str('. Extracted spike using cl40439 conserved domain from S2')
        spike_rec.seq = spike_seq
        spike_list.append(spike_rec)

    print('Spike found: {} out of {}'.format(found, total_seqs))

    SeqIO.write(spike_list, output_path, "fasta")

    
def save_accessions_from_fasta(input_path, output_path):
    record_iterator = SeqIO.parse(input_path,"fasta")
    # for seq_record in SeqIO.parse(kan_path,"fasta"):  # Another option for parsing
    total_seqs = 0
    acc_list = []
#     for i in range(5):  # for testing and debugging
    while True:
        try:
            seq_record = next(record_iterator)
        except:
            break
        total_seqs += 1
        this_seq = [seq_record.id, len(seq_record)]
        acc_list.append(this_seq)

    print('Total records: {}'.format(total_seqs))

    return np.savetxt(output_path, acc_list, fmt ='% s')


def fasta_to_df(fasta_path):
    """ Creates a dataframe from a fasta file with Accession and Sequence columns
    """
    record_iterator = SeqIO.parse(fasta_path,"fasta")
    Accession = []
    Sequence =[]
    while True:
        try:
            seq_record = next(record_iterator)
        except:
            break
        Accession.append(seq_record.id)
        Sequence.append(seq_record.seq)
    return pd.DataFrame({'Accession': Accession, 'Sequence': Sequence})


def create_SARS1_dataset(to_remove=None):
    """ Merges manually curated host data with NCBI Entrez nucleotide database for sarbecovirus
    excluding SARS-CoV-2 from 07/26/2021. Creates column 'Host_agg' that aggregates the different host into:
    [Homo sapiens, Vero cells, Patent human, Bat, Civet, Mus musculus, Pangolin, Ferret, Other]. The final
    data set consists of 561 SARS1 records
    
    Parameters
    ----------
    to_remove: list
        List of host to be removed from the data set. Default=None, 
        Consider using to_remove = ['Vero cells', 'Patent human']
    Returns
    -------
    final_host: DataFrame
        Final data frame with columns: ['Accession', 'Sequence', 'Species', 'Virus', 'Note', 'Host',
       'Host_agg', 'Human']
    """

    def remove_cov19():
        """ Remove SARS-CoV-2 from sarbecovirus dataset downloaded from NCBI 07/27/2021
        with 993,959 records for a total of 1,597 records to obtain the field 'Host' when
        available"""
        hosts_path ='./sarbecovirus_hosts_993959.csv'
        host = pd.read_csv(hosts_path)
        covid19 = 'Severe acute respiratory syndrome coronavirus 2'
        host['COV2'] = 0
        host.loc[host.GenBank_Title.str.contains(covid19,regex=False, na=False),'COV2']=1
        noCov2 = host[host.COV2!=1].copy()
        noCov2.drop(columns='COV2', inplace=True)
        noCov2['Virus'] = noCov2.GenBank_Title.str.split(',').str[0]  # only keep first element
        noCov2.drop(columns='GenBank_Title', inplace=True)
        return noCov2

    
    def remove_hosts(df, host_list):
        initial = df.shape[0]
        for host in host_list:
            to_remove = df[df.Host == host].shape[0]
            df = df[df.Host != host]
            print('Removed {} {}'.format(to_remove, host))
        print('\n')
        removed = initial - df.shape[0]
        humans = df.Host[df.Host == 'Homo sapiens'].shape[0]
        print('There are {} humans'.format(humans))
        rest = df.shape[0] - humans - df.Host.isna().sum()
        print('There are {} non-humans\n'.format(rest))
        return df
    
    
    def aggregate_host(df):
        """ Aggregates the Hosts into only 8 different host
        Returns
        -------
        df: DataFrame
            With new feature 'Host_agg'
        """
        df['Host_agg'] = df['Host'].copy()
        Civets = ['Paradoxurus hermaphroditus', 'Palm civet', 'Paguma larvata', 'Viverridae', 'Civet']
        Pangolins = ['Pholidota', 'Manis javanica' ]
        Ferrets = [ 'Ferret', 'Melogale moschata']
        Mice = ['Mus musculus','BALB/c mice']
        Others = ['Raccoon dog', 'Chlorocebus aethiops']
        for c in Civets:
            df.Host_agg.replace(c, 'Civet', inplace=True)
        for p in Pangolins:
            df.Host_agg.replace(p, 'Pangolin', inplace=True)
        for f in Ferrets:
            df.Host_agg.replace(f, 'Ferret', inplace=True)
        for m in Mice:
            df.Host_agg.replace(m, 'Mouse', inplace=True)
        for o in Others:
            df.Host_agg.replace(o, 'Other', inplace=True)
        not_this = ['Homo sapiens', 'Civet', 'Pangolin', 'Ferret', 'Mouse', 'Other',
                   'Vero cells', 'Patent human']
        df.loc[df.Host_agg.isin(not_this) == False, 'Host_agg'] = 'Bat'
        return df


    # Retrieve sarbecovirus data from NCBI to get available host and remove Covid19 (07/27/2021)
    dfhost = remove_cov19()   # 1597 records, SARS1_hosts_NCBI_072721
    # Sequence data from NCBI Entrez nucleotide database
    SARS1_spike = '../SARS1/sarbe_noCov2/sars1_spike_26jul21_560.fasta'   # Accesion includes version
    
    # Extract the spike using S2 conserved domain. 
    # Data extracted using Entrez Nucleotide database. Accession|GenBank_title and Sequence
    SARS1 = '../SARS1/sarbe_noCov2/sars1_1544_26jul21.fasta'
    # create_spike_fasta(SARS1, SARS1_spike)   # uncomment if you do not have the file
    dfspike = fasta_to_df(SARS1_spike)
    dfspike['Accession'] = dfspike['Accession'].str.split('.').str[0]
    
    # Merge both files to add host information to the spike
    spike_host = dfspike.merge(dfhost, on='Accession', how='left')  # 560 most missing the host
    
    # Manually curated the host field for spike_host records. Retrieve as host_man
    SARS1_host = '../SARS1/sarbe_noCov2/sars1_spike_560_26jul21_host_manual.csv'  # 560 host manually annotated
    host_man = pd.read_csv(SARS1_host, sep=',')
    host_man.Accession = host_man.Accession.str.split('.').str[0]
    
    # Merge file with sequence for spike (and some hosts) with manually curated hosts
    final_host = spike_host.merge(host_man, on='Accession', how='left')
    
    # Combine manually curated hosts with oficial hosts
    final_host['Host_comb'] = final_host['Host_manual'].copy()
    host_nonan = ~final_host.Host.isna()
    final_host.loc[host_nonan,'Host_comb'] = final_host.loc[host_nonan,'Host'].copy()
    final_host.drop(columns=['Host_manual', 'Host'], inplace=True)
    final_host.rename(columns={'Host_comb':'Host'}, inplace=True)
    
    # Drop rows with missing values for host
    final_host.drop(final_host[final_host.Host.isna()].index, inplace=True)
   
    final_host = aggregate_host(final_host)
    
    # Create target feature: 'Human'
    final_host['Human'] = final_host['Host_agg'].copy()
    final_host.loc[final_host.Human != 'Homo sapiens', 'Human'] = 0
    final_host.loc[final_host.Human == 'Homo sapiens', 'Human'] = 1
    
    if to_remove is not None:
        final_host = remove_hosts(final_host, to_remove)
    
    print('\n Total records in dataset: ', final_host.shape[0], '\n')
    print('Available features: ', final_host.columns, '\n')
    print(final_host.Host_agg.value_counts(),'\n')
    print(final_host.Human.value_counts(),'\n')
    
    return final_host


def selected_sequences():
    """ For SARS1 (sars1_spike_26jul21_560.fasta), comparing against these records was enough
    to remove all duplicates. From duplicates, either the most relevant was kept or the one with the
    smallest number in the sequence for the strain. If different host, kept original host 
    (as opposed to new host infected with it)
    """
    # Vero cells labeled human host
    keep_refseqs_human_vero = {'AY278741': 'Urbani', 'AB257344': 'Frankfurt 1'}
    keep_seqs_human = {'AY278554': 'CUHK-W1',
                'AY394997': 'ZS-A', 'AY278489': 'GD01', 'AY274119': 'Tor2',
                'AY310120': 'FRA', 'AY278490': 'BJ03', 
                'AY525636': 'GD03T0013', 'AY394981': 'HGZ8L1-A', 'AY304495': 'GZ50', 
                'AY304491': 'GZ60', 'AY394999': 'LC2', 'DQ412623': 'mut S nonfunctional',
                'DQ412622': 'mut S nonfunctional', 
                'DQ412619': 'mut S nonfunctional', 'DQ412615': 'CUHKtc39TI', 'DQ412614': 'CUHKtc38NP', 
                'GU553363': 'HKU-39849 isolate TCVSP-HARROD-00001', 'AY394987': 'HZS2-Fb',
                'AY279354': 'BJ04', 'EU371560': 'BJ182a', 'AY559081': 'Sin842',
                'AY394986': 'HSZ-Cb','AY394985': 'HSZ-Bb', 'CS123510': 'Patent human',
                'CS244437': 'Patent human', 'DJ059769': 'Patent human', 'AY463059': 'ShanghaiQXC1'}
    keep_seqs_civet = {'AY572038': 'civet020', 'AY686863': 'A022', 'AY304488': 'SZ16',
                       'AY304486': 'SZ3', 'AY687372': 'C029', 'AY687371': 'C028', 'AY687365': 'C013',
                       'AY687361': 'B029', 'AY687359': 'B012', 'AY687355': 'A013', 'AY687354': 'A001',
                       'AY572034': 'civet007', 'AY572037': 'civet019', 'AY613950': 'PC4-227'}
    keep_other = {'FJ882957': 'Mice MA15', 'NC_014470': 'Bat BM48-31/BGR/2008', 'MN996532': 'Bat RaTG13',
               'KF367457': 'Bat WIV1', 'DQ022305': 'HKU3-1', 'GQ153544': 'HKU3-9',
               'KC881005': 'Bat RsSHC014', 'MT799521': 'Pangolin cDNA8-S',
               'JX163927': 'Ferret Tor2/FP1-10851', 'HQ890538': 'Mouse MA15 ExoN1 isolate d2om5',
               'HQ890534': 'Mouse MA15 ExoN1 isolate d2om1',
               'FJ882931': 'Vero cells ExoN1 isolate P3pp12',
               'FJ882951': 'Vero cells ExoN1 isolate P3pp3'}
    
    #===== Below this line for documentation purposes only: not enforced to keep ====
    # Mutated sequences from Urbani (human) to become lethal for BALB/c mice
    seqs_MA15_vero_mut = {'FJ882942': 'BALB/c mice', 'FJ882951': 'BALB/c mice'}
    # Vero cells labeled human host
    seqs_human_vero = {'KF514407': 'ExoN1 c5.7P20', 'JX162087': 'ExoN1 c5P10',
                        'FJ882960': 'ExoN1 P3pp34', 'FJ882956': 'ExoN1 P3pp53',
                        'FJ882950': 'ExoN1 P3pp60', 'FJ882949': 'wtic-MB P3pp23',
                        'FJ882944': 'ExoN1 P3pp23', 'FJ882941': 'ExoN1 P3pp8',
                        'FJ882940': 'ExoN1 P3pp37', 'FJ882929': 'ExoN1 P3pp1',
                        'FJ882931': 'ExoN1 P3pp12'}
    # Other sequences of interest. Duplicates with different host. 
    do_not_include = {'MT308984': 'host=mouse, mut from bat KC881005'}
    
    return keep_seqs_human | keep_seqs_civet | keep_other | keep_refseqs_human_vero


def remove_selected_duplicates(acc_list, df):
    """This functions will go through a list of accession numbers and remove all records in the
    dataframe matching the exact same sequence

    Parameters
    ----------
    acc_list: List
        List of accession numbers to match. This one will stay in the final dataframe
    df: DataFrame
        Dataframe from were duplicates will be removed
    
    Return
    ------
    final: DataFrame
        Final dataframe with no duplicates for records from the accession list
    """
    df_copy = df.copy()
    print('Original dataset: {} sequences'.format(df.shape[0]))
    keep = df[df.Accession.isin(acc_list)]
    for s in acc_list:
        dups = df_copy[df_copy.Sequence.isin(df_copy.Sequence[df_copy.Accession == s])].index
        df_copy.drop(dups, inplace = True)
    final = pd.concat([df_copy, keep], axis=0)
    print('After removing duplicates: {} sequences'.format(final.shape[0]))
    return final


if __name__=='__main__':
    """ Function create_SARS1_dataset() is reading files from:
          1) hosts_path ='./sarbecovirus_hosts_993959.csv' with host info
          2) SARS1_spike = '../SARS1/sarbe_noCov2/sars1_spike_26jul21_560.fasta' with sequences
          3) SARS1_host = '../SARS1/sarbe_noCov2/sars1_spike_560_26jul21_host_manual.csv' curated hosts

        Sanity check for extract_spike function
          1) All sequences from kan paper were found.
          2) All 78 of SARS1 records from original Hostexp DS are present in our new DS with spike
          3) SARS1 (noCov-2): Found 560 of 1544 (36% found).
               Found: len(seq_prot) >=1000 except for 9 (79, 534, 568, 577, 577, 839, 885, 926, 930)
               Not found: 88 have 534 or more residues, the rest are smaller (896)

        Kan dataset
        + 23 unique records used as selected accessions for final dataset
        
        Duplicates
        Selected sequences were used as a reference to remove duplicates. They have been documented on 
        selected_sequences()

        Vero cells
        Urbani and other 15 vero cells were assigned to their corresponding host after removing
        duplicates. Some of them include MA15 from BALB/c mice, but most were Homo sapiens. They
        have been documented on selected_sequences()

        Patents
        After removing duplicates, only 5 unique patents remained. Given that there was no information
        related to the host, they were removed from the final dataset
    """
    
    SARS1 = create_SARS1_dataset()
    sars1_path = '../SARS1/sars1_552_final_dups.csv'
    SARS1.to_csv(sars1_path,index=False)
    keep_seqs_sars1 = selected_sequences()
    nodups = remove_selected_duplicates(keep_seqs_sars1,SARS1)
    print('Patents in final dataset to be removed:')
    nodups.loc[nodups.Host == 'Patent human', ['Accession', 'Virus', 'Host']]
    nodups = nodups[nodups.Host != 'Patent human']
    print('Total after removing patents: ', nodups.shape[0])
    sars1_nodups_path = '../SARS1/sars1_195_final.csv'
    nodups.to_csv(sars1_nodups_path,index=False)
    nodups.Host_agg.value_counts()
    