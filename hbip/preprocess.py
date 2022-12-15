import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
import gensim
from gensim.models import Word2Vec
import os
import pickle


def save_model_here(model_name, model_object, object_name):
    directory = './models/' + model_name
    if not os.path.isdir(directory):
        os.makedirs(directory)
    obj_path = directory + '/' + object_name + '.pkl'
    with open(obj_path, 'wb') as f:
        pickle.dump(model_object, f)


def train_test(df, strat='Species_agg', test_frac=0.15):
    """Split into train and test after shuffling stratifying by Virus

    Parameters
    ----------
    df :
    strat : str
        Name of the feature used for stratification
    test_frac : float
        Fraction to create the test set
    Returns
    -------
    train, test : DataFrame    
    """
    random_indices = np.random.RandomState(seed=42).permutation(df.shape[0])
    print('Shuffling data with seed')
    print('Seed working OK if indices are always the same: ', random_indices[:3])
    df_shuffled_raw = df.iloc[random_indices].copy()
    print('Test set fraction: {}%'.format(test_frac*100))
    train, test = train_test_split(df_shuffled_raw, test_size=test_frac, 
                                   stratify=df_shuffled_raw[strat], random_state=42)
    train = train.reset_index(drop=True)
    test = test.reset_index(drop=True)
    print('\nTrain data: ', train.shape[0])
    print(train[strat].value_counts())
    print('\nTest data: ', test.shape[0])
    print(test[strat].value_counts())
    return train, test


def continuous_trimer(seqs):
    """Creates all possible trimers for a sequence padding the beginning and end
    (e.g. 'ABCD' = '--A, -AB, ABC, BCD, CD-, D--'

    Parameters
    ----------
    seqs: Series

    Returns
    -------
    list
     """
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
        cont3_seq = continuous_trimer_one(s)
        trigram.append(cont3_seq)
    return trigram


def trigram_3slides(df, label='Human'):
    """ Creates trigrams for all sequences. It triplicates the total number of sequences
        (data-augmentation) by sliding the word by one character.
        Parameters
       ----------
       df : DataFrame
           Dataframe with: Sequence, Accession and Human columns

       Returns
       -------
       trigram : list
           List with the trigrams for all sequences
       target : DataFrame
           A dataframe with details about the target variable after data augmentation (triplicates)
           The new column "Slide" will identify the 3 possible trigrams (slide1, slide2, slide3)
       """
       
    def trigram_3slides_oneseq(seq):
        """Produces the 3 possible (slides) trigrams from one sequence"""
        n = len(seq)
        word_lists = [[], [seq[0]], [seq[0:2]]]
        for aa in range(n):
            word_lists[aa % 3].append(seq[aa:aa + 3])
        return word_lists
    trigram = []
    accession = []
    slide = []
    lab = []
    for i in range(df.shape[0]):
        trigram_seq = trigram_3slides_oneseq(df.iloc[i].Sequence)
        acc_triplicate = [df.iloc[i].Accession]*3
        label_triplicate = [df.iloc[i][label]] * 3
        trigram.extend(trigram_seq)
        accession.extend(acc_triplicate)
        lab.extend(label_triplicate)
        slide.extend(['slide1', 'slide2', 'slide3'])
    df_target = pd.DataFrame({'Accession': accession, 'Slide': slide, label: lab})
    return trigram, df_target


def create_word_vectors(_wordsof3, _vector_size=100, _context_size=25):
    """ Convert trigram sequences into vectors using word2vec skip-gram with negative sampling

    Parameters
    ----------
    _wordsof3 : series
        Series with trigrams from each sequence
    _vector_size : int
        Desired final dimensionality of the word vectors
    _context_size : int
        Context size, that is number of words used to left and right of current word

    Returns
    -------
    _word_vectors : numpy.ndarray
        Values from word vectors (embedding matrix)
    _word2id : dict
        Dictionary of available words in the embedding matrix with correspondent index in _word_vectors
    _unknown_vector : numpy.ndarray
        Mean vector to impute new/unknown words

    """
    _model = gensim.models.Word2Vec(_wordsof3, min_count=1, vector_size=_vector_size,
                                    window=_context_size, seed=42, negative=1, sg=1, workers=1)
    _word_vectors = _model.wv.vectors
    _word2id = dict(_model.wv.key_to_index)
    _unknown_vector = _word_vectors.mean(0)
    return _word_vectors, _word2id, _unknown_vector


def create_seq_embedding(_trigram, _embedding_matrix, _word2id, _unknown_vector, _word_vector_size=100):
    """ A single embedding vector for the sequence with the sum of all word vectors
    Parameters
    ----------
    _trigram : list
        trigram from sequence
    _embedding_matrix : list
        Matrix with word vectors
    _word2id : dict
        Dictionary with vocabulary and index to corresponding vector
    _unknown_vector :
        Vector used to input unknown words
    _word_vector_size: int
        Desired size for embedding
    
    Returns
    -------
    numpy.ndarray
    """

    print('\nCreating sequence embeddings... \n')
    _seq_embedding_list = []
    _count_new_words = 0
    for seq in _trigram:
        _seq_embedding = [0] * _word_vector_size
        for wof3 in seq:
            if wof3 in _word2id:
                _seq_embedding += _embedding_matrix[_word2id.get(wof3)]
            else:
                _count_new_words += 1
                _seq_embedding += _unknown_vector
        _seq_embedding_list.append(_seq_embedding)
    print('There were :', _count_new_words, 'new words in the vocabulary')
    return np.array(_seq_embedding_list)


def embedding(embed_name, train, test, label='Human', save_embed=False):
    """Creates trigrams from the training set and feeds them to a skip-gram model with
    negative sampling (word2vec) to create the sequence embedding. The test set uses
    the resulting model to create an embedding. Unknown words are replaced with a mean vector

    Parameters
    ----------
    embed_name: str
        Name for the model. It will be used to save the binary file with the embedding
    train, test : DataFrame
        Dataframes for the train and test set with columns [ 'Sequence', 'Accession', label]
    label: str
        Name of target variable
    save_embed: Boolean
        It will save the model as binary when True
    Returns
    --------
    train_triple, test_triple: np.ndarray
        Embedding for the Sequence after data augmentation (triplicate)
    train_target, test_target : DataFrame
        Details about the target and sequence after data augmentation: ['Accession', 'Slide', 'Source', label].
        Source will be 'Train' for all at the train set and 'Test' at the test set
    """
    sars1_3w_train, train_target = trigram_3slides(train, label=label)
    train_target['Source'] = ['Train'] * (train_target.shape[0])
    sars1_3w_test, test_target = trigram_3slides(test, label=label)
    test_target['Source'] = ['Test'] * (test_target.shape[0])
    mat_vect, vocab, mean_vec = create_word_vectors(sars1_3w_train)
    print('Training set:')
    train_triple = create_seq_embedding(sars1_3w_train, mat_vect, vocab, mean_vec)
    if save_embed:
        save_model_here(embed_name, mat_vect, 'mat_vect')
        save_model_here(embed_name, vocab, 'vocab')
        save_model_here(embed_name, mean_vec, 'mean_vec')
    print('\nTest set:')
    test_triple = create_seq_embedding(sars1_3w_test, mat_vect, vocab, mean_vec)
    return train_triple, train_target, test_triple, test_target


def embedding_continuous(embed_name, train_df, test_df, label='Human', save_embed=False):
    """Creates trigrams from the training set and feeds them to a skip-gram model with
    negative sampling (word2vec) to create the sequence embedding. The test set uses
    the resulting model to create an embedding. Unknown words are replaced with a mean vector

    Parameters
    ----------
    embed_name: str
        Name for the model. It will be used to save the binary file with the embedding
    train_df, test_df : DataFrame
        Dataframes for the train and test set with columns [ 'Sequence', 'Accession', label]
    label: string
        Name of the target variable
    save_embed: Boolean
        It will save the model as binary when True
    Returns
    --------
    train_embedding, test_embedding: np.ndarray
        Embedding for the Sequence
    train_target, test_target : DataFrame
        Details about the target and sequence: ['Accession', 'Source', label].
        Source will be 'Train' for all at the train set and 'Test' at the test set
    """
    df_3w_train = continuous_trimer(train_df.Sequence)
    train_target = train_df.loc[:, ['Accession', label]].copy()
    train_target['Source'] = ['Train'] * (train_target.shape[0])
    df_3w_test = continuous_trimer(test_df.Sequence)
    test_target = test_df.loc[:, ['Accession', label]].copy()
    test_target['Source'] = ['Test'] * (test_target.shape[0])
    mat_vect, vocab, mean_vec = create_word_vectors(df_3w_train)
    print('Training set:')
    train_embedding = create_seq_embedding(df_3w_train, mat_vect, vocab, mean_vec)
    if save_embed:
        save_model_here(embed_name, mat_vect, 'mat_vect')
        save_model_here(embed_name, vocab, 'vocab')
        save_model_here(embed_name, mean_vec, 'mean_vec')
    print('\nTest set:')
    test_embedding = create_seq_embedding(df_3w_test, mat_vect, vocab, mean_vec)
    return train_embedding, train_target, test_embedding, test_target


def single_embedding(model_name, sequence):
    """
    Generate the embedding for a single sequence from the trained model from word2vec

    :param model_name: str
    :param sequence: Series
    :return: ndarray
    """
    directory = './models/' + model_name
    obj_path = directory + '/mat_vect.pkl'
    with open(obj_path, 'rb') as f:
        mat_vect = pickle.load(f)
    obj_path = directory + '/vocab.pkl'
    with open(obj_path, 'rb') as f:
        vocab = pickle.load(f)
    obj_path = directory + '/mean_vec.pkl'
    with open(obj_path, 'rb') as f:
        mean_vec = pickle.load(f)
    seq_3w = continuous_trimer(sequence)
    return create_seq_embedding(seq_3w, mat_vect, vocab, mean_vec)


# if __name__ == '__main__':
#     Urbani = ['MFIFLLFLTLTSGSDLDRCTTFDDVQAPNYTQHTSSMRGVYYPDEIFRSDTLYLTQDLFLPFYSNVTGFHTINHTFGNPVIPFKDGIYFAATEKSNVV\
#     RGWVFGSTMNNKSQSVIIINNSTNVVIRACNFELCDNPFFAVSKPMGTQTHTMIFDNAFNCTFEYISDAFSLDVSEKSGNFKHLREFVFKNKDGFLYVYKGYQPIDVVR\
#     DLPSGFNTLKPIFKLPLGINITNFRAILTAFSPAQDIWGTSAAAYFVGYLKPTTFMLKYDENGTITDAVDCSQNPLAELKCSVKSFEIDKGIYQTSNFRVVPSGDVVRF\
#     PNITNLCPFGEVFNATKFPSVYAWERKKISNCVADYSVLYNSTFFSTFKCYGVSATKLNDLCFSNVYADSFVVKGDDVRQIAPGQTGVIADYNYKLPDDFMGCVLAWNT\
#     RNIDATSTGNYNYKYRYLRHGKLRPFERDISNVPFSPDGKPCTPPALNCYWPLNDYGFYTTTGIGYQPYRVVVLSFELLNAPATVCGPKLSTDLIKNQCVNFNFNGLTG\
#     TGVLTPSSKRFQPFQQFGRDVSDFTDSVRDPKTSEILDISPCSFGGVSVITPGTNASSEVAVLYQDVNCTDVSTAIHADQLTPAWRIYSTGNNVFQTQAGCLIGAEHVD\
#     TSYECDIPIGAGICASYHTVSLLRSTSQKSIVAYTMSLGADSSIAYSNNTIAIPTNFSISITTEVMPVSMAKTSVDCNMYICGDSTECANLLLQYGSFCTQLNRALSGI\
#     AAEQDRNTREVFAQVKQMYKTPTLKYFGGFNFSQILPDPLKPTKRSFIEDLLFNKVTLADAGFMKQYGECLGDINARDLICAQKFNGLTVLPPLLTDDMIAAYTAALVS\
#     GTATAGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKQIANQFNKAISQIQESLTTTSTALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVE\
#     AEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQAAPHGVVFLHVTYVPSQERNFTTAPAICHEGKAYFPREGVFVF\
#     NGTSWFITQRNFFSPQIITTDNTFVSGNCDVVIGIINNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKY\
#     EQYIKWPWYVWLGFIAGLIAIVMVTILLCCMTSCCSCLKGACSCGSCCKFDEDDSEPVLKGVKLHYT']
#
#     print(single_embedding('ab_no_sars2', Urbani))
