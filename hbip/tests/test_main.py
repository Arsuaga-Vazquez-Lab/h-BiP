import unittest
import numpy as np
from numpy.linalg import norm
import os
import pickle

import hbip.predict
import hbip_reproduce as mn
from hbip import preprocess as pp


def base_model_path(name):
    return os.path.join("./models", name)


def model_files(name):
    return [
        os.path.join(base_model_path(name), filename)
        for filename in "LR85.pkl	mat_vect.pkl	mean_vec.pkl	scale.pkl	vocab.pkl".split()
    ]


def unlink_model(name):
    path = base_model_path(name)
    if not os.path.isdir(path):
        return
    for filename in model_files(name):
        os.unlink(filename)
    os.rmdir(path)


class TestMain(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        unlink_model("rbd_189")
        unlink_model("sars1_195")

    def test_rbd_189_binds(self):
        embedding, final_scores = mn.main("./hbip/tests/testdata/rbd_189_config.yml")

        # verify that embeddings are as expected
        expected = np.load("./hbip/tests/fixture/rbd_189_embedding.npy")
        actual = embedding
        self.assertAlmostEqual(0, norm(expected - actual))

        # verify that highest probability has expected index
        probs = final_scores.Binds_prob.values
        expected = 183
        actual = np.argmax(probs)
        self.assertAlmostEqual(0, norm(expected - actual))

        # verify that probability distribution is as expected
        expected = np.load("./hbip/tests/fixture/rbd_189_final_scores_binds_prob.npy")
        actual = probs
        self.assertAlmostEqual(0, norm(expected - actual))

    def test_sars1_195(self):
        embedding, final_scores = mn.main("./hbip/tests/testdata/sars1_195_config.yml")

        # verify that embeddings are as expected
        expected = np.load("./hbip/tests/fixture/sars1_195_embedding.npy")
        actual = embedding
        self.assertAlmostEqual(0, norm(expected - actual))

        # verify that highest probability has expected index
        probs = final_scores.Human_prob.values
        expected = 53
        actual = np.argmax(probs)
        self.assertAlmostEqual(0, norm(expected - actual))

        # verify that probability distribution is as expected
        expected = np.load("./hbip/tests/fixture/sars1_195_final_scores_binds_prob.npy")
        actual = probs
        self.assertAlmostEqual(0, norm(expected - actual), places=6)

    def test_model_save_locations(self):
        # verify that model is saved with expected file names in expected directory
        mn.main("./hbip/tests/testdata/rbd_189_config.yml", save_model=True)
        base_path = base_model_path("rbd_189")
        filesnames = model_files("rbd_189")
        self.assertTrue(os.path.isdir(base_path))
        for filename in filesnames:
            self.assertTrue(os.path.exists(filename), msg=filename)

    def test_restore_model(self):
        # train and save model
        mn.main("./hbip/tests/testdata/rbd_189_config.yml", save_model=True)

        # reload model
        with open("./models/rbd_189/LR85.pkl", "rb") as f:
            LR85 = pickle.load(f)
        seq = np.load("./hbip/tests/fixture/rbd_189_Xall_norm.npy")
        actual = hbip.predict.predict(LR85, seq)
        expected = np.load("./hbip/tests/fixture/rbd_189_final_scores_binds_prob.npy")
        self.assertAlmostEqual(0, norm(expected - actual))

    def test_single_embedding(self):
        # train and save model
        mn.main("./hbip/tests/testdata/sars1_195_config.yml", save_model=True)

        # load sequences
        filename = "./hbip/tests/testdata/virus_predict.fasta"
        seq = hbip.predict.read_seqs(filename).iloc[:, 1]

        actual = pp.single_embedding("sars1_195", seq)
        expected = np.load("./hbip/tests/fixture/sars1_195_single_embedding.npy")
        self.assertAlmostEqual(0, norm(expected - actual))

    def test_predict_new(self):
        # train and save model
        mn.main("./hbip/tests/testdata/sars1_195_config.yml", save_model=True)

        # load sequences
        filename = "./hbip/tests/testdata/virus_predict.fasta"

        seq = hbip.predict.read_seqs(filename).iloc[:, 1]

        actual = hbip.predict.predict_new("sars1_195", seq, models_dir="./models/")
        expected = np.array(
            [
                0.002782201706697373,
                7.63459742402316e-09,
                4.858718049397064e-05,
                0.22100315897583508,
                9.807818834934887e-08,
            ]
        )
        self.assertAlmostEqual(0, norm(expected - actual))
