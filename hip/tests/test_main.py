import unittest
import numpy as np
from numpy.linalg import norm

import main as mn


class TestMain(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_rbd_189_binds(self):
        embedding, final_scores = mn.main("./hip/tests/testdata/rbd_189_config.yml")

        # verify that embeddings are as expected
        expected = np.load("./hip/tests/fixture/rbd_189_embedding.npy")
        actual = embedding
        self.assertAlmostEqual(0, norm(expected - actual))

        # verify that highest probability has expected index
        probs = final_scores.Binds_prob.values
        expected = 183
        actual = np.argmax(probs)
        self.assertAlmostEqual(0, norm(expected - actual))

        # verify that probability distribution is as expected
        expected = np.load("./hip/tests/fixture/rbd_189_final_scores_binds_prob.npy")
        actual = probs
        self.assertAlmostEqual(0, norm(expected - actual))
    #
    # def test_model_save_locations(self):
    #     self.fail("incomplete")
    #
    # def test_restore_model(self):
    #     self.fail("incomplete")
