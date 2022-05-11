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
        expected = np.load("./hip/tests/fixture/rbd_189_embedding.npy")
        actual = embedding
        self.assertAlmostEqual(0, norm(expected - actual))

        expected = np.load("./hip/tests/fixture/rbd_189_final_scores_binds_prob.npy")
        actual = final_scores.Binds_prob.values
        self.assertAlmostEqual(0, norm(expected - actual))
