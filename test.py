#!/usr/bin/env python
import unittest

from hbip.tests.test_example import TestExample
from hbip.tests.test_load_test_data import TestLoadTestData
from hbip.tests.test_main import TestMain


class CountSuite(object):
    def __init__(self):
        self.count = 0
        self.s = unittest.TestSuite()

    def add(self, tests):
        self.count += 1
        print("%d: %s" % (self.count, tests.__name__))
        self.s.addTest(unittest.makeSuite(tests))


def suite():
    s = CountSuite()

    s.add(TestExample)
    s.add(TestLoadTestData)
    s.add(TestMain)

    return s.s


if __name__ == "__main__":
    runner = unittest.TextTestRunner()
    runner.run(suite())
