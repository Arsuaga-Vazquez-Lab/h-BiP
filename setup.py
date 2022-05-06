import sys
from setuptools import setup, find_packages


def forbid_publish():
    argv = sys.argv
    blacklist = ["register", "upload"]

    for command in blacklist:
        if command in argv:
            values = {"command": command}
            print('Command "%(command)s" has been blacklisted, exiting...' % values)
            sys.exit(2)


forbid_publish()

setup(
    name="hip",
    version="0.0.1",
    author="Georgina Gonzalez",
    install_requires=[],
    packages=find_packages(),
    include_package_data=True,
)
