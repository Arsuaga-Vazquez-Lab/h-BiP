import unittest

import yaml
import pandas as pd

from hbip.predict import read_seqs


def load_data(base_path, data_name):
    config_path = base_path + data_name + "_config.yml"
    with open(config_path, "r") as file:
        config = yaml.safe_load(file)

    data_path = config["path"]
    return pd.read_csv(data_path)


class TestLoadTestData(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_load_rbd_189_binds(self):
        data = load_data("./hbip/tests/testdata/", "rbd_189")

        expected = (189, 8)
        actual = data.shape
        self.assertEqual(expected, actual)

        expected = [
            "Accession",
            "Species",
            "Virus",
            "Note",
            "Host_agg",
            "Human",
            "Sequence",
            "Binds",
        ]
        actual = list(data.columns)
        self.assertEqual(expected, actual)

        expected = {
            "Accession": "MT121216",
            "Species": "Pangolin coronavirus",
            "Virus": "Pangolin coronavirus isolate MP789",
            "Host_agg": "Pangolin",
            "Human": 0,
            "Sequence": "ITNLCPFGEVFNATTFASVYAWNRKRISNCVADYSVLYNSTSFSTFKCYGVSPTKLNDLCFTNVYADSFVVRGDEVRQIAPGQTGRIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFHPTNGVGYQPYRVVVLSFELLKAPATV",
            "Binds": 1.0,
        }
        actual = data.iloc[0].to_dict()
        del actual["Note"]
        for key in actual:
            self.assertEqual(expected[key], actual[key], msg=key)
        self.assertEqual(expected, actual)

    def test_load_sars1_195(self):
        data = load_data("./hbip/tests/testdata/", "sars1_195")

        expected = (195, 8)
        actual = data.shape
        self.assertEqual(expected, actual)

        expected = [
            "Accession",
            "Sequence",
            "Species",
            "Virus",
            "Note",
            "Host",
            "Host_agg",
            "Human",
        ]
        actual = list(data.columns)
        self.assertEqual(expected, actual)

        expected = {
            "Accession": "MT121216",
            "Sequence": "MLFFFFLHFALVNSQCVNLTGRAAIQPSFTNSSQRGVYYPDTIFRSNTLVLSQGYFLPFYSNVSWYYALTKTNSAEKRVDNPVLDFKDGIYFAATEKSNIVRGWIFGTTLDNTSQSLLIVNNATNVIIKVCNFQFCYDPYLSGYYHNNKTWSTREFAVYSSYANCTFEYVSKSFMLDIAGKSGLFDTLREFVFRNVDGYFKIYSKYTPVNVNSNLPIGFSALEPLVEIPAGINITKFRTLLTIHRGDPMPNNGWTVFSAAYYVGYLAPRTFMLNYNENGTITDAVDCALDPLSEAKCTLKSLTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATTFASVYAWNRKRISNCVADYSVLYNSTSFSTFKCYGVSPTKLNDLCFTNVYADSFVVRGDEVRQIAPGQTGRIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFHPTNGVGYQPYRVVVLSFELLKAPATVCGPKQSTNLVKNKCVNFNFNGLTGTGVLTESSKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNTYECDIPIGAGICASYQTQTNSRSVSSQAIIAYTMSLGAENSVAYANNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSIECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPSQEKNFTTTPAICHEGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGSCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIIMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT",
            "Species": "Pangolin coronavirus",
            "Virus": "Pangolin coronavirus isolate MP789",
            "Host": "Manis javanica",
            "Host_agg": "Pangolin",
            "Human": 0,
        }
        actual = data.iloc[0].to_dict()
        del actual["Note"]
        for key in actual:
            self.assertEqual(expected[key], actual[key], msg=key)
        self.assertEqual(expected, actual)

    def test_read_seq(self):
        filename = "./hbip/tests/testdata/virus_predict.fasta"
        df = read_seqs(filename)
        expected = ["Description", "Sequence"]
        actual = list(df.columns)
        self.assertEqual(expected, actual)

        expected = (5, 2)
        actual = df.shape
        self.assertEqual(expected, actual)

        expected = {
            "Description": "MT799521,Pangolin coronavirus, Manis javanica, Pangolin coronavirus isolate cDNA8-S surface glycoprotein (S) gene",
            "Sequence": "MLFFFFLHFALVNSQCVNLTGRAAIQPSFTNSSQRGVYYPDTIFRSNTLVLSQGYFLPFYSNVSWYYALTKTNSAEKRVDNPVLDFKDGIYFAATEKSNIVRGWIFGTTLDNTSQSLLIVNNATNVIIKVCNFQFCYDPYLSGYYHNNKTWSTREFAVYSSYANCTFEYVSKSFMLDIAGKSGLFDTLREFVFRNVDGYFKIYSKYTPVNVNSNLPIGFSALEPLVEIPAGINITKFRTLLTIHRGDPMPNNGWTVFSAAYYVGYLAPRTFMLNYNENGTITDAVDCALDPLSEAKCTLKSLTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATTFASVYAWNRKRISNCVADYSVLYNSTSFSTFKCYGVSPTKLNDLCFTNVYADSFVVRGDEVRQIAPGQTGRIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFHPTNGVGYQPYRVVVLSFELLNAPATVCGPKQSTNLVKNKCVNFNFNGLTGTGVLTESSKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNTYECDIPIGAGICASYQTQTNSRSVSSQAIIAYTMSLGAENSVAYANNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSIECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPSQEKNFTTTPAICHEGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGSCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIIMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT",
        }
        actual = df.iloc[0].to_dict()
        self.assertEqual(expected, actual)
