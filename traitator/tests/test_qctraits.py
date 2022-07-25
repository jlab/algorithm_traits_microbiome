from unittest import TestCase, main
import pandas as pd


def _cols_check(table: pd.DataFrame, cols: [str]):
    assert all([c in table.columns for c in cols]), \
        "some columns in 'cols' are not in your table!"


def _find_nan_rows(table: pd.DataFrame, cols: [str]) -> pd.DataFrame:
    _cols_check(table, cols)

    errors = []
    for col in cols:
        errors.append(table[~pd.notnull(table[col])])

    return pd.concat(errors)


class QCTests(TestCase):
    def setUp(self):
        self.fp_guitar = 'traitator/tests/data/trait_data_n10906.csv'
        self.fp_tresor = 'data/Traits.csv'
        self.cols_binomial = ['Genus', 'Species']

    def tearDown(self):
        pass

    def test_qctresor(self):
        tresor = pd.read_csv(self.fp_tresor, sep=",")

        # test dimension of data
        self.assertEqual(tresor.shape, (370479, 5))

        obs = _find_nan_rows(tresor, self.cols_binomial)
        self.assertTrue(
            obs.shape[0] == 0,
            msg="all Genus and Species names should be non-emptry:\n%s" % obs)


if __name__ == '__main__':
    main()
