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


def _find_padded_rows(table: pd.DataFrame, cols: [str]) -> pd.DataFrame:
    _cols_check(table, cols)

    _filter = None
    for col in cols:
        err = table[col].apply(lambda x: x != x.strip())
        if _filter is None:
            _filter = err
        else:
            _filter |= err

    # enclose names in "" to indicate whitespaces
    errors = table[_filter].copy()
    for col in cols:
        errors[col] = errors[col].apply(lambda x: '"%s"' % x)

    return errors


def _find_multiword_rows(table: pd.DataFrame, cols: [str]) -> pd.DataFrame:
    _cols_check(table, cols)

    _filter = None
    for col in cols:
        err = table[col].apply(lambda x: len(x.split()) > 1)
        if _filter is None:
            _filter = err
        else:
            _filter |= err

    return table[_filter]


class QCTests(TestCase):
    def setUp(self):
        self.fp_guitar = 'traitator/tests/data/trait_data_n10906.csv'
        self.fp_tresor = 'data/traitdata_tresor.csv'
        self.cols_binomial = ['Genus', 'Species']

    def tearDown(self):
        pass

    def test_qctresor(self):
        tresor = pd.read_csv(self.fp_tresor, sep="\t")

        # test dimension of data
        self.assertEqual(tresor.shape, (370479, 5))

        # ensure all Genus/Species are not np.nan
        obs = _find_nan_rows(tresor, self.cols_binomial)
        self.assertTrue(
            obs.shape[0] == 0,
            msg="all Genus and Species names should be non-emptry:\n%s" % obs)

        # ensure all Genus/Species names have no leading or trailing
        # whitespaces
        obs = _find_padded_rows(tresor, self.cols_binomial)
        self.assertTrue(
            obs.shape[0] == 0,
            msg="you have whitespaces in your names:\n%s" % obs)

        # ensure Genus/Species names are exactly one word
        obs = _find_multiword_rows(tresor, self.cols_binomial)
        self.assertTrue(
            obs.shape[0] == 0,
            msg="Genus or Species consists of more than one word!:\n%s" % obs)

        self.assertTrue(
            any(map(lambda x: 'GW2011' in x,
                    tresor[
                        tresor['Genus'] == 'archaeon']['Species'].unique())),
            "Species name of Genus 'archaeon' should contain GW2011")

        # test critical name resolution
        for (genus, species) in [
                ('Anaerosalibacter', 'Anaerosalibacter sp.'),
                ('Aphanizomenon', 'flos-aquae'),
                ('Bacillus', 'safensis'),
                ('Brevundimonas', 'Brevundimonas sp.'),
                ('Dolichospermum', 'flos-aquae'),
                ('Finegoldia', 'Finegoldia sp.'),
                ('Halomonas', 'denitrificans'),
                ('Methanothermobacter', 'Methanothermobacter sp.'),
                ('Mucilaginibacter rigui', 'rigui'),
                ('Mycobacterium gordonae', 'paragordonae'),
                ('Oscillatoria', 'nigro-viridis'),
                ('Paulownia', 'witches-broom'),
                ('Plasticicumulans', 'lactativoran'),
                ('Pseudoclavibacter', 'Pseudoclavibacter sp.'),
                ('Ruania', 'albidiflava'),
                ('Selenomonas', 'Selenomonas sp.'),
                ('Sphingobacterium', 'composti'),
                ('Thalassospira', 'A40-3'),
                ('Thermovibrio', 'ammonificans'),
                ('archaeon', 'GW2011_AR10'),
                ('archaeon', 'GW2011_AR20'),
                ('haloarchaeon', '3A1-DGR'),
                ('olei', 'IMMIBHF-1T')]:
            self.AssertTrue(
                tresor[(tresor['Genus'] == genus) &
                       (tresor['Species'] == species)].shape[0] > 0)


if __name__ == '__main__':
    main()
