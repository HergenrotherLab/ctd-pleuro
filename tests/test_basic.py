from .context import calc_props

import unittest
from rdkit import Chem


class PropTestSuite(unittest.TestCase):
    """Basic test cases using ampiciliin and benzyl penicillin"""

    def setUp(self):
        self.amp = Chem.MolFromSmiles(
            'O=C(O)[C@@H]2N3C(=O)[C@@H](NC(=O)[C@@H](c1ccccc1)N)[C@H]3SC2(C)C')
        self.amp.SetProp('_Name', 'amp')
        self.pen = Chem.MolFromSmiles(
            'CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)Cc3ccccc3)C(=O)O)C')


    def test_absolute_truth_and_meaning(self):
        assert True


    def test_builtin_props(self):
        amp = calc_props.calc_builtin_props(self.amp)
        props = amp.GetPropNames()
        self.assertEqual(len(props), 3)
        self.assertAlmostEqual(float(amp.GetProp('FractionCSP3')), 0.4375)
        self.assertAlmostEqual(float(amp.GetProp('MolWt')), 349.412)
        self.assertEqual(int(amp.GetProp('RingCount')), 3)


    def test_csv_header(self):
        amp = calc_props.calc_builtin_props(self.amp)
        header = calc_props.csv_header(amp)
        self.assertEqual(header, 'name,smiles,FractionCSP3,MolWt,RingCount')


    def test_mol_to_csv(self):
        amp = calc_props.calc_builtin_props(self.amp)
        line = calc_props.mol_props_to_csv(amp)
        self.assertEqual(line, 'amp,CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccccc3)C(=O)N2[C@H]1C(=O)O,0.4375,349.4120000000001,3')


    def test_ring_props(self):
        amp = self.amp
        calc_props.calcRingDescriptors(amp)
        props = amp.GetPropNames()
        self.assertEqual(len(props), 9)


if __name__ == '__main__':
    unittest.main()
