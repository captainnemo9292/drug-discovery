import * as React from 'react';
import { Component, useState } from 'react';
import { DataTable } from 'react-native-paper';
import {
  TouchableOpacity,
  Alert,
  StyleSheet,
  Text,
  TextInput,
  View,
  FlatList,
  ScrollView,
  ActivityIndicator,
  Image,
} from 'react-native';

export default class MoleculeGeneration extends React.Component {


  state = {
    input_n: '',
    show_input: true,
    show_data: false,
    show_loading: false,
    data_arr: '',
  };

  async GenerateMolecule(input_n) {
    this.setState({ show_input: false });
    this.setState({ show_loading: true });
    try {
      let response = await fetch(
        'https://molecule-generation-dv4wcr3seq-an.a.run.app/generate',
        {
          method: 'POST',
          headers: {
            Accept: 'application/json',
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({
            n_sample: input_n,
          }),
        }
      );
      let json = await response.json();
      if (json['compound_data'] == 'failed') {
        this.setState({ show_loading: false });
        this.setState({ show_input: true });
        this.FailedAlert();
      } else {
        this.setState({
          data_arr: json['compound_data']['PropertyTable']['Properties'],
        });
        this.setState({ show_loading: false });
        this.setState({ show_data: true });
      }
    } catch (error) {
      console.error(error);
    }
  }

  Back() {
    this.setState({ show_loading: false });
    this.setState({ show_data: false });
    this.setState({ show_input: true });
  }

  DataProcessing(input_n) {
    var n_sample = Number(input_n);
    console.log(n_sample);
    this.GenerateMolecule(n_sample);
  }

  DataBase({ molecule_data, navigation }) {

    var image_uri =
      'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/' +
      molecule_data['InChIKey'] +
      '/PNG';
    console.log('data:');
    console.log(JSON.stringify(molecule_data));
    console.log(image_uri);
    return (
      <View>
      <DataTable>
        <DataTable.Header style={{ height: 80 }}>
          <DataTable.Title>Molecule</DataTable.Title>
        </DataTable.Header>

        <DataTable.Row>
          <Image
            style={{
              width: 300,
              height: 300,
              resizeMode: 'contain',
            }}
            source={{
              uri: image_uri,
            }}
          />
        </DataTable.Row>
        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert(
              'MolecularFormula',
              String(molecule_data['MolecularFormula'])
            )
          }>
          <DataTable.Cell>MolecularFormula</DataTable.Cell>
          <DataTable.Cell>
            {String(molecule_data['MolecularFormula'])}
          </DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert(
              'MolecularWeight',
              String(molecule_data['MolecularWeight'])
            )
          }>
          <DataTable.Cell>MolecularWeight</DataTable.Cell>
          <DataTable.Cell>
            {String(molecule_data['MolecularWeight'])}
          </DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert(
              'CanonicalSMILES',
              String(molecule_data['CanonicalSMILES'])
            )
          }>
          <DataTable.Cell>CanonicalSMILES</DataTable.Cell>
          <DataTable.Cell>
            {String(molecule_data['CanonicalSMILES'])}
          </DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert(
              'IsomericSMILES',
              String(molecule_data['IsomericSMILES'])
            )
          }>
          <DataTable.Cell>IsomericSMILES</DataTable.Cell>
          <DataTable.Cell>
            {String(molecule_data['IsomericSMILES'])}
          </DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() => Alert.alert('InChI', String(molecule_data['InChI']))}>
          <DataTable.Cell>InChI</DataTable.Cell>
          <DataTable.Cell>{String(molecule_data['InChI'])}</DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert('InChIKey', String(molecule_data['InChIKey']))
          }>
          <DataTable.Cell>InChIKey</DataTable.Cell>
          <DataTable.Cell>{String(molecule_data['InChIKey'])}</DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert('IUPACName', String(molecule_data['IUPACName']))
          }>
          <DataTable.Cell>IUPACName</DataTable.Cell>
          <DataTable.Cell>{String(molecule_data['IUPACName'])}</DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() => Alert.alert('XLogP', String(molecule_data['XLogP']))}>
          <DataTable.Cell>XLogP</DataTable.Cell>
          <DataTable.Cell>{String(molecule_data['XLogP'])}</DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert('ExactMass', String(molecule_data['ExactMass']))
          }>
          <DataTable.Cell>ExactMass</DataTable.Cell>
          <DataTable.Cell>{String(molecule_data['ExactMass'])}</DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert(
              'MonoisotopicMass',
              String(molecule_data['MonoisotopicMass'])
            )
          }>
          <DataTable.Cell>MonoisotopicMass</DataTable.Cell>
          <DataTable.Cell>
            {String(molecule_data['MonoisotopicMass'])}
          </DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() => Alert.alert('TPSA', String(molecule_data['TPSA']))}>
          <DataTable.Cell>TPSA</DataTable.Cell>
          <DataTable.Cell>{String(molecule_data['TPSA'])}</DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert('Complexity', String(molecule_data['Complexity']))
          }>
          <DataTable.Cell>Complexity</DataTable.Cell>
          <DataTable.Cell>{String(molecule_data['Complexity'])}</DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert('Charge', String(molecule_data['Charge']))
          }>
          <DataTable.Cell>Charge</DataTable.Cell>
          <DataTable.Cell>{String(molecule_data['Charge'])}</DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert(
              'HBondDonorCount',
              String(molecule_data['HBondDonorCount'])
            )
          }>
          <DataTable.Cell>HBondDonorCount</DataTable.Cell>
          <DataTable.Cell>
            {String(molecule_data['HBondDonorCount'])}
          </DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert(
              'HBondAcceptorCount',
              String(molecule_data['HBondAcceptorCount'])
            )
          }>
          <DataTable.Cell>HBondAcceptorCount</DataTable.Cell>
          <DataTable.Cell>
            {String(molecule_data['HBondAcceptorCount'])}
          </DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert(
              'RotatableBondCount',
              String(molecule_data['RotatableBondCount'])
            )
          }>
          <DataTable.Cell>RotatableBondCount</DataTable.Cell>
          <DataTable.Cell>
            {String(molecule_data['RotatableBondCount'])}
          </DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert(
              'HeavyAtomCount',
              String(molecule_data['HeavyAtomCount'])
            )
          }>
          <DataTable.Cell>HeavyAtomCount</DataTable.Cell>
          <DataTable.Cell>
            {String(molecule_data['HeavyAtomCount'])}
          </DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert(
              'IsotopeAtomCount',
              String(molecule_data['IsotopeAtomCount'])
            )
          }>
          <DataTable.Cell>IsotopeAtomCount</DataTable.Cell>
          <DataTable.Cell>
            {String(molecule_data['IsotopeAtomCount'])}
          </DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert(
              'AtomStereoCount',
              String(molecule_data['AtomStereoCount'])
            )
          }>
          <DataTable.Cell>AtomStereoCount</DataTable.Cell>
          <DataTable.Cell>
            {String(molecule_data['AtomStereoCount'])}
          </DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert(
              'DefinedAtomStereoCount',
              String(molecule_data['DefinedAtomStereoCount'])
            )
          }>
          <DataTable.Cell>DefinedAtomStereoCount</DataTable.Cell>
          <DataTable.Cell>
            {String(molecule_data['DefinedAtomStereoCount'])}
          </DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert(
              'UndefinedAtomStereoCount',
              String(molecule_data['UndefinedAtomStereoCount'])
            )
          }>
          <DataTable.Cell>UndefinedAtomStereoCount</DataTable.Cell>
          <DataTable.Cell>
            {String(molecule_data['UndefinedAtomStereoCount'])}
          </DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert(
              'BondStereoCount',
              String(molecule_data['BondStereoCount'])
            )
          }>
          <DataTable.Cell>BondStereoCount</DataTable.Cell>
          <DataTable.Cell>
            {String(molecule_data['BondStereoCount'])}
          </DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert(
              'DefinedBondStereoCount',
              String(molecule_data['DefinedBondStereoCount'])
            )
          }>
          <DataTable.Cell>DefinedBondStereoCount</DataTable.Cell>
          <DataTable.Cell>
            {String(molecule_data['DefinedBondStereoCount'])}
          </DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert(
              'UndefinedBondStereoCount',
              String(molecule_data['UndefinedBondStereoCount'])
            )
          }>
          <DataTable.Cell>UndefinedBondStereoCount</DataTable.Cell>
          <DataTable.Cell>
            {String(molecule_data['UndefinedBondStereoCount'])}
          </DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert(
              'CovalentUnitCount',
              String(molecule_data['CovalentUnitCount'])
            )
          }>
          <DataTable.Cell>CovalentUnitCount</DataTable.Cell>
          <DataTable.Cell>
            {String(molecule_data['CovalentUnitCount'])}
          </DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert('Volume3D', String(molecule_data['Volume3D']))
          }>
          <DataTable.Cell>Volume3D</DataTable.Cell>
          <DataTable.Cell>{String(molecule_data['Volume3D'])}</DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert(
              'XStericQuadrupole3D',
              String(molecule_data['XStericQuadrupole3D'])
            )
          }>
          <DataTable.Cell>XStericQuadrupole3D</DataTable.Cell>
          <DataTable.Cell>
            {String(molecule_data['XStericQuadrupole3D'])}
          </DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert(
              'YStericQuadrupole3D',
              String(molecule_data['YStericQuadrupole3D'])
            )
          }>
          <DataTable.Cell>YStericQuadrupole3D</DataTable.Cell>
          <DataTable.Cell>
            {String(molecule_data['YStericQuadrupole3D'])}
          </DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert(
              'ZStericQuadrupole3D',
              String(molecule_data['ZStericQuadrupole3D'])
            )
          }>
          <DataTable.Cell>ZStericQuadrupole3D</DataTable.Cell>
          <DataTable.Cell>
            {String(molecule_data['ZStericQuadrupole3D'])}
          </DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert(
              'FeatureCount3D',
              String(molecule_data['FeatureCount3D'])
            )
          }>
          <DataTable.Cell>FeatureCount3D</DataTable.Cell>
          <DataTable.Cell>
            {String(molecule_data['FeatureCount3D'])}
          </DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert(
              'FeatureAcceptorCount3D',
              String(molecule_data['FeatureAcceptorCount3D'])
            )
          }>
          <DataTable.Cell>FeatureAcceptorCount3D</DataTable.Cell>
          <DataTable.Cell>
            {String(molecule_data['FeatureAcceptorCount3D'])}
          </DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert(
              'FeatureDonorCount3D',
              String(molecule_data['FeatureDonorCount3D'])
            )
          }>
          <DataTable.Cell>FeatureDonorCount3D</DataTable.Cell>
          <DataTable.Cell>
            {String(molecule_data['FeatureDonorCount3D'])}
          </DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert(
              'FeatureAnionCount3D',
              String(molecule_data['FeatureAnionCount3D'])
            )
          }>
          <DataTable.Cell>FeatureAnionCount3D</DataTable.Cell>
          <DataTable.Cell>
            {String(molecule_data['FeatureAnionCount3D'])}
          </DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert(
              'FeatureCationCount3D',
              molecule_data['FeatureCationCount3D']
            )
          }>
          <DataTable.Cell>FeatureCationCount3D</DataTable.Cell>
          <DataTable.Cell>
            {String(molecule_data['FeatureCationCount3D'])}
          </DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert(
              'FeatureRingCount3D',
              String(molecule_data['FeatureRingCount3D'])
            )
          }>
          <DataTable.Cell>FeatureRingCount3D</DataTable.Cell>
          <DataTable.Cell>
            {String(molecule_data['FeatureRingCount3D'])}
          </DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert(
              'FeatureHydrophobeCount3D',
              String(molecule_data['FeatureHydrophobeCount3D'])
            )
          }>
          <DataTable.Cell>FeatureHydrophobeCount3D</DataTable.Cell>
          <DataTable.Cell>
            {String(molecule_data['FeatureHydrophobeCount3D'])}
          </DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert(
              'ConformerModelRMSD3D',
              String(molecule_data['ConformerModelRMSD3D'])
            )
          }>
          <DataTable.Cell>ConformerModelRMSD3D</DataTable.Cell>
          <DataTable.Cell>
            {String(molecule_data['ConformerModelRMSD3D'])}
          </DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert(
              'EffectiveRotorCount3D',
              String(molecule_data['EffectiveRotorCount3D'])
            )
          }>
          <DataTable.Cell>EffectiveRotorCount3D</DataTable.Cell>
          <DataTable.Cell>
            {String(molecule_data['EffectiveRotorCount3D'])}
          </DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert(
              'ConformerCount3D',
              String(molecule_data['ConformerCount3D'])
            )
          }>
          <DataTable.Cell>ConformerCount3D</DataTable.Cell>
          <DataTable.Cell>
            {String(molecule_data['ConformerCount3D'])}
          </DataTable.Cell>
        </DataTable.Row>

        <DataTable.Row
          style={{ height: 70 }}
          onPress={() =>
            Alert.alert('Fingerprint2D', String(molecule_data['Fingerprint2D']))
          }>
          <DataTable.Cell>Fingerprint2D</DataTable.Cell>
          <DataTable.Cell>
            {String(molecule_data['Fingerprint2D'])}
          </DataTable.Cell>
        </DataTable.Row>
      </DataTable>

      <TouchableOpacity
        style={styles.button}
        onPress={() => navigation.navigate('Virtual Screening',
        {
            data: molecule_data['CanonicalSMILES']
        })}>
        <Text style={styles.buttonText}>해당 데이터로 Virtual Screening 진행</Text>
      </TouchableOpacity>
      </View>
    );
  }

  FailedAlert() {
    Alert.alert(
      'Error',
      'API 통신 중 문제가 발생했습니다. 다시 시도해 주십시오...'
    );
  }

  render() {
    const { navigation } = this.props;
    return (
      <View style={styles.container}>
        {this.state.show_input && (
          <View>
            <TextInput
              style={[styles.input]}
              placeholder="생성할 신약 후보 물질의 개수를 입력하세요."
              onChangeText={input_n => this.setState({ input_n })}
            />
            <TouchableOpacity
              onPress={() => this.DataProcessing(this.state.input_n)}
              style={styles.button}>
              <Text style={styles.buttonText}>제출</Text>
            </TouchableOpacity>
          </View>
        )}
        {this.state.show_loading && (
          <View style={[styles.container]}>
            <ActivityIndicator size="large" color="#0DA900" />
          </View>
        )}
        {this.state.show_data && (
          <ScrollView>
            <View style={{ marginTop: 20 }} />
            <TouchableOpacity onPress={() => this.Back()} style={styles.button}>
              <Text style={styles.buttonText}>새로운 신약 후보 물질 생성</Text>
            </TouchableOpacity>
            <DataTable style={{ marginTop: 20 }}>
              <DataTable.Header style={{ height: 80 }}>
                <DataTable.Title>Molecule Data Format</DataTable.Title>
              </DataTable.Header>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert('MolecularFormula', 'Molecular formula')
                }>
                <DataTable.Cell>MolecularFormula</DataTable.Cell>
                <DataTable.Cell>Molecular formula</DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'MolecularWeight',
                    'The molecular weight is the sum of all atomic weights of the constituent atoms in a compound, measured in g/mol. In the absence of explicit isotope labelling, averaged natural abundance is assumed. If an atom bears an explicit isotope label, 100% isotopic purity is assumed at this location'
                  )
                }>
                <DataTable.Cell>MolecularWeight</DataTable.Cell>
                <DataTable.Cell>
                  The molecular weight is the sum of all atomic weights of the
                  constituent atoms in a compound, measured in g/mol. In the
                  absence of explicit isotope labelling, averaged natural
                  abundance is assumed. If an atom bears an explicit isotope
                  label, 100% isotopic purity is assumed at this location
                </DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'CanonicalSMILES',
                    'Canonical SMILES (Simplified Molecular Input Line Entry System) string.  It is a unique SMILES string of a compound, generated by a “canonicalization” algorithm'
                  )
                }>
                <DataTable.Cell>CanonicalSMILES</DataTable.Cell>
                <DataTable.Cell>
                  Canonical SMILES (Simplified Molecular Input Line Entry
                  System) string. It is a unique SMILES string of a compound,
                  generated by a “canonicalization” algorithm
                </DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'IsomericSMILES',
                    'Isomeric SMILES string.  It is a SMILES string with stereochemical and isotopic specifications'
                  )
                }>
                <DataTable.Cell>IsomericSMILES</DataTable.Cell>
                <DataTable.Cell>
                  Isomeric SMILES string. It is a SMILES string with
                  stereochemical and isotopic specifications
                </DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'InChI',
                    'Standard IUPAC International Chemical Identifier (InChI).  It does not allow for user selectable options in dealing with the stereochemistry and tautomer layers of the InChI string'
                  )
                }>
                <DataTable.Cell>InChI</DataTable.Cell>
                <DataTable.Cell>
                  Standard IUPAC International Chemical Identifier (InChI). It
                  does not allow for user selectable options in dealing with the
                  stereochemistry and tautomer layers of the InChI string
                </DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'InChIKey',
                    'Hashed version of the full standard InChI, consisting of 27 characters'
                  )
                }>
                <DataTable.Cell>InChIKey</DataTable.Cell>
                <DataTable.Cell>
                  Hashed version of the full standard InChI, consisting of 27
                  characters
                </DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'IUPACName',
                    'Chemical name systematically determined according to the IUPAC nomenclatures'
                  )
                }>
                <DataTable.Cell>IUPACName</DataTable.Cell>
                <DataTable.Cell>
                  Chemical name systematically determined according to the IUPAC
                  nomenclatures
                </DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'XLogP',
                    'Computationally generated octanol-water partition coefficient or distribution coefficient. XLogP is used as a measure of hydrophilicity or hydrophobicity of a molecule'
                  )
                }>
                <DataTable.Cell>XLogP</DataTable.Cell>
                <DataTable.Cell>
                  Computationally generated octanol-water partition coefficient
                  or distribution coefficient. XLogP is used as a measure of
                  hydrophilicity or hydrophobicity of a molecule
                </DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'ExactMass',
                    'The mass of the most likely isotopic composition for a single molecule, corresponding to the most intense ion/molecule peak in a mass spectrum'
                  )
                }>
                <DataTable.Cell>ExactMass</DataTable.Cell>
                <DataTable.Cell>
                  The mass of the most likely isotopic composition for a single
                  molecule, corresponding to the most intense ion/molecule peak
                  in a mass spectrum
                </DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'MonoisotopicMass',
                    'The mass of a molecule, calculated using the mass of the most abundant isotope of each element'
                  )
                }>
                <DataTable.Cell>MonoisotopicMass</DataTable.Cell>
                <DataTable.Cell>
                  The mass of a molecule, calculated using the mass of the most
                  abundant isotope of each element
                </DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'TPSA',
                    'Topological polar surface area, computed by the algorithm described in the paper by Ertl et al'
                  )
                }>
                <DataTable.Cell>TPSA</DataTable.Cell>
                <DataTable.Cell>
                  Topological polar surface area, computed by the algorithm
                  described in the paper by Ertl et al
                </DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'Complexity',
                    'The molecular complexity rating of a compound, computed using the Bertz/Hendrickson/Ihlenfeldt formula'
                  )
                }>
                <DataTable.Cell>Complexity</DataTable.Cell>
                <DataTable.Cell>
                  The molecular complexity rating of a compound, computed using
                  the Bertz/Hendrickson/Ihlenfeldt formula
                </DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'Charge',
                    'The total (or net) charge of a molecule'
                  )
                }>
                <DataTable.Cell>Charge</DataTable.Cell>
                <DataTable.Cell>
                  The total (or net) charge of a molecule
                </DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'HBondDonorCount',
                    'Number of hydrogen-bond donors in the structure'
                  )
                }>
                <DataTable.Cell>HBondDonorCount</DataTable.Cell>
                <DataTable.Cell>
                  Number of hydrogen-bond donors in the structure
                </DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'HBondAcceptorCount',
                    'Number of hydrogen-bond acceptors in the structure'
                  )
                }>
                <DataTable.Cell>HBondAcceptorCount</DataTable.Cell>
                <DataTable.Cell>
                  Number of hydrogen-bond acceptors in the structure
                </DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'RotatableBondCount',
                    'Number of rotatable bonds.'
                  )
                }>
                <DataTable.Cell>RotatableBondCount</DataTable.Cell>
                <DataTable.Cell>Number of rotatable bonds</DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert('HeavyAtomCount', 'Number of non-hydrogen atoms')
                }>
                <DataTable.Cell>HeavyAtomCount</DataTable.Cell>
                <DataTable.Cell>Number of non-hydrogen atoms</DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'IsotopeAtomCount',
                    'Number of atoms with enriched isotope(s)'
                  )
                }>
                <DataTable.Cell>IsotopeAtomCount</DataTable.Cell>
                <DataTable.Cell>
                  Number of atoms with enriched isotope(s)
                </DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'AtomStereoCount',
                    'Total number of atoms with tetrahedral (sp3) stereo [e.g., (R)- or (S)-configuration]'
                  )
                }>
                <DataTable.Cell>AtomStereoCount</DataTable.Cell>
                <DataTable.Cell>
                  Total number of atoms with tetrahedral (sp3) stereo [e.g.,
                  (R)- or (S)-configuration]
                </DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'DefinedAtomStereoCount',
                    'Number of atoms with defined tetrahedral (sp3) stereo'
                  )
                }>
                <DataTable.Cell>DefinedAtomStereoCount</DataTable.Cell>
                <DataTable.Cell>
                  Number of atoms with defined tetrahedral (sp3) stereo
                </DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'UndefinedAtomStereoCount',
                    'Number of atoms with undefined tetrahedral (sp3) stereo'
                  )
                }>
                <DataTable.Cell>UndefinedAtomStereoCount</DataTable.Cell>
                <DataTable.Cell>
                  Number of atoms with undefined tetrahedral (sp3) stereo
                </DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'BondStereoCount',
                    'Total number of bonds with planar (sp2) stereo [e.g., (E)- or (Z)-configuration]'
                  )
                }>
                <DataTable.Cell>BondStereoCount</DataTable.Cell>
                <DataTable.Cell>
                  Total number of bonds with planar (sp2) stereo [e.g., (E)- or
                  (Z)-configuration]
                </DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'DefinedBondStereoCount',
                    'Number of atoms with defined planar (sp2) stereo'
                  )
                }>
                <DataTable.Cell>DefinedBondStereoCount</DataTable.Cell>
                <DataTable.Cell>
                  Number of atoms with defined planar (sp2) stereo
                </DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'UndefinedBondStereoCount',
                    'Number of atoms with undefined planar (sp2) stereo'
                  )
                }>
                <DataTable.Cell>UndefinedBondStereoCount</DataTable.Cell>
                <DataTable.Cell>
                  Number of atoms with undefined planar (sp2) stereo
                </DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'CovalentUnitCount',
                    'Number of covalently bound units'
                  )
                }>
                <DataTable.Cell>CovalentUnitCount</DataTable.Cell>
                <DataTable.Cell>
                  Number of covalently bound units
                </DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'Volume3D',
                    'Analytic volume of the first diverse conformer (default conformer) for a compound'
                  )
                }>
                <DataTable.Cell>Volume3D</DataTable.Cell>
                <DataTable.Cell>
                  Analytic volume of the first diverse conformer (default
                  conformer) for a compound
                </DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'XStericQuadrupole3D',
                    'The x component of the quadrupole moment (Qx) of the first diverse conformer (default conformer) for a compound'
                  )
                }>
                <DataTable.Cell>XStericQuadrupole3D</DataTable.Cell>
                <DataTable.Cell>
                  The x component of the quadrupole moment (Qx) of the first
                  diverse conformer (default conformer) for a compound
                </DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'YStericQuadrupole3D',
                    'The y component of the quadrupole moment (Qy) of the first diverse conformer (default conformer) for a compound'
                  )
                }>
                <DataTable.Cell>YStericQuadrupole3D</DataTable.Cell>
                <DataTable.Cell>
                  The y component of the quadrupole moment (Qy) of the first
                  diverse conformer (default conformer) for a compound
                </DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'ZStericQuadrupole3D',
                    'The z component of the quadrupole moment (Qz) of the first diverse conformer (default conformer) for a compound'
                  )
                }>
                <DataTable.Cell>ZStericQuadrupole3D</DataTable.Cell>
                <DataTable.Cell>
                  The z component of the quadrupole moment (Qz) of the first
                  diverse conformer (default conformer) for a compound
                </DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'FeatureCount3D',
                    'Total number of 3D features (the sum of FeatureAcceptorCount3D, FeatureDonorCount3D, FeatureAnionCount3D, FeatureCationCount3D, FeatureRingCount3D and FeatureHydrophobeCount3D)'
                  )
                }>
                <DataTable.Cell>FeatureCount3D</DataTable.Cell>
                <DataTable.Cell>
                  Total number of 3D features (the sum of
                  FeatureAcceptorCount3D, FeatureDonorCount3D,
                  FeatureAnionCount3D, FeatureCationCount3D, FeatureRingCount3D
                  and FeatureHydrophobeCount3D)
                </DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'FeatureAcceptorCount3D',
                    'Number of hydrogen-bond acceptors of a conformer'
                  )
                }>
                <DataTable.Cell>FeatureAcceptorCount3D</DataTable.Cell>
                <DataTable.Cell>
                  Number of hydrogen-bond acceptors of a conformer
                </DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'FeatureDonorCount3D',
                    'Number of hydrogen-bond donors of a conformer'
                  )
                }>
                <DataTable.Cell>FeatureDonorCount3D</DataTable.Cell>
                <DataTable.Cell>
                  Number of hydrogen-bond donors of a conformer
                </DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'FeatureAnionCount3D',
                    'Number of anionic centers (at pH 7) of a conformer'
                  )
                }>
                <DataTable.Cell>FeatureAnionCount3D</DataTable.Cell>
                <DataTable.Cell>
                  Number of anionic centers (at pH 7) of a conformer
                </DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'FeatureCationCount3D',
                    'Number of cationic centers (at pH 7) of a conformer'
                  )
                }>
                <DataTable.Cell>FeatureCationCount3D</DataTable.Cell>
                <DataTable.Cell>
                  Number of cationic centers (at pH 7) of a conformer
                </DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'FeatureRingCount3D',
                    'Number of rings of a conformer'
                  )
                }>
                <DataTable.Cell>FeatureRingCount3D</DataTable.Cell>
                <DataTable.Cell>Number of rings of a conformer</DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'FeatureHydrophobeCount3D',
                    'Number of hydrophobes of a conformer'
                  )
                }>
                <DataTable.Cell>FeatureHydrophobeCount3D</DataTable.Cell>
                <DataTable.Cell>
                  Number of hydrophobes of a conformer
                </DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'ConformerModelRMSD3D',
                    'Conformer sampling RMSD in Å'
                  )
                }>
                <DataTable.Cell>ConformerModelRMSD3D</DataTable.Cell>
                <DataTable.Cell>Conformer sampling RMSD in Å</DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'EffectiveRotorCount3D',
                    'Total number of 3D features (the sum of FeatureAcceptorCount3D, FeatureDonorCount3D, FeatureAnionCount3D, FeatureCationCount3D, FeatureRingCount3D and FeatureHydrophobeCount3D)'
                  )
                }>
                <DataTable.Cell>EffectiveRotorCount3D</DataTable.Cell>
                <DataTable.Cell>
                  Total number of 3D features (the sum of
                  FeatureAcceptorCount3D, FeatureDonorCount3D,
                  FeatureAnionCount3D, FeatureCationCount3D, FeatureRingCount3D
                  and FeatureHydrophobeCount3D)
                </DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'ConformerCount3D',
                    'The number of conformers in the conformer model for a compound'
                  )
                }>
                <DataTable.Cell>ConformerCount3D</DataTable.Cell>
                <DataTable.Cell>
                  The number of conformers in the conformer model for a compound
                </DataTable.Cell>
              </DataTable.Row>

              <DataTable.Row
                style={{ height: 70 }}
                onPress={() =>
                  Alert.alert(
                    'Fingerprint2D',
                    'Base64-encoded PubChem Substructure Fingerprint of a molecule'
                  )
                }>
                <DataTable.Cell>Fingerprint2D</DataTable.Cell>
                <DataTable.Cell>
                  Base64-encoded PubChem Substructure Fingerprint of a molecule
                </DataTable.Cell>
              </DataTable.Row>
            </DataTable>
            <FlatList
              data={this.state.data_arr}
              renderItem={({ item }) => <this.DataBase molecule_data={item} navigation={navigation} />}
              keyExtractor={item => item['InChIKey']}
              style={{ marginBottom: 10 }}
            />
          </ScrollView>
        )}
      </View>
    );
  }
}

const styles = StyleSheet.create({
  button: {
    backgroundColor: '#0DA900',
    padding: 15,
    borderRadius: 5,
    marginBottom: 10,
  },
  buttonText: {
    fontSize: 15,
    color: '#fff',
  },
  input: {
    height: 50,
    borderColor: '#0DA900',
    borderWidth: 1,
    padding: 10,
    borderRadius: 5,
    marginBottom: 10,
    justifyContent: 'center',
    alignItems: 'center',
  },
  container: {
    flex: 1,
    justifyContent: 'center',
    backgroundColor: '#f3f3f3',
    padding: 8,
  },
});
