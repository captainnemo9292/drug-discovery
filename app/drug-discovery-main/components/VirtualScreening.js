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

export default class VirtualScreening extends React.Component {

    state = {
      ligand_smiles: '',
      protein_smiles: '',
      show_input: true,
      show_data: false,
      show_loading: false,
      binding_affinity: '',
      toxicity: '',
    };

    async BindingAffinity(ligand_smiles, protein_smiles) {

      try {
        let response = await fetch(
          'https://virtual-screening-dv4wcr3seq-an.a.run.app/binding_affinity',
          {
            method: 'POST',
            headers: {
              Accept: 'application/json',
              'Content-Type': 'application/json',
            },
            body: JSON.stringify({
              ligand_smiles: ligand_smiles,
              protein_smiles: protein_smiles
            }),
          }
        );
        let json = await response.json();
        this.setState({
            binding_affinity: json['binding_affinity']
        });
      } catch (error) {
        console.error(error);
      }
    };

    async Toxicity(ligand_smiles) {
      try {
        let response = await fetch(
          'https://virtual-screening-dv4wcr3seq-an.a.run.app/toxicity',
          {
            method: 'POST',
            headers: {
              Accept: 'application/json',
              'Content-Type': 'application/json',
            },
            body: JSON.stringify({
              smiles: ligand_smiles,
            }),
          }
        );
        let json = await response.json();
        this.setState({
            toxicity: json
        });

      } catch (error) {
        console.error(error);
        this.setState({ show_loading: false });
        this.setState({ show_input: true });
        this.FailedAlert();
      }
    }

    FailedAlert() {
      Alert.alert(
        'Error',
        'API 통신 중 문제가 발생했습니다. 다시 시도해 주십시오...'
      );
    }

    DataProcessing(ligand_smiles, protein_smiles) {
      console.log(ligand_smiles, protein_smiles);
      this.Toxicity(ligand_smiles);
      this.BindingAffinity(ligand_smiles, protein_smiles);
      this.setState({ show_loading: false });
      this.setState({ show_input: false });
      this.setState({ show_data: true });
    }

    Back() {
      this.setState({ show_loading: false });
      this.setState({ show_data: false });
      this.setState({ show_input: true });
    }

  render() {
    const { navigation, route } = this.props;
    const params = route.params;

    console.log(params);
    return(
      <View style={styles.container}>
        {this.state.show_input && (
          <View>
            <TextInput
              style={[styles.input]}
              defaultValue={params.data}
              placeholder={"Ligand 물질을 SMILES 형식으로 입력하세요."}
              onChangeText={ligand_smiles => this.setState({ ligand_smiles })}
            />
            <TextInput
              style={[styles.input]}
              placeholder="Protein 물질을 SMILES 형식으로 입력하세요."
              onChangeText={protein_smiles => this.setState({ protein_smiles })}
            />
            <TouchableOpacity
              onPress={() => this.DataProcessing(this.state.ligand_smiles, this.state.protein_smiles)}
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
        <ScrollView style={{marginTop:20}}>
        <TouchableOpacity onPress={() => this.Back()} style={styles.button}>
          <Text style={styles.buttonText}>Virtual Screening 다시 진행</Text>
        </TouchableOpacity>
          <DataTable>
            <DataTable.Header style={{ height: 80 }}>
              <DataTable.Title>Binding Affinity 데이터</DataTable.Title>
            </DataTable.Header>
            <DataTable.Row
              style={{ height: 70 }}
              onPress={() =>
                Alert.alert(
                  'Binding Affinity (-log Ki)',
                  String(Number(this.state.binding_affinity).toFixed(2))
                )
              }>
              <DataTable.Cell>Binding Affinity (-log Ki)</DataTable.Cell>
              <DataTable.Cell>
                {String(Number(this.state.binding_affinity).toFixed(2))}
              </DataTable.Cell>
            </DataTable.Row>
          </DataTable>

          <DataTable>
            <DataTable.Header style={{ height: 80 }}>
              <DataTable.Title>Toxicity 데이터</DataTable.Title>
            </DataTable.Header>
            <DataTable.Row
              style={{ height: 70 }}
              onPress={() =>
                Alert.alert(
                  'estrogen receptor alpha, LBD (ER, LBD)',
                  String(Number(this.state.toxicity['estrogen receptor alpha, LBD (ER, LBD)']).toFixed(2))
                )
              }>
              <DataTable.Cell>estrogen receptor alpha, LBD (ER, LBD)</DataTable.Cell>
              <DataTable.Cell>
                {String(Number(this.state.toxicity['estrogen receptor alpha, LBD (ER, LBD)']).toFixed(2))}
              </DataTable.Cell>
            </DataTable.Row>
            <DataTable.Row
              style={{ height: 70 }}
              onPress={() =>
                Alert.alert(
                  'estrogen receptor alpha, full (ER, full)',
                  String(Number(this.state.toxicity['estrogen receptor alpha, full (ER, full)']).toFixed(2))
                )
              }>
              <DataTable.Cell>estrogen receptor alpha, full (ER, full)</DataTable.Cell>
              <DataTable.Cell>
                {String(Number(this.state.toxicity['estrogen receptor alpha, full (ER, full)']).toFixed(2))}
              </DataTable.Cell>
            </DataTable.Row>
            <DataTable.Row
              style={{ height: 70 }}
              onPress={() =>
                Alert.alert(
                  'estrogen receptor alpha, LBD (ER, LBD)',
                  String(Number(this.state.toxicity['aromatase']).toFixed(2))
                )
              }>
              <DataTable.Cell>aromatase</DataTable.Cell>
              <DataTable.Cell>
                {String(Number(this.state.toxicity['aromatase']).toFixed(2))}
              </DataTable.Cell>
            </DataTable.Row>
            <DataTable.Row
              style={{ height: 70 }}
              onPress={() =>
                Alert.alert(
                  'aryl hydrocarbon receptor (AhR)',
                  String(Number(this.state.toxicity['aryl hydrocarbon receptor (AhR)']).toFixed(2))
                )
              }>
              <DataTable.Cell>aryl hydrocarbon receptor (AhR)</DataTable.Cell>
              <DataTable.Cell>
                {String(Number(this.state.toxicity['aryl hydrocarbon receptor (AhR)']).toFixed(2))}
              </DataTable.Cell>
            </DataTable.Row>
            <DataTable.Row
              style={{ height: 70 }}
              onPress={() =>
                Alert.alert(
                  'androgen receptor, full (AR, full)',
                  String(Number(this.state.toxicity['androgen receptor, full (AR, full)']).toFixed(2))
                )
              }>
              <DataTable.Cell>androgen receptor, full (AR, full)</DataTable.Cell>
              <DataTable.Cell>
                {String(Number(this.state.toxicity['androgen receptor, full (AR, full)']).toFixed(2))}
              </DataTable.Cell>
            </DataTable.Row>
            <DataTable.Row
              style={{ height: 70 }}
              onPress={() =>
                Alert.alert(
                  'androgen receptor, LBD (AR, LBD)',
                  String(Number(this.state.toxicity['androgen receptor, LBD (AR, LBD)']).toFixed(2))
                )
              }>
              <DataTable.Cell>androgen receptor, LBD (AR, LBD)</DataTable.Cell>
              <DataTable.Cell>
                {String(Number(this.state.toxicity['androgen receptor, LBD (AR, LBD)']).toFixed(2))}
              </DataTable.Cell>
            </DataTable.Row>
            <DataTable.Row
              style={{ height: 70 }}
              onPress={() =>
                Alert.alert(
                  'peroxisome proliferator-activated receptor gamma (PPAR-gamma)',
                  String(Number(this.state.toxicity['peroxisome proliferator-activated receptor gamma (PPAR-gamma)']).toFixed(2))
                )
              }>
              <DataTable.Cell>peroxisome proliferator-activated receptor gamma (PPAR-gamma)</DataTable.Cell>
              <DataTable.Cell>
                {String(Number(this.state.toxicity['peroxisome proliferator-activated receptor gamma (PPAR-gamma)']).toFixed(2))}
              </DataTable.Cell>
            </DataTable.Row>
            <DataTable.Row
              style={{ height: 70 }}
              onPress={() =>
                Alert.alert(
                  'nuclear factor (erythroid-derived 2)-like 2/antioxidant responsive element (Nrf2/ARE)',
                  String(Number(this.state.toxicity['nuclear factor (erythroid-derived 2)-like 2/antioxidant responsive element (Nrf2/ARE)']).toFixed(2))
                )
              }>
              <DataTable.Cell>nuclear factor (erythroid-derived 2)-like 2/antioxidant responsive element (Nrf2/ARE)</DataTable.Cell>
              <DataTable.Cell>
                {String(Number(this.state.toxicity['nuclear factor (erythroid-derived 2)-like 2/antioxidant responsive element (Nrf2/ARE)']).toFixed(2))}
              </DataTable.Cell>
            </DataTable.Row>
            <DataTable.Row
              style={{ height: 70 }}
              onPress={() =>
                Alert.alert(
                  'heat shock factor response element (HSE)',
                  String(Number(this.state.toxicity['heat shock factor response element (HSE)']).toFixed(2))
                )
              }>
              <DataTable.Cell>heat shock factor response element (HSE)</DataTable.Cell>
              <DataTable.Cell>
                {String(Number(this.state.toxicity['heat shock factor response element (HSE)']).toFixed(2))}
              </DataTable.Cell>
            </DataTable.Row>
            <DataTable.Row
              style={{ height: 70 }}
              onPress={() =>
                Alert.alert(
                  'ATAD5',
                  String(Number(this.state.toxicity['ATAD5']).toFixed(2))
                )
              }>
              <DataTable.Cell>ATAD5</DataTable.Cell>
              <DataTable.Cell>
                {String(Number(this.state.toxicity['ATAD5']).toFixed(2))}
              </DataTable.Cell>
            </DataTable.Row>
            <DataTable.Row
              style={{ height: 70 }}
              onPress={() =>
                Alert.alert(
                  'mitochondrial membrane potential (MMP)',
                  String(Number(this.state.toxicity['mitochondrial membrane potential (MMP)']).toFixed(2))
                )
              }>
              <DataTable.Cell>mitochondrial membrane potential (MMP)</DataTable.Cell>
              <DataTable.Cell>
                {String(Number(this.state.toxicity['mitochondrial membrane potential (MMP)']).toFixed(2))}
              </DataTable.Cell>
            </DataTable.Row>
            <DataTable.Row
              style={{ height: 70 }}
              onPress={() =>
                Alert.alert(
                  'p53',
                  String(Number(this.state.toxicity['p53']).toFixed(2))
                )
              }>
              <DataTable.Cell>p53</DataTable.Cell>
              <DataTable.Cell>
                {String(Number(this.state.toxicity['p53']).toFixed(2))}
              </DataTable.Cell>
            </DataTable.Row>
          </DataTable>
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
