import * as React from 'react';
import { TouchableOpacity,Text, View, StyleSheet } from 'react-native';
import { NavigationContainer } from '@react-navigation/native';
import { createStackNavigator } from '@react-navigation/stack';
import Constants from 'expo-constants';

import MoleculeGeneration from './components/MoleculeGeneration.js';
import VirtualScreening from './components/VirtualScreening.js';

const Stack = createStackNavigator();

function Home(props) {
  return(
  <View style={styles.container}>
  <TouchableOpacity
    style={styles.button}
    onPress={() => props.navigation.navigate('Molecule Generation')}>
    <Text style={styles.buttonText}>신약 후보 물질 생성</Text>
  </TouchableOpacity>
  <TouchableOpacity
    style={styles.button}
    onPress={() => props.navigation.navigate('Virtual Screening', {
        data: ''
    })}>
    <Text style={styles.buttonText}>Virtual Screening</Text>
  </TouchableOpacity>
  </View>
  );
}

export default function App() {
  
  return (
    <NavigationContainer>
      <Stack.Navigator>
        <Stack.Screen name="Home" component={Home} />
        <Stack.Screen name="Molecule Generation" component={MoleculeGeneration} initialRouteName='Home'
        screenOptions={{
          gestureEnabled: true
        }} />
        <Stack.Screen name="Virtual Screening" component={VirtualScreening} initialRouteName='Home'
        screenOptions={{
          gestureEnabled: true
        }}/>
      </Stack.Navigator>
    </NavigationContainer>
  );
}

const styles = StyleSheet.create({
  container: {
    flex: 1,
    justifyContent: 'center',
    backgroundColor: '#f3f3f3',
    padding: 8,
  },
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
});
