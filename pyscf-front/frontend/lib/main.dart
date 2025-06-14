import 'package:flutter/material.dart';
import 'screens/home_screen.dart';

void main() {
  runApp(const PySCFApp());
}

class PySCFApp extends StatelessWidget {
  const PySCFApp({super.key});

  @override
  Widget build(BuildContext context) {
    return MaterialApp(
      title: 'PySCF Frontend',
      theme: ThemeData(
        colorScheme: ColorScheme.fromSeed(seedColor: Colors.blue),
        useMaterial3: true,
      ),
      home: const HomeScreen(),
      debugShowCheckedModeBanner: false,
    );
  }
}

