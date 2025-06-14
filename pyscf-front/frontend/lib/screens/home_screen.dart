import 'package:flutter/material.dart';
import '../widgets/molecule_input_widget.dart';
import '../widgets/calculation_settings_widget.dart';
import '../widgets/results_display_widget.dart';
import '../widgets/molecule_visualizer_widget.dart';
import '../models/molecule.dart';
import '../models/calculation_result.dart';
import '../services/grpc_service.dart';
import '../generated/calculation.pb.dart';

class HomeScreen extends StatefulWidget {
  const HomeScreen({super.key});

  @override
  State<HomeScreen> createState() => _HomeScreenState();
}

class _HomeScreenState extends State<HomeScreen> {
  MoleculeData? _currentMolecule;
  CalculationResult? _calculationResult;
  bool _isCalculating = false;
  String _status = 'Ready';

  // Calculation settings
  String _selectedMethod = 'HF';
  String _selectedBasisSet = 'sto-3g';
  int _charge = 0;
  int _multiplicity = 1;

  @override
  void initState() {
    super.initState();
    _initializeGrpc();
  }

  Future<void> _initializeGrpc() async {
    try {
      await grpcService.initialize();
      setState(() {
        _status = 'Connected to PySCF backend';
      });
    } catch (e) {
      setState(() {
        _status = 'Failed to connect to backend: $e';
      });
    }
  }

  void _onMoleculeChanged(MoleculeData molecule) {
    setState(() {
      _currentMolecule = molecule;
      _calculationResult = null; // Clear previous results
      _status = 'Molecule loaded: ${molecule.name}';
    });
  }

  void _onCalculationSettingsChanged({
    String? method,
    String? basisSet,
    int? charge,
    int? multiplicity,
  }) {
    setState(() {
      if (method != null) _selectedMethod = method;
      if (basisSet != null) _selectedBasisSet = basisSet;
      if (charge != null) _charge = charge;
      if (multiplicity != null) _multiplicity = multiplicity;
    });
  }

  Future<void> _startCalculation() async {
    if (_currentMolecule == null) {
      ScaffoldMessenger.of(context).showSnackBar(
        const SnackBar(content: Text('Please load a molecule first')),
      );
      return;
    }

    setState(() {
      _isCalculating = true;
      _status = 'Starting calculation...';
    });

    try {
      // First create the molecule in the backend
      final atoms = _currentMolecule!.atoms.map((atom) {
        return Atom()
          ..symbol = atom.symbol
          ..x = atom.x
          ..y = atom.y
          ..z = atom.z;
      }).toList();

      final moleculeResponse = await grpcService.createMolecule(
        name: _currentMolecule!.name,
        formula: _currentMolecule!.formula,
        atoms: atoms,
        charge: _charge,
        multiplicity: _multiplicity,
      );

      if (!moleculeResponse.success || !moleculeResponse.hasMolecule()) {
        throw Exception('Failed to create molecule: ${moleculeResponse.message}');
      }

      setState(() {
        _status = 'Creating calculation instance...';
      });

      // Create calculation instance
      final instanceResponse = await grpcService.createInstance(
        name: '${_currentMolecule!.name} ${_selectedMethod}/${_selectedBasisSet}',
        description: 'Calculation with method=$_selectedMethod, basis=$_selectedBasisSet',
        moleculeId: moleculeResponse.molecule?.id ?? '',
      );

      if (!instanceResponse.success || !instanceResponse.hasInstance()) {
        throw Exception('Failed to create instance: ${instanceResponse.message}');
      }

      setState(() {
        _status = 'Running $_selectedMethod calculation...';
      });

      // Start the calculation
      final calculationStream = grpcService.startCalculation(
        instanceId: instanceResponse.instance?.id ?? '',
        method: _selectedMethod,
        basisSet: _selectedBasisSet,
        priority: 1,
      );

      String lastMessage = '';
      await for (final progress in calculationStream) {
        setState(() {
          _status = '${progress.currentStep}: ${progress.progressPercentage.toStringAsFixed(1)}%';
          lastMessage = progress.message;
        });

        // Check if calculation is completed
        if (progress.status == CalculationStatus.CALCULATION_STATUS_COMPLETED) {
          setState(() {
            _calculationResult = CalculationResult(
              success: true,
              message: 'Calculation completed successfully',
              energy: 0.0, // Will be filled from actual results
              calculationTime: 0, // Will be filled from actual results
            );
            _status = 'Calculation completed successfully';
          });
          break;
        } else if (progress.status == CalculationStatus.CALCULATION_STATUS_FAILED) {
          throw Exception('Calculation failed: $lastMessage');
        }
      }
    } catch (e) {
      setState(() {
        _calculationResult = CalculationResult(
          success: false,
          message: 'Calculation failed: $e',
          energy: 0.0,
          calculationTime: 0,
        );
        _status = 'Calculation failed';
      });
      
      ScaffoldMessenger.of(context).showSnackBar(
        SnackBar(content: Text('Calculation failed: $e')),
      );
    } finally {
      setState(() {
        _isCalculating = false;
      });
    }
  }

  @override
  Widget build(BuildContext context) {
    return Scaffold(
      appBar: AppBar(
        title: const Text('PySCF Frontend'),
        backgroundColor: Theme.of(context).colorScheme.inversePrimary,
        actions: [
          IconButton(
            icon: const Icon(Icons.info_outline),
            onPressed: () {
              showDialog(
                context: context,
                builder: (context) => AlertDialog(
                  title: const Text('About'),
                  content: const Text(
                    'PySCF Frontend - Quantum Chemistry Calculations\n\n'
                    'This application provides a user-friendly interface for '
                    'running quantum chemistry calculations using PySCF.',
                  ),
                  actions: [
                    TextButton(
                      onPressed: () => Navigator.of(context).pop(),
                      child: const Text('OK'),
                    ),
                  ],
                ),
              );
            },
          ),
        ],
      ),
      body: Column(
        children: [
          // Status bar
          Container(
            width: double.infinity,
            padding: const EdgeInsets.all(8.0),
            color: _status.contains('Failed') || _status.contains('failed')
                ? Colors.red.withOpacity(0.1)
                : _status.contains('completed')
                    ? Colors.green.withOpacity(0.1)
                    : Colors.blue.withOpacity(0.1),
            child: Text(
              _status,
              style: TextStyle(
                color: _status.contains('Failed') || _status.contains('failed')
                    ? Colors.red.shade700
                    : _status.contains('completed')
                        ? Colors.green.shade700
                        : Colors.blue.shade700,
                fontWeight: FontWeight.w500,
              ),
            ),
          ),
          
          // Main content
          Expanded(
            child: Row(
              children: [
                // Left panel - Input and settings
                Expanded(
                  flex: 1,
                  child: Column(
                    children: [
                      // Molecule input
                      Expanded(
                        flex: 2,
                        child: Card(
                          margin: const EdgeInsets.all(8.0),
                          child: Padding(
                            padding: const EdgeInsets.all(16.0),
                            child: MoleculeInputWidget(
                              onMoleculeChanged: _onMoleculeChanged,
                            ),
                          ),
                        ),
                      ),
                      
                      // Calculation settings
                      Expanded(
                        flex: 1,
                        child: Card(
                          margin: const EdgeInsets.all(8.0),
                          child: Padding(
                            padding: const EdgeInsets.all(16.0),
                            child: CalculationSettingsWidget(
                              selectedMethod: _selectedMethod,
                              selectedBasisSet: _selectedBasisSet,
                              charge: _charge,
                              multiplicity: _multiplicity,
                              onSettingsChanged: _onCalculationSettingsChanged,
                            ),
                          ),
                        ),
                      ),
                      
                      // Run calculation button
                      Padding(
                        padding: const EdgeInsets.all(16.0),
                        child: SizedBox(
                          width: double.infinity,
                          height: 48,
                          child: ElevatedButton.icon(
                            onPressed: _isCalculating || _currentMolecule == null
                                ? null
                                : _startCalculation,
                            icon: _isCalculating
                                ? const SizedBox(
                                    width: 20,
                                    height: 20,
                                    child: CircularProgressIndicator(
                                      strokeWidth: 2,
                                    ),
                                  )
                                : const Icon(Icons.play_arrow),
                            label: Text(_isCalculating 
                                ? 'Calculating...' 
                                : 'Run Calculation'),
                            style: ElevatedButton.styleFrom(
                              backgroundColor: Theme.of(context).primaryColor,
                              foregroundColor: Colors.white,
                            ),
                          ),
                        ),
                      ),
                    ],
                  ),
                ),
                
                // Right panel - Visualization and results
                Expanded(
                  flex: 2,
                  child: Column(
                    children: [
                      // 3D Molecule viewer
                      Expanded(
                        flex: 1,
                        child: Card(
                          margin: const EdgeInsets.all(8.0),
                          child: MoleculeVisualizerWidget(
                            molecule: _currentMolecule,
                          ),
                        ),
                      ),
                      
                      // Results display
                      if (_calculationResult != null)
                        Expanded(
                          flex: 1,
                          child: Card(
                            margin: const EdgeInsets.all(8.0),
                            child: Padding(
                              padding: const EdgeInsets.all(16.0),
                              child: ResultsDisplayWidget(
                                result: _calculationResult!,
                                molecule: _currentMolecule,
                              ),
                            ),
                          ),
                        ),
                    ],
                  ),
                ),
              ],
            ),
          ),
        ],
      ),
    );
  }

  @override
  void dispose() {
    grpcService.close();
    super.dispose();
  }
}