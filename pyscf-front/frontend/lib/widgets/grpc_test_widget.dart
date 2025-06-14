import 'package:flutter/material.dart';
import '../services/grpc_service.dart';
import '../generated/calculation.pb.dart';

class GrpcTestWidget extends StatefulWidget {
  const GrpcTestWidget({super.key});

  @override
  State<GrpcTestWidget> createState() => _GrpcTestWidgetState();
}

class _GrpcTestWidgetState extends State<GrpcTestWidget> {
  String _status = 'Not connected';
  String _systemInfo = '';
  String _testResults = '';
  bool _isLoading = false;

  @override
  void initState() {
    super.initState();
    _initializeGrpc();
  }

  Future<void> _initializeGrpc() async {
    setState(() {
      _isLoading = true;
      _status = 'Connecting...';
    });

    try {
      await grpcService.initialize();
      setState(() {
        _status = 'Connected to gRPC server';
      });
    } catch (e) {
      setState(() {
        _status = 'Connection failed: $e';
      });
    } finally {
      setState(() {
        _isLoading = false;
      });
    }
  }

  Future<void> _testHealthCheck() async {
    setState(() {
      _isLoading = true;
      _testResults = 'Testing health check...';
    });

    try {
      final response = await grpcService.healthCheck();
      setState(() {
        _testResults = 'Health Check: ${response.healthy ? "Healthy" : "Unhealthy"}\n'
            'Status: ${response.status}\n'
            'Details: ${response.details}';
      });
    } catch (e) {
      setState(() {
        _testResults = 'Health check failed: $e';
      });
    } finally {
      setState(() {
        _isLoading = false;
      });
    }
  }

  Future<void> _testSystemInfo() async {
    setState(() {
      _isLoading = true;
      _systemInfo = 'Getting system info...';
    });

    try {
      final response = await grpcService.getSystemInfo();
      if (response.success && response.hasSystemInfo()) {
        final info = response.systemInfo;
        setState(() {
          _systemInfo = 'System Info:\n'
              'Version: ${info.version}\n'
              'Python: ${info.pythonVersion}\n'
              'PySCF: ${info.pyscfVersion}\n'
              'GPU Available: ${info.gpuAvailable}\n'
              'CPU Cores: ${info.resources.cpuCores}\n'
              'Memory: ${info.resources.totalMemoryMb} MB\n'
              'Available Methods: ${info.availableMethods.join(", ")}\n'
              'Basis Sets: ${info.availableBasisSets.join(", ")}';
        });
      } else {
        setState(() {
          _systemInfo = 'Failed to get system info: ${response.message}';
        });
      }
    } catch (e) {
      setState(() {
        _systemInfo = 'System info failed: $e';
      });
    } finally {
      setState(() {
        _isLoading = false;
      });
    }
  }

  Future<void> _testCreateMolecule() async {
    setState(() {
      _isLoading = true;
      _testResults = 'Creating water molecule...';
    });

    try {
      // Create water molecule (H2O)
      final atoms = [
        Atom()
          ..symbol = 'O'
          ..x = 0.0
          ..y = 0.0
          ..z = 0.0,
        Atom()
          ..symbol = 'H'
          ..x = 0.757
          ..y = 0.586
          ..z = 0.0,
        Atom()
          ..symbol = 'H'
          ..x = -0.757
          ..y = 0.586
          ..z = 0.0,
      ];

      final response = await grpcService.createMolecule(
        name: 'Water',
        formula: 'H2O',
        atoms: atoms,
        charge: 0,
        multiplicity: 1,
      );

      if (response.success && response.hasMolecule()) {
        final molecule = response.molecule;
        setState(() {
          _testResults = 'Molecule Created Successfully!\n'
              'ID: ${molecule.id}\n'
              'Name: ${molecule.name}\n'
              'Formula: ${molecule.formula}\n'
              'Molecular Weight: ${molecule.molecularWeight.toStringAsFixed(2)}\n'
              'Atoms: ${molecule.atoms.length}\n'
              'Charge: ${molecule.charge}\n'
              'Multiplicity: ${molecule.multiplicity}';
        });
      } else {
        setState(() {
          _testResults = 'Failed to create molecule: ${response.message}';
        });
      }
    } catch (e) {
      setState(() {
        _testResults = 'Create molecule failed: $e';
      });
    } finally {
      setState(() {
        _isLoading = false;
      });
    }
  }

  @override
  Widget build(BuildContext context) {
    return Scaffold(
      appBar: AppBar(
        title: const Text('gRPC Communication Test'),
        backgroundColor: Theme.of(context).colorScheme.inversePrimary,
      ),
      body: Padding(
        padding: const EdgeInsets.all(16.0),
        child: Column(
          crossAxisAlignment: CrossAxisAlignment.start,
          children: [
            Card(
              child: Padding(
                padding: const EdgeInsets.all(16.0),
                child: Column(
                  crossAxisAlignment: CrossAxisAlignment.start,
                  children: [
                    Text(
                      'Connection Status',
                      style: Theme.of(context).textTheme.titleMedium,
                    ),
                    const SizedBox(height: 8),
                    Text(
                      _status,
                      style: TextStyle(
                        color: _status.contains('Connected')
                            ? Colors.green
                            : _status.contains('failed')
                                ? Colors.red
                                : Colors.orange,
                      ),
                    ),
                  ],
                ),
              ),
            ),
            const SizedBox(height: 16),
            Row(
              children: [
                Expanded(
                  child: ElevatedButton(
                    onPressed: _isLoading ? null : _testHealthCheck,
                    child: const Text('Health Check'),
                  ),
                ),
                const SizedBox(width: 8),
                Expanded(
                  child: ElevatedButton(
                    onPressed: _isLoading ? null : _testSystemInfo,
                    child: const Text('System Info'),
                  ),
                ),
                const SizedBox(width: 8),
                Expanded(
                  child: ElevatedButton(
                    onPressed: _isLoading ? null : _testCreateMolecule,
                    child: const Text('Create H2O'),
                  ),
                ),
              ],
            ),
            const SizedBox(height: 16),
            if (_isLoading)
              const Center(
                child: CircularProgressIndicator(),
              ),
            if (_systemInfo.isNotEmpty) ...[
              Card(
                child: Padding(
                  padding: const EdgeInsets.all(16.0),
                  child: Column(
                    crossAxisAlignment: CrossAxisAlignment.start,
                    children: [
                      Text(
                        'System Information',
                        style: Theme.of(context).textTheme.titleMedium,
                      ),
                      const SizedBox(height: 8),
                      Text(
                        _systemInfo,
                        style: const TextStyle(fontFamily: 'monospace'),
                      ),
                    ],
                  ),
                ),
              ),
              const SizedBox(height: 16),
            ],
            if (_testResults.isNotEmpty) ...[
              Card(
                child: Padding(
                  padding: const EdgeInsets.all(16.0),
                  child: Column(
                    crossAxisAlignment: CrossAxisAlignment.start,
                    children: [
                      Text(
                        'Test Results',
                        style: Theme.of(context).textTheme.titleMedium,
                      ),
                      const SizedBox(height: 8),
                      Text(
                        _testResults,
                        style: const TextStyle(fontFamily: 'monospace'),
                      ),
                    ],
                  ),
                ),
              ),
            ],
          ],
        ),
      ),
    );
  }

  @override
  void dispose() {
    grpcService.close();
    super.dispose();
  }
}