import 'package:flutter/material.dart';

class CalculationSettingsWidget extends StatefulWidget {
  final String selectedMethod;
  final String selectedBasisSet;
  final int charge;
  final int multiplicity;
  final Function({
    String? method,
    String? basisSet,
    int? charge,
    int? multiplicity,
  }) onSettingsChanged;

  const CalculationSettingsWidget({
    super.key,
    required this.selectedMethod,
    required this.selectedBasisSet,
    required this.charge,
    required this.multiplicity,
    required this.onSettingsChanged,
  });

  @override
  State<CalculationSettingsWidget> createState() => _CalculationSettingsWidgetState();
}

class _CalculationSettingsWidgetState extends State<CalculationSettingsWidget> {
  late TextEditingController _chargeController;
  late TextEditingController _multiplicityController;
  bool _showAdvancedSettings = false;

  // Available calculation methods organized by category
  static const Map<String, List<String>> _methodCategories = {
    'Hartree-Fock': ['HF', 'RHF', 'UHF', 'ROHF'],
    'DFT (GGA)': ['PBE', 'BLYP', 'PW91'],
    'DFT (Hybrid)': ['B3LYP', 'PBE0', 'HSE06'],
    'DFT (Meta-GGA)': ['TPSS', 'M06-L'],
    'DFT (Range-separated)': ['CAM-B3LYP', 'ωB97X-D'],
    'Post-HF': ['MP2', 'CCSD', 'CCSD(T)'],
  };

  // Available basis sets organized by category
  static const Map<String, List<String>> _basisSetCategories = {
    'Minimal': ['sto-3g'],
    'Split-valence': ['3-21g', '6-31g', '6-311g'],
    'Polarized': ['6-31g*', '6-31g**', '6-311g*', '6-311g**'],
    'Diffuse': ['6-31+g', '6-31++g', '6-311+g', '6-311++g'],
    'Correlation-consistent': ['cc-pvdz', 'cc-pvtz', 'cc-pvqz'],
    'Dunning (aug)': ['aug-cc-pvdz', 'aug-cc-pvtz'],
  };

  // Method descriptions for educational purposes
  static const Map<String, String> _methodDescriptions = {
    'HF': 'Hartree-Fock: Mean-field approximation, fast but no correlation',
    'B3LYP': 'B3LYP: Popular hybrid DFT functional, good balance of accuracy and speed',
    'PBE': 'PBE: GGA functional, good for geometry optimization',
    'MP2': 'MP2: Second-order perturbation theory, includes electron correlation',
    'CCSD': 'CCSD: Coupled cluster singles and doubles, high accuracy',
    'CCSD(T)': 'CCSD(T): Gold standard for small molecules, very expensive',
  };

  static const Map<String, String> _basisSetDescriptions = {
    'sto-3g': 'Minimal basis set, very fast but inaccurate',
    '6-31g': 'Split-valence basis, good balance for initial calculations',
    '6-31g*': 'Adds polarization functions, better for bonding',
    'cc-pvdz': 'Correlation-consistent, systematically improvable',
    'aug-cc-pvdz': 'Adds diffuse functions, important for anions and weak interactions',
  };

  @override
  void initState() {
    super.initState();
    _chargeController = TextEditingController(text: widget.charge.toString());
    _multiplicityController = TextEditingController(text: widget.multiplicity.toString());
  }

  @override
  void dispose() {
    _chargeController.dispose();
    _multiplicityController.dispose();
    super.dispose();
  }

  void _updateCharge(String value) {
    final charge = int.tryParse(value) ?? 0;
    widget.onSettingsChanged(charge: charge);
  }

  void _updateMultiplicity(String value) {
    final multiplicity = int.tryParse(value) ?? 1;
    if (multiplicity < 1) return;
    widget.onSettingsChanged(multiplicity: multiplicity);
  }

  Widget _buildMethodSelector() {
    return Column(
      crossAxisAlignment: CrossAxisAlignment.start,
      children: [
        Text(
          'Calculation Method',
          style: Theme.of(context).textTheme.titleSmall,
        ),
        const SizedBox(height: 8),
        
        Container(
          width: double.infinity,
          padding: const EdgeInsets.symmetric(horizontal: 12),
          decoration: BoxDecoration(
            border: Border.all(color: Theme.of(context).dividerColor),
            borderRadius: BorderRadius.circular(4),
          ),
          child: DropdownButtonHideUnderline(
            child: DropdownButton<String>(
              value: widget.selectedMethod,
              isExpanded: true,
              onChanged: (String? newValue) {
                if (newValue != null) {
                  widget.onSettingsChanged(method: newValue);
                }
              },
              items: _methodCategories.entries.expand((category) {
                return [
                  // Category header
                  DropdownMenuItem<String>(
                    enabled: false,
                    value: null,
                    child: Text(
                      category.key,
                      style: Theme.of(context).textTheme.labelSmall?.copyWith(
                        fontWeight: FontWeight.bold,
                        color: Theme.of(context).primaryColor,
                      ),
                    ),
                  ),
                  // Methods in category
                  ...category.value.map((method) => DropdownMenuItem<String>(
                    value: method,
                    child: Padding(
                      padding: const EdgeInsets.only(left: 16),
                      child: Text(method),
                    ),
                  )),
                ];
              }).toList(),
            ),
          ),
        ),
        
        // Method description
        if (_methodDescriptions.containsKey(widget.selectedMethod))
          Container(
            width: double.infinity,
            padding: const EdgeInsets.all(8),
            margin: const EdgeInsets.only(top: 8),
            decoration: BoxDecoration(
              color: Theme.of(context).primaryColor.withOpacity(0.1),
              borderRadius: BorderRadius.circular(4),
            ),
            child: Text(
              _methodDescriptions[widget.selectedMethod]!,
              style: Theme.of(context).textTheme.bodySmall,
            ),
          ),
      ],
    );
  }

  Widget _buildBasisSetSelector() {
    return Column(
      crossAxisAlignment: CrossAxisAlignment.start,
      children: [
        Text(
          'Basis Set',
          style: Theme.of(context).textTheme.titleSmall,
        ),
        const SizedBox(height: 8),
        
        Container(
          width: double.infinity,
          padding: const EdgeInsets.symmetric(horizontal: 12),
          decoration: BoxDecoration(
            border: Border.all(color: Theme.of(context).dividerColor),
            borderRadius: BorderRadius.circular(4),
          ),
          child: DropdownButtonHideUnderline(
            child: DropdownButton<String>(
              value: widget.selectedBasisSet,
              isExpanded: true,
              onChanged: (String? newValue) {
                if (newValue != null) {
                  widget.onSettingsChanged(basisSet: newValue);
                }
              },
              items: _basisSetCategories.entries.expand((category) {
                return [
                  // Category header
                  DropdownMenuItem<String>(
                    enabled: false,
                    value: null,
                    child: Text(
                      category.key,
                      style: Theme.of(context).textTheme.labelSmall?.copyWith(
                        fontWeight: FontWeight.bold,
                        color: Theme.of(context).primaryColor,
                      ),
                    ),
                  ),
                  // Basis sets in category
                  ...category.value.map((basisSet) => DropdownMenuItem<String>(
                    value: basisSet,
                    child: Padding(
                      padding: const EdgeInsets.only(left: 16),
                      child: Text(basisSet),
                    ),
                  )),
                ];
              }).toList(),
            ),
          ),
        ),
        
        // Basis set description
        if (_basisSetDescriptions.containsKey(widget.selectedBasisSet))
          Container(
            width: double.infinity,
            padding: const EdgeInsets.all(8),
            margin: const EdgeInsets.only(top: 8),
            decoration: BoxDecoration(
              color: Theme.of(context).primaryColor.withOpacity(0.1),
              borderRadius: BorderRadius.circular(4),
            ),
            child: Text(
              _basisSetDescriptions[widget.selectedBasisSet]!,
              style: Theme.of(context).textTheme.bodySmall,
            ),
          ),
      ],
    );
  }

  Widget _buildBasicSettings() {
    return Column(
      children: [
        Row(
          children: [
            Expanded(
              child: Column(
                crossAxisAlignment: CrossAxisAlignment.start,
                children: [
                  Text(
                    'Charge',
                    style: Theme.of(context).textTheme.titleSmall,
                  ),
                  const SizedBox(height: 8),
                  TextField(
                    controller: _chargeController,
                    decoration: const InputDecoration(
                      border: OutlineInputBorder(),
                      hintText: '0',
                    ),
                    keyboardType: TextInputType.number,
                    onChanged: _updateCharge,
                  ),
                ],
              ),
            ),
            const SizedBox(width: 16),
            Expanded(
              child: Column(
                crossAxisAlignment: CrossAxisAlignment.start,
                children: [
                  Text(
                    'Multiplicity',
                    style: Theme.of(context).textTheme.titleSmall,
                  ),
                  const SizedBox(height: 8),
                  TextField(
                    controller: _multiplicityController,
                    decoration: const InputDecoration(
                      border: OutlineInputBorder(),
                      hintText: '1',
                    ),
                    keyboardType: TextInputType.number,
                    onChanged: _updateMultiplicity,
                  ),
                ],
              ),
            ),
          ],
        ),
        
        const SizedBox(height: 16),
        
        // Spin state helper
        Container(
          width: double.infinity,
          padding: const EdgeInsets.all(12),
          decoration: BoxDecoration(
            color: Colors.blue.withOpacity(0.1),
            borderRadius: BorderRadius.circular(4),
          ),
          child: Column(
            crossAxisAlignment: CrossAxisAlignment.start,
            children: [
              Text(
                'Spin State Guide:',
                style: Theme.of(context).textTheme.labelMedium?.copyWith(
                  fontWeight: FontWeight.bold,
                ),
              ),
              const SizedBox(height: 4),
              Text(
                '• Multiplicity = 2S + 1 (S = total spin)',
                style: Theme.of(context).textTheme.bodySmall,
              ),
              Text(
                '• Singlet: M=1 (closed shell), Doublet: M=2 (1 unpaired e⁻)',
                style: Theme.of(context).textTheme.bodySmall,
              ),
              Text(
                '• Triplet: M=3 (2 unpaired e⁻), Quartet: M=4 (3 unpaired e⁻)',
                style: Theme.of(context).textTheme.bodySmall,
              ),
            ],
          ),
        ),
      ],
    );
  }

  Widget _buildAdvancedSettings() {
    return Column(
      crossAxisAlignment: CrossAxisAlignment.start,
      children: [
        Text(
          'Advanced Settings',
          style: Theme.of(context).textTheme.titleSmall,
        ),
        const SizedBox(height: 16),
        
        // Convergence criteria
        Row(
          children: [
            Expanded(
              child: Column(
                crossAxisAlignment: CrossAxisAlignment.start,
                children: [
                  Text(
                    'Energy Convergence',
                    style: Theme.of(context).textTheme.labelMedium,
                  ),
                  const SizedBox(height: 8),
                  const TextField(
                    decoration: InputDecoration(
                      border: OutlineInputBorder(),
                      hintText: '1e-9',
                      enabled: false, // Disabled for MVP
                    ),
                  ),
                ],
              ),
            ),
            const SizedBox(width: 16),
            Expanded(
              child: Column(
                crossAxisAlignment: CrossAxisAlignment.start,
                children: [
                  Text(
                    'Max Iterations',
                    style: Theme.of(context).textTheme.labelMedium,
                  ),
                  const SizedBox(height: 8),
                  const TextField(
                    decoration: InputDecoration(
                      border: OutlineInputBorder(),
                      hintText: '100',
                      enabled: false, // Disabled for MVP
                    ),
                  ),
                ],
              ),
            ),
          ],
        ),
        
        const SizedBox(height: 16),
        
        // Memory and parallelization
        Row(
          children: [
            Expanded(
              child: Column(
                crossAxisAlignment: CrossAxisAlignment.start,
                children: [
                  Text(
                    'Memory (GB)',
                    style: Theme.of(context).textTheme.labelMedium,
                  ),
                  const SizedBox(height: 8),
                  const TextField(
                    decoration: InputDecoration(
                      border: OutlineInputBorder(),
                      hintText: '4',
                      enabled: false, // Disabled for MVP
                    ),
                  ),
                ],
              ),
            ),
            const SizedBox(width: 16),
            Expanded(
              child: Column(
                crossAxisAlignment: CrossAxisAlignment.start,
                children: [
                  Text(
                    'CPU Cores',
                    style: Theme.of(context).textTheme.labelMedium,
                  ),
                  const SizedBox(height: 8),
                  const TextField(
                    decoration: InputDecoration(
                      border: OutlineInputBorder(),
                      hintText: 'Auto',
                      enabled: false, // Disabled for MVP
                    ),
                  ),
                ],
              ),
            ),
          ],
        ),
        
        const SizedBox(height: 16),
        
        // Note about advanced settings
        Container(
          width: double.infinity,
          padding: const EdgeInsets.all(12),
          decoration: BoxDecoration(
            color: Colors.orange.withOpacity(0.1),
            borderRadius: BorderRadius.circular(4),
          ),
          child: Text(
            'Note: Advanced settings will be available in Phase 2. Default values are optimized for most calculations.',
            style: Theme.of(context).textTheme.bodySmall,
          ),
        ),
      ],
    );
  }

  Widget _buildRecommendations() {
    return Container(
      width: double.infinity,
      padding: const EdgeInsets.all(12),
      decoration: BoxDecoration(
        color: Colors.green.withOpacity(0.1),
        borderRadius: BorderRadius.circular(4),
      ),
      child: Column(
        crossAxisAlignment: CrossAxisAlignment.start,
        children: [
          Row(
            children: [
              Icon(
                Icons.lightbulb_outline,
                color: Colors.green.shade700,
                size: 16,
              ),
              const SizedBox(width: 8),
              Text(
                'Recommendations',
                style: Theme.of(context).textTheme.labelMedium?.copyWith(
                  fontWeight: FontWeight.bold,
                  color: Colors.green.shade700,
                ),
              ),
            ],
          ),
          const SizedBox(height: 8),
          
          // Method-specific recommendations
          if (widget.selectedMethod == 'HF')
            Text(
              '• HF: Fast, good starting point for geometry optimization\n'
              '• Try 6-31G* for better results with minimal cost increase',
              style: Theme.of(context).textTheme.bodySmall,
            )
          else if (widget.selectedMethod == 'B3LYP')
            Text(
              '• B3LYP: Excellent for most organic molecules\n'
              '• 6-31G* or cc-pVDZ recommended for production calculations',
              style: Theme.of(context).textTheme.bodySmall,
            )
          else if (widget.selectedMethod == 'MP2')
            Text(
              '• MP2: Includes correlation, more expensive than DFT\n'
              '• Consider cc-pVDZ or larger basis sets',
              style: Theme.of(context).textTheme.bodySmall,
            )
          else
            Text(
              '• Start with HF/6-31G for quick tests\n'
              '• Use B3LYP/6-31G* for production calculations',
              style: Theme.of(context).textTheme.bodySmall,
            ),
        ],
      ),
    );
  }

  @override
  Widget build(BuildContext context) {
    return Column(
      crossAxisAlignment: CrossAxisAlignment.start,
      children: [
        // Header
        Text(
          'Calculation Settings',
          style: Theme.of(context).textTheme.titleLarge,
        ),
        
        const SizedBox(height: 16),
        
        Expanded(
          child: SingleChildScrollView(
            child: Column(
              crossAxisAlignment: CrossAxisAlignment.start,
              children: [
                // Method selection
                _buildMethodSelector(),
                
                const SizedBox(height: 16),
                
                // Basis set selection
                _buildBasisSetSelector(),
                
                const SizedBox(height: 16),
                
                // Basic settings
                _buildBasicSettings(),
                
                const SizedBox(height: 16),
                
                // Advanced settings toggle
                Row(
                  children: [
                    Text(
                      'Advanced Settings',
                      style: Theme.of(context).textTheme.titleSmall,
                    ),
                    const Spacer(),
                    Switch(
                      value: _showAdvancedSettings,
                      onChanged: (value) {
                        setState(() {
                          _showAdvancedSettings = value;
                        });
                      },
                    ),
                  ],
                ),
                
                if (_showAdvancedSettings) ...[
                  const SizedBox(height: 16),
                  _buildAdvancedSettings(),
                ],
                
                const SizedBox(height: 16),
                
                // Recommendations
                _buildRecommendations(),
              ],
            ),
          ),
        ),
      ],
    );
  }
}