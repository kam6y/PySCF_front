import 'package:flutter/material.dart';
import 'package:fl_chart/fl_chart.dart';
import '../models/calculation_result.dart';
import '../models/molecule.dart';

class ResultsDisplayWidget extends StatefulWidget {
  final CalculationResult result;
  final MoleculeData? molecule;

  const ResultsDisplayWidget({
    super.key,
    required this.result,
    this.molecule,
  });

  @override
  State<ResultsDisplayWidget> createState() => _ResultsDisplayWidgetState();
}

class _ResultsDisplayWidgetState extends State<ResultsDisplayWidget>
    with SingleTickerProviderStateMixin {
  late TabController _tabController;

  @override
  void initState() {
    super.initState();
    _tabController = TabController(
      length: widget.result.success ? 3 : 1,
      vsync: this,
    );
  }

  @override
  void dispose() {
    _tabController.dispose();
    super.dispose();
  }

  Widget _buildSummaryTab() {
    return Padding(
      padding: const EdgeInsets.all(16.0),
      child: SingleChildScrollView(
        child: Column(
          crossAxisAlignment: CrossAxisAlignment.start,
          children: [
            // Status indicator
            Container(
              width: double.infinity,
              padding: const EdgeInsets.all(16),
              decoration: BoxDecoration(
                color: widget.result.success
                    ? Colors.green.withOpacity(0.1)
                    : Colors.red.withOpacity(0.1),
                border: Border.all(
                  color: widget.result.success ? Colors.green : Colors.red,
                ),
                borderRadius: BorderRadius.circular(8),
              ),
              child: Column(
                crossAxisAlignment: CrossAxisAlignment.start,
                children: [
                  Row(
                    children: [
                      Icon(
                        widget.result.success
                            ? Icons.check_circle
                            : Icons.error,
                        color: widget.result.success ? Colors.green : Colors.red,
                        size: 20,
                      ),
                      const SizedBox(width: 8),
                      Text(
                        widget.result.success
                            ? 'Calculation Completed'
                            : 'Calculation Failed',
                        style: TextStyle(
                          color: widget.result.success ? Colors.green : Colors.red,
                          fontWeight: FontWeight.bold,
                          fontSize: 16,
                        ),
                      ),
                    ],
                  ),
                  const SizedBox(height: 8),
                  Text(
                    widget.result.message,
                    style: TextStyle(
                      color: widget.result.success ? Colors.green : Colors.red,
                    ),
                  ),
                ],
              ),
            ),

            if (!widget.result.success) ...[
              const SizedBox(height: 16),
              // Error details
              Card(
                child: Padding(
                  padding: const EdgeInsets.all(16),
                  child: Column(
                    crossAxisAlignment: CrossAxisAlignment.start,
                    children: [
                      Text(
                        'Error Details',
                        style: Theme.of(context).textTheme.titleMedium,
                      ),
                      const SizedBox(height: 8),
                      Text(widget.result.message),
                      if (widget.result.additionalData != null) ...[
                        const SizedBox(height: 16),
                        Text(
                          'Additional Information:',
                          style: Theme.of(context).textTheme.titleSmall,
                        ),
                        const SizedBox(height: 4),
                        Text(
                          widget.result.additionalData.toString(),
                          style: const TextStyle(fontFamily: 'monospace'),
                        ),
                      ],
                    ],
                  ),
                ),
              ),
            ] else ...[
              const SizedBox(height: 16),

              // Main results
              Card(
                child: Padding(
                  padding: const EdgeInsets.all(16),
                  child: Column(
                    crossAxisAlignment: CrossAxisAlignment.start,
                    children: [
                      Text(
                        'Calculation Results',
                        style: Theme.of(context).textTheme.titleMedium,
                      ),
                      const SizedBox(height: 16),

                      // Energy information
                      _buildEnergySection(),

                      const SizedBox(height: 16),

                      // Timing information
                      Row(
                        children: [
                          Icon(
                            Icons.timer,
                            size: 16,
                            color: Theme.of(context).hintColor,
                          ),
                          const SizedBox(width: 8),
                          Text(
                            'Calculation Time: ${widget.result.formattedCalculationTime}',
                            style: Theme.of(context).textTheme.bodyMedium,
                          ),
                        ],
                      ),

                      if (widget.result.completedAt != null) ...[
                        const SizedBox(height: 8),
                        Row(
                          children: [
                            Icon(
                              Icons.schedule,
                              size: 16,
                              color: Theme.of(context).hintColor,
                            ),
                            const SizedBox(width: 8),
                            Text(
                              'Completed: ${widget.result.completedAt!.toLocal().toString().split('.')[0]}',
                              style: Theme.of(context).textTheme.bodyMedium,
                            ),
                          ],
                        ),
                      ],
                    ],
                  ),
                ),
              ),

              // Molecule information
              if (widget.molecule != null) ...[
                const SizedBox(height: 16),
                Card(
                  child: Padding(
                    padding: const EdgeInsets.all(16),
                    child: Column(
                      crossAxisAlignment: CrossAxisAlignment.start,
                      children: [
                        Text(
                          'Molecule Information',
                          style: Theme.of(context).textTheme.titleMedium,
                        ),
                        const SizedBox(height: 16),
                        _buildMoleculeInfo(),
                      ],
                    ),
                  ),
                ),
              ],
            ],
          ],
        ),
      ),
    );
  }

  Widget _buildEnergySection() {
    if (!widget.result.hasEnergyData) {
      return Column(
        crossAxisAlignment: CrossAxisAlignment.start,
        children: [
          Row(
            crossAxisAlignment: CrossAxisAlignment.baseline,
            textBaseline: TextBaseline.alphabetic,
            children: [
              Text(
                '${widget.result.energy.toStringAsFixed(6)}',
                style: Theme.of(context).textTheme.headlineSmall?.copyWith(
                  fontWeight: FontWeight.bold,
                  fontFamily: 'monospace',
                ),
              ),
              const SizedBox(width: 8),
              Text(
                'Hartree',
                style: Theme.of(context).textTheme.bodyMedium,
              ),
            ],
          ),
          const SizedBox(height: 8),
          Text(
            'Total Energy',
            style: Theme.of(context).textTheme.titleSmall,
          ),
        ],
      );
    }

    final energyData = widget.result.energyData!;
    return Column(
      crossAxisAlignment: CrossAxisAlignment.start,
      children: [
        // Total energy
        Row(
          crossAxisAlignment: CrossAxisAlignment.baseline,
          textBaseline: TextBaseline.alphabetic,
          children: [
            Text(
              '${energyData.totalEnergy.toStringAsFixed(6)}',
              style: Theme.of(context).textTheme.headlineSmall?.copyWith(
                fontWeight: FontWeight.bold,
                fontFamily: 'monospace',
              ),
            ),
            const SizedBox(width: 8),
            Text(
              'Hartree',
              style: Theme.of(context).textTheme.bodyMedium,
            ),
          ],
        ),
        const SizedBox(height: 4),
        Text(
          'Total Energy',
          style: Theme.of(context).textTheme.titleSmall,
        ),

        const SizedBox(height: 16),

        // Energy in different units
        Row(
          children: [
            Expanded(
              child: Column(
                crossAxisAlignment: CrossAxisAlignment.start,
                children: [
                  Text(
                    '${energyData.totalEnergyEv.toStringAsFixed(2)} eV',
                    style: const TextStyle(fontFamily: 'monospace'),
                  ),
                  Text(
                    'Electron Volts',
                    style: Theme.of(context).textTheme.bodySmall,
                  ),
                ],
              ),
            ),
            Expanded(
              child: Column(
                crossAxisAlignment: CrossAxisAlignment.start,
                children: [
                  Text(
                    '${energyData.totalEnergyKcalMol.toStringAsFixed(1)} kcal/mol',
                    style: const TextStyle(fontFamily: 'monospace'),
                  ),
                  Text(
                    'kcal/mol',
                    style: Theme.of(context).textTheme.bodySmall,
                  ),
                ],
              ),
            ),
          ],
        ),

        // HOMO-LUMO gap if available
        if (energyData.homoLumoGap != null) ...[
          const SizedBox(height: 16),
          Container(
            width: double.infinity,
            padding: const EdgeInsets.all(12),
            decoration: BoxDecoration(
              color: Theme.of(context).primaryColor.withOpacity(0.1),
              borderRadius: BorderRadius.circular(4),
            ),
            child: Column(
              crossAxisAlignment: CrossAxisAlignment.start,
              children: [
                Text(
                  'HOMO-LUMO Gap',
                  style: Theme.of(context).textTheme.titleSmall,
                ),
                const SizedBox(height: 4),
                Text(
                  '${energyData.homoLumoGap!.toStringAsFixed(4)} Hartree',
                  style: const TextStyle(fontFamily: 'monospace'),
                ),
                Text(
                  '${energyData.homoLumoGapEv!.toStringAsFixed(2)} eV',
                  style: const TextStyle(fontFamily: 'monospace'),
                ),
              ],
            ),
          ),
        ],
      ],
    );
  }

  Widget _buildMoleculeInfo() {
    final molecule = widget.molecule!;
    return Column(
      crossAxisAlignment: CrossAxisAlignment.start,
      children: [
        Row(
          children: [
            Expanded(
              child: Column(
                crossAxisAlignment: CrossAxisAlignment.start,
                children: [
                  Text(
                    'Name',
                    style: Theme.of(context).textTheme.bodySmall,
                  ),
                  Text(
                    molecule.name,
                    style: Theme.of(context).textTheme.bodyMedium?.copyWith(
                      fontWeight: FontWeight.w500,
                    ),
                  ),
                ],
              ),
            ),
            Expanded(
              child: Column(
                crossAxisAlignment: CrossAxisAlignment.start,
                children: [
                  Text(
                    'Formula',
                    style: Theme.of(context).textTheme.bodySmall,
                  ),
                  Text(
                    molecule.formula,
                    style: Theme.of(context).textTheme.bodyMedium?.copyWith(
                      fontWeight: FontWeight.w500,
                    ),
                  ),
                ],
              ),
            ),
          ],
        ),
        const SizedBox(height: 12),
        Row(
          children: [
            Expanded(
              child: Column(
                crossAxisAlignment: CrossAxisAlignment.start,
                children: [
                  Text(
                    'Charge',
                    style: Theme.of(context).textTheme.bodySmall,
                  ),
                  Text(
                    '${molecule.charge}',
                    style: Theme.of(context).textTheme.bodyMedium?.copyWith(
                      fontWeight: FontWeight.w500,
                    ),
                  ),
                ],
              ),
            ),
            Expanded(
              child: Column(
                crossAxisAlignment: CrossAxisAlignment.start,
                children: [
                  Text(
                    'Multiplicity',
                    style: Theme.of(context).textTheme.bodySmall,
                  ),
                  Text(
                    '${molecule.multiplicity}',
                    style: Theme.of(context).textTheme.bodyMedium?.copyWith(
                      fontWeight: FontWeight.w500,
                    ),
                  ),
                ],
              ),
            ),
            Expanded(
              child: Column(
                crossAxisAlignment: CrossAxisAlignment.start,
                children: [
                  Text(
                    'Atoms',
                    style: Theme.of(context).textTheme.bodySmall,
                  ),
                  Text(
                    '${molecule.atoms.length}',
                    style: Theme.of(context).textTheme.bodyMedium?.copyWith(
                      fontWeight: FontWeight.w500,
                    ),
                  ),
                ],
              ),
            ),
          ],
        ),
      ],
    );
  }

  Widget _buildEnergyTab() {
    if (!widget.result.hasOrbitalData) {
      return const Center(
        child: Text('No orbital data available'),
      );
    }

    return Padding(
      padding: const EdgeInsets.all(16.0),
      child: Column(
        crossAxisAlignment: CrossAxisAlignment.start,
        children: [
          Text(
            'Molecular Orbitals',
            style: Theme.of(context).textTheme.titleLarge,
          ),
          const SizedBox(height: 16),
          
          Expanded(
            child: _buildOrbitalChart(),
          ),
        ],
      ),
    );
  }

  Widget _buildOrbitalChart() {
    final orbitals = widget.result.orbitals!;
    final occupiedOrbitals = widget.result.occupiedOrbitals;
    final virtualOrbitals = widget.result.virtualOrbitals.take(10).toList(); // Show first 10 virtual

    return Column(
      children: [
        // HOMO-LUMO information
        if (widget.result.homoOrbital != null && widget.result.lumoOrbital != null)
          Container(
            width: double.infinity,
            padding: const EdgeInsets.all(12),
            margin: const EdgeInsets.only(bottom: 16),
            decoration: BoxDecoration(
              color: Theme.of(context).primaryColor.withOpacity(0.1),
              borderRadius: BorderRadius.circular(8),
            ),
            child: Row(
              children: [
                Expanded(
                  child: Column(
                    crossAxisAlignment: CrossAxisAlignment.start,
                    children: [
                      Text(
                        'HOMO',
                        style: Theme.of(context).textTheme.titleSmall,
                      ),
                      Text(
                        '${widget.result.homoOrbital!.energy.toStringAsFixed(4)} Ha',
                        style: const TextStyle(fontFamily: 'monospace'),
                      ),
                    ],
                  ),
                ),
                Expanded(
                  child: Column(
                    crossAxisAlignment: CrossAxisAlignment.start,
                    children: [
                      Text(
                        'LUMO',
                        style: Theme.of(context).textTheme.titleSmall,
                      ),
                      Text(
                        '${widget.result.lumoOrbital!.energy.toStringAsFixed(4)} Ha',
                        style: const TextStyle(fontFamily: 'monospace'),
                      ),
                    ],
                  ),
                ),
                if (widget.result.energyData?.homoLumoGap != null)
                  Expanded(
                    child: Column(
                      crossAxisAlignment: CrossAxisAlignment.start,
                      children: [
                        Text(
                          'Gap',
                          style: Theme.of(context).textTheme.titleSmall,
                        ),
                        Text(
                          '${widget.result.energyData!.homoLumoGapEv!.toStringAsFixed(2)} eV',
                          style: const TextStyle(fontFamily: 'monospace'),
                        ),
                      ],
                    ),
                  ),
              ],
            ),
          ),

        // Orbital energy diagram (simplified)
        Expanded(
          child: ListView.builder(
            itemCount: occupiedOrbitals.length + virtualOrbitals.length,
            itemBuilder: (context, index) {
              late OrbitalData orbital;
              bool isOccupied;
              
              if (index < occupiedOrbitals.length) {
                orbital = occupiedOrbitals[occupiedOrbitals.length - 1 - index]; // Reverse order
                isOccupied = true;
              } else {
                orbital = virtualOrbitals[index - occupiedOrbitals.length];
                isOccupied = false;
              }

              final isHomo = orbital.type == 'homo';
              final isLumo = orbital.type == 'lumo';

              return Container(
                margin: const EdgeInsets.only(bottom: 4),
                child: Row(
                  children: [
                    // Orbital index
                    SizedBox(
                      width: 40,
                      child: Text(
                        '${orbital.index}',
                        style: const TextStyle(fontFamily: 'monospace'),
                      ),
                    ),
                    
                    // Energy bar
                    Expanded(
                      child: Container(
                        height: 24,
                        decoration: BoxDecoration(
                          color: isHomo
                              ? Colors.red.withOpacity(0.7)
                              : isLumo
                                  ? Colors.blue.withOpacity(0.7)
                                  : isOccupied
                                      ? Colors.green.withOpacity(0.5)
                                      : Colors.grey.withOpacity(0.3),
                          borderRadius: BorderRadius.circular(4),
                        ),
                        child: Center(
                          child: Text(
                            isHomo
                                ? 'HOMO'
                                : isLumo
                                    ? 'LUMO'
                                    : isOccupied
                                        ? 'Occ'
                                        : 'Virt',
                            style: const TextStyle(
                              fontSize: 10,
                              fontWeight: FontWeight.bold,
                              color: Colors.white,
                            ),
                          ),
                        ),
                      ),
                    ),
                    
                    const SizedBox(width: 8),
                    
                    // Energy value
                    SizedBox(
                      width: 100,
                      child: Text(
                        '${orbital.energy.toStringAsFixed(4)}',
                        style: const TextStyle(fontFamily: 'monospace'),
                        textAlign: TextAlign.right,
                      ),
                    ),
                  ],
                ),
              );
            },
          ),
        ),
      ],
    );
  }

  Widget _buildPropertiesTab() {
    return Padding(
      padding: const EdgeInsets.all(16.0),
      child: SingleChildScrollView(
        child: Column(
          crossAxisAlignment: CrossAxisAlignment.start,
          children: [
            Text(
              'Molecular Properties',
              style: Theme.of(context).textTheme.titleLarge,
            ),
            const SizedBox(height: 16),
            
            // Placeholder for molecular properties
            Card(
              child: Padding(
                padding: const EdgeInsets.all(16),
                child: Column(
                  crossAxisAlignment: CrossAxisAlignment.start,
                  children: [
                    Text(
                      'Electronic Properties',
                      style: Theme.of(context).textTheme.titleMedium,
                    ),
                    const SizedBox(height: 16),
                    
                    // Properties coming in Phase 2
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
                            'Coming in Phase 2:',
                            style: Theme.of(context).textTheme.titleSmall?.copyWith(
                              fontWeight: FontWeight.bold,
                            ),
                          ),
                          const SizedBox(height: 8),
                          Text(
                            '• Mulliken population analysis\n'
                            '• Dipole moment\n'
                            '• Vibrational frequencies\n'
                            '• IR spectrum\n'
                            '• UV-Vis spectrum\n'
                            '• Thermodynamic properties',
                            style: Theme.of(context).textTheme.bodySmall,
                          ),
                        ],
                      ),
                    ),
                  ],
                ),
              ),
            ),
          ],
        ),
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
          'Calculation Results',
          style: Theme.of(context).textTheme.titleLarge,
        ),
        
        const SizedBox(height: 16),
        
        // Tab bar (only show if calculation succeeded)
        if (widget.result.success)
          TabBar(
            controller: _tabController,
            tabs: const [
              Tab(text: 'Summary'),
              Tab(text: 'Orbitals'),
              Tab(text: 'Properties'),
            ],
          ),
        
        // Tab content
        Expanded(
          child: widget.result.success
              ? TabBarView(
                  controller: _tabController,
                  children: [
                    _buildSummaryTab(),
                    _buildEnergyTab(),
                    _buildPropertiesTab(),
                  ],
                )
              : _buildSummaryTab(),
        ),
      ],
    );
  }
}