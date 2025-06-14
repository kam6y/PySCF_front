import 'package:flutter/material.dart';
import 'package:file_picker/file_picker.dart';
import '../models/molecule.dart';

class MoleculeInputWidget extends StatefulWidget {
  final Function(MoleculeData) onMoleculeChanged;

  const MoleculeInputWidget({
    super.key,
    required this.onMoleculeChanged,
  });

  @override
  State<MoleculeInputWidget> createState() => _MoleculeInputWidgetState();
}

class _MoleculeInputWidgetState extends State<MoleculeInputWidget>
    with SingleTickerProviderStateMixin {
  late TabController _tabController;
  final TextEditingController _nameController = TextEditingController();
  final TextEditingController _xyzController = TextEditingController();
  String? _errorMessage;
  MoleculeData? _currentMolecule;

  @override
  void initState() {
    super.initState();
    _tabController = TabController(length: 3, vsync: this);
    _loadSampleMolecule();
  }

  @override
  void dispose() {
    _tabController.dispose();
    _nameController.dispose();
    _xyzController.dispose();
    super.dispose();
  }

  void _loadSampleMolecule() {
    // Load a water molecule as default
    const waterXyz = '''2
Water molecule
O    0.000000    0.000000    0.000000
H    0.757000    0.586000    0.000000
H   -0.757000    0.586000    0.000000''';
    
    _nameController.text = 'Water';
    _xyzController.text = waterXyz;
    _parseMolecule();
  }

  void _parseMolecule() {
    setState(() {
      _errorMessage = null;
    });

    try {
      final xyzContent = _xyzController.text.trim();
      if (xyzContent.isEmpty) {
        setState(() {
          _currentMolecule = null;
        });
        return;
      }

      final molecule = MoleculeData.fromXyz(
        xyzContent,
        name: _nameController.text.isNotEmpty ? _nameController.text : null,
      );

      setState(() {
        _currentMolecule = molecule;
        _errorMessage = null;
      });

      widget.onMoleculeChanged(molecule);
    } catch (e) {
      setState(() {
        _errorMessage = 'Error parsing molecule: $e';
        _currentMolecule = null;
      });
    }
  }

  Future<void> _loadFromFile() async {
    try {
      FilePickerResult? result = await FilePicker.platform.pickFiles(
        type: FileType.custom,
        allowedExtensions: ['xyz'],
        withData: true,
      );

      if (result != null && result.files.single.bytes != null) {
        final content = String.fromCharCodes(result.files.single.bytes!);
        final fileName = result.files.single.name;
        
        setState(() {
          _xyzController.text = content;
          _nameController.text = fileName.replaceAll('.xyz', '');
        });
        
        _parseMolecule();
      }
    } catch (e) {
      setState(() {
        _errorMessage = 'Error loading file: $e';
      });
    }
  }

  void _addAtom() {
    final currentText = _xyzController.text.trim();
    final lines = currentText.split('\n');
    
    if (lines.length < 2) {
      // Initialize with basic structure
      _xyzController.text = '''1
New molecule
H    0.000000    0.000000    0.000000''';
    } else {
      // Parse current atom count and increment
      try {
        final atomCount = int.parse(lines[0]);
        final newAtomCount = atomCount + 1;
        
        // Update atom count
        lines[0] = newAtomCount.toString();
        
        // Add new hydrogen atom
        lines.add('H    0.000000    0.000000    0.000000');
        
        _xyzController.text = lines.join('\n');
      } catch (e) {
        setState(() {
          _errorMessage = 'Error adding atom: $e';
        });
        return;
      }
    }
    
    _parseMolecule();
  }

  void _clearMolecule() {
    setState(() {
      _nameController.clear();
      _xyzController.clear();
      _currentMolecule = null;
      _errorMessage = null;
    });
  }

  Widget _buildFileInputTab() {
    return Padding(
      padding: const EdgeInsets.all(16.0),
      child: Column(
        crossAxisAlignment: CrossAxisAlignment.start,
        children: [
          Text(
            'File Input',
            style: Theme.of(context).textTheme.titleMedium,
          ),
          const SizedBox(height: 16),
          
          // File picker button
          SizedBox(
            width: double.infinity,
            child: ElevatedButton.icon(
              onPressed: _loadFromFile,
              icon: const Icon(Icons.folder_open),
              label: const Text('Load XYZ File'),
            ),
          ),
          
          const SizedBox(height: 16),
          
          // Drag and drop area (placeholder)
          Container(
            width: double.infinity,
            height: 120,
            decoration: BoxDecoration(
              border: Border.all(
                color: Theme.of(context).dividerColor,
                style: BorderStyle.solid,
                width: 2,
              ),
              borderRadius: BorderRadius.circular(8),
            ),
            child: Column(
              mainAxisAlignment: MainAxisAlignment.center,
              children: [
                Icon(
                  Icons.cloud_upload_outlined,
                  size: 48,
                  color: Theme.of(context).hintColor,
                ),
                const SizedBox(height: 8),
                Text(
                  'Drag & Drop XYZ files here',
                  style: TextStyle(
                    color: Theme.of(context).hintColor,
                    fontSize: 16,
                  ),
                ),
                Text(
                  '(Feature coming soon)',
                  style: TextStyle(
                    color: Theme.of(context).hintColor,
                    fontSize: 12,
                  ),
                ),
              ],
            ),
          ),
          
          const SizedBox(height: 16),
          
          // Sample molecules
          Text(
            'Sample Molecules',
            style: Theme.of(context).textTheme.titleSmall,
          ),
          const SizedBox(height: 8),
          Wrap(
            spacing: 8,
            children: [
              _buildSampleButton('Water', '''2
Water molecule
O    0.000000    0.000000    0.000000
H    0.757000    0.586000    0.000000
H   -0.757000    0.586000    0.000000'''),
              _buildSampleButton('Methane', '''5
Methane molecule
C    0.000000    0.000000    0.000000
H    1.089000    0.000000    0.000000
H   -0.363000    1.026000    0.000000
H   -0.363000   -0.513000    0.889000
H   -0.363000   -0.513000   -0.889000'''),
              _buildSampleButton('Ammonia', '''4
Ammonia molecule
N    0.000000    0.000000    0.000000
H    0.000000    1.012000    0.000000
H    0.876000   -0.506000    0.000000
H   -0.876000   -0.506000    0.000000'''),
            ],
          ),
        ],
      ),
    );
  }

  Widget _buildSampleButton(String name, String xyzContent) {
    return ElevatedButton(
      onPressed: () {
        setState(() {
          _nameController.text = name;
          _xyzController.text = xyzContent;
        });
        _parseMolecule();
      },
      child: Text(name),
    );
  }

  Widget _buildDirectInputTab() {
    return Padding(
      padding: const EdgeInsets.all(16.0),
      child: Column(
        crossAxisAlignment: CrossAxisAlignment.start,
        children: [
          // Molecule name
          TextField(
            controller: _nameController,
            decoration: const InputDecoration(
              labelText: 'Molecule Name',
              border: OutlineInputBorder(),
            ),
            onChanged: (_) => _parseMolecule(),
          ),
          
          const SizedBox(height: 16),
          
          // XYZ coordinates input
          Text(
            'XYZ Coordinates',
            style: Theme.of(context).textTheme.titleSmall,
          ),
          const SizedBox(height: 8),
          
          Expanded(
            child: TextField(
              controller: _xyzController,
              decoration: const InputDecoration(
                border: OutlineInputBorder(),
                hintText: '''2
Molecule name
O  0.000  0.000  0.000
H  0.757  0.586  0.000''',
                contentPadding: EdgeInsets.all(12),
              ),
              style: const TextStyle(fontFamily: 'monospace'),
              maxLines: null,
              expands: true,
              textAlignVertical: TextAlignVertical.top,
              onChanged: (_) => _parseMolecule(),
            ),
          ),
          
          const SizedBox(height: 16),
          
          // Action buttons
          Row(
            children: [
              ElevatedButton.icon(
                onPressed: _addAtom,
                icon: const Icon(Icons.add),
                label: const Text('Add Atom'),
              ),
              const SizedBox(width: 8),
              ElevatedButton.icon(
                onPressed: _clearMolecule,
                icon: const Icon(Icons.clear),
                label: const Text('Clear'),
              ),
              const Spacer(),
              ElevatedButton(
                onPressed: _parseMolecule,
                child: const Text('Parse'),
              ),
            ],
          ),
        ],
      ),
    );
  }

  Widget _buildDatabaseTab() {
    return Padding(
      padding: const EdgeInsets.all(16.0),
      child: Column(
        crossAxisAlignment: CrossAxisAlignment.start,
        children: [
          Text(
            'Database Search',
            style: Theme.of(context).textTheme.titleMedium,
          ),
          const SizedBox(height: 16),
          
          // Search field
          TextField(
            decoration: const InputDecoration(
              labelText: 'Search molecule by name or formula',
              border: OutlineInputBorder(),
              prefixIcon: Icon(Icons.search),
              hintText: 'e.g., water, H2O, benzene, C6H6',
            ),
            onSubmitted: (value) {
              // TODO: Implement database search
            },
          ),
          
          const SizedBox(height: 16),
          
          // Placeholder for search results
          Expanded(
            child: Container(
              width: double.infinity,
              decoration: BoxDecoration(
                border: Border.all(color: Theme.of(context).dividerColor),
                borderRadius: BorderRadius.circular(8),
              ),
              child: Column(
                mainAxisAlignment: MainAxisAlignment.center,
                children: [
                  Icon(
                    Icons.storage_outlined,
                    size: 64,
                    color: Theme.of(context).hintColor,
                  ),
                  const SizedBox(height: 16),
                  Text(
                    'Database Search',
                    style: Theme.of(context).textTheme.titleMedium?.copyWith(
                      color: Theme.of(context).hintColor,
                    ),
                  ),
                  const SizedBox(height: 8),
                  Text(
                    'PubChem integration coming soon',
                    style: TextStyle(
                      color: Theme.of(context).hintColor,
                    ),
                  ),
                ],
              ),
            ),
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
          'Molecule Input',
          style: Theme.of(context).textTheme.titleLarge,
        ),
        
        const SizedBox(height: 16),
        
        // Tab bar
        TabBar(
          controller: _tabController,
          tabs: const [
            Tab(text: 'File'),
            Tab(text: 'Direct'),
            Tab(text: 'Database'),
          ],
        ),
        
        const SizedBox(height: 8),
        
        // Tab content with fixed height
        Expanded(
          child: TabBarView(
            controller: _tabController,
            children: [
              _buildFileInputTab(),
              _buildDirectInputTab(),
              _buildDatabaseTab(),
            ],
          ),
        ),
        
        // Status messages at bottom with constrained height
        Container(
          constraints: const BoxConstraints(maxHeight: 120),
          child: SingleChildScrollView(
            child: Column(
              children: [
                // Error message
                if (_errorMessage != null)
                  Container(
                    width: double.infinity,
                    padding: const EdgeInsets.all(8),
                    margin: const EdgeInsets.only(top: 4),
                    decoration: BoxDecoration(
                      color: Colors.red.withOpacity(0.1),
                      border: Border.all(color: Colors.red),
                      borderRadius: BorderRadius.circular(4),
                    ),
                    child: Text(
                      _errorMessage!,
                      style: const TextStyle(color: Colors.red, fontSize: 12),
                      maxLines: 2,
                      overflow: TextOverflow.ellipsis,
                    ),
                  ),
                
                // Molecule info
                if (_currentMolecule != null)
                  Container(
                    width: double.infinity,
                    padding: const EdgeInsets.all(8),
                    margin: const EdgeInsets.only(top: 4),
                    decoration: BoxDecoration(
                      color: Colors.green.withOpacity(0.1),
                      border: Border.all(color: Colors.green),
                      borderRadius: BorderRadius.circular(4),
                    ),
                    child: Column(
                      crossAxisAlignment: CrossAxisAlignment.start,
                      children: [
                        Text(
                          'Loaded: ${_currentMolecule!.name}',
                          style: const TextStyle(
                            color: Colors.green,
                            fontWeight: FontWeight.bold,
                            fontSize: 12,
                          ),
                        ),
                        Text(
                          'Formula: ${_currentMolecule!.formula}',
                          style: const TextStyle(color: Colors.green, fontSize: 11),
                        ),
                        Text(
                          'Atoms: ${_currentMolecule!.atoms.length}',
                          style: const TextStyle(color: Colors.green, fontSize: 11),
                        ),
                      ],
                    ),
                  ),
              ],
            ),
          ),
        ),
      ],
    );
  }
}