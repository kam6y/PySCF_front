import 'package:flutter/material.dart';
import 'dart:math' as math;
import '../models/molecule.dart';

class MoleculeVisualizerWidget extends StatefulWidget {
  final MoleculeData? molecule;

  const MoleculeVisualizerWidget({
    super.key,
    this.molecule,
  });

  @override
  State<MoleculeVisualizerWidget> createState() => _MoleculeVisualizerWidgetState();
}

class _MoleculeVisualizerWidgetState extends State<MoleculeVisualizerWidget>
    with TickerProviderStateMixin {
  late AnimationController _rotationController;
  double _rotationX = 0.0;
  double _rotationY = 0.0;
  double _zoom = 1.0;
  bool _showBonds = true;
  String _renderMode = 'ball_stick'; // 'ball_stick', 'space_filling', 'wireframe'

  // Atom colors (CPK coloring scheme)
  static const Map<String, Color> _atomColors = {
    'H': Color(0xFFFFFFFF), // White
    'C': Color(0xFF909090), // Gray
    'N': Color(0xFF3050F8), // Blue
    'O': Color(0xFFFF0D0D), // Red
    'F': Color(0xFF90E050), // Green
    'P': Color(0xFFFF8000), // Orange
    'S': Color(0xFFFFFF30), // Yellow
    'Cl': Color(0xFF1FF01F), // Green
    'Br': Color(0xFFA62929), // Brown
    'I': Color(0xFF940094), // Purple
    'He': Color(0xFFD9FFFF), // Cyan
    'Li': Color(0xFFCC80FF), // Violet
    'Be': Color(0xFFC2FF00), // Green
    'B': Color(0xFFFFB5B5), // Pink
    'Ne': Color(0xFFB3E3F5), // Light blue
    'Na': Color(0xFFAB5CF2), // Violet
    'Mg': Color(0xFF8AFF00), // Green
    'Al': Color(0xFFBFA6A6), // Gray
    'Si': Color(0xFFF0C8A0), // Tan
    'Ar': Color(0xFF80D1E3), // Cyan
    'K': Color(0xFF8F40D4), // Violet
    'Ca': Color(0xFF3DFF00), // Green
  };

  // Van der Waals radii in Angstroms
  static const Map<String, double> _atomRadii = {
    'H': 1.20,
    'C': 1.70,
    'N': 1.55,
    'O': 1.52,
    'F': 1.47,
    'P': 1.80,
    'S': 1.80,
    'Cl': 1.75,
    'Br': 1.85,
    'I': 1.98,
    'He': 1.40,
    'Li': 1.82,
    'Be': 1.53,
    'B': 1.92,
    'Ne': 1.54,
    'Na': 2.27,
    'Mg': 1.73,
    'Al': 1.84,
    'Si': 2.10,
    'Ar': 1.88,
    'K': 2.75,
    'Ca': 2.31,
  };

  @override
  void initState() {
    super.initState();
    _rotationController = AnimationController(
      duration: const Duration(seconds: 10),
      vsync: this,
    );
    
    // Auto-rotation (can be toggled)
    // _rotationController.repeat();
  }

  @override
  void dispose() {
    _rotationController.dispose();
    super.dispose();
  }

  Color _getAtomColor(String symbol) {
    return _atomColors[symbol] ?? const Color(0xFFFF1493); // Default: deep pink
  }

  double _getAtomRadius(String symbol) {
    return _atomRadii[symbol] ?? 1.5; // Default radius
  }

  // Calculate bond distance between two atoms
  double _calculateDistance(AtomData atom1, AtomData atom2) {
    final dx = atom1.x - atom2.x;
    final dy = atom1.y - atom2.y;
    final dz = atom1.z - atom2.z;
    return math.sqrt(dx * dx + dy * dy + dz * dz);
  }

  // Simple bond detection based on distance
  List<List<int>> _detectBonds(List<AtomData> atoms) {
    final bonds = <List<int>>[];
    
    for (int i = 0; i < atoms.length; i++) {
      for (int j = i + 1; j < atoms.length; j++) {
        final distance = _calculateDistance(atoms[i], atoms[j]);
        final atom1Radius = _getAtomRadius(atoms[i].symbol);
        final atom2Radius = _getAtomRadius(atoms[j].symbol);
        
        // Bond if distance is less than sum of covalent radii * 1.3
        final bondThreshold = (atom1Radius + atom2Radius) * 0.6; // Adjusted for covalent radii
        
        if (distance < bondThreshold) {
          bonds.add([i, j]);
        }
      }
    }
    
    return bonds;
  }

  // Transform 3D coordinates to 2D screen coordinates
  Offset _project3DTo2D(AtomData atom, Size size) {
    // Apply rotations
    final cosX = math.cos(_rotationX);
    final sinX = math.sin(_rotationX);
    final cosY = math.cos(_rotationY);
    final sinY = math.sin(_rotationY);
    
    // Rotate around Y axis
    final x1 = atom.x * cosY - atom.z * sinY;
    final z1 = atom.x * sinY + atom.z * cosY;
    
    // Rotate around X axis
    final y2 = atom.y * cosX - z1 * sinX;
    final z2 = atom.y * sinX + z1 * cosX;
    
    // Apply zoom and center
    final scale = _zoom * math.min(size.width, size.height) * 0.3;
    final centerX = size.width / 2;
    final centerY = size.height / 2;
    
    return Offset(
      centerX + x1 * scale,
      centerY - y2 * scale, // Flip Y axis
    );
  }

  Widget _buildMoleculeViewer() {
    if (widget.molecule == null || widget.molecule!.atoms.isEmpty) {
      return _buildEmptyState();
    }

    return GestureDetector(
      onPanUpdate: (details) {
        setState(() {
          _rotationY += details.delta.dx * 0.01;
          _rotationX += details.delta.dy * 0.01;
        });
      },
      child: CustomPaint(
        painter: MoleculePainter(
          molecule: widget.molecule!,
          rotationX: _rotationX,
          rotationY: _rotationY,
          zoom: _zoom,
          showBonds: _showBonds,
          renderMode: _renderMode,
        ),
        size: Size.infinite,
      ),
    );
  }

  Widget _buildEmptyState() {
    return Column(
      mainAxisAlignment: MainAxisAlignment.center,
      children: [
        Icon(
          Icons.science_outlined,
          size: 64,
          color: Theme.of(context).hintColor,
        ),
        const SizedBox(height: 16),
        Text(
          '3D Molecule Viewer',
          style: Theme.of(context).textTheme.titleMedium?.copyWith(
            color: Theme.of(context).hintColor,
          ),
        ),
        const SizedBox(height: 8),
        Text(
          'Load a molecule to view its 3D structure',
          style: TextStyle(
            color: Theme.of(context).hintColor,
          ),
        ),
        const SizedBox(height: 16),
        Container(
          padding: const EdgeInsets.all(12),
          decoration: BoxDecoration(
            color: Colors.blue.withOpacity(0.1),
            borderRadius: BorderRadius.circular(8),
          ),
          child: Column(
            children: [
              Text(
                'Controls:',
                style: Theme.of(context).textTheme.titleSmall?.copyWith(
                  fontWeight: FontWeight.bold,
                ),
              ),
              const SizedBox(height: 8),
              Text(
                '• Drag to rotate molecule\n'
                '• Scroll to zoom in/out\n'
                '• Use controls to change view mode',
                style: Theme.of(context).textTheme.bodySmall,
                textAlign: TextAlign.center,
              ),
            ],
          ),
        ),
      ],
    );
  }

  Widget _buildControls() {
    return Container(
      padding: const EdgeInsets.all(8),
      decoration: BoxDecoration(
        color: Theme.of(context).cardColor,
        border: Border.all(color: Theme.of(context).dividerColor),
        borderRadius: BorderRadius.circular(8),
      ),
      child: Column(
        mainAxisSize: MainAxisSize.min,
        children: [
          // Zoom controls
          Row(
            mainAxisAlignment: MainAxisAlignment.spaceEvenly,
            children: [
              IconButton(
                onPressed: () {
                  setState(() {
                    _zoom = (_zoom * 0.8).clamp(0.1, 5.0);
                  });
                },
                icon: const Icon(Icons.zoom_out),
                tooltip: 'Zoom Out',
              ),
              Text(
                '${(_zoom * 100).round()}%',
                style: Theme.of(context).textTheme.bodySmall,
              ),
              IconButton(
                onPressed: () {
                  setState(() {
                    _zoom = (_zoom * 1.25).clamp(0.1, 5.0);
                  });
                },
                icon: const Icon(Icons.zoom_in),
                tooltip: 'Zoom In',
              ),
            ],
          ),
          
          const Divider(height: 1),
          
          // Render mode
          Row(
            mainAxisAlignment: MainAxisAlignment.spaceEvenly,
            children: [
              IconButton(
                onPressed: () {
                  setState(() {
                    _renderMode = 'ball_stick';
                  });
                },
                icon: Icon(
                  Icons.scatter_plot,
                  color: _renderMode == 'ball_stick'
                      ? Theme.of(context).primaryColor
                      : null,
                ),
                tooltip: 'Ball & Stick',
              ),
              IconButton(
                onPressed: () {
                  setState(() {
                    _renderMode = 'space_filling';
                  });
                },
                icon: Icon(
                  Icons.circle,
                  color: _renderMode == 'space_filling'
                      ? Theme.of(context).primaryColor
                      : null,
                ),
                tooltip: 'Space Filling',
              ),
              IconButton(
                onPressed: () {
                  setState(() {
                    _renderMode = 'wireframe';
                  });
                },
                icon: Icon(
                  Icons.account_tree,
                  color: _renderMode == 'wireframe'
                      ? Theme.of(context).primaryColor
                      : null,
                ),
                tooltip: 'Wireframe',
              ),
            ],
          ),
          
          const Divider(height: 1),
          
          // Additional controls
          Row(
            mainAxisAlignment: MainAxisAlignment.spaceEvenly,
            children: [
              IconButton(
                onPressed: () {
                  setState(() {
                    _showBonds = !_showBonds;
                  });
                },
                icon: Icon(
                  _showBonds ? Icons.link : Icons.link_off,
                  color: _showBonds ? Theme.of(context).primaryColor : null,
                ),
                tooltip: 'Toggle Bonds',
              ),
              IconButton(
                onPressed: () {
                  setState(() {
                    _rotationX = 0.0;
                    _rotationY = 0.0;
                    _zoom = 1.0;
                  });
                },
                icon: const Icon(Icons.refresh),
                tooltip: 'Reset View',
              ),
              IconButton(
                onPressed: () {
                  if (_rotationController.isAnimating) {
                    _rotationController.stop();
                  } else {
                    _rotationController.repeat();
                  }
                },
                icon: Icon(
                  _rotationController.isAnimating
                      ? Icons.pause
                      : Icons.play_arrow,
                ),
                tooltip: 'Auto Rotate',
              ),
            ],
          ),
        ],
      ),
    );
  }

  @override
  Widget build(BuildContext context) {
    return Column(
      children: [
        // Header
        Row(
          children: [
            Text(
              '3D Molecule Viewer',
              style: Theme.of(context).textTheme.titleMedium,
            ),
            const Spacer(),
            if (widget.molecule != null)
              Text(
                widget.molecule!.name,
                style: Theme.of(context).textTheme.bodyMedium?.copyWith(
                  fontWeight: FontWeight.w500,
                ),
              ),
          ],
        ),
        
        const SizedBox(height: 8),
        
        // 3D viewer
        Expanded(
          child: Container(
            decoration: BoxDecoration(
              color: Colors.black.withOpacity(0.05),
              borderRadius: BorderRadius.circular(8),
              border: Border.all(color: Theme.of(context).dividerColor),
            ),
            child: ClipRRect(
              borderRadius: BorderRadius.circular(8),
              child: _buildMoleculeViewer(),
            ),
          ),
        ),
        
        const SizedBox(height: 8),
        
        // Controls
        _buildControls(),
      ],
    );
  }
}

class MoleculePainter extends CustomPainter {
  final MoleculeData molecule;
  final double rotationX;
  final double rotationY;
  final double zoom;
  final bool showBonds;
  final String renderMode;

  MoleculePainter({
    required this.molecule,
    required this.rotationX,
    required this.rotationY,
    required this.zoom,
    required this.showBonds,
    required this.renderMode,
  });

  @override
  void paint(Canvas canvas, Size size) {
    if (molecule.atoms.isEmpty) return;

    final paint = Paint()..isAntiAlias = true;

    // Project all atoms to 2D
    final projectedAtoms = molecule.atoms.map((atom) {
      return _project3DTo2D(atom, size);
    }).toList();

    // Draw bonds first (behind atoms)
    if (showBonds && renderMode != 'space_filling') {
      _drawBonds(canvas, projectedAtoms, paint);
    }

    // Draw atoms
    _drawAtoms(canvas, projectedAtoms, paint);

    // Draw atom labels if wireframe mode
    if (renderMode == 'wireframe') {
      _drawAtomLabels(canvas, projectedAtoms, size);
    }
  }

  Offset _project3DTo2D(AtomData atom, Size size) {
    // Apply rotations
    final cosX = math.cos(rotationX);
    final sinX = math.sin(rotationX);
    final cosY = math.cos(rotationY);
    final sinY = math.sin(rotationY);
    
    // Rotate around Y axis
    final x1 = atom.x * cosY - atom.z * sinY;
    final z1 = atom.x * sinY + atom.z * cosY;
    
    // Rotate around X axis
    final y2 = atom.y * cosX - z1 * sinX;
    final z2 = atom.y * sinX + z1 * cosX;
    
    // Apply zoom and center
    final scale = zoom * math.min(size.width, size.height) * 0.3;
    final centerX = size.width / 2;
    final centerY = size.height / 2;
    
    return Offset(
      centerX + x1 * scale,
      centerY - y2 * scale, // Flip Y axis
    );
  }

  void _drawBonds(Canvas canvas, List<Offset> projectedAtoms, Paint paint) {
    final bonds = _detectBonds(molecule.atoms);
    
    for (final bond in bonds) {
      final start = projectedAtoms[bond[0]];
      final end = projectedAtoms[bond[1]];
      
      paint.color = Colors.grey.shade600;
      paint.strokeWidth = renderMode == 'wireframe' ? 1.0 : 3.0;
      paint.style = PaintingStyle.stroke;
      
      canvas.drawLine(start, end, paint);
    }
  }

  void _drawAtoms(Canvas canvas, List<Offset> projectedAtoms, Paint paint) {
    paint.style = PaintingStyle.fill;
    
    for (int i = 0; i < molecule.atoms.length; i++) {
      final atom = molecule.atoms[i];
      final position = projectedAtoms[i];
      
      double radius;
      switch (renderMode) {
        case 'space_filling':
          radius = _getAtomRadius(atom.symbol) * zoom * 20;
          break;
        case 'wireframe':
          radius = 2.0;
          break;
        default: // ball_stick
          radius = math.max(8.0, _getAtomRadius(atom.symbol) * zoom * 8);
      }
      
      paint.color = _getAtomColor(atom.symbol);
      
      // Add shadow for depth
      if (renderMode != 'wireframe') {
        final shadowPaint = Paint()
          ..color = Colors.black.withOpacity(0.3)
          ..style = PaintingStyle.fill;
        canvas.drawCircle(
          position + const Offset(2, 2),
          radius,
          shadowPaint,
        );
      }
      
      canvas.drawCircle(position, radius, paint);
      
      // Add highlight for 3D effect
      if (renderMode != 'wireframe') {
        final highlightPaint = Paint()
          ..color = Colors.white.withOpacity(0.6)
          ..style = PaintingStyle.fill;
        canvas.drawCircle(
          position + Offset(-radius * 0.3, -radius * 0.3),
          radius * 0.3,
          highlightPaint,
        );
      }
    }
  }

  void _drawAtomLabels(Canvas canvas, List<Offset> projectedAtoms, Size size) {
    for (int i = 0; i < molecule.atoms.length; i++) {
      final atom = molecule.atoms[i];
      final position = projectedAtoms[i];
      
      final textPainter = TextPainter(
        text: TextSpan(
          text: atom.symbol,
          style: const TextStyle(
            color: Colors.black,
            fontSize: 12,
            fontWeight: FontWeight.bold,
          ),
        ),
        textDirection: TextDirection.ltr,
      );
      
      textPainter.layout();
      textPainter.paint(
        canvas,
        position - Offset(textPainter.width / 2, textPainter.height / 2),
      );
    }
  }

  List<List<int>> _detectBonds(List<AtomData> atoms) {
    final bonds = <List<int>>[];
    
    for (int i = 0; i < atoms.length; i++) {
      for (int j = i + 1; j < atoms.length; j++) {
        final distance = _calculateDistance(atoms[i], atoms[j]);
        final atom1Radius = _getAtomRadius(atoms[i].symbol);
        final atom2Radius = _getAtomRadius(atoms[j].symbol);
        
        // Bond if distance is less than sum of covalent radii * 1.3
        final bondThreshold = (atom1Radius + atom2Radius) * 0.6;
        
        if (distance < bondThreshold) {
          bonds.add([i, j]);
        }
      }
    }
    
    return bonds;
  }

  double _calculateDistance(AtomData atom1, AtomData atom2) {
    final dx = atom1.x - atom2.x;
    final dy = atom1.y - atom2.y;
    final dz = atom1.z - atom2.z;
    return math.sqrt(dx * dx + dy * dy + dz * dz);
  }

  Color _getAtomColor(String symbol) {
    const atomColors = {
      'H': Color(0xFFFFFFFF), // White
      'C': Color(0xFF909090), // Gray
      'N': Color(0xFF3050F8), // Blue
      'O': Color(0xFFFF0D0D), // Red
      'F': Color(0xFF90E050), // Green
      'P': Color(0xFFFF8000), // Orange
      'S': Color(0xFFFFFF30), // Yellow
      'Cl': Color(0xFF1FF01F), // Green
    };
    return atomColors[symbol] ?? const Color(0xFFFF1493);
  }

  double _getAtomRadius(String symbol) {
    const atomRadii = {
      'H': 1.20,
      'C': 1.70,
      'N': 1.55,
      'O': 1.52,
      'F': 1.47,
      'P': 1.80,
      'S': 1.80,
      'Cl': 1.75,
    };
    return atomRadii[symbol] ?? 1.5;
  }

  @override
  bool shouldRepaint(covariant CustomPainter oldDelegate) {
    return true;
  }
}