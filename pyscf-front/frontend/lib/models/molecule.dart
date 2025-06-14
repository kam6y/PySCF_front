class AtomData {
  final String symbol;
  final double x;
  final double y;
  final double z;

  const AtomData({
    required this.symbol,
    required this.x,
    required this.y,
    required this.z,
  });

  factory AtomData.fromJson(Map<String, dynamic> json) {
    return AtomData(
      symbol: json['symbol'] as String,
      x: (json['x'] as num).toDouble(),
      y: (json['y'] as num).toDouble(),
      z: (json['z'] as num).toDouble(),
    );
  }

  Map<String, dynamic> toJson() {
    return {
      'symbol': symbol,
      'x': x,
      'y': y,
      'z': z,
    };
  }

  AtomData copyWith({
    String? symbol,
    double? x,
    double? y,
    double? z,
  }) {
    return AtomData(
      symbol: symbol ?? this.symbol,
      x: x ?? this.x,
      y: y ?? this.y,
      z: z ?? this.z,
    );
  }

  @override
  String toString() {
    return 'AtomData(symbol: $symbol, x: $x, y: $y, z: $z)';
  }

  @override
  bool operator ==(Object other) {
    if (identical(this, other)) return true;
    return other is AtomData &&
        other.symbol == symbol &&
        other.x == x &&
        other.y == y &&
        other.z == z;
  }

  @override
  int get hashCode {
    return Object.hash(symbol, x, y, z);
  }
}

class MoleculeData {
  final String id;
  final String name;
  final String formula;
  final List<AtomData> atoms;
  final int charge;
  final int multiplicity;
  final double molecularWeight;

  const MoleculeData({
    required this.id,
    required this.name,
    required this.formula,
    required this.atoms,
    this.charge = 0,
    this.multiplicity = 1,
    this.molecularWeight = 0.0,
  });

  factory MoleculeData.fromJson(Map<String, dynamic> json) {
    return MoleculeData(
      id: json['id'] as String,
      name: json['name'] as String,
      formula: json['formula'] as String,
      atoms: (json['atoms'] as List)
          .map((atomJson) => AtomData.fromJson(atomJson as Map<String, dynamic>))
          .toList(),
      charge: json['charge'] as int? ?? 0,
      multiplicity: json['multiplicity'] as int? ?? 1,
      molecularWeight: (json['molecular_weight'] as num?)?.toDouble() ?? 0.0,
    );
  }

  Map<String, dynamic> toJson() {
    return {
      'id': id,
      'name': name,
      'formula': formula,
      'atoms': atoms.map((atom) => atom.toJson()).toList(),
      'charge': charge,
      'multiplicity': multiplicity,
      'molecular_weight': molecularWeight,
    };
  }

  factory MoleculeData.fromXyz(String xyzContent, {String? name}) {
    final lines = xyzContent.trim().split('\n');
    if (lines.length < 2) {
      throw ArgumentError('Invalid XYZ format');
    }

    final atomCount = int.parse(lines[0]);
    final comment = lines[1];
    
    final atoms = <AtomData>[];
    for (int i = 2; i < 2 + atomCount; i++) {
      if (i >= lines.length) break;
      
      final parts = lines[i].trim().split(RegExp(r'\s+'));
      if (parts.length >= 4) {
        atoms.add(AtomData(
          symbol: parts[0],
          x: double.parse(parts[1]),
          y: double.parse(parts[2]),
          z: double.parse(parts[3]),
        ));
      }
    }

    // Simple formula generation
    final symbolCounts = <String, int>{};
    for (final atom in atoms) {
      symbolCounts[atom.symbol] = (symbolCounts[atom.symbol] ?? 0) + 1;
    }
    
    final formula = symbolCounts.entries
        .map((entry) => entry.value == 1 ? entry.key : '${entry.key}${entry.value}')
        .join('');

    return MoleculeData(
      id: DateTime.now().millisecondsSinceEpoch.toString(),
      name: name ?? 'Molecule',
      formula: formula,
      atoms: atoms,
    );
  }

  String toXyz() {
    final buffer = StringBuffer();
    buffer.writeln(atoms.length);
    buffer.writeln(name);
    
    for (final atom in atoms) {
      buffer.writeln('${atom.symbol.padRight(4)} ${atom.x.toStringAsFixed(6)} ${atom.y.toStringAsFixed(6)} ${atom.z.toStringAsFixed(6)}');
    }
    
    return buffer.toString();
  }

  MoleculeData copyWith({
    String? id,
    String? name,
    String? formula,
    List<AtomData>? atoms,
    int? charge,
    int? multiplicity,
    double? molecularWeight,
  }) {
    return MoleculeData(
      id: id ?? this.id,
      name: name ?? this.name,
      formula: formula ?? this.formula,
      atoms: atoms ?? this.atoms,
      charge: charge ?? this.charge,
      multiplicity: multiplicity ?? this.multiplicity,
      molecularWeight: molecularWeight ?? this.molecularWeight,
    );
  }

  @override
  String toString() {
    return 'MoleculeData(id: $id, name: $name, formula: $formula, atoms: ${atoms.length})';
  }

  @override
  bool operator ==(Object other) {
    if (identical(this, other)) return true;
    return other is MoleculeData &&
        other.id == id &&
        other.name == name &&
        other.formula == formula &&
        _listEquals(other.atoms, atoms) &&
        other.charge == charge &&
        other.multiplicity == multiplicity &&
        other.molecularWeight == molecularWeight;
  }

  @override
  int get hashCode {
    return Object.hash(
      id,
      name,
      formula,
      Object.hashAll(atoms),
      charge,
      multiplicity,
      molecularWeight,
    );
  }

  bool _listEquals<T>(List<T> list1, List<T> list2) {
    if (list1.length != list2.length) return false;
    for (int i = 0; i < list1.length; i++) {
      if (list1[i] != list2[i]) return false;
    }
    return true;
  }
}