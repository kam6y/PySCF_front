class OrbitalData {
  final int index;
  final double energy;
  final double occupancy;
  final String type; // 'occupied', 'virtual', 'homo', 'lumo'

  const OrbitalData({
    required this.index,
    required this.energy,
    required this.occupancy,
    required this.type,
  });

  factory OrbitalData.fromJson(Map<String, dynamic> json) {
    return OrbitalData(
      index: json['index'] as int,
      energy: (json['energy'] as num).toDouble(),
      occupancy: (json['occupancy'] as num).toDouble(),
      type: json['type'] as String,
    );
  }

  Map<String, dynamic> toJson() {
    return {
      'index': index,
      'energy': energy,
      'occupancy': occupancy,
      'type': type,
    };
  }

  @override
  String toString() {
    return 'OrbitalData(index: $index, energy: $energy, occupancy: $occupancy, type: $type)';
  }
}

class EnergyData {
  final double totalEnergy; // in Hartree
  final double? homoEnergy;
  final double? lumoEnergy;
  final double? homoLumoGap;
  final double? nuclearRepulsion;
  final double? electronicEnergy;

  const EnergyData({
    required this.totalEnergy,
    this.homoEnergy,
    this.lumoEnergy,
    this.homoLumoGap,
    this.nuclearRepulsion,
    this.electronicEnergy,
  });

  factory EnergyData.fromJson(Map<String, dynamic> json) {
    return EnergyData(
      totalEnergy: (json['total_energy'] as num).toDouble(),
      homoEnergy: (json['homo_energy'] as num?)?.toDouble(),
      lumoEnergy: (json['lumo_energy'] as num?)?.toDouble(),
      homoLumoGap: (json['homo_lumo_gap'] as num?)?.toDouble(),
      nuclearRepulsion: (json['nuclear_repulsion'] as num?)?.toDouble(),
      electronicEnergy: (json['electronic_energy'] as num?)?.toDouble(),
    );
  }

  Map<String, dynamic> toJson() {
    return {
      'total_energy': totalEnergy,
      'homo_energy': homoEnergy,
      'lumo_energy': lumoEnergy,
      'homo_lumo_gap': homoLumoGap,
      'nuclear_repulsion': nuclearRepulsion,
      'electronic_energy': electronicEnergy,
    };
  }

  // Energy unit conversions
  double get totalEnergyEv => totalEnergy * 27.211386245988; // Hartree to eV
  double get totalEnergyKcalMol => totalEnergy * 627.5094738; // Hartree to kcal/mol
  double? get homoLumoGapEv => homoLumoGap != null ? homoLumoGap! * 27.211386245988 : null;

  @override
  String toString() {
    return 'EnergyData(totalEnergy: $totalEnergy Ha, homoLumoGap: $homoLumoGap Ha)';
  }
}

class CalculationResult {
  final bool success;
  final String message;
  final double energy;
  final int calculationTime; // in seconds
  final EnergyData? energyData;
  final List<OrbitalData>? orbitals;
  final Map<String, dynamic>? additionalData;
  final DateTime? completedAt;

  const CalculationResult({
    required this.success,
    required this.message,
    required this.energy,
    required this.calculationTime,
    this.energyData,
    this.orbitals,
    this.additionalData,
    this.completedAt,
  });

  factory CalculationResult.fromJson(Map<String, dynamic> json) {
    return CalculationResult(
      success: json['success'] as bool,
      message: json['message'] as String,
      energy: (json['energy'] as num).toDouble(),
      calculationTime: json['calculation_time'] as int,
      energyData: json['energy_data'] != null
          ? EnergyData.fromJson(json['energy_data'] as Map<String, dynamic>)
          : null,
      orbitals: json['orbitals'] != null
          ? (json['orbitals'] as List)
              .map((orbitalJson) => OrbitalData.fromJson(orbitalJson as Map<String, dynamic>))
              .toList()
          : null,
      additionalData: json['additional_data'] as Map<String, dynamic>?,
      completedAt: json['completed_at'] != null
          ? DateTime.parse(json['completed_at'] as String)
          : null,
    );
  }

  Map<String, dynamic> toJson() {
    return {
      'success': success,
      'message': message,
      'energy': energy,
      'calculation_time': calculationTime,
      'energy_data': energyData?.toJson(),
      'orbitals': orbitals?.map((orbital) => orbital.toJson()).toList(),
      'additional_data': additionalData,
      'completed_at': completedAt?.toIso8601String(),
    };
  }

  factory CalculationResult.success({
    required double energy,
    required int calculationTime,
    EnergyData? energyData,
    List<OrbitalData>? orbitals,
    Map<String, dynamic>? additionalData,
    String? message,
  }) {
    return CalculationResult(
      success: true,
      message: message ?? 'Calculation completed successfully',
      energy: energy,
      calculationTime: calculationTime,
      energyData: energyData,
      orbitals: orbitals,
      additionalData: additionalData,
      completedAt: DateTime.now(),
    );
  }

  factory CalculationResult.failure({
    required String message,
    int calculationTime = 0,
    Map<String, dynamic>? additionalData,
  }) {
    return CalculationResult(
      success: false,
      message: message,
      energy: 0.0,
      calculationTime: calculationTime,
      additionalData: additionalData,
      completedAt: DateTime.now(),
    );
  }

  // Helper methods for result analysis
  bool get hasEnergyData => energyData != null;
  bool get hasOrbitalData => orbitals != null && orbitals!.isNotEmpty;
  
  OrbitalData? get homoOrbital {
    if (!hasOrbitalData) return null;
    try {
      return orbitals!.firstWhere((orbital) => orbital.type == 'homo');
    } catch (e) {
      return null;
    }
  }

  OrbitalData? get lumoOrbital {
    if (!hasOrbitalData) return null;
    try {
      return orbitals!.firstWhere((orbital) => orbital.type == 'lumo');
    } catch (e) {
      return null;
    }
  }

  List<OrbitalData> get occupiedOrbitals {
    if (!hasOrbitalData) return [];
    return orbitals!.where((orbital) => orbital.occupancy > 0).toList();
  }

  List<OrbitalData> get virtualOrbitals {
    if (!hasOrbitalData) return [];
    return orbitals!.where((orbital) => orbital.occupancy == 0).toList();
  }

  String get formattedCalculationTime {
    if (calculationTime < 60) {
      return '${calculationTime}s';
    } else if (calculationTime < 3600) {
      final minutes = calculationTime ~/ 60;
      final seconds = calculationTime % 60;
      return '${minutes}m ${seconds}s';
    } else {
      final hours = calculationTime ~/ 3600;
      final minutes = (calculationTime % 3600) ~/ 60;
      final seconds = calculationTime % 60;
      return '${hours}h ${minutes}m ${seconds}s';
    }
  }

  CalculationResult copyWith({
    bool? success,
    String? message,
    double? energy,
    int? calculationTime,
    EnergyData? energyData,
    List<OrbitalData>? orbitals,
    Map<String, dynamic>? additionalData,
    DateTime? completedAt,
  }) {
    return CalculationResult(
      success: success ?? this.success,
      message: message ?? this.message,
      energy: energy ?? this.energy,
      calculationTime: calculationTime ?? this.calculationTime,
      energyData: energyData ?? this.energyData,
      orbitals: orbitals ?? this.orbitals,
      additionalData: additionalData ?? this.additionalData,
      completedAt: completedAt ?? this.completedAt,
    );
  }

  @override
  String toString() {
    return 'CalculationResult(success: $success, energy: $energy Ha, time: $formattedCalculationTime)';
  }

  @override
  bool operator ==(Object other) {
    if (identical(this, other)) return true;
    return other is CalculationResult &&
        other.success == success &&
        other.message == message &&
        other.energy == energy &&
        other.calculationTime == calculationTime;
  }

  @override
  int get hashCode {
    return Object.hash(success, message, energy, calculationTime);
  }
}