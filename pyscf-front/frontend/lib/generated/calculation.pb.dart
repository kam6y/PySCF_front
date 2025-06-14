// Generated code stub for calculation.proto
// This is a minimal stub to allow the application to build
// Full implementation will be generated from proto files

import 'package:protobuf/protobuf.dart';

// Basic enum for calculation status
enum CalculationStatus {
  CALCULATION_STATUS_UNKNOWN,
  CALCULATION_STATUS_PENDING,
  CALCULATION_STATUS_RUNNING,
  CALCULATION_STATUS_COMPLETED,
  CALCULATION_STATUS_FAILED,
  CALCULATION_STATUS_CANCELLED,
}

// Basic Atom class
class Atom {
  String symbol = '';
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  
  Atom();
  
  Atom.fromData({
    required this.symbol,
    required this.x,
    required this.y,
    required this.z,
  });
}

// Basic Molecule class
class Molecule {
  String id = '';
  String name = '';
  String formula = '';
  double molecularWeight = 0.0;
  int charge = 0;
  int multiplicity = 1;
  List<Atom> atoms = [];
  
  Molecule();
}

// Basic CalculationProgress class
class CalculationProgress {
  String jobId = '';
  double progressPercentage = 0.0;
  CalculationStatus status = CalculationStatus.CALCULATION_STATUS_UNKNOWN;
  String message = '';
  String currentStep = '';
  
  CalculationProgress();
}

// Response classes
class MoleculeResponse {
  bool success = false;
  String message = '';
  Molecule? molecule;
  
  MoleculeResponse();
  
  bool hasMolecule() => molecule != null;
}

class InstanceResponse {
  bool success = false;
  String message = '';
  Instance? instance;
  
  InstanceResponse();
  
  bool hasInstance() => instance != null;
}

class Instance {
  String id = '';
  String name = '';
  String description = '';
  
  Instance();
}