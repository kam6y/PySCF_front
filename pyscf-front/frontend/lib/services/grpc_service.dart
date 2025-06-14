import 'package:grpc/grpc.dart';
import '../generated/calculation.pbgrpc.dart';

class GrpcService {
  late ClientChannel _channel;
  late CalculationServiceClient _client;
  
  static const String defaultHost = 'localhost';
  static const int defaultPort = 50051;
  
  /// Initialize the gRPC connection
  Future<void> initialize({
    String host = defaultHost,
    int port = defaultPort,
  }) async {
    _channel = ClientChannel(
      host,
      port: port,
      options: const ChannelOptions(
        credentials: ChannelCredentials.insecure(),
      ),
    );
    
    _client = CalculationServiceClient(_channel);
  }
  
  /// Close the gRPC connection
  Future<void> close() async {
    await _channel.shutdown();
  }
  
  /// Get the calculation service client
  CalculationServiceClient get client => _client;
  
  /// Health check
  Future<HealthCheckResponse> healthCheck() async {
    final request = HealthCheckRequest();
    return await _client.healthCheck(request);
  }
  
  /// Get system information
  Future<SystemInfoResponse> getSystemInfo() async {
    final request = GetSystemInfoRequest();
    return await _client.getSystemInfo(request);
  }
  
  /// Create a new molecule
  Future<MoleculeResponse> createMolecule({
    required String name,
    required String formula,
    required List<Atom> atoms,
    int charge = 0,
    int multiplicity = 1,
    String symmetry = '',
    GeometryType geometryType = GeometryType.GEOMETRY_TYPE_XYZ,
  }) async {
    final request = CreateMoleculeRequest()
      ..name = name
      ..formula = formula
      ..atoms.addAll(atoms)
      ..charge = charge
      ..multiplicity = multiplicity
      ..symmetry = symmetry
      ..geometryType = geometryType;
    
    return await _client.createMolecule(request);
  }
  
  /// Get molecule by ID
  Future<MoleculeResponse> getMolecule(String moleculeId) async {
    final request = GetMoleculeRequest()..moleculeId = moleculeId;
    return await _client.getMolecule(request);
  }
  
  /// Create a new calculation instance
  Future<InstanceResponse> createInstance({
    required String name,
    required String description,
    required String moleculeId,
    String? projectId,
  }) async {
    final request = CreateInstanceRequest()
      ..name = name
      ..description = description
      ..moleculeId = moleculeId;
    
    if (projectId != null) {
      request.projectId = projectId;
    }
    
    return await _client.createInstance(request);
  }
  
  /// Get instance by ID
  Future<InstanceResponse> getInstance(String instanceId) async {
    final request = GetInstanceRequest()..instanceId = instanceId;
    return await _client.getInstance(request);
  }
  
  /// List instances with optional filtering
  Future<ListInstancesResponse> listInstances({
    String? projectId,
    InstanceStatus? status,
    int page = 1,
    int pageSize = 10,
  }) async {
    final request = ListInstancesRequest()
      ..page = page
      ..pageSize = pageSize;
    
    if (projectId != null) {
      request.projectId = projectId;
    }
    
    if (status != null) {
      request.status = status;
    }
    
    return await _client.listInstances(request);
  }
  
  /// Start a calculation and listen to progress updates
  Stream<CalculationProgress> startCalculation({
    required String instanceId,
    required String method,
    required String basisSet,
    Map<String, String>? parameters,
    ConvergenceCriteria? convergenceCriteria,
    int? maxIterations,
    int? priority,
  }) {
    final request = StartCalculationRequest()
      ..instanceId = instanceId
      ..method = method
      ..basisSet = basisSet;
    
    if (parameters != null) {
      request.parameters.addAll(parameters);
    }
    
    if (convergenceCriteria != null) {
      request.convergenceCriteria = convergenceCriteria;
    }
    
    if (maxIterations != null) {
      request.maxIterations = maxIterations;
    }
    
    if (priority != null) {
      request.priority = priority;
    }
    
    return _client.startCalculation(request);
  }
  
  /// Get calculation results
  Future<ResultsResponse> getResults({
    required String calculationId,
    List<String>? resultTypes,
  }) async {
    final request = GetResultsRequest()..calculationId = calculationId;
    
    if (resultTypes != null) {
      request.resultTypes.addAll(resultTypes);
    }
    
    return await _client.getResults(request);
  }
  
  /// Cancel a running calculation
  Future<CancelCalculationResponse> cancelCalculation(String calculationId) async {
    final request = CancelCalculationRequest()..calculationId = calculationId;
    return await _client.cancelCalculation(request);
  }
}

/// Singleton instance for global access
final grpcService = GrpcService();