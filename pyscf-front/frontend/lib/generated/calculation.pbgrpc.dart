//
//  Generated code. Do not modify.
//  source: calculation.proto
//
// @dart = 3.3

// ignore_for_file: annotate_overrides, camel_case_types, comment_references
// ignore_for_file: constant_identifier_names
// ignore_for_file: curly_braces_in_flow_control_structures
// ignore_for_file: deprecated_member_use_from_same_package, library_prefixes
// ignore_for_file: non_constant_identifier_names

import 'dart:async' as $async;
import 'dart:core' as $core;

import 'package:grpc/service_api.dart' as $grpc;
import 'package:protobuf/protobuf.dart' as $pb;

import 'calculation.pb.dart' as $0;

export 'calculation.pb.dart';

/// CalculationService は量子化学計算を管理するメインサービス
@$pb.GrpcServiceName('pyscf_front.CalculationService')
class CalculationServiceClient extends $grpc.Client {
  /// The hostname for this service.
  static const $core.String defaultHost = '';

  /// OAuth scopes needed for the client.
  static const $core.List<$core.String> oauthScopes = [
    '',
  ];

  static final _$createMolecule = $grpc.ClientMethod<$0.CreateMoleculeRequest, $0.MoleculeResponse>(
      '/pyscf_front.CalculationService/CreateMolecule',
      ($0.CreateMoleculeRequest value) => value.writeToBuffer(),
      ($core.List<$core.int> value) => $0.MoleculeResponse.fromBuffer(value));
  static final _$getMolecule = $grpc.ClientMethod<$0.GetMoleculeRequest, $0.MoleculeResponse>(
      '/pyscf_front.CalculationService/GetMolecule',
      ($0.GetMoleculeRequest value) => value.writeToBuffer(),
      ($core.List<$core.int> value) => $0.MoleculeResponse.fromBuffer(value));
  static final _$updateMolecule = $grpc.ClientMethod<$0.UpdateMoleculeRequest, $0.MoleculeResponse>(
      '/pyscf_front.CalculationService/UpdateMolecule',
      ($0.UpdateMoleculeRequest value) => value.writeToBuffer(),
      ($core.List<$core.int> value) => $0.MoleculeResponse.fromBuffer(value));
  static final _$deleteMolecule = $grpc.ClientMethod<$0.DeleteMoleculeRequest, $0.DeleteMoleculeResponse>(
      '/pyscf_front.CalculationService/DeleteMolecule',
      ($0.DeleteMoleculeRequest value) => value.writeToBuffer(),
      ($core.List<$core.int> value) => $0.DeleteMoleculeResponse.fromBuffer(value));
  static final _$createInstance = $grpc.ClientMethod<$0.CreateInstanceRequest, $0.InstanceResponse>(
      '/pyscf_front.CalculationService/CreateInstance',
      ($0.CreateInstanceRequest value) => value.writeToBuffer(),
      ($core.List<$core.int> value) => $0.InstanceResponse.fromBuffer(value));
  static final _$getInstance = $grpc.ClientMethod<$0.GetInstanceRequest, $0.InstanceResponse>(
      '/pyscf_front.CalculationService/GetInstance',
      ($0.GetInstanceRequest value) => value.writeToBuffer(),
      ($core.List<$core.int> value) => $0.InstanceResponse.fromBuffer(value));
  static final _$listInstances = $grpc.ClientMethod<$0.ListInstancesRequest, $0.ListInstancesResponse>(
      '/pyscf_front.CalculationService/ListInstances',
      ($0.ListInstancesRequest value) => value.writeToBuffer(),
      ($core.List<$core.int> value) => $0.ListInstancesResponse.fromBuffer(value));
  static final _$updateInstance = $grpc.ClientMethod<$0.UpdateInstanceRequest, $0.InstanceResponse>(
      '/pyscf_front.CalculationService/UpdateInstance',
      ($0.UpdateInstanceRequest value) => value.writeToBuffer(),
      ($core.List<$core.int> value) => $0.InstanceResponse.fromBuffer(value));
  static final _$deleteInstance = $grpc.ClientMethod<$0.DeleteInstanceRequest, $0.DeleteInstanceResponse>(
      '/pyscf_front.CalculationService/DeleteInstance',
      ($0.DeleteInstanceRequest value) => value.writeToBuffer(),
      ($core.List<$core.int> value) => $0.DeleteInstanceResponse.fromBuffer(value));
  static final _$startCalculation = $grpc.ClientMethod<$0.StartCalculationRequest, $0.CalculationProgress>(
      '/pyscf_front.CalculationService/StartCalculation',
      ($0.StartCalculationRequest value) => value.writeToBuffer(),
      ($core.List<$core.int> value) => $0.CalculationProgress.fromBuffer(value));
  static final _$cancelCalculation = $grpc.ClientMethod<$0.CancelCalculationRequest, $0.CancelCalculationResponse>(
      '/pyscf_front.CalculationService/CancelCalculation',
      ($0.CancelCalculationRequest value) => value.writeToBuffer(),
      ($core.List<$core.int> value) => $0.CancelCalculationResponse.fromBuffer(value));
  static final _$getCalculationStatus = $grpc.ClientMethod<$0.GetCalculationStatusRequest, $0.CalculationStatusResponse>(
      '/pyscf_front.CalculationService/GetCalculationStatus',
      ($0.GetCalculationStatusRequest value) => value.writeToBuffer(),
      ($core.List<$core.int> value) => $0.CalculationStatusResponse.fromBuffer(value));
  static final _$getResults = $grpc.ClientMethod<$0.GetResultsRequest, $0.ResultsResponse>(
      '/pyscf_front.CalculationService/GetResults',
      ($0.GetResultsRequest value) => value.writeToBuffer(),
      ($core.List<$core.int> value) => $0.ResultsResponse.fromBuffer(value));
  static final _$exportResults = $grpc.ClientMethod<$0.ExportResultsRequest, $0.ExportResultsResponse>(
      '/pyscf_front.CalculationService/ExportResults',
      ($0.ExportResultsRequest value) => value.writeToBuffer(),
      ($core.List<$core.int> value) => $0.ExportResultsResponse.fromBuffer(value));
  static final _$getJobQueue = $grpc.ClientMethod<$0.GetJobQueueRequest, $0.JobQueueResponse>(
      '/pyscf_front.CalculationService/GetJobQueue',
      ($0.GetJobQueueRequest value) => value.writeToBuffer(),
      ($core.List<$core.int> value) => $0.JobQueueResponse.fromBuffer(value));
  static final _$updateJobPriority = $grpc.ClientMethod<$0.UpdateJobPriorityRequest, $0.UpdateJobPriorityResponse>(
      '/pyscf_front.CalculationService/UpdateJobPriority',
      ($0.UpdateJobPriorityRequest value) => value.writeToBuffer(),
      ($core.List<$core.int> value) => $0.UpdateJobPriorityResponse.fromBuffer(value));
  static final _$getSystemInfo = $grpc.ClientMethod<$0.GetSystemInfoRequest, $0.SystemInfoResponse>(
      '/pyscf_front.CalculationService/GetSystemInfo',
      ($0.GetSystemInfoRequest value) => value.writeToBuffer(),
      ($core.List<$core.int> value) => $0.SystemInfoResponse.fromBuffer(value));
  static final _$healthCheck = $grpc.ClientMethod<$0.HealthCheckRequest, $0.HealthCheckResponse>(
      '/pyscf_front.CalculationService/HealthCheck',
      ($0.HealthCheckRequest value) => value.writeToBuffer(),
      ($core.List<$core.int> value) => $0.HealthCheckResponse.fromBuffer(value));

  CalculationServiceClient(super.channel, {super.options, super.interceptors});

  /// 分子操作
  $grpc.ResponseFuture<$0.MoleculeResponse> createMolecule($0.CreateMoleculeRequest request, {$grpc.CallOptions? options}) {
    return $createUnaryCall(_$createMolecule, request, options: options);
  }

  $grpc.ResponseFuture<$0.MoleculeResponse> getMolecule($0.GetMoleculeRequest request, {$grpc.CallOptions? options}) {
    return $createUnaryCall(_$getMolecule, request, options: options);
  }

  $grpc.ResponseFuture<$0.MoleculeResponse> updateMolecule($0.UpdateMoleculeRequest request, {$grpc.CallOptions? options}) {
    return $createUnaryCall(_$updateMolecule, request, options: options);
  }

  $grpc.ResponseFuture<$0.DeleteMoleculeResponse> deleteMolecule($0.DeleteMoleculeRequest request, {$grpc.CallOptions? options}) {
    return $createUnaryCall(_$deleteMolecule, request, options: options);
  }

  /// インスタンス操作
  $grpc.ResponseFuture<$0.InstanceResponse> createInstance($0.CreateInstanceRequest request, {$grpc.CallOptions? options}) {
    return $createUnaryCall(_$createInstance, request, options: options);
  }

  $grpc.ResponseFuture<$0.InstanceResponse> getInstance($0.GetInstanceRequest request, {$grpc.CallOptions? options}) {
    return $createUnaryCall(_$getInstance, request, options: options);
  }

  $grpc.ResponseFuture<$0.ListInstancesResponse> listInstances($0.ListInstancesRequest request, {$grpc.CallOptions? options}) {
    return $createUnaryCall(_$listInstances, request, options: options);
  }

  $grpc.ResponseFuture<$0.InstanceResponse> updateInstance($0.UpdateInstanceRequest request, {$grpc.CallOptions? options}) {
    return $createUnaryCall(_$updateInstance, request, options: options);
  }

  $grpc.ResponseFuture<$0.DeleteInstanceResponse> deleteInstance($0.DeleteInstanceRequest request, {$grpc.CallOptions? options}) {
    return $createUnaryCall(_$deleteInstance, request, options: options);
  }

  /// 計算実行
  $grpc.ResponseStream<$0.CalculationProgress> startCalculation($0.StartCalculationRequest request, {$grpc.CallOptions? options}) {
    return $createStreamingCall(_$startCalculation, $async.Stream.fromIterable([request]), options: options);
  }

  $grpc.ResponseFuture<$0.CancelCalculationResponse> cancelCalculation($0.CancelCalculationRequest request, {$grpc.CallOptions? options}) {
    return $createUnaryCall(_$cancelCalculation, request, options: options);
  }

  $grpc.ResponseFuture<$0.CalculationStatusResponse> getCalculationStatus($0.GetCalculationStatusRequest request, {$grpc.CallOptions? options}) {
    return $createUnaryCall(_$getCalculationStatus, request, options: options);
  }

  /// 結果取得
  $grpc.ResponseFuture<$0.ResultsResponse> getResults($0.GetResultsRequest request, {$grpc.CallOptions? options}) {
    return $createUnaryCall(_$getResults, request, options: options);
  }

  $grpc.ResponseFuture<$0.ExportResultsResponse> exportResults($0.ExportResultsRequest request, {$grpc.CallOptions? options}) {
    return $createUnaryCall(_$exportResults, request, options: options);
  }

  /// ジョブ管理
  $grpc.ResponseFuture<$0.JobQueueResponse> getJobQueue($0.GetJobQueueRequest request, {$grpc.CallOptions? options}) {
    return $createUnaryCall(_$getJobQueue, request, options: options);
  }

  $grpc.ResponseFuture<$0.UpdateJobPriorityResponse> updateJobPriority($0.UpdateJobPriorityRequest request, {$grpc.CallOptions? options}) {
    return $createUnaryCall(_$updateJobPriority, request, options: options);
  }

  /// システム情報
  $grpc.ResponseFuture<$0.SystemInfoResponse> getSystemInfo($0.GetSystemInfoRequest request, {$grpc.CallOptions? options}) {
    return $createUnaryCall(_$getSystemInfo, request, options: options);
  }

  $grpc.ResponseFuture<$0.HealthCheckResponse> healthCheck($0.HealthCheckRequest request, {$grpc.CallOptions? options}) {
    return $createUnaryCall(_$healthCheck, request, options: options);
  }
}

@$pb.GrpcServiceName('pyscf_front.CalculationService')
abstract class CalculationServiceBase extends $grpc.Service {
  $core.String get $name => 'pyscf_front.CalculationService';

  CalculationServiceBase() {
    $addMethod($grpc.ServiceMethod<$0.CreateMoleculeRequest, $0.MoleculeResponse>(
        'CreateMolecule',
        createMolecule_Pre,
        false,
        false,
        ($core.List<$core.int> value) => $0.CreateMoleculeRequest.fromBuffer(value),
        ($0.MoleculeResponse value) => value.writeToBuffer()));
    $addMethod($grpc.ServiceMethod<$0.GetMoleculeRequest, $0.MoleculeResponse>(
        'GetMolecule',
        getMolecule_Pre,
        false,
        false,
        ($core.List<$core.int> value) => $0.GetMoleculeRequest.fromBuffer(value),
        ($0.MoleculeResponse value) => value.writeToBuffer()));
    $addMethod($grpc.ServiceMethod<$0.UpdateMoleculeRequest, $0.MoleculeResponse>(
        'UpdateMolecule',
        updateMolecule_Pre,
        false,
        false,
        ($core.List<$core.int> value) => $0.UpdateMoleculeRequest.fromBuffer(value),
        ($0.MoleculeResponse value) => value.writeToBuffer()));
    $addMethod($grpc.ServiceMethod<$0.DeleteMoleculeRequest, $0.DeleteMoleculeResponse>(
        'DeleteMolecule',
        deleteMolecule_Pre,
        false,
        false,
        ($core.List<$core.int> value) => $0.DeleteMoleculeRequest.fromBuffer(value),
        ($0.DeleteMoleculeResponse value) => value.writeToBuffer()));
    $addMethod($grpc.ServiceMethod<$0.CreateInstanceRequest, $0.InstanceResponse>(
        'CreateInstance',
        createInstance_Pre,
        false,
        false,
        ($core.List<$core.int> value) => $0.CreateInstanceRequest.fromBuffer(value),
        ($0.InstanceResponse value) => value.writeToBuffer()));
    $addMethod($grpc.ServiceMethod<$0.GetInstanceRequest, $0.InstanceResponse>(
        'GetInstance',
        getInstance_Pre,
        false,
        false,
        ($core.List<$core.int> value) => $0.GetInstanceRequest.fromBuffer(value),
        ($0.InstanceResponse value) => value.writeToBuffer()));
    $addMethod($grpc.ServiceMethod<$0.ListInstancesRequest, $0.ListInstancesResponse>(
        'ListInstances',
        listInstances_Pre,
        false,
        false,
        ($core.List<$core.int> value) => $0.ListInstancesRequest.fromBuffer(value),
        ($0.ListInstancesResponse value) => value.writeToBuffer()));
    $addMethod($grpc.ServiceMethod<$0.UpdateInstanceRequest, $0.InstanceResponse>(
        'UpdateInstance',
        updateInstance_Pre,
        false,
        false,
        ($core.List<$core.int> value) => $0.UpdateInstanceRequest.fromBuffer(value),
        ($0.InstanceResponse value) => value.writeToBuffer()));
    $addMethod($grpc.ServiceMethod<$0.DeleteInstanceRequest, $0.DeleteInstanceResponse>(
        'DeleteInstance',
        deleteInstance_Pre,
        false,
        false,
        ($core.List<$core.int> value) => $0.DeleteInstanceRequest.fromBuffer(value),
        ($0.DeleteInstanceResponse value) => value.writeToBuffer()));
    $addMethod($grpc.ServiceMethod<$0.StartCalculationRequest, $0.CalculationProgress>(
        'StartCalculation',
        startCalculation_Pre,
        false,
        true,
        ($core.List<$core.int> value) => $0.StartCalculationRequest.fromBuffer(value),
        ($0.CalculationProgress value) => value.writeToBuffer()));
    $addMethod($grpc.ServiceMethod<$0.CancelCalculationRequest, $0.CancelCalculationResponse>(
        'CancelCalculation',
        cancelCalculation_Pre,
        false,
        false,
        ($core.List<$core.int> value) => $0.CancelCalculationRequest.fromBuffer(value),
        ($0.CancelCalculationResponse value) => value.writeToBuffer()));
    $addMethod($grpc.ServiceMethod<$0.GetCalculationStatusRequest, $0.CalculationStatusResponse>(
        'GetCalculationStatus',
        getCalculationStatus_Pre,
        false,
        false,
        ($core.List<$core.int> value) => $0.GetCalculationStatusRequest.fromBuffer(value),
        ($0.CalculationStatusResponse value) => value.writeToBuffer()));
    $addMethod($grpc.ServiceMethod<$0.GetResultsRequest, $0.ResultsResponse>(
        'GetResults',
        getResults_Pre,
        false,
        false,
        ($core.List<$core.int> value) => $0.GetResultsRequest.fromBuffer(value),
        ($0.ResultsResponse value) => value.writeToBuffer()));
    $addMethod($grpc.ServiceMethod<$0.ExportResultsRequest, $0.ExportResultsResponse>(
        'ExportResults',
        exportResults_Pre,
        false,
        false,
        ($core.List<$core.int> value) => $0.ExportResultsRequest.fromBuffer(value),
        ($0.ExportResultsResponse value) => value.writeToBuffer()));
    $addMethod($grpc.ServiceMethod<$0.GetJobQueueRequest, $0.JobQueueResponse>(
        'GetJobQueue',
        getJobQueue_Pre,
        false,
        false,
        ($core.List<$core.int> value) => $0.GetJobQueueRequest.fromBuffer(value),
        ($0.JobQueueResponse value) => value.writeToBuffer()));
    $addMethod($grpc.ServiceMethod<$0.UpdateJobPriorityRequest, $0.UpdateJobPriorityResponse>(
        'UpdateJobPriority',
        updateJobPriority_Pre,
        false,
        false,
        ($core.List<$core.int> value) => $0.UpdateJobPriorityRequest.fromBuffer(value),
        ($0.UpdateJobPriorityResponse value) => value.writeToBuffer()));
    $addMethod($grpc.ServiceMethod<$0.GetSystemInfoRequest, $0.SystemInfoResponse>(
        'GetSystemInfo',
        getSystemInfo_Pre,
        false,
        false,
        ($core.List<$core.int> value) => $0.GetSystemInfoRequest.fromBuffer(value),
        ($0.SystemInfoResponse value) => value.writeToBuffer()));
    $addMethod($grpc.ServiceMethod<$0.HealthCheckRequest, $0.HealthCheckResponse>(
        'HealthCheck',
        healthCheck_Pre,
        false,
        false,
        ($core.List<$core.int> value) => $0.HealthCheckRequest.fromBuffer(value),
        ($0.HealthCheckResponse value) => value.writeToBuffer()));
  }

  $async.Future<$0.MoleculeResponse> createMolecule_Pre($grpc.ServiceCall $call, $async.Future<$0.CreateMoleculeRequest> $request) async {
    return createMolecule($call, await $request);
  }

  $async.Future<$0.MoleculeResponse> getMolecule_Pre($grpc.ServiceCall $call, $async.Future<$0.GetMoleculeRequest> $request) async {
    return getMolecule($call, await $request);
  }

  $async.Future<$0.MoleculeResponse> updateMolecule_Pre($grpc.ServiceCall $call, $async.Future<$0.UpdateMoleculeRequest> $request) async {
    return updateMolecule($call, await $request);
  }

  $async.Future<$0.DeleteMoleculeResponse> deleteMolecule_Pre($grpc.ServiceCall $call, $async.Future<$0.DeleteMoleculeRequest> $request) async {
    return deleteMolecule($call, await $request);
  }

  $async.Future<$0.InstanceResponse> createInstance_Pre($grpc.ServiceCall $call, $async.Future<$0.CreateInstanceRequest> $request) async {
    return createInstance($call, await $request);
  }

  $async.Future<$0.InstanceResponse> getInstance_Pre($grpc.ServiceCall $call, $async.Future<$0.GetInstanceRequest> $request) async {
    return getInstance($call, await $request);
  }

  $async.Future<$0.ListInstancesResponse> listInstances_Pre($grpc.ServiceCall $call, $async.Future<$0.ListInstancesRequest> $request) async {
    return listInstances($call, await $request);
  }

  $async.Future<$0.InstanceResponse> updateInstance_Pre($grpc.ServiceCall $call, $async.Future<$0.UpdateInstanceRequest> $request) async {
    return updateInstance($call, await $request);
  }

  $async.Future<$0.DeleteInstanceResponse> deleteInstance_Pre($grpc.ServiceCall $call, $async.Future<$0.DeleteInstanceRequest> $request) async {
    return deleteInstance($call, await $request);
  }

  $async.Stream<$0.CalculationProgress> startCalculation_Pre($grpc.ServiceCall $call, $async.Future<$0.StartCalculationRequest> $request) async* {
    yield* startCalculation($call, await $request);
  }

  $async.Future<$0.CancelCalculationResponse> cancelCalculation_Pre($grpc.ServiceCall $call, $async.Future<$0.CancelCalculationRequest> $request) async {
    return cancelCalculation($call, await $request);
  }

  $async.Future<$0.CalculationStatusResponse> getCalculationStatus_Pre($grpc.ServiceCall $call, $async.Future<$0.GetCalculationStatusRequest> $request) async {
    return getCalculationStatus($call, await $request);
  }

  $async.Future<$0.ResultsResponse> getResults_Pre($grpc.ServiceCall $call, $async.Future<$0.GetResultsRequest> $request) async {
    return getResults($call, await $request);
  }

  $async.Future<$0.ExportResultsResponse> exportResults_Pre($grpc.ServiceCall $call, $async.Future<$0.ExportResultsRequest> $request) async {
    return exportResults($call, await $request);
  }

  $async.Future<$0.JobQueueResponse> getJobQueue_Pre($grpc.ServiceCall $call, $async.Future<$0.GetJobQueueRequest> $request) async {
    return getJobQueue($call, await $request);
  }

  $async.Future<$0.UpdateJobPriorityResponse> updateJobPriority_Pre($grpc.ServiceCall $call, $async.Future<$0.UpdateJobPriorityRequest> $request) async {
    return updateJobPriority($call, await $request);
  }

  $async.Future<$0.SystemInfoResponse> getSystemInfo_Pre($grpc.ServiceCall $call, $async.Future<$0.GetSystemInfoRequest> $request) async {
    return getSystemInfo($call, await $request);
  }

  $async.Future<$0.HealthCheckResponse> healthCheck_Pre($grpc.ServiceCall $call, $async.Future<$0.HealthCheckRequest> $request) async {
    return healthCheck($call, await $request);
  }

  $async.Future<$0.MoleculeResponse> createMolecule($grpc.ServiceCall call, $0.CreateMoleculeRequest request);
  $async.Future<$0.MoleculeResponse> getMolecule($grpc.ServiceCall call, $0.GetMoleculeRequest request);
  $async.Future<$0.MoleculeResponse> updateMolecule($grpc.ServiceCall call, $0.UpdateMoleculeRequest request);
  $async.Future<$0.DeleteMoleculeResponse> deleteMolecule($grpc.ServiceCall call, $0.DeleteMoleculeRequest request);
  $async.Future<$0.InstanceResponse> createInstance($grpc.ServiceCall call, $0.CreateInstanceRequest request);
  $async.Future<$0.InstanceResponse> getInstance($grpc.ServiceCall call, $0.GetInstanceRequest request);
  $async.Future<$0.ListInstancesResponse> listInstances($grpc.ServiceCall call, $0.ListInstancesRequest request);
  $async.Future<$0.InstanceResponse> updateInstance($grpc.ServiceCall call, $0.UpdateInstanceRequest request);
  $async.Future<$0.DeleteInstanceResponse> deleteInstance($grpc.ServiceCall call, $0.DeleteInstanceRequest request);
  $async.Stream<$0.CalculationProgress> startCalculation($grpc.ServiceCall call, $0.StartCalculationRequest request);
  $async.Future<$0.CancelCalculationResponse> cancelCalculation($grpc.ServiceCall call, $0.CancelCalculationRequest request);
  $async.Future<$0.CalculationStatusResponse> getCalculationStatus($grpc.ServiceCall call, $0.GetCalculationStatusRequest request);
  $async.Future<$0.ResultsResponse> getResults($grpc.ServiceCall call, $0.GetResultsRequest request);
  $async.Future<$0.ExportResultsResponse> exportResults($grpc.ServiceCall call, $0.ExportResultsRequest request);
  $async.Future<$0.JobQueueResponse> getJobQueue($grpc.ServiceCall call, $0.GetJobQueueRequest request);
  $async.Future<$0.UpdateJobPriorityResponse> updateJobPriority($grpc.ServiceCall call, $0.UpdateJobPriorityRequest request);
  $async.Future<$0.SystemInfoResponse> getSystemInfo($grpc.ServiceCall call, $0.GetSystemInfoRequest request);
  $async.Future<$0.HealthCheckResponse> healthCheck($grpc.ServiceCall call, $0.HealthCheckRequest request);
}
