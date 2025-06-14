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

import 'dart:convert' as $convert;
import 'dart:core' as $core;
import 'dart:typed_data' as $typed_data;

@$core.Deprecated('Use geometryTypeDescriptor instead')
const GeometryType$json = {
  '1': 'GeometryType',
  '2': [
    {'1': 'GEOMETRY_TYPE_UNSPECIFIED', '2': 0},
    {'1': 'GEOMETRY_TYPE_XYZ', '2': 1},
    {'1': 'GEOMETRY_TYPE_ZMATRIX', '2': 2},
  ],
};

/// Descriptor for `GeometryType`. Decode as a `google.protobuf.EnumDescriptorProto`.
final $typed_data.Uint8List geometryTypeDescriptor = $convert.base64Decode(
    'CgxHZW9tZXRyeVR5cGUSHQoZR0VPTUVUUllfVFlQRV9VTlNQRUNJRklFRBAAEhUKEUdFT01FVF'
    'JZX1RZUEVfWFlaEAESGQoVR0VPTUVUUllfVFlQRV9aTUFUUklYEAI=');

@$core.Deprecated('Use instanceStatusDescriptor instead')
const InstanceStatus$json = {
  '1': 'InstanceStatus',
  '2': [
    {'1': 'INSTANCE_STATUS_UNSPECIFIED', '2': 0},
    {'1': 'INSTANCE_STATUS_DRAFT', '2': 1},
    {'1': 'INSTANCE_STATUS_READY', '2': 2},
    {'1': 'INSTANCE_STATUS_RUNNING', '2': 3},
    {'1': 'INSTANCE_STATUS_COMPLETED', '2': 4},
    {'1': 'INSTANCE_STATUS_ERROR', '2': 5},
  ],
};

/// Descriptor for `InstanceStatus`. Decode as a `google.protobuf.EnumDescriptorProto`.
final $typed_data.Uint8List instanceStatusDescriptor = $convert.base64Decode(
    'Cg5JbnN0YW5jZVN0YXR1cxIfChtJTlNUQU5DRV9TVEFUVVNfVU5TUEVDSUZJRUQQABIZChVJTl'
    'NUQU5DRV9TVEFUVVNfRFJBRlQQARIZChVJTlNUQU5DRV9TVEFUVVNfUkVBRFkQAhIbChdJTlNU'
    'QU5DRV9TVEFUVVNfUlVOTklORxADEh0KGUlOU1RBTkNFX1NUQVRVU19DT01QTEVURUQQBBIZCh'
    'VJTlNUQU5DRV9TVEFUVVNfRVJST1IQBQ==');

@$core.Deprecated('Use calculationStatusDescriptor instead')
const CalculationStatus$json = {
  '1': 'CalculationStatus',
  '2': [
    {'1': 'CALCULATION_STATUS_UNSPECIFIED', '2': 0},
    {'1': 'CALCULATION_STATUS_PENDING', '2': 1},
    {'1': 'CALCULATION_STATUS_RUNNING', '2': 2},
    {'1': 'CALCULATION_STATUS_COMPLETED', '2': 3},
    {'1': 'CALCULATION_STATUS_FAILED', '2': 4},
    {'1': 'CALCULATION_STATUS_CANCELLED', '2': 5},
  ],
};

/// Descriptor for `CalculationStatus`. Decode as a `google.protobuf.EnumDescriptorProto`.
final $typed_data.Uint8List calculationStatusDescriptor = $convert.base64Decode(
    'ChFDYWxjdWxhdGlvblN0YXR1cxIiCh5DQUxDVUxBVElPTl9TVEFUVVNfVU5TUEVDSUZJRUQQAB'
    'IeChpDQUxDVUxBVElPTl9TVEFUVVNfUEVORElORxABEh4KGkNBTENVTEFUSU9OX1NUQVRVU19S'
    'VU5OSU5HEAISIAocQ0FMQ1VMQVRJT05fU1RBVFVTX0NPTVBMRVRFRBADEh0KGUNBTENVTEFUSU'
    '9OX1NUQVRVU19GQUlMRUQQBBIgChxDQUxDVUxBVElPTl9TVEFUVVNfQ0FOQ0VMTEVEEAU=');

@$core.Deprecated('Use exportFormatDescriptor instead')
const ExportFormat$json = {
  '1': 'ExportFormat',
  '2': [
    {'1': 'EXPORT_FORMAT_UNSPECIFIED', '2': 0},
    {'1': 'EXPORT_FORMAT_JSON', '2': 1},
    {'1': 'EXPORT_FORMAT_XML', '2': 2},
    {'1': 'EXPORT_FORMAT_HDF5', '2': 3},
    {'1': 'EXPORT_FORMAT_PDF', '2': 4},
  ],
};

/// Descriptor for `ExportFormat`. Decode as a `google.protobuf.EnumDescriptorProto`.
final $typed_data.Uint8List exportFormatDescriptor = $convert.base64Decode(
    'CgxFeHBvcnRGb3JtYXQSHQoZRVhQT1JUX0ZPUk1BVF9VTlNQRUNJRklFRBAAEhYKEkVYUE9SVF'
    '9GT1JNQVRfSlNPThABEhUKEUVYUE9SVF9GT1JNQVRfWE1MEAISFgoSRVhQT1JUX0ZPUk1BVF9I'
    'REY1EAMSFQoRRVhQT1JUX0ZPUk1BVF9QREYQBA==');

@$core.Deprecated('Use jobStatusDescriptor instead')
const JobStatus$json = {
  '1': 'JobStatus',
  '2': [
    {'1': 'JOB_STATUS_UNSPECIFIED', '2': 0},
    {'1': 'JOB_STATUS_WAITING', '2': 1},
    {'1': 'JOB_STATUS_RUNNING', '2': 2},
    {'1': 'JOB_STATUS_COMPLETED', '2': 3},
    {'1': 'JOB_STATUS_FAILED', '2': 4},
  ],
};

/// Descriptor for `JobStatus`. Decode as a `google.protobuf.EnumDescriptorProto`.
final $typed_data.Uint8List jobStatusDescriptor = $convert.base64Decode(
    'CglKb2JTdGF0dXMSGgoWSk9CX1NUQVRVU19VTlNQRUNJRklFRBAAEhYKEkpPQl9TVEFUVVNfV0'
    'FJVElORxABEhYKEkpPQl9TVEFUVVNfUlVOTklORxACEhgKFEpPQl9TVEFUVVNfQ09NUExFVEVE'
    'EAMSFQoRSk9CX1NUQVRVU19GQUlMRUQQBA==');

@$core.Deprecated('Use createMoleculeRequestDescriptor instead')
const CreateMoleculeRequest$json = {
  '1': 'CreateMoleculeRequest',
  '2': [
    {'1': 'name', '3': 1, '4': 1, '5': 9, '10': 'name'},
    {'1': 'formula', '3': 2, '4': 1, '5': 9, '10': 'formula'},
    {'1': 'atoms', '3': 3, '4': 3, '5': 11, '6': '.pyscf_front.Atom', '10': 'atoms'},
    {'1': 'charge', '3': 4, '4': 1, '5': 5, '10': 'charge'},
    {'1': 'multiplicity', '3': 5, '4': 1, '5': 5, '10': 'multiplicity'},
    {'1': 'symmetry', '3': 6, '4': 1, '5': 9, '10': 'symmetry'},
    {'1': 'geometry_type', '3': 7, '4': 1, '5': 14, '6': '.pyscf_front.GeometryType', '10': 'geometryType'},
  ],
};

/// Descriptor for `CreateMoleculeRequest`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List createMoleculeRequestDescriptor = $convert.base64Decode(
    'ChVDcmVhdGVNb2xlY3VsZVJlcXVlc3QSEgoEbmFtZRgBIAEoCVIEbmFtZRIYCgdmb3JtdWxhGA'
    'IgASgJUgdmb3JtdWxhEicKBWF0b21zGAMgAygLMhEucHlzY2ZfZnJvbnQuQXRvbVIFYXRvbXMS'
    'FgoGY2hhcmdlGAQgASgFUgZjaGFyZ2USIgoMbXVsdGlwbGljaXR5GAUgASgFUgxtdWx0aXBsaW'
    'NpdHkSGgoIc3ltbWV0cnkYBiABKAlSCHN5bW1ldHJ5Ej4KDWdlb21ldHJ5X3R5cGUYByABKA4y'
    'GS5weXNjZl9mcm9udC5HZW9tZXRyeVR5cGVSDGdlb21ldHJ5VHlwZQ==');

@$core.Deprecated('Use getMoleculeRequestDescriptor instead')
const GetMoleculeRequest$json = {
  '1': 'GetMoleculeRequest',
  '2': [
    {'1': 'molecule_id', '3': 1, '4': 1, '5': 9, '10': 'moleculeId'},
  ],
};

/// Descriptor for `GetMoleculeRequest`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List getMoleculeRequestDescriptor = $convert.base64Decode(
    'ChJHZXRNb2xlY3VsZVJlcXVlc3QSHwoLbW9sZWN1bGVfaWQYASABKAlSCm1vbGVjdWxlSWQ=');

@$core.Deprecated('Use updateMoleculeRequestDescriptor instead')
const UpdateMoleculeRequest$json = {
  '1': 'UpdateMoleculeRequest',
  '2': [
    {'1': 'molecule_id', '3': 1, '4': 1, '5': 9, '10': 'moleculeId'},
    {'1': 'name', '3': 2, '4': 1, '5': 9, '9': 0, '10': 'name', '17': true},
    {'1': 'formula', '3': 3, '4': 1, '5': 9, '9': 1, '10': 'formula', '17': true},
    {'1': 'atoms', '3': 4, '4': 3, '5': 11, '6': '.pyscf_front.Atom', '10': 'atoms'},
    {'1': 'charge', '3': 5, '4': 1, '5': 5, '9': 2, '10': 'charge', '17': true},
    {'1': 'multiplicity', '3': 6, '4': 1, '5': 5, '9': 3, '10': 'multiplicity', '17': true},
    {'1': 'symmetry', '3': 7, '4': 1, '5': 9, '9': 4, '10': 'symmetry', '17': true},
  ],
  '8': [
    {'1': '_name'},
    {'1': '_formula'},
    {'1': '_charge'},
    {'1': '_multiplicity'},
    {'1': '_symmetry'},
  ],
};

/// Descriptor for `UpdateMoleculeRequest`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List updateMoleculeRequestDescriptor = $convert.base64Decode(
    'ChVVcGRhdGVNb2xlY3VsZVJlcXVlc3QSHwoLbW9sZWN1bGVfaWQYASABKAlSCm1vbGVjdWxlSW'
    'QSFwoEbmFtZRgCIAEoCUgAUgRuYW1liAEBEh0KB2Zvcm11bGEYAyABKAlIAVIHZm9ybXVsYYgB'
    'ARInCgVhdG9tcxgEIAMoCzIRLnB5c2NmX2Zyb250LkF0b21SBWF0b21zEhsKBmNoYXJnZRgFIA'
    'EoBUgCUgZjaGFyZ2WIAQESJwoMbXVsdGlwbGljaXR5GAYgASgFSANSDG11bHRpcGxpY2l0eYgB'
    'ARIfCghzeW1tZXRyeRgHIAEoCUgEUghzeW1tZXRyeYgBAUIHCgVfbmFtZUIKCghfZm9ybXVsYU'
    'IJCgdfY2hhcmdlQg8KDV9tdWx0aXBsaWNpdHlCCwoJX3N5bW1ldHJ5');

@$core.Deprecated('Use deleteMoleculeRequestDescriptor instead')
const DeleteMoleculeRequest$json = {
  '1': 'DeleteMoleculeRequest',
  '2': [
    {'1': 'molecule_id', '3': 1, '4': 1, '5': 9, '10': 'moleculeId'},
  ],
};

/// Descriptor for `DeleteMoleculeRequest`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List deleteMoleculeRequestDescriptor = $convert.base64Decode(
    'ChVEZWxldGVNb2xlY3VsZVJlcXVlc3QSHwoLbW9sZWN1bGVfaWQYASABKAlSCm1vbGVjdWxlSW'
    'Q=');

@$core.Deprecated('Use deleteMoleculeResponseDescriptor instead')
const DeleteMoleculeResponse$json = {
  '1': 'DeleteMoleculeResponse',
  '2': [
    {'1': 'success', '3': 1, '4': 1, '5': 8, '10': 'success'},
    {'1': 'message', '3': 2, '4': 1, '5': 9, '10': 'message'},
  ],
};

/// Descriptor for `DeleteMoleculeResponse`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List deleteMoleculeResponseDescriptor = $convert.base64Decode(
    'ChZEZWxldGVNb2xlY3VsZVJlc3BvbnNlEhgKB3N1Y2Nlc3MYASABKAhSB3N1Y2Nlc3MSGAoHbW'
    'Vzc2FnZRgCIAEoCVIHbWVzc2FnZQ==');

@$core.Deprecated('Use moleculeResponseDescriptor instead')
const MoleculeResponse$json = {
  '1': 'MoleculeResponse',
  '2': [
    {'1': 'success', '3': 1, '4': 1, '5': 8, '10': 'success'},
    {'1': 'message', '3': 2, '4': 1, '5': 9, '10': 'message'},
    {'1': 'molecule', '3': 3, '4': 1, '5': 11, '6': '.pyscf_front.Molecule', '9': 0, '10': 'molecule', '17': true},
  ],
  '8': [
    {'1': '_molecule'},
  ],
};

/// Descriptor for `MoleculeResponse`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List moleculeResponseDescriptor = $convert.base64Decode(
    'ChBNb2xlY3VsZVJlc3BvbnNlEhgKB3N1Y2Nlc3MYASABKAhSB3N1Y2Nlc3MSGAoHbWVzc2FnZR'
    'gCIAEoCVIHbWVzc2FnZRI2Cghtb2xlY3VsZRgDIAEoCzIVLnB5c2NmX2Zyb250Lk1vbGVjdWxl'
    'SABSCG1vbGVjdWxliAEBQgsKCV9tb2xlY3VsZQ==');

@$core.Deprecated('Use moleculeDescriptor instead')
const Molecule$json = {
  '1': 'Molecule',
  '2': [
    {'1': 'id', '3': 1, '4': 1, '5': 9, '10': 'id'},
    {'1': 'name', '3': 2, '4': 1, '5': 9, '10': 'name'},
    {'1': 'formula', '3': 3, '4': 1, '5': 9, '10': 'formula'},
    {'1': 'molecular_weight', '3': 4, '4': 1, '5': 1, '10': 'molecularWeight'},
    {'1': 'atoms', '3': 5, '4': 3, '5': 11, '6': '.pyscf_front.Atom', '10': 'atoms'},
    {'1': 'charge', '3': 6, '4': 1, '5': 5, '10': 'charge'},
    {'1': 'multiplicity', '3': 7, '4': 1, '5': 5, '10': 'multiplicity'},
    {'1': 'symmetry', '3': 8, '4': 1, '5': 9, '10': 'symmetry'},
    {'1': 'geometry_type', '3': 9, '4': 1, '5': 14, '6': '.pyscf_front.GeometryType', '10': 'geometryType'},
    {'1': 'created_at', '3': 10, '4': 1, '5': 3, '10': 'createdAt'},
  ],
};

/// Descriptor for `Molecule`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List moleculeDescriptor = $convert.base64Decode(
    'CghNb2xlY3VsZRIOCgJpZBgBIAEoCVICaWQSEgoEbmFtZRgCIAEoCVIEbmFtZRIYCgdmb3JtdW'
    'xhGAMgASgJUgdmb3JtdWxhEikKEG1vbGVjdWxhcl93ZWlnaHQYBCABKAFSD21vbGVjdWxhcldl'
    'aWdodBInCgVhdG9tcxgFIAMoCzIRLnB5c2NmX2Zyb250LkF0b21SBWF0b21zEhYKBmNoYXJnZR'
    'gGIAEoBVIGY2hhcmdlEiIKDG11bHRpcGxpY2l0eRgHIAEoBVIMbXVsdGlwbGljaXR5EhoKCHN5'
    'bW1ldHJ5GAggASgJUghzeW1tZXRyeRI+Cg1nZW9tZXRyeV90eXBlGAkgASgOMhkucHlzY2ZfZn'
    'JvbnQuR2VvbWV0cnlUeXBlUgxnZW9tZXRyeVR5cGUSHQoKY3JlYXRlZF9hdBgKIAEoA1IJY3Jl'
    'YXRlZEF0');

@$core.Deprecated('Use atomDescriptor instead')
const Atom$json = {
  '1': 'Atom',
  '2': [
    {'1': 'symbol', '3': 1, '4': 1, '5': 9, '10': 'symbol'},
    {'1': 'x', '3': 2, '4': 1, '5': 1, '10': 'x'},
    {'1': 'y', '3': 3, '4': 1, '5': 1, '10': 'y'},
    {'1': 'z', '3': 4, '4': 1, '5': 1, '10': 'z'},
  ],
};

/// Descriptor for `Atom`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List atomDescriptor = $convert.base64Decode(
    'CgRBdG9tEhYKBnN5bWJvbBgBIAEoCVIGc3ltYm9sEgwKAXgYAiABKAFSAXgSDAoBeRgDIAEoAV'
    'IBeRIMCgF6GAQgASgBUgF6');

@$core.Deprecated('Use createInstanceRequestDescriptor instead')
const CreateInstanceRequest$json = {
  '1': 'CreateInstanceRequest',
  '2': [
    {'1': 'name', '3': 1, '4': 1, '5': 9, '10': 'name'},
    {'1': 'description', '3': 2, '4': 1, '5': 9, '10': 'description'},
    {'1': 'molecule_id', '3': 3, '4': 1, '5': 9, '10': 'moleculeId'},
    {'1': 'project_id', '3': 4, '4': 1, '5': 9, '9': 0, '10': 'projectId', '17': true},
  ],
  '8': [
    {'1': '_project_id'},
  ],
};

/// Descriptor for `CreateInstanceRequest`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List createInstanceRequestDescriptor = $convert.base64Decode(
    'ChVDcmVhdGVJbnN0YW5jZVJlcXVlc3QSEgoEbmFtZRgBIAEoCVIEbmFtZRIgCgtkZXNjcmlwdG'
    'lvbhgCIAEoCVILZGVzY3JpcHRpb24SHwoLbW9sZWN1bGVfaWQYAyABKAlSCm1vbGVjdWxlSWQS'
    'IgoKcHJvamVjdF9pZBgEIAEoCUgAUglwcm9qZWN0SWSIAQFCDQoLX3Byb2plY3RfaWQ=');

@$core.Deprecated('Use getInstanceRequestDescriptor instead')
const GetInstanceRequest$json = {
  '1': 'GetInstanceRequest',
  '2': [
    {'1': 'instance_id', '3': 1, '4': 1, '5': 9, '10': 'instanceId'},
  ],
};

/// Descriptor for `GetInstanceRequest`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List getInstanceRequestDescriptor = $convert.base64Decode(
    'ChJHZXRJbnN0YW5jZVJlcXVlc3QSHwoLaW5zdGFuY2VfaWQYASABKAlSCmluc3RhbmNlSWQ=');

@$core.Deprecated('Use listInstancesRequestDescriptor instead')
const ListInstancesRequest$json = {
  '1': 'ListInstancesRequest',
  '2': [
    {'1': 'project_id', '3': 1, '4': 1, '5': 9, '9': 0, '10': 'projectId', '17': true},
    {'1': 'status', '3': 2, '4': 1, '5': 14, '6': '.pyscf_front.InstanceStatus', '9': 1, '10': 'status', '17': true},
    {'1': 'page', '3': 3, '4': 1, '5': 5, '10': 'page'},
    {'1': 'page_size', '3': 4, '4': 1, '5': 5, '10': 'pageSize'},
  ],
  '8': [
    {'1': '_project_id'},
    {'1': '_status'},
  ],
};

/// Descriptor for `ListInstancesRequest`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List listInstancesRequestDescriptor = $convert.base64Decode(
    'ChRMaXN0SW5zdGFuY2VzUmVxdWVzdBIiCgpwcm9qZWN0X2lkGAEgASgJSABSCXByb2plY3RJZI'
    'gBARI4CgZzdGF0dXMYAiABKA4yGy5weXNjZl9mcm9udC5JbnN0YW5jZVN0YXR1c0gBUgZzdGF0'
    'dXOIAQESEgoEcGFnZRgDIAEoBVIEcGFnZRIbCglwYWdlX3NpemUYBCABKAVSCHBhZ2VTaXplQg'
    '0KC19wcm9qZWN0X2lkQgkKB19zdGF0dXM=');

@$core.Deprecated('Use listInstancesResponseDescriptor instead')
const ListInstancesResponse$json = {
  '1': 'ListInstancesResponse',
  '2': [
    {'1': 'success', '3': 1, '4': 1, '5': 8, '10': 'success'},
    {'1': 'message', '3': 2, '4': 1, '5': 9, '10': 'message'},
    {'1': 'instances', '3': 3, '4': 3, '5': 11, '6': '.pyscf_front.Instance', '10': 'instances'},
    {'1': 'total_count', '3': 4, '4': 1, '5': 5, '10': 'totalCount'},
    {'1': 'page', '3': 5, '4': 1, '5': 5, '10': 'page'},
    {'1': 'page_size', '3': 6, '4': 1, '5': 5, '10': 'pageSize'},
  ],
};

/// Descriptor for `ListInstancesResponse`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List listInstancesResponseDescriptor = $convert.base64Decode(
    'ChVMaXN0SW5zdGFuY2VzUmVzcG9uc2USGAoHc3VjY2VzcxgBIAEoCFIHc3VjY2VzcxIYCgdtZX'
    'NzYWdlGAIgASgJUgdtZXNzYWdlEjMKCWluc3RhbmNlcxgDIAMoCzIVLnB5c2NmX2Zyb250Lklu'
    'c3RhbmNlUglpbnN0YW5jZXMSHwoLdG90YWxfY291bnQYBCABKAVSCnRvdGFsQ291bnQSEgoEcG'
    'FnZRgFIAEoBVIEcGFnZRIbCglwYWdlX3NpemUYBiABKAVSCHBhZ2VTaXpl');

@$core.Deprecated('Use updateInstanceRequestDescriptor instead')
const UpdateInstanceRequest$json = {
  '1': 'UpdateInstanceRequest',
  '2': [
    {'1': 'instance_id', '3': 1, '4': 1, '5': 9, '10': 'instanceId'},
    {'1': 'name', '3': 2, '4': 1, '5': 9, '9': 0, '10': 'name', '17': true},
    {'1': 'description', '3': 3, '4': 1, '5': 9, '9': 1, '10': 'description', '17': true},
    {'1': 'status', '3': 4, '4': 1, '5': 14, '6': '.pyscf_front.InstanceStatus', '9': 2, '10': 'status', '17': true},
  ],
  '8': [
    {'1': '_name'},
    {'1': '_description'},
    {'1': '_status'},
  ],
};

/// Descriptor for `UpdateInstanceRequest`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List updateInstanceRequestDescriptor = $convert.base64Decode(
    'ChVVcGRhdGVJbnN0YW5jZVJlcXVlc3QSHwoLaW5zdGFuY2VfaWQYASABKAlSCmluc3RhbmNlSW'
    'QSFwoEbmFtZRgCIAEoCUgAUgRuYW1liAEBEiUKC2Rlc2NyaXB0aW9uGAMgASgJSAFSC2Rlc2Ny'
    'aXB0aW9uiAEBEjgKBnN0YXR1cxgEIAEoDjIbLnB5c2NmX2Zyb250Lkluc3RhbmNlU3RhdHVzSA'
    'JSBnN0YXR1c4gBAUIHCgVfbmFtZUIOCgxfZGVzY3JpcHRpb25CCQoHX3N0YXR1cw==');

@$core.Deprecated('Use deleteInstanceRequestDescriptor instead')
const DeleteInstanceRequest$json = {
  '1': 'DeleteInstanceRequest',
  '2': [
    {'1': 'instance_id', '3': 1, '4': 1, '5': 9, '10': 'instanceId'},
  ],
};

/// Descriptor for `DeleteInstanceRequest`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List deleteInstanceRequestDescriptor = $convert.base64Decode(
    'ChVEZWxldGVJbnN0YW5jZVJlcXVlc3QSHwoLaW5zdGFuY2VfaWQYASABKAlSCmluc3RhbmNlSW'
    'Q=');

@$core.Deprecated('Use deleteInstanceResponseDescriptor instead')
const DeleteInstanceResponse$json = {
  '1': 'DeleteInstanceResponse',
  '2': [
    {'1': 'success', '3': 1, '4': 1, '5': 8, '10': 'success'},
    {'1': 'message', '3': 2, '4': 1, '5': 9, '10': 'message'},
  ],
};

/// Descriptor for `DeleteInstanceResponse`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List deleteInstanceResponseDescriptor = $convert.base64Decode(
    'ChZEZWxldGVJbnN0YW5jZVJlc3BvbnNlEhgKB3N1Y2Nlc3MYASABKAhSB3N1Y2Nlc3MSGAoHbW'
    'Vzc2FnZRgCIAEoCVIHbWVzc2FnZQ==');

@$core.Deprecated('Use instanceResponseDescriptor instead')
const InstanceResponse$json = {
  '1': 'InstanceResponse',
  '2': [
    {'1': 'success', '3': 1, '4': 1, '5': 8, '10': 'success'},
    {'1': 'message', '3': 2, '4': 1, '5': 9, '10': 'message'},
    {'1': 'instance', '3': 3, '4': 1, '5': 11, '6': '.pyscf_front.Instance', '9': 0, '10': 'instance', '17': true},
  ],
  '8': [
    {'1': '_instance'},
  ],
};

/// Descriptor for `InstanceResponse`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List instanceResponseDescriptor = $convert.base64Decode(
    'ChBJbnN0YW5jZVJlc3BvbnNlEhgKB3N1Y2Nlc3MYASABKAhSB3N1Y2Nlc3MSGAoHbWVzc2FnZR'
    'gCIAEoCVIHbWVzc2FnZRI2CghpbnN0YW5jZRgDIAEoCzIVLnB5c2NmX2Zyb250Lkluc3RhbmNl'
    'SABSCGluc3RhbmNliAEBQgsKCV9pbnN0YW5jZQ==');

@$core.Deprecated('Use instanceDescriptor instead')
const Instance$json = {
  '1': 'Instance',
  '2': [
    {'1': 'id', '3': 1, '4': 1, '5': 9, '10': 'id'},
    {'1': 'name', '3': 2, '4': 1, '5': 9, '10': 'name'},
    {'1': 'description', '3': 3, '4': 1, '5': 9, '10': 'description'},
    {'1': 'status', '3': 4, '4': 1, '5': 14, '6': '.pyscf_front.InstanceStatus', '10': 'status'},
    {'1': 'user_id', '3': 5, '4': 1, '5': 9, '9': 0, '10': 'userId', '17': true},
    {'1': 'project_id', '3': 6, '4': 1, '5': 9, '9': 1, '10': 'projectId', '17': true},
    {'1': 'created_at', '3': 7, '4': 1, '5': 3, '10': 'createdAt'},
    {'1': 'updated_at', '3': 8, '4': 1, '5': 3, '10': 'updatedAt'},
    {'1': 'molecule', '3': 9, '4': 1, '5': 11, '6': '.pyscf_front.Molecule', '9': 2, '10': 'molecule', '17': true},
    {'1': 'calculations', '3': 10, '4': 3, '5': 11, '6': '.pyscf_front.Calculation', '10': 'calculations'},
  ],
  '8': [
    {'1': '_user_id'},
    {'1': '_project_id'},
    {'1': '_molecule'},
  ],
};

/// Descriptor for `Instance`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List instanceDescriptor = $convert.base64Decode(
    'CghJbnN0YW5jZRIOCgJpZBgBIAEoCVICaWQSEgoEbmFtZRgCIAEoCVIEbmFtZRIgCgtkZXNjcm'
    'lwdGlvbhgDIAEoCVILZGVzY3JpcHRpb24SMwoGc3RhdHVzGAQgASgOMhsucHlzY2ZfZnJvbnQu'
    'SW5zdGFuY2VTdGF0dXNSBnN0YXR1cxIcCgd1c2VyX2lkGAUgASgJSABSBnVzZXJJZIgBARIiCg'
    'pwcm9qZWN0X2lkGAYgASgJSAFSCXByb2plY3RJZIgBARIdCgpjcmVhdGVkX2F0GAcgASgDUglj'
    'cmVhdGVkQXQSHQoKdXBkYXRlZF9hdBgIIAEoA1IJdXBkYXRlZEF0EjYKCG1vbGVjdWxlGAkgAS'
    'gLMhUucHlzY2ZfZnJvbnQuTW9sZWN1bGVIAlIIbW9sZWN1bGWIAQESPAoMY2FsY3VsYXRpb25z'
    'GAogAygLMhgucHlzY2ZfZnJvbnQuQ2FsY3VsYXRpb25SDGNhbGN1bGF0aW9uc0IKCghfdXNlcl'
    '9pZEINCgtfcHJvamVjdF9pZEILCglfbW9sZWN1bGU=');

@$core.Deprecated('Use startCalculationRequestDescriptor instead')
const StartCalculationRequest$json = {
  '1': 'StartCalculationRequest',
  '2': [
    {'1': 'instance_id', '3': 1, '4': 1, '5': 9, '10': 'instanceId'},
    {'1': 'method', '3': 2, '4': 1, '5': 9, '10': 'method'},
    {'1': 'basis_set', '3': 3, '4': 1, '5': 9, '10': 'basisSet'},
    {'1': 'parameters', '3': 4, '4': 3, '5': 11, '6': '.pyscf_front.StartCalculationRequest.ParametersEntry', '10': 'parameters'},
    {'1': 'convergence_criteria', '3': 5, '4': 1, '5': 11, '6': '.pyscf_front.ConvergenceCriteria', '9': 0, '10': 'convergenceCriteria', '17': true},
    {'1': 'max_iterations', '3': 6, '4': 1, '5': 5, '9': 1, '10': 'maxIterations', '17': true},
    {'1': 'priority', '3': 7, '4': 1, '5': 5, '9': 2, '10': 'priority', '17': true},
  ],
  '3': [StartCalculationRequest_ParametersEntry$json],
  '8': [
    {'1': '_convergence_criteria'},
    {'1': '_max_iterations'},
    {'1': '_priority'},
  ],
};

@$core.Deprecated('Use startCalculationRequestDescriptor instead')
const StartCalculationRequest_ParametersEntry$json = {
  '1': 'ParametersEntry',
  '2': [
    {'1': 'key', '3': 1, '4': 1, '5': 9, '10': 'key'},
    {'1': 'value', '3': 2, '4': 1, '5': 9, '10': 'value'},
  ],
  '7': {'7': true},
};

/// Descriptor for `StartCalculationRequest`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List startCalculationRequestDescriptor = $convert.base64Decode(
    'ChdTdGFydENhbGN1bGF0aW9uUmVxdWVzdBIfCgtpbnN0YW5jZV9pZBgBIAEoCVIKaW5zdGFuY2'
    'VJZBIWCgZtZXRob2QYAiABKAlSBm1ldGhvZBIbCgliYXNpc19zZXQYAyABKAlSCGJhc2lzU2V0'
    'ElQKCnBhcmFtZXRlcnMYBCADKAsyNC5weXNjZl9mcm9udC5TdGFydENhbGN1bGF0aW9uUmVxdW'
    'VzdC5QYXJhbWV0ZXJzRW50cnlSCnBhcmFtZXRlcnMSWAoUY29udmVyZ2VuY2VfY3JpdGVyaWEY'
    'BSABKAsyIC5weXNjZl9mcm9udC5Db252ZXJnZW5jZUNyaXRlcmlhSABSE2NvbnZlcmdlbmNlQ3'
    'JpdGVyaWGIAQESKgoObWF4X2l0ZXJhdGlvbnMYBiABKAVIAVINbWF4SXRlcmF0aW9uc4gBARIf'
    'Cghwcmlvcml0eRgHIAEoBUgCUghwcmlvcml0eYgBARo9Cg9QYXJhbWV0ZXJzRW50cnkSEAoDa2'
    'V5GAEgASgJUgNrZXkSFAoFdmFsdWUYAiABKAlSBXZhbHVlOgI4AUIXChVfY29udmVyZ2VuY2Vf'
    'Y3JpdGVyaWFCEQoPX21heF9pdGVyYXRpb25zQgsKCV9wcmlvcml0eQ==');

@$core.Deprecated('Use cancelCalculationRequestDescriptor instead')
const CancelCalculationRequest$json = {
  '1': 'CancelCalculationRequest',
  '2': [
    {'1': 'calculation_id', '3': 1, '4': 1, '5': 9, '10': 'calculationId'},
  ],
};

/// Descriptor for `CancelCalculationRequest`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List cancelCalculationRequestDescriptor = $convert.base64Decode(
    'ChhDYW5jZWxDYWxjdWxhdGlvblJlcXVlc3QSJQoOY2FsY3VsYXRpb25faWQYASABKAlSDWNhbG'
    'N1bGF0aW9uSWQ=');

@$core.Deprecated('Use cancelCalculationResponseDescriptor instead')
const CancelCalculationResponse$json = {
  '1': 'CancelCalculationResponse',
  '2': [
    {'1': 'success', '3': 1, '4': 1, '5': 8, '10': 'success'},
    {'1': 'message', '3': 2, '4': 1, '5': 9, '10': 'message'},
  ],
};

/// Descriptor for `CancelCalculationResponse`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List cancelCalculationResponseDescriptor = $convert.base64Decode(
    'ChlDYW5jZWxDYWxjdWxhdGlvblJlc3BvbnNlEhgKB3N1Y2Nlc3MYASABKAhSB3N1Y2Nlc3MSGA'
    'oHbWVzc2FnZRgCIAEoCVIHbWVzc2FnZQ==');

@$core.Deprecated('Use getCalculationStatusRequestDescriptor instead')
const GetCalculationStatusRequest$json = {
  '1': 'GetCalculationStatusRequest',
  '2': [
    {'1': 'calculation_id', '3': 1, '4': 1, '5': 9, '10': 'calculationId'},
  ],
};

/// Descriptor for `GetCalculationStatusRequest`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List getCalculationStatusRequestDescriptor = $convert.base64Decode(
    'ChtHZXRDYWxjdWxhdGlvblN0YXR1c1JlcXVlc3QSJQoOY2FsY3VsYXRpb25faWQYASABKAlSDW'
    'NhbGN1bGF0aW9uSWQ=');

@$core.Deprecated('Use calculationStatusResponseDescriptor instead')
const CalculationStatusResponse$json = {
  '1': 'CalculationStatusResponse',
  '2': [
    {'1': 'success', '3': 1, '4': 1, '5': 8, '10': 'success'},
    {'1': 'message', '3': 2, '4': 1, '5': 9, '10': 'message'},
    {'1': 'status', '3': 3, '4': 1, '5': 14, '6': '.pyscf_front.CalculationStatus', '9': 0, '10': 'status', '17': true},
  ],
  '8': [
    {'1': '_status'},
  ],
};

/// Descriptor for `CalculationStatusResponse`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List calculationStatusResponseDescriptor = $convert.base64Decode(
    'ChlDYWxjdWxhdGlvblN0YXR1c1Jlc3BvbnNlEhgKB3N1Y2Nlc3MYASABKAhSB3N1Y2Nlc3MSGA'
    'oHbWVzc2FnZRgCIAEoCVIHbWVzc2FnZRI7CgZzdGF0dXMYAyABKA4yHi5weXNjZl9mcm9udC5D'
    'YWxjdWxhdGlvblN0YXR1c0gAUgZzdGF0dXOIAQFCCQoHX3N0YXR1cw==');

@$core.Deprecated('Use calculationProgressDescriptor instead')
const CalculationProgress$json = {
  '1': 'CalculationProgress',
  '2': [
    {'1': 'calculation_id', '3': 1, '4': 1, '5': 9, '10': 'calculationId'},
    {'1': 'job_id', '3': 2, '4': 1, '5': 9, '10': 'jobId'},
    {'1': 'status', '3': 3, '4': 1, '5': 14, '6': '.pyscf_front.CalculationStatus', '10': 'status'},
    {'1': 'progress_percentage', '3': 4, '4': 1, '5': 1, '10': 'progressPercentage'},
    {'1': 'current_step', '3': 5, '4': 1, '5': 9, '10': 'currentStep'},
    {'1': 'message', '3': 6, '4': 1, '5': 9, '10': 'message'},
    {'1': 'timestamp', '3': 7, '4': 1, '5': 3, '10': 'timestamp'},
    {'1': 'estimated_time_remaining', '3': 8, '4': 1, '5': 9, '9': 0, '10': 'estimatedTimeRemaining', '17': true},
  ],
  '8': [
    {'1': '_estimated_time_remaining'},
  ],
};

/// Descriptor for `CalculationProgress`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List calculationProgressDescriptor = $convert.base64Decode(
    'ChNDYWxjdWxhdGlvblByb2dyZXNzEiUKDmNhbGN1bGF0aW9uX2lkGAEgASgJUg1jYWxjdWxhdG'
    'lvbklkEhUKBmpvYl9pZBgCIAEoCVIFam9iSWQSNgoGc3RhdHVzGAMgASgOMh4ucHlzY2ZfZnJv'
    'bnQuQ2FsY3VsYXRpb25TdGF0dXNSBnN0YXR1cxIvChNwcm9ncmVzc19wZXJjZW50YWdlGAQgAS'
    'gBUhJwcm9ncmVzc1BlcmNlbnRhZ2USIQoMY3VycmVudF9zdGVwGAUgASgJUgtjdXJyZW50U3Rl'
    'cBIYCgdtZXNzYWdlGAYgASgJUgdtZXNzYWdlEhwKCXRpbWVzdGFtcBgHIAEoA1IJdGltZXN0YW'
    '1wEj0KGGVzdGltYXRlZF90aW1lX3JlbWFpbmluZxgIIAEoCUgAUhZlc3RpbWF0ZWRUaW1lUmVt'
    'YWluaW5niAEBQhsKGV9lc3RpbWF0ZWRfdGltZV9yZW1haW5pbmc=');

@$core.Deprecated('Use calculationDescriptor instead')
const Calculation$json = {
  '1': 'Calculation',
  '2': [
    {'1': 'id', '3': 1, '4': 1, '5': 9, '10': 'id'},
    {'1': 'instance_id', '3': 2, '4': 1, '5': 9, '10': 'instanceId'},
    {'1': 'method', '3': 3, '4': 1, '5': 9, '10': 'method'},
    {'1': 'basis_set', '3': 4, '4': 1, '5': 9, '10': 'basisSet'},
    {'1': 'parameters', '3': 5, '4': 3, '5': 11, '6': '.pyscf_front.Calculation.ParametersEntry', '10': 'parameters'},
    {'1': 'convergence_criteria', '3': 6, '4': 1, '5': 11, '6': '.pyscf_front.ConvergenceCriteria', '10': 'convergenceCriteria'},
    {'1': 'max_iterations', '3': 7, '4': 1, '5': 5, '10': 'maxIterations'},
    {'1': 'start_time', '3': 8, '4': 1, '5': 3, '9': 0, '10': 'startTime', '17': true},
    {'1': 'end_time', '3': 9, '4': 1, '5': 3, '9': 1, '10': 'endTime', '17': true},
    {'1': 'status', '3': 10, '4': 1, '5': 14, '6': '.pyscf_front.CalculationStatus', '10': 'status'},
    {'1': 'error_message', '3': 11, '4': 1, '5': 9, '9': 2, '10': 'errorMessage', '17': true},
    {'1': 'created_at', '3': 12, '4': 1, '5': 3, '10': 'createdAt'},
  ],
  '3': [Calculation_ParametersEntry$json],
  '8': [
    {'1': '_start_time'},
    {'1': '_end_time'},
    {'1': '_error_message'},
  ],
};

@$core.Deprecated('Use calculationDescriptor instead')
const Calculation_ParametersEntry$json = {
  '1': 'ParametersEntry',
  '2': [
    {'1': 'key', '3': 1, '4': 1, '5': 9, '10': 'key'},
    {'1': 'value', '3': 2, '4': 1, '5': 9, '10': 'value'},
  ],
  '7': {'7': true},
};

/// Descriptor for `Calculation`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List calculationDescriptor = $convert.base64Decode(
    'CgtDYWxjdWxhdGlvbhIOCgJpZBgBIAEoCVICaWQSHwoLaW5zdGFuY2VfaWQYAiABKAlSCmluc3'
    'RhbmNlSWQSFgoGbWV0aG9kGAMgASgJUgZtZXRob2QSGwoJYmFzaXNfc2V0GAQgASgJUghiYXNp'
    'c1NldBJICgpwYXJhbWV0ZXJzGAUgAygLMigucHlzY2ZfZnJvbnQuQ2FsY3VsYXRpb24uUGFyYW'
    '1ldGVyc0VudHJ5UgpwYXJhbWV0ZXJzElMKFGNvbnZlcmdlbmNlX2NyaXRlcmlhGAYgASgLMiAu'
    'cHlzY2ZfZnJvbnQuQ29udmVyZ2VuY2VDcml0ZXJpYVITY29udmVyZ2VuY2VDcml0ZXJpYRIlCg'
    '5tYXhfaXRlcmF0aW9ucxgHIAEoBVINbWF4SXRlcmF0aW9ucxIiCgpzdGFydF90aW1lGAggASgD'
    'SABSCXN0YXJ0VGltZYgBARIeCghlbmRfdGltZRgJIAEoA0gBUgdlbmRUaW1liAEBEjYKBnN0YX'
    'R1cxgKIAEoDjIeLnB5c2NmX2Zyb250LkNhbGN1bGF0aW9uU3RhdHVzUgZzdGF0dXMSKAoNZXJy'
    'b3JfbWVzc2FnZRgLIAEoCUgCUgxlcnJvck1lc3NhZ2WIAQESHQoKY3JlYXRlZF9hdBgMIAEoA1'
    'IJY3JlYXRlZEF0Gj0KD1BhcmFtZXRlcnNFbnRyeRIQCgNrZXkYASABKAlSA2tleRIUCgV2YWx1'
    'ZRgCIAEoCVIFdmFsdWU6AjgBQg0KC19zdGFydF90aW1lQgsKCV9lbmRfdGltZUIQCg5fZXJyb3'
    'JfbWVzc2FnZQ==');

@$core.Deprecated('Use convergenceCriteriaDescriptor instead')
const ConvergenceCriteria$json = {
  '1': 'ConvergenceCriteria',
  '2': [
    {'1': 'energy_threshold', '3': 1, '4': 1, '5': 1, '10': 'energyThreshold'},
    {'1': 'density_threshold', '3': 2, '4': 1, '5': 1, '10': 'densityThreshold'},
    {'1': 'gradient_threshold', '3': 3, '4': 1, '5': 1, '10': 'gradientThreshold'},
  ],
};

/// Descriptor for `ConvergenceCriteria`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List convergenceCriteriaDescriptor = $convert.base64Decode(
    'ChNDb252ZXJnZW5jZUNyaXRlcmlhEikKEGVuZXJneV90aHJlc2hvbGQYASABKAFSD2VuZXJneV'
    'RocmVzaG9sZBIrChFkZW5zaXR5X3RocmVzaG9sZBgCIAEoAVIQZGVuc2l0eVRocmVzaG9sZBIt'
    'ChJncmFkaWVudF90aHJlc2hvbGQYAyABKAFSEWdyYWRpZW50VGhyZXNob2xk');

@$core.Deprecated('Use getResultsRequestDescriptor instead')
const GetResultsRequest$json = {
  '1': 'GetResultsRequest',
  '2': [
    {'1': 'calculation_id', '3': 1, '4': 1, '5': 9, '10': 'calculationId'},
    {'1': 'result_types', '3': 2, '4': 3, '5': 9, '10': 'resultTypes'},
  ],
};

/// Descriptor for `GetResultsRequest`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List getResultsRequestDescriptor = $convert.base64Decode(
    'ChFHZXRSZXN1bHRzUmVxdWVzdBIlCg5jYWxjdWxhdGlvbl9pZBgBIAEoCVINY2FsY3VsYXRpb2'
    '5JZBIhCgxyZXN1bHRfdHlwZXMYAiADKAlSC3Jlc3VsdFR5cGVz');

@$core.Deprecated('Use exportResultsRequestDescriptor instead')
const ExportResultsRequest$json = {
  '1': 'ExportResultsRequest',
  '2': [
    {'1': 'calculation_id', '3': 1, '4': 1, '5': 9, '10': 'calculationId'},
    {'1': 'format', '3': 2, '4': 1, '5': 14, '6': '.pyscf_front.ExportFormat', '10': 'format'},
    {'1': 'result_types', '3': 3, '4': 3, '5': 9, '10': 'resultTypes'},
  ],
};

/// Descriptor for `ExportResultsRequest`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List exportResultsRequestDescriptor = $convert.base64Decode(
    'ChRFeHBvcnRSZXN1bHRzUmVxdWVzdBIlCg5jYWxjdWxhdGlvbl9pZBgBIAEoCVINY2FsY3VsYX'
    'Rpb25JZBIxCgZmb3JtYXQYAiABKA4yGS5weXNjZl9mcm9udC5FeHBvcnRGb3JtYXRSBmZvcm1h'
    'dBIhCgxyZXN1bHRfdHlwZXMYAyADKAlSC3Jlc3VsdFR5cGVz');

@$core.Deprecated('Use exportResultsResponseDescriptor instead')
const ExportResultsResponse$json = {
  '1': 'ExportResultsResponse',
  '2': [
    {'1': 'success', '3': 1, '4': 1, '5': 8, '10': 'success'},
    {'1': 'message', '3': 2, '4': 1, '5': 9, '10': 'message'},
    {'1': 'file_path', '3': 3, '4': 1, '5': 9, '9': 0, '10': 'filePath', '17': true},
    {'1': 'file_data', '3': 4, '4': 1, '5': 12, '9': 1, '10': 'fileData', '17': true},
  ],
  '8': [
    {'1': '_file_path'},
    {'1': '_file_data'},
  ],
};

/// Descriptor for `ExportResultsResponse`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List exportResultsResponseDescriptor = $convert.base64Decode(
    'ChVFeHBvcnRSZXN1bHRzUmVzcG9uc2USGAoHc3VjY2VzcxgBIAEoCFIHc3VjY2VzcxIYCgdtZX'
    'NzYWdlGAIgASgJUgdtZXNzYWdlEiAKCWZpbGVfcGF0aBgDIAEoCUgAUghmaWxlUGF0aIgBARIg'
    'CglmaWxlX2RhdGEYBCABKAxIAVIIZmlsZURhdGGIAQFCDAoKX2ZpbGVfcGF0aEIMCgpfZmlsZV'
    '9kYXRh');

@$core.Deprecated('Use resultsResponseDescriptor instead')
const ResultsResponse$json = {
  '1': 'ResultsResponse',
  '2': [
    {'1': 'success', '3': 1, '4': 1, '5': 8, '10': 'success'},
    {'1': 'message', '3': 2, '4': 1, '5': 9, '10': 'message'},
    {'1': 'results', '3': 3, '4': 1, '5': 11, '6': '.pyscf_front.CalculationResults', '9': 0, '10': 'results', '17': true},
  ],
  '8': [
    {'1': '_results'},
  ],
};

/// Descriptor for `ResultsResponse`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List resultsResponseDescriptor = $convert.base64Decode(
    'Cg9SZXN1bHRzUmVzcG9uc2USGAoHc3VjY2VzcxgBIAEoCFIHc3VjY2VzcxIYCgdtZXNzYWdlGA'
    'IgASgJUgdtZXNzYWdlEj4KB3Jlc3VsdHMYAyABKAsyHy5weXNjZl9mcm9udC5DYWxjdWxhdGlv'
    'blJlc3VsdHNIAFIHcmVzdWx0c4gBAUIKCghfcmVzdWx0cw==');

@$core.Deprecated('Use calculationResultsDescriptor instead')
const CalculationResults$json = {
  '1': 'CalculationResults',
  '2': [
    {'1': 'calculation_id', '3': 1, '4': 1, '5': 9, '10': 'calculationId'},
    {'1': 'converged', '3': 2, '4': 1, '5': 8, '10': 'converged'},
    {'1': 'total_energy', '3': 3, '4': 1, '5': 1, '10': 'totalEnergy'},
    {'1': 'energy_unit', '3': 4, '4': 1, '5': 9, '10': 'energyUnit'},
    {'1': 'calculation_time_seconds', '3': 5, '4': 1, '5': 5, '10': 'calculationTimeSeconds'},
    {'1': 'energy_results', '3': 6, '4': 1, '5': 11, '6': '.pyscf_front.EnergyResults', '9': 0, '10': 'energyResults', '17': true},
    {'1': 'orbital_results', '3': 7, '4': 1, '5': 11, '6': '.pyscf_front.OrbitalResults', '9': 1, '10': 'orbitalResults', '17': true},
    {'1': 'frequency_results', '3': 8, '4': 1, '5': 11, '6': '.pyscf_front.FrequencyResults', '9': 2, '10': 'frequencyResults', '17': true},
    {'1': 'property_results', '3': 9, '4': 1, '5': 11, '6': '.pyscf_front.PropertyResults', '9': 3, '10': 'propertyResults', '17': true},
  ],
  '8': [
    {'1': '_energy_results'},
    {'1': '_orbital_results'},
    {'1': '_frequency_results'},
    {'1': '_property_results'},
  ],
};

/// Descriptor for `CalculationResults`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List calculationResultsDescriptor = $convert.base64Decode(
    'ChJDYWxjdWxhdGlvblJlc3VsdHMSJQoOY2FsY3VsYXRpb25faWQYASABKAlSDWNhbGN1bGF0aW'
    '9uSWQSHAoJY29udmVyZ2VkGAIgASgIUgljb252ZXJnZWQSIQoMdG90YWxfZW5lcmd5GAMgASgB'
    'Ugt0b3RhbEVuZXJneRIfCgtlbmVyZ3lfdW5pdBgEIAEoCVIKZW5lcmd5VW5pdBI4ChhjYWxjdW'
    'xhdGlvbl90aW1lX3NlY29uZHMYBSABKAVSFmNhbGN1bGF0aW9uVGltZVNlY29uZHMSRgoOZW5l'
    'cmd5X3Jlc3VsdHMYBiABKAsyGi5weXNjZl9mcm9udC5FbmVyZ3lSZXN1bHRzSABSDWVuZXJneV'
    'Jlc3VsdHOIAQESSQoPb3JiaXRhbF9yZXN1bHRzGAcgASgLMhsucHlzY2ZfZnJvbnQuT3JiaXRh'
    'bFJlc3VsdHNIAVIOb3JiaXRhbFJlc3VsdHOIAQESTwoRZnJlcXVlbmN5X3Jlc3VsdHMYCCABKA'
    'syHS5weXNjZl9mcm9udC5GcmVxdWVuY3lSZXN1bHRzSAJSEGZyZXF1ZW5jeVJlc3VsdHOIAQES'
    'TAoQcHJvcGVydHlfcmVzdWx0cxgJIAEoCzIcLnB5c2NmX2Zyb250LlByb3BlcnR5UmVzdWx0c0'
    'gDUg9wcm9wZXJ0eVJlc3VsdHOIAQFCEQoPX2VuZXJneV9yZXN1bHRzQhIKEF9vcmJpdGFsX3Jl'
    'c3VsdHNCFAoSX2ZyZXF1ZW5jeV9yZXN1bHRzQhMKEV9wcm9wZXJ0eV9yZXN1bHRz');

@$core.Deprecated('Use energyResultsDescriptor instead')
const EnergyResults$json = {
  '1': 'EnergyResults',
  '2': [
    {'1': 'total_energy', '3': 1, '4': 1, '5': 1, '10': 'totalEnergy'},
    {'1': 'nuclear_repulsion_energy', '3': 2, '4': 1, '5': 1, '10': 'nuclearRepulsionEnergy'},
    {'1': 'electronic_energy', '3': 3, '4': 1, '5': 1, '10': 'electronicEnergy'},
    {'1': 'correlation_energy', '3': 4, '4': 1, '5': 1, '9': 0, '10': 'correlationEnergy', '17': true},
    {'1': 'components', '3': 5, '4': 3, '5': 11, '6': '.pyscf_front.EnergyComponent', '10': 'components'},
  ],
  '8': [
    {'1': '_correlation_energy'},
  ],
};

/// Descriptor for `EnergyResults`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List energyResultsDescriptor = $convert.base64Decode(
    'Cg1FbmVyZ3lSZXN1bHRzEiEKDHRvdGFsX2VuZXJneRgBIAEoAVILdG90YWxFbmVyZ3kSOAoYbn'
    'VjbGVhcl9yZXB1bHNpb25fZW5lcmd5GAIgASgBUhZudWNsZWFyUmVwdWxzaW9uRW5lcmd5EisK'
    'EWVsZWN0cm9uaWNfZW5lcmd5GAMgASgBUhBlbGVjdHJvbmljRW5lcmd5EjIKEmNvcnJlbGF0aW'
    '9uX2VuZXJneRgEIAEoAUgAUhFjb3JyZWxhdGlvbkVuZXJneYgBARI8Cgpjb21wb25lbnRzGAUg'
    'AygLMhwucHlzY2ZfZnJvbnQuRW5lcmd5Q29tcG9uZW50Ugpjb21wb25lbnRzQhUKE19jb3JyZW'
    'xhdGlvbl9lbmVyZ3k=');

@$core.Deprecated('Use energyComponentDescriptor instead')
const EnergyComponent$json = {
  '1': 'EnergyComponent',
  '2': [
    {'1': 'name', '3': 1, '4': 1, '5': 9, '10': 'name'},
    {'1': 'value', '3': 2, '4': 1, '5': 1, '10': 'value'},
    {'1': 'unit', '3': 3, '4': 1, '5': 9, '10': 'unit'},
  ],
};

/// Descriptor for `EnergyComponent`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List energyComponentDescriptor = $convert.base64Decode(
    'Cg9FbmVyZ3lDb21wb25lbnQSEgoEbmFtZRgBIAEoCVIEbmFtZRIUCgV2YWx1ZRgCIAEoAVIFdm'
    'FsdWUSEgoEdW5pdBgDIAEoCVIEdW5pdA==');

@$core.Deprecated('Use orbitalResultsDescriptor instead')
const OrbitalResults$json = {
  '1': 'OrbitalResults',
  '2': [
    {'1': 'homo_energy', '3': 1, '4': 1, '5': 1, '10': 'homoEnergy'},
    {'1': 'lumo_energy', '3': 2, '4': 1, '5': 1, '10': 'lumoEnergy'},
    {'1': 'homo_lumo_gap', '3': 3, '4': 1, '5': 1, '10': 'homoLumoGap'},
    {'1': 'orbital_energies', '3': 4, '4': 3, '5': 1, '10': 'orbitalEnergies'},
    {'1': 'occupations', '3': 5, '4': 3, '5': 1, '10': 'occupations'},
    {'1': 'molecular_orbitals', '3': 6, '4': 1, '5': 12, '9': 0, '10': 'molecularOrbitals', '17': true},
  ],
  '8': [
    {'1': '_molecular_orbitals'},
  ],
};

/// Descriptor for `OrbitalResults`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List orbitalResultsDescriptor = $convert.base64Decode(
    'Cg5PcmJpdGFsUmVzdWx0cxIfCgtob21vX2VuZXJneRgBIAEoAVIKaG9tb0VuZXJneRIfCgtsdW'
    '1vX2VuZXJneRgCIAEoAVIKbHVtb0VuZXJneRIiCg1ob21vX2x1bW9fZ2FwGAMgASgBUgtob21v'
    'THVtb0dhcBIpChBvcmJpdGFsX2VuZXJnaWVzGAQgAygBUg9vcmJpdGFsRW5lcmdpZXMSIAoLb2'
    'NjdXBhdGlvbnMYBSADKAFSC29jY3VwYXRpb25zEjIKEm1vbGVjdWxhcl9vcmJpdGFscxgGIAEo'
    'DEgAUhFtb2xlY3VsYXJPcmJpdGFsc4gBAUIVChNfbW9sZWN1bGFyX29yYml0YWxz');

@$core.Deprecated('Use frequencyResultsDescriptor instead')
const FrequencyResults$json = {
  '1': 'FrequencyResults',
  '2': [
    {'1': 'frequencies', '3': 1, '4': 3, '5': 1, '10': 'frequencies'},
    {'1': 'intensities', '3': 2, '4': 3, '5': 1, '10': 'intensities'},
    {'1': 'normal_modes', '3': 3, '4': 1, '5': 12, '9': 0, '10': 'normalModes', '17': true},
  ],
  '8': [
    {'1': '_normal_modes'},
  ],
};

/// Descriptor for `FrequencyResults`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List frequencyResultsDescriptor = $convert.base64Decode(
    'ChBGcmVxdWVuY3lSZXN1bHRzEiAKC2ZyZXF1ZW5jaWVzGAEgAygBUgtmcmVxdWVuY2llcxIgCg'
    'tpbnRlbnNpdGllcxgCIAMoAVILaW50ZW5zaXRpZXMSJgoMbm9ybWFsX21vZGVzGAMgASgMSABS'
    'C25vcm1hbE1vZGVziAEBQg8KDV9ub3JtYWxfbW9kZXM=');

@$core.Deprecated('Use propertyResultsDescriptor instead')
const PropertyResults$json = {
  '1': 'PropertyResults',
  '2': [
    {'1': 'dipole_moment', '3': 1, '4': 1, '5': 11, '6': '.pyscf_front.DipoleProperty', '9': 0, '10': 'dipoleMoment', '17': true},
    {'1': 'quadrupole_moment', '3': 2, '4': 1, '5': 11, '6': '.pyscf_front.QuadrupoleProperty', '9': 1, '10': 'quadrupoleMoment', '17': true},
    {'1': 'atomic_charges', '3': 3, '4': 3, '5': 11, '6': '.pyscf_front.ChargeProperty', '10': 'atomicCharges'},
  ],
  '8': [
    {'1': '_dipole_moment'},
    {'1': '_quadrupole_moment'},
  ],
};

/// Descriptor for `PropertyResults`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List propertyResultsDescriptor = $convert.base64Decode(
    'Cg9Qcm9wZXJ0eVJlc3VsdHMSRQoNZGlwb2xlX21vbWVudBgBIAEoCzIbLnB5c2NmX2Zyb250Lk'
    'RpcG9sZVByb3BlcnR5SABSDGRpcG9sZU1vbWVudIgBARJRChFxdWFkcnVwb2xlX21vbWVudBgC'
    'IAEoCzIfLnB5c2NmX2Zyb250LlF1YWRydXBvbGVQcm9wZXJ0eUgBUhBxdWFkcnVwb2xlTW9tZW'
    '50iAEBEkIKDmF0b21pY19jaGFyZ2VzGAMgAygLMhsucHlzY2ZfZnJvbnQuQ2hhcmdlUHJvcGVy'
    'dHlSDWF0b21pY0NoYXJnZXNCEAoOX2RpcG9sZV9tb21lbnRCFAoSX3F1YWRydXBvbGVfbW9tZW'
    '50');

@$core.Deprecated('Use dipolePropertyDescriptor instead')
const DipoleProperty$json = {
  '1': 'DipoleProperty',
  '2': [
    {'1': 'magnitude', '3': 1, '4': 1, '5': 1, '10': 'magnitude'},
    {'1': 'x', '3': 2, '4': 1, '5': 1, '10': 'x'},
    {'1': 'y', '3': 3, '4': 1, '5': 1, '10': 'y'},
    {'1': 'z', '3': 4, '4': 1, '5': 1, '10': 'z'},
    {'1': 'unit', '3': 5, '4': 1, '5': 9, '10': 'unit'},
  ],
};

/// Descriptor for `DipoleProperty`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List dipolePropertyDescriptor = $convert.base64Decode(
    'Cg5EaXBvbGVQcm9wZXJ0eRIcCgltYWduaXR1ZGUYASABKAFSCW1hZ25pdHVkZRIMCgF4GAIgAS'
    'gBUgF4EgwKAXkYAyABKAFSAXkSDAoBehgEIAEoAVIBehISCgR1bml0GAUgASgJUgR1bml0');

@$core.Deprecated('Use quadrupolePropertyDescriptor instead')
const QuadrupoleProperty$json = {
  '1': 'QuadrupoleProperty',
  '2': [
    {'1': 'tensor', '3': 1, '4': 3, '5': 1, '10': 'tensor'},
    {'1': 'unit', '3': 2, '4': 1, '5': 9, '10': 'unit'},
  ],
};

/// Descriptor for `QuadrupoleProperty`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List quadrupolePropertyDescriptor = $convert.base64Decode(
    'ChJRdWFkcnVwb2xlUHJvcGVydHkSFgoGdGVuc29yGAEgAygBUgZ0ZW5zb3ISEgoEdW5pdBgCIA'
    'EoCVIEdW5pdA==');

@$core.Deprecated('Use chargePropertyDescriptor instead')
const ChargeProperty$json = {
  '1': 'ChargeProperty',
  '2': [
    {'1': 'atom_index', '3': 1, '4': 1, '5': 5, '10': 'atomIndex'},
    {'1': 'method', '3': 2, '4': 1, '5': 9, '10': 'method'},
    {'1': 'charge', '3': 3, '4': 1, '5': 1, '10': 'charge'},
  ],
};

/// Descriptor for `ChargeProperty`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List chargePropertyDescriptor = $convert.base64Decode(
    'Cg5DaGFyZ2VQcm9wZXJ0eRIdCgphdG9tX2luZGV4GAEgASgFUglhdG9tSW5kZXgSFgoGbWV0aG'
    '9kGAIgASgJUgZtZXRob2QSFgoGY2hhcmdlGAMgASgBUgZjaGFyZ2U=');

@$core.Deprecated('Use getJobQueueRequestDescriptor instead')
const GetJobQueueRequest$json = {
  '1': 'GetJobQueueRequest',
  '2': [
    {'1': 'status_filter', '3': 1, '4': 1, '5': 14, '6': '.pyscf_front.JobStatus', '9': 0, '10': 'statusFilter', '17': true},
  ],
  '8': [
    {'1': '_status_filter'},
  ],
};

/// Descriptor for `GetJobQueueRequest`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List getJobQueueRequestDescriptor = $convert.base64Decode(
    'ChJHZXRKb2JRdWV1ZVJlcXVlc3QSQAoNc3RhdHVzX2ZpbHRlchgBIAEoDjIWLnB5c2NmX2Zyb2'
    '50LkpvYlN0YXR1c0gAUgxzdGF0dXNGaWx0ZXKIAQFCEAoOX3N0YXR1c19maWx0ZXI=');

@$core.Deprecated('Use updateJobPriorityRequestDescriptor instead')
const UpdateJobPriorityRequest$json = {
  '1': 'UpdateJobPriorityRequest',
  '2': [
    {'1': 'job_id', '3': 1, '4': 1, '5': 9, '10': 'jobId'},
    {'1': 'new_priority', '3': 2, '4': 1, '5': 5, '10': 'newPriority'},
  ],
};

/// Descriptor for `UpdateJobPriorityRequest`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List updateJobPriorityRequestDescriptor = $convert.base64Decode(
    'ChhVcGRhdGVKb2JQcmlvcml0eVJlcXVlc3QSFQoGam9iX2lkGAEgASgJUgVqb2JJZBIhCgxuZX'
    'dfcHJpb3JpdHkYAiABKAVSC25ld1ByaW9yaXR5');

@$core.Deprecated('Use updateJobPriorityResponseDescriptor instead')
const UpdateJobPriorityResponse$json = {
  '1': 'UpdateJobPriorityResponse',
  '2': [
    {'1': 'success', '3': 1, '4': 1, '5': 8, '10': 'success'},
    {'1': 'message', '3': 2, '4': 1, '5': 9, '10': 'message'},
  ],
};

/// Descriptor for `UpdateJobPriorityResponse`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List updateJobPriorityResponseDescriptor = $convert.base64Decode(
    'ChlVcGRhdGVKb2JQcmlvcml0eVJlc3BvbnNlEhgKB3N1Y2Nlc3MYASABKAhSB3N1Y2Nlc3MSGA'
    'oHbWVzc2FnZRgCIAEoCVIHbWVzc2FnZQ==');

@$core.Deprecated('Use jobQueueResponseDescriptor instead')
const JobQueueResponse$json = {
  '1': 'JobQueueResponse',
  '2': [
    {'1': 'success', '3': 1, '4': 1, '5': 8, '10': 'success'},
    {'1': 'message', '3': 2, '4': 1, '5': 9, '10': 'message'},
    {'1': 'jobs', '3': 3, '4': 3, '5': 11, '6': '.pyscf_front.JobInfo', '10': 'jobs'},
  ],
};

/// Descriptor for `JobQueueResponse`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List jobQueueResponseDescriptor = $convert.base64Decode(
    'ChBKb2JRdWV1ZVJlc3BvbnNlEhgKB3N1Y2Nlc3MYASABKAhSB3N1Y2Nlc3MSGAoHbWVzc2FnZR'
    'gCIAEoCVIHbWVzc2FnZRIoCgRqb2JzGAMgAygLMhQucHlzY2ZfZnJvbnQuSm9iSW5mb1IEam9i'
    'cw==');

@$core.Deprecated('Use jobInfoDescriptor instead')
const JobInfo$json = {
  '1': 'JobInfo',
  '2': [
    {'1': 'id', '3': 1, '4': 1, '5': 9, '10': 'id'},
    {'1': 'calculation_id', '3': 2, '4': 1, '5': 9, '10': 'calculationId'},
    {'1': 'priority', '3': 3, '4': 1, '5': 5, '10': 'priority'},
    {'1': 'status', '3': 4, '4': 1, '5': 14, '6': '.pyscf_front.JobStatus', '10': 'status'},
    {'1': 'assigned_worker', '3': 5, '4': 1, '5': 9, '9': 0, '10': 'assignedWorker', '17': true},
    {'1': 'created_at', '3': 6, '4': 1, '5': 3, '10': 'createdAt'},
    {'1': 'started_at', '3': 7, '4': 1, '5': 3, '9': 1, '10': 'startedAt', '17': true},
    {'1': 'completed_at', '3': 8, '4': 1, '5': 3, '9': 2, '10': 'completedAt', '17': true},
  ],
  '8': [
    {'1': '_assigned_worker'},
    {'1': '_started_at'},
    {'1': '_completed_at'},
  ],
};

/// Descriptor for `JobInfo`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List jobInfoDescriptor = $convert.base64Decode(
    'CgdKb2JJbmZvEg4KAmlkGAEgASgJUgJpZBIlCg5jYWxjdWxhdGlvbl9pZBgCIAEoCVINY2FsY3'
    'VsYXRpb25JZBIaCghwcmlvcml0eRgDIAEoBVIIcHJpb3JpdHkSLgoGc3RhdHVzGAQgASgOMhYu'
    'cHlzY2ZfZnJvbnQuSm9iU3RhdHVzUgZzdGF0dXMSLAoPYXNzaWduZWRfd29ya2VyGAUgASgJSA'
    'BSDmFzc2lnbmVkV29ya2VyiAEBEh0KCmNyZWF0ZWRfYXQYBiABKANSCWNyZWF0ZWRBdBIiCgpz'
    'dGFydGVkX2F0GAcgASgDSAFSCXN0YXJ0ZWRBdIgBARImCgxjb21wbGV0ZWRfYXQYCCABKANIAl'
    'ILY29tcGxldGVkQXSIAQFCEgoQX2Fzc2lnbmVkX3dvcmtlckINCgtfc3RhcnRlZF9hdEIPCg1f'
    'Y29tcGxldGVkX2F0');

@$core.Deprecated('Use getSystemInfoRequestDescriptor instead')
const GetSystemInfoRequest$json = {
  '1': 'GetSystemInfoRequest',
};

/// Descriptor for `GetSystemInfoRequest`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List getSystemInfoRequestDescriptor = $convert.base64Decode(
    'ChRHZXRTeXN0ZW1JbmZvUmVxdWVzdA==');

@$core.Deprecated('Use systemInfoResponseDescriptor instead')
const SystemInfoResponse$json = {
  '1': 'SystemInfoResponse',
  '2': [
    {'1': 'success', '3': 1, '4': 1, '5': 8, '10': 'success'},
    {'1': 'message', '3': 2, '4': 1, '5': 9, '10': 'message'},
    {'1': 'system_info', '3': 3, '4': 1, '5': 11, '6': '.pyscf_front.SystemInfo', '9': 0, '10': 'systemInfo', '17': true},
  ],
  '8': [
    {'1': '_system_info'},
  ],
};

/// Descriptor for `SystemInfoResponse`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List systemInfoResponseDescriptor = $convert.base64Decode(
    'ChJTeXN0ZW1JbmZvUmVzcG9uc2USGAoHc3VjY2VzcxgBIAEoCFIHc3VjY2VzcxIYCgdtZXNzYW'
    'dlGAIgASgJUgdtZXNzYWdlEj0KC3N5c3RlbV9pbmZvGAMgASgLMhcucHlzY2ZfZnJvbnQuU3lz'
    'dGVtSW5mb0gAUgpzeXN0ZW1JbmZviAEBQg4KDF9zeXN0ZW1faW5mbw==');

@$core.Deprecated('Use systemInfoDescriptor instead')
const SystemInfo$json = {
  '1': 'SystemInfo',
  '2': [
    {'1': 'version', '3': 1, '4': 1, '5': 9, '10': 'version'},
    {'1': 'python_version', '3': 2, '4': 1, '5': 9, '10': 'pythonVersion'},
    {'1': 'pyscf_version', '3': 3, '4': 1, '5': 9, '10': 'pyscfVersion'},
    {'1': 'gpu_available', '3': 4, '4': 1, '5': 8, '10': 'gpuAvailable'},
    {'1': 'available_methods', '3': 5, '4': 3, '5': 9, '10': 'availableMethods'},
    {'1': 'available_basis_sets', '3': 6, '4': 3, '5': 9, '10': 'availableBasisSets'},
    {'1': 'resources', '3': 7, '4': 1, '5': 11, '6': '.pyscf_front.SystemResources', '10': 'resources'},
  ],
};

/// Descriptor for `SystemInfo`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List systemInfoDescriptor = $convert.base64Decode(
    'CgpTeXN0ZW1JbmZvEhgKB3ZlcnNpb24YASABKAlSB3ZlcnNpb24SJQoOcHl0aG9uX3ZlcnNpb2'
    '4YAiABKAlSDXB5dGhvblZlcnNpb24SIwoNcHlzY2ZfdmVyc2lvbhgDIAEoCVIMcHlzY2ZWZXJz'
    'aW9uEiMKDWdwdV9hdmFpbGFibGUYBCABKAhSDGdwdUF2YWlsYWJsZRIrChFhdmFpbGFibGVfbW'
    'V0aG9kcxgFIAMoCVIQYXZhaWxhYmxlTWV0aG9kcxIwChRhdmFpbGFibGVfYmFzaXNfc2V0cxgG'
    'IAMoCVISYXZhaWxhYmxlQmFzaXNTZXRzEjoKCXJlc291cmNlcxgHIAEoCzIcLnB5c2NmX2Zyb2'
    '50LlN5c3RlbVJlc291cmNlc1IJcmVzb3VyY2Vz');

@$core.Deprecated('Use systemResourcesDescriptor instead')
const SystemResources$json = {
  '1': 'SystemResources',
  '2': [
    {'1': 'cpu_cores', '3': 1, '4': 1, '5': 5, '10': 'cpuCores'},
    {'1': 'total_memory_mb', '3': 2, '4': 1, '5': 3, '10': 'totalMemoryMb'},
    {'1': 'available_memory_mb', '3': 3, '4': 1, '5': 3, '10': 'availableMemoryMb'},
    {'1': 'disk_space_mb', '3': 4, '4': 1, '5': 3, '10': 'diskSpaceMb'},
    {'1': 'gpus', '3': 5, '4': 3, '5': 11, '6': '.pyscf_front.GpuInfo', '10': 'gpus'},
  ],
};

/// Descriptor for `SystemResources`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List systemResourcesDescriptor = $convert.base64Decode(
    'Cg9TeXN0ZW1SZXNvdXJjZXMSGwoJY3B1X2NvcmVzGAEgASgFUghjcHVDb3JlcxImCg90b3RhbF'
    '9tZW1vcnlfbWIYAiABKANSDXRvdGFsTWVtb3J5TWISLgoTYXZhaWxhYmxlX21lbW9yeV9tYhgD'
    'IAEoA1IRYXZhaWxhYmxlTWVtb3J5TWISIgoNZGlza19zcGFjZV9tYhgEIAEoA1ILZGlza1NwYW'
    'NlTWISKAoEZ3B1cxgFIAMoCzIULnB5c2NmX2Zyb250LkdwdUluZm9SBGdwdXM=');

@$core.Deprecated('Use gpuInfoDescriptor instead')
const GpuInfo$json = {
  '1': 'GpuInfo',
  '2': [
    {'1': 'name', '3': 1, '4': 1, '5': 9, '10': 'name'},
    {'1': 'memory_mb', '3': 2, '4': 1, '5': 3, '10': 'memoryMb'},
    {'1': 'driver_version', '3': 3, '4': 1, '5': 9, '10': 'driverVersion'},
  ],
};

/// Descriptor for `GpuInfo`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List gpuInfoDescriptor = $convert.base64Decode(
    'CgdHcHVJbmZvEhIKBG5hbWUYASABKAlSBG5hbWUSGwoJbWVtb3J5X21iGAIgASgDUghtZW1vcn'
    'lNYhIlCg5kcml2ZXJfdmVyc2lvbhgDIAEoCVINZHJpdmVyVmVyc2lvbg==');

@$core.Deprecated('Use healthCheckRequestDescriptor instead')
const HealthCheckRequest$json = {
  '1': 'HealthCheckRequest',
};

/// Descriptor for `HealthCheckRequest`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List healthCheckRequestDescriptor = $convert.base64Decode(
    'ChJIZWFsdGhDaGVja1JlcXVlc3Q=');

@$core.Deprecated('Use healthCheckResponseDescriptor instead')
const HealthCheckResponse$json = {
  '1': 'HealthCheckResponse',
  '2': [
    {'1': 'healthy', '3': 1, '4': 1, '5': 8, '10': 'healthy'},
    {'1': 'status', '3': 2, '4': 1, '5': 9, '10': 'status'},
    {'1': 'details', '3': 3, '4': 3, '5': 11, '6': '.pyscf_front.HealthCheckResponse.DetailsEntry', '10': 'details'},
    {'1': 'timestamp', '3': 4, '4': 1, '5': 3, '10': 'timestamp'},
  ],
  '3': [HealthCheckResponse_DetailsEntry$json],
};

@$core.Deprecated('Use healthCheckResponseDescriptor instead')
const HealthCheckResponse_DetailsEntry$json = {
  '1': 'DetailsEntry',
  '2': [
    {'1': 'key', '3': 1, '4': 1, '5': 9, '10': 'key'},
    {'1': 'value', '3': 2, '4': 1, '5': 9, '10': 'value'},
  ],
  '7': {'7': true},
};

/// Descriptor for `HealthCheckResponse`. Decode as a `google.protobuf.DescriptorProto`.
final $typed_data.Uint8List healthCheckResponseDescriptor = $convert.base64Decode(
    'ChNIZWFsdGhDaGVja1Jlc3BvbnNlEhgKB2hlYWx0aHkYASABKAhSB2hlYWx0aHkSFgoGc3RhdH'
    'VzGAIgASgJUgZzdGF0dXMSRwoHZGV0YWlscxgDIAMoCzItLnB5c2NmX2Zyb250LkhlYWx0aENo'
    'ZWNrUmVzcG9uc2UuRGV0YWlsc0VudHJ5UgdkZXRhaWxzEhwKCXRpbWVzdGFtcBgEIAEoA1IJdG'
    'ltZXN0YW1wGjoKDERldGFpbHNFbnRyeRIQCgNrZXkYASABKAlSA2tleRIUCgV2YWx1ZRgCIAEo'
    'CVIFdmFsdWU6AjgB');

