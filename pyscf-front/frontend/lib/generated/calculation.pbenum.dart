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

import 'dart:core' as $core;

import 'package:protobuf/protobuf.dart' as $pb;

class GeometryType extends $pb.ProtobufEnum {
  static const GeometryType GEOMETRY_TYPE_UNSPECIFIED = GeometryType._(0, _omitEnumNames ? '' : 'GEOMETRY_TYPE_UNSPECIFIED');
  static const GeometryType GEOMETRY_TYPE_XYZ = GeometryType._(1, _omitEnumNames ? '' : 'GEOMETRY_TYPE_XYZ');
  static const GeometryType GEOMETRY_TYPE_ZMATRIX = GeometryType._(2, _omitEnumNames ? '' : 'GEOMETRY_TYPE_ZMATRIX');

  static const $core.List<GeometryType> values = <GeometryType> [
    GEOMETRY_TYPE_UNSPECIFIED,
    GEOMETRY_TYPE_XYZ,
    GEOMETRY_TYPE_ZMATRIX,
  ];

  static final $core.List<GeometryType?> _byValue = $pb.ProtobufEnum.$_initByValueList(values, 2);
  static GeometryType? valueOf($core.int value) =>  value < 0 || value >= _byValue.length ? null : _byValue[value];

  const GeometryType._(super.value, super.name);
}

class InstanceStatus extends $pb.ProtobufEnum {
  static const InstanceStatus INSTANCE_STATUS_UNSPECIFIED = InstanceStatus._(0, _omitEnumNames ? '' : 'INSTANCE_STATUS_UNSPECIFIED');
  static const InstanceStatus INSTANCE_STATUS_DRAFT = InstanceStatus._(1, _omitEnumNames ? '' : 'INSTANCE_STATUS_DRAFT');
  static const InstanceStatus INSTANCE_STATUS_READY = InstanceStatus._(2, _omitEnumNames ? '' : 'INSTANCE_STATUS_READY');
  static const InstanceStatus INSTANCE_STATUS_RUNNING = InstanceStatus._(3, _omitEnumNames ? '' : 'INSTANCE_STATUS_RUNNING');
  static const InstanceStatus INSTANCE_STATUS_COMPLETED = InstanceStatus._(4, _omitEnumNames ? '' : 'INSTANCE_STATUS_COMPLETED');
  static const InstanceStatus INSTANCE_STATUS_ERROR = InstanceStatus._(5, _omitEnumNames ? '' : 'INSTANCE_STATUS_ERROR');

  static const $core.List<InstanceStatus> values = <InstanceStatus> [
    INSTANCE_STATUS_UNSPECIFIED,
    INSTANCE_STATUS_DRAFT,
    INSTANCE_STATUS_READY,
    INSTANCE_STATUS_RUNNING,
    INSTANCE_STATUS_COMPLETED,
    INSTANCE_STATUS_ERROR,
  ];

  static final $core.List<InstanceStatus?> _byValue = $pb.ProtobufEnum.$_initByValueList(values, 5);
  static InstanceStatus? valueOf($core.int value) =>  value < 0 || value >= _byValue.length ? null : _byValue[value];

  const InstanceStatus._(super.value, super.name);
}

class CalculationStatus extends $pb.ProtobufEnum {
  static const CalculationStatus CALCULATION_STATUS_UNSPECIFIED = CalculationStatus._(0, _omitEnumNames ? '' : 'CALCULATION_STATUS_UNSPECIFIED');
  static const CalculationStatus CALCULATION_STATUS_PENDING = CalculationStatus._(1, _omitEnumNames ? '' : 'CALCULATION_STATUS_PENDING');
  static const CalculationStatus CALCULATION_STATUS_RUNNING = CalculationStatus._(2, _omitEnumNames ? '' : 'CALCULATION_STATUS_RUNNING');
  static const CalculationStatus CALCULATION_STATUS_COMPLETED = CalculationStatus._(3, _omitEnumNames ? '' : 'CALCULATION_STATUS_COMPLETED');
  static const CalculationStatus CALCULATION_STATUS_FAILED = CalculationStatus._(4, _omitEnumNames ? '' : 'CALCULATION_STATUS_FAILED');
  static const CalculationStatus CALCULATION_STATUS_CANCELLED = CalculationStatus._(5, _omitEnumNames ? '' : 'CALCULATION_STATUS_CANCELLED');

  static const $core.List<CalculationStatus> values = <CalculationStatus> [
    CALCULATION_STATUS_UNSPECIFIED,
    CALCULATION_STATUS_PENDING,
    CALCULATION_STATUS_RUNNING,
    CALCULATION_STATUS_COMPLETED,
    CALCULATION_STATUS_FAILED,
    CALCULATION_STATUS_CANCELLED,
  ];

  static final $core.List<CalculationStatus?> _byValue = $pb.ProtobufEnum.$_initByValueList(values, 5);
  static CalculationStatus? valueOf($core.int value) =>  value < 0 || value >= _byValue.length ? null : _byValue[value];

  const CalculationStatus._(super.value, super.name);
}

class ExportFormat extends $pb.ProtobufEnum {
  static const ExportFormat EXPORT_FORMAT_UNSPECIFIED = ExportFormat._(0, _omitEnumNames ? '' : 'EXPORT_FORMAT_UNSPECIFIED');
  static const ExportFormat EXPORT_FORMAT_JSON = ExportFormat._(1, _omitEnumNames ? '' : 'EXPORT_FORMAT_JSON');
  static const ExportFormat EXPORT_FORMAT_XML = ExportFormat._(2, _omitEnumNames ? '' : 'EXPORT_FORMAT_XML');
  static const ExportFormat EXPORT_FORMAT_HDF5 = ExportFormat._(3, _omitEnumNames ? '' : 'EXPORT_FORMAT_HDF5');
  static const ExportFormat EXPORT_FORMAT_PDF = ExportFormat._(4, _omitEnumNames ? '' : 'EXPORT_FORMAT_PDF');

  static const $core.List<ExportFormat> values = <ExportFormat> [
    EXPORT_FORMAT_UNSPECIFIED,
    EXPORT_FORMAT_JSON,
    EXPORT_FORMAT_XML,
    EXPORT_FORMAT_HDF5,
    EXPORT_FORMAT_PDF,
  ];

  static final $core.List<ExportFormat?> _byValue = $pb.ProtobufEnum.$_initByValueList(values, 4);
  static ExportFormat? valueOf($core.int value) =>  value < 0 || value >= _byValue.length ? null : _byValue[value];

  const ExportFormat._(super.value, super.name);
}

class JobStatus extends $pb.ProtobufEnum {
  static const JobStatus JOB_STATUS_UNSPECIFIED = JobStatus._(0, _omitEnumNames ? '' : 'JOB_STATUS_UNSPECIFIED');
  static const JobStatus JOB_STATUS_WAITING = JobStatus._(1, _omitEnumNames ? '' : 'JOB_STATUS_WAITING');
  static const JobStatus JOB_STATUS_RUNNING = JobStatus._(2, _omitEnumNames ? '' : 'JOB_STATUS_RUNNING');
  static const JobStatus JOB_STATUS_COMPLETED = JobStatus._(3, _omitEnumNames ? '' : 'JOB_STATUS_COMPLETED');
  static const JobStatus JOB_STATUS_FAILED = JobStatus._(4, _omitEnumNames ? '' : 'JOB_STATUS_FAILED');

  static const $core.List<JobStatus> values = <JobStatus> [
    JOB_STATUS_UNSPECIFIED,
    JOB_STATUS_WAITING,
    JOB_STATUS_RUNNING,
    JOB_STATUS_COMPLETED,
    JOB_STATUS_FAILED,
  ];

  static final $core.List<JobStatus?> _byValue = $pb.ProtobufEnum.$_initByValueList(values, 4);
  static JobStatus? valueOf($core.int value) =>  value < 0 || value >= _byValue.length ? null : _byValue[value];

  const JobStatus._(super.value, super.name);
}


const $core.bool _omitEnumNames = $core.bool.fromEnvironment('protobuf.omit_enum_names');
