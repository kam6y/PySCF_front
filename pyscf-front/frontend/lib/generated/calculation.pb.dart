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

import 'package:fixnum/fixnum.dart' as $fixnum;
import 'package:protobuf/protobuf.dart' as $pb;

import 'calculation.pbenum.dart';

export 'package:protobuf/protobuf.dart' show GeneratedMessageGenericExtensions;

export 'calculation.pbenum.dart';

/// 分子関連メッセージ
class CreateMoleculeRequest extends $pb.GeneratedMessage {
  factory CreateMoleculeRequest({
    $core.String? name,
    $core.String? formula,
    $core.Iterable<Atom>? atoms,
    $core.int? charge,
    $core.int? multiplicity,
    $core.String? symmetry,
    GeometryType? geometryType,
  }) {
    final result = create();
    if (name != null) result.name = name;
    if (formula != null) result.formula = formula;
    if (atoms != null) result.atoms.addAll(atoms);
    if (charge != null) result.charge = charge;
    if (multiplicity != null) result.multiplicity = multiplicity;
    if (symmetry != null) result.symmetry = symmetry;
    if (geometryType != null) result.geometryType = geometryType;
    return result;
  }

  CreateMoleculeRequest._();

  factory CreateMoleculeRequest.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory CreateMoleculeRequest.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'CreateMoleculeRequest', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOS(1, _omitFieldNames ? '' : 'name')
    ..aOS(2, _omitFieldNames ? '' : 'formula')
    ..pc<Atom>(3, _omitFieldNames ? '' : 'atoms', $pb.PbFieldType.PM, subBuilder: Atom.create)
    ..a<$core.int>(4, _omitFieldNames ? '' : 'charge', $pb.PbFieldType.O3)
    ..a<$core.int>(5, _omitFieldNames ? '' : 'multiplicity', $pb.PbFieldType.O3)
    ..aOS(6, _omitFieldNames ? '' : 'symmetry')
    ..e<GeometryType>(7, _omitFieldNames ? '' : 'geometryType', $pb.PbFieldType.OE, defaultOrMaker: GeometryType.GEOMETRY_TYPE_UNSPECIFIED, valueOf: GeometryType.valueOf, enumValues: GeometryType.values)
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  CreateMoleculeRequest clone() => CreateMoleculeRequest()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  CreateMoleculeRequest copyWith(void Function(CreateMoleculeRequest) updates) => super.copyWith((message) => updates(message as CreateMoleculeRequest)) as CreateMoleculeRequest;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static CreateMoleculeRequest create() => CreateMoleculeRequest._();
  @$core.override
  CreateMoleculeRequest createEmptyInstance() => create();
  static $pb.PbList<CreateMoleculeRequest> createRepeated() => $pb.PbList<CreateMoleculeRequest>();
  @$core.pragma('dart2js:noInline')
  static CreateMoleculeRequest getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<CreateMoleculeRequest>(create);
  static CreateMoleculeRequest? _defaultInstance;

  @$pb.TagNumber(1)
  $core.String get name => $_getSZ(0);
  @$pb.TagNumber(1)
  set name($core.String value) => $_setString(0, value);
  @$pb.TagNumber(1)
  $core.bool hasName() => $_has(0);
  @$pb.TagNumber(1)
  void clearName() => $_clearField(1);

  @$pb.TagNumber(2)
  $core.String get formula => $_getSZ(1);
  @$pb.TagNumber(2)
  set formula($core.String value) => $_setString(1, value);
  @$pb.TagNumber(2)
  $core.bool hasFormula() => $_has(1);
  @$pb.TagNumber(2)
  void clearFormula() => $_clearField(2);

  @$pb.TagNumber(3)
  $pb.PbList<Atom> get atoms => $_getList(2);

  @$pb.TagNumber(4)
  $core.int get charge => $_getIZ(3);
  @$pb.TagNumber(4)
  set charge($core.int value) => $_setSignedInt32(3, value);
  @$pb.TagNumber(4)
  $core.bool hasCharge() => $_has(3);
  @$pb.TagNumber(4)
  void clearCharge() => $_clearField(4);

  @$pb.TagNumber(5)
  $core.int get multiplicity => $_getIZ(4);
  @$pb.TagNumber(5)
  set multiplicity($core.int value) => $_setSignedInt32(4, value);
  @$pb.TagNumber(5)
  $core.bool hasMultiplicity() => $_has(4);
  @$pb.TagNumber(5)
  void clearMultiplicity() => $_clearField(5);

  @$pb.TagNumber(6)
  $core.String get symmetry => $_getSZ(5);
  @$pb.TagNumber(6)
  set symmetry($core.String value) => $_setString(5, value);
  @$pb.TagNumber(6)
  $core.bool hasSymmetry() => $_has(5);
  @$pb.TagNumber(6)
  void clearSymmetry() => $_clearField(6);

  @$pb.TagNumber(7)
  GeometryType get geometryType => $_getN(6);
  @$pb.TagNumber(7)
  set geometryType(GeometryType value) => $_setField(7, value);
  @$pb.TagNumber(7)
  $core.bool hasGeometryType() => $_has(6);
  @$pb.TagNumber(7)
  void clearGeometryType() => $_clearField(7);
}

class GetMoleculeRequest extends $pb.GeneratedMessage {
  factory GetMoleculeRequest({
    $core.String? moleculeId,
  }) {
    final result = create();
    if (moleculeId != null) result.moleculeId = moleculeId;
    return result;
  }

  GetMoleculeRequest._();

  factory GetMoleculeRequest.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory GetMoleculeRequest.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'GetMoleculeRequest', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOS(1, _omitFieldNames ? '' : 'moleculeId')
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  GetMoleculeRequest clone() => GetMoleculeRequest()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  GetMoleculeRequest copyWith(void Function(GetMoleculeRequest) updates) => super.copyWith((message) => updates(message as GetMoleculeRequest)) as GetMoleculeRequest;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static GetMoleculeRequest create() => GetMoleculeRequest._();
  @$core.override
  GetMoleculeRequest createEmptyInstance() => create();
  static $pb.PbList<GetMoleculeRequest> createRepeated() => $pb.PbList<GetMoleculeRequest>();
  @$core.pragma('dart2js:noInline')
  static GetMoleculeRequest getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<GetMoleculeRequest>(create);
  static GetMoleculeRequest? _defaultInstance;

  @$pb.TagNumber(1)
  $core.String get moleculeId => $_getSZ(0);
  @$pb.TagNumber(1)
  set moleculeId($core.String value) => $_setString(0, value);
  @$pb.TagNumber(1)
  $core.bool hasMoleculeId() => $_has(0);
  @$pb.TagNumber(1)
  void clearMoleculeId() => $_clearField(1);
}

class UpdateMoleculeRequest extends $pb.GeneratedMessage {
  factory UpdateMoleculeRequest({
    $core.String? moleculeId,
    $core.String? name,
    $core.String? formula,
    $core.Iterable<Atom>? atoms,
    $core.int? charge,
    $core.int? multiplicity,
    $core.String? symmetry,
  }) {
    final result = create();
    if (moleculeId != null) result.moleculeId = moleculeId;
    if (name != null) result.name = name;
    if (formula != null) result.formula = formula;
    if (atoms != null) result.atoms.addAll(atoms);
    if (charge != null) result.charge = charge;
    if (multiplicity != null) result.multiplicity = multiplicity;
    if (symmetry != null) result.symmetry = symmetry;
    return result;
  }

  UpdateMoleculeRequest._();

  factory UpdateMoleculeRequest.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory UpdateMoleculeRequest.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'UpdateMoleculeRequest', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOS(1, _omitFieldNames ? '' : 'moleculeId')
    ..aOS(2, _omitFieldNames ? '' : 'name')
    ..aOS(3, _omitFieldNames ? '' : 'formula')
    ..pc<Atom>(4, _omitFieldNames ? '' : 'atoms', $pb.PbFieldType.PM, subBuilder: Atom.create)
    ..a<$core.int>(5, _omitFieldNames ? '' : 'charge', $pb.PbFieldType.O3)
    ..a<$core.int>(6, _omitFieldNames ? '' : 'multiplicity', $pb.PbFieldType.O3)
    ..aOS(7, _omitFieldNames ? '' : 'symmetry')
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  UpdateMoleculeRequest clone() => UpdateMoleculeRequest()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  UpdateMoleculeRequest copyWith(void Function(UpdateMoleculeRequest) updates) => super.copyWith((message) => updates(message as UpdateMoleculeRequest)) as UpdateMoleculeRequest;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static UpdateMoleculeRequest create() => UpdateMoleculeRequest._();
  @$core.override
  UpdateMoleculeRequest createEmptyInstance() => create();
  static $pb.PbList<UpdateMoleculeRequest> createRepeated() => $pb.PbList<UpdateMoleculeRequest>();
  @$core.pragma('dart2js:noInline')
  static UpdateMoleculeRequest getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<UpdateMoleculeRequest>(create);
  static UpdateMoleculeRequest? _defaultInstance;

  @$pb.TagNumber(1)
  $core.String get moleculeId => $_getSZ(0);
  @$pb.TagNumber(1)
  set moleculeId($core.String value) => $_setString(0, value);
  @$pb.TagNumber(1)
  $core.bool hasMoleculeId() => $_has(0);
  @$pb.TagNumber(1)
  void clearMoleculeId() => $_clearField(1);

  @$pb.TagNumber(2)
  $core.String get name => $_getSZ(1);
  @$pb.TagNumber(2)
  set name($core.String value) => $_setString(1, value);
  @$pb.TagNumber(2)
  $core.bool hasName() => $_has(1);
  @$pb.TagNumber(2)
  void clearName() => $_clearField(2);

  @$pb.TagNumber(3)
  $core.String get formula => $_getSZ(2);
  @$pb.TagNumber(3)
  set formula($core.String value) => $_setString(2, value);
  @$pb.TagNumber(3)
  $core.bool hasFormula() => $_has(2);
  @$pb.TagNumber(3)
  void clearFormula() => $_clearField(3);

  @$pb.TagNumber(4)
  $pb.PbList<Atom> get atoms => $_getList(3);

  @$pb.TagNumber(5)
  $core.int get charge => $_getIZ(4);
  @$pb.TagNumber(5)
  set charge($core.int value) => $_setSignedInt32(4, value);
  @$pb.TagNumber(5)
  $core.bool hasCharge() => $_has(4);
  @$pb.TagNumber(5)
  void clearCharge() => $_clearField(5);

  @$pb.TagNumber(6)
  $core.int get multiplicity => $_getIZ(5);
  @$pb.TagNumber(6)
  set multiplicity($core.int value) => $_setSignedInt32(5, value);
  @$pb.TagNumber(6)
  $core.bool hasMultiplicity() => $_has(5);
  @$pb.TagNumber(6)
  void clearMultiplicity() => $_clearField(6);

  @$pb.TagNumber(7)
  $core.String get symmetry => $_getSZ(6);
  @$pb.TagNumber(7)
  set symmetry($core.String value) => $_setString(6, value);
  @$pb.TagNumber(7)
  $core.bool hasSymmetry() => $_has(6);
  @$pb.TagNumber(7)
  void clearSymmetry() => $_clearField(7);
}

class DeleteMoleculeRequest extends $pb.GeneratedMessage {
  factory DeleteMoleculeRequest({
    $core.String? moleculeId,
  }) {
    final result = create();
    if (moleculeId != null) result.moleculeId = moleculeId;
    return result;
  }

  DeleteMoleculeRequest._();

  factory DeleteMoleculeRequest.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory DeleteMoleculeRequest.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'DeleteMoleculeRequest', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOS(1, _omitFieldNames ? '' : 'moleculeId')
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  DeleteMoleculeRequest clone() => DeleteMoleculeRequest()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  DeleteMoleculeRequest copyWith(void Function(DeleteMoleculeRequest) updates) => super.copyWith((message) => updates(message as DeleteMoleculeRequest)) as DeleteMoleculeRequest;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static DeleteMoleculeRequest create() => DeleteMoleculeRequest._();
  @$core.override
  DeleteMoleculeRequest createEmptyInstance() => create();
  static $pb.PbList<DeleteMoleculeRequest> createRepeated() => $pb.PbList<DeleteMoleculeRequest>();
  @$core.pragma('dart2js:noInline')
  static DeleteMoleculeRequest getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<DeleteMoleculeRequest>(create);
  static DeleteMoleculeRequest? _defaultInstance;

  @$pb.TagNumber(1)
  $core.String get moleculeId => $_getSZ(0);
  @$pb.TagNumber(1)
  set moleculeId($core.String value) => $_setString(0, value);
  @$pb.TagNumber(1)
  $core.bool hasMoleculeId() => $_has(0);
  @$pb.TagNumber(1)
  void clearMoleculeId() => $_clearField(1);
}

class DeleteMoleculeResponse extends $pb.GeneratedMessage {
  factory DeleteMoleculeResponse({
    $core.bool? success,
    $core.String? message,
  }) {
    final result = create();
    if (success != null) result.success = success;
    if (message != null) result.message = message;
    return result;
  }

  DeleteMoleculeResponse._();

  factory DeleteMoleculeResponse.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory DeleteMoleculeResponse.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'DeleteMoleculeResponse', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOB(1, _omitFieldNames ? '' : 'success')
    ..aOS(2, _omitFieldNames ? '' : 'message')
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  DeleteMoleculeResponse clone() => DeleteMoleculeResponse()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  DeleteMoleculeResponse copyWith(void Function(DeleteMoleculeResponse) updates) => super.copyWith((message) => updates(message as DeleteMoleculeResponse)) as DeleteMoleculeResponse;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static DeleteMoleculeResponse create() => DeleteMoleculeResponse._();
  @$core.override
  DeleteMoleculeResponse createEmptyInstance() => create();
  static $pb.PbList<DeleteMoleculeResponse> createRepeated() => $pb.PbList<DeleteMoleculeResponse>();
  @$core.pragma('dart2js:noInline')
  static DeleteMoleculeResponse getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<DeleteMoleculeResponse>(create);
  static DeleteMoleculeResponse? _defaultInstance;

  @$pb.TagNumber(1)
  $core.bool get success => $_getBF(0);
  @$pb.TagNumber(1)
  set success($core.bool value) => $_setBool(0, value);
  @$pb.TagNumber(1)
  $core.bool hasSuccess() => $_has(0);
  @$pb.TagNumber(1)
  void clearSuccess() => $_clearField(1);

  @$pb.TagNumber(2)
  $core.String get message => $_getSZ(1);
  @$pb.TagNumber(2)
  set message($core.String value) => $_setString(1, value);
  @$pb.TagNumber(2)
  $core.bool hasMessage() => $_has(1);
  @$pb.TagNumber(2)
  void clearMessage() => $_clearField(2);
}

class MoleculeResponse extends $pb.GeneratedMessage {
  factory MoleculeResponse({
    $core.bool? success,
    $core.String? message,
    Molecule? molecule,
  }) {
    final result = create();
    if (success != null) result.success = success;
    if (message != null) result.message = message;
    if (molecule != null) result.molecule = molecule;
    return result;
  }

  MoleculeResponse._();

  factory MoleculeResponse.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory MoleculeResponse.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'MoleculeResponse', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOB(1, _omitFieldNames ? '' : 'success')
    ..aOS(2, _omitFieldNames ? '' : 'message')
    ..aOM<Molecule>(3, _omitFieldNames ? '' : 'molecule', subBuilder: Molecule.create)
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  MoleculeResponse clone() => MoleculeResponse()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  MoleculeResponse copyWith(void Function(MoleculeResponse) updates) => super.copyWith((message) => updates(message as MoleculeResponse)) as MoleculeResponse;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static MoleculeResponse create() => MoleculeResponse._();
  @$core.override
  MoleculeResponse createEmptyInstance() => create();
  static $pb.PbList<MoleculeResponse> createRepeated() => $pb.PbList<MoleculeResponse>();
  @$core.pragma('dart2js:noInline')
  static MoleculeResponse getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<MoleculeResponse>(create);
  static MoleculeResponse? _defaultInstance;

  @$pb.TagNumber(1)
  $core.bool get success => $_getBF(0);
  @$pb.TagNumber(1)
  set success($core.bool value) => $_setBool(0, value);
  @$pb.TagNumber(1)
  $core.bool hasSuccess() => $_has(0);
  @$pb.TagNumber(1)
  void clearSuccess() => $_clearField(1);

  @$pb.TagNumber(2)
  $core.String get message => $_getSZ(1);
  @$pb.TagNumber(2)
  set message($core.String value) => $_setString(1, value);
  @$pb.TagNumber(2)
  $core.bool hasMessage() => $_has(1);
  @$pb.TagNumber(2)
  void clearMessage() => $_clearField(2);

  @$pb.TagNumber(3)
  Molecule get molecule => $_getN(2);
  @$pb.TagNumber(3)
  set molecule(Molecule value) => $_setField(3, value);
  @$pb.TagNumber(3)
  $core.bool hasMolecule() => $_has(2);
  @$pb.TagNumber(3)
  void clearMolecule() => $_clearField(3);
  @$pb.TagNumber(3)
  Molecule ensureMolecule() => $_ensure(2);
}

class Molecule extends $pb.GeneratedMessage {
  factory Molecule({
    $core.String? id,
    $core.String? name,
    $core.String? formula,
    $core.double? molecularWeight,
    $core.Iterable<Atom>? atoms,
    $core.int? charge,
    $core.int? multiplicity,
    $core.String? symmetry,
    GeometryType? geometryType,
    $fixnum.Int64? createdAt,
  }) {
    final result = create();
    if (id != null) result.id = id;
    if (name != null) result.name = name;
    if (formula != null) result.formula = formula;
    if (molecularWeight != null) result.molecularWeight = molecularWeight;
    if (atoms != null) result.atoms.addAll(atoms);
    if (charge != null) result.charge = charge;
    if (multiplicity != null) result.multiplicity = multiplicity;
    if (symmetry != null) result.symmetry = symmetry;
    if (geometryType != null) result.geometryType = geometryType;
    if (createdAt != null) result.createdAt = createdAt;
    return result;
  }

  Molecule._();

  factory Molecule.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory Molecule.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'Molecule', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOS(1, _omitFieldNames ? '' : 'id')
    ..aOS(2, _omitFieldNames ? '' : 'name')
    ..aOS(3, _omitFieldNames ? '' : 'formula')
    ..a<$core.double>(4, _omitFieldNames ? '' : 'molecularWeight', $pb.PbFieldType.OD)
    ..pc<Atom>(5, _omitFieldNames ? '' : 'atoms', $pb.PbFieldType.PM, subBuilder: Atom.create)
    ..a<$core.int>(6, _omitFieldNames ? '' : 'charge', $pb.PbFieldType.O3)
    ..a<$core.int>(7, _omitFieldNames ? '' : 'multiplicity', $pb.PbFieldType.O3)
    ..aOS(8, _omitFieldNames ? '' : 'symmetry')
    ..e<GeometryType>(9, _omitFieldNames ? '' : 'geometryType', $pb.PbFieldType.OE, defaultOrMaker: GeometryType.GEOMETRY_TYPE_UNSPECIFIED, valueOf: GeometryType.valueOf, enumValues: GeometryType.values)
    ..aInt64(10, _omitFieldNames ? '' : 'createdAt')
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  Molecule clone() => Molecule()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  Molecule copyWith(void Function(Molecule) updates) => super.copyWith((message) => updates(message as Molecule)) as Molecule;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static Molecule create() => Molecule._();
  @$core.override
  Molecule createEmptyInstance() => create();
  static $pb.PbList<Molecule> createRepeated() => $pb.PbList<Molecule>();
  @$core.pragma('dart2js:noInline')
  static Molecule getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<Molecule>(create);
  static Molecule? _defaultInstance;

  @$pb.TagNumber(1)
  $core.String get id => $_getSZ(0);
  @$pb.TagNumber(1)
  set id($core.String value) => $_setString(0, value);
  @$pb.TagNumber(1)
  $core.bool hasId() => $_has(0);
  @$pb.TagNumber(1)
  void clearId() => $_clearField(1);

  @$pb.TagNumber(2)
  $core.String get name => $_getSZ(1);
  @$pb.TagNumber(2)
  set name($core.String value) => $_setString(1, value);
  @$pb.TagNumber(2)
  $core.bool hasName() => $_has(1);
  @$pb.TagNumber(2)
  void clearName() => $_clearField(2);

  @$pb.TagNumber(3)
  $core.String get formula => $_getSZ(2);
  @$pb.TagNumber(3)
  set formula($core.String value) => $_setString(2, value);
  @$pb.TagNumber(3)
  $core.bool hasFormula() => $_has(2);
  @$pb.TagNumber(3)
  void clearFormula() => $_clearField(3);

  @$pb.TagNumber(4)
  $core.double get molecularWeight => $_getN(3);
  @$pb.TagNumber(4)
  set molecularWeight($core.double value) => $_setDouble(3, value);
  @$pb.TagNumber(4)
  $core.bool hasMolecularWeight() => $_has(3);
  @$pb.TagNumber(4)
  void clearMolecularWeight() => $_clearField(4);

  @$pb.TagNumber(5)
  $pb.PbList<Atom> get atoms => $_getList(4);

  @$pb.TagNumber(6)
  $core.int get charge => $_getIZ(5);
  @$pb.TagNumber(6)
  set charge($core.int value) => $_setSignedInt32(5, value);
  @$pb.TagNumber(6)
  $core.bool hasCharge() => $_has(5);
  @$pb.TagNumber(6)
  void clearCharge() => $_clearField(6);

  @$pb.TagNumber(7)
  $core.int get multiplicity => $_getIZ(6);
  @$pb.TagNumber(7)
  set multiplicity($core.int value) => $_setSignedInt32(6, value);
  @$pb.TagNumber(7)
  $core.bool hasMultiplicity() => $_has(6);
  @$pb.TagNumber(7)
  void clearMultiplicity() => $_clearField(7);

  @$pb.TagNumber(8)
  $core.String get symmetry => $_getSZ(7);
  @$pb.TagNumber(8)
  set symmetry($core.String value) => $_setString(7, value);
  @$pb.TagNumber(8)
  $core.bool hasSymmetry() => $_has(7);
  @$pb.TagNumber(8)
  void clearSymmetry() => $_clearField(8);

  @$pb.TagNumber(9)
  GeometryType get geometryType => $_getN(8);
  @$pb.TagNumber(9)
  set geometryType(GeometryType value) => $_setField(9, value);
  @$pb.TagNumber(9)
  $core.bool hasGeometryType() => $_has(8);
  @$pb.TagNumber(9)
  void clearGeometryType() => $_clearField(9);

  @$pb.TagNumber(10)
  $fixnum.Int64 get createdAt => $_getI64(9);
  @$pb.TagNumber(10)
  set createdAt($fixnum.Int64 value) => $_setInt64(9, value);
  @$pb.TagNumber(10)
  $core.bool hasCreatedAt() => $_has(9);
  @$pb.TagNumber(10)
  void clearCreatedAt() => $_clearField(10);
}

class Atom extends $pb.GeneratedMessage {
  factory Atom({
    $core.String? symbol,
    $core.double? x,
    $core.double? y,
    $core.double? z,
  }) {
    final result = create();
    if (symbol != null) result.symbol = symbol;
    if (x != null) result.x = x;
    if (y != null) result.y = y;
    if (z != null) result.z = z;
    return result;
  }

  Atom._();

  factory Atom.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory Atom.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'Atom', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOS(1, _omitFieldNames ? '' : 'symbol')
    ..a<$core.double>(2, _omitFieldNames ? '' : 'x', $pb.PbFieldType.OD)
    ..a<$core.double>(3, _omitFieldNames ? '' : 'y', $pb.PbFieldType.OD)
    ..a<$core.double>(4, _omitFieldNames ? '' : 'z', $pb.PbFieldType.OD)
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  Atom clone() => Atom()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  Atom copyWith(void Function(Atom) updates) => super.copyWith((message) => updates(message as Atom)) as Atom;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static Atom create() => Atom._();
  @$core.override
  Atom createEmptyInstance() => create();
  static $pb.PbList<Atom> createRepeated() => $pb.PbList<Atom>();
  @$core.pragma('dart2js:noInline')
  static Atom getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<Atom>(create);
  static Atom? _defaultInstance;

  @$pb.TagNumber(1)
  $core.String get symbol => $_getSZ(0);
  @$pb.TagNumber(1)
  set symbol($core.String value) => $_setString(0, value);
  @$pb.TagNumber(1)
  $core.bool hasSymbol() => $_has(0);
  @$pb.TagNumber(1)
  void clearSymbol() => $_clearField(1);

  @$pb.TagNumber(2)
  $core.double get x => $_getN(1);
  @$pb.TagNumber(2)
  set x($core.double value) => $_setDouble(1, value);
  @$pb.TagNumber(2)
  $core.bool hasX() => $_has(1);
  @$pb.TagNumber(2)
  void clearX() => $_clearField(2);

  @$pb.TagNumber(3)
  $core.double get y => $_getN(2);
  @$pb.TagNumber(3)
  set y($core.double value) => $_setDouble(2, value);
  @$pb.TagNumber(3)
  $core.bool hasY() => $_has(2);
  @$pb.TagNumber(3)
  void clearY() => $_clearField(3);

  @$pb.TagNumber(4)
  $core.double get z => $_getN(3);
  @$pb.TagNumber(4)
  set z($core.double value) => $_setDouble(3, value);
  @$pb.TagNumber(4)
  $core.bool hasZ() => $_has(3);
  @$pb.TagNumber(4)
  void clearZ() => $_clearField(4);
}

/// インスタンス関連メッセージ
class CreateInstanceRequest extends $pb.GeneratedMessage {
  factory CreateInstanceRequest({
    $core.String? name,
    $core.String? description,
    $core.String? moleculeId,
    $core.String? projectId,
  }) {
    final result = create();
    if (name != null) result.name = name;
    if (description != null) result.description = description;
    if (moleculeId != null) result.moleculeId = moleculeId;
    if (projectId != null) result.projectId = projectId;
    return result;
  }

  CreateInstanceRequest._();

  factory CreateInstanceRequest.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory CreateInstanceRequest.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'CreateInstanceRequest', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOS(1, _omitFieldNames ? '' : 'name')
    ..aOS(2, _omitFieldNames ? '' : 'description')
    ..aOS(3, _omitFieldNames ? '' : 'moleculeId')
    ..aOS(4, _omitFieldNames ? '' : 'projectId')
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  CreateInstanceRequest clone() => CreateInstanceRequest()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  CreateInstanceRequest copyWith(void Function(CreateInstanceRequest) updates) => super.copyWith((message) => updates(message as CreateInstanceRequest)) as CreateInstanceRequest;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static CreateInstanceRequest create() => CreateInstanceRequest._();
  @$core.override
  CreateInstanceRequest createEmptyInstance() => create();
  static $pb.PbList<CreateInstanceRequest> createRepeated() => $pb.PbList<CreateInstanceRequest>();
  @$core.pragma('dart2js:noInline')
  static CreateInstanceRequest getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<CreateInstanceRequest>(create);
  static CreateInstanceRequest? _defaultInstance;

  @$pb.TagNumber(1)
  $core.String get name => $_getSZ(0);
  @$pb.TagNumber(1)
  set name($core.String value) => $_setString(0, value);
  @$pb.TagNumber(1)
  $core.bool hasName() => $_has(0);
  @$pb.TagNumber(1)
  void clearName() => $_clearField(1);

  @$pb.TagNumber(2)
  $core.String get description => $_getSZ(1);
  @$pb.TagNumber(2)
  set description($core.String value) => $_setString(1, value);
  @$pb.TagNumber(2)
  $core.bool hasDescription() => $_has(1);
  @$pb.TagNumber(2)
  void clearDescription() => $_clearField(2);

  @$pb.TagNumber(3)
  $core.String get moleculeId => $_getSZ(2);
  @$pb.TagNumber(3)
  set moleculeId($core.String value) => $_setString(2, value);
  @$pb.TagNumber(3)
  $core.bool hasMoleculeId() => $_has(2);
  @$pb.TagNumber(3)
  void clearMoleculeId() => $_clearField(3);

  @$pb.TagNumber(4)
  $core.String get projectId => $_getSZ(3);
  @$pb.TagNumber(4)
  set projectId($core.String value) => $_setString(3, value);
  @$pb.TagNumber(4)
  $core.bool hasProjectId() => $_has(3);
  @$pb.TagNumber(4)
  void clearProjectId() => $_clearField(4);
}

class GetInstanceRequest extends $pb.GeneratedMessage {
  factory GetInstanceRequest({
    $core.String? instanceId,
  }) {
    final result = create();
    if (instanceId != null) result.instanceId = instanceId;
    return result;
  }

  GetInstanceRequest._();

  factory GetInstanceRequest.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory GetInstanceRequest.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'GetInstanceRequest', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOS(1, _omitFieldNames ? '' : 'instanceId')
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  GetInstanceRequest clone() => GetInstanceRequest()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  GetInstanceRequest copyWith(void Function(GetInstanceRequest) updates) => super.copyWith((message) => updates(message as GetInstanceRequest)) as GetInstanceRequest;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static GetInstanceRequest create() => GetInstanceRequest._();
  @$core.override
  GetInstanceRequest createEmptyInstance() => create();
  static $pb.PbList<GetInstanceRequest> createRepeated() => $pb.PbList<GetInstanceRequest>();
  @$core.pragma('dart2js:noInline')
  static GetInstanceRequest getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<GetInstanceRequest>(create);
  static GetInstanceRequest? _defaultInstance;

  @$pb.TagNumber(1)
  $core.String get instanceId => $_getSZ(0);
  @$pb.TagNumber(1)
  set instanceId($core.String value) => $_setString(0, value);
  @$pb.TagNumber(1)
  $core.bool hasInstanceId() => $_has(0);
  @$pb.TagNumber(1)
  void clearInstanceId() => $_clearField(1);
}

class ListInstancesRequest extends $pb.GeneratedMessage {
  factory ListInstancesRequest({
    $core.String? projectId,
    InstanceStatus? status,
    $core.int? page,
    $core.int? pageSize,
  }) {
    final result = create();
    if (projectId != null) result.projectId = projectId;
    if (status != null) result.status = status;
    if (page != null) result.page = page;
    if (pageSize != null) result.pageSize = pageSize;
    return result;
  }

  ListInstancesRequest._();

  factory ListInstancesRequest.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory ListInstancesRequest.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'ListInstancesRequest', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOS(1, _omitFieldNames ? '' : 'projectId')
    ..e<InstanceStatus>(2, _omitFieldNames ? '' : 'status', $pb.PbFieldType.OE, defaultOrMaker: InstanceStatus.INSTANCE_STATUS_UNSPECIFIED, valueOf: InstanceStatus.valueOf, enumValues: InstanceStatus.values)
    ..a<$core.int>(3, _omitFieldNames ? '' : 'page', $pb.PbFieldType.O3)
    ..a<$core.int>(4, _omitFieldNames ? '' : 'pageSize', $pb.PbFieldType.O3)
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  ListInstancesRequest clone() => ListInstancesRequest()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  ListInstancesRequest copyWith(void Function(ListInstancesRequest) updates) => super.copyWith((message) => updates(message as ListInstancesRequest)) as ListInstancesRequest;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static ListInstancesRequest create() => ListInstancesRequest._();
  @$core.override
  ListInstancesRequest createEmptyInstance() => create();
  static $pb.PbList<ListInstancesRequest> createRepeated() => $pb.PbList<ListInstancesRequest>();
  @$core.pragma('dart2js:noInline')
  static ListInstancesRequest getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<ListInstancesRequest>(create);
  static ListInstancesRequest? _defaultInstance;

  @$pb.TagNumber(1)
  $core.String get projectId => $_getSZ(0);
  @$pb.TagNumber(1)
  set projectId($core.String value) => $_setString(0, value);
  @$pb.TagNumber(1)
  $core.bool hasProjectId() => $_has(0);
  @$pb.TagNumber(1)
  void clearProjectId() => $_clearField(1);

  @$pb.TagNumber(2)
  InstanceStatus get status => $_getN(1);
  @$pb.TagNumber(2)
  set status(InstanceStatus value) => $_setField(2, value);
  @$pb.TagNumber(2)
  $core.bool hasStatus() => $_has(1);
  @$pb.TagNumber(2)
  void clearStatus() => $_clearField(2);

  @$pb.TagNumber(3)
  $core.int get page => $_getIZ(2);
  @$pb.TagNumber(3)
  set page($core.int value) => $_setSignedInt32(2, value);
  @$pb.TagNumber(3)
  $core.bool hasPage() => $_has(2);
  @$pb.TagNumber(3)
  void clearPage() => $_clearField(3);

  @$pb.TagNumber(4)
  $core.int get pageSize => $_getIZ(3);
  @$pb.TagNumber(4)
  set pageSize($core.int value) => $_setSignedInt32(3, value);
  @$pb.TagNumber(4)
  $core.bool hasPageSize() => $_has(3);
  @$pb.TagNumber(4)
  void clearPageSize() => $_clearField(4);
}

class ListInstancesResponse extends $pb.GeneratedMessage {
  factory ListInstancesResponse({
    $core.bool? success,
    $core.String? message,
    $core.Iterable<Instance>? instances,
    $core.int? totalCount,
    $core.int? page,
    $core.int? pageSize,
  }) {
    final result = create();
    if (success != null) result.success = success;
    if (message != null) result.message = message;
    if (instances != null) result.instances.addAll(instances);
    if (totalCount != null) result.totalCount = totalCount;
    if (page != null) result.page = page;
    if (pageSize != null) result.pageSize = pageSize;
    return result;
  }

  ListInstancesResponse._();

  factory ListInstancesResponse.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory ListInstancesResponse.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'ListInstancesResponse', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOB(1, _omitFieldNames ? '' : 'success')
    ..aOS(2, _omitFieldNames ? '' : 'message')
    ..pc<Instance>(3, _omitFieldNames ? '' : 'instances', $pb.PbFieldType.PM, subBuilder: Instance.create)
    ..a<$core.int>(4, _omitFieldNames ? '' : 'totalCount', $pb.PbFieldType.O3)
    ..a<$core.int>(5, _omitFieldNames ? '' : 'page', $pb.PbFieldType.O3)
    ..a<$core.int>(6, _omitFieldNames ? '' : 'pageSize', $pb.PbFieldType.O3)
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  ListInstancesResponse clone() => ListInstancesResponse()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  ListInstancesResponse copyWith(void Function(ListInstancesResponse) updates) => super.copyWith((message) => updates(message as ListInstancesResponse)) as ListInstancesResponse;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static ListInstancesResponse create() => ListInstancesResponse._();
  @$core.override
  ListInstancesResponse createEmptyInstance() => create();
  static $pb.PbList<ListInstancesResponse> createRepeated() => $pb.PbList<ListInstancesResponse>();
  @$core.pragma('dart2js:noInline')
  static ListInstancesResponse getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<ListInstancesResponse>(create);
  static ListInstancesResponse? _defaultInstance;

  @$pb.TagNumber(1)
  $core.bool get success => $_getBF(0);
  @$pb.TagNumber(1)
  set success($core.bool value) => $_setBool(0, value);
  @$pb.TagNumber(1)
  $core.bool hasSuccess() => $_has(0);
  @$pb.TagNumber(1)
  void clearSuccess() => $_clearField(1);

  @$pb.TagNumber(2)
  $core.String get message => $_getSZ(1);
  @$pb.TagNumber(2)
  set message($core.String value) => $_setString(1, value);
  @$pb.TagNumber(2)
  $core.bool hasMessage() => $_has(1);
  @$pb.TagNumber(2)
  void clearMessage() => $_clearField(2);

  @$pb.TagNumber(3)
  $pb.PbList<Instance> get instances => $_getList(2);

  @$pb.TagNumber(4)
  $core.int get totalCount => $_getIZ(3);
  @$pb.TagNumber(4)
  set totalCount($core.int value) => $_setSignedInt32(3, value);
  @$pb.TagNumber(4)
  $core.bool hasTotalCount() => $_has(3);
  @$pb.TagNumber(4)
  void clearTotalCount() => $_clearField(4);

  @$pb.TagNumber(5)
  $core.int get page => $_getIZ(4);
  @$pb.TagNumber(5)
  set page($core.int value) => $_setSignedInt32(4, value);
  @$pb.TagNumber(5)
  $core.bool hasPage() => $_has(4);
  @$pb.TagNumber(5)
  void clearPage() => $_clearField(5);

  @$pb.TagNumber(6)
  $core.int get pageSize => $_getIZ(5);
  @$pb.TagNumber(6)
  set pageSize($core.int value) => $_setSignedInt32(5, value);
  @$pb.TagNumber(6)
  $core.bool hasPageSize() => $_has(5);
  @$pb.TagNumber(6)
  void clearPageSize() => $_clearField(6);
}

class UpdateInstanceRequest extends $pb.GeneratedMessage {
  factory UpdateInstanceRequest({
    $core.String? instanceId,
    $core.String? name,
    $core.String? description,
    InstanceStatus? status,
  }) {
    final result = create();
    if (instanceId != null) result.instanceId = instanceId;
    if (name != null) result.name = name;
    if (description != null) result.description = description;
    if (status != null) result.status = status;
    return result;
  }

  UpdateInstanceRequest._();

  factory UpdateInstanceRequest.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory UpdateInstanceRequest.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'UpdateInstanceRequest', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOS(1, _omitFieldNames ? '' : 'instanceId')
    ..aOS(2, _omitFieldNames ? '' : 'name')
    ..aOS(3, _omitFieldNames ? '' : 'description')
    ..e<InstanceStatus>(4, _omitFieldNames ? '' : 'status', $pb.PbFieldType.OE, defaultOrMaker: InstanceStatus.INSTANCE_STATUS_UNSPECIFIED, valueOf: InstanceStatus.valueOf, enumValues: InstanceStatus.values)
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  UpdateInstanceRequest clone() => UpdateInstanceRequest()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  UpdateInstanceRequest copyWith(void Function(UpdateInstanceRequest) updates) => super.copyWith((message) => updates(message as UpdateInstanceRequest)) as UpdateInstanceRequest;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static UpdateInstanceRequest create() => UpdateInstanceRequest._();
  @$core.override
  UpdateInstanceRequest createEmptyInstance() => create();
  static $pb.PbList<UpdateInstanceRequest> createRepeated() => $pb.PbList<UpdateInstanceRequest>();
  @$core.pragma('dart2js:noInline')
  static UpdateInstanceRequest getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<UpdateInstanceRequest>(create);
  static UpdateInstanceRequest? _defaultInstance;

  @$pb.TagNumber(1)
  $core.String get instanceId => $_getSZ(0);
  @$pb.TagNumber(1)
  set instanceId($core.String value) => $_setString(0, value);
  @$pb.TagNumber(1)
  $core.bool hasInstanceId() => $_has(0);
  @$pb.TagNumber(1)
  void clearInstanceId() => $_clearField(1);

  @$pb.TagNumber(2)
  $core.String get name => $_getSZ(1);
  @$pb.TagNumber(2)
  set name($core.String value) => $_setString(1, value);
  @$pb.TagNumber(2)
  $core.bool hasName() => $_has(1);
  @$pb.TagNumber(2)
  void clearName() => $_clearField(2);

  @$pb.TagNumber(3)
  $core.String get description => $_getSZ(2);
  @$pb.TagNumber(3)
  set description($core.String value) => $_setString(2, value);
  @$pb.TagNumber(3)
  $core.bool hasDescription() => $_has(2);
  @$pb.TagNumber(3)
  void clearDescription() => $_clearField(3);

  @$pb.TagNumber(4)
  InstanceStatus get status => $_getN(3);
  @$pb.TagNumber(4)
  set status(InstanceStatus value) => $_setField(4, value);
  @$pb.TagNumber(4)
  $core.bool hasStatus() => $_has(3);
  @$pb.TagNumber(4)
  void clearStatus() => $_clearField(4);
}

class DeleteInstanceRequest extends $pb.GeneratedMessage {
  factory DeleteInstanceRequest({
    $core.String? instanceId,
  }) {
    final result = create();
    if (instanceId != null) result.instanceId = instanceId;
    return result;
  }

  DeleteInstanceRequest._();

  factory DeleteInstanceRequest.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory DeleteInstanceRequest.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'DeleteInstanceRequest', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOS(1, _omitFieldNames ? '' : 'instanceId')
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  DeleteInstanceRequest clone() => DeleteInstanceRequest()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  DeleteInstanceRequest copyWith(void Function(DeleteInstanceRequest) updates) => super.copyWith((message) => updates(message as DeleteInstanceRequest)) as DeleteInstanceRequest;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static DeleteInstanceRequest create() => DeleteInstanceRequest._();
  @$core.override
  DeleteInstanceRequest createEmptyInstance() => create();
  static $pb.PbList<DeleteInstanceRequest> createRepeated() => $pb.PbList<DeleteInstanceRequest>();
  @$core.pragma('dart2js:noInline')
  static DeleteInstanceRequest getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<DeleteInstanceRequest>(create);
  static DeleteInstanceRequest? _defaultInstance;

  @$pb.TagNumber(1)
  $core.String get instanceId => $_getSZ(0);
  @$pb.TagNumber(1)
  set instanceId($core.String value) => $_setString(0, value);
  @$pb.TagNumber(1)
  $core.bool hasInstanceId() => $_has(0);
  @$pb.TagNumber(1)
  void clearInstanceId() => $_clearField(1);
}

class DeleteInstanceResponse extends $pb.GeneratedMessage {
  factory DeleteInstanceResponse({
    $core.bool? success,
    $core.String? message,
  }) {
    final result = create();
    if (success != null) result.success = success;
    if (message != null) result.message = message;
    return result;
  }

  DeleteInstanceResponse._();

  factory DeleteInstanceResponse.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory DeleteInstanceResponse.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'DeleteInstanceResponse', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOB(1, _omitFieldNames ? '' : 'success')
    ..aOS(2, _omitFieldNames ? '' : 'message')
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  DeleteInstanceResponse clone() => DeleteInstanceResponse()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  DeleteInstanceResponse copyWith(void Function(DeleteInstanceResponse) updates) => super.copyWith((message) => updates(message as DeleteInstanceResponse)) as DeleteInstanceResponse;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static DeleteInstanceResponse create() => DeleteInstanceResponse._();
  @$core.override
  DeleteInstanceResponse createEmptyInstance() => create();
  static $pb.PbList<DeleteInstanceResponse> createRepeated() => $pb.PbList<DeleteInstanceResponse>();
  @$core.pragma('dart2js:noInline')
  static DeleteInstanceResponse getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<DeleteInstanceResponse>(create);
  static DeleteInstanceResponse? _defaultInstance;

  @$pb.TagNumber(1)
  $core.bool get success => $_getBF(0);
  @$pb.TagNumber(1)
  set success($core.bool value) => $_setBool(0, value);
  @$pb.TagNumber(1)
  $core.bool hasSuccess() => $_has(0);
  @$pb.TagNumber(1)
  void clearSuccess() => $_clearField(1);

  @$pb.TagNumber(2)
  $core.String get message => $_getSZ(1);
  @$pb.TagNumber(2)
  set message($core.String value) => $_setString(1, value);
  @$pb.TagNumber(2)
  $core.bool hasMessage() => $_has(1);
  @$pb.TagNumber(2)
  void clearMessage() => $_clearField(2);
}

class InstanceResponse extends $pb.GeneratedMessage {
  factory InstanceResponse({
    $core.bool? success,
    $core.String? message,
    Instance? instance,
  }) {
    final result = create();
    if (success != null) result.success = success;
    if (message != null) result.message = message;
    if (instance != null) result.instance = instance;
    return result;
  }

  InstanceResponse._();

  factory InstanceResponse.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory InstanceResponse.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'InstanceResponse', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOB(1, _omitFieldNames ? '' : 'success')
    ..aOS(2, _omitFieldNames ? '' : 'message')
    ..aOM<Instance>(3, _omitFieldNames ? '' : 'instance', subBuilder: Instance.create)
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  InstanceResponse clone() => InstanceResponse()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  InstanceResponse copyWith(void Function(InstanceResponse) updates) => super.copyWith((message) => updates(message as InstanceResponse)) as InstanceResponse;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static InstanceResponse create() => InstanceResponse._();
  @$core.override
  InstanceResponse createEmptyInstance() => create();
  static $pb.PbList<InstanceResponse> createRepeated() => $pb.PbList<InstanceResponse>();
  @$core.pragma('dart2js:noInline')
  static InstanceResponse getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<InstanceResponse>(create);
  static InstanceResponse? _defaultInstance;

  @$pb.TagNumber(1)
  $core.bool get success => $_getBF(0);
  @$pb.TagNumber(1)
  set success($core.bool value) => $_setBool(0, value);
  @$pb.TagNumber(1)
  $core.bool hasSuccess() => $_has(0);
  @$pb.TagNumber(1)
  void clearSuccess() => $_clearField(1);

  @$pb.TagNumber(2)
  $core.String get message => $_getSZ(1);
  @$pb.TagNumber(2)
  set message($core.String value) => $_setString(1, value);
  @$pb.TagNumber(2)
  $core.bool hasMessage() => $_has(1);
  @$pb.TagNumber(2)
  void clearMessage() => $_clearField(2);

  @$pb.TagNumber(3)
  Instance get instance => $_getN(2);
  @$pb.TagNumber(3)
  set instance(Instance value) => $_setField(3, value);
  @$pb.TagNumber(3)
  $core.bool hasInstance() => $_has(2);
  @$pb.TagNumber(3)
  void clearInstance() => $_clearField(3);
  @$pb.TagNumber(3)
  Instance ensureInstance() => $_ensure(2);
}

class Instance extends $pb.GeneratedMessage {
  factory Instance({
    $core.String? id,
    $core.String? name,
    $core.String? description,
    InstanceStatus? status,
    $core.String? userId,
    $core.String? projectId,
    $fixnum.Int64? createdAt,
    $fixnum.Int64? updatedAt,
    Molecule? molecule,
    $core.Iterable<Calculation>? calculations,
  }) {
    final result = create();
    if (id != null) result.id = id;
    if (name != null) result.name = name;
    if (description != null) result.description = description;
    if (status != null) result.status = status;
    if (userId != null) result.userId = userId;
    if (projectId != null) result.projectId = projectId;
    if (createdAt != null) result.createdAt = createdAt;
    if (updatedAt != null) result.updatedAt = updatedAt;
    if (molecule != null) result.molecule = molecule;
    if (calculations != null) result.calculations.addAll(calculations);
    return result;
  }

  Instance._();

  factory Instance.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory Instance.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'Instance', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOS(1, _omitFieldNames ? '' : 'id')
    ..aOS(2, _omitFieldNames ? '' : 'name')
    ..aOS(3, _omitFieldNames ? '' : 'description')
    ..e<InstanceStatus>(4, _omitFieldNames ? '' : 'status', $pb.PbFieldType.OE, defaultOrMaker: InstanceStatus.INSTANCE_STATUS_UNSPECIFIED, valueOf: InstanceStatus.valueOf, enumValues: InstanceStatus.values)
    ..aOS(5, _omitFieldNames ? '' : 'userId')
    ..aOS(6, _omitFieldNames ? '' : 'projectId')
    ..aInt64(7, _omitFieldNames ? '' : 'createdAt')
    ..aInt64(8, _omitFieldNames ? '' : 'updatedAt')
    ..aOM<Molecule>(9, _omitFieldNames ? '' : 'molecule', subBuilder: Molecule.create)
    ..pc<Calculation>(10, _omitFieldNames ? '' : 'calculations', $pb.PbFieldType.PM, subBuilder: Calculation.create)
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  Instance clone() => Instance()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  Instance copyWith(void Function(Instance) updates) => super.copyWith((message) => updates(message as Instance)) as Instance;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static Instance create() => Instance._();
  @$core.override
  Instance createEmptyInstance() => create();
  static $pb.PbList<Instance> createRepeated() => $pb.PbList<Instance>();
  @$core.pragma('dart2js:noInline')
  static Instance getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<Instance>(create);
  static Instance? _defaultInstance;

  @$pb.TagNumber(1)
  $core.String get id => $_getSZ(0);
  @$pb.TagNumber(1)
  set id($core.String value) => $_setString(0, value);
  @$pb.TagNumber(1)
  $core.bool hasId() => $_has(0);
  @$pb.TagNumber(1)
  void clearId() => $_clearField(1);

  @$pb.TagNumber(2)
  $core.String get name => $_getSZ(1);
  @$pb.TagNumber(2)
  set name($core.String value) => $_setString(1, value);
  @$pb.TagNumber(2)
  $core.bool hasName() => $_has(1);
  @$pb.TagNumber(2)
  void clearName() => $_clearField(2);

  @$pb.TagNumber(3)
  $core.String get description => $_getSZ(2);
  @$pb.TagNumber(3)
  set description($core.String value) => $_setString(2, value);
  @$pb.TagNumber(3)
  $core.bool hasDescription() => $_has(2);
  @$pb.TagNumber(3)
  void clearDescription() => $_clearField(3);

  @$pb.TagNumber(4)
  InstanceStatus get status => $_getN(3);
  @$pb.TagNumber(4)
  set status(InstanceStatus value) => $_setField(4, value);
  @$pb.TagNumber(4)
  $core.bool hasStatus() => $_has(3);
  @$pb.TagNumber(4)
  void clearStatus() => $_clearField(4);

  @$pb.TagNumber(5)
  $core.String get userId => $_getSZ(4);
  @$pb.TagNumber(5)
  set userId($core.String value) => $_setString(4, value);
  @$pb.TagNumber(5)
  $core.bool hasUserId() => $_has(4);
  @$pb.TagNumber(5)
  void clearUserId() => $_clearField(5);

  @$pb.TagNumber(6)
  $core.String get projectId => $_getSZ(5);
  @$pb.TagNumber(6)
  set projectId($core.String value) => $_setString(5, value);
  @$pb.TagNumber(6)
  $core.bool hasProjectId() => $_has(5);
  @$pb.TagNumber(6)
  void clearProjectId() => $_clearField(6);

  @$pb.TagNumber(7)
  $fixnum.Int64 get createdAt => $_getI64(6);
  @$pb.TagNumber(7)
  set createdAt($fixnum.Int64 value) => $_setInt64(6, value);
  @$pb.TagNumber(7)
  $core.bool hasCreatedAt() => $_has(6);
  @$pb.TagNumber(7)
  void clearCreatedAt() => $_clearField(7);

  @$pb.TagNumber(8)
  $fixnum.Int64 get updatedAt => $_getI64(7);
  @$pb.TagNumber(8)
  set updatedAt($fixnum.Int64 value) => $_setInt64(7, value);
  @$pb.TagNumber(8)
  $core.bool hasUpdatedAt() => $_has(7);
  @$pb.TagNumber(8)
  void clearUpdatedAt() => $_clearField(8);

  @$pb.TagNumber(9)
  Molecule get molecule => $_getN(8);
  @$pb.TagNumber(9)
  set molecule(Molecule value) => $_setField(9, value);
  @$pb.TagNumber(9)
  $core.bool hasMolecule() => $_has(8);
  @$pb.TagNumber(9)
  void clearMolecule() => $_clearField(9);
  @$pb.TagNumber(9)
  Molecule ensureMolecule() => $_ensure(8);

  @$pb.TagNumber(10)
  $pb.PbList<Calculation> get calculations => $_getList(9);
}

/// 計算関連メッセージ
class StartCalculationRequest extends $pb.GeneratedMessage {
  factory StartCalculationRequest({
    $core.String? instanceId,
    $core.String? method,
    $core.String? basisSet,
    $core.Iterable<$core.MapEntry<$core.String, $core.String>>? parameters,
    ConvergenceCriteria? convergenceCriteria,
    $core.int? maxIterations,
    $core.int? priority,
  }) {
    final result = create();
    if (instanceId != null) result.instanceId = instanceId;
    if (method != null) result.method = method;
    if (basisSet != null) result.basisSet = basisSet;
    if (parameters != null) result.parameters.addEntries(parameters);
    if (convergenceCriteria != null) result.convergenceCriteria = convergenceCriteria;
    if (maxIterations != null) result.maxIterations = maxIterations;
    if (priority != null) result.priority = priority;
    return result;
  }

  StartCalculationRequest._();

  factory StartCalculationRequest.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory StartCalculationRequest.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'StartCalculationRequest', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOS(1, _omitFieldNames ? '' : 'instanceId')
    ..aOS(2, _omitFieldNames ? '' : 'method')
    ..aOS(3, _omitFieldNames ? '' : 'basisSet')
    ..m<$core.String, $core.String>(4, _omitFieldNames ? '' : 'parameters', entryClassName: 'StartCalculationRequest.ParametersEntry', keyFieldType: $pb.PbFieldType.OS, valueFieldType: $pb.PbFieldType.OS, packageName: const $pb.PackageName('pyscf_front'))
    ..aOM<ConvergenceCriteria>(5, _omitFieldNames ? '' : 'convergenceCriteria', subBuilder: ConvergenceCriteria.create)
    ..a<$core.int>(6, _omitFieldNames ? '' : 'maxIterations', $pb.PbFieldType.O3)
    ..a<$core.int>(7, _omitFieldNames ? '' : 'priority', $pb.PbFieldType.O3)
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  StartCalculationRequest clone() => StartCalculationRequest()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  StartCalculationRequest copyWith(void Function(StartCalculationRequest) updates) => super.copyWith((message) => updates(message as StartCalculationRequest)) as StartCalculationRequest;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static StartCalculationRequest create() => StartCalculationRequest._();
  @$core.override
  StartCalculationRequest createEmptyInstance() => create();
  static $pb.PbList<StartCalculationRequest> createRepeated() => $pb.PbList<StartCalculationRequest>();
  @$core.pragma('dart2js:noInline')
  static StartCalculationRequest getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<StartCalculationRequest>(create);
  static StartCalculationRequest? _defaultInstance;

  @$pb.TagNumber(1)
  $core.String get instanceId => $_getSZ(0);
  @$pb.TagNumber(1)
  set instanceId($core.String value) => $_setString(0, value);
  @$pb.TagNumber(1)
  $core.bool hasInstanceId() => $_has(0);
  @$pb.TagNumber(1)
  void clearInstanceId() => $_clearField(1);

  @$pb.TagNumber(2)
  $core.String get method => $_getSZ(1);
  @$pb.TagNumber(2)
  set method($core.String value) => $_setString(1, value);
  @$pb.TagNumber(2)
  $core.bool hasMethod() => $_has(1);
  @$pb.TagNumber(2)
  void clearMethod() => $_clearField(2);

  @$pb.TagNumber(3)
  $core.String get basisSet => $_getSZ(2);
  @$pb.TagNumber(3)
  set basisSet($core.String value) => $_setString(2, value);
  @$pb.TagNumber(3)
  $core.bool hasBasisSet() => $_has(2);
  @$pb.TagNumber(3)
  void clearBasisSet() => $_clearField(3);

  @$pb.TagNumber(4)
  $pb.PbMap<$core.String, $core.String> get parameters => $_getMap(3);

  @$pb.TagNumber(5)
  ConvergenceCriteria get convergenceCriteria => $_getN(4);
  @$pb.TagNumber(5)
  set convergenceCriteria(ConvergenceCriteria value) => $_setField(5, value);
  @$pb.TagNumber(5)
  $core.bool hasConvergenceCriteria() => $_has(4);
  @$pb.TagNumber(5)
  void clearConvergenceCriteria() => $_clearField(5);
  @$pb.TagNumber(5)
  ConvergenceCriteria ensureConvergenceCriteria() => $_ensure(4);

  @$pb.TagNumber(6)
  $core.int get maxIterations => $_getIZ(5);
  @$pb.TagNumber(6)
  set maxIterations($core.int value) => $_setSignedInt32(5, value);
  @$pb.TagNumber(6)
  $core.bool hasMaxIterations() => $_has(5);
  @$pb.TagNumber(6)
  void clearMaxIterations() => $_clearField(6);

  @$pb.TagNumber(7)
  $core.int get priority => $_getIZ(6);
  @$pb.TagNumber(7)
  set priority($core.int value) => $_setSignedInt32(6, value);
  @$pb.TagNumber(7)
  $core.bool hasPriority() => $_has(6);
  @$pb.TagNumber(7)
  void clearPriority() => $_clearField(7);
}

class CancelCalculationRequest extends $pb.GeneratedMessage {
  factory CancelCalculationRequest({
    $core.String? calculationId,
  }) {
    final result = create();
    if (calculationId != null) result.calculationId = calculationId;
    return result;
  }

  CancelCalculationRequest._();

  factory CancelCalculationRequest.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory CancelCalculationRequest.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'CancelCalculationRequest', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOS(1, _omitFieldNames ? '' : 'calculationId')
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  CancelCalculationRequest clone() => CancelCalculationRequest()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  CancelCalculationRequest copyWith(void Function(CancelCalculationRequest) updates) => super.copyWith((message) => updates(message as CancelCalculationRequest)) as CancelCalculationRequest;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static CancelCalculationRequest create() => CancelCalculationRequest._();
  @$core.override
  CancelCalculationRequest createEmptyInstance() => create();
  static $pb.PbList<CancelCalculationRequest> createRepeated() => $pb.PbList<CancelCalculationRequest>();
  @$core.pragma('dart2js:noInline')
  static CancelCalculationRequest getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<CancelCalculationRequest>(create);
  static CancelCalculationRequest? _defaultInstance;

  @$pb.TagNumber(1)
  $core.String get calculationId => $_getSZ(0);
  @$pb.TagNumber(1)
  set calculationId($core.String value) => $_setString(0, value);
  @$pb.TagNumber(1)
  $core.bool hasCalculationId() => $_has(0);
  @$pb.TagNumber(1)
  void clearCalculationId() => $_clearField(1);
}

class CancelCalculationResponse extends $pb.GeneratedMessage {
  factory CancelCalculationResponse({
    $core.bool? success,
    $core.String? message,
  }) {
    final result = create();
    if (success != null) result.success = success;
    if (message != null) result.message = message;
    return result;
  }

  CancelCalculationResponse._();

  factory CancelCalculationResponse.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory CancelCalculationResponse.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'CancelCalculationResponse', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOB(1, _omitFieldNames ? '' : 'success')
    ..aOS(2, _omitFieldNames ? '' : 'message')
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  CancelCalculationResponse clone() => CancelCalculationResponse()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  CancelCalculationResponse copyWith(void Function(CancelCalculationResponse) updates) => super.copyWith((message) => updates(message as CancelCalculationResponse)) as CancelCalculationResponse;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static CancelCalculationResponse create() => CancelCalculationResponse._();
  @$core.override
  CancelCalculationResponse createEmptyInstance() => create();
  static $pb.PbList<CancelCalculationResponse> createRepeated() => $pb.PbList<CancelCalculationResponse>();
  @$core.pragma('dart2js:noInline')
  static CancelCalculationResponse getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<CancelCalculationResponse>(create);
  static CancelCalculationResponse? _defaultInstance;

  @$pb.TagNumber(1)
  $core.bool get success => $_getBF(0);
  @$pb.TagNumber(1)
  set success($core.bool value) => $_setBool(0, value);
  @$pb.TagNumber(1)
  $core.bool hasSuccess() => $_has(0);
  @$pb.TagNumber(1)
  void clearSuccess() => $_clearField(1);

  @$pb.TagNumber(2)
  $core.String get message => $_getSZ(1);
  @$pb.TagNumber(2)
  set message($core.String value) => $_setString(1, value);
  @$pb.TagNumber(2)
  $core.bool hasMessage() => $_has(1);
  @$pb.TagNumber(2)
  void clearMessage() => $_clearField(2);
}

class GetCalculationStatusRequest extends $pb.GeneratedMessage {
  factory GetCalculationStatusRequest({
    $core.String? calculationId,
  }) {
    final result = create();
    if (calculationId != null) result.calculationId = calculationId;
    return result;
  }

  GetCalculationStatusRequest._();

  factory GetCalculationStatusRequest.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory GetCalculationStatusRequest.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'GetCalculationStatusRequest', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOS(1, _omitFieldNames ? '' : 'calculationId')
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  GetCalculationStatusRequest clone() => GetCalculationStatusRequest()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  GetCalculationStatusRequest copyWith(void Function(GetCalculationStatusRequest) updates) => super.copyWith((message) => updates(message as GetCalculationStatusRequest)) as GetCalculationStatusRequest;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static GetCalculationStatusRequest create() => GetCalculationStatusRequest._();
  @$core.override
  GetCalculationStatusRequest createEmptyInstance() => create();
  static $pb.PbList<GetCalculationStatusRequest> createRepeated() => $pb.PbList<GetCalculationStatusRequest>();
  @$core.pragma('dart2js:noInline')
  static GetCalculationStatusRequest getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<GetCalculationStatusRequest>(create);
  static GetCalculationStatusRequest? _defaultInstance;

  @$pb.TagNumber(1)
  $core.String get calculationId => $_getSZ(0);
  @$pb.TagNumber(1)
  set calculationId($core.String value) => $_setString(0, value);
  @$pb.TagNumber(1)
  $core.bool hasCalculationId() => $_has(0);
  @$pb.TagNumber(1)
  void clearCalculationId() => $_clearField(1);
}

class CalculationStatusResponse extends $pb.GeneratedMessage {
  factory CalculationStatusResponse({
    $core.bool? success,
    $core.String? message,
    CalculationStatus? status,
  }) {
    final result = create();
    if (success != null) result.success = success;
    if (message != null) result.message = message;
    if (status != null) result.status = status;
    return result;
  }

  CalculationStatusResponse._();

  factory CalculationStatusResponse.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory CalculationStatusResponse.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'CalculationStatusResponse', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOB(1, _omitFieldNames ? '' : 'success')
    ..aOS(2, _omitFieldNames ? '' : 'message')
    ..e<CalculationStatus>(3, _omitFieldNames ? '' : 'status', $pb.PbFieldType.OE, defaultOrMaker: CalculationStatus.CALCULATION_STATUS_UNSPECIFIED, valueOf: CalculationStatus.valueOf, enumValues: CalculationStatus.values)
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  CalculationStatusResponse clone() => CalculationStatusResponse()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  CalculationStatusResponse copyWith(void Function(CalculationStatusResponse) updates) => super.copyWith((message) => updates(message as CalculationStatusResponse)) as CalculationStatusResponse;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static CalculationStatusResponse create() => CalculationStatusResponse._();
  @$core.override
  CalculationStatusResponse createEmptyInstance() => create();
  static $pb.PbList<CalculationStatusResponse> createRepeated() => $pb.PbList<CalculationStatusResponse>();
  @$core.pragma('dart2js:noInline')
  static CalculationStatusResponse getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<CalculationStatusResponse>(create);
  static CalculationStatusResponse? _defaultInstance;

  @$pb.TagNumber(1)
  $core.bool get success => $_getBF(0);
  @$pb.TagNumber(1)
  set success($core.bool value) => $_setBool(0, value);
  @$pb.TagNumber(1)
  $core.bool hasSuccess() => $_has(0);
  @$pb.TagNumber(1)
  void clearSuccess() => $_clearField(1);

  @$pb.TagNumber(2)
  $core.String get message => $_getSZ(1);
  @$pb.TagNumber(2)
  set message($core.String value) => $_setString(1, value);
  @$pb.TagNumber(2)
  $core.bool hasMessage() => $_has(1);
  @$pb.TagNumber(2)
  void clearMessage() => $_clearField(2);

  @$pb.TagNumber(3)
  CalculationStatus get status => $_getN(2);
  @$pb.TagNumber(3)
  set status(CalculationStatus value) => $_setField(3, value);
  @$pb.TagNumber(3)
  $core.bool hasStatus() => $_has(2);
  @$pb.TagNumber(3)
  void clearStatus() => $_clearField(3);
}

class CalculationProgress extends $pb.GeneratedMessage {
  factory CalculationProgress({
    $core.String? calculationId,
    $core.String? jobId,
    CalculationStatus? status,
    $core.double? progressPercentage,
    $core.String? currentStep,
    $core.String? message,
    $fixnum.Int64? timestamp,
    $core.String? estimatedTimeRemaining,
  }) {
    final result = create();
    if (calculationId != null) result.calculationId = calculationId;
    if (jobId != null) result.jobId = jobId;
    if (status != null) result.status = status;
    if (progressPercentage != null) result.progressPercentage = progressPercentage;
    if (currentStep != null) result.currentStep = currentStep;
    if (message != null) result.message = message;
    if (timestamp != null) result.timestamp = timestamp;
    if (estimatedTimeRemaining != null) result.estimatedTimeRemaining = estimatedTimeRemaining;
    return result;
  }

  CalculationProgress._();

  factory CalculationProgress.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory CalculationProgress.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'CalculationProgress', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOS(1, _omitFieldNames ? '' : 'calculationId')
    ..aOS(2, _omitFieldNames ? '' : 'jobId')
    ..e<CalculationStatus>(3, _omitFieldNames ? '' : 'status', $pb.PbFieldType.OE, defaultOrMaker: CalculationStatus.CALCULATION_STATUS_UNSPECIFIED, valueOf: CalculationStatus.valueOf, enumValues: CalculationStatus.values)
    ..a<$core.double>(4, _omitFieldNames ? '' : 'progressPercentage', $pb.PbFieldType.OD)
    ..aOS(5, _omitFieldNames ? '' : 'currentStep')
    ..aOS(6, _omitFieldNames ? '' : 'message')
    ..aInt64(7, _omitFieldNames ? '' : 'timestamp')
    ..aOS(8, _omitFieldNames ? '' : 'estimatedTimeRemaining')
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  CalculationProgress clone() => CalculationProgress()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  CalculationProgress copyWith(void Function(CalculationProgress) updates) => super.copyWith((message) => updates(message as CalculationProgress)) as CalculationProgress;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static CalculationProgress create() => CalculationProgress._();
  @$core.override
  CalculationProgress createEmptyInstance() => create();
  static $pb.PbList<CalculationProgress> createRepeated() => $pb.PbList<CalculationProgress>();
  @$core.pragma('dart2js:noInline')
  static CalculationProgress getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<CalculationProgress>(create);
  static CalculationProgress? _defaultInstance;

  @$pb.TagNumber(1)
  $core.String get calculationId => $_getSZ(0);
  @$pb.TagNumber(1)
  set calculationId($core.String value) => $_setString(0, value);
  @$pb.TagNumber(1)
  $core.bool hasCalculationId() => $_has(0);
  @$pb.TagNumber(1)
  void clearCalculationId() => $_clearField(1);

  @$pb.TagNumber(2)
  $core.String get jobId => $_getSZ(1);
  @$pb.TagNumber(2)
  set jobId($core.String value) => $_setString(1, value);
  @$pb.TagNumber(2)
  $core.bool hasJobId() => $_has(1);
  @$pb.TagNumber(2)
  void clearJobId() => $_clearField(2);

  @$pb.TagNumber(3)
  CalculationStatus get status => $_getN(2);
  @$pb.TagNumber(3)
  set status(CalculationStatus value) => $_setField(3, value);
  @$pb.TagNumber(3)
  $core.bool hasStatus() => $_has(2);
  @$pb.TagNumber(3)
  void clearStatus() => $_clearField(3);

  @$pb.TagNumber(4)
  $core.double get progressPercentage => $_getN(3);
  @$pb.TagNumber(4)
  set progressPercentage($core.double value) => $_setDouble(3, value);
  @$pb.TagNumber(4)
  $core.bool hasProgressPercentage() => $_has(3);
  @$pb.TagNumber(4)
  void clearProgressPercentage() => $_clearField(4);

  @$pb.TagNumber(5)
  $core.String get currentStep => $_getSZ(4);
  @$pb.TagNumber(5)
  set currentStep($core.String value) => $_setString(4, value);
  @$pb.TagNumber(5)
  $core.bool hasCurrentStep() => $_has(4);
  @$pb.TagNumber(5)
  void clearCurrentStep() => $_clearField(5);

  @$pb.TagNumber(6)
  $core.String get message => $_getSZ(5);
  @$pb.TagNumber(6)
  set message($core.String value) => $_setString(5, value);
  @$pb.TagNumber(6)
  $core.bool hasMessage() => $_has(5);
  @$pb.TagNumber(6)
  void clearMessage() => $_clearField(6);

  @$pb.TagNumber(7)
  $fixnum.Int64 get timestamp => $_getI64(6);
  @$pb.TagNumber(7)
  set timestamp($fixnum.Int64 value) => $_setInt64(6, value);
  @$pb.TagNumber(7)
  $core.bool hasTimestamp() => $_has(6);
  @$pb.TagNumber(7)
  void clearTimestamp() => $_clearField(7);

  @$pb.TagNumber(8)
  $core.String get estimatedTimeRemaining => $_getSZ(7);
  @$pb.TagNumber(8)
  set estimatedTimeRemaining($core.String value) => $_setString(7, value);
  @$pb.TagNumber(8)
  $core.bool hasEstimatedTimeRemaining() => $_has(7);
  @$pb.TagNumber(8)
  void clearEstimatedTimeRemaining() => $_clearField(8);
}

class Calculation extends $pb.GeneratedMessage {
  factory Calculation({
    $core.String? id,
    $core.String? instanceId,
    $core.String? method,
    $core.String? basisSet,
    $core.Iterable<$core.MapEntry<$core.String, $core.String>>? parameters,
    ConvergenceCriteria? convergenceCriteria,
    $core.int? maxIterations,
    $fixnum.Int64? startTime,
    $fixnum.Int64? endTime,
    CalculationStatus? status,
    $core.String? errorMessage,
    $fixnum.Int64? createdAt,
  }) {
    final result = create();
    if (id != null) result.id = id;
    if (instanceId != null) result.instanceId = instanceId;
    if (method != null) result.method = method;
    if (basisSet != null) result.basisSet = basisSet;
    if (parameters != null) result.parameters.addEntries(parameters);
    if (convergenceCriteria != null) result.convergenceCriteria = convergenceCriteria;
    if (maxIterations != null) result.maxIterations = maxIterations;
    if (startTime != null) result.startTime = startTime;
    if (endTime != null) result.endTime = endTime;
    if (status != null) result.status = status;
    if (errorMessage != null) result.errorMessage = errorMessage;
    if (createdAt != null) result.createdAt = createdAt;
    return result;
  }

  Calculation._();

  factory Calculation.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory Calculation.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'Calculation', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOS(1, _omitFieldNames ? '' : 'id')
    ..aOS(2, _omitFieldNames ? '' : 'instanceId')
    ..aOS(3, _omitFieldNames ? '' : 'method')
    ..aOS(4, _omitFieldNames ? '' : 'basisSet')
    ..m<$core.String, $core.String>(5, _omitFieldNames ? '' : 'parameters', entryClassName: 'Calculation.ParametersEntry', keyFieldType: $pb.PbFieldType.OS, valueFieldType: $pb.PbFieldType.OS, packageName: const $pb.PackageName('pyscf_front'))
    ..aOM<ConvergenceCriteria>(6, _omitFieldNames ? '' : 'convergenceCriteria', subBuilder: ConvergenceCriteria.create)
    ..a<$core.int>(7, _omitFieldNames ? '' : 'maxIterations', $pb.PbFieldType.O3)
    ..aInt64(8, _omitFieldNames ? '' : 'startTime')
    ..aInt64(9, _omitFieldNames ? '' : 'endTime')
    ..e<CalculationStatus>(10, _omitFieldNames ? '' : 'status', $pb.PbFieldType.OE, defaultOrMaker: CalculationStatus.CALCULATION_STATUS_UNSPECIFIED, valueOf: CalculationStatus.valueOf, enumValues: CalculationStatus.values)
    ..aOS(11, _omitFieldNames ? '' : 'errorMessage')
    ..aInt64(12, _omitFieldNames ? '' : 'createdAt')
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  Calculation clone() => Calculation()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  Calculation copyWith(void Function(Calculation) updates) => super.copyWith((message) => updates(message as Calculation)) as Calculation;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static Calculation create() => Calculation._();
  @$core.override
  Calculation createEmptyInstance() => create();
  static $pb.PbList<Calculation> createRepeated() => $pb.PbList<Calculation>();
  @$core.pragma('dart2js:noInline')
  static Calculation getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<Calculation>(create);
  static Calculation? _defaultInstance;

  @$pb.TagNumber(1)
  $core.String get id => $_getSZ(0);
  @$pb.TagNumber(1)
  set id($core.String value) => $_setString(0, value);
  @$pb.TagNumber(1)
  $core.bool hasId() => $_has(0);
  @$pb.TagNumber(1)
  void clearId() => $_clearField(1);

  @$pb.TagNumber(2)
  $core.String get instanceId => $_getSZ(1);
  @$pb.TagNumber(2)
  set instanceId($core.String value) => $_setString(1, value);
  @$pb.TagNumber(2)
  $core.bool hasInstanceId() => $_has(1);
  @$pb.TagNumber(2)
  void clearInstanceId() => $_clearField(2);

  @$pb.TagNumber(3)
  $core.String get method => $_getSZ(2);
  @$pb.TagNumber(3)
  set method($core.String value) => $_setString(2, value);
  @$pb.TagNumber(3)
  $core.bool hasMethod() => $_has(2);
  @$pb.TagNumber(3)
  void clearMethod() => $_clearField(3);

  @$pb.TagNumber(4)
  $core.String get basisSet => $_getSZ(3);
  @$pb.TagNumber(4)
  set basisSet($core.String value) => $_setString(3, value);
  @$pb.TagNumber(4)
  $core.bool hasBasisSet() => $_has(3);
  @$pb.TagNumber(4)
  void clearBasisSet() => $_clearField(4);

  @$pb.TagNumber(5)
  $pb.PbMap<$core.String, $core.String> get parameters => $_getMap(4);

  @$pb.TagNumber(6)
  ConvergenceCriteria get convergenceCriteria => $_getN(5);
  @$pb.TagNumber(6)
  set convergenceCriteria(ConvergenceCriteria value) => $_setField(6, value);
  @$pb.TagNumber(6)
  $core.bool hasConvergenceCriteria() => $_has(5);
  @$pb.TagNumber(6)
  void clearConvergenceCriteria() => $_clearField(6);
  @$pb.TagNumber(6)
  ConvergenceCriteria ensureConvergenceCriteria() => $_ensure(5);

  @$pb.TagNumber(7)
  $core.int get maxIterations => $_getIZ(6);
  @$pb.TagNumber(7)
  set maxIterations($core.int value) => $_setSignedInt32(6, value);
  @$pb.TagNumber(7)
  $core.bool hasMaxIterations() => $_has(6);
  @$pb.TagNumber(7)
  void clearMaxIterations() => $_clearField(7);

  @$pb.TagNumber(8)
  $fixnum.Int64 get startTime => $_getI64(7);
  @$pb.TagNumber(8)
  set startTime($fixnum.Int64 value) => $_setInt64(7, value);
  @$pb.TagNumber(8)
  $core.bool hasStartTime() => $_has(7);
  @$pb.TagNumber(8)
  void clearStartTime() => $_clearField(8);

  @$pb.TagNumber(9)
  $fixnum.Int64 get endTime => $_getI64(8);
  @$pb.TagNumber(9)
  set endTime($fixnum.Int64 value) => $_setInt64(8, value);
  @$pb.TagNumber(9)
  $core.bool hasEndTime() => $_has(8);
  @$pb.TagNumber(9)
  void clearEndTime() => $_clearField(9);

  @$pb.TagNumber(10)
  CalculationStatus get status => $_getN(9);
  @$pb.TagNumber(10)
  set status(CalculationStatus value) => $_setField(10, value);
  @$pb.TagNumber(10)
  $core.bool hasStatus() => $_has(9);
  @$pb.TagNumber(10)
  void clearStatus() => $_clearField(10);

  @$pb.TagNumber(11)
  $core.String get errorMessage => $_getSZ(10);
  @$pb.TagNumber(11)
  set errorMessage($core.String value) => $_setString(10, value);
  @$pb.TagNumber(11)
  $core.bool hasErrorMessage() => $_has(10);
  @$pb.TagNumber(11)
  void clearErrorMessage() => $_clearField(11);

  @$pb.TagNumber(12)
  $fixnum.Int64 get createdAt => $_getI64(11);
  @$pb.TagNumber(12)
  set createdAt($fixnum.Int64 value) => $_setInt64(11, value);
  @$pb.TagNumber(12)
  $core.bool hasCreatedAt() => $_has(11);
  @$pb.TagNumber(12)
  void clearCreatedAt() => $_clearField(12);
}

class ConvergenceCriteria extends $pb.GeneratedMessage {
  factory ConvergenceCriteria({
    $core.double? energyThreshold,
    $core.double? densityThreshold,
    $core.double? gradientThreshold,
  }) {
    final result = create();
    if (energyThreshold != null) result.energyThreshold = energyThreshold;
    if (densityThreshold != null) result.densityThreshold = densityThreshold;
    if (gradientThreshold != null) result.gradientThreshold = gradientThreshold;
    return result;
  }

  ConvergenceCriteria._();

  factory ConvergenceCriteria.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory ConvergenceCriteria.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'ConvergenceCriteria', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..a<$core.double>(1, _omitFieldNames ? '' : 'energyThreshold', $pb.PbFieldType.OD)
    ..a<$core.double>(2, _omitFieldNames ? '' : 'densityThreshold', $pb.PbFieldType.OD)
    ..a<$core.double>(3, _omitFieldNames ? '' : 'gradientThreshold', $pb.PbFieldType.OD)
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  ConvergenceCriteria clone() => ConvergenceCriteria()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  ConvergenceCriteria copyWith(void Function(ConvergenceCriteria) updates) => super.copyWith((message) => updates(message as ConvergenceCriteria)) as ConvergenceCriteria;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static ConvergenceCriteria create() => ConvergenceCriteria._();
  @$core.override
  ConvergenceCriteria createEmptyInstance() => create();
  static $pb.PbList<ConvergenceCriteria> createRepeated() => $pb.PbList<ConvergenceCriteria>();
  @$core.pragma('dart2js:noInline')
  static ConvergenceCriteria getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<ConvergenceCriteria>(create);
  static ConvergenceCriteria? _defaultInstance;

  @$pb.TagNumber(1)
  $core.double get energyThreshold => $_getN(0);
  @$pb.TagNumber(1)
  set energyThreshold($core.double value) => $_setDouble(0, value);
  @$pb.TagNumber(1)
  $core.bool hasEnergyThreshold() => $_has(0);
  @$pb.TagNumber(1)
  void clearEnergyThreshold() => $_clearField(1);

  @$pb.TagNumber(2)
  $core.double get densityThreshold => $_getN(1);
  @$pb.TagNumber(2)
  set densityThreshold($core.double value) => $_setDouble(1, value);
  @$pb.TagNumber(2)
  $core.bool hasDensityThreshold() => $_has(1);
  @$pb.TagNumber(2)
  void clearDensityThreshold() => $_clearField(2);

  @$pb.TagNumber(3)
  $core.double get gradientThreshold => $_getN(2);
  @$pb.TagNumber(3)
  set gradientThreshold($core.double value) => $_setDouble(2, value);
  @$pb.TagNumber(3)
  $core.bool hasGradientThreshold() => $_has(2);
  @$pb.TagNumber(3)
  void clearGradientThreshold() => $_clearField(3);
}

/// 結果関連メッセージ
class GetResultsRequest extends $pb.GeneratedMessage {
  factory GetResultsRequest({
    $core.String? calculationId,
    $core.Iterable<$core.String>? resultTypes,
  }) {
    final result = create();
    if (calculationId != null) result.calculationId = calculationId;
    if (resultTypes != null) result.resultTypes.addAll(resultTypes);
    return result;
  }

  GetResultsRequest._();

  factory GetResultsRequest.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory GetResultsRequest.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'GetResultsRequest', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOS(1, _omitFieldNames ? '' : 'calculationId')
    ..pPS(2, _omitFieldNames ? '' : 'resultTypes')
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  GetResultsRequest clone() => GetResultsRequest()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  GetResultsRequest copyWith(void Function(GetResultsRequest) updates) => super.copyWith((message) => updates(message as GetResultsRequest)) as GetResultsRequest;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static GetResultsRequest create() => GetResultsRequest._();
  @$core.override
  GetResultsRequest createEmptyInstance() => create();
  static $pb.PbList<GetResultsRequest> createRepeated() => $pb.PbList<GetResultsRequest>();
  @$core.pragma('dart2js:noInline')
  static GetResultsRequest getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<GetResultsRequest>(create);
  static GetResultsRequest? _defaultInstance;

  @$pb.TagNumber(1)
  $core.String get calculationId => $_getSZ(0);
  @$pb.TagNumber(1)
  set calculationId($core.String value) => $_setString(0, value);
  @$pb.TagNumber(1)
  $core.bool hasCalculationId() => $_has(0);
  @$pb.TagNumber(1)
  void clearCalculationId() => $_clearField(1);

  @$pb.TagNumber(2)
  $pb.PbList<$core.String> get resultTypes => $_getList(1);
}

class ExportResultsRequest extends $pb.GeneratedMessage {
  factory ExportResultsRequest({
    $core.String? calculationId,
    ExportFormat? format,
    $core.Iterable<$core.String>? resultTypes,
  }) {
    final result = create();
    if (calculationId != null) result.calculationId = calculationId;
    if (format != null) result.format = format;
    if (resultTypes != null) result.resultTypes.addAll(resultTypes);
    return result;
  }

  ExportResultsRequest._();

  factory ExportResultsRequest.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory ExportResultsRequest.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'ExportResultsRequest', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOS(1, _omitFieldNames ? '' : 'calculationId')
    ..e<ExportFormat>(2, _omitFieldNames ? '' : 'format', $pb.PbFieldType.OE, defaultOrMaker: ExportFormat.EXPORT_FORMAT_UNSPECIFIED, valueOf: ExportFormat.valueOf, enumValues: ExportFormat.values)
    ..pPS(3, _omitFieldNames ? '' : 'resultTypes')
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  ExportResultsRequest clone() => ExportResultsRequest()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  ExportResultsRequest copyWith(void Function(ExportResultsRequest) updates) => super.copyWith((message) => updates(message as ExportResultsRequest)) as ExportResultsRequest;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static ExportResultsRequest create() => ExportResultsRequest._();
  @$core.override
  ExportResultsRequest createEmptyInstance() => create();
  static $pb.PbList<ExportResultsRequest> createRepeated() => $pb.PbList<ExportResultsRequest>();
  @$core.pragma('dart2js:noInline')
  static ExportResultsRequest getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<ExportResultsRequest>(create);
  static ExportResultsRequest? _defaultInstance;

  @$pb.TagNumber(1)
  $core.String get calculationId => $_getSZ(0);
  @$pb.TagNumber(1)
  set calculationId($core.String value) => $_setString(0, value);
  @$pb.TagNumber(1)
  $core.bool hasCalculationId() => $_has(0);
  @$pb.TagNumber(1)
  void clearCalculationId() => $_clearField(1);

  @$pb.TagNumber(2)
  ExportFormat get format => $_getN(1);
  @$pb.TagNumber(2)
  set format(ExportFormat value) => $_setField(2, value);
  @$pb.TagNumber(2)
  $core.bool hasFormat() => $_has(1);
  @$pb.TagNumber(2)
  void clearFormat() => $_clearField(2);

  @$pb.TagNumber(3)
  $pb.PbList<$core.String> get resultTypes => $_getList(2);
}

class ExportResultsResponse extends $pb.GeneratedMessage {
  factory ExportResultsResponse({
    $core.bool? success,
    $core.String? message,
    $core.String? filePath,
    $core.List<$core.int>? fileData,
  }) {
    final result = create();
    if (success != null) result.success = success;
    if (message != null) result.message = message;
    if (filePath != null) result.filePath = filePath;
    if (fileData != null) result.fileData = fileData;
    return result;
  }

  ExportResultsResponse._();

  factory ExportResultsResponse.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory ExportResultsResponse.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'ExportResultsResponse', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOB(1, _omitFieldNames ? '' : 'success')
    ..aOS(2, _omitFieldNames ? '' : 'message')
    ..aOS(3, _omitFieldNames ? '' : 'filePath')
    ..a<$core.List<$core.int>>(4, _omitFieldNames ? '' : 'fileData', $pb.PbFieldType.OY)
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  ExportResultsResponse clone() => ExportResultsResponse()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  ExportResultsResponse copyWith(void Function(ExportResultsResponse) updates) => super.copyWith((message) => updates(message as ExportResultsResponse)) as ExportResultsResponse;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static ExportResultsResponse create() => ExportResultsResponse._();
  @$core.override
  ExportResultsResponse createEmptyInstance() => create();
  static $pb.PbList<ExportResultsResponse> createRepeated() => $pb.PbList<ExportResultsResponse>();
  @$core.pragma('dart2js:noInline')
  static ExportResultsResponse getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<ExportResultsResponse>(create);
  static ExportResultsResponse? _defaultInstance;

  @$pb.TagNumber(1)
  $core.bool get success => $_getBF(0);
  @$pb.TagNumber(1)
  set success($core.bool value) => $_setBool(0, value);
  @$pb.TagNumber(1)
  $core.bool hasSuccess() => $_has(0);
  @$pb.TagNumber(1)
  void clearSuccess() => $_clearField(1);

  @$pb.TagNumber(2)
  $core.String get message => $_getSZ(1);
  @$pb.TagNumber(2)
  set message($core.String value) => $_setString(1, value);
  @$pb.TagNumber(2)
  $core.bool hasMessage() => $_has(1);
  @$pb.TagNumber(2)
  void clearMessage() => $_clearField(2);

  @$pb.TagNumber(3)
  $core.String get filePath => $_getSZ(2);
  @$pb.TagNumber(3)
  set filePath($core.String value) => $_setString(2, value);
  @$pb.TagNumber(3)
  $core.bool hasFilePath() => $_has(2);
  @$pb.TagNumber(3)
  void clearFilePath() => $_clearField(3);

  @$pb.TagNumber(4)
  $core.List<$core.int> get fileData => $_getN(3);
  @$pb.TagNumber(4)
  set fileData($core.List<$core.int> value) => $_setBytes(3, value);
  @$pb.TagNumber(4)
  $core.bool hasFileData() => $_has(3);
  @$pb.TagNumber(4)
  void clearFileData() => $_clearField(4);
}

class ResultsResponse extends $pb.GeneratedMessage {
  factory ResultsResponse({
    $core.bool? success,
    $core.String? message,
    CalculationResults? results,
  }) {
    final result = create();
    if (success != null) result.success = success;
    if (message != null) result.message = message;
    if (results != null) result.results = results;
    return result;
  }

  ResultsResponse._();

  factory ResultsResponse.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory ResultsResponse.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'ResultsResponse', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOB(1, _omitFieldNames ? '' : 'success')
    ..aOS(2, _omitFieldNames ? '' : 'message')
    ..aOM<CalculationResults>(3, _omitFieldNames ? '' : 'results', subBuilder: CalculationResults.create)
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  ResultsResponse clone() => ResultsResponse()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  ResultsResponse copyWith(void Function(ResultsResponse) updates) => super.copyWith((message) => updates(message as ResultsResponse)) as ResultsResponse;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static ResultsResponse create() => ResultsResponse._();
  @$core.override
  ResultsResponse createEmptyInstance() => create();
  static $pb.PbList<ResultsResponse> createRepeated() => $pb.PbList<ResultsResponse>();
  @$core.pragma('dart2js:noInline')
  static ResultsResponse getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<ResultsResponse>(create);
  static ResultsResponse? _defaultInstance;

  @$pb.TagNumber(1)
  $core.bool get success => $_getBF(0);
  @$pb.TagNumber(1)
  set success($core.bool value) => $_setBool(0, value);
  @$pb.TagNumber(1)
  $core.bool hasSuccess() => $_has(0);
  @$pb.TagNumber(1)
  void clearSuccess() => $_clearField(1);

  @$pb.TagNumber(2)
  $core.String get message => $_getSZ(1);
  @$pb.TagNumber(2)
  set message($core.String value) => $_setString(1, value);
  @$pb.TagNumber(2)
  $core.bool hasMessage() => $_has(1);
  @$pb.TagNumber(2)
  void clearMessage() => $_clearField(2);

  @$pb.TagNumber(3)
  CalculationResults get results => $_getN(2);
  @$pb.TagNumber(3)
  set results(CalculationResults value) => $_setField(3, value);
  @$pb.TagNumber(3)
  $core.bool hasResults() => $_has(2);
  @$pb.TagNumber(3)
  void clearResults() => $_clearField(3);
  @$pb.TagNumber(3)
  CalculationResults ensureResults() => $_ensure(2);
}

class CalculationResults extends $pb.GeneratedMessage {
  factory CalculationResults({
    $core.String? calculationId,
    $core.bool? converged,
    $core.double? totalEnergy,
    $core.String? energyUnit,
    $core.int? calculationTimeSeconds,
    EnergyResults? energyResults,
    OrbitalResults? orbitalResults,
    FrequencyResults? frequencyResults,
    PropertyResults? propertyResults,
  }) {
    final result = create();
    if (calculationId != null) result.calculationId = calculationId;
    if (converged != null) result.converged = converged;
    if (totalEnergy != null) result.totalEnergy = totalEnergy;
    if (energyUnit != null) result.energyUnit = energyUnit;
    if (calculationTimeSeconds != null) result.calculationTimeSeconds = calculationTimeSeconds;
    if (energyResults != null) result.energyResults = energyResults;
    if (orbitalResults != null) result.orbitalResults = orbitalResults;
    if (frequencyResults != null) result.frequencyResults = frequencyResults;
    if (propertyResults != null) result.propertyResults = propertyResults;
    return result;
  }

  CalculationResults._();

  factory CalculationResults.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory CalculationResults.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'CalculationResults', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOS(1, _omitFieldNames ? '' : 'calculationId')
    ..aOB(2, _omitFieldNames ? '' : 'converged')
    ..a<$core.double>(3, _omitFieldNames ? '' : 'totalEnergy', $pb.PbFieldType.OD)
    ..aOS(4, _omitFieldNames ? '' : 'energyUnit')
    ..a<$core.int>(5, _omitFieldNames ? '' : 'calculationTimeSeconds', $pb.PbFieldType.O3)
    ..aOM<EnergyResults>(6, _omitFieldNames ? '' : 'energyResults', subBuilder: EnergyResults.create)
    ..aOM<OrbitalResults>(7, _omitFieldNames ? '' : 'orbitalResults', subBuilder: OrbitalResults.create)
    ..aOM<FrequencyResults>(8, _omitFieldNames ? '' : 'frequencyResults', subBuilder: FrequencyResults.create)
    ..aOM<PropertyResults>(9, _omitFieldNames ? '' : 'propertyResults', subBuilder: PropertyResults.create)
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  CalculationResults clone() => CalculationResults()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  CalculationResults copyWith(void Function(CalculationResults) updates) => super.copyWith((message) => updates(message as CalculationResults)) as CalculationResults;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static CalculationResults create() => CalculationResults._();
  @$core.override
  CalculationResults createEmptyInstance() => create();
  static $pb.PbList<CalculationResults> createRepeated() => $pb.PbList<CalculationResults>();
  @$core.pragma('dart2js:noInline')
  static CalculationResults getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<CalculationResults>(create);
  static CalculationResults? _defaultInstance;

  @$pb.TagNumber(1)
  $core.String get calculationId => $_getSZ(0);
  @$pb.TagNumber(1)
  set calculationId($core.String value) => $_setString(0, value);
  @$pb.TagNumber(1)
  $core.bool hasCalculationId() => $_has(0);
  @$pb.TagNumber(1)
  void clearCalculationId() => $_clearField(1);

  @$pb.TagNumber(2)
  $core.bool get converged => $_getBF(1);
  @$pb.TagNumber(2)
  set converged($core.bool value) => $_setBool(1, value);
  @$pb.TagNumber(2)
  $core.bool hasConverged() => $_has(1);
  @$pb.TagNumber(2)
  void clearConverged() => $_clearField(2);

  @$pb.TagNumber(3)
  $core.double get totalEnergy => $_getN(2);
  @$pb.TagNumber(3)
  set totalEnergy($core.double value) => $_setDouble(2, value);
  @$pb.TagNumber(3)
  $core.bool hasTotalEnergy() => $_has(2);
  @$pb.TagNumber(3)
  void clearTotalEnergy() => $_clearField(3);

  @$pb.TagNumber(4)
  $core.String get energyUnit => $_getSZ(3);
  @$pb.TagNumber(4)
  set energyUnit($core.String value) => $_setString(3, value);
  @$pb.TagNumber(4)
  $core.bool hasEnergyUnit() => $_has(3);
  @$pb.TagNumber(4)
  void clearEnergyUnit() => $_clearField(4);

  @$pb.TagNumber(5)
  $core.int get calculationTimeSeconds => $_getIZ(4);
  @$pb.TagNumber(5)
  set calculationTimeSeconds($core.int value) => $_setSignedInt32(4, value);
  @$pb.TagNumber(5)
  $core.bool hasCalculationTimeSeconds() => $_has(4);
  @$pb.TagNumber(5)
  void clearCalculationTimeSeconds() => $_clearField(5);

  @$pb.TagNumber(6)
  EnergyResults get energyResults => $_getN(5);
  @$pb.TagNumber(6)
  set energyResults(EnergyResults value) => $_setField(6, value);
  @$pb.TagNumber(6)
  $core.bool hasEnergyResults() => $_has(5);
  @$pb.TagNumber(6)
  void clearEnergyResults() => $_clearField(6);
  @$pb.TagNumber(6)
  EnergyResults ensureEnergyResults() => $_ensure(5);

  @$pb.TagNumber(7)
  OrbitalResults get orbitalResults => $_getN(6);
  @$pb.TagNumber(7)
  set orbitalResults(OrbitalResults value) => $_setField(7, value);
  @$pb.TagNumber(7)
  $core.bool hasOrbitalResults() => $_has(6);
  @$pb.TagNumber(7)
  void clearOrbitalResults() => $_clearField(7);
  @$pb.TagNumber(7)
  OrbitalResults ensureOrbitalResults() => $_ensure(6);

  @$pb.TagNumber(8)
  FrequencyResults get frequencyResults => $_getN(7);
  @$pb.TagNumber(8)
  set frequencyResults(FrequencyResults value) => $_setField(8, value);
  @$pb.TagNumber(8)
  $core.bool hasFrequencyResults() => $_has(7);
  @$pb.TagNumber(8)
  void clearFrequencyResults() => $_clearField(8);
  @$pb.TagNumber(8)
  FrequencyResults ensureFrequencyResults() => $_ensure(7);

  @$pb.TagNumber(9)
  PropertyResults get propertyResults => $_getN(8);
  @$pb.TagNumber(9)
  set propertyResults(PropertyResults value) => $_setField(9, value);
  @$pb.TagNumber(9)
  $core.bool hasPropertyResults() => $_has(8);
  @$pb.TagNumber(9)
  void clearPropertyResults() => $_clearField(9);
  @$pb.TagNumber(9)
  PropertyResults ensurePropertyResults() => $_ensure(8);
}

class EnergyResults extends $pb.GeneratedMessage {
  factory EnergyResults({
    $core.double? totalEnergy,
    $core.double? nuclearRepulsionEnergy,
    $core.double? electronicEnergy,
    $core.double? correlationEnergy,
    $core.Iterable<EnergyComponent>? components,
  }) {
    final result = create();
    if (totalEnergy != null) result.totalEnergy = totalEnergy;
    if (nuclearRepulsionEnergy != null) result.nuclearRepulsionEnergy = nuclearRepulsionEnergy;
    if (electronicEnergy != null) result.electronicEnergy = electronicEnergy;
    if (correlationEnergy != null) result.correlationEnergy = correlationEnergy;
    if (components != null) result.components.addAll(components);
    return result;
  }

  EnergyResults._();

  factory EnergyResults.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory EnergyResults.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'EnergyResults', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..a<$core.double>(1, _omitFieldNames ? '' : 'totalEnergy', $pb.PbFieldType.OD)
    ..a<$core.double>(2, _omitFieldNames ? '' : 'nuclearRepulsionEnergy', $pb.PbFieldType.OD)
    ..a<$core.double>(3, _omitFieldNames ? '' : 'electronicEnergy', $pb.PbFieldType.OD)
    ..a<$core.double>(4, _omitFieldNames ? '' : 'correlationEnergy', $pb.PbFieldType.OD)
    ..pc<EnergyComponent>(5, _omitFieldNames ? '' : 'components', $pb.PbFieldType.PM, subBuilder: EnergyComponent.create)
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  EnergyResults clone() => EnergyResults()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  EnergyResults copyWith(void Function(EnergyResults) updates) => super.copyWith((message) => updates(message as EnergyResults)) as EnergyResults;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static EnergyResults create() => EnergyResults._();
  @$core.override
  EnergyResults createEmptyInstance() => create();
  static $pb.PbList<EnergyResults> createRepeated() => $pb.PbList<EnergyResults>();
  @$core.pragma('dart2js:noInline')
  static EnergyResults getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<EnergyResults>(create);
  static EnergyResults? _defaultInstance;

  @$pb.TagNumber(1)
  $core.double get totalEnergy => $_getN(0);
  @$pb.TagNumber(1)
  set totalEnergy($core.double value) => $_setDouble(0, value);
  @$pb.TagNumber(1)
  $core.bool hasTotalEnergy() => $_has(0);
  @$pb.TagNumber(1)
  void clearTotalEnergy() => $_clearField(1);

  @$pb.TagNumber(2)
  $core.double get nuclearRepulsionEnergy => $_getN(1);
  @$pb.TagNumber(2)
  set nuclearRepulsionEnergy($core.double value) => $_setDouble(1, value);
  @$pb.TagNumber(2)
  $core.bool hasNuclearRepulsionEnergy() => $_has(1);
  @$pb.TagNumber(2)
  void clearNuclearRepulsionEnergy() => $_clearField(2);

  @$pb.TagNumber(3)
  $core.double get electronicEnergy => $_getN(2);
  @$pb.TagNumber(3)
  set electronicEnergy($core.double value) => $_setDouble(2, value);
  @$pb.TagNumber(3)
  $core.bool hasElectronicEnergy() => $_has(2);
  @$pb.TagNumber(3)
  void clearElectronicEnergy() => $_clearField(3);

  @$pb.TagNumber(4)
  $core.double get correlationEnergy => $_getN(3);
  @$pb.TagNumber(4)
  set correlationEnergy($core.double value) => $_setDouble(3, value);
  @$pb.TagNumber(4)
  $core.bool hasCorrelationEnergy() => $_has(3);
  @$pb.TagNumber(4)
  void clearCorrelationEnergy() => $_clearField(4);

  @$pb.TagNumber(5)
  $pb.PbList<EnergyComponent> get components => $_getList(4);
}

class EnergyComponent extends $pb.GeneratedMessage {
  factory EnergyComponent({
    $core.String? name,
    $core.double? value,
    $core.String? unit,
  }) {
    final result = create();
    if (name != null) result.name = name;
    if (value != null) result.value = value;
    if (unit != null) result.unit = unit;
    return result;
  }

  EnergyComponent._();

  factory EnergyComponent.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory EnergyComponent.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'EnergyComponent', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOS(1, _omitFieldNames ? '' : 'name')
    ..a<$core.double>(2, _omitFieldNames ? '' : 'value', $pb.PbFieldType.OD)
    ..aOS(3, _omitFieldNames ? '' : 'unit')
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  EnergyComponent clone() => EnergyComponent()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  EnergyComponent copyWith(void Function(EnergyComponent) updates) => super.copyWith((message) => updates(message as EnergyComponent)) as EnergyComponent;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static EnergyComponent create() => EnergyComponent._();
  @$core.override
  EnergyComponent createEmptyInstance() => create();
  static $pb.PbList<EnergyComponent> createRepeated() => $pb.PbList<EnergyComponent>();
  @$core.pragma('dart2js:noInline')
  static EnergyComponent getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<EnergyComponent>(create);
  static EnergyComponent? _defaultInstance;

  @$pb.TagNumber(1)
  $core.String get name => $_getSZ(0);
  @$pb.TagNumber(1)
  set name($core.String value) => $_setString(0, value);
  @$pb.TagNumber(1)
  $core.bool hasName() => $_has(0);
  @$pb.TagNumber(1)
  void clearName() => $_clearField(1);

  @$pb.TagNumber(2)
  $core.double get value => $_getN(1);
  @$pb.TagNumber(2)
  set value($core.double value) => $_setDouble(1, value);
  @$pb.TagNumber(2)
  $core.bool hasValue() => $_has(1);
  @$pb.TagNumber(2)
  void clearValue() => $_clearField(2);

  @$pb.TagNumber(3)
  $core.String get unit => $_getSZ(2);
  @$pb.TagNumber(3)
  set unit($core.String value) => $_setString(2, value);
  @$pb.TagNumber(3)
  $core.bool hasUnit() => $_has(2);
  @$pb.TagNumber(3)
  void clearUnit() => $_clearField(3);
}

class OrbitalResults extends $pb.GeneratedMessage {
  factory OrbitalResults({
    $core.double? homoEnergy,
    $core.double? lumoEnergy,
    $core.double? homoLumoGap,
    $core.Iterable<$core.double>? orbitalEnergies,
    $core.Iterable<$core.double>? occupations,
    $core.List<$core.int>? molecularOrbitals,
  }) {
    final result = create();
    if (homoEnergy != null) result.homoEnergy = homoEnergy;
    if (lumoEnergy != null) result.lumoEnergy = lumoEnergy;
    if (homoLumoGap != null) result.homoLumoGap = homoLumoGap;
    if (orbitalEnergies != null) result.orbitalEnergies.addAll(orbitalEnergies);
    if (occupations != null) result.occupations.addAll(occupations);
    if (molecularOrbitals != null) result.molecularOrbitals = molecularOrbitals;
    return result;
  }

  OrbitalResults._();

  factory OrbitalResults.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory OrbitalResults.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'OrbitalResults', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..a<$core.double>(1, _omitFieldNames ? '' : 'homoEnergy', $pb.PbFieldType.OD)
    ..a<$core.double>(2, _omitFieldNames ? '' : 'lumoEnergy', $pb.PbFieldType.OD)
    ..a<$core.double>(3, _omitFieldNames ? '' : 'homoLumoGap', $pb.PbFieldType.OD)
    ..p<$core.double>(4, _omitFieldNames ? '' : 'orbitalEnergies', $pb.PbFieldType.KD)
    ..p<$core.double>(5, _omitFieldNames ? '' : 'occupations', $pb.PbFieldType.KD)
    ..a<$core.List<$core.int>>(6, _omitFieldNames ? '' : 'molecularOrbitals', $pb.PbFieldType.OY)
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  OrbitalResults clone() => OrbitalResults()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  OrbitalResults copyWith(void Function(OrbitalResults) updates) => super.copyWith((message) => updates(message as OrbitalResults)) as OrbitalResults;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static OrbitalResults create() => OrbitalResults._();
  @$core.override
  OrbitalResults createEmptyInstance() => create();
  static $pb.PbList<OrbitalResults> createRepeated() => $pb.PbList<OrbitalResults>();
  @$core.pragma('dart2js:noInline')
  static OrbitalResults getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<OrbitalResults>(create);
  static OrbitalResults? _defaultInstance;

  @$pb.TagNumber(1)
  $core.double get homoEnergy => $_getN(0);
  @$pb.TagNumber(1)
  set homoEnergy($core.double value) => $_setDouble(0, value);
  @$pb.TagNumber(1)
  $core.bool hasHomoEnergy() => $_has(0);
  @$pb.TagNumber(1)
  void clearHomoEnergy() => $_clearField(1);

  @$pb.TagNumber(2)
  $core.double get lumoEnergy => $_getN(1);
  @$pb.TagNumber(2)
  set lumoEnergy($core.double value) => $_setDouble(1, value);
  @$pb.TagNumber(2)
  $core.bool hasLumoEnergy() => $_has(1);
  @$pb.TagNumber(2)
  void clearLumoEnergy() => $_clearField(2);

  @$pb.TagNumber(3)
  $core.double get homoLumoGap => $_getN(2);
  @$pb.TagNumber(3)
  set homoLumoGap($core.double value) => $_setDouble(2, value);
  @$pb.TagNumber(3)
  $core.bool hasHomoLumoGap() => $_has(2);
  @$pb.TagNumber(3)
  void clearHomoLumoGap() => $_clearField(3);

  @$pb.TagNumber(4)
  $pb.PbList<$core.double> get orbitalEnergies => $_getList(3);

  @$pb.TagNumber(5)
  $pb.PbList<$core.double> get occupations => $_getList(4);

  @$pb.TagNumber(6)
  $core.List<$core.int> get molecularOrbitals => $_getN(5);
  @$pb.TagNumber(6)
  set molecularOrbitals($core.List<$core.int> value) => $_setBytes(5, value);
  @$pb.TagNumber(6)
  $core.bool hasMolecularOrbitals() => $_has(5);
  @$pb.TagNumber(6)
  void clearMolecularOrbitals() => $_clearField(6);
}

class FrequencyResults extends $pb.GeneratedMessage {
  factory FrequencyResults({
    $core.Iterable<$core.double>? frequencies,
    $core.Iterable<$core.double>? intensities,
    $core.List<$core.int>? normalModes,
  }) {
    final result = create();
    if (frequencies != null) result.frequencies.addAll(frequencies);
    if (intensities != null) result.intensities.addAll(intensities);
    if (normalModes != null) result.normalModes = normalModes;
    return result;
  }

  FrequencyResults._();

  factory FrequencyResults.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory FrequencyResults.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'FrequencyResults', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..p<$core.double>(1, _omitFieldNames ? '' : 'frequencies', $pb.PbFieldType.KD)
    ..p<$core.double>(2, _omitFieldNames ? '' : 'intensities', $pb.PbFieldType.KD)
    ..a<$core.List<$core.int>>(3, _omitFieldNames ? '' : 'normalModes', $pb.PbFieldType.OY)
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  FrequencyResults clone() => FrequencyResults()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  FrequencyResults copyWith(void Function(FrequencyResults) updates) => super.copyWith((message) => updates(message as FrequencyResults)) as FrequencyResults;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static FrequencyResults create() => FrequencyResults._();
  @$core.override
  FrequencyResults createEmptyInstance() => create();
  static $pb.PbList<FrequencyResults> createRepeated() => $pb.PbList<FrequencyResults>();
  @$core.pragma('dart2js:noInline')
  static FrequencyResults getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<FrequencyResults>(create);
  static FrequencyResults? _defaultInstance;

  @$pb.TagNumber(1)
  $pb.PbList<$core.double> get frequencies => $_getList(0);

  @$pb.TagNumber(2)
  $pb.PbList<$core.double> get intensities => $_getList(1);

  @$pb.TagNumber(3)
  $core.List<$core.int> get normalModes => $_getN(2);
  @$pb.TagNumber(3)
  set normalModes($core.List<$core.int> value) => $_setBytes(2, value);
  @$pb.TagNumber(3)
  $core.bool hasNormalModes() => $_has(2);
  @$pb.TagNumber(3)
  void clearNormalModes() => $_clearField(3);
}

class PropertyResults extends $pb.GeneratedMessage {
  factory PropertyResults({
    DipoleProperty? dipoleMoment,
    QuadrupoleProperty? quadrupoleMoment,
    $core.Iterable<ChargeProperty>? atomicCharges,
  }) {
    final result = create();
    if (dipoleMoment != null) result.dipoleMoment = dipoleMoment;
    if (quadrupoleMoment != null) result.quadrupoleMoment = quadrupoleMoment;
    if (atomicCharges != null) result.atomicCharges.addAll(atomicCharges);
    return result;
  }

  PropertyResults._();

  factory PropertyResults.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory PropertyResults.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'PropertyResults', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOM<DipoleProperty>(1, _omitFieldNames ? '' : 'dipoleMoment', subBuilder: DipoleProperty.create)
    ..aOM<QuadrupoleProperty>(2, _omitFieldNames ? '' : 'quadrupoleMoment', subBuilder: QuadrupoleProperty.create)
    ..pc<ChargeProperty>(3, _omitFieldNames ? '' : 'atomicCharges', $pb.PbFieldType.PM, subBuilder: ChargeProperty.create)
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  PropertyResults clone() => PropertyResults()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  PropertyResults copyWith(void Function(PropertyResults) updates) => super.copyWith((message) => updates(message as PropertyResults)) as PropertyResults;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static PropertyResults create() => PropertyResults._();
  @$core.override
  PropertyResults createEmptyInstance() => create();
  static $pb.PbList<PropertyResults> createRepeated() => $pb.PbList<PropertyResults>();
  @$core.pragma('dart2js:noInline')
  static PropertyResults getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<PropertyResults>(create);
  static PropertyResults? _defaultInstance;

  @$pb.TagNumber(1)
  DipoleProperty get dipoleMoment => $_getN(0);
  @$pb.TagNumber(1)
  set dipoleMoment(DipoleProperty value) => $_setField(1, value);
  @$pb.TagNumber(1)
  $core.bool hasDipoleMoment() => $_has(0);
  @$pb.TagNumber(1)
  void clearDipoleMoment() => $_clearField(1);
  @$pb.TagNumber(1)
  DipoleProperty ensureDipoleMoment() => $_ensure(0);

  @$pb.TagNumber(2)
  QuadrupoleProperty get quadrupoleMoment => $_getN(1);
  @$pb.TagNumber(2)
  set quadrupoleMoment(QuadrupoleProperty value) => $_setField(2, value);
  @$pb.TagNumber(2)
  $core.bool hasQuadrupoleMoment() => $_has(1);
  @$pb.TagNumber(2)
  void clearQuadrupoleMoment() => $_clearField(2);
  @$pb.TagNumber(2)
  QuadrupoleProperty ensureQuadrupoleMoment() => $_ensure(1);

  @$pb.TagNumber(3)
  $pb.PbList<ChargeProperty> get atomicCharges => $_getList(2);
}

class DipoleProperty extends $pb.GeneratedMessage {
  factory DipoleProperty({
    $core.double? magnitude,
    $core.double? x,
    $core.double? y,
    $core.double? z,
    $core.String? unit,
  }) {
    final result = create();
    if (magnitude != null) result.magnitude = magnitude;
    if (x != null) result.x = x;
    if (y != null) result.y = y;
    if (z != null) result.z = z;
    if (unit != null) result.unit = unit;
    return result;
  }

  DipoleProperty._();

  factory DipoleProperty.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory DipoleProperty.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'DipoleProperty', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..a<$core.double>(1, _omitFieldNames ? '' : 'magnitude', $pb.PbFieldType.OD)
    ..a<$core.double>(2, _omitFieldNames ? '' : 'x', $pb.PbFieldType.OD)
    ..a<$core.double>(3, _omitFieldNames ? '' : 'y', $pb.PbFieldType.OD)
    ..a<$core.double>(4, _omitFieldNames ? '' : 'z', $pb.PbFieldType.OD)
    ..aOS(5, _omitFieldNames ? '' : 'unit')
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  DipoleProperty clone() => DipoleProperty()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  DipoleProperty copyWith(void Function(DipoleProperty) updates) => super.copyWith((message) => updates(message as DipoleProperty)) as DipoleProperty;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static DipoleProperty create() => DipoleProperty._();
  @$core.override
  DipoleProperty createEmptyInstance() => create();
  static $pb.PbList<DipoleProperty> createRepeated() => $pb.PbList<DipoleProperty>();
  @$core.pragma('dart2js:noInline')
  static DipoleProperty getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<DipoleProperty>(create);
  static DipoleProperty? _defaultInstance;

  @$pb.TagNumber(1)
  $core.double get magnitude => $_getN(0);
  @$pb.TagNumber(1)
  set magnitude($core.double value) => $_setDouble(0, value);
  @$pb.TagNumber(1)
  $core.bool hasMagnitude() => $_has(0);
  @$pb.TagNumber(1)
  void clearMagnitude() => $_clearField(1);

  @$pb.TagNumber(2)
  $core.double get x => $_getN(1);
  @$pb.TagNumber(2)
  set x($core.double value) => $_setDouble(1, value);
  @$pb.TagNumber(2)
  $core.bool hasX() => $_has(1);
  @$pb.TagNumber(2)
  void clearX() => $_clearField(2);

  @$pb.TagNumber(3)
  $core.double get y => $_getN(2);
  @$pb.TagNumber(3)
  set y($core.double value) => $_setDouble(2, value);
  @$pb.TagNumber(3)
  $core.bool hasY() => $_has(2);
  @$pb.TagNumber(3)
  void clearY() => $_clearField(3);

  @$pb.TagNumber(4)
  $core.double get z => $_getN(3);
  @$pb.TagNumber(4)
  set z($core.double value) => $_setDouble(3, value);
  @$pb.TagNumber(4)
  $core.bool hasZ() => $_has(3);
  @$pb.TagNumber(4)
  void clearZ() => $_clearField(4);

  @$pb.TagNumber(5)
  $core.String get unit => $_getSZ(4);
  @$pb.TagNumber(5)
  set unit($core.String value) => $_setString(4, value);
  @$pb.TagNumber(5)
  $core.bool hasUnit() => $_has(4);
  @$pb.TagNumber(5)
  void clearUnit() => $_clearField(5);
}

class QuadrupoleProperty extends $pb.GeneratedMessage {
  factory QuadrupoleProperty({
    $core.Iterable<$core.double>? tensor,
    $core.String? unit,
  }) {
    final result = create();
    if (tensor != null) result.tensor.addAll(tensor);
    if (unit != null) result.unit = unit;
    return result;
  }

  QuadrupoleProperty._();

  factory QuadrupoleProperty.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory QuadrupoleProperty.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'QuadrupoleProperty', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..p<$core.double>(1, _omitFieldNames ? '' : 'tensor', $pb.PbFieldType.KD)
    ..aOS(2, _omitFieldNames ? '' : 'unit')
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  QuadrupoleProperty clone() => QuadrupoleProperty()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  QuadrupoleProperty copyWith(void Function(QuadrupoleProperty) updates) => super.copyWith((message) => updates(message as QuadrupoleProperty)) as QuadrupoleProperty;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static QuadrupoleProperty create() => QuadrupoleProperty._();
  @$core.override
  QuadrupoleProperty createEmptyInstance() => create();
  static $pb.PbList<QuadrupoleProperty> createRepeated() => $pb.PbList<QuadrupoleProperty>();
  @$core.pragma('dart2js:noInline')
  static QuadrupoleProperty getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<QuadrupoleProperty>(create);
  static QuadrupoleProperty? _defaultInstance;

  @$pb.TagNumber(1)
  $pb.PbList<$core.double> get tensor => $_getList(0);

  @$pb.TagNumber(2)
  $core.String get unit => $_getSZ(1);
  @$pb.TagNumber(2)
  set unit($core.String value) => $_setString(1, value);
  @$pb.TagNumber(2)
  $core.bool hasUnit() => $_has(1);
  @$pb.TagNumber(2)
  void clearUnit() => $_clearField(2);
}

class ChargeProperty extends $pb.GeneratedMessage {
  factory ChargeProperty({
    $core.int? atomIndex,
    $core.String? method,
    $core.double? charge,
  }) {
    final result = create();
    if (atomIndex != null) result.atomIndex = atomIndex;
    if (method != null) result.method = method;
    if (charge != null) result.charge = charge;
    return result;
  }

  ChargeProperty._();

  factory ChargeProperty.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory ChargeProperty.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'ChargeProperty', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..a<$core.int>(1, _omitFieldNames ? '' : 'atomIndex', $pb.PbFieldType.O3)
    ..aOS(2, _omitFieldNames ? '' : 'method')
    ..a<$core.double>(3, _omitFieldNames ? '' : 'charge', $pb.PbFieldType.OD)
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  ChargeProperty clone() => ChargeProperty()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  ChargeProperty copyWith(void Function(ChargeProperty) updates) => super.copyWith((message) => updates(message as ChargeProperty)) as ChargeProperty;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static ChargeProperty create() => ChargeProperty._();
  @$core.override
  ChargeProperty createEmptyInstance() => create();
  static $pb.PbList<ChargeProperty> createRepeated() => $pb.PbList<ChargeProperty>();
  @$core.pragma('dart2js:noInline')
  static ChargeProperty getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<ChargeProperty>(create);
  static ChargeProperty? _defaultInstance;

  @$pb.TagNumber(1)
  $core.int get atomIndex => $_getIZ(0);
  @$pb.TagNumber(1)
  set atomIndex($core.int value) => $_setSignedInt32(0, value);
  @$pb.TagNumber(1)
  $core.bool hasAtomIndex() => $_has(0);
  @$pb.TagNumber(1)
  void clearAtomIndex() => $_clearField(1);

  @$pb.TagNumber(2)
  $core.String get method => $_getSZ(1);
  @$pb.TagNumber(2)
  set method($core.String value) => $_setString(1, value);
  @$pb.TagNumber(2)
  $core.bool hasMethod() => $_has(1);
  @$pb.TagNumber(2)
  void clearMethod() => $_clearField(2);

  @$pb.TagNumber(3)
  $core.double get charge => $_getN(2);
  @$pb.TagNumber(3)
  set charge($core.double value) => $_setDouble(2, value);
  @$pb.TagNumber(3)
  $core.bool hasCharge() => $_has(2);
  @$pb.TagNumber(3)
  void clearCharge() => $_clearField(3);
}

/// ジョブキュー関連メッセージ
class GetJobQueueRequest extends $pb.GeneratedMessage {
  factory GetJobQueueRequest({
    JobStatus? statusFilter,
  }) {
    final result = create();
    if (statusFilter != null) result.statusFilter = statusFilter;
    return result;
  }

  GetJobQueueRequest._();

  factory GetJobQueueRequest.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory GetJobQueueRequest.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'GetJobQueueRequest', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..e<JobStatus>(1, _omitFieldNames ? '' : 'statusFilter', $pb.PbFieldType.OE, defaultOrMaker: JobStatus.JOB_STATUS_UNSPECIFIED, valueOf: JobStatus.valueOf, enumValues: JobStatus.values)
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  GetJobQueueRequest clone() => GetJobQueueRequest()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  GetJobQueueRequest copyWith(void Function(GetJobQueueRequest) updates) => super.copyWith((message) => updates(message as GetJobQueueRequest)) as GetJobQueueRequest;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static GetJobQueueRequest create() => GetJobQueueRequest._();
  @$core.override
  GetJobQueueRequest createEmptyInstance() => create();
  static $pb.PbList<GetJobQueueRequest> createRepeated() => $pb.PbList<GetJobQueueRequest>();
  @$core.pragma('dart2js:noInline')
  static GetJobQueueRequest getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<GetJobQueueRequest>(create);
  static GetJobQueueRequest? _defaultInstance;

  @$pb.TagNumber(1)
  JobStatus get statusFilter => $_getN(0);
  @$pb.TagNumber(1)
  set statusFilter(JobStatus value) => $_setField(1, value);
  @$pb.TagNumber(1)
  $core.bool hasStatusFilter() => $_has(0);
  @$pb.TagNumber(1)
  void clearStatusFilter() => $_clearField(1);
}

class UpdateJobPriorityRequest extends $pb.GeneratedMessage {
  factory UpdateJobPriorityRequest({
    $core.String? jobId,
    $core.int? newPriority,
  }) {
    final result = create();
    if (jobId != null) result.jobId = jobId;
    if (newPriority != null) result.newPriority = newPriority;
    return result;
  }

  UpdateJobPriorityRequest._();

  factory UpdateJobPriorityRequest.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory UpdateJobPriorityRequest.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'UpdateJobPriorityRequest', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOS(1, _omitFieldNames ? '' : 'jobId')
    ..a<$core.int>(2, _omitFieldNames ? '' : 'newPriority', $pb.PbFieldType.O3)
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  UpdateJobPriorityRequest clone() => UpdateJobPriorityRequest()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  UpdateJobPriorityRequest copyWith(void Function(UpdateJobPriorityRequest) updates) => super.copyWith((message) => updates(message as UpdateJobPriorityRequest)) as UpdateJobPriorityRequest;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static UpdateJobPriorityRequest create() => UpdateJobPriorityRequest._();
  @$core.override
  UpdateJobPriorityRequest createEmptyInstance() => create();
  static $pb.PbList<UpdateJobPriorityRequest> createRepeated() => $pb.PbList<UpdateJobPriorityRequest>();
  @$core.pragma('dart2js:noInline')
  static UpdateJobPriorityRequest getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<UpdateJobPriorityRequest>(create);
  static UpdateJobPriorityRequest? _defaultInstance;

  @$pb.TagNumber(1)
  $core.String get jobId => $_getSZ(0);
  @$pb.TagNumber(1)
  set jobId($core.String value) => $_setString(0, value);
  @$pb.TagNumber(1)
  $core.bool hasJobId() => $_has(0);
  @$pb.TagNumber(1)
  void clearJobId() => $_clearField(1);

  @$pb.TagNumber(2)
  $core.int get newPriority => $_getIZ(1);
  @$pb.TagNumber(2)
  set newPriority($core.int value) => $_setSignedInt32(1, value);
  @$pb.TagNumber(2)
  $core.bool hasNewPriority() => $_has(1);
  @$pb.TagNumber(2)
  void clearNewPriority() => $_clearField(2);
}

class UpdateJobPriorityResponse extends $pb.GeneratedMessage {
  factory UpdateJobPriorityResponse({
    $core.bool? success,
    $core.String? message,
  }) {
    final result = create();
    if (success != null) result.success = success;
    if (message != null) result.message = message;
    return result;
  }

  UpdateJobPriorityResponse._();

  factory UpdateJobPriorityResponse.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory UpdateJobPriorityResponse.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'UpdateJobPriorityResponse', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOB(1, _omitFieldNames ? '' : 'success')
    ..aOS(2, _omitFieldNames ? '' : 'message')
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  UpdateJobPriorityResponse clone() => UpdateJobPriorityResponse()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  UpdateJobPriorityResponse copyWith(void Function(UpdateJobPriorityResponse) updates) => super.copyWith((message) => updates(message as UpdateJobPriorityResponse)) as UpdateJobPriorityResponse;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static UpdateJobPriorityResponse create() => UpdateJobPriorityResponse._();
  @$core.override
  UpdateJobPriorityResponse createEmptyInstance() => create();
  static $pb.PbList<UpdateJobPriorityResponse> createRepeated() => $pb.PbList<UpdateJobPriorityResponse>();
  @$core.pragma('dart2js:noInline')
  static UpdateJobPriorityResponse getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<UpdateJobPriorityResponse>(create);
  static UpdateJobPriorityResponse? _defaultInstance;

  @$pb.TagNumber(1)
  $core.bool get success => $_getBF(0);
  @$pb.TagNumber(1)
  set success($core.bool value) => $_setBool(0, value);
  @$pb.TagNumber(1)
  $core.bool hasSuccess() => $_has(0);
  @$pb.TagNumber(1)
  void clearSuccess() => $_clearField(1);

  @$pb.TagNumber(2)
  $core.String get message => $_getSZ(1);
  @$pb.TagNumber(2)
  set message($core.String value) => $_setString(1, value);
  @$pb.TagNumber(2)
  $core.bool hasMessage() => $_has(1);
  @$pb.TagNumber(2)
  void clearMessage() => $_clearField(2);
}

class JobQueueResponse extends $pb.GeneratedMessage {
  factory JobQueueResponse({
    $core.bool? success,
    $core.String? message,
    $core.Iterable<JobInfo>? jobs,
  }) {
    final result = create();
    if (success != null) result.success = success;
    if (message != null) result.message = message;
    if (jobs != null) result.jobs.addAll(jobs);
    return result;
  }

  JobQueueResponse._();

  factory JobQueueResponse.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory JobQueueResponse.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'JobQueueResponse', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOB(1, _omitFieldNames ? '' : 'success')
    ..aOS(2, _omitFieldNames ? '' : 'message')
    ..pc<JobInfo>(3, _omitFieldNames ? '' : 'jobs', $pb.PbFieldType.PM, subBuilder: JobInfo.create)
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  JobQueueResponse clone() => JobQueueResponse()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  JobQueueResponse copyWith(void Function(JobQueueResponse) updates) => super.copyWith((message) => updates(message as JobQueueResponse)) as JobQueueResponse;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static JobQueueResponse create() => JobQueueResponse._();
  @$core.override
  JobQueueResponse createEmptyInstance() => create();
  static $pb.PbList<JobQueueResponse> createRepeated() => $pb.PbList<JobQueueResponse>();
  @$core.pragma('dart2js:noInline')
  static JobQueueResponse getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<JobQueueResponse>(create);
  static JobQueueResponse? _defaultInstance;

  @$pb.TagNumber(1)
  $core.bool get success => $_getBF(0);
  @$pb.TagNumber(1)
  set success($core.bool value) => $_setBool(0, value);
  @$pb.TagNumber(1)
  $core.bool hasSuccess() => $_has(0);
  @$pb.TagNumber(1)
  void clearSuccess() => $_clearField(1);

  @$pb.TagNumber(2)
  $core.String get message => $_getSZ(1);
  @$pb.TagNumber(2)
  set message($core.String value) => $_setString(1, value);
  @$pb.TagNumber(2)
  $core.bool hasMessage() => $_has(1);
  @$pb.TagNumber(2)
  void clearMessage() => $_clearField(2);

  @$pb.TagNumber(3)
  $pb.PbList<JobInfo> get jobs => $_getList(2);
}

class JobInfo extends $pb.GeneratedMessage {
  factory JobInfo({
    $core.String? id,
    $core.String? calculationId,
    $core.int? priority,
    JobStatus? status,
    $core.String? assignedWorker,
    $fixnum.Int64? createdAt,
    $fixnum.Int64? startedAt,
    $fixnum.Int64? completedAt,
  }) {
    final result = create();
    if (id != null) result.id = id;
    if (calculationId != null) result.calculationId = calculationId;
    if (priority != null) result.priority = priority;
    if (status != null) result.status = status;
    if (assignedWorker != null) result.assignedWorker = assignedWorker;
    if (createdAt != null) result.createdAt = createdAt;
    if (startedAt != null) result.startedAt = startedAt;
    if (completedAt != null) result.completedAt = completedAt;
    return result;
  }

  JobInfo._();

  factory JobInfo.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory JobInfo.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'JobInfo', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOS(1, _omitFieldNames ? '' : 'id')
    ..aOS(2, _omitFieldNames ? '' : 'calculationId')
    ..a<$core.int>(3, _omitFieldNames ? '' : 'priority', $pb.PbFieldType.O3)
    ..e<JobStatus>(4, _omitFieldNames ? '' : 'status', $pb.PbFieldType.OE, defaultOrMaker: JobStatus.JOB_STATUS_UNSPECIFIED, valueOf: JobStatus.valueOf, enumValues: JobStatus.values)
    ..aOS(5, _omitFieldNames ? '' : 'assignedWorker')
    ..aInt64(6, _omitFieldNames ? '' : 'createdAt')
    ..aInt64(7, _omitFieldNames ? '' : 'startedAt')
    ..aInt64(8, _omitFieldNames ? '' : 'completedAt')
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  JobInfo clone() => JobInfo()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  JobInfo copyWith(void Function(JobInfo) updates) => super.copyWith((message) => updates(message as JobInfo)) as JobInfo;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static JobInfo create() => JobInfo._();
  @$core.override
  JobInfo createEmptyInstance() => create();
  static $pb.PbList<JobInfo> createRepeated() => $pb.PbList<JobInfo>();
  @$core.pragma('dart2js:noInline')
  static JobInfo getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<JobInfo>(create);
  static JobInfo? _defaultInstance;

  @$pb.TagNumber(1)
  $core.String get id => $_getSZ(0);
  @$pb.TagNumber(1)
  set id($core.String value) => $_setString(0, value);
  @$pb.TagNumber(1)
  $core.bool hasId() => $_has(0);
  @$pb.TagNumber(1)
  void clearId() => $_clearField(1);

  @$pb.TagNumber(2)
  $core.String get calculationId => $_getSZ(1);
  @$pb.TagNumber(2)
  set calculationId($core.String value) => $_setString(1, value);
  @$pb.TagNumber(2)
  $core.bool hasCalculationId() => $_has(1);
  @$pb.TagNumber(2)
  void clearCalculationId() => $_clearField(2);

  @$pb.TagNumber(3)
  $core.int get priority => $_getIZ(2);
  @$pb.TagNumber(3)
  set priority($core.int value) => $_setSignedInt32(2, value);
  @$pb.TagNumber(3)
  $core.bool hasPriority() => $_has(2);
  @$pb.TagNumber(3)
  void clearPriority() => $_clearField(3);

  @$pb.TagNumber(4)
  JobStatus get status => $_getN(3);
  @$pb.TagNumber(4)
  set status(JobStatus value) => $_setField(4, value);
  @$pb.TagNumber(4)
  $core.bool hasStatus() => $_has(3);
  @$pb.TagNumber(4)
  void clearStatus() => $_clearField(4);

  @$pb.TagNumber(5)
  $core.String get assignedWorker => $_getSZ(4);
  @$pb.TagNumber(5)
  set assignedWorker($core.String value) => $_setString(4, value);
  @$pb.TagNumber(5)
  $core.bool hasAssignedWorker() => $_has(4);
  @$pb.TagNumber(5)
  void clearAssignedWorker() => $_clearField(5);

  @$pb.TagNumber(6)
  $fixnum.Int64 get createdAt => $_getI64(5);
  @$pb.TagNumber(6)
  set createdAt($fixnum.Int64 value) => $_setInt64(5, value);
  @$pb.TagNumber(6)
  $core.bool hasCreatedAt() => $_has(5);
  @$pb.TagNumber(6)
  void clearCreatedAt() => $_clearField(6);

  @$pb.TagNumber(7)
  $fixnum.Int64 get startedAt => $_getI64(6);
  @$pb.TagNumber(7)
  set startedAt($fixnum.Int64 value) => $_setInt64(6, value);
  @$pb.TagNumber(7)
  $core.bool hasStartedAt() => $_has(6);
  @$pb.TagNumber(7)
  void clearStartedAt() => $_clearField(7);

  @$pb.TagNumber(8)
  $fixnum.Int64 get completedAt => $_getI64(7);
  @$pb.TagNumber(8)
  set completedAt($fixnum.Int64 value) => $_setInt64(7, value);
  @$pb.TagNumber(8)
  $core.bool hasCompletedAt() => $_has(7);
  @$pb.TagNumber(8)
  void clearCompletedAt() => $_clearField(8);
}

/// システム情報関連メッセージ
class GetSystemInfoRequest extends $pb.GeneratedMessage {
  factory GetSystemInfoRequest() => create();

  GetSystemInfoRequest._();

  factory GetSystemInfoRequest.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory GetSystemInfoRequest.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'GetSystemInfoRequest', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  GetSystemInfoRequest clone() => GetSystemInfoRequest()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  GetSystemInfoRequest copyWith(void Function(GetSystemInfoRequest) updates) => super.copyWith((message) => updates(message as GetSystemInfoRequest)) as GetSystemInfoRequest;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static GetSystemInfoRequest create() => GetSystemInfoRequest._();
  @$core.override
  GetSystemInfoRequest createEmptyInstance() => create();
  static $pb.PbList<GetSystemInfoRequest> createRepeated() => $pb.PbList<GetSystemInfoRequest>();
  @$core.pragma('dart2js:noInline')
  static GetSystemInfoRequest getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<GetSystemInfoRequest>(create);
  static GetSystemInfoRequest? _defaultInstance;
}

class SystemInfoResponse extends $pb.GeneratedMessage {
  factory SystemInfoResponse({
    $core.bool? success,
    $core.String? message,
    SystemInfo? systemInfo,
  }) {
    final result = create();
    if (success != null) result.success = success;
    if (message != null) result.message = message;
    if (systemInfo != null) result.systemInfo = systemInfo;
    return result;
  }

  SystemInfoResponse._();

  factory SystemInfoResponse.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory SystemInfoResponse.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'SystemInfoResponse', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOB(1, _omitFieldNames ? '' : 'success')
    ..aOS(2, _omitFieldNames ? '' : 'message')
    ..aOM<SystemInfo>(3, _omitFieldNames ? '' : 'systemInfo', subBuilder: SystemInfo.create)
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  SystemInfoResponse clone() => SystemInfoResponse()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  SystemInfoResponse copyWith(void Function(SystemInfoResponse) updates) => super.copyWith((message) => updates(message as SystemInfoResponse)) as SystemInfoResponse;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static SystemInfoResponse create() => SystemInfoResponse._();
  @$core.override
  SystemInfoResponse createEmptyInstance() => create();
  static $pb.PbList<SystemInfoResponse> createRepeated() => $pb.PbList<SystemInfoResponse>();
  @$core.pragma('dart2js:noInline')
  static SystemInfoResponse getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<SystemInfoResponse>(create);
  static SystemInfoResponse? _defaultInstance;

  @$pb.TagNumber(1)
  $core.bool get success => $_getBF(0);
  @$pb.TagNumber(1)
  set success($core.bool value) => $_setBool(0, value);
  @$pb.TagNumber(1)
  $core.bool hasSuccess() => $_has(0);
  @$pb.TagNumber(1)
  void clearSuccess() => $_clearField(1);

  @$pb.TagNumber(2)
  $core.String get message => $_getSZ(1);
  @$pb.TagNumber(2)
  set message($core.String value) => $_setString(1, value);
  @$pb.TagNumber(2)
  $core.bool hasMessage() => $_has(1);
  @$pb.TagNumber(2)
  void clearMessage() => $_clearField(2);

  @$pb.TagNumber(3)
  SystemInfo get systemInfo => $_getN(2);
  @$pb.TagNumber(3)
  set systemInfo(SystemInfo value) => $_setField(3, value);
  @$pb.TagNumber(3)
  $core.bool hasSystemInfo() => $_has(2);
  @$pb.TagNumber(3)
  void clearSystemInfo() => $_clearField(3);
  @$pb.TagNumber(3)
  SystemInfo ensureSystemInfo() => $_ensure(2);
}

class SystemInfo extends $pb.GeneratedMessage {
  factory SystemInfo({
    $core.String? version,
    $core.String? pythonVersion,
    $core.String? pyscfVersion,
    $core.bool? gpuAvailable,
    $core.Iterable<$core.String>? availableMethods,
    $core.Iterable<$core.String>? availableBasisSets,
    SystemResources? resources,
  }) {
    final result = create();
    if (version != null) result.version = version;
    if (pythonVersion != null) result.pythonVersion = pythonVersion;
    if (pyscfVersion != null) result.pyscfVersion = pyscfVersion;
    if (gpuAvailable != null) result.gpuAvailable = gpuAvailable;
    if (availableMethods != null) result.availableMethods.addAll(availableMethods);
    if (availableBasisSets != null) result.availableBasisSets.addAll(availableBasisSets);
    if (resources != null) result.resources = resources;
    return result;
  }

  SystemInfo._();

  factory SystemInfo.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory SystemInfo.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'SystemInfo', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOS(1, _omitFieldNames ? '' : 'version')
    ..aOS(2, _omitFieldNames ? '' : 'pythonVersion')
    ..aOS(3, _omitFieldNames ? '' : 'pyscfVersion')
    ..aOB(4, _omitFieldNames ? '' : 'gpuAvailable')
    ..pPS(5, _omitFieldNames ? '' : 'availableMethods')
    ..pPS(6, _omitFieldNames ? '' : 'availableBasisSets')
    ..aOM<SystemResources>(7, _omitFieldNames ? '' : 'resources', subBuilder: SystemResources.create)
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  SystemInfo clone() => SystemInfo()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  SystemInfo copyWith(void Function(SystemInfo) updates) => super.copyWith((message) => updates(message as SystemInfo)) as SystemInfo;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static SystemInfo create() => SystemInfo._();
  @$core.override
  SystemInfo createEmptyInstance() => create();
  static $pb.PbList<SystemInfo> createRepeated() => $pb.PbList<SystemInfo>();
  @$core.pragma('dart2js:noInline')
  static SystemInfo getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<SystemInfo>(create);
  static SystemInfo? _defaultInstance;

  @$pb.TagNumber(1)
  $core.String get version => $_getSZ(0);
  @$pb.TagNumber(1)
  set version($core.String value) => $_setString(0, value);
  @$pb.TagNumber(1)
  $core.bool hasVersion() => $_has(0);
  @$pb.TagNumber(1)
  void clearVersion() => $_clearField(1);

  @$pb.TagNumber(2)
  $core.String get pythonVersion => $_getSZ(1);
  @$pb.TagNumber(2)
  set pythonVersion($core.String value) => $_setString(1, value);
  @$pb.TagNumber(2)
  $core.bool hasPythonVersion() => $_has(1);
  @$pb.TagNumber(2)
  void clearPythonVersion() => $_clearField(2);

  @$pb.TagNumber(3)
  $core.String get pyscfVersion => $_getSZ(2);
  @$pb.TagNumber(3)
  set pyscfVersion($core.String value) => $_setString(2, value);
  @$pb.TagNumber(3)
  $core.bool hasPyscfVersion() => $_has(2);
  @$pb.TagNumber(3)
  void clearPyscfVersion() => $_clearField(3);

  @$pb.TagNumber(4)
  $core.bool get gpuAvailable => $_getBF(3);
  @$pb.TagNumber(4)
  set gpuAvailable($core.bool value) => $_setBool(3, value);
  @$pb.TagNumber(4)
  $core.bool hasGpuAvailable() => $_has(3);
  @$pb.TagNumber(4)
  void clearGpuAvailable() => $_clearField(4);

  @$pb.TagNumber(5)
  $pb.PbList<$core.String> get availableMethods => $_getList(4);

  @$pb.TagNumber(6)
  $pb.PbList<$core.String> get availableBasisSets => $_getList(5);

  @$pb.TagNumber(7)
  SystemResources get resources => $_getN(6);
  @$pb.TagNumber(7)
  set resources(SystemResources value) => $_setField(7, value);
  @$pb.TagNumber(7)
  $core.bool hasResources() => $_has(6);
  @$pb.TagNumber(7)
  void clearResources() => $_clearField(7);
  @$pb.TagNumber(7)
  SystemResources ensureResources() => $_ensure(6);
}

class SystemResources extends $pb.GeneratedMessage {
  factory SystemResources({
    $core.int? cpuCores,
    $fixnum.Int64? totalMemoryMb,
    $fixnum.Int64? availableMemoryMb,
    $fixnum.Int64? diskSpaceMb,
    $core.Iterable<GpuInfo>? gpus,
  }) {
    final result = create();
    if (cpuCores != null) result.cpuCores = cpuCores;
    if (totalMemoryMb != null) result.totalMemoryMb = totalMemoryMb;
    if (availableMemoryMb != null) result.availableMemoryMb = availableMemoryMb;
    if (diskSpaceMb != null) result.diskSpaceMb = diskSpaceMb;
    if (gpus != null) result.gpus.addAll(gpus);
    return result;
  }

  SystemResources._();

  factory SystemResources.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory SystemResources.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'SystemResources', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..a<$core.int>(1, _omitFieldNames ? '' : 'cpuCores', $pb.PbFieldType.O3)
    ..aInt64(2, _omitFieldNames ? '' : 'totalMemoryMb')
    ..aInt64(3, _omitFieldNames ? '' : 'availableMemoryMb')
    ..aInt64(4, _omitFieldNames ? '' : 'diskSpaceMb')
    ..pc<GpuInfo>(5, _omitFieldNames ? '' : 'gpus', $pb.PbFieldType.PM, subBuilder: GpuInfo.create)
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  SystemResources clone() => SystemResources()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  SystemResources copyWith(void Function(SystemResources) updates) => super.copyWith((message) => updates(message as SystemResources)) as SystemResources;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static SystemResources create() => SystemResources._();
  @$core.override
  SystemResources createEmptyInstance() => create();
  static $pb.PbList<SystemResources> createRepeated() => $pb.PbList<SystemResources>();
  @$core.pragma('dart2js:noInline')
  static SystemResources getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<SystemResources>(create);
  static SystemResources? _defaultInstance;

  @$pb.TagNumber(1)
  $core.int get cpuCores => $_getIZ(0);
  @$pb.TagNumber(1)
  set cpuCores($core.int value) => $_setSignedInt32(0, value);
  @$pb.TagNumber(1)
  $core.bool hasCpuCores() => $_has(0);
  @$pb.TagNumber(1)
  void clearCpuCores() => $_clearField(1);

  @$pb.TagNumber(2)
  $fixnum.Int64 get totalMemoryMb => $_getI64(1);
  @$pb.TagNumber(2)
  set totalMemoryMb($fixnum.Int64 value) => $_setInt64(1, value);
  @$pb.TagNumber(2)
  $core.bool hasTotalMemoryMb() => $_has(1);
  @$pb.TagNumber(2)
  void clearTotalMemoryMb() => $_clearField(2);

  @$pb.TagNumber(3)
  $fixnum.Int64 get availableMemoryMb => $_getI64(2);
  @$pb.TagNumber(3)
  set availableMemoryMb($fixnum.Int64 value) => $_setInt64(2, value);
  @$pb.TagNumber(3)
  $core.bool hasAvailableMemoryMb() => $_has(2);
  @$pb.TagNumber(3)
  void clearAvailableMemoryMb() => $_clearField(3);

  @$pb.TagNumber(4)
  $fixnum.Int64 get diskSpaceMb => $_getI64(3);
  @$pb.TagNumber(4)
  set diskSpaceMb($fixnum.Int64 value) => $_setInt64(3, value);
  @$pb.TagNumber(4)
  $core.bool hasDiskSpaceMb() => $_has(3);
  @$pb.TagNumber(4)
  void clearDiskSpaceMb() => $_clearField(4);

  @$pb.TagNumber(5)
  $pb.PbList<GpuInfo> get gpus => $_getList(4);
}

class GpuInfo extends $pb.GeneratedMessage {
  factory GpuInfo({
    $core.String? name,
    $fixnum.Int64? memoryMb,
    $core.String? driverVersion,
  }) {
    final result = create();
    if (name != null) result.name = name;
    if (memoryMb != null) result.memoryMb = memoryMb;
    if (driverVersion != null) result.driverVersion = driverVersion;
    return result;
  }

  GpuInfo._();

  factory GpuInfo.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory GpuInfo.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'GpuInfo', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOS(1, _omitFieldNames ? '' : 'name')
    ..aInt64(2, _omitFieldNames ? '' : 'memoryMb')
    ..aOS(3, _omitFieldNames ? '' : 'driverVersion')
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  GpuInfo clone() => GpuInfo()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  GpuInfo copyWith(void Function(GpuInfo) updates) => super.copyWith((message) => updates(message as GpuInfo)) as GpuInfo;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static GpuInfo create() => GpuInfo._();
  @$core.override
  GpuInfo createEmptyInstance() => create();
  static $pb.PbList<GpuInfo> createRepeated() => $pb.PbList<GpuInfo>();
  @$core.pragma('dart2js:noInline')
  static GpuInfo getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<GpuInfo>(create);
  static GpuInfo? _defaultInstance;

  @$pb.TagNumber(1)
  $core.String get name => $_getSZ(0);
  @$pb.TagNumber(1)
  set name($core.String value) => $_setString(0, value);
  @$pb.TagNumber(1)
  $core.bool hasName() => $_has(0);
  @$pb.TagNumber(1)
  void clearName() => $_clearField(1);

  @$pb.TagNumber(2)
  $fixnum.Int64 get memoryMb => $_getI64(1);
  @$pb.TagNumber(2)
  set memoryMb($fixnum.Int64 value) => $_setInt64(1, value);
  @$pb.TagNumber(2)
  $core.bool hasMemoryMb() => $_has(1);
  @$pb.TagNumber(2)
  void clearMemoryMb() => $_clearField(2);

  @$pb.TagNumber(3)
  $core.String get driverVersion => $_getSZ(2);
  @$pb.TagNumber(3)
  set driverVersion($core.String value) => $_setString(2, value);
  @$pb.TagNumber(3)
  $core.bool hasDriverVersion() => $_has(2);
  @$pb.TagNumber(3)
  void clearDriverVersion() => $_clearField(3);
}

/// ヘルスチェック関連メッセージ
class HealthCheckRequest extends $pb.GeneratedMessage {
  factory HealthCheckRequest() => create();

  HealthCheckRequest._();

  factory HealthCheckRequest.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory HealthCheckRequest.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'HealthCheckRequest', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  HealthCheckRequest clone() => HealthCheckRequest()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  HealthCheckRequest copyWith(void Function(HealthCheckRequest) updates) => super.copyWith((message) => updates(message as HealthCheckRequest)) as HealthCheckRequest;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static HealthCheckRequest create() => HealthCheckRequest._();
  @$core.override
  HealthCheckRequest createEmptyInstance() => create();
  static $pb.PbList<HealthCheckRequest> createRepeated() => $pb.PbList<HealthCheckRequest>();
  @$core.pragma('dart2js:noInline')
  static HealthCheckRequest getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<HealthCheckRequest>(create);
  static HealthCheckRequest? _defaultInstance;
}

class HealthCheckResponse extends $pb.GeneratedMessage {
  factory HealthCheckResponse({
    $core.bool? healthy,
    $core.String? status,
    $core.Iterable<$core.MapEntry<$core.String, $core.String>>? details,
    $fixnum.Int64? timestamp,
  }) {
    final result = create();
    if (healthy != null) result.healthy = healthy;
    if (status != null) result.status = status;
    if (details != null) result.details.addEntries(details);
    if (timestamp != null) result.timestamp = timestamp;
    return result;
  }

  HealthCheckResponse._();

  factory HealthCheckResponse.fromBuffer($core.List<$core.int> data, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromBuffer(data, registry);
  factory HealthCheckResponse.fromJson($core.String json, [$pb.ExtensionRegistry registry = $pb.ExtensionRegistry.EMPTY]) => create()..mergeFromJson(json, registry);

  static final $pb.BuilderInfo _i = $pb.BuilderInfo(_omitMessageNames ? '' : 'HealthCheckResponse', package: const $pb.PackageName(_omitMessageNames ? '' : 'pyscf_front'), createEmptyInstance: create)
    ..aOB(1, _omitFieldNames ? '' : 'healthy')
    ..aOS(2, _omitFieldNames ? '' : 'status')
    ..m<$core.String, $core.String>(3, _omitFieldNames ? '' : 'details', entryClassName: 'HealthCheckResponse.DetailsEntry', keyFieldType: $pb.PbFieldType.OS, valueFieldType: $pb.PbFieldType.OS, packageName: const $pb.PackageName('pyscf_front'))
    ..aInt64(4, _omitFieldNames ? '' : 'timestamp')
    ..hasRequiredFields = false
  ;

  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  HealthCheckResponse clone() => HealthCheckResponse()..mergeFromMessage(this);
  @$core.Deprecated('See https://github.com/google/protobuf.dart/issues/998.')
  HealthCheckResponse copyWith(void Function(HealthCheckResponse) updates) => super.copyWith((message) => updates(message as HealthCheckResponse)) as HealthCheckResponse;

  @$core.override
  $pb.BuilderInfo get info_ => _i;

  @$core.pragma('dart2js:noInline')
  static HealthCheckResponse create() => HealthCheckResponse._();
  @$core.override
  HealthCheckResponse createEmptyInstance() => create();
  static $pb.PbList<HealthCheckResponse> createRepeated() => $pb.PbList<HealthCheckResponse>();
  @$core.pragma('dart2js:noInline')
  static HealthCheckResponse getDefault() => _defaultInstance ??= $pb.GeneratedMessage.$_defaultFor<HealthCheckResponse>(create);
  static HealthCheckResponse? _defaultInstance;

  @$pb.TagNumber(1)
  $core.bool get healthy => $_getBF(0);
  @$pb.TagNumber(1)
  set healthy($core.bool value) => $_setBool(0, value);
  @$pb.TagNumber(1)
  $core.bool hasHealthy() => $_has(0);
  @$pb.TagNumber(1)
  void clearHealthy() => $_clearField(1);

  @$pb.TagNumber(2)
  $core.String get status => $_getSZ(1);
  @$pb.TagNumber(2)
  set status($core.String value) => $_setString(1, value);
  @$pb.TagNumber(2)
  $core.bool hasStatus() => $_has(1);
  @$pb.TagNumber(2)
  void clearStatus() => $_clearField(2);

  @$pb.TagNumber(3)
  $pb.PbMap<$core.String, $core.String> get details => $_getMap(2);

  @$pb.TagNumber(4)
  $fixnum.Int64 get timestamp => $_getI64(3);
  @$pb.TagNumber(4)
  set timestamp($fixnum.Int64 value) => $_setInt64(3, value);
  @$pb.TagNumber(4)
  $core.bool hasTimestamp() => $_has(3);
  @$pb.TagNumber(4)
  void clearTimestamp() => $_clearField(4);
}


const $core.bool _omitFieldNames = $core.bool.fromEnvironment('protobuf.omit_field_names');
const $core.bool _omitMessageNames = $core.bool.fromEnvironment('protobuf.omit_message_names');
