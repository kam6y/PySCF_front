"""
OpenAPI contract tests.

These tests ensure that public API routes and OpenAPI definitions stay aligned:
- HTTP method + path contracts
- Query parameter names
"""

from __future__ import annotations

import ast
import re
from functools import lru_cache
from pathlib import Path

from fastapi.routing import APIRoute
import yaml
from app import create_fastapi_app
from quantum_calc.method_defaults import PARAMETER_CONSTRAINTS


HTTP_METHODS = {"GET", "POST", "PUT", "PATCH", "DELETE"}
PENDING_FASTAPI_MIGRATION_ROUTES = {
    ("DELETE", "/api/quantum/calculations/{}"),
    ("DELETE", "/api/quantum/calculations/{}/orbitals/cube-files"),
    ("GET", "/api/quantum/calculations"),
    ("GET", "/api/quantum/calculations/{}"),
    ("GET", "/api/quantum/calculations/{}/ir-spectrum"),
    ("GET", "/api/quantum/calculations/{}/orbitals"),
    ("GET", "/api/quantum/calculations/{}/orbitals/{}/cube"),
    ("GET", "/api/quantum/calculations/{}/orbitals/cube-files"),
    ("GET", "/api/quantum/status"),
    ("GET", "/api/quantum/supported-parameters"),
    ("GET", "/api/system/gpu4pyscf-status"),
    ("GET", "/api/system/resource-status"),
    ("POST", "/api/agent/chat"),
    ("POST", "/api/quantum/calculate"),
    ("POST", "/api/quantum/calculations/{}/pause"),
    ("POST", "/api/quantum/calculations/{}/resume"),
    ("POST", "/api/system/gpu4pyscf-install"),
    ("PUT", "/api/quantum/calculations/{}"),
}

PYTHON_DIR = Path(__file__).resolve().parents[3]
OPENAPI_PATH = PYTHON_DIR.parent / "api-spec" / "openapi.yaml"
ORBITAL_GENERATOR_PATH = PYTHON_DIR / "quantum_calc" / "orbital_generator.py"

_OPENAPI_PATH_PARAM_PATTERN = re.compile(r"\{[^}]+\}")


def _is_public_contract_path(path: str) -> bool:
    if path == "/health":
        return True

    if not path.startswith("/api/"):
        return False

    # Debug endpoints are intentionally internal and not part of OpenAPI contract.
    return not path.startswith("/api/debug/")


def _normalize_impl_path(path: str) -> str:
    return _OPENAPI_PATH_PARAM_PATTERN.sub("{}", path)


def _normalize_openapi_path(path: str) -> str:
    return _OPENAPI_PATH_PARAM_PATTERN.sub("{}", path)


def _extract_route_query_params(route: APIRoute) -> set[str]:
    query_params: set[str] = set()
    for field in route.dependant.query_params:
        name = getattr(field, "alias", None) or getattr(field, "name", None)
        if isinstance(name, str):
            query_params.add(name)
    return query_params


@lru_cache(maxsize=1)
def _extract_implementation_contract() -> tuple[set[tuple[str, str]], dict[tuple[str, str], set[str]]]:
    routes: set[tuple[str, str]] = set()
    query_params_by_route: dict[tuple[str, str], set[str]] = {}

    app = create_fastapi_app(server_port=5000, test_config={"TESTING": True})
    for route in app.routes:
        if not isinstance(route, APIRoute):
            continue

        if not _is_public_contract_path(route.path):
            continue

        normalized_path = _normalize_impl_path(route.path)
        query_params = _extract_route_query_params(route)
        for method in (route.methods or set()) & HTTP_METHODS:
            key = (method, normalized_path)
            routes.add(key)
            query_params_by_route[key] = query_params

    return routes, query_params_by_route


@lru_cache(maxsize=1)
def _extract_openapi_contract() -> tuple[set[tuple[str, str]], dict[tuple[str, str], set[str]]]:
    spec = yaml.safe_load(OPENAPI_PATH.read_text(encoding="utf-8"))
    path_items = spec.get("paths", {})

    routes: set[tuple[str, str]] = set()
    query_params_by_route: dict[tuple[str, str], set[str]] = {}

    for raw_path, path_item in path_items.items():
        if not _is_public_contract_path(raw_path):
            continue

        normalized_path = _normalize_openapi_path(raw_path)
        path_level_parameters = path_item.get("parameters", [])

        for method, operation in path_item.items():
            upper_method = method.upper()
            if upper_method not in HTTP_METHODS:
                continue

            key = (upper_method, normalized_path)
            routes.add(key)

            all_parameters = []
            if isinstance(path_level_parameters, list):
                all_parameters.extend(path_level_parameters)
            if isinstance(operation, dict) and isinstance(operation.get("parameters"), list):
                all_parameters.extend(operation["parameters"])

            query_names = {
                parameter["name"]
                for parameter in all_parameters
                if isinstance(parameter, dict)
                and parameter.get("in") == "query"
                and isinstance(parameter.get("name"), str)
            }
            query_params_by_route[key] = query_names

    return routes, query_params_by_route


def _format_route(route: tuple[str, str]) -> str:
    method, path = route
    return f"{method} {path}"


def _resolve_schema(spec: dict, schema: dict) -> dict:
    ref = schema.get("$ref")
    if not isinstance(ref, str):
        return schema

    prefix = "#/components/schemas/"
    assert ref.startswith(prefix), f"Unsupported schema reference: {ref}"
    schema_name = ref.removeprefix(prefix)
    return spec["components"]["schemas"][schema_name]


def _extract_string_assignments(path: Path, variable_name: str) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    values: set[str] = set()

    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign):
            continue
        if not isinstance(node.value, ast.Constant) or not isinstance(node.value.value, str):
            continue
        if any(isinstance(target, ast.Name) and target.id == variable_name for target in node.targets):
            values.add(node.value.value)

    return values


def test_openapi_and_implementation_have_same_public_routes() -> None:
    impl_routes, _ = _extract_implementation_contract()
    openapi_routes, _ = _extract_openapi_contract()

    only_in_impl = sorted(impl_routes - openapi_routes)
    only_in_openapi = sorted(
        openapi_routes - impl_routes - PENDING_FASTAPI_MIGRATION_ROUTES
    )

    issues = []
    if only_in_impl:
        formatted = ", ".join(_format_route(route) for route in only_in_impl)
        issues.append(f"Implemented but missing in OpenAPI: {formatted}")
    if only_in_openapi:
        formatted = ", ".join(_format_route(route) for route in only_in_openapi)
        issues.append(f"Defined in OpenAPI but missing in implementation: {formatted}")

    assert not issues, "\n".join(issues)


def test_pending_fastapi_migration_routes_are_defined_in_openapi() -> None:
    openapi_routes, _ = _extract_openapi_contract()

    missing_from_openapi = sorted(PENDING_FASTAPI_MIGRATION_ROUTES - openapi_routes)

    assert not missing_from_openapi, (
        "Pending migration allowlist routes must exist in OpenAPI: "
        + ", ".join(_format_route(route) for route in missing_from_openapi)
    )


def test_pending_fastapi_migration_routes_are_not_registered() -> None:
    impl_routes, _ = _extract_implementation_contract()

    registered_routes = sorted(PENDING_FASTAPI_MIGRATION_ROUTES & impl_routes)

    assert not registered_routes, (
        "Registered routes must be removed from PENDING_FASTAPI_MIGRATION_ROUTES: "
        + ", ".join(_format_route(route) for route in registered_routes)
    )


def test_openapi_and_implementation_have_same_query_parameter_names() -> None:
    impl_routes, impl_query_params = _extract_implementation_contract()
    openapi_routes, openapi_query_params = _extract_openapi_contract()

    mismatches = []
    for route in sorted(impl_routes & openapi_routes):
        impl_query = impl_query_params.get(route, set())
        openapi_query = openapi_query_params.get(route, set())

        only_in_impl = sorted(impl_query - openapi_query)
        only_in_openapi = sorted(openapi_query - impl_query)

        if only_in_impl or only_in_openapi:
            mismatch_lines = [f"{_format_route(route)}:"]
            if only_in_impl:
                mismatch_lines.append(
                    f"  - Implemented but undocumented query params: {', '.join(only_in_impl)}"
                )
            if only_in_openapi:
                mismatch_lines.append(
                    f"  - Documented but unused query params: {', '.join(only_in_openapi)}"
                )
            mismatches.append("\n".join(mismatch_lines))

    assert not mismatches, "\n".join(mismatches)


def test_openapi_active_space_limits_match_runtime_constraints() -> None:
    """CASCI/CASSCF OpenAPI limits must match runtime validation constraints."""
    spec = yaml.safe_load(OPENAPI_PATH.read_text(encoding="utf-8"))
    schemas = spec["components"]["schemas"]

    for schema_name in ("CASCICalculationRequest", "CASSCFCalculationRequest"):
        properties = schemas[schema_name]["allOf"][1]["properties"]
        for param_name in ("ncas", "nelecas", "max_cycle_micro"):
            assert properties[param_name]["maximum"] == PARAMETER_CONSTRAINTS[param_name]["max"]


def test_list_calculations_method_filter_matches_calculation_method_schema() -> None:
    """The list endpoint method filter must allow every supported calculation method."""
    spec = yaml.safe_load(OPENAPI_PATH.read_text(encoding="utf-8"))
    operation = spec["paths"]["/api/quantum/calculations"]["get"]
    method_parameter = next(
        parameter
        for parameter in operation["parameters"]
        if parameter["name"] == "calculation_method"
    )

    parameter_schema = _resolve_schema(spec, method_parameter["schema"])
    calculation_method_schema = spec["components"]["schemas"]["CalculationMethod"]

    assert parameter_schema["enum"] == calculation_method_schema["enum"]


def test_openapi_orbital_type_enum_allows_runtime_orbital_generator_values() -> None:
    """OpenAPI must allow every orbital_type emitted by MolecularOrbitalGenerator."""
    spec = yaml.safe_load(OPENAPI_PATH.read_text(encoding="utf-8"))
    orbital_type_schema = spec["components"]["schemas"]["OrbitalInfo"]["properties"][
        "orbital_type"
    ]

    runtime_orbital_types = _extract_string_assignments(
        ORBITAL_GENERATOR_PATH,
        "orbital_type",
    )

    assert runtime_orbital_types <= set(orbital_type_schema["enum"])


def test_delete_calculation_documents_validation_error_response() -> None:
    """DELETE calculation must document the 400 response used for non-terminal statuses."""
    spec = yaml.safe_load(OPENAPI_PATH.read_text(encoding="utf-8"))
    responses = spec["paths"]["/api/quantum/calculations/{calculationId}"]["delete"]["responses"]

    assert responses["400"]["content"]["application/json"]["schema"] == {
        "$ref": "#/components/schemas/ErrorResponse"
    }
