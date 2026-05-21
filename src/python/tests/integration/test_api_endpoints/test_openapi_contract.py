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

import yaml
from quantum_calc.method_defaults import PARAMETER_CONSTRAINTS


HTTP_METHODS = {"GET", "POST", "PUT", "PATCH", "DELETE"}
FASTAPI_DECORATOR_METHODS = {
    "get": "GET",
    "post": "POST",
    "put": "PUT",
    "patch": "PATCH",
    "delete": "DELETE",
}

PYTHON_DIR = Path(__file__).resolve().parents[3]
API_DIR = PYTHON_DIR / "api"
OPENAPI_PATH = PYTHON_DIR.parent / "api-spec" / "openapi.yaml"
ORBITAL_GENERATOR_PATH = PYTHON_DIR / "quantum_calc" / "orbital_generator.py"

_IMPL_PATH_PARAM_PATTERN = re.compile(r"<[^>]+>")
_OPENAPI_PATH_PARAM_PATTERN = re.compile(r"\{[^}]+\}")


def _is_public_contract_path(path: str) -> bool:
    if path == "/health":
        return True

    if not path.startswith("/api/"):
        return False

    # Debug endpoints are intentionally internal and not part of OpenAPI contract.
    return not path.startswith("/api/debug/")


def _normalize_impl_path(path: str) -> str:
    path = _IMPL_PATH_PARAM_PATTERN.sub("{}", path)
    return _OPENAPI_PATH_PARAM_PATTERN.sub("{}", path)


def _normalize_openapi_path(path: str) -> str:
    return _OPENAPI_PATH_PARAM_PATTERN.sub("{}", path)


def _constant_string(node: ast.expr) -> str | None:
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    return None


def _join_paths(prefix: str, path: str) -> str:
    if not prefix:
        return path or "/"
    if not path:
        return prefix
    return f"{prefix.rstrip('/')}/{path.lstrip('/')}"


def _is_apirouter_call(node: ast.expr) -> bool:
    if not isinstance(node, ast.Call):
        return False
    if isinstance(node.func, ast.Name):
        return node.func.id == "APIRouter"
    if isinstance(node.func, ast.Attribute):
        return node.func.attr == "APIRouter"
    return False


def _extract_fastapi_router_prefixes(tree: ast.Module) -> dict[str, str]:
    prefixes: dict[str, str] = {}

    for node in tree.body:
        if not isinstance(node, ast.Assign) or not _is_apirouter_call(node.value):
            continue

        prefix = ""
        for keyword in node.value.keywords:
            if keyword.arg == "prefix":
                prefix = _constant_string(keyword.value) or ""

        for target in node.targets:
            if isinstance(target, ast.Name):
                prefixes[target.id] = prefix

    return prefixes


def _get_fastapi_decorator_path(decorator: ast.Call) -> str | None:
    if decorator.args:
        path = _constant_string(decorator.args[0])
        if path is not None:
            return path

    for keyword in decorator.keywords:
        if keyword.arg == "path":
            return _constant_string(keyword.value)

    return None


def _parse_flask_route_decorator(
    decorator: ast.expr,
) -> tuple[str, set[str]] | None:
    if not isinstance(decorator, ast.Call):
        return None

    if not isinstance(decorator.func, ast.Attribute) or decorator.func.attr != "route":
        return None

    if not decorator.args:
        return None

    path = _constant_string(decorator.args[0])
    if path is None:
        return None

    methods: set[str] = {"GET"}
    for keyword in decorator.keywords:
        if keyword.arg != "methods":
            continue

        if not isinstance(keyword.value, (ast.List, ast.Tuple)):
            continue

        parsed_methods = set()
        for item in keyword.value.elts:
            if isinstance(item, ast.Constant) and isinstance(item.value, str):
                upper_method = item.value.upper()
                if upper_method in HTTP_METHODS:
                    parsed_methods.add(upper_method)

        if parsed_methods:
            methods = parsed_methods

    return path, methods


def _parse_fastapi_route_decorator(
    decorator: ast.expr,
    router_prefixes: dict[str, str],
) -> tuple[str, set[str]] | None:
    if not isinstance(decorator, ast.Call):
        return None
    if not isinstance(decorator.func, ast.Attribute):
        return None
    if decorator.func.attr not in FASTAPI_DECORATOR_METHODS:
        return None
    if not isinstance(decorator.func.value, ast.Name):
        return None

    router_name = decorator.func.value.id
    if router_name not in router_prefixes:
        return None

    path = _get_fastapi_decorator_path(decorator)
    if path is None:
        return None

    method = FASTAPI_DECORATOR_METHODS[decorator.func.attr]
    return _join_paths(router_prefixes[router_name], path), {method}


def _parse_route_decorator(
    decorator: ast.expr,
    router_prefixes: dict[str, str],
) -> tuple[str, set[str]] | None:
    return (
        _parse_flask_route_decorator(decorator)
        or _parse_fastapi_route_decorator(decorator, router_prefixes)
    )


class _RequestArgsVisitor(ast.NodeVisitor):
    def __init__(self) -> None:
        self.query_names: set[str] = set()

    def visit_Call(self, node: ast.Call) -> None:
        if (
            isinstance(node.func, ast.Attribute)
            and node.func.attr == "get"
            and isinstance(node.func.value, ast.Attribute)
            and node.func.value.attr == "args"
            and isinstance(node.func.value.value, ast.Name)
            and node.func.value.value.id == "request"
            and node.args
            and isinstance(node.args[0], ast.Constant)
            and isinstance(node.args[0].value, str)
        ):
            self.query_names.add(node.args[0].value)

        self.generic_visit(node)


def _is_fastapi_query_default(node: ast.expr) -> bool:
    if not isinstance(node, ast.Call):
        return False
    if isinstance(node.func, ast.Name):
        return node.func.id == "Query"
    if isinstance(node.func, ast.Attribute):
        return node.func.attr == "Query"
    return False


def _extract_fastapi_query_params(node: ast.FunctionDef) -> set[str]:
    args = node.args.args
    defaults = node.args.defaults
    default_offset = len(args) - len(defaults)

    query_params: set[str] = set()
    for index, default in enumerate(defaults, start=default_offset):
        if _is_fastapi_query_default(default):
            query_params.add(args[index].arg)

    return query_params


@lru_cache(maxsize=1)
def _extract_implementation_contract() -> tuple[set[tuple[str, str]], dict[tuple[str, str], set[str]]]:
    routes: set[tuple[str, str]] = set()
    query_params_by_route: dict[tuple[str, str], set[str]] = {}

    for api_file in sorted(API_DIR.glob("*.py")):
        if api_file.name == "__init__.py":
            continue

        tree = ast.parse(api_file.read_text(encoding="utf-8"), filename=str(api_file))
        router_prefixes = _extract_fastapi_router_prefixes(tree)

        for node in tree.body:
            if not isinstance(node, ast.FunctionDef):
                continue

            route_specs: list[tuple[str, set[str]]] = []
            for decorator in node.decorator_list:
                route_spec = _parse_route_decorator(decorator, router_prefixes)
                if route_spec:
                    route_specs.append(route_spec)

            if not route_specs:
                continue

            visitor = _RequestArgsVisitor()
            visitor.visit(node)
            fastapi_query_params = _extract_fastapi_query_params(node)

            for raw_path, methods in route_specs:
                if not _is_public_contract_path(raw_path):
                    continue

                normalized_path = _normalize_impl_path(raw_path)
                for method in methods:
                    key = (method, normalized_path)
                    routes.add(key)
                    query_params_by_route.setdefault(key, set()).update(
                        visitor.query_names | fastapi_query_params
                    )

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
    only_in_openapi = sorted(openapi_routes - impl_routes)

    issues = []
    if only_in_impl:
        formatted = ", ".join(_format_route(route) for route in only_in_impl)
        issues.append(f"Implemented but missing in OpenAPI: {formatted}")
    if only_in_openapi:
        formatted = ", ".join(_format_route(route) for route in only_in_openapi)
        issues.append(f"Defined in OpenAPI but missing in implementation: {formatted}")

    assert not issues, "\n".join(issues)


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
