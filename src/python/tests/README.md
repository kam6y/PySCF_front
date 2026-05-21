# PySCF Front Backend Test Suite

This directory contains backend tests for the FastAPI + ASGI runtime.

## Running Tests

```bash
conda activate pyscf-env
cd src/python

python -m pytest tests/ -v
python -m pytest tests/integration/test_api_endpoints -v
python -m pytest tests/integration/test_websocket_handlers.py tests/integration/test_socketio_asgi.py -v
```

## Fixtures

Key fixtures are defined in `conftest.py`:

- `app`: FastAPI application configured for tests.
- `client`: FastAPI `TestClient` for HTTP endpoint tests.
- `asgi_server`: real Uvicorn server for Socket.IO ASGI smoke tests.
- `dummy_executor`: synchronous executor for deterministic workflow tests.
- `sample_h2_xyz`, `sample_water_xyz`, `valid_dft_params`, `valid_hf_params`: common chemistry payloads.

## Examples

```python
def test_api_endpoint(client):
    response = client.get("/api/quantum/calculations")
    assert response.status_code == 200
    data = response.json()
    assert data["success"] is True
```

```python
def test_with_sample_data(sample_h2_xyz, valid_dft_params):
    assert "H" in sample_h2_xyz
    assert valid_dft_params["calculation_method"] == "DFT"
```

## Test Strategy

- Unit tests isolate service, calculation, and utility logic with mocks.
- HTTP integration tests use FastAPI `TestClient`.
- Socket.IO handler tests call registered async handlers directly for deterministic edge-case coverage.
- `test_socketio_asgi.py` starts the real ASGI app with Uvicorn and verifies token auth, room joins, and `calculation_update` delivery over WebSocket transport.
- Workflow tests use `DummyExecutor` where needed to keep process-manager behavior deterministic.

## Notes

- API contract source is `src/api-spec/openapi.yaml`.
- Generated models should be refreshed through `npm run codegen`; do not hand-edit generated files.
- Prefer the conda environment Python for backend tests:

```bash
~/miniforge3/envs/pyscf-env/bin/python -m pytest tests/ -v
```

## References

- [Pytest documentation](https://docs.pytest.org/)
- [FastAPI testing](https://fastapi.tiangolo.com/tutorial/testing/)
- [python-socketio documentation](https://python-socketio.readthedocs.io/)
