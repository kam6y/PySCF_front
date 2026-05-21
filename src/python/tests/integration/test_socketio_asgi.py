"""End-to-end Socket.IO smoke tests against the ASGI app."""

import asyncio
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Callable

import engineio.async_client
import pytest
import socketio
import websockets
from socketio.exceptions import ConnectionError as SocketIOConnectionError


AUTH_TOKEN = "socket-token"
ORIGIN = "http://127.0.0.1:3000"
CALCULATION_ID = "calc-123"
UPDATE_TIMEOUT_SECONDS = 2.0
POLL_INTERVAL_SECONDS = 0.05


class _ClientWSTimeout:
    def __init__(self, **kwargs: Any) -> None:
        self.kwargs = kwargs


class _ClientExceptions:
    class WSServerHandshakeError(Exception):
        pass

    class ServerConnectionError(Exception):
        pass

    class ClientConnectionError(Exception):
        pass

    class ServerDisconnectedError(Exception):
        pass


class _CookieJar:
    def update_cookies(self, cookies: dict[str, str]) -> None:
        return None


class _WebSocketResponse:
    def __init__(self, websocket: Any) -> None:
        self._websocket = websocket

    async def receive(self) -> SimpleNamespace:
        try:
            return SimpleNamespace(data=await self._websocket.recv(), type=None)
        except websockets.exceptions.ConnectionClosed as exc:
            raise _ClientExceptions.ServerDisconnectedError(str(exc)) from exc

    async def send_str(self, data: str) -> None:
        try:
            await self._websocket.send(data)
        except websockets.exceptions.ConnectionClosed as exc:
            raise _ClientExceptions.ServerDisconnectedError(str(exc)) from exc

    async def send_bytes(self, data: bytes) -> None:
        try:
            await self._websocket.send(data)
        except websockets.exceptions.ConnectionClosed as exc:
            raise _ClientExceptions.ServerDisconnectedError(str(exc)) from exc

    async def close(self) -> None:
        await self._websocket.close()


class _ClientSession:
    def __init__(self) -> None:
        self.closed = False
        self.cookie_jar = _CookieJar()
        self._websockets: list[Any] = []

    async def ws_connect(self, url: str, **kwargs: Any) -> _WebSocketResponse:
        headers = dict(kwargs.get("headers") or {})
        origin = headers.pop("Origin", None) or headers.pop("origin", None)
        try:
            websocket = await websockets.connect(
                url,
                origin=origin,
                additional_headers=headers,
                proxy=None,
            )
        except websockets.exceptions.InvalidStatus as exc:
            raise _ClientExceptions.WSServerHandshakeError(str(exc)) from exc
        except OSError as exc:
            raise _ClientExceptions.ClientConnectionError(str(exc)) from exc
        self._websockets.append(websocket)
        return _WebSocketResponse(websocket)

    async def close(self) -> None:
        self.closed = True
        for websocket in self._websockets:
            await websocket.close()


class _AiohttpWebsocketsAdapter:
    ClientError = Exception
    ClientSession = _ClientSession
    ClientTimeout = _ClientWSTimeout
    ClientWSTimeout = _ClientWSTimeout
    WSMsgType = SimpleNamespace(CLOSE="close", CLOSING="closing")
    client_exceptions = _ClientExceptions


@pytest.fixture(autouse=True)
def patch_engineio_aiohttp_when_missing(monkeypatch: pytest.MonkeyPatch) -> None:
    if importlib.util.find_spec("aiohttp") is not None:
        return

    monkeypatch.setattr(
        engineio.async_client,
        "aiohttp",
        _AiohttpWebsocketsAdapter,
    )


def _run(coro: Any) -> Any:
    return asyncio.run(coro)


async def _connect(base_url: str, token: str | None) -> socketio.AsyncClient:
    client = socketio.AsyncClient(logger=False, engineio_logger=False)
    auth = {"token": token} if token is not None else None
    await client.connect(
        base_url,
        auth=auth,
        headers={"Origin": ORIGIN},
        transports=["websocket"],
        wait_timeout=2,
    )
    return client


async def _wait_for_update(
    updates: list[dict[str, Any]],
    predicate: Callable[[dict[str, Any]], bool],
) -> dict[str, Any]:
    deadline = asyncio.get_running_loop().time() + UPDATE_TIMEOUT_SECONDS
    while asyncio.get_running_loop().time() < deadline:
        for update in updates:
            if predicate(update):
                return update
        await asyncio.sleep(POLL_INTERVAL_SECONDS)
    raise AssertionError("Timed out waiting for calculation_update")


def _write_calculation(calculations_dir: Path) -> None:
    calc_dir = calculations_dir / CALCULATION_ID
    calc_dir.mkdir(parents=True, exist_ok=True)
    parameters = {
        "name": "ASGI Smoke Calculation",
        "calculation_method": "HF",
        "basis_function": "sto-3g",
        "created_at": "2024-01-01T00:00:00",
    }
    (calc_dir / "parameters.json").write_text(json.dumps(parameters), encoding="utf-8")
    (calc_dir / "status.json").write_text(
        json.dumps({"status": "running"}),
        encoding="utf-8",
    )


def test_socketio_rejects_missing_token(asgi_server: str) -> None:
    async def main() -> None:
        client = socketio.AsyncClient(logger=False, engineio_logger=False)
        try:
            with pytest.raises(SocketIOConnectionError):
                await client.connect(
                    asgi_server,
                    headers={"Origin": ORIGIN},
                    transports=["websocket"],
                    wait_timeout=2,
                )
        finally:
            if client.connected:
                await client.disconnect()

    _run(main())


def test_socketio_rejects_wrong_token(asgi_server: str) -> None:
    async def main() -> None:
        client = socketio.AsyncClient(logger=False, engineio_logger=False)
        try:
            with pytest.raises(SocketIOConnectionError):
                await client.connect(
                    asgi_server,
                    auth={"token": "wrong-token"},
                    headers={"Origin": ORIGIN},
                    transports=["websocket"],
                    wait_timeout=2,
                )
        finally:
            if client.connected:
                await client.disconnect()

    _run(main())


def test_socketio_accepts_token_and_receives_calculation_update(
    asgi_server: str,
) -> None:
    from quantum_calc import get_current_settings
    from services.notification_service import get_notification_service

    calculations_dir = Path(get_current_settings().calculations_directory)
    _write_calculation(calculations_dir)

    async def main() -> None:
        updates: list[dict[str, Any]] = []
        client = await _connect(asgi_server, AUTH_TOKEN)

        @client.on("calculation_update")
        async def on_calculation_update(payload: dict[str, Any]) -> None:
            updates.append(payload)

        try:
            await client.call("join_global_updates", timeout=2)
            await client.call(
                "join_calculation",
                {"calculation_id": CALCULATION_ID},
                timeout=2,
            )

            initial_update = await _wait_for_update(
                updates,
                lambda update: update.get("id") == CALCULATION_ID,
            )
            assert initial_update["status"] == "running"

            updates.clear()
            get_notification_service().send_calculation_update(
                CALCULATION_ID,
                "completed",
            )
            completed_update = await _wait_for_update(
                updates,
                lambda update: (
                    update.get("id") == CALCULATION_ID
                    and update.get("status") == "completed"
                ),
            )

            assert completed_update["id"] == CALCULATION_ID
            assert completed_update["status"] == "completed"
            assert completed_update["updatedAt"]
        finally:
            if client.connected:
                await client.disconnect()

    _run(main())
