import asyncio
import inspect
import threading
from concurrent.futures import Future

from websocket.event_loop_bridge import (
    bind_event_loop,
    clear_event_loop,
    schedule_coroutine,
)


def test_schedule_coroutine_from_background_thread():
    async def main():
        loop = asyncio.get_running_loop()
        bind_event_loop(loop)
        results = []

        async def append_value(value):
            results.append(value)
            return value

        def worker():
            future = schedule_coroutine(append_value("sent"))
            assert future.result(timeout=2) == "sent"

        thread = threading.Thread(target=worker)
        thread.start()
        await asyncio.to_thread(thread.join, 2)

        assert not thread.is_alive()
        assert results == ["sent"]
        clear_event_loop()

    asyncio.run(main())


def test_schedule_coroutine_without_bound_loop_returns_none(caplog):
    clear_event_loop()

    async def noop():
        return None

    future = schedule_coroutine(noop())

    assert future is None
    assert "ASGI event loop is not bound" in caplog.text


def test_schedule_coroutine_runtime_error_returns_none(mocker, caplog):
    async def main():
        loop = asyncio.get_running_loop()
        bind_event_loop(loop)

        async def noop():
            return None

        mocker.patch(
            "asyncio.run_coroutine_threadsafe",
            side_effect=RuntimeError("loop unavailable"),
        )

        coro = noop()
        future = schedule_coroutine(coro)

        assert future is None
        assert inspect.getcoroutinestate(coro) == inspect.CORO_CLOSED
        assert "Failed to schedule ASGI coroutine" in caplog.text
        clear_event_loop()

    asyncio.run(main())


def test_schedule_coroutine_uses_bound_loop_snapshot(mocker):
    class ClearingLoop:
        def is_closed(self):
            clear_event_loop()
            return False

    async def noop():
        return None

    loop = ClearingLoop()
    scheduled_future = Future()
    scheduled_future.set_result(None)
    run_threadsafe = mocker.patch(
        "asyncio.run_coroutine_threadsafe",
        return_value=scheduled_future,
    )

    bind_event_loop(loop)
    coro = noop()
    future = schedule_coroutine(coro)

    assert future is scheduled_future
    run_threadsafe.assert_called_once_with(coro, loop)
    coro.close()
    clear_event_loop()
