import asyncio
import threading

from websocket.event_loop_bridge import bind_event_loop, clear_event_loop, schedule_coroutine


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
