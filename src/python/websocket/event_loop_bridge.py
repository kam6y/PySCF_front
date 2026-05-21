import asyncio
import logging
from concurrent.futures import Future
from typing import Coroutine, TypeVar

logger = logging.getLogger(__name__)
T = TypeVar("T")

_loop: asyncio.AbstractEventLoop | None = None


def bind_event_loop(loop: asyncio.AbstractEventLoop) -> None:
    global _loop
    _loop = loop
    logger.info("ASGI event loop bound for background notifications")


def clear_event_loop() -> None:
    global _loop
    _loop = None
    logger.info("ASGI event loop cleared for background notifications")


def schedule_coroutine(coro: Coroutine[object, object, T]) -> Future[T] | None:
    if _loop is None or _loop.is_closed():
        logger.warning("ASGI event loop is not bound; dropping scheduled coroutine")
        coro.close()
        return None

    future = asyncio.run_coroutine_threadsafe(coro, _loop)

    def log_failure(done_future: Future[T]) -> None:
        try:
            done_future.result()
        except Exception:
            logger.exception("Scheduled ASGI coroutine failed")

    future.add_done_callback(log_failure)
    return future
