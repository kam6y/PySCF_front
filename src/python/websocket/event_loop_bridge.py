import asyncio
import logging

logger = logging.getLogger(__name__)
_loop: asyncio.AbstractEventLoop | None = None


def bind_event_loop(loop: asyncio.AbstractEventLoop) -> None:
    global _loop
    _loop = loop
    logger.info("ASGI event loop bound for background notifications")


def clear_event_loop() -> None:
    global _loop
    _loop = None
    logger.info("ASGI event loop cleared for background notifications")
