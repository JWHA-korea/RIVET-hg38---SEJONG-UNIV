from __future__ import annotations
import logging

_DEF_FMT = "[%(levelname)s] %(message)s"

def get_logger(name: str = "rivet-gs", level: int = logging.INFO) -> logging.Logger:
    log = logging.getLogger(name)
    if not log.handlers:
        h = logging.StreamHandler()
        h.setFormatter(logging.Formatter(_DEF_FMT))
        log.addHandler(h)
    log.setLevel(level)
    return log
