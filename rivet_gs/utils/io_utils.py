from __future__ import annotations
import os, hashlib, json
from typing import Dict, Optional

def ensure_dir(p: str) -> None:
    os.makedirs(p, exist_ok=True)

def md5sum(path: str, chunk_size: int = 1 << 20) -> Optional[str]:
    if not path or not os.path.exists(path):
        return None
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()

def write_json(path: str, data: Dict) -> None:
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
