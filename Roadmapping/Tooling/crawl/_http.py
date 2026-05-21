"""Tiny urllib-based HTTP helper shared by all crawlers.

Provides retry-on-429/5xx and JSON parsing without pulling in `requests`.
"""

import json
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

USER_AGENT = "PyPhysics-crawl/0.1 (mailto:morris.trey.j@gmail.com)"
DEFAULT_TIMEOUT = 30
MAX_RETRIES = 3
RETRY_DELAY = 2.0  # seconds, exponentially backed off


def get_json(url: str, params: dict | None = None, accept: str = "application/json",
             timeout: int = DEFAULT_TIMEOUT, max_retries: int = MAX_RETRIES) -> dict | list:
    """GET a JSON resource with retries on 429/5xx. Raises on persistent failure."""
    if params:
        url = f"{url}?{urlencode(params)}"
    req = Request(url, headers={"User-Agent": USER_AGENT, "Accept": accept})
    last_err: Exception | None = None
    for attempt in range(max_retries):
        try:
            with urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except HTTPError as e:
            last_err = e
            if e.code in (429, 500, 502, 503, 504) and attempt < max_retries - 1:
                time.sleep(RETRY_DELAY * (2 ** attempt))
                continue
            raise
        except (URLError, TimeoutError) as e:
            last_err = e
            if attempt < max_retries - 1:
                time.sleep(RETRY_DELAY * (2 ** attempt))
                continue
            raise
    if last_err:
        raise last_err
    raise RuntimeError(f"unexpected fallthrough for GET {url}")


def get_text(url: str, params: dict | None = None, accept: str = "*/*",
             timeout: int = DEFAULT_TIMEOUT, max_retries: int = MAX_RETRIES) -> str:
    """GET a text resource (XML, plain) with retries. Returns the decoded body."""
    if params:
        url = f"{url}?{urlencode(params)}"
    req = Request(url, headers={"User-Agent": USER_AGENT, "Accept": accept})
    last_err: Exception | None = None
    for attempt in range(max_retries):
        try:
            with urlopen(req, timeout=timeout) as resp:
                return resp.read().decode("utf-8", errors="replace")
        except HTTPError as e:
            last_err = e
            if e.code in (429, 500, 502, 503, 504) and attempt < max_retries - 1:
                time.sleep(RETRY_DELAY * (2 ** attempt))
                continue
            raise
        except (URLError, TimeoutError) as e:
            last_err = e
            if attempt < max_retries - 1:
                time.sleep(RETRY_DELAY * (2 ** attempt))
                continue
            raise
    if last_err:
        raise last_err
    raise RuntimeError(f"unexpected fallthrough for GET {url}")
