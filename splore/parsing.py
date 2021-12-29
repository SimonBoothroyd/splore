import base64
import json
import re
from typing import Optional

from splore.models import Cursor, Page, SortBy


def encode_base64(value: str) -> str:
    return base64.urlsafe_b64encode(value.encode()).decode().rstrip("=")


def parse_base64(value: str) -> str:

    value_padding = 4 - (len(value) % 4)
    value = value + ("=" * value_padding)

    return base64.urlsafe_b64decode(value).decode()


def encode_cursor(value: Cursor) -> str:
    return encode_base64(json.dumps(value))


def parse_cursor(value: str) -> Cursor:
    cursor = json.loads(parse_base64(value))
    return tuple(cursor) if isinstance(cursor, list) else cursor


def encode_page(value: Page) -> str:
    encoded_cursor = encode_cursor(value[0])
    return f"{value[1]}({encoded_cursor})"


def parse_page(page: Optional[str]) -> Page:

    if page is None:
        return None, "next"

    direction, encoded_cursor = re.match(r"(.+)\((.+)\)", page).groups()
    cursor = parse_cursor(encoded_cursor)

    return cursor, direction


def parse_sort_by(sort_by: Optional[str]) -> Optional[SortBy]:

    if sort_by is None:
        return None

    direction, column = re.match(r"(.+)\((.+)\)", sort_by).groups()
    return column, direction


def encode_sort_by(sort_by: SortBy) -> str:
    return f"{sort_by[1]}({sort_by[0]})"
