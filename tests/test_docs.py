"""The principles document links to functions by line range; keep the
ranges in step with the code."""

import ast
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LINK = re.compile(r"\[`([\w.]+)`\]\(\.\./nmrtools/(\w+\.py)#L(\d+)-L(\d+)\)")


def definitions(path: Path) -> dict[str, tuple[int, int]]:
    """Name -> (first line incl. decorators, last line) of every function and method."""
    found: dict[str, tuple[int, int]] = {}

    def add(node, prefix=""):
        for child in node.body:
            if isinstance(child, (ast.FunctionDef, ast.ClassDef)):
                start = min([child.lineno, *(d.lineno for d in child.decorator_list)])
                found[prefix + child.name] = (start, child.end_lineno)
                if isinstance(child, ast.ClassDef):
                    add(child, prefix + child.name + ".")

    add(ast.parse(path.read_text()))
    return found


def expected(name: str, module: str) -> tuple[int, int]:
    defs = definitions(ROOT / "nmrtools" / module)
    return defs[name] if name in defs else defs["NmrShell." + name]


def test_principles_links_point_at_the_named_functions():
    links = LINK.findall((ROOT / "docs" / "principles.ko.md").read_text())
    assert len(links) > 20
    for name, module, first, last in links:
        assert (int(first), int(last)) == expected(name, module), f"{name} in {module}"
