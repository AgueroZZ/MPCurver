"""Audit the curated MPCurver 0.3.3 package site."""
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit
import json
import re

root = Path("docs").resolve()


class Page(HTMLParser):
    def __init__(self, text):
        super().__init__()
        self.links = []
        self.ids = set()
        self.feed(text)

    def handle_starttag(self, tag, attrs):
        data = dict(attrs)
        if "id" in data:
            self.ids.add(data["id"])
        for key in ("href", "src"):
            if data.get(key):
                self.links.append(data[key])


pages = {path: Page(path.read_text()) for path in root.rglob("*.html")}
errors = []
checked = 0
for path, page in pages.items():
    for link in page.links:
        url = urlsplit(link)
        if url.scheme or url.netloc:
            continue
        name = unquote(url.path)
        if re.search(
            r"AGENTS|CLAUDE|PACKAGE_OVERVIEW_FOR_AGENTS|important_derivations|"
            r"(^|/)(internal|log|plan)(/|$)|(^|/)Intro[.]",
            name,
        ):
            errors.append(f"{path.relative_to(root)}: internal link {link}")
        if name.startswith("/MPCurver/"):
            target = root / name[len("/MPCurver/") :]
        elif name.startswith("/"):
            continue
        else:
            target = (path.parent / name).resolve() if name else path
        if target.is_dir():
            target = target / "index.html"
        checked += 1
        if not target.exists():
            errors.append(f"{path.relative_to(root)}: missing {link}")
        elif url.fragment and target in pages and unquote(url.fragment) not in pages[target].ids:
            errors.append(f"{path.relative_to(root)}: missing anchor {link}")

allowed = {"index.html", "404.html", "authors.html", "LICENSE-text.html"}
for path in pages:
    relative = path.relative_to(root)
    if len(relative.parts) == 1 and relative.name not in allowed:
        errors.append(f"Unexpected root page: {relative}")
    if relative.parts[0] == "articles" and relative.name not in {
        "index.html",
        "mpcurve_intro.html",
        "partition.html",
        "fitness.html",
    }:
        errors.append(f"Unexpected article: {relative}")

required_text = {
    "index.html": ["0.3.3", 'intrinsic_dim = "auto"'],
    "reference/fit_mpcurve.html": [
        "0.3.3",
        "similarity_min_cluster_size",
        "spline_r2",
    ],
    "articles/partition.html": [
        "0.3.3",
        'intrinsic_dim = "auto"',
        "mean silhouette",
    ],
}
for name, required in required_text.items():
    text = (root / name).read_text()
    for phrase in required:
        if phrase not in text:
            errors.append(f"Missing {phrase!r} in {name}")

home = (root / "index.html").read_text()
if "https://github.com/AgueroZZ/MPCurver/blob/main/NEWS.md" not in home:
    errors.append("Homepage Changelog link missing")

search_entries = json.loads((root / "search.json").read_text())
for entry in search_entries:
    if not entry["path"]:
        continue
    path = urlsplit(entry["path"]).path.removeprefix("/MPCurver/")
    if not (root / path).exists():
        errors.append(f"Missing search target: {path}")

result = {
    "pages": len(pages),
    "local_targets_checked": checked,
    "search_entries_checked": len(search_entries),
    "errors": errors,
}
Path("experiments/v033_auto_dimension_release/site_audit.json").write_text(
    json.dumps(result, indent=2) + "\n"
)
print(json.dumps(result, indent=2))
raise SystemExit(bool(errors))
