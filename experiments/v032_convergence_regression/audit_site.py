"""Check curated package pages for missing local targets and internal links."""
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit
import json
import re

root = Path('docs').resolve()
class Page(HTMLParser):
    def __init__(self, text):
        super().__init__()
        self.links = []
        self.ids = set()
        self.feed(text)
    def handle_starttag(self, tag, attrs):
        data = dict(attrs)
        if 'id' in data:
            self.ids.add(data['id'])
        for key in ('href', 'src'):
            if data.get(key):
                self.links.append(data[key])

pages = {path: Page(path.read_text()) for path in root.rglob('*.html')}
errors = []
checked = 0
for path, page in pages.items():
    for link in page.links:
        url = urlsplit(link)
        if url.scheme or url.netloc:
            continue
        name = unquote(url.path)
        if re.search(r'AGENTS|CLAUDE|PACKAGE_OVERVIEW_FOR_AGENTS|important_derivations|(^|/)(internal|log|plan)(/|$)|(^|/)Intro[.]', name):
            errors.append(f'{path.relative_to(root)}: internal link {link}')
        if name.startswith('/MPCurver/'):
            target = root / name[len('/MPCurver/'):]
        elif name.startswith('/'):
            continue
        else:
            target = (path.parent / name).resolve() if name else path
        if target.is_dir():
            target = target / 'index.html'
        checked += 1
        if not target.exists():
            errors.append(f'{path.relative_to(root)}: missing {link}')
        elif url.fragment and target in pages and unquote(url.fragment) not in pages[target].ids:
            errors.append(f'{path.relative_to(root)}: missing anchor {link}')

allowed = {'index.html','404.html','authors.html','LICENSE-text.html'}
for path in pages:
    rel = path.relative_to(root)
    if len(rel.parts) == 1 and rel.name not in allowed:
        errors.append(f'Unexpected root page: {rel}')
    if rel.parts[0] == 'articles' and rel.name not in {'index.html','mpcurve_intro.html','partition.html','fitness.html'}:
        errors.append(f'Unexpected article: {rel}')

home = (root / 'index.html').read_text()
if 'convergence' in pages[root / 'index.html'].ids:
    errors.append('Detailed convergence section remains on homepage')
if '0.3.2' not in home or 'https://github.com/AgueroZZ/MPCurver/blob/main/NEWS.md' not in home:
    errors.append('Homepage version or Changelog link missing')

for file in ['reference/fit_mpcurve.html','reference/do_mpcurve.html','articles/mpcurve_intro.html']:
    text = (root / file).read_text()
    required = 'N * D' if file == 'articles/mpcurve_intro.html' else 'normalized'
    if '0.3.2' not in text or required not in text:
        errors.append(f'New version or convergence documentation missing: {file}')

search_entries = json.loads((root / 'search.json').read_text())
for entry in search_entries:
    if not entry['path']:
        continue
    path = urlsplit(entry['path']).path.removeprefix('/MPCurver/')
    if not (root / path).exists():
        errors.append(f'Missing search target: {path}')
result = dict(pages=len(pages), local_targets_checked=checked,
              search_entries_checked=len(search_entries), errors=errors)
Path('experiments/v032_convergence_regression/site_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
raise SystemExit(bool(errors))
