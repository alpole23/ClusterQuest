"""One page per gene cluster family, linked from the report's master table.

The report used to inline every family's detail: an 8-column representative gene table
and a consensus gene table for each of 18 families, 998 KB between them — 46% of a
2.2 MB file, all of it collapsed behind a click. Moving the detail to one page per
family leaves the report a table of families, each row linking to the family it names.

The page is assembled from HTML the report already knows how to build — the card body
from viz.clustering and the consensus table from viz.report_sections — so a family's
detail looks the same wherever it is read, and there is one implementation of each.
"""
import html as _html
import json
from pathlib import Path

_PAGE_CSS = """
    body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
           margin: 0; background: #f7f8f9; color: #212529; }
    .wrap { max-width: 1180px; margin: 0 auto; padding: 24px 20px 64px; }
    h1 { font-size: 1.55rem; margin: 0 0 4px; }
    .sub { color: #666; margin: 0 0 18px; font-size: .95rem; }
    .crumbs { font-size: .9rem; margin-bottom: 14px; }
    .crumbs a { color: #2c5aa0; text-decoration: none; }
    .crumbs a:hover { text-decoration: underline; }
    .facts { display: flex; flex-wrap: wrap; gap: 10px; margin: 0 0 22px; }
    .fact { background: #fff; border: 1px solid #e3e6e8; border-radius: 6px;
            padding: 10px 14px; min-width: 110px; }
    .fact .k { font-size: .72rem; letter-spacing: .06em; text-transform: uppercase; color: #6c757d; }
    .fact .v { font-size: 1.15rem; font-weight: 600; font-variant-numeric: tabular-nums; }
    .panel { background: #fff; border: 1px solid #e3e6e8; border-radius: 8px;
             padding: 18px 20px; margin-bottom: 20px; }
    .panel > h2 { font-size: 1.05rem; margin: 0 0 12px; }
    .gcf-card { border: 0 !important; }
    .gcf-header { display: none !important; }
    .gcf-content { display: block !important; }
    .gene-diagram svg { max-width: 100%; height: auto; }
    .gene-legend { font-size: .8rem; color: #555; margin: 10px 0; }
    table { border-collapse: collapse; width: 100%; }
    details > summary { cursor: pointer; }
    @media (prefers-color-scheme: dark) {
      body { background: #15181b; color: #e6e8ea; }
      .panel, .fact { background: #1c2024; border-color: #2c3237; }
      .sub, .fact .k { color: #9aa3ab; }
      .crumbs a { color: #7fb8d8; }
    }
"""


def create_gcf_pages(outdir, taxon, cards, consensus_blocks=None, meta=None):
    """Write gcf/GCF-<id>.html for each family; return {family_id: relative href}.

    cards             {family_id: representative-cluster HTML} from viz.clustering
    consensus_blocks  {family_id: consensus-table HTML} from viz.report_sections
    meta              {family_id: {members, genomes, product, organism, coupling_class,
                                   headgroup, priority_rank, kcb_hit, ...}}
    """
    consensus_blocks = consensus_blocks or {}
    meta = meta or {}
    gcf_dir = Path(outdir) / 'gcf'
    families = sorted(set(cards) | set(consensus_blocks), key=lambda f: (len(f), f))
    if not families:
        return {}
    gcf_dir.mkdir(parents=True, exist_ok=True)

    hrefs = {}
    for fam in families:
        m = meta.get(fam, {})
        facts = []
        for key, label in (('members', 'BGCs'), ('genomes', 'Genomes'), ('genera', 'Genera'),
                           ('coupling_class', 'Coupling class'), ('headgroup', 'Headgroup'),
                           ('priority_rank', 'Priority rank'), ('isolation', 'Isolation')):
            v = m.get(key)
            if v not in (None, '', '-'):
                facts.append(f'<div class="fact"><div class="k">{label}</div>'
                             f'<div class="v">{_html.escape(str(v))}</div></div>')

        known = m.get('kcb_hit') or ''
        product = m.get('product') or 'phosphonate'
        blocks = []
        if cards.get(fam):
            blocks.append('<div class="panel"><h2>Representative cluster</h2>'
                          f'{cards[fam]}</div>')
        if consensus_blocks.get(fam):
            blocks.append('<div class="panel"><h2>Consensus gene content across the family'
                          '</h2><p class="sub">Assembled from every member rather than one '
                          'representative. Prevalence is the column to read: a gene at 1.00 '
                          'is in every member; one at 0.24 is accessory.</p>'
                          f'{consensus_blocks[fam]}</div>')
        if not blocks:
            blocks.append('<div class="panel"><p class="sub">No detail available for this '
                          'family.</p></div>')

        page = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>GCF-{_html.escape(fam)} — {_html.escape(taxon)}</title>
<style>{_PAGE_CSS}</style></head>
<body><div class="wrap">
  <div class="crumbs"><a href="../bgc_report.html">← {_html.escape(taxon)} report</a></div>
  <h1>GCF-{_html.escape(fam)}</h1>
  <p class="sub">{_html.escape(product)}
     {' · ' + _html.escape(str(known)) if known else ' · no known-cluster match'}</p>
  <div class="facts">{''.join(facts)}</div>
  {''.join(blocks)}
  <div class="crumbs"><a href="../bgc_report.html">← back to the report</a></div>
</div></body></html>
"""
        (gcf_dir / f'GCF-{fam}.html').write_text(page, encoding='utf-8')
        hrefs[fam] = f'gcf/GCF-{fam}.html'
    return hrefs


def family_meta(gcf_data_file=None, novelty_path=None, transfer_summary=None,
                headgroup_path=None):
    """Per-family facts for the page headers, from whatever files the run produced."""
    meta = {}

    def slot(fam):
        return meta.setdefault(str(fam), {})

    if gcf_data_file and Path(gcf_data_file).exists():
        try:
            data = json.loads(Path(gcf_data_file).read_text())
            for gcf in data.get('gcfs', []):
                s = slot(gcf.get('family_id', '?'))
                s['members'] = gcf.get('member_count')
                s['product'] = gcf.get('product')
                s['organism'] = gcf.get('organism')
                kcb = gcf.get('kcb_hit')
                s['kcb_hit'] = '' if isinstance(kcb, float) else (kcb or '')
        except Exception as exc:
            print(f'Warning: GCF metadata unavailable: {exc}')

    for path, cols in ((novelty_path, {'gcf': 'gcf', 'rank': 'priority_rank',
                                       'isolation': 'isolation', 'genomes': 'genomes',
                                       'genera': 'genera', 'coupling_class': 'coupling_class'}),
                       (headgroup_path, {'gcf': 'gcf', 'headgroup': 'headgroup'})):
        if not path or not Path(path).exists():
            continue
        try:
            lines = [l.rstrip('\n').split('\t') for l in Path(path).open() if l.strip()]
            head = lines[0]
            for row in lines[1:]:
                r = dict(zip(head, row))
                fam = r.get('gcf') or r.get('family') or r.get('family_id')
                if not fam:
                    continue
                s = slot(str(fam).replace('GCF-', ''))
                for src, dest in cols.items():
                    if src in r and r[src] not in ('', '-'):
                        s[dest] = r[src]
        except Exception as exc:
            print(f'Warning: could not read {path}: {exc}')
    return meta
