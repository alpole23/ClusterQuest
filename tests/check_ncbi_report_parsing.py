#!/usr/bin/env python3
"""Check that the NCBI assembly report reader survives NCBI's two key spellings.

WHY THIS EXISTS
---------------
NCBI's `datasets` ships the same records under the same filename with two
different key spellings. `datasets download genome` puts an
assembly_data_report.jsonl inside its archive using protobuf's lowerCamelCase --
taxId, organismName, infraspecificNames. `datasets summary genome
... --as-json-lines` writes snake_case -- tax_id, organism_name,
infraspecific_names.

NCBI_FETCH_METADATA moved from the first to the second in f884c0e (2026-10-27),
so that each date shard could retry independently, and kept the filename.
extract_taxonomy.py went on asking for `taxId`, got None from all 2,863 records,
skipped every one and wrote `{}`. Three runs then published a taxonomy tree
reading `total_genomes: 0` -- Brenneria (2026-11-02), Pantoea ananatis
(2026-11-02) and the Erwiniaceae regression gate (2026-11-03) -- and Nextflow
reported every one of them as a success, because a reader that finds nothing is
not a reader that failed.

So the reader is held to both spellings, and to refusing an empty result from a
non-empty report. The second half is the part that would have caught this on the
day it broke: a format change the reader does not understand now stops the run
instead of publishing a zero.
"""
import json
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'scripts' / 'taxonomy'))

import extract_taxonomy as et                                  # noqa: E402

FAIL = []


def check(ok, label, detail=''):
    print(f'  {"PASS" if ok else "FAIL"}  {label}' + (f' ({detail})' if detail else ''))
    if not ok:
        FAIL.append(label)


def report(records):
    fh = tempfile.NamedTemporaryFile('w', suffix='.jsonl', delete=False)
    for rec in records:
        fh.write(json.dumps(rec) + '\n')
    fh.close()
    return fh.name


SNAKE = {
    'accession': 'GCA_003710245.1',
    'organism': {
        'tax_id': 221276,
        'organism_name': '[Curtobacterium] plantarum',
        'infraspecific_names': {'strain': 'LMG 16222'},
    },
}

CAMEL = {
    'accession': 'GCA_000007565.2',
    'organism': {
        'taxId': 208964,
        'organismName': 'Pseudomonas aeruginosa',
        'infraspecificNames': {'strain': 'PAO1'},
    },
}


def main():
    # 1. snake_case -- what `datasets summary --as-json-lines` writes today, and
    #    what the pipeline has actually been feeding this reader since 2026-10-27.
    asm, tax = et.read_assembly_report(report([SNAKE]))
    got = asm.get('GCA_003710245.1', {})
    check(got.get('taxId') == 221276 and tax == {221276},
          'snake_case tax_id is read', f'got {got.get("taxId")!r}')
    check(got.get('organism') == '[Curtobacterium] plantarum',
          'snake_case organism_name is read', repr(got.get('organism')))
    check(got.get('strain') == 'LMG 16222',
          'snake_case infraspecific_names.strain is read', repr(got.get('strain')))

    # 2. lowerCamelCase -- what the archive from `datasets download` carries. Still
    #    accepted, because a report staged from either tool reaches the same reader.
    asm, tax = et.read_assembly_report(report([CAMEL]))
    got = asm.get('GCA_000007565.2', {})
    check(got.get('taxId') == 208964 and tax == {208964},
          'lowerCamelCase taxId is still read', f'got {got.get("taxId")!r}')
    check(got.get('organism') == 'Pseudomonas aeruginosa'
          and got.get('strain') == 'PAO1',
          'lowerCamelCase organismName and strain are still read')

    # 3. Mixed in one file. Shards are concatenated, so one run can in principle
    #    stage records from both tools; neither spelling may shadow the other.
    asm, _tax = et.read_assembly_report(report([SNAKE, CAMEL]))
    check(len(asm) == 2, 'both spellings parse from one file', f'{len(asm)} of 2')

    # 4. A report with no records is a real condition -- a taxon with no genomes
    #    at the requested assembly level -- and must stay a quiet empty result.
    asm, _tax = et.read_assembly_report(report([]))
    check(asm == {}, 'an empty report yields an empty map without complaint')

    # 5. The guard. Records in, nothing out, means the reader does not understand
    #    the format; that must stop the run rather than publish a zero.
    raised = ''
    try:
        et.read_assembly_report(report([{'accession': 'GCA_1', 'organism':
                                        {'taxonId': 1, 'name': 'x'}}]))
    except ValueError as exc:
        raised = str(exc)
    check(bool(raised), 'a non-empty report that yields nothing raises')
    check('tax_id/taxId' in raised,
          'the error names the keys it tried, so the next rename is one edit')

    return 1 if FAIL else 0


if __name__ == '__main__':
    sys.exit(main())
