#!/usr/bin/env python3
"""Build a Nextflow params file by asking, then optionally launch the run.

    python scripts/configure_run.py              # ask the parameters that matter
    python scripts/configure_run.py --advanced   # also ask the tuning knobs
    python scripts/configure_run.py -o my.json   # choose the output path

WHY THIS IS A WRAPPER AND NOT A PROMPT INSIDE main.nf
-----------------------------------------------------

The obvious design is to prompt from `main.nf` itself. It cannot work, and it fails
*silently*, which is worse than failing. Measured on Nextflow 26.04.3:

    println "console: ${System.console()}"   -> java.io.ProxyingConsole@...   (a console IS attached)
    params.taxon                             -> Erwiniaceae
    params.taxon = 'MutatedTaxon'            -> no exception
    params.taxon                             -> Erwiniaceae                  (unchanged)

So a script CAN read stdin under `nextflow run`, and the answer CANNOT be applied:
`params` is immutable and assigning to it is a no-op that raises nothing. Every gate in
the workflow body (`if (params.pepm_prescreen)`) and every process directive reads
`params`, so a prompt there would collect an answer and then ignore it. Three further
reasons the wrapper is the right shape even if that were fixed:

  * `validateParameters()` runs before the workflow body, so prompts would land after
    validation had already rejected or accepted whatever arrived on the command line.
  * `-resume` would re-prompt on every resume, and an answer that differed from the
    cached run would change the task hashes without saying so.
  * `sbatch submit_slurm.sh` and any CI invocation have no terminal at all.

Writing a params file instead has a property worth more than the convenience: the file
is a **record of exactly what was run**, which can be committed beside a result and
replayed. `docs/comparisons/` already holds several (`run_erwin_screen_off.json` and
friends) for precisely that reason.

It also sidesteps the CLI boolean trap documented in CLAUDE.md: `--run_gtdbtk false`
arrives as the *string* `"false"`, which is truthy in Groovy, so it ENABLES what it
appears to disable. A params file carries a real JSON boolean. Everything this script
writes goes through `-params-file`, never `--flag value`.

WHAT IT ASKS
------------

`workflow` first, then only the parameters that mode actually consumes -- a `download`
run is never asked about BiG-SCAPE. Within that, questions are gated on their own
answers: BiG-SCAPE parameters appear only if clustering is on, partition parameters only
if partitioning is on, and so forth.

Types, defaults and help text are read from `nextflow_schema.json` at runtime rather
than restated here, so this script cannot drift from the schema the pipeline validates
against. `tests/check_configure_run.py` asserts every name mentioned below exists there.

The curation is the one thing that IS hand-written: 61 parameters is far too many to ask
about, so each is tagged `ASK` (always), `ADV` (only with `--advanced`), or left out of
the questionnaire entirely and allowed to keep its default. Infrastructure and tuning --
batch sizes, database pins, publish mode -- are deliberately in the last group; they are
set from measurement in `nextflow.config`, and a prompt invites a worse value.
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCHEMA = ROOT / 'nextflow_schema.json'

ASK, ADV = 'ask', 'advanced'

# Choice lists. The schema declares no enums, so these are restated from the validation
# in main.nf -- validateParams() rejects anything else, and a prompt that offers an
# invalid option would send the user into that error instead of away from it. Keep in
# step with main.nf; tests/check_configure_run.py checks these against it.
CHOICES = {
    'workflow':                ['download', 'bgc_analysis', 'full'],
    'clustering':              ['bigscape', 'none'],
    # 'all' is deliberately NOT offered. It is both GenBank and RefSeq, which is most
    # genomes twice -- a GCA_ and a GCF_ accession for one assembly. Measured on
    # Brenneria: GenBank 60, RefSeq 58, all 119, of which 59 are the same assembly under
    # both accessions. Both twins survive renaming, so both are screened, both reach
    # antiSMASH and both land in BiG-SCAPE, inflating family member counts in a way that
    # reads as biological. main.nf rejects it too.
    'assembly_source':         ['GenBank', 'RefSeq'],
    'bigscape_alignment_mode': ['auto', 'global', 'glocal'],
    'bigscape_classify':       ['category', 'class', 'legacy', ''],
    'publish_mode':            ['link', 'copy', 'symlink'],
}

# Questions, in the order they are asked. Each entry is
#   (param, tier, modes, predicate_or_None)
# `modes` is which values of `workflow` consume it. `predicate` takes the answers so far
# and decides whether the question still applies -- that is what makes the interview
# branch rather than march through everything.
QUESTIONS = [
    # ---- what to run, and where it goes ---------------------------------------
    ('outdir',                   ASK, {'download', 'bgc_analysis', 'full'}, None),

    # ---- input -----------------------------------------------------------------
    ('taxon',                    ASK, {'download', 'bgc_analysis', 'full'}, None),
    ('input_genomes',            ASK, {'bgc_analysis'}, None),
    ('assembly_source',          ASK, {'download', 'full'}, None),
    ('assembly_reference_only',  ASK, {'download', 'full'}, None),
    ('assembly_level',           ADV, {'download', 'full'},
     lambda a: not a.get('assembly_reference_only')),
    ('download_batch_size',      ADV, {'download', 'full'}, None),
    ('download_max_forks',       ADV, {'download', 'full'}, None),
    ('metadata_shard_threshold', ADV, {'download', 'full'}, None),

    # ---- the pepM pre-screen ---------------------------------------------------
    ('pepm_prescreen',           ASK, {'download', 'bgc_analysis', 'full'}, None),
    ('pepm_prescreen_bitscore',  ADV, {'download', 'bgc_analysis', 'full'},
     lambda a: a.get('pepm_prescreen')),
    ('pepm_prescreen_min_density', ADV, {'download', 'bgc_analysis', 'full'},
     lambda a: a.get('pepm_prescreen')),
    ('pepm_prescreen_batch_size', ADV, {'download', 'bgc_analysis', 'full'},
     lambda a: a.get('pepm_prescreen')),

    # ---- antiSMASH -------------------------------------------------------------
    ('recover_orfs',             ASK, {'bgc_analysis', 'full'}, None),
    ('annotation_transfer',      ADV, {'bgc_analysis', 'full'},
     lambda a: a.get('recover_orfs')),
    ('transfer_min_identity',    ADV, {'bgc_analysis', 'full'},
     lambda a: a.get('annotation_transfer')),
    ('transfer_min_coverage',    ADV, {'bgc_analysis', 'full'},
     lambda a: a.get('annotation_transfer')),
    ('antismash_minimal',        ASK, {'bgc_analysis', 'full'}, None),
    # These three are ignored in minimal mode -- nextflow.config says so on each line --
    # so asking about them there would be asking about nothing.
    ('antismash_cb_general',     ADV, {'bgc_analysis', 'full'},
     lambda a: not a.get('antismash_minimal')),
    ('antismash_cc_mibig',       ADV, {'bgc_analysis', 'full'},
     lambda a: not a.get('antismash_minimal')),
    ('antismash_smcog_trees',    ADV, {'bgc_analysis', 'full'},
     lambda a: not a.get('antismash_minimal')),
    ('antismash_phosphonate_neighbourhood', ADV, {'bgc_analysis', 'full'}, None),
    ('antismash_batch_size',     ADV, {'bgc_analysis', 'full'}, None),
    ('antismash_html',           ADV, {'bgc_analysis', 'full'}, None),
    ('antismash_summary_gbk',    ADV, {'bgc_analysis', 'full'}, None),
    ('reuse_antismash_from',     ADV, {'bgc_analysis', 'full'}, None),

    # ---- clustering ------------------------------------------------------------
    ('clustering',               ASK, {'bgc_analysis', 'full'}, None),
    ('bigscape_cutoffs',         ASK, {'bgc_analysis', 'full'},
     lambda a: a.get('clustering') == 'bigscape'),
    ('bigscape_min_bgc_length',  ADV, {'bgc_analysis', 'full'},
     lambda a: a.get('clustering') == 'bigscape'),
    ('bigscape_partition',       ASK, {'bgc_analysis', 'full'},
     lambda a: a.get('clustering') == 'bigscape'),
    ('bigscape_partition_identity', ADV, {'bgc_analysis', 'full'},
     lambda a: a.get('bigscape_partition')),
    ('bigscape_partition_threshold', ADV, {'bgc_analysis', 'full'},
     lambda a: a.get('bigscape_partition')),
    ('bigscape_partition_max_size', ADV, {'bgc_analysis', 'full'},
     lambda a: a.get('bigscape_partition')),
    ('bigscape_partition_memory_margin', ADV, {'bgc_analysis', 'full'},
     lambda a: a.get('bigscape_partition')),
    ('bigscape_alignment_mode',  ADV, {'bgc_analysis', 'full'},
     lambda a: a.get('clustering') == 'bigscape'),
    ('bigscape_classify',        ADV, {'bgc_analysis', 'full'},
     lambda a: a.get('clustering') == 'bigscape'),
    ('bigscape_include_singletons', ADV, {'bgc_analysis', 'full'},
     lambda a: a.get('clustering') == 'bigscape'),
    ('bigscape_mix',             ADV, {'bgc_analysis', 'full'},
     lambda a: a.get('clustering') == 'bigscape'),
    ('bigscape_mibig_version',   ADV, {'bgc_analysis', 'full'},
     lambda a: a.get('clustering') == 'bigscape'),

    # ---- phylogeny -------------------------------------------------------------
    ('run_gtdbtk',               ASK, {'bgc_analysis', 'full'}, None),
    ('gtdbtk_bgc_genomes_only',  ADV, {'bgc_analysis', 'full'},
     lambda a: a.get('run_gtdbtk')),
    ('gtdbtk_skip_ani_screen',   ADV, {'bgc_analysis', 'full'},
     lambda a: a.get('run_gtdbtk')),
    ('gtdbtk_cpus',              ADV, {'bgc_analysis', 'full'},
     lambda a: a.get('run_gtdbtk')),
    ('gtdbtk_pplacer_cpus',      ADV, {'bgc_analysis', 'full'},
     lambda a: a.get('run_gtdbtk')),
    ('gtdbtk_shard_size',        ADV, {'bgc_analysis', 'full'},
     lambda a: a.get('run_gtdbtk')),
    ('gtdbtk_min_perc_aa',       ADV, {'bgc_analysis', 'full'},
     lambda a: a.get('run_gtdbtk')),
    ('reuse_gtdbtk_from',        ADV, {'bgc_analysis', 'full'},
     lambda a: a.get('run_gtdbtk')),

    # ---- report ----------------------------------------------------------------
    ('run_analysis',             ADV, {'bgc_analysis', 'full'}, None),
    ('skip_tree',                ADV, {'bgc_analysis', 'full'}, None),
    ('count_per_contig',         ADV, {'bgc_analysis', 'full'}, None),
    ('split_hybrids',            ADV, {'bgc_analysis', 'full'}, None),

    # ---- infrastructure --------------------------------------------------------
    ('publish_mode',             ADV, {'download', 'bgc_analysis', 'full'}, None),
    ('task_batch_size',          ADV, {'download', 'bgc_analysis', 'full'}, None),
    ('check_db_updates',         ADV, {'download', 'bgc_analysis', 'full'}, None),
]

# Warnings printed when an answer commits the user to something expensive or
# irreversible. Keyed by (param, value); the text is shown immediately after the answer,
# while it can still be changed, rather than in a summary nobody reads.
WARNINGS = {
    ('run_gtdbtk', True): (
        'GTDB-Tk needs ~140 GB of reference database (downloaded once into\n'
        '    <outdir>/databases/gtdbtk) and ~95 GB of scratch per concurrent task.\n'
        '    On a workstation this is usually the largest single cost of a run.'),
    ('pepm_prescreen', False): (
        'With the screen off, antiSMASH sees EVERY genome. On Erwiniaceae that is\n'
        '    5.5x the antiSMASH CPU for an identical BGC count; on a dilute taxon it is\n'
        '    far more. Turn it off only to measure what it saves.'),
    ('bigscape_partition', True): (
        'Partitioning does nothing below bigscape_partition_threshold BGCs (default\n'
        '    10,000) and no run to date has been large enough to need it. It is a memory\n'
        '    backstop for a dense taxon, and it costs more total compute.'),
    ('antismash_minimal', True): (
        'Minimal mode skips --clusterhmmer and --tigrfam, so regions carry no Pfam\n'
        '    domain annotation. Faster, and it changes what the report can show.'),
}


def load_schema():
    if not SCHEMA.exists():
        sys.exit(f'error: {SCHEMA} not found -- run this from the pipeline checkout')
    return json.loads(SCHEMA.read_text())['properties']


# --------------------------------------------------------------------------- asking

def _types(spec):
    """The schema's `type` as a set of names; it may be a list for nullable params."""
    t = spec.get('type')
    return set(t) if isinstance(t, list) else {t}


def _fmt(value):
    """Render a default the way the user will have to type it back."""
    if value is None:
        return 'null'
    if isinstance(value, bool):
        return 'yes' if value else 'no'
    if value == '':
        return '(empty)'
    return str(value)


def _wrap(text, width=76, indent='    '):
    out, line = [], ''
    for word in re.sub(r'\s+', ' ', text).strip().split():
        if line and len(line) + 1 + len(word) > width:
            out.append(indent + line)
            line = word
        else:
            line = f'{line} {word}'.strip()
    if line:
        out.append(indent + line)
    return '\n'.join(out)


def coerce(raw, spec, param):
    """Turn typed text into a JSON value of the schema's type, or raise ValueError.

    Nullable params accept the literal `null`, because that is the documented "off"
    value for several of them and `antismash_phosphonate_neighbourhood` in particular
    means "restore antiSMASH's own default" rather than "zero".
    """
    types = _types(spec)
    raw = raw.strip()

    if 'null' in types and raw.lower() in ('null', 'none', '-'):
        return None
    if 'boolean' in types:
        if raw.lower() in ('y', 'yes', 'true', 't', '1'):
            return True
        if raw.lower() in ('n', 'no', 'false', 'f', '0'):
            return False
        raise ValueError('answer yes or no')
    if param in CHOICES:
        if raw in CHOICES[param]:
            return raw
        raise ValueError(f'choose one of: {", ".join(c or "(empty)" for c in CHOICES[param])}')
    if 'integer' in types:
        try:
            return int(raw)
        except ValueError:
            raise ValueError('whole number expected')
    if 'number' in types:
        try:
            return float(raw)
        except ValueError:
            raise ValueError('number expected')
    return raw


def ask(param, spec, default):
    """Prompt until the answer parses. Empty input accepts the default."""
    print()
    print(f'  \033[1m{param}\033[0m')
    desc = spec.get('description') or ''
    if desc:
        print(_wrap(desc))
        # The schema generator caps descriptions at 240 characters, so the longest ones
        # stop mid-sentence. Say where the rest is rather than leaving a dangling clause.
        if len(desc) >= 240:
            print(f'    \033[2m(truncated; full note on {param} is in nextflow.config)\033[0m')
    if param in CHOICES:
        print(_wrap('options: ' + ', '.join(c or '(empty)' for c in CHOICES[param])))

    hint = _fmt(default)
    if 'boolean' in _types(spec):
        hint = f'{hint} [yes/no]'
    while True:
        try:
            raw = input(f'    > [{hint}] ').strip()
        except EOFError:
            print()
            sys.exit('aborted: stdin closed')
        except KeyboardInterrupt:
            print()
            sys.exit('aborted')
        if not raw:
            return default
        try:
            return coerce(raw, spec, param)
        except ValueError as exc:
            print(f'    \033[31m{exc}\033[0m')


def validate_answer(param, value):
    """Post-hoc checks the schema cannot express. Returns a complaint or None."""
    if param == 'input_genomes':
        if not value:
            return 'bgc_analysis needs a genome directory'
        p = Path(str(value)).expanduser()
        if not p.is_dir():
            return f'not a directory: {p}'
        if not any(p.glob('*.gbff')):
            return (f'no *.gbff files in {p} -- main.nf reads '
                    '"${params.input_genomes}/*.gbff"')
    if param == 'outdir' and not str(value).strip():
        return 'outdir cannot be empty'
    if param == 'taxon' and not str(value).strip():
        return 'taxon cannot be empty'
    return None


# --------------------------------------------------------------------------- driving

def interview(schema, advanced):
    answers = {}

    print('\n\033[1mClusterQuest run configuration\033[0m')
    print(_wrap('Press Enter to accept the value in brackets. Answers are written to a '
                'params file; nothing runs until you confirm.', indent='  '))

    wf_spec = schema['workflow']
    answers['workflow'] = ask('workflow', wf_spec, wf_spec.get('default', 'full'))
    mode = answers['workflow']

    for param, tier, modes, predicate in QUESTIONS:
        if param == 'workflow' or mode not in modes:
            continue
        if tier == ADV and not advanced:
            continue
        if predicate and not predicate(answers):
            continue
        spec = schema.get(param)
        if spec is None:          # cannot happen while the test below passes
            continue
        while True:
            value = ask(param, spec, spec.get('default'))
            complaint = validate_answer(param, value)
            if complaint is None:
                break
            print(f'    \033[31m{complaint}\033[0m')
        answers[param] = value
        warn = WARNINGS.get((param, value))
        if warn:
            print(f'    \033[33mnote:\033[0m {warn}')

    return answers


def summarise(answers, schema):
    """Show what will be written, marking the values that differ from the default."""
    print('\n\033[1mParameters\033[0m')
    for k, v in answers.items():
        default = schema.get(k, {}).get('default')
        changed = ' \033[36m<- changed\033[0m' if v != default else ''
        print(f'  {k:38s} {_fmt(v)}{changed}')


def default_outfile(answers):
    slug = re.sub(r'[^A-Za-z0-9]+', '_', str(answers.get('taxon', 'run'))).strip('_').lower()
    stamp = datetime.now().strftime('%Y%m%d_%H%M')
    return Path(f'run_{slug}_{answers["workflow"]}_{stamp}.json')


def main():
    ap = argparse.ArgumentParser(
        description='Build a Nextflow params file by asking, then optionally launch.')
    ap.add_argument('--advanced', action='store_true',
                    help='also ask the tuning parameters, not just the common ones')
    ap.add_argument('-o', '--out', type=Path,
                    help='where to write the params file (default: run_<taxon>_<mode>_<stamp>.json)')
    ap.add_argument('--print-only', action='store_true',
                    help='write the file and print the command; never launch')
    ap.add_argument('--list', action='store_true',
                    help='list which parameters each mode asks about, and exit')
    args = ap.parse_args()

    schema = load_schema()

    if args.list:
        for mode in CHOICES['workflow']:
            common = [p for p, t, m, _ in QUESTIONS if mode in m and t == ASK]
            adv = [p for p, t, m, _ in QUESTIONS if mode in m and t == ADV]
            print(f'\n{mode}:')
            print(f'  asked:    workflow, {", ".join(common)}')
            print(f'  advanced: {", ".join(adv)}')
        return 0

    if not sys.stdin.isatty():
        sys.exit('error: no terminal on stdin. This script is interactive by design; '
                 'for a non-interactive run write the params file directly and pass it '
                 'with -params-file.')

    answers = interview(schema, args.advanced)
    summarise(answers, schema)

    out = args.out or default_outfile(answers)
    if out.exists():
        reply = input(f'\n{out} exists. Overwrite? [no] ').strip().lower()
        if reply not in ('y', 'yes'):
            sys.exit('aborted: nothing written')
    out.write_text(json.dumps(answers, indent=2) + '\n')
    print(f'\nwrote \033[1m{out}\033[0m')

    # -profile is a Nextflow option, not a pipeline parameter, so it cannot live in the
    # params file and has to go on the command line. It is always included: a run
    # without it gets the `standard` profile, which does NOT carry the local memory
    # caps in nextflow.config, so BiG-SCAPE would ask for 32 GB on a 56 GB box.
    profile = 'local'
    if not args.print_only:
        profile = input('\nprofile [local] ').strip() or 'local'
    cmd = ['nextflow', 'run', str(ROOT / 'main.nf'),
           '-profile', profile, '-params-file', str(out)]

    print('\n\033[1mcommand\033[0m')
    print('  ' + ' '.join(cmd))

    if args.print_only:
        return 0
    if input('\nrun it now? [no] ').strip().lower() not in ('y', 'yes'):
        print('not launched. The command above reproduces this configuration.')
        return 0

    if not shutil.which('nextflow'):
        sys.exit('error: nextflow is not on PATH. Activate its environment first:\n'
                 '  export JAVA_HOME="$HOME/miniconda3/envs/nextflow"\n'
                 '  export PATH="$PATH:$HOME/miniconda3/envs/nextflow/bin"')
    if not os.environ.get('JAVA_HOME'):
        print('\033[33mnote:\033[0m JAVA_HOME is unset; Nextflow may not find a JVM.')
    return subprocess.call(cmd, cwd=ROOT)


if __name__ == '__main__':
    sys.exit(main())
