#!/usr/bin/env python3
"""The antiSMASH detection rule names must have ONE value, spelled in two places.

`Utils.ANTISMASH_RULE_NAMES` is what the report prints in its header; the literal in
`modules/analysis/antismash.nf` is what antiSMASH is actually invoked with. A reader
has to be able to trust that the header describes the run, so the two must agree.

The obvious way to guarantee that is to interpolate the constant into the module. That
was done on 2026-11-05 and reverted the same day: **Nextflow hashes the unevaluated
script block**, so editing that one line changed every ANTISMASH task hash while the
rendered command stayed byte-identical, and `-resume` began recomputing ~3 h of
antiSMASH and GTDB-Tk output it already had on disk. antiSMASH is the most expensive
process in the pipeline and the one most likely to be resumed into.

So the module keeps a literal and this check holds the two spellings together. It costs
nothing at runtime and invalidates nothing.

Note `Utils.antismashParamsHash` separately hardcodes `hmmdetection_rules=phosphonate`,
which is what governs `--reuse_antismash_from` ACROSS runs. That is a third spelling of
the same value, so it is checked here too.
"""

import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent

UTILS = ROOT / 'lib' / 'Utils.groovy'
MODULE = ROOT / 'modules' / 'analysis' / 'antismash.nf'

# static final String ANTISMASH_RULE_NAMES = 'phosphonate'
_CONST = re.compile(
    r"\bANTISMASH_RULE_NAMES\s*=\s*(['\"])(?P<value>[^'\"]+)\1")

# def hmmdetection_flag = '--hmmdetection-limit-to-rule-names phosphonate'
_FLAG = re.compile(
    r"--hmmdetection-limit-to-rule-names\s+(?P<value>[^'\"\s}]+)")

# "hmmdetection_rules=phosphonate"
_PARAMS_HASH = re.compile(
    r"hmmdetection_rules=(?P<value>[^'\"\s;]+)")


def _one(pattern, text, what, path):
    hits = [m.group('value') for m in pattern.finditer(text)]
    if not hits:
        return None, f"{path.name}: no {what} found -- has it been renamed?"
    if len(set(hits)) > 1:
        return None, f"{path.name}: {what} disagrees with itself: {sorted(set(hits))}"
    return hits[0], None


def main():
    problems = []

    utils_text = UTILS.read_text(encoding='utf-8')
    module_text = MODULE.read_text(encoding='utf-8')

    const, err = _one(_CONST, utils_text, 'ANTISMASH_RULE_NAMES', UTILS)
    if err:
        problems.append(err)
    flag, err = _one(_FLAG, module_text, '--hmmdetection-limit-to-rule-names value', MODULE)
    if err:
        problems.append(err)
    phash, err = _one(_PARAMS_HASH, utils_text, 'hmmdetection_rules in antismashParamsHash', UTILS)
    if err:
        problems.append(err)

    # An interpolation is how this drifted into invalidating the cache; reject it by
    # name so the revert cannot be quietly undone.
    if flag and ('$' in flag or '{' in flag):
        problems.append(
            f"{MODULE.name}: the rule-names flag interpolates ({flag!r}) instead of "
            "carrying a literal. Nextflow hashes the unevaluated script block, so "
            "editing that line re-runs EVERY cached antiSMASH task even when the "
            "rendered command does not change. Keep the literal; this check is what "
            "keeps it honest.")

    spellings = {'Utils.ANTISMASH_RULE_NAMES': const,
                 f'{MODULE.name} script block': flag,
                 'Utils.antismashParamsHash': phash}
    present = {k: v for k, v in spellings.items() if v is not None}
    if len(set(present.values())) > 1:
        problems.append(
            "the detection rule names disagree between: "
            + ", ".join(f"{k} = {v!r}" for k, v in sorted(present.items()))
            + ". The report header would then describe a scope antiSMASH did not run.")

    if problems:
        for p in problems:
            print(f"FAIL: {p}", file=sys.stderr)
        return 1

    print(f"ok: detection rule names agree in {len(present)} places ({const!r})")
    return 0


if __name__ == '__main__':
    sys.exit(main())
