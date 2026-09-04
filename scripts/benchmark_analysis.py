#!/usr/bin/env python3
"""
ClusterQuest Pipeline Benchmark Analysis

Reads the Nextflow trace file from a completed run and generates a self-contained
HTML report covering:
  - Run summary (genome count, CPU hours, success rate)
  - Bottleneck analysis — Pareto chart of CPU hours by process
  - Per-process statistics (timing, memory, CPU efficiency)
  - Wall-time distributions for high-count processes
  - Scaling extrapolation to larger genome sets with cost estimates

Usage:
    python scripts/benchmark_analysis.py \\
        --trace results/pipeline_info/pipeline_trace.tsv \\
        --outdir results/benchmark \\
        [--n_genomes 2548] \\
        [--targets 10000,100000,500000,1000000,5000000] \\
        [--max_parallel 200] \\
        [--cpu_cost 0.048]
"""

import argparse
import base64
import io
import math
import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from statistics import mean, median, stdev

import matplotlib
matplotlib.use('Agg')
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

sys.path.insert(0, str(Path(__file__).parent))
from utils.parsers import (
    parse_trace_file, parse_duration, parse_memory,
    format_bytes, format_duration_str,
)

# ── Scaling classification ──────────────────────────────────────────────────

# BiG-SCAPE's all-vs-all comparison scales O(n²) in BGC count.
QUADRATIC_PROCESSES = {'BIGSCAPE'}

# Processes that run once per pipeline invocation regardless of genome count.
FIXED_PROCESSES = {
    'NCBI_DATASETS_DOWNLOAD', 'CREATE_NAME_MAP', 'EXTRACT_TAXONOMY',
    'GET_ANTISMASH_VERSION', 'DOWNLOAD_ANTISMASH_DBS', 'DOWNLOAD_TAXONKIT_DB',
    'DOWNLOAD_PFAM', 'DOWNLOAD_GTDBTK_DB', 'AGGREGATE_TAXONOMY',
    'TABULATE_REGIONS', 'COUNT_REGIONS', 'COLLECT_VERSIONS',
    'VISUALIZE_RESULTS', 'GCF_BIOSYNTHETIC_TREE', 'COUPLING_ENZYME_TREE',
    'CLUSTERING_STATS', 'EXTRACT_GCF_REPRESENTATIVES',
    'GTDBTK_CLASSIFY',
}

SCALE_COLORS = {
    'per_genome': '#2c7bb6',
    'quadratic':  '#d7191c',
    'fixed':      '#74add1',
}

SCALE_LABELS = {
    'per_genome': 'per-genome (linear)',
    'quadratic':  'quadratic (O·n²)',
    'fixed':      'fixed / once-per-run',
}

# ── Helpers ─────────────────────────────────────────────────────────────────

def extract_process_name(full_name):
    """'BGC_ANALYSIS:ANTISMASH_ANALYSIS:ANTISMASH (genome)' → 'ANTISMASH'"""
    name = re.sub(r'\s*\(.*\)$', '', full_name).strip()
    return name.split(':')[-1]


def extract_subworkflow(full_name):
    """'BGC_ANALYSIS:ANTISMASH_ANALYSIS:ANTISMASH (genome)' → 'BGC_ANALYSIS:ANTISMASH_ANALYSIS'"""
    name = re.sub(r'\s*\(.*\)$', '', full_name).strip()
    parts = name.split(':')
    return ':'.join(parts[:-1]) if len(parts) > 1 else 'root'


def classify_scaling(proc_name, n_tasks, n_genomes):
    """Return 'per_genome', 'quadratic', or 'fixed'."""
    if proc_name in QUADRATIC_PROCESSES:
        return 'quadratic'
    if proc_name in FIXED_PROCESSES:
        return 'fixed'
    if n_genomes and n_tasks >= n_genomes * 0.4:
        return 'per_genome'
    return 'fixed'


def fig_to_b64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


# ── Data loading ─────────────────────────────────────────────────────────────

def load_trace(trace_path):
    """Parse trace TSV and add derived numeric columns."""
    raw = parse_trace_file(trace_path)
    tasks = []
    for t in raw:
        if t.get('status') not in ('COMPLETED', 'FAILED', 'CACHED'):
            continue
        task = dict(t)
        task['process']       = extract_process_name(t['name'])
        task['subworkflow']   = extract_subworkflow(t['name'])
        task['realtime_s']    = parse_duration(t.get('realtime', '0'))
        task['peak_rss_b']    = parse_memory(t.get('peak_rss', '0'))
        task['memory_b']      = parse_memory(t.get('memory', '0'))
        task['cpus_int']      = int(t.get('cpus', 1) or 1)
        task['rchar_b']       = parse_memory(t.get('rchar', '0'))
        task['wchar_b']       = parse_memory(t.get('wchar', '0'))
        try:
            task['cpu_pct'] = float(str(t.get('%cpu', '0')).replace('%', ''))
        except ValueError:
            task['cpu_pct'] = 0.0
        # Actual CPU-seconds consumed
        task['cpu_s'] = task['realtime_s'] * task['cpus_int'] * task['cpu_pct'] / 100
        tasks.append(task)
    return tasks


def infer_n_genomes(tasks):
    """Infer genome count from per-genome process task counts."""
    for proc in ('ANTISMASH', 'CHECK_ANTISMASH_REUSE', 'COPY_ANTISMASH_RESULT', 'RENAME_GENOMES'):
        count = sum(1 for t in tasks if t['process'] == proc
                    and t['status'] in ('COMPLETED', 'FAILED'))
        if count > 10:
            return count
    return None


# ── Statistics ───────────────────────────────────────────────────────────────

def compute_process_stats(tasks):
    """Return {process_name: stats_dict} for all processes."""
    groups = defaultdict(list)
    for t in tasks:
        groups[t['process']].append(t)

    stats = {}
    for proc, ptasks in groups.items():
        done    = [t for t in ptasks if t['status'] == 'COMPLETED']
        failed  = [t for t in ptasks if t['status'] == 'FAILED']
        cached  = [t for t in ptasks if t['status'] == 'CACHED']
        times   = [t['realtime_s']  for t in done if t['realtime_s'] > 0]
        mems    = [t['peak_rss_b']  for t in done if t['peak_rss_b'] > 0]
        cpus    = [t['cpu_pct']     for t in done if t['cpu_pct']    > 0]
        cpu_s   = [t['cpu_s']       for t in done]

        def pct(lst, p):
            if not lst:
                return 0
            s = sorted(lst)
            idx = max(0, int(len(s) * p / 100) - 1)
            return s[idx]

        stats[proc] = {
            'n_completed':    len(done),
            'n_failed':       len(failed),
            'n_cached':       len(cached),
            'subworkflow':    (done or ptasks)[0]['subworkflow'],
            'cpus_allocated': (done or ptasks)[0]['cpus_int'],
            'mean_time_s':    mean(times)    if times else 0,
            'median_time_s':  median(times)  if times else 0,
            'p95_time_s':     pct(times, 95) if times else 0,
            'max_time_s':     max(times)     if times else 0,
            'total_time_s':   sum(times),
            'mean_peak_rss':  mean(mems)     if mems else 0,
            'max_peak_rss':   max(mems)      if mems else 0,
            'mean_cpu_pct':   mean(cpus)     if cpus else 0,
            'total_cpu_h':    sum(cpu_s) / 3600,
            'times':          times,
        }
    return stats


# ── Extrapolation ────────────────────────────────────────────────────────────

def extrapolate(proc_stats, n_current, targets, max_parallel):
    """
    Estimate wall time and CPU hours for each target genome count.
    Returns {n_target: {'wall_h': float, 'cpu_h': float, 'by_proc': {proc: wall_h}}}
    """
    results = {}
    for n_target in targets:
        by_proc = {}
        per_genome_walls = []
        bigscape_wall_h  = 0.0
        fixed_wall_h     = 0.0
        total_cpu_h      = 0.0

        for proc, s in proc_stats.items():
            scale = classify_scaling(proc, s['n_completed'], n_current)
            t_s   = s['total_time_s']
            cpu_h = s['total_cpu_h']

            if scale == 'per_genome' and s['mean_time_s'] > 0:
                # Wall time limited by parallelism
                parallel  = min(max_parallel, n_target)
                wall_s    = math.ceil(n_target / parallel) * s['mean_time_s']
                ratio     = n_target / n_current if n_current else 1
                wall_h    = wall_s / 3600
                est_cpu_h = cpu_h * ratio
                by_proc[proc]  = wall_h
                per_genome_walls.append(wall_h)
                total_cpu_h   += est_cpu_h

            elif scale == 'quadratic' and t_s > 0 and n_current:
                # O(n²) — BiG-SCAPE scales with BGC count; use genome count as proxy
                factor    = (n_target / n_current) ** 2
                wall_h    = t_s * factor / 3600
                by_proc[proc]  = wall_h
                bigscape_wall_h = wall_h
                total_cpu_h    += cpu_h * factor

            else:  # fixed
                by_proc[proc] = t_s / 3600
                fixed_wall_h += t_s / 3600
                total_cpu_h  += cpu_h

        # Critical-path wall time:
        # sequential fixed overhead + dominant per-genome bottleneck + BiG-SCAPE
        # (GTDB-Tk and antiSMASH overlap in practice, so use max not sum)
        per_genome_bottleneck = max(per_genome_walls, default=0)
        # Fixed sequential overhead: only a few small once-per-run processes are truly
        # on the critical path; cap at 2 hours as a conservative estimate for download + setup
        fixed_critical = min(fixed_wall_h, 2.0)
        wall_h_total = fixed_critical + per_genome_bottleneck + bigscape_wall_h

        results[n_target] = {
            'wall_h':  wall_h_total,
            'cpu_h':   total_cpu_h,
            'by_proc': by_proc,
        }
    return results


# ── Charts ───────────────────────────────────────────────────────────────────

def plot_pareto(proc_stats, n_genomes):
    """Horizontal Pareto chart of total CPU hours by process."""
    items = sorted(proc_stats.items(), key=lambda x: x[1]['total_cpu_h'], reverse=True)[:15]
    if not items:
        return None

    names     = [p for p, _ in items]
    cpu_hours = [s['total_cpu_h'] for _, s in items]
    total     = sum(cpu_hours) or 1
    cum_pct   = []
    running   = 0
    for h in cpu_hours:
        running += h
        cum_pct.append(running / total * 100)

    colors = [SCALE_COLORS[classify_scaling(p, proc_stats[p]['n_completed'], n_genomes)]
              for p in names]

    fig, ax1 = plt.subplots(figsize=(11, max(4, len(names) * 0.45 + 1)))
    bars = ax1.barh(range(len(names)), cpu_hours, color=colors, alpha=0.85, edgecolor='white')
    ax1.set_yticks(range(len(names)))
    ax1.set_yticklabels(names, fontsize=9)
    ax1.set_xlabel('Total CPU Hours', fontsize=10)
    ax1.invert_yaxis()

    max_h = max(cpu_hours) if cpu_hours else 1
    for bar, val in zip(bars, cpu_hours):
        ax1.text(bar.get_width() + max_h * 0.01,
                 bar.get_y() + bar.get_height() / 2,
                 f'{val:.1f}h', va='center', fontsize=8)

    ax2 = ax1.twiny()
    ax2.plot(cum_pct, range(len(names)), 'ko-', markersize=4, linewidth=1.5, alpha=0.7)
    ax2.axvline(80, color='#888', linestyle='--', linewidth=1, alpha=0.6)
    ax2.set_xlabel('Cumulative %', fontsize=10)
    ax2.set_xlim(0, 115)

    legend_handles = [mpatches.Patch(color=c, label=SCALE_LABELS[k])
                      for k, c in SCALE_COLORS.items()]
    ax1.legend(handles=legend_handles, loc='lower right', fontsize=8)
    ax1.set_title('CPU Hours by Process (Pareto)', fontsize=11, fontweight='bold', pad=12)
    fig.tight_layout()
    return fig_to_b64(fig)


def plot_timing_distribution(proc_stats):
    """Box plots of per-task wall times for variable processes."""
    candidates = [(p, s) for p, s in proc_stats.items() if len(s['times']) >= 10]
    candidates.sort(key=lambda x: x[1]['median_time_s'], reverse=True)
    candidates = candidates[:12]
    if not candidates:
        return None

    data_min = [[t / 60 for t in s['times']] for _, s in candidates]
    labels   = [p for p, _ in candidates]

    fig, ax = plt.subplots(figsize=(11, max(4, len(candidates) * 0.5 + 1.5)))
    bp = ax.boxplot(data_min, vert=False, patch_artist=True,
                    flierprops=dict(marker='o', markersize=2, alpha=0.35, color='#666'))
    for patch in bp['boxes']:
        patch.set_facecolor('#abd9e9')
        patch.set_alpha(0.75)

    ax.set_yticks(range(1, len(labels) + 1))
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel('Wall Time (minutes)', fontsize=10)
    ax.set_title('Per-Task Wall Time Distribution', fontsize=11, fontweight='bold')
    ax.grid(axis='x', alpha=0.3)
    fig.tight_layout()
    return fig_to_b64(fig)


def plot_memory(proc_stats):
    """Horizontal bar chart of peak RSS memory by process."""
    candidates = [(p, s) for p, s in proc_stats.items() if s['max_peak_rss'] > 0]
    candidates.sort(key=lambda x: x[1]['max_peak_rss'], reverse=True)
    candidates = candidates[:15]
    if not candidates:
        return None

    names    = [p for p, _ in candidates]
    mean_gb  = [s['mean_peak_rss'] / 1024**3 for _, s in candidates]
    max_gb   = [s['max_peak_rss']  / 1024**3 for _, s in candidates]

    fig, ax = plt.subplots(figsize=(11, max(4, len(names) * 0.45 + 1)))
    y = range(len(names))
    ax.barh(y, max_gb,  color='#fdae61', alpha=0.6, label='Peak RSS (max)')
    ax.barh(y, mean_gb, color='#d7191c', alpha=0.85, label='Peak RSS (mean)')
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=9)
    ax.set_xlabel('Memory (GB)', fontsize=10)
    ax.invert_yaxis()
    ax.legend(fontsize=9)
    ax.set_title('Peak Memory Usage by Process', fontsize=11, fontweight='bold')
    ax.grid(axis='x', alpha=0.3)
    fig.tight_layout()
    return fig_to_b64(fig)


def plot_extrapolation(proc_stats, n_current, targets, max_parallel, cpu_cost):
    """Two-panel: wall time and CPU hours vs genome count (log-log)."""
    # Dense targets for smooth curves, always including n_current
    dense = sorted(set(targets + [n_current]))
    extrap = extrapolate(proc_stats, n_current, dense, max_parallel)

    wall_h   = [extrap[n]['wall_h']  for n in dense]
    cpu_h    = [extrap[n]['cpu_h']   for n in dense]
    cur_wall = extrap[n_current]['wall_h']
    cur_cpu  = extrap[n_current]['cpu_h']

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # — Wall time panel —
    ax1.plot(dense, wall_h, 'b-o', lw=2, ms=4, label='Estimated wall time')
    ax1.scatter([n_current], [cur_wall], color='red', zorder=6, s=90,
                label=f'This run ({n_current:,} genomes)')
    ax1.axvline(n_current, color='#aaa', lw=1, linestyle='--')
    ax1.set_xscale('log'); ax1.set_yscale('log')
    ax1.set_xlabel('Genomes', fontsize=10)
    ax1.set_ylabel('Wall Time (hours)', fontsize=10)
    ax1.set_title('Estimated Wall Time', fontsize=11, fontweight='bold')
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.25)
    ax1.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f'{int(x):,}'))

    # — CPU hours + cost panel —
    ax2.plot(dense, cpu_h, 'g-o', lw=2, ms=4, label='CPU hours')
    ax2.scatter([n_current], [cur_cpu], color='red', zorder=6, s=90)
    ax2.axvline(n_current, color='#aaa', lw=1, linestyle='--')
    if cpu_cost > 0:
        ax2b = ax2.twinx()
        costs = [h * cpu_cost for h in cpu_h]
        ax2b.plot(dense, costs, 'r--s', lw=1.5, ms=4, alpha=0.7,
                  label=f'Cost (${cpu_cost}/CPU-hr)')
        ax2b.set_ylabel(f'Est. Cost (USD @ ${cpu_cost}/CPU-hr)', fontsize=9, color='#c0392b')
        ax2b.tick_params(axis='y', labelcolor='#c0392b')
        ax2b.set_yscale('log')
        ax2b.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f'${x:,.0f}'))
        ax2b.legend(fontsize=9, loc='upper left')
    ax2.set_xscale('log'); ax2.set_yscale('log')
    ax2.set_xlabel('Genomes', fontsize=10)
    ax2.set_ylabel('Total CPU Hours', fontsize=10)
    ax2.set_title('CPU Hours & Cost Estimate', fontsize=11, fontweight='bold')
    ax2.legend(fontsize=9, loc='lower right')
    ax2.grid(True, alpha=0.25)
    ax2.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f'{int(x):,}'))

    fig.tight_layout()
    return fig_to_b64(fig)


# ── HTML report ───────────────────────────────────────────────────────────────

def generate_report(tasks, proc_stats, n_genomes, targets, max_parallel,
                    cpu_cost, outdir, trace_path):

    done   = [t for t in tasks if t['status'] == 'COMPLETED']
    failed = [t for t in tasks if t['status'] == 'FAILED']
    cached = [t for t in tasks if t['status'] == 'CACHED']

    total_cpu_h  = sum(s['total_cpu_h'] for s in proc_stats.values())
    success_rate = (len(done) / (len(done) + len(failed)) * 100
                    if (done or failed) else 100.0)

    # Charts
    pareto_b64  = plot_pareto(proc_stats, n_genomes)
    timing_b64  = plot_timing_distribution(proc_stats)
    memory_b64  = plot_memory(proc_stats)
    extrap_b64  = plot_extrapolation(proc_stats, n_genomes, targets, max_parallel, cpu_cost)

    # Extrapolation table
    extrap_data  = extrapolate(proc_stats, n_genomes, targets, max_parallel)

    # ── Process stats table rows ──
    sorted_procs = sorted(proc_stats.items(), key=lambda x: x[1]['total_cpu_h'], reverse=True)
    proc_rows = ''
    for proc, s in sorted_procs:
        scale = classify_scaling(proc, s['n_completed'], n_genomes)
        color = SCALE_COLORS[scale]
        label = SCALE_LABELS[scale]
        fail_badge = (f' <span style="color:#e74c3c;font-size:0.82em;">({s["n_failed"]} failed)</span>'
                      if s['n_failed'] else '')
        proc_rows += f'''
            <tr>
                <td style="font-family:monospace;font-size:0.85em;">{proc}{fail_badge}</td>
                <td style="color:{color};font-size:0.82em;white-space:nowrap;">{label}</td>
                <td style="text-align:right;">{s["n_completed"]:,}</td>
                <td style="text-align:right;">{format_duration_str(s["mean_time_s"])}</td>
                <td style="text-align:right;">{format_duration_str(s["p95_time_s"])}</td>
                <td style="text-align:right;">{format_duration_str(s["max_time_s"])}</td>
                <td style="text-align:right;">{format_bytes(s["mean_peak_rss"])}</td>
                <td style="text-align:right;">{s["mean_cpu_pct"]:.0f}%</td>
                <td style="text-align:right;">{s["total_cpu_h"]:.1f}h</td>
            </tr>'''

    # ── Extrapolation table rows ──
    extrap_rows  = ''
    cost_header  = f'<th>Est. Cost (@ ${cpu_cost}/CPU-hr)</th>' if cpu_cost > 0 else ''
    for n_target in targets:
        row      = extrap_data[n_target]
        wall_str = format_duration_str(row['wall_h'] * 3600)
        cpu_str  = f'{row["cpu_h"]:,.0f}'
        cost_str = (f'<td style="text-align:right;">${row["cpu_h"] * cpu_cost:,.0f}</td>'
                    if cpu_cost > 0 else '')
        highlight = (' style="background:#fffbea;font-weight:600;"'
                     if n_target == n_genomes else '')
        extrap_rows += f'''
            <tr{highlight}>
                <td style="text-align:right;">{n_target:,}</td>
                <td style="text-align:right;">{wall_str}</td>
                <td style="text-align:right;">{cpu_str}</td>
                {cost_str}
            </tr>'''

    # ── Optional sections ──
    timing_html = (f'<img src="data:image/png;base64,{timing_b64}" style="max-width:100%;">'
                   if timing_b64 else
                   '<p style="color:#888;">No per-genome processes with ≥10 tasks for distribution analysis.</p>')
    memory_html = (f'<img src="data:image/png;base64,{memory_b64}" style="max-width:100%;">'
                   if memory_b64 else '')
    bigscape_warning = (
        '<div class="warning-box">⚠️  <strong>BiG-SCAPE scales quadratically</strong> with BGC count. '
        'At large genome counts this will dominate runtime and may be impractical without sharding '
        'or replacing with a linear-time alternative (e.g. MASH-based pre-clustering).</div>'
        if 'BIGSCAPE' in proc_stats else ''
    )

    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>ClusterQuest Benchmark Report</title>
    <style>
        * {{ box-sizing: border-box; }}
        body {{ font-family: "Segoe UI", Tahoma, sans-serif; margin: 0; padding: 24px 44px;
               background: #f5f7fa; color: #2d3436; }}
        h1   {{ color: #2c3e50; margin-bottom: 4px; font-size: 1.8em; }}
        h2   {{ color: #2c5aa0; border-bottom: 2px solid #2c5aa0; padding-bottom: 6px;
               margin-top: 44px; font-size: 1.2em; }}
        .subtitle {{ color: #888; font-size: 0.92em; margin-bottom: 28px; }}
        .cards {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
                 gap: 14px; margin: 18px 0 10px; }}
        .card {{ background: white; border-radius: 10px; padding: 18px 14px;
                box-shadow: 0 2px 8px rgba(0,0,0,0.07); text-align: center; }}
        .card .val {{ font-size: 1.9em; font-weight: 700; color: #2c5aa0; }}
        .card .lbl {{ color: #666; font-size: 0.83em; margin-top: 5px; }}
        .card .sub {{ color: #aaa; font-size: 0.76em; margin-top: 2px; }}
        .plot {{ background: white; border-radius: 10px; padding: 20px;
                box-shadow: 0 2px 8px rgba(0,0,0,0.07); margin: 18px 0; }}
        table {{ width: 100%; border-collapse: collapse; background: white;
                border-radius: 10px; overflow: hidden;
                box-shadow: 0 2px 8px rgba(0,0,0,0.07); margin: 16px 0; }}
        th {{ background: #2c5aa0; color: white; padding: 10px 12px;
             text-align: left; font-size: 0.86em; white-space: nowrap; }}
        td {{ padding: 7px 12px; border-bottom: 1px solid #f0f0f0; font-size: 0.86em; }}
        tr:last-child td {{ border-bottom: none; }}
        tr:hover {{ background: #f8f9ff; }}
        .info-box    {{ background: #e8f4fd; border-left: 4px solid #2c5aa0;
                       padding: 10px 14px; border-radius: 4px; margin: 12px 0;
                       font-size: 0.9em; line-height: 1.5; }}
        .warning-box {{ background: #fff8e1; border-left: 4px solid #f39c12;
                       padding: 10px 14px; border-radius: 4px; margin: 12px 0;
                       font-size: 0.9em; line-height: 1.5; }}
        code {{ background: #f0f0f0; padding: 1px 5px; border-radius: 3px;
               font-size: 0.9em; }}
    </style>
</head>
<body>
    <h1>ClusterQuest Pipeline Benchmark Report</h1>
    <p class="subtitle">
        Trace: <code>{trace_path}</code>
        &nbsp;|&nbsp; Generated: <span id="ts"></span>
    </p>
    <script>document.getElementById('ts').textContent = new Date().toLocaleString();</script>

    <!-- ── Run Summary ── -->
    <h2>Run Summary</h2>
    <div class="cards">
        <div class="card">
            <div class="val">{n_genomes:,}</div>
            <div class="lbl">Genomes Analyzed</div>
        </div>
        <div class="card">
            <div class="val">{len(done):,}</div>
            <div class="lbl">Tasks Completed</div>
            <div class="sub">{len(cached)} cached &nbsp;·&nbsp; {len(failed)} failed</div>
        </div>
        <div class="card">
            <div class="val">{success_rate:.1f}%</div>
            <div class="lbl">Success Rate</div>
        </div>
        <div class="card">
            <div class="val">{total_cpu_h:,.0f}h</div>
            <div class="lbl">Total CPU Hours</div>
        </div>
        <div class="card">
            <div class="val">{total_cpu_h / n_genomes:.2f}h</div>
            <div class="lbl">CPU Hours / Genome</div>
        </div>
        <div class="card">
            <div class="val">{max_parallel}</div>
            <div class="lbl">Max Parallel Jobs</div>
            <div class="sub">used for extrapolation</div>
        </div>
    </div>

    <!-- ── Bottleneck Analysis ── -->
    <h2>Bottleneck Analysis</h2>
    <div class="info-box">
        Processes are classified by scaling behavior:
        <strong style="color:{SCALE_COLORS["per_genome"]}">per-genome (linear)</strong> — scales linearly, parallelizable across SLURM jobs;
        <strong style="color:{SCALE_COLORS["quadratic"]}">quadratic</strong> — BiG-SCAPE all-vs-all comparison, scales O(n²) in BGC count;
        <strong style="color:{SCALE_COLORS["fixed"]}">fixed / once-per-run</strong> — constant cost regardless of genome count.
        The 80% cumulative line (dashed) marks where engineering effort has the most leverage.
    </div>
    <div class="plot">
        <img src="data:image/png;base64,{pareto_b64}" style="max-width:100%;">
    </div>

    <!-- ── Per-Process Statistics ── -->
    <h2>Per-Process Statistics</h2>
    <table>
        <thead>
            <tr>
                <th>Process</th>
                <th>Scaling</th>
                <th>Tasks</th>
                <th>Mean Time</th>
                <th>p95 Time</th>
                <th>Max Time</th>
                <th>Mean Peak RSS</th>
                <th>CPU Efficiency</th>
                <th>Total CPU Hours</th>
            </tr>
        </thead>
        <tbody>{proc_rows}</tbody>
    </table>

    <!-- ── Wall Time Distributions ── -->
    <h2>Wall Time Distributions (Per-Genome Processes)</h2>
    <div class="plot">{timing_html}</div>

    <!-- ── Memory Usage ── -->
    <h2>Memory Usage</h2>
    <div class="plot">{memory_html}</div>

    <!-- ── Scaling Extrapolation ── -->
    <h2>Scaling Extrapolation</h2>
    <div class="info-box">
        Wall time estimates assume <strong>{max_parallel} maximum parallel SLURM jobs</strong>.
        Per-genome wall time = ⌈n / {max_parallel}⌉ × mean task time.
        BiG-SCAPE uses O(n²) scaling with genome count as a proxy for BGC count.
        These are optimistic lower bounds — they assume perfect parallelism and no I/O contention.
    </div>
    {bigscape_warning}
    <div class="plot">
        <img src="data:image/png;base64,{extrap_b64}" style="max-width:100%;">
    </div>
    <table>
        <thead>
            <tr>
                <th>Genomes</th>
                <th>Est. Wall Time</th>
                <th>Est. Total CPU Hours</th>
                {cost_header}
            </tr>
        </thead>
        <tbody>{extrap_rows}</tbody>
    </table>
    <p style="color:#aaa;font-size:0.82em;margin-top:6px;">
        Highlighted row = this run.
        Cost uses ${cpu_cost}/CPU-hr (approximate AWS c5 on-demand).
        Actual cloud costs vary by instance type, spot pricing, and storage I/O.
    </p>
</body>
</html>'''

    os.makedirs(outdir, exist_ok=True)
    out_path = os.path.join(outdir, 'benchmark_report.html')
    with open(out_path, 'w') as f:
        f.write(html)
    return out_path


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='ClusterQuest Pipeline Benchmark Analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('--trace',        required=True,
                        help='Nextflow trace TSV (results/pipeline_info/pipeline_trace.tsv)')
    parser.add_argument('--outdir',       default='results/benchmark',
                        help='Output directory (default: results/benchmark)')
    parser.add_argument('--n_genomes',    type=int, default=None,
                        help='Genome count (inferred from trace if omitted)')
    parser.add_argument('--targets',
                        default='10000,50000,100000,500000,1000000,5000000',
                        help='Comma-separated genome counts for extrapolation table/chart')
    parser.add_argument('--max_parallel', type=int, default=200,
                        help='Max concurrent SLURM jobs (default: 200)')
    parser.add_argument('--cpu_cost',     type=float, default=0.048,
                        help='USD per CPU-hour for cost estimates (0 to omit, default: 0.048)')
    args = parser.parse_args()

    target_list = sorted(set(int(x.strip()) for x in args.targets.split(',')))

    print(f'Loading trace: {args.trace}')
    tasks = load_trace(args.trace)
    if not tasks:
        print('ERROR: No completed/failed/cached tasks found in trace file.')
        sys.exit(1)
    print(f'  {len(tasks)} tasks loaded')

    n_genomes = args.n_genomes or infer_n_genomes(tasks)
    if not n_genomes:
        print('ERROR: Could not infer genome count. Pass --n_genomes explicitly.')
        sys.exit(1)
    print(f'  Genome count: {n_genomes:,}')

    proc_stats = compute_process_stats(tasks)
    print(f'  {len(proc_stats)} distinct processes\n')

    print(f'{"Process":<45} {"CPU Hours":>10}  {"Mean Time":>10}  {"Tasks":>6}')
    print('-' * 75)
    for proc, s in sorted(proc_stats.items(), key=lambda x: x[1]['total_cpu_h'], reverse=True)[:10]:
        print(f'  {proc:<43} {s["total_cpu_h"]:>10.1f}  '
              f'{format_duration_str(s["mean_time_s"]):>10}  {s["n_completed"]:>6,}')

    print(f'\nGenerating report...')
    out_path = generate_report(
        tasks, proc_stats, n_genomes, target_list,
        args.max_parallel, args.cpu_cost, args.outdir, args.trace,
    )
    print(f'Saved: {out_path}')
    print('Done.')


if __name__ == '__main__':
    main()
