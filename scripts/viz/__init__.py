"""Visualization modules for BGC analysis pipeline."""

from .charts import get_bgc_color, plot_kcb_identification_chart
from .tables import (
    get_genome_count,
    generate_genome_table_html,
    calculate_summary_statistics,
    create_bgc_distribution_table
)
from .clustering import (
    generate_bigscape_stats_html,
    generate_gcf_visualization_html
)
from .taxonomy import generate_taxonomy_tree_html
from .distribution import (
    extract_assembly_id_from_genome_name,
    build_gcf_taxonomy_distribution,
    generate_bgc_distribution_html
)
from .genome_pages import create_genome_metadata_pages
from .rarefaction import generate_rarefaction_curve
from .report_assets import REPORT_CSS, REPORT_JS
from .report_sections import build_coupling_table_rows
# Trace/resource rendering lives in utils.trace (single implementation); re-exported
# here so `viz.generate_resource_usage_html` keeps working.
from utils.trace import (
    aggregate_trace_by_process,
    generate_gantt_chart_html,
    generate_resource_usage_html
)

__all__ = [
    # charts
    'get_bgc_color',
    'plot_kcb_identification_chart',
    # tree_viz
    # tables
    'get_genome_count',
    'generate_genome_table_html',
    'calculate_summary_statistics',
    'create_bgc_distribution_table',
    # clustering
    'generate_bigscape_stats_html',
    'generate_gcf_visualization_html',
    # taxonomy
    'generate_taxonomy_tree_html',
    # resources
    'aggregate_trace_by_process',
    'generate_gantt_chart_html',
    'generate_resource_usage_html',
    # distribution / genome pages / rarefaction
    'extract_assembly_id_from_genome_name',
    'build_gcf_taxonomy_distribution',
    'generate_bgc_distribution_html',
    'create_genome_metadata_pages',
    'generate_rarefaction_curve',
    # report assembly
    'REPORT_CSS',
    'REPORT_JS',
    'build_coupling_table_rows',
]
