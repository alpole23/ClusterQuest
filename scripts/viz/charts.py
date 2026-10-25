#!/usr/bin/env python3
"""Chart generation functions for BGC visualization."""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.constants import BGC_COLORS


def get_bgc_color(bgc_type):
    """Get color for a BGC type, handling hybrids and unknown types."""
    if bgc_type in BGC_COLORS:
        return BGC_COLORS[bgc_type]

    # Check if it's a hybrid (contains + or -)
    if '+' in bgc_type or '-' in bgc_type:
        # Try to match first component
        parts = bgc_type.replace('+', '-').split('-')
        for part in parts:
            if part in BGC_COLORS:
                return BGC_COLORS[part]
        return '#9b59b6'  # Default purple for hybrids

    # Check for partial matches
    bgc_lower = bgc_type.lower()
    for key, color in BGC_COLORS.items():
        if key.lower() in bgc_lower or bgc_lower in key.lower():
            return color

    return '#95a5a6'  # Default gray for unknown
