"""Static CSS and JavaScript for the BGC HTML report.

Plain module-level strings (not f-strings) so braces in the CSS/JS need no escaping.
generate_html_report() inlines them into the page, keeping the report self-contained.
"""


REPORT_CSS = """\
        * { box-sizing: border-box; }
        body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 0; padding: 20px 40px; background: #f8f9fa; color: #333; }
        h1 { color: #333; text-align: center; margin-bottom: 5px; }
        h2 { color: #333; border-bottom: 2px solid #5b8ac5; padding-bottom: 10px; margin-top: 30px; }
        h3 { color: #444; margin-top: 25px; }
        .subtitle { text-align: center; color: #666; font-size: 1.1em; margin-bottom: 20px; }

        /* Tab Styles */
        .tabs {
            margin-top: 20px;
        }
        .tabs input[type="radio"] {
            display: none;
        }
        .tabs label {
            display: inline-block;
            padding: 12px 24px;
            background: #e9ecef;
            color: #666;
            cursor: pointer;
            border-radius: 8px 8px 0 0;
            margin-right: 4px;
            font-weight: 500;
            transition: all 0.2s;
            border: 1px solid #ddd;
            border-bottom: none;
            white-space: nowrap;
        }
        .tabs label:hover {
            background: #dee2e6;
            color: #333;
        }
        .tabs input[type="radio"]:checked + label {
            background: white;
            color: #2c5aa0;
            border-color: #5b8ac5;
            font-weight: 600;
        }
        .tab-content {
            display: none;
            background: white;
            padding: 30px;
            border: 1px solid #ddd;
            border-radius: 0 8px 8px 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.08);
            min-height: 400px;
        }
        #tab1:checked ~ #content1,
        #tab2:checked ~ #content2,
        #tab3:checked ~ #content3,
        #tab4:checked ~ #content4,
        #tab5:checked ~ #content5,
        #tab6:checked ~ #content6,
        #tab7:checked ~ #content7 {
            display: block;
        }

        /* Collapsible details block (Pipeline Info in Overview) */
        details.pipeline-info {
            margin-top: 30px;
            border: 1px solid #dee2e6;
            border-radius: 8px;
            overflow: hidden;
        }
        details.pipeline-info > summary {
            padding: 12px 18px;
            background: #f8f9fa;
            cursor: pointer;
            font-weight: 600;
            font-size: 1em;
            color: #444;
            list-style: none;
            display: flex;
            align-items: center;
            gap: 8px;
        }
        details.pipeline-info > summary::-webkit-details-marker { display: none; }
        details.pipeline-info > summary::before {
            content: '▶';
            font-size: 0.75em;
            transition: transform 0.2s;
            display: inline-block;
        }
        details.pipeline-info[open] > summary::before { transform: rotate(90deg); }
        details.pipeline-info > .details-body {
            padding: 20px 24px;
        }

        /* Section divider within merged tabs */
        .tab-section-divider {
            border: none;
            border-top: 2px solid #e9ecef;
            margin: 32px 0 28px;
        }

        /* Stats Dashboard */
        .stats-container {
            display: grid;
            grid-template-columns: repeat(8, 1fr);
            gap: 8px;
            margin: 15px 0;
        }
        @media (max-width: 1200px) {
            .stats-container {
                grid-template-columns: repeat(4, 1fr);
            }
        }
        @media (max-width: 768px) {
            .stats-container {
                grid-template-columns: repeat(2, 1fr);
            }
        }
        .stat-box {
            background: #f8f9fa;
            padding: 12px 8px;
            border-radius: 6px;
            border: 1px solid #e9ecef;
            text-align: center;
        }
        .stat-box.highlight {
            /* Same as regular stat-box */
        }
        .stat-value {
            font-size: 1.4em;
            font-weight: bold;
            color: #2c5aa0;
        }
        .stat-label {
            color: #666;
            margin-top: 3px;
            font-size: 0.8em;
        }

        /* Plots */
        .plot {
            margin: 20px 0;
            text-align: center;
            background: #fafafa;
            padding: 20px;
            border-radius: 8px;
            border: 1px solid #eee;
        }
        .plot img {
            max-width: 100%;
            height: auto;
            border-radius: 4px;
        }
        .plot-row {
            display: flex;
            gap: 20px;
            flex-wrap: wrap;
            justify-content: center;
        }
        .plot-row .plot {
            flex: 1;
            min-width: 300px;
            max-width: 500px;
        }

        /* Tables */
        table {
            border-collapse: collapse;
            width: 100%;
            background: white;
            font-size: 0.9em;
        }
        th, td {
            border: 1px solid #ddd;
            padding: 10px 12px;
            text-align: left;
        }
        th {
            background: #2c5aa0;
            color: white;
            font-weight: 600;
            position: sticky;
            top: 0;
        }
        tr:nth-child(even) { background-color: #f8f9fa; }
        tr:hover { background-color: #e3f2fd; }
        .table-container {
            max-height: 500px;
            overflow-y: auto;
            border: 1px solid #ddd;
            border-radius: 8px;
        }

        /* Search Box */
        .search-box {
            margin-bottom: 15px;
        }
        .search-box input {
            width: 100%;
            max-width: 400px;
            padding: 10px 15px;
            font-size: 1em;
            border: 2px solid #ddd;
            border-radius: 6px;
            outline: none;
            transition: border-color 0.2s;
        }
        .search-box input:focus {
            border-color: #5b8ac5;
        }

        /* Info Box */
        .info-box {
            background: #e8f4f8;
            border-left: 4px solid #2c5aa0;
            padding: 20px;
            margin: 20px 0;
            border-radius: 0 8px 8px 0;
        }
        .info-box.warning {
            background: #fff8e1;
            border-left-color: #ffc107;
        }

        /* Links & Buttons */
        a { color: #2c5aa0; text-decoration: none; }
        a:hover { color: #5b8ac5; text-decoration: underline; }
        .btn {
            display: inline-block;
            padding: 10px 20px;
            background: #2c5aa0;
            color: white;
            border-radius: 6px;
            font-weight: 500;
            text-decoration: none;
            transition: background 0.2s;
        }
        .btn:hover { background: #5b8ac5; color: white; text-decoration: none; }
        code {
            background: #f4f4f4;
            padding: 3px 8px;
            border-radius: 4px;
            font-family: 'Consolas', monospace;
            font-size: 0.9em;
            border: 1px solid #e0e0e0;
        }

        /* Clustering sections */
        .clustering-section {
            margin-bottom: 30px;
            padding-bottom: 30px;
            border-bottom: 1px solid #eee;
        }
        .clustering-section:last-child {
            border-bottom: none;
        }
"""

REPORT_JS = """\
            const input = document.getElementById('genomeSearch');
            const filter = searchNorm(input.value);
            const tbody = document.getElementById('genomeTableBody');
            const rows = tbody.getElementsByTagName('tr');

            for (let i = 0; i < rows.length; i++) {
                const cells = rows[i].getElementsByTagName('td');
                let found = false;
                for (let j = 0; j < cells.length; j++) {
                    if (searchMatches(cells[j].textContent, filter)) {
                        found = true;
                        break;
                    }
                }
                rows[i].style.display = found ? '' : 'none';
            }
        }

        // Genome names use underscores (Pantoea_ananatis_LMG_5342) but people type
        // spaces. Collapse both to a single space on each side so "LMG 5342",
        // "LMG_5342" and "lmg  5342" all match the same row.
        // Filters the Genomes tab. This was previously missing entirely: the search
        // box called filterGenomes(), which was never defined, so typing there threw
        // a ReferenceError and filtered nothing.
        // ---- Genome table: rendered from JSON, not from a 1,735-row DOM ----------
        // The fully-rendered table was 612 KB, the largest single element in the report,
        // parsed and painted on load though almost nobody scrolls past the first screen.
        // Rows now live in a JSON island and are rendered on demand; search runs over the
        // array, so it still covers every genome rather than only the rendered ones.
        const GENOME_RENDER_CAP = 500;
        let _genomeRows = null;

        function genomeRows() {
            if (_genomeRows === null) {
                const el = document.getElementById('genomeData');
                try { _genomeRows = el ? JSON.parse(el.textContent) : []; }
                catch (e) { _genomeRows = []; }
            }
            return _genomeRows;
        }

        function esc(v) {
            return String(v === null || v === undefined ? '' : v)
                .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
                .replace(/"/g, '&quot;');
        }

        // Mirrors _genome_row_html() in viz/tables.py — change both together.
        function renderGenomeRows(rows, cap) {
            const tbody = document.getElementById('genomeTableBody');
            if (!tbody) return 0;
            const limit = cap === null ? rows.length : Math.min(rows.length, cap);
            const html = [];
            for (let i = 0; i < limit; i++) {
                const r = rows[i];
                html.push('<tr><td><a href="genomes/' + esc(r[0]) + '.html">' + esc(r[0]) + '</a></td>' +
                          '<td>' + esc(r[1]) + '</td>' +
                          '<td title="' + esc(r[2]) + '">' + esc(r[2]) + '</td>' +
                          '<td>' + esc(r[3]) + '</td>' +
                          '<td>' + esc(r[4]) + '</td>' +
                          '<td>' + esc(r[5]) + '</td></tr>');
            }
            tbody.innerHTML = html.join('');
            return limit;
        }

        function setGenomeStatus(shown, matched, total, filtered) {
            const el = document.getElementById('genomeTableStatus');
            if (!el) return;
            let msg;
            if (filtered) {
                msg = 'Showing ' + shown + ' of ' + matched + ' matching genomes'
                    + (matched < total ? ' (' + total + ' total)' : '') + '.';
            } else {
                msg = 'Showing ' + shown + ' of ' + total + ' genomes.';
            }
            const more = shown < (filtered ? matched : total);
            el.innerHTML = msg + (more
                ? ' <button type="button" onclick="showAllGenomes()" style="background:none;border:none;'
                  + 'color:#2c5aa0;cursor:pointer;padding:0;font:inherit;text-decoration:underline;">Show all</button>'
                : '');
        }

        function showAllGenomes() {
            const input = document.getElementById('genomeSearch');
            const filter = input ? searchNorm(input.value) : '';
            const rows = genomeRows();
            const matched = filter
                ? rows.filter(function (r) { return r.some(function (c) { return searchMatches(c, filter); }); })
                : rows;
            const shown = renderGenomeRows(matched, null);
            setGenomeStatus(shown, matched.length, rows.length, !!filter);
        }

        function _filterGenomes() {
            const input = document.getElementById('genomeSearch');
            if (!input) return;
            const filter = searchNorm(input.value);
            const rows = genomeRows();
            if (!rows.length) return;          // no JSON island: leave the server-rendered rows alone
            const matched = filter
                ? rows.filter(function (r) { return r.some(function (c) { return searchMatches(c, filter); }); })
                : rows;
            const shown = renderGenomeRows(matched, GENOME_RENDER_CAP);
            setGenomeStatus(shown, matched.length, rows.length, !!filter);
        }

        // The genome table can hold thousands of rows (1,735 on the Pantoea genus) and
        // each keystroke scans every cell — ~10k reads. Debouncing keeps typing
        // responsive; the public names are unchanged so the onkeyup markup still works.
        const _filterTimers = {};
        function _debounce(key, fn, ms) {
            clearTimeout(_filterTimers[key]);
            _filterTimers[key] = setTimeout(fn, ms);
        }
        function filterGenomes()   { _debounce('genomes', _filterGenomes,   150); }
        function filterNovelBGCs() { _debounce('novel',   _filterNovelBGCs, 150); }
        function filterKCBHits()   { _debounce('kcb',     _filterKCBHits,   150); }

        function searchNorm(text) {
            return text.toLowerCase().replace(/[_\\s]+/g, ' ').trim();
        }

        // Whether a cell matches the query. Plain substring, except that a query
        // ending in a digit will not match inside a longer number: "GCF-1" finds
        // GCF-1 but not GCF-10..GCF-13, and "1" does not hit the member count 215.
        // Prefix search on words is unaffected — "Panto" still finds Pantoea.
        function searchMatches(cellText, needle) {
            if (!needle) return true;
            const hay = searchNorm(cellText);
            if (!/\\d$/.test(needle)) return hay.indexOf(needle) > -1;
            let from = 0, i;
            while ((i = hay.indexOf(needle, from)) > -1) {
                const before = i > 0 ? hay.charAt(i - 1) : '';
                const after  = hay.charAt(i + needle.length);
                if (!/\\d/.test(before) && !/\\d/.test(after)) return true;
                from = i + 1;
            }
            return false;
        }

        function _filterNovelBGCs() {
            const input = document.getElementById('novelSearch');
            if (!input) return;
            const filter = searchNorm(input.value);
            const tbody = document.getElementById('novelTableBody');
            if (!tbody) return;   // table absent, e.g. the empty-hits KCB tab
            const rows = tbody.getElementsByTagName('tr');

            for (let i = 0; i < rows.length; i++) {
                const cells = rows[i].getElementsByTagName('td');
                let found = false;
                for (let j = 0; j < cells.length; j++) {
                    if (searchMatches(cells[j].textContent, filter)) {
                        found = true;
                        break;
                    }
                }
                rows[i].style.display = found ? '' : 'none';
            }
        }

        function _filterKCBHits() {
            const input = document.getElementById('kcbSearch');
            if (!input) return;
            const filter = searchNorm(input.value);
            const tbody = document.getElementById('kcbTableBody');
            if (!tbody) return;   // table absent, e.g. the empty-hits KCB tab
            const rows = tbody.getElementsByTagName('tr');

            for (let i = 0; i < rows.length; i++) {
                const cells = rows[i].getElementsByTagName('td');
                let found = false;
                for (let j = 0; j < cells.length; j++) {
                    if (searchMatches(cells[j].textContent, filter)) {
                        found = true;
                        break;
                    }
                }
                rows[i].style.display = found ? '' : 'none';
            }
        }

        function toggleNode(nodeId) {
            const element = document.getElementById(nodeId);
            const header = element.previousElementSibling;
            const icon = header.querySelector('.toggle-icon');
            if (element.style.display === 'none') {
                element.style.display = 'block';
                icon.innerHTML = '&#9660;';
            } else {
                element.style.display = 'none';
                icon.innerHTML = '&#9654;';
            }
        }

        // Initialize all tree nodes as collapsed except root
        document.addEventListener('DOMContentLoaded', function() {
            const nodeChildren = document.querySelectorAll('.node-children');
            nodeChildren.forEach(function(node, index) {
                if (index > 0) {
                    node.style.display = 'none';
                    const header = node.previousElementSibling;
                    if (header) {
                        const icon = header.querySelector('.toggle-icon');
                        if (icon) icon.innerHTML = '&#9654;';
                    }
                }
            });
        });
"""
