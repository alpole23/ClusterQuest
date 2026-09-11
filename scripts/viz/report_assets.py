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

        /* ---- Sidebar navigation -------------------------------------------
           Still pure CSS: the radio inputs drive `#tabN:checked ~ #contentN`,
           which does not care whether the labels sit above or beside the panes.
           Group headings are static labels, so sub-sections need no mechanism of
           their own — they are simply more radios under a heading. */
        .tabs {
            margin-top: 20px;
            display: grid;
            grid-template-columns: 224px minmax(0, 1fr);
            gap: 0 26px;
            align-items: start;
        }
        .tabs input[type="radio"] { display: none; }
        .tabs label {
            grid-column: 1;
            display: block;
            padding: 9px 14px;
            color: #555;
            cursor: pointer;
            border-left: 3px solid transparent;
            border-radius: 0 5px 5px 0;
            font-size: 0.93em;
            font-weight: 500;
            transition: background 0.15s, color 0.15s, border-color 0.15s;
        }
        .tabs label.sub { padding-left: 28px; font-size: 0.9em; }
        .tabs label:hover { background: #eef2f6; color: #1f3b5f; }
        .tabs input[type="radio"]:checked + label {
            background: #e8eff7;
            border-left-color: #2c5aa0;
            color: #2c5aa0;
            font-weight: 600;
        }
        .nav-group {
            grid-column: 1;
            padding: 16px 14px 5px;
            font-size: 0.7em;
            font-weight: 700;
            letter-spacing: 0.09em;
            text-transform: uppercase;
            color: #97a3b0;
        }
        .tab-content {
            grid-column: 2;
            grid-row: 1 / 99;
            display: none;
            background: white;
            padding: 28px 30px;
            border: 1px solid #ddd;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.08);
            min-height: 400px;
        }
        #tab1:checked ~ #content1,
        #tab2:checked ~ #content2,
        #tab3:checked ~ #content3,
        #tab4:checked ~ #content4,
        #tab5:checked ~ #content5,
        #tab6:checked ~ #content6,
        #tab7:checked ~ #content7,
        #tab8:checked ~ #content8 {
            display: block;
        }
        /* Below this width a 224px rail costs more than it gives, so the nav
           stacks above the pane and the labels sit inline. */
        @media (max-width: 760px) {
            .tabs { grid-template-columns: 1fr; gap: 0; }
            .tabs label, .nav-group, .tab-content { grid-column: 1; }
            .tab-content { grid-row: auto; margin-top: 14px; }
            .tabs label {
                display: inline-block;
                border-left: none;
                border-bottom: 3px solid transparent;
                border-radius: 5px 5px 0 0;
            }
            .tabs label.sub { padding-left: 14px; }
            .tabs input[type="radio"]:checked + label {
                border-left-color: transparent;
                border-bottom-color: #2c5aa0;
            }
            .nav-group { padding: 12px 4px 2px; }
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
        // Jump from the BGC Novelty priority table to a family's representative card.
        // Defined here, not in viz/clustering.py, because the caller and the target
        // live in different sections: clustering.py's <script> is only emitted when
        // there are GCF cards to render, so a handler defined there is undefined
        // whenever a run has a novelty ranking and no representatives. The report
        // linter caught exactly that.
        //
        // The tabs are CSS radio buttons, so <a href="#gcf_7"> would scroll to an
        // element that is display:none and appear to do nothing — the radio has to be
        // checked first.
        function showGCF(familyId) {
            const card = document.getElementById('gcf_' + familyId);
            if (!card) return;                      // no representatives in this run
            const tab = document.getElementById('tab3');   // Gene Cluster Families > Analysis
            if (tab) tab.checked = true;
            const content = document.getElementById('gcf_' + familyId + '_content');
            const toggle = document.getElementById('gcf_' + familyId + '_toggle');
            if (content && content.style.display === 'none') {
                content.style.display = 'block';
                if (toggle) toggle.textContent = '-';
            }
            card.scrollIntoView({behavior: 'smooth', block: 'center'});
            card.style.transition = 'box-shadow .3s';
            card.style.boxShadow = '0 0 0 3px #2c5aa0';
            setTimeout(function () { card.style.boxShadow = ''; }, 1600);
        }

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

        // ---- Taxonomy tree: species genome lists rendered on first expand -------
        // Inlined, these were a second copy of all 1,735 genomes and ~90% of the tree's
        // 657 KB. Colour thresholds mirror get_green_bg_color / get_red_font_color in
        // viz/taxonomy.py — change both together.
        let _taxGenomes = null;

        function taxGenomes() {
            if (_taxGenomes === null) {
                const el = document.getElementById('taxonomyGenomeData');
                try { _taxGenomes = el ? JSON.parse(el.textContent) : {}; }
                catch (e) { _taxGenomes = {}; }
            }
            return _taxGenomes;
        }

        function greenBg(count, max) {
            if (!count) return '';
            const i = count / max;
            if (i <= 0.2) return '#e8f5e9';
            if (i <= 0.4) return '#a5d6a7';
            if (i <= 0.6) return '#66bb6a';
            if (i <= 0.8) return '#43a047';
            return '#2e7d32';
        }

        function redFont(count, max) {
            if (!count) return '';
            const i = count / max;
            if (i <= 0.2) return '#ffcdd2';
            if (i <= 0.4) return '#ef5350';
            if (i <= 0.6) return '#e53935';
            if (i <= 0.8) return '#c62828';
            return '#b71c1c';
        }

        function renderTaxonomyGenomes(container) {
            const nodeId = container.dataset.node;
            const rows = taxGenomes()[nodeId];
            if (!rows) { container.innerHTML = ''; return; }

            const cols = [];
            rows.forEach(function (r) {
                Object.keys(r[2] || {}).forEach(function (k) {
                    if (cols.indexOf(k) === -1) cols.push(k);
                });
            });
            cols.sort();

            let typeMax = 1, totalMax = 1;
            rows.forEach(function (r) {
                if (r[1] > totalMax) totalMax = r[1];
                Object.keys(r[2] || {}).forEach(function (k) {
                    if (r[2][k] > typeMax) typeMax = r[2][k];
                });
            });

            const head = ['<th style="text-align:left;">Genome</th>',
                          '<th>Total BGCs</th>'].concat(
                          cols.map(function (c) { return '<th>' + esc(c) + '</th>'; })).join('');
            const body = rows.map(function (r) {
                const cells = cols.map(function (c) {
                    const v = (r[2] || {})[c] || 0;
                    const col = redFont(v, typeMax);
                    return '<td style="text-align:center;' + (col ? 'color:' + col + ';font-weight:600;' : '') + '">'
                           + (v || '') + '</td>';
                }).join('');
                const bg = greenBg(r[1], totalMax);
                return '<tr><td><a href="genomes/' + esc(r[0]) + '.html">' + esc(r[0]) + '</a></td>'
                     + '<td style="text-align:center;' + (bg ? 'background:' + bg + ';' : '') + '">' + r[1] + '</td>'
                     + cells + '</tr>';
            }).join('');

            container.innerHTML = '<div class="genome-list"><table style="width:100%;border-collapse:collapse;font-size:0.9em;">'
                                + '<thead><tr>' + head + '</tr></thead><tbody>' + body + '</tbody></table></div>';
            container.dataset.rendered = '1';
        }

        function toggleNode(nodeId) {
            const element = document.getElementById(nodeId);
            const header = element.previousElementSibling;
            const icon = header.querySelector('.toggle-icon');
            if (element.style.display === 'none') {
                element.style.display = 'block';
                icon.innerHTML = '&#9660;';
                // render any deferred genome list the first time this node opens
                element.querySelectorAll('.genome-list-lazy:not([data-rendered])')
                       .forEach(renderTaxonomyGenomes);
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
