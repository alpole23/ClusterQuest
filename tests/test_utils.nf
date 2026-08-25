/*
 * Unit checks for lib/Utils.groovy helpers that guard optional pipeline inputs.
 * Prints PASS/FAIL lines; run via tests/run_tests.sh.
 */

nextflow.enable.dsl=2

def check(String name, actual, expected) {
    def ok = (actual == expected)
    println "${ok ? 'PASS' : 'FAIL'} ${name}  (got '${actual}', expected '${expected}')"
    if (!ok) {
        System.exit(1)
    }
}

workflow {
    def real = file("${params.fixtures}/name_map.json")
    def ph   = file('NO_COUNTS')

    check('optArg real file',   Utils.optArg('--counts', real),  "--counts ${real}")
    check('optArg placeholder', Utils.optArg('--counts', ph),    '')
    check('optArg list real',   Utils.optArg('--db', [real]),    "--db ${real}")
    check('optArg list ph',     Utils.optArg('--db', [ph]),      '')
    check('optArg empty list',  Utils.optArg('--db', []),        '')
    check('optArg null',        Utils.optArg('--db', null),      '')

    check('isValidInput real',  Utils.isValidInput(real),        true)
    check('isValidInput ph',    Utils.isValidInput(ph),          false)

    check('sanitizeTaxon',      Utils.sanitizeTaxon('Pantoea ananatis'), 'Pantoea_ananatis')
    check('sanitizeTaxon punct', Utils.sanitizeTaxon('E. coli (K-12)'),  'E_coli_K_12')

    // batchSize() must coerce: params from the command line arrive as strings and
    // collate() silently fails to dispatch on anything but a real Integer.
    check('batchSize is Integer', (params.task_batch_size.toString().toInteger()) instanceof Integer, true)
}
