#!/bin/bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <input_prefix> <output_directory>"
    exit 1
fi

PREFIX="$1"
OUTPUT_DIR="$2"

# --- Auto-detect omic_type from path ---
OMIC_TYPE=$(echo "$PREFIX" | grep -oP 'tensorqtl_outputs/\K[^/]+')
if [[ -z "$OMIC_TYPE" ]]; then
    echo "ERROR: Could not detect omic_type from path."
    exit 1
fi

# --- Auto-detect celltype and assay from filename ---
BASENAME=$(basename "$PREFIX")
NAME_PART="${BASENAME%.cis_qtl_pairs}"
CELLTYPE="${NAME_PART%%_*}"
ASSAY="${NAME_PART#*_}"

echo "=== Configuration ==="
echo "  Prefix:     $PREFIX"
echo "  Output dir: $OUTPUT_DIR"
echo "  Omic type:  $OMIC_TYPE"
echo "  Cell type:  $CELLTYPE"
echo "  Assay:      $ASSAY"
echo ""

# --- Find chromosome files ---
mapfile -t FILES < <(ls "${PREFIX}."*.csv 2>/dev/null | sort -t'.' -k2 -n)

if [[ ${#FILES[@]} -eq 0 ]]; then
    echo "ERROR: No files found matching ${PREFIX}.*.csv"
    exit 1
fi

echo "  ${#FILES[@]} chromosome files found"
echo "  First: $(basename "${FILES[0]}")"
echo "  Last:  $(basename "${FILES[-1]}")"
echo ""

# --- Check dependencies ---
if command -v pigz &>/dev/null; then
    COMPRESS="pigz -p 2"
else
    echo "WARNING: pigz not found, using gzip"
    COMPRESS="gzip"
fi

# --- Read header from first file ---
HEADER=$(head -1 "${FILES[0]}")

# --- Create output directory ---
mkdir -p "${OUTPUT_DIR}"
OUT_BASE="${OUTPUT_DIR%/}/${CELLTYPE}_${OMIC_TYPE}_${ASSAY}"

# --- Process via Python ---
python3 - "${HEADER}" "${OUT_BASE}" "${COMPRESS}" "${FILES[@]}" << 'PYEOF'
import sys, re, subprocess, signal, os, tempfile

signal.signal(signal.SIGPIPE, signal.SIG_DFL)

args = sys.argv[1:]
header = args[0].split(',')
out_base = args[1]
compress_cmd = args[2]
files = args[3:]

print(f"  Header: {len(header)} columns", file=sys.stderr)

# --- Column index mapping (0-based) ---
col_map = {c: i for i, c in enumerate(header)}

# --- Required common columns ---
common_names = ['phenotype_id', 'variant_id', 'af', 'ma_samples', 'ma_count']
for c in common_names:
    if c not in col_map:
        print(f"  ERROR: Required column '{c}' not found in header!", file=sys.stderr)
        sys.exit(1)
common_idx = [col_map[c] for c in common_names]

# --- Required main genotype effect columns ---
main_names = ['pval_g', 'b_g', 'b_g_se']
for c in main_names:
    if c not in col_map:
        print(f"  ERROR: Required column '{c}' not found in header!", file=sys.stderr)
        sys.exit(1)
main_idx = [col_map[c] for c in main_names]

# --- Discover modules from header ---
modules = {}
skip_set = set(common_names + main_names + ['start_distance', 'end_distance'])

for i, c in enumerate(header):
    if c in skip_set:
        continue

    m = re.search(r'(ME[a-zA-Z0-9]+)', c)
    if not m:
        continue
    mod = m.group(1)
    if mod not in modules:
        modules[mod] = {}

    if c.startswith('pval_g-'):
        modules[mod]['pval_int'] = i
    elif c.startswith('b_g-') and c.endswith('_se'):
        modules[mod]['b_int_se'] = i
    elif c.startswith('b_g-'):
        modules[mod]['b_int'] = i
    elif c.startswith('pval_'):
        modules[mod]['pval_me'] = i
    elif c.startswith('b_') and c.endswith('_se'):
        modules[mod]['b_me_se'] = i
    elif c.startswith('b_'):
        modules[mod]['b_me'] = i

# --- Validate modules ---
expected_keys = {'pval_me', 'b_me', 'b_me_se', 'pval_int', 'b_int', 'b_int_se'}
valid_modules = {}
for mod, cols in modules.items():
    if set(cols.keys()) == expected_keys:
        valid_modules[mod] = cols
    else:
        missing = expected_keys - set(cols.keys())
        print(f"  WARNING: Module {mod} missing {missing}. Skipping.", file=sys.stderr)

modules = valid_modules
mod_names_sorted = sorted(modules.keys())
print(f"  {len(modules)} valid modules detected", file=sys.stderr)
if len(modules) <= 30:
    print(f"  Modules: {' '.join(mod_names_sorted)}", file=sys.stderr)

if len(modules) == 0:
    print("  ERROR: No valid modules found!", file=sys.stderr)
    sys.exit(1)

# --- Build source stream command ---
parts = [f'cat "{files[0]}"']
for f in files[1:]:
    parts.append(f'tail -n +2 "{f}"')
src_cmd = '{ ' + '; '.join(parts) + '; }'

# --- Define output headers ---
out_header_main = r'phenotype_id\tvariant_id\taf\tma_samples\tma_count\tpval_g\tb_g\tb_g_se'
out_header_mod = r'phenotype_id\tvariant_id\taf\tma_samples\tma_count\tpval_g\tb_g\tb_g_se\tpval_ME\tb_ME\tb_ME_se\tpval\tb_int\tb_int_se'

# --- Create awk scripts and launch processes ---
temp_files = []
procs = []

def make_awk_script(field_indices, out_header):
    print_expr = ','.join(f'${x+1}' for x in field_indices)
    f = tempfile.NamedTemporaryFile(mode='w', suffix='.awk', delete=False)
    f.write('BEGIN { OFS="\\t" }\n')
    f.write(f'NR==1 {{ print "{out_header}"; next }}\n')
    f.write(f'{{ print {print_expr} }}\n')
    f.close()
    temp_files.append(f.name)
    return f.name

# Main effect file
main_fields = common_idx + main_idx
main_awk = make_awk_script(main_fields, out_header_main)
main_outfile = f'{out_base}_main_effect_trans_qtl.tsv.gz'
cmd = f'{src_cmd} | awk -F"," -f "{main_awk}" | {compress_cmd} > "{main_outfile}"'
procs.append(('main_effect', main_outfile, subprocess.Popen(cmd, shell=True, executable='/bin/bash')))

# Module files
for mod_name in mod_names_sorted:
    mod = modules[mod_name]
    mod_fields = common_idx + main_idx + [
        mod['pval_me'], mod['b_me'], mod['b_me_se'],
        mod['pval_int'], mod['b_int'], mod['b_int_se']
    ]
    mod_awk = make_awk_script(mod_fields, out_header_mod)
    mod_outfile = f'{out_base}_{mod_name}_trans_qtl.tsv.gz'
    cmd = f'{src_cmd} | awk -F"," -f "{mod_awk}" | {compress_cmd} > "{mod_outfile}"'
    procs.append((mod_name, mod_outfile, subprocess.Popen(cmd, shell=True, executable='/bin/bash')))

# --- Wait for all processes ---
print(f"  Launched {len(procs)} parallel pipelines...", file=sys.stderr)
print("", file=sys.stderr)

failures = 0
for label, fname, p in procs:
    rc = p.wait()
    bname = os.path.basename(fname)
    if rc != 0:
        print(f"  FAILED [{label}]: {bname} (exit code {rc})", file=sys.stderr)
        failures += 1
    else:
        try:
            size_mb = os.path.getsize(fname) / (1024 * 1024)
            print(f"  OK [{label}]: {bname} ({size_mb:.1f} MB)", file=sys.stderr)
        except OSError:
            print(f"  OK [{label}]: {bname}", file=sys.stderr)

# --- Clean up ---
for tf in temp_files:
    try:
        os.unlink(tf)
    except OSError:
        pass

print("", file=sys.stderr)
if failures > 0:
    print(f"  ERROR: {failures} process(es) failed!", file=sys.stderr)
    sys.exit(1)

print(f"  SUCCESS: All {len(procs)} files written.", file=sys.stderr)
PYEOF

echo ""
echo "Done: ${CELLTYPE} / ${ASSAY} / ${OMIC_TYPE}"
echo "Output: ${OUT_BASE}_*_trans_qtl.tsv.gz"