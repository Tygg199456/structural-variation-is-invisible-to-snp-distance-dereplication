#!/usr/bin/env bash
# Replace compute-host absolute paths with a placeholder so the archive is portable.
# Run on an archive tree:  bash scripts/scrub_paths.sh <tree>
#   SCRUB_HOSTS=<alias>,<alias>  -> also replace those host aliases with <compute-host>
set -euo pipefail
TREE="${1:?usage: scrub_paths.sh <tree>}"
cd "$TREE"
SCRUB_HOSTS="${SCRUB_HOSTS:-}" python3 - <<'PY'
import os, re
HOSTS = [h.strip() for h in os.environ.get("SCRUB_HOSTS", "").split(",") if h.strip()]
root = os.getcwd()
for dirpath, dirnames, filenames in os.walk(root):
    if '.git' in dirpath.split(os.sep):
        continue
    for fn in filenames:
        p = os.path.join(dirpath, fn)
        if fn.endswith(('.npy', '.npz', '.png', '.pdf', '.gif', '.jpg')):
            continue
        try:
            s = open(p, encoding='utf-8', errors='strict').read()
        except (UnicodeDecodeError, OSError):
            continue
        o = s
        s = re.sub(r'ROOT\s*=\s*"/home/[^"]*/work/SY/pipeline/plague-decay"',
                   'ROOT = os.environ.get("PLAGUE_DECAY_ROOT",\n'
                   '                       os.path.dirname(os.path.dirname(os.path.abspath(__file__))))', s)
        s = re.sub(r'/home/[A-Za-z0-9_.-]+/envs/([A-Za-z0-9_.-]+)', r'<CONDA_ENV_PREFIX>/\1', s)
        s = re.sub(r'/home/[A-Za-z0-9_.-]+', '<HOME>', s)
        s = re.sub(r'/Users/[A-Za-z0-9_.-]+', '<HOME>', s)
        # compute-host aliases (supply your own: SCRUB_HOSTS=a,b,c) and the project dir
        for h in HOSTS:
            if h:
                s = s.replace(h, '<compute-host>')
        s = re.sub(r'~/Desktop/SY01|\$HOME/Desktop/SY01', '<project-dir>', s)
        if s != o:
            open(p, 'w', encoding='utf-8').write(s)
            print('scrubbed', os.path.relpath(p, root))
PY
echo "residual absolute home paths: $(grep -rIl '/Users/\|/home/' . 2>/dev/null | grep -v '^./.git/' | wc -l)"
