#!/bin/bash

echo "=== SNMAKE RULES BY .SMK FILE ==="
for smk_file in ./workflow/rules/*.smk; do
    if [[ -f "$smk_file" ]]; then
        echo
        echo "📁 $smk_file:"
        # Extract rule names using grep (look for 'rule ' followed by word)
        grep -oE 'rule [a-zA-Z0-9_]+:' "$smk_file" | sed 's/:$//' | sort
    fi
done

echo
echo "=== SCRIPTS IN 'scripts/' DIRECTORY ==="
ls -1 scripts/
