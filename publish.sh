#!/usr/bin/env bash
# ============================================================
#  One-click publish to GitHub
#  Edit the variables below, then run:  bash publish.sh
# ============================================================
set -e

REPO_NAME="plasma-immunome-phenome-atlas"
VISIBILITY="public"          # public | private
RELEASE_TAG="v1.0.0"

cd "$(dirname "$0")"
echo
echo "=== Publishing \"$REPO_NAME\" ($VISIBILITY) ==="
echo

# 1. init if needed
if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    echo "Initializing git repository..."
    git init -b main
fi

# 2. stage + commit
git add -A
if ! git diff --cached --quiet; then
    git commit -m "Release v1.0.0: open genetics-anchored plasma immunome atlas"
else
    echo "Nothing new to commit."
fi

# 3. create remote + push + release
if command -v gh >/dev/null 2>&1; then
    if ! gh auth status >/dev/null 2>&1; then
        echo "You are not logged into gh. Run:  gh auth login"; exit 1
    fi
    if git remote get-url origin >/dev/null 2>&1; then
        git push -u origin main
    else
        gh repo create "$REPO_NAME" --"$VISIBILITY" --source=. --remote=origin --push
    fi
    echo "Creating release $RELEASE_TAG..."
    gh release create "$RELEASE_TAG" \
        --title "$REPO_NAME $RELEASE_TAG" \
        --notes "Initial public release: four-layer causal-evidence plasma immunome atlas (9 main + 70 supplementary figures, evidence-tiered targets, manuscript)."
    echo
    echo "DONE. Repository and release published."
else
    echo "GitHub CLI (gh) not found."
    if git remote get-url origin >/dev/null 2>&1; then
        echo "Pushing to existing origin remote..."
        git push -u origin main
        echo "DONE. Pushed to origin."
    else
        cat <<EOF

No 'origin' remote set. To finish, either:
  A) Install GitHub CLI (https://cli.github.com), then re-run this script, OR
  B) Create an empty repo on github.com and run:
       git remote add origin https://github.com/<you>/$REPO_NAME.git
       git push -u origin main
EOF
    fi
fi
