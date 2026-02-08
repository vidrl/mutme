#!/usr/bin/env sh
set -eu

FILE="pyproject.toml"
MSG_FILE=".git/COMMIT_EDITMSG"

# If commit message doesn't exist yet (rare, but happens), do nothing
[ -f "$MSG_FILE" ] || exit 0

HEADER="$(head -n1 "$MSG_FILE" | tr -d '\r')"

# Decide bump type
BUMP=""
echo "$HEADER" | grep -Eq '^(fix|hotfix)(\(.+\))?:' && BUMP="patch"
echo "$HEADER" | grep -Eq '^feat(\(.+\))?:' && BUMP="minor"

[ -n "$BUMP" ] || exit 0

# Extract version from [project]
OLD_VERSION="$(
  awk '
    BEGIN { in_project=0 }
    /^\[project\][[:space:]]*$/ { in_project=1; next }
    /^\[[^]]+\][[:space:]]*$/ { in_project=0 }
    in_project && $0 ~ /^[[:space:]]*version[[:space:]]*=/ {
      if (match($0, /"([0-9]+\.[0-9]+\.[0-9]+)"/, m)) {
        print m[1]; exit
      }
    }
  ' "$FILE"
)"

if [ -z "$OLD_VERSION" ]; then
  echo "ERROR: Could not find [project].version in $FILE" >&2
  exit 1
fi

MAJOR=$(echo "$OLD_VERSION" | cut -d. -f1)
MINOR=$(echo "$OLD_VERSION" | cut -d. -f2)
PATCH=$(echo "$OLD_VERSION" | cut -d. -f3)

case "$BUMP" in
  patch)
    PATCH=$((PATCH + 1))
    ;;
  minor)
    MINOR=$((MINOR + 1))
    PATCH=0
    ;;
esac

NEW_VERSION="${MAJOR}.${MINOR}.${PATCH}"

# Rewrite only the [project] version
TMP="$(mktemp)"
awk -v old="$OLD_VERSION" -v new="$NEW_VERSION" '
  BEGIN { in_project=0; replaced=0 }
  /^\[project\][[:space:]]*$/ { in_project=1; print; next }
  /^\[[^]]+\][[:space:]]*$/ { in_project=0; print; next }

  in_project && !replaced && $0 ~ /^[[:space:]]*version[[:space:]]*=/ {
    gsub("\"" old "\"", "\"" new "\"")
    replaced=1
    print
    next
  }

  { print }
' "$FILE" > "$TMP"

mv "$TMP" "$FILE"

# Ensure the bumped version is included in the commit
git add "$FILE"

echo "Pre-commit: bumped version $OLD_VERSION -> $NEW_VERSION ($BUMP)"
