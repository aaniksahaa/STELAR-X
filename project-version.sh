#!/usr/bin/env bash
# Print the STELAR-X project version declared in pom.xml.
# All repository build and launcher scripts use this helper so pom.xml remains
# the single source of truth for release versions and artifact names.
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
POM_FILE="$PROJECT_ROOT/pom.xml"

if [[ ! -f "$POM_FILE" ]]; then
  echo "Error: Maven project file not found at $POM_FILE" >&2
  exit 1
fi

VERSION="$(sed -nE 's/^[[:space:]]*<version>([^<]+)<\/version>[[:space:]]*$/\1/p' "$POM_FILE" | head -n 1)"

if [[ -z "$VERSION" || "$VERSION" == *[!0-9A-Za-z._-]* ]]; then
  echo "Error: Could not read a valid project version from $POM_FILE" >&2
  exit 1
fi

printf '%s\n' "$VERSION"
