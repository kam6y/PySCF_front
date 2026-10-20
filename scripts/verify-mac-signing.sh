#!/bin/bash
# Verify that the packaged macOS app is correctly code-signed and, for release
# builds, notarized and stapled. Run it after `electron-builder --mac`.
#
# Usage:
#   scripts/verify-mac-signing.sh [--allow-development] [path/to/App.app ...]
#
# Without explicit paths every .app under release/mac*/ is checked.
#
# --allow-development   Accept an "Apple Development" signing identity and skip
#                       the Gatekeeper / notarization-ticket checks. Only for
#                       local tests on a machine without a Developer ID
#                       certificate; a build that needs this flag must not be
#                       distributed.
set -euo pipefail

ALLOW_DEVELOPMENT=0
APPS=()
for arg in "$@"; do
  case "$arg" in
    --allow-development) ALLOW_DEVELOPMENT=1 ;;
    -h|--help) sed -n '2,15p' "$0"; exit 0 ;;
    *) APPS+=("$arg") ;;
  esac
done

if [ ${#APPS[@]} -eq 0 ]; then
  while IFS= read -r -d '' app; do
    APPS+=("$app")
  done < <(find release -maxdepth 2 -name '*.app' -print0 2>/dev/null)
fi

if [ ${#APPS[@]} -eq 0 ]; then
  echo "✗ No .app bundle found under release/. Run 'npx electron-builder --mac' first." >&2
  exit 1
fi

fail() { echo "✗ $*" >&2; exit 1; }
ok()   { echo "✓ $*"; }

for APP in "${APPS[@]}"; do
  echo "=== Verifying $APP ==="
  [ -d "$APP" ] || fail "Not a directory: $APP"

  # 1. Whole-bundle signature, deep + strict (what Gatekeeper evaluates)
  codesign --verify --deep --strict --verbose=2 "$APP" 2>&1 | sed 's/^/    /' || true
  codesign --verify --deep --strict "$APP" || fail "codesign --verify failed for $APP"
  ok "Bundle signature is valid"

  # 2. Signing identity, Hardened Runtime and secure timestamp on the main executable
  details=$(codesign -dvv "$APP" 2>&1)
  echo "$details" | grep -E '^(Authority|TeamIdentifier|Timestamp)=' | sed 's/^/    /'
  if echo "$details" | grep -q '^Authority=Developer ID Application:'; then
    ok "Signed with a Developer ID Application certificate"
  elif [ "$ALLOW_DEVELOPMENT" -eq 1 ] && echo "$details" | grep -q '^Authority=Apple Development:'; then
    echo "⚠ Signed with an Apple Development certificate (accepted by --allow-development; NOT distributable)"
  else
    fail "Main executable is not signed with a Developer ID Application certificate"
  fi
  echo "$details" | grep -q 'flags=.*runtime' || fail "Hardened Runtime flag is missing on the main executable"
  ok "Hardened Runtime is enabled on the main executable"
  echo "$details" | grep -q '^Timestamp=' || fail "Secure timestamp is missing (notarization requires it)"
  ok "Secure timestamp present"

  # 3. Every Mach-O inside the bundle (Electron helpers, native modules and the
  #    whole conda environment under Resources/conda_env) must carry a real
  #    signature with the Hardened Runtime; notarization rejects the upload otherwise.
  echo "Scanning Mach-O binaries inside the bundle (includes conda_env)..."
  # `file` prints one extra line per slice for universal binaries
  # ("... (for architecture x86_64): application/x-mach-binary"); strip that
  # suffix and de-duplicate so every path is checked exactly once.
  machos=$(find "$APP" -type f -print0 \
    | xargs -0 -n 500 file --mime-type -- 2>/dev/null \
    | grep 'application/x-mach-binary$' \
    | sed -E 's/( \(for architecture [^)]*\))?:[[:space:]]*application\/x-mach-binary$//' \
    | sort -u || true)
  total=0
  bad=0
  while IFS= read -r bin; do
    [ -n "$bin" ] || continue
    total=$((total + 1))
    if ! info=$(codesign -dv "$bin" 2>&1); then
      echo "  ✗ unsigned:            ${bin#"$APP"/}"
      bad=$((bad + 1))
      continue
    fi
    if echo "$info" | grep -q '^Signature=adhoc'; then
      echo "  ✗ ad-hoc signature:    ${bin#"$APP"/}"
      bad=$((bad + 1))
      continue
    fi
    if ! echo "$info" | grep -q 'flags=.*runtime'; then
      echo "  ✗ no hardened runtime: ${bin#"$APP"/}"
      bad=$((bad + 1))
    fi
  done <<< "$machos"
  [ "$total" -gt 0 ] || fail "No Mach-O binaries found in $APP (is 'file' available?)"
  [ "$bad" -eq 0 ] || fail "$bad of $total Mach-O binaries are not properly signed"
  ok "All $total Mach-O binaries are signed with the Hardened Runtime"

  # 4. Notarization ticket and Gatekeeper assessment
  if [ "$ALLOW_DEVELOPMENT" -eq 1 ]; then
    echo "⚠ Skipping stapler validation and Gatekeeper assessment (--allow-development)"
  else
    xcrun stapler validate "$APP" || fail "No notarization ticket is stapled to $APP"
    ok "Notarization ticket is stapled"
    spctl --assess --type execute --verbose=4 "$APP" 2>&1 | sed 's/^/    /'
    spctl --assess --type execute "$APP" || fail "Gatekeeper rejected $APP"
    ok "Gatekeeper accepts the app"
  fi
done

echo "All macOS signing checks passed."
