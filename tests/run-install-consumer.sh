#!/usr/bin/env bash
# Runner for tests/install_consumer (PR #73 review B-Test1).
#
# Installs Exasim from $EXASIM_BUILD into $INSTALL_PREFIX, configures
# the install_consumer against that prefix, builds it, runs it. Exit
# code 0 means the install layout (B1) and public-header ODR (B2)
# are both intact.
#
# Usage:  bash tests/run-install-consumer.sh
# Env:    EXASIM_BUILD     in-tree build dir            (default: $REPO/build)
#         KOKKOS_DIR       Kokkos install/build dir     (default: $REPO/kokkos/buildserial)
#         INSTALL_PREFIX   where to install Exasim      (default: /tmp/exasim_install)
#         CONSUMER_BUILD   where to build the consumer  (default: /tmp/install_consumer_build)

set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
EXASIM_BUILD="${EXASIM_BUILD:-$REPO/build}"
KOKKOS_DIR="${KOKKOS_DIR:-$REPO/kokkos/buildserial}"
INSTALL_PREFIX="${INSTALL_PREFIX:-/tmp/exasim_install}"
CONSUMER_BUILD="${CONSUMER_BUILD:-/tmp/install_consumer_build}"

echo "[B-Test1] Repo            $REPO"
echo "[B-Test1] Exasim build    $EXASIM_BUILD"
echo "[B-Test1] Install prefix  $INSTALL_PREFIX"
echo "[B-Test1] Consumer build  $CONSUMER_BUILD"
echo "[B-Test1] Kokkos          $KOKKOS_DIR"

if [ ! -d "$EXASIM_BUILD" ]; then
    echo "[B-Test1] FAIL: Exasim build dir does not exist: $EXASIM_BUILD" >&2
    exit 2
fi

echo "[B-Test1] Installing Exasim..."
rm -rf "$INSTALL_PREFIX"
cmake --install "$EXASIM_BUILD" --prefix "$INSTALL_PREFIX" > /dev/null

echo "[B-Test1] Configuring consumer..."
rm -rf "$CONSUMER_BUILD"
cmake -S "$REPO/tests/install_consumer" -B "$CONSUMER_BUILD" \
      -D "CMAKE_PREFIX_PATH=$INSTALL_PREFIX;$KOKKOS_DIR" > /dev/null

echo "[B-Test1] Building consumer..."
cmake --build "$CONSUMER_BUILD" > /dev/null

echo "[B-Test1] Running consumer..."
"$CONSUMER_BUILD/install_consumer"

echo "[B-Test1] PASS"
