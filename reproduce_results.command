#!/usr/bin/env bash
DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$DIR"
chmod +x ./reproduce_results.sh
./reproduce_results.sh
