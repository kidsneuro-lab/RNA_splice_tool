#!/usr/bin/env bash

set -euo pipefail

Rscript -e "devtools::test()"
