#!/bin/bash
set -euo pipefail

Rscript .devcontainer/setup.R
R --version