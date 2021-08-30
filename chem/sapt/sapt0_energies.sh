#!/usr/bin/env bash
find . -maxdepth 2 -name "*log" | xargs grep "Total sSAPT0" | sort -k8n | tr -s ' ' | cut -d ' ' -f1,8,9 | column -t
