#!/usr/bin/env bash

## autobuild.sh: automatically rebuild mansucript outputs and the webpage when content changes
## Depends on watchdog https://github.com/gorakhargosh/watchdog

watchmedo shell-command \
  --wait \
  --command='bash paper/build/build.sh && python paper/build/webpage.py' \
  content
