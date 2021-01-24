#!/bin/bash

for dirs in $(ls --color=never -l | grep "^d" | awk '{print $9}'); do du -hs $dirs;done
