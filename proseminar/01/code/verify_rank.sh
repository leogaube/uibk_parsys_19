#!/bin/bash

/bin/hostname
cat /proc/self/stat | awk '{print $39}'
