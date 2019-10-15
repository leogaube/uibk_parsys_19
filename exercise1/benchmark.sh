#!/bin/bash

qsub latency_1h1s.script
qsub latency_1h2s.script
qsub latency_2h.script

qsub bandwidth_1h1s.script
qsub bandwidth_1h2s.script
qsub bandwidth_2h.script
