COMPILER=chpl

mat_mult: mat_mult.chpl
	$(COMPILER) mat_mult.chpl -o mat_mult

.PHONEY: run
run: mat_mult
	qsub 8_threads.script
	qsub 4_threads.script
	qsub 2_threads.script
	qsub 1_thread.script


