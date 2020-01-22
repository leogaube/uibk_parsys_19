COMPILER=chpl

mat_mult: mat_mult.chpl
	$(COMPILER) $< -o $%

.PHONEY: run
run: mat_mult
	qsub 8_threads.script
	qsub 4_threads.script
	qsub 2_threads.script
	qsub 1_thread.script


