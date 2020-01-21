COMPILER=chpl

SOURCE=pi.chpl

pi: pi.chpl
	$(COMPILER) $(SOURCE) -o pi

.PHONEY: run
run: pi
	qsub 1_thread.script
	qsub 2_threads.script
	qsub 4_threads.script
	qsub 8_threads.script


