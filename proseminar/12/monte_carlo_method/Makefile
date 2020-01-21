COMPILER=chpl

SOURCE=pi.chpl

pi: pi.chpl
	$(COMPILER) $(SOURCE) -o pi

.PHONEY: run
run: pi
	qsub 1_thread.script
	qsub 2_thread.script
	qsub 4_thread.script
	qsub 8_thread.script


