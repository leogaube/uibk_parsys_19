CC=gcc
CC_FLAGS=-O2 -std=c99 -Wall

# Verification
.PHONEY: verify
verify: verification.c
	@$(CC) $(CC_FLAGS) $< -o verification
	@./verification


# clean
.PHONEY: clean
clean:
	@rm -f verification
