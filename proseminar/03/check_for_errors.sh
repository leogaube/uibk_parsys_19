#!/bin/bash

if [ "$(ls -A $1)" ]; then
	for filename in "$1/*.err"; do
		cat filename;
	done
else
	echo "$1 is empty!"
fi