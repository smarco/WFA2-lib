#!/bin/bash
awk '{if($0~/^>/){printf("\n>")}else{printf($0)}}' | tail -n +2
