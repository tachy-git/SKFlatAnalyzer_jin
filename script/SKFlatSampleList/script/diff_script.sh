#!/bin/sh 
diffout="$(diff $1 <(find $2 -type f|sort -V))"
echo $diffout >> $4/tmp$3.txt
