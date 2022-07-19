#!/bin/sh 
diffout="$(diff $1 $2)"
echo $diffout
echo $diffout > $3/tmp.txt
