#!/usr/bin/env bash
# -*- coding:  utf-8 -*-
# Author: Gaurav Sablok
# date: 2023-10-19

echo " updating the genome assembly"
echo "if you have an already assembled genome then use the 
        update option or else use the genome assembly to update"

read -r -p "please select the choice of the option:": option
read -r -p "please select the choice of the assembler:": assembler