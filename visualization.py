#!/usr/bin/env python
# a streamlit application for the genome visualization
# -*- coding:  utf-8 -*-
# Author: Gaurav Sablok
# date: 2023-10-20
import streamlit as st
import jbrowse 
import pysam 

class Browsertracts:
    def __init__(self, bam, )