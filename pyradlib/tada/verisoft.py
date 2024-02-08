# -*- coding: utf-8 -*-
""" Verisoft patient-specific quality assurance report reader.

This module provides parsing functionality for patient-specific quality assurance reports at the Royal Brisbane & Women's Hospital.
"""

# authorship information
__author__ = 'Scott Crowe'
__email__ = 'sb.crowe@gmail.com'
__license__ = "GPL3"

# import required code
import os
import pandas as pd
import pdfplumber
import re


def analyse_verisoft_report(path: str):
    """ Read contents of a PTW Verisoft report generated at the RBWH. 

    Args: 
        path (str): The path of the PTW report file to be opened.
    """
    _non_decimal = re.compile(r'[^\d.]+')
    with pdfplumber.open(path) as report:            
        text = report.pages[0].extract_text().split('\n')
        if len(text) > 17:
            patient_id = text[3].replace('PatientID ','')
            patient_name = text[4].replace('Patient Name ','').rstrip(',').replace(', ','^').rstrip('^')
            comment = text[5].replace('Comment ','')
            comment_number = re.findall(r'\d+', comment)
            if len(comment_number) > 0:
                # if number is found, assume it is for an arc or beam
                if 'arc' in comment.lower():
                    beam_name = comment_number[0].zfill(2) + 'ARC'
                    plan_name = re.sub(r'\d+','',comment).lower().replace(':','').replace('arc','').replace('  ',' ').strip().upper()
                else:
                    beam_name = comment_number[0].zfill(2) + 'BEAM'
                    plan_name = re.sub(r'\d+','',comment).lower().replace(':','').replace('arc','').replace('  ',' ').strip().upper()
            else:
                beam_name = 'TOTAL'
                plan_name = comment.lower().replace(':','').replace('total','').replace('  ',' ').strip().upper()
            gamma_dta = float(_non_decimal.sub('',text[7]))
            gamma_dd = float(_non_decimal.sub('',text[8]).rstrip('.'))
            if text[16].startswith('Result'):
                gamma_pass_rate = float(_non_decimal.sub('',text[16]))
            elif text[17].startswith('Result'):
                gamma_pass_rate = float(_non_decimal.sub('',text[17]))
            gamma_pass_rate = min(gamma_pass_rate, 100)
            return patient_id, patient_name, plan_name, beam_name, gamma_dd, gamma_dta, gamma_pass_rate


def analyse_verisoft_report_archive(path: str):
    """ Read contents of multiple RBWH PTW Verisoft reports at the RBWH. 

    Args: 
        path (str): The path of the folder to be walked through.
    """
    qa_results = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if '.pdf' in name.lower():
                qa_results.append(analyse_verisoft_report(os.path.join(root,name)))
    return pd.DataFrame(qa_results, columns = ['Patient ID', 'Patient Name', 'Plan Name', 'Beam Name', 'Dose Criteria', 'DTA Criteria', 'Gamma Pass Rate'])