"""
Report Writer Agent Module

This module provides the Report Writer agent for generating comprehensive
scientific reports from quantum chemistry calculations and literature searches.
"""

from .report_writer_worker import create_report_writer

__all__ = ['create_report_writer']
