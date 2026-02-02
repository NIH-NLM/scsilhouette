# src/scsilhouette/logging_config.py

import logging
import sys

def setup_logger(name: str = "scsilhouette") -> logging.Logger:
    """Configure logger for scsilhouette"""
    
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    
    # Avoid duplicate handlers
    if logger.handlers:
        return logger
    
    # Console handler
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.INFO)
    
    # Simple format matching NSForest style
    formatter = logging.Formatter('%(message)s')
    handler.setFormatter(formatter)
    
    logger.addHandler(handler)
    
    return logger
