import os
import shutil
import logging

def setup_logger(name: str = "CFDSolver") -> logging.Logger:
    """
    Configures and returns a basic logger for the application.
    """
    logger = logging.getLogger(name)
    if not logger.handlers:
        logger.setLevel(logging.INFO)
        handler = logging.StreamHandler()
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%H:%M:%S')
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    return logger

def clean_results_dir(results_dir: str, logger: logging.Logger = None) -> None:
    """
    Cleans the results directory to ensure no legacy data overlaps with new runs.
    
    Args:
        results_dir (str): Path to the results directory.
        logger (logging.Logger): Optional logger to log the cleaning process.
    """
    if os.path.exists(results_dir):
        if logger:
            logger.info(f"Cleaning results directory: {results_dir}")
        for filename in os.listdir(results_dir):
            file_path = os.path.join(results_dir, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                if logger:
                    logger.error(f"Failed to delete {file_path}. Reason: {e}")
    else:
        if logger:
            logger.info(f"Creating results directory: {results_dir}")
        os.makedirs(results_dir, exist_ok=True)
