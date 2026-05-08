"""Standardized logging for KwanTube v3.5.1."""
import logging
from datetime import datetime, timezone

def get_auditor_logger(name: str) -> logging.Logger:
    """Returns a logger with ISO-8601 UTC timestamp format.
    
    Args:
        name: Logger name (typically script filename)
    
    Returns:
        Configured logger instance with UTC timestamps
    """
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler()
        formatter = logging.Formatter(
            "[%(asctime)s] [%(filename)s] %(message)s",
            datefmt="%Y-%m-%dT%H:%M:%S%z"
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
    return logger
