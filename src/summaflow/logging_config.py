# Built-in imports
import logging
import sys

# Logging configuration module
def setup_logger(name=__name__, level=logging.INFO):
    """Set up and return a configured logger."""
    logger = logging.getLogger(name)
    logger.setLevel(level)
 
    # Prevent adding handlers multiple times
    if not logger.handlers:
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )

        # Console handler
        ch = logging.StreamHandler(sys.stdout)
        ch.setFormatter(formatter)
        logger.addHandler(ch)

    return logger
