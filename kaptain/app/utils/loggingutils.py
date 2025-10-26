import logging


logger = logging.getLogger('KAPTAIN')


def initialize_logging(activate_log: bool = False) -> None:
    """
    Initializes the logging.
    :return: None
    """
    formatter = logging.Formatter('%(asctime)s - %(module)15s - %(levelname)7s - %(message)s')

    # Remove all existing handlers
    if logger.hasHandlers():
        logger.handlers.clear()

    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    console_handler.setLevel(logging.DEBUG)
    console_handler.name = 'console'
    logger.addHandler(console_handler)

    if activate_log:
        file_handler = logging.FileHandler('kt.log')
        file_handler.setFormatter(formatter)
        file_handler.setLevel(logging.DEBUG)
        logger.addHandler(file_handler)

    # General logging level
    logger.setLevel(logging.DEBUG)
    logger.propagate = False